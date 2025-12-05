#!/usr/bin/env python3
"""
BLAST+ Genome Search Script
Searches a nucleotide or protein query FASTA file against a genome sequence from a specific accession number

This script follows the BLAST+ workflow:
1. Downloads genome FASTA for the accession (using coral/ncbi/download.py)
2. Creates BLAST+ nucleotide database from genome sequence (cached)
3. Performs BLAST search (blastn for DNA queries, tblastn for protein queries) using the created database
4. Outputs results to specified file
5. Cleans up temporary files while preserving the BLAST database

Usage: python blast-genome-search.py <query_fasta> <accession> <output_file>
"""

import argparse
import os
import sys
import subprocess
import shutil
import tempfile
import uuid
from pathlib import Path

from needle.ncbi import download_and_extract_by_accession
from needle.blast import Results


def detect_query_type(fasta_file):
    """
    Detect if the query FASTA file contains DNA or protein sequences.
    
    Args:
        fasta_file: Path to query FASTA file
    
    Returns:
        str: 'protein' or 'dna'
    """
    from Bio import SeqIO
    from Bio.SeqUtils import IUPACData
    
    protein_letters = set(IUPACData.protein_letters)
    dna_letters = {'A', 'T', 'G', 'C', 'N'}
    
    protein_count = 0
    dna_count = 0
    
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq).upper()
            if not seq:
                continue
            
            # Check first 100 characters to determine type
            sample = seq[:100]
            sample_letters = set(sample)
            
            # If we see protein-specific letters (not in DNA), it's protein
            if sample_letters & (protein_letters - dna_letters):
                protein_count += 1
            else:
                # Check if it looks like DNA (mostly A, T, G, C, N)
                non_dna = sample_letters - dna_letters
                if len(non_dna) / len(sample_letters) < 0.1 if sample_letters else True:
                    dna_count += 1
                else:
                    protein_count += 1
            
            # Only check first sequence for speed
            break
    except Exception:
        # Default to DNA if we can't determine
        return 'dna'
    
    return 'protein' if protein_count > dna_count else 'dna'

def run_docker_command(cmd_args, work_dir):
    """
    Run a Docker command with NCBI BLAST+.
    
    Args:
        cmd_args: List of command arguments for BLAST+
        work_dir: Working directory to mount as /blast/blastdb_custom
    
    Returns:
        subprocess.CompletedProcess: Result of the command execution
    """
    # Convert to absolute path for Docker volume mounting
    abs_work_dir = os.path.abspath(work_dir)
    
    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{abs_work_dir}:/blast/blastdb_custom",
        "ncbi/blast"
    ] + cmd_args
    
    print(f"Running: {' '.join(docker_cmd)}")
    result = subprocess.run(docker_cmd, capture_output=True, text=True, cwd=work_dir)
    
    if result.returncode != 0:
        print(f"Error running BLAST+ command: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, docker_cmd, result.stdout, result.stderr)
    
    return result

def check_blast_db_exists(db_path):
    """
    Check if a BLAST+ database already exists.
    
    Args:
        db_path: Path to the database (without extension)
    
    Returns:
        bool: True if database exists, False otherwise
    """
    # BLAST+ creates multiple files with different extensions
    required_files = [f"{db_path}.nhr", f"{db_path}.nin", f"{db_path}.nsq"]
    
    # Check if core database files exist
    return all(os.path.exists(f) for f in required_files)

def create_blast_database(fasta_path, db_path, work_dir):
    """
    Create a BLAST+ database from a FASTA file.
    
    Args:
        fasta_path: Path to input FASTA file (relative to work_dir)
        db_path: Path for output database (relative to work_dir, without extension)
        work_dir: Working directory for Docker operations
    """
    # print(f"Creating BLAST+ database from {fasta_path}...")
    
    # Paths need /blast/blastdb_custom/ prefix since that's how they appear in the Docker container
    cmd_args = ["makeblastdb", "-in", f"/blast/blastdb_custom/{fasta_path}", "-dbtype", "nucl", "-out", f"/blast/blastdb_custom/{db_path}"]
    run_docker_command(cmd_args, work_dir)
    
    # Verify database was created
    db_files = [f"{db_path}.nhr", f"{db_path}.nin", f"{db_path}.nsq"]
    for db_file in db_files:
        full_path = os.path.join(work_dir, db_file)
        if os.path.exists(full_path):
            # print(f"✓ Database file created: {db_file}")
            pass
        else:
            print(f"Warning: Database file missing: {db_file}")
    
    # print(f"✓ BLAST+ database created: {db_path}")

def run_blast_search(query_fasta, db_path, output_file, work_dir, evalue, query_type, min_word_size):
    """
    Run BLAST+ search (blastn for DNA, tblastn for protein).
    
    Args:
        query_fasta: Path to query FASTA file (relative to work_dir)
        db_path: Path to target database (relative to work_dir, without extension) - this is REUSED
        output_file: Path for output results (relative to work_dir)
        work_dir: Working directory for Docker operations
        evalue: E-value threshold for BLAST search
        query_type: 'dna' or 'protein' to determine which BLAST program to use
        min_word_size: Minimum word size for small matches (None for default)
    """
    is_protein = (query_type == 'protein')
    search_type = "tblastn" if is_protein else "blastn"
    
    # print(f"Running BLAST+ {search_type} search...")
    
    # Create a temporary directory for the search workflow
    temp_dir_name = f"blast_tmp_{uuid.uuid4().hex[:8]}"
    temp_dir_path = os.path.join(work_dir, temp_dir_name)
    os.makedirs(temp_dir_path, exist_ok=True)
    
    try:
        # Copy query FASTA to the temp directory for Docker operations
        query_fasta_dest = os.path.join(temp_dir_path, "query.fasta")
        shutil.copy2(os.path.join(work_dir, query_fasta), query_fasta_dest)
        
        # Build BLAST command
        cmd_args = [
            search_type,
            "-query", f"/blast/blastdb_custom/{temp_dir_name}/query.fasta",
            "-db", f"/blast/blastdb_custom/{db_path}",
            "-out", f"/blast/blastdb_custom/{output_file}",
            "-outfmt", "6 qseqid sseqid evalue pident qstart qend stitle sseq sstart send",
            "-evalue", str(evalue)
        ]
        
        # Add word size parameter for small matches
        if min_word_size is not None:
            if is_protein:
                # For tblastn, word_size is in amino acids (default 3)
                # Allow down to 2 amino acids (6 bp when translated)
                word_size = max(2, int(min_word_size))
                cmd_args.extend(["-word_size", str(word_size)])
            else:
                # For blastn, word_size is in nucleotides (default 11)
                # Allow down to 4 nucleotides for very small matches
                word_size = max(4, int(min_word_size))
                cmd_args.extend(["-word_size", str(word_size)])

        run_docker_command(cmd_args, work_dir)
        
        # print(f"✓ BLAST+ {search_type} search completed successfully")
        
    finally:
        # Clean up temporary files
        # print("Cleaning up search temporary files...")
        try:
            if os.path.exists(temp_dir_path):
                shutil.rmtree(temp_dir_path, ignore_errors=True)
                # print(f"  ✓ Removed temp directory: {temp_dir_name}")
        except Exception as e:
            print(f"  Warning: Could not remove temp directory {temp_dir_name}: {e}")

def format_output_file(raw_results, output_file):
    """
    Format the raw BLAST+ results with a proper header.
    
    Args:
        raw_results: Path to raw results file
        output_file: Path for formatted output file
    """
    # print("Formatting output file...")
    import csv

    # Use canonical NCBI-style headers from Results and explicit raw->header mapping
    header_columns = Results.PRODUCER_HEADER
    header = "\t".join(header_columns) + "\n"

    with open(raw_results, 'r') as infile, open(output_file, 'w') as outfile:
        outfile.write(header)

        # Raw BLAST outfmt (configured in run_blast_search); use constants
        reader = csv.DictReader(infile, delimiter='\t', fieldnames=Results.RAW_OUTFMT_FIELDS)

        for row in reader:
            # For each header column, fetch the value from its corresponding raw field
            ordered_values = []
            for header_name in header_columns:
                raw_key = Results.HEADER_TO_RAW.get(header_name)
                ordered_values.append(row.get(raw_key, '') if raw_key else '')
            outfile.write("\t".join(ordered_values) + "\n")

    # print(f"✓ Results written to: {output_file}")

def cleanup_blast_temp_files(work_dir, pattern_prefixes=None):
    """
    Clean up any remaining BLAST+ temporary files that might have been missed.
    This function ONLY removes temporary blast_tmp_* files, NOT the target blast_db files.
    
    Args:
        work_dir: Working directory to clean
        pattern_prefixes: List of prefixes to look for (default: common BLAST+ temp prefixes)
    """
    if pattern_prefixes is None:
        pattern_prefixes = ['blast_tmp_']
    
    # print("Performing final cleanup of any remaining temporary files...")
    # print("Note: Only removing temporary blast_tmp_* files, preserving target BLAST database files")
    
    files_removed = 0
    for prefix in pattern_prefixes:
        pattern = f"{prefix}*"
        matches = list(Path(work_dir).glob(pattern))
        for match in matches:
            # Double-check: only remove files that start with our temporary prefixes
            if match.name.startswith(tuple(pattern_prefixes)):
                try:
                    if match.is_dir():
                        shutil.rmtree(match, ignore_errors=True)
                    else:
                        match.unlink()
                    # print(f"  ✓ Removed: {match.name}")
                    files_removed += 1
                except Exception as e:
                    print(f"  Warning: Could not remove {match.name}: {e}")
    
    if files_removed > 0:
        # print(f"✓ Final cleanup completed: {files_removed} files removed")
        pass
    else:
        # print("✓ No additional temporary files found")
        pass

def get_sequence_count(fasta_file):
    """Get the number of sequences in a FASTA file."""
    try:
        count = 0
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count
    except Exception:
        return "Unknown"

def main():
    """Main function to handle command line arguments and execute the BLAST+ nucleotide workflow."""
    parser = argparse.ArgumentParser(
        description="Search a nucleotide or protein query FASTA against a genome sequence using BLAST+",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python blast-genome-search.py query.fasta GCA_000507305.1 results.tsv
  python blast-genome-search.py dna_query.fasta GCA_001939145.1 search_results.tsv
  python blast-genome-search.py protein_query.faa GCA_000507305.1 results.tsv --query-type protein
  python blast-genome-search.py query.fasta GCA_000507305.1 results.tsv --evalue 1e-5
  python blast-genome-search.py protein.faa GCA_000507305.1 results.tsv --min-word-size 2
        """
    )
    
    parser.add_argument("query_fasta", help="Path to query FASTA file (DNA or protein)")
    parser.add_argument("accession", help="NCBI genome accession (e.g., GCA_000507305.1)")
    parser.add_argument("output_file", help="Path for output results file")
    parser.add_argument("--force-rebuild", action="store_true",
                       help="Force rebuilding of BLAST+ database")
    parser.add_argument("--evalue", type=float, default=1e-3,
                       help="E-value threshold for BLAST search (default: 1e-3)")
    parser.add_argument("--query-type", choices=['auto', 'dna', 'protein'], default='auto',
                       help="Query sequence type: 'auto' to detect, 'dna' for nucleotide, 'protein' for amino acid (default: auto)")
    parser.add_argument("--min-word-size", type=int, default=None,
                       help="Minimum word size for small matches. For DNA: minimum 4 nucleotides (default: 11). For protein: minimum 2 amino acids (default: 3). Use this flag to find very small matches (e.g., 2-3 amino acids = 6-9 bp)")

    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.query_fasta):
        print(f"Error: Query FASTA file '{args.query_fasta}' not found")
        sys.exit(1)
    
    if not args.accession:
        print("Error: Accession cannot be empty")
        sys.exit(1)
    
    # Detect or use specified query type
    if args.query_type == 'auto':
        # print("Detecting query sequence type...")
        query_type = detect_query_type(args.query_fasta)
        # print(f"✓ Detected query type: {query_type}")
    else:
        query_type = args.query_type
    
    is_protein = (query_type == 'protein')
    search_type = "tblastn" if is_protein else "blastn"
   
    """ 
    print("=== BLAST+ Genome Search ===")
    print(f"Query FASTA: {args.query_fasta}")
    print(f"Genome Accession: {args.accession}")
    print(f"Output file: {args.output_file}")
    print(f"Cache directory: {args.cache_dir}")
    print(f"Query type: {query_type}")
    print(f"Search type: {search_type}")
    print(f"E-value threshold: {args.evalue}")
    if args.min_word_size is not None:
        print(f"Minimum word size: {args.min_word_size} ({'amino acids' if is_protein else 'nucleotides'})")
    print()
    """ 
    
    # Step 1: Download genome FASTA using coral/ncbi/download.py
    # print("Step 1: Downloading genome FASTA...")
    try:
        cache_dir = download_and_extract_by_accession(args.accession, args.cache_dir)
        # print(f"✓ Genome data cached in: {cache_dir}")
    except Exception as e:
        print(f"Error downloading genome: {e}")
        sys.exit(1)
    
    # Find the genome FASTA file in the cache
    genome_fasta = None
    for pattern in ["*.fna", "*.fasta"]:
        matches = list(Path(cache_dir).glob(pattern))
        if matches:
            genome_fasta = str(matches[0])
            break
    
    if not genome_fasta:
        print("Error: Could not find genome FASTA file in cache")
        sys.exit(1)
    
    # print(f"✓ Genome FASTA found: {genome_fasta}")
    # print()
    
    # Step 2: Check if BLAST+ database already exists in cache
    db_path = "blast_db"  # Relative path in the mounted directory
    
    if not args.force_rebuild and check_blast_db_exists(os.path.join(cache_dir, db_path)):
        # print("Step 2: Using existing BLAST+ database from cache...")
        # print(f"✓ Database found: {db_path}")
        pass
    else:
        # print("Step 2: Creating BLAST+ database...")
        # Use the genome FASTA filename relative to the cache directory
        genome_fasta_rel = os.path.basename(genome_fasta)
        
        create_blast_database(genome_fasta_rel, db_path, cache_dir)
    
    # print()

    # Step 3: Perform BLAST+ search
    # print("Step 3: Running BLAST+ search...")
    # Copy query FASTA to the cache directory for Docker operations
    query_fasta_rel = "query.fasta"
    query_fasta_dest = os.path.join(cache_dir, query_fasta_rel)
    shutil.copy2(args.query_fasta, query_fasta_dest)
    
    raw_results = "blast_results.out"  # Relative path in the mounted directory
    raw_results_full = os.path.join(cache_dir, raw_results)
    
    try:
        run_blast_search(query_fasta_rel, db_path, raw_results, cache_dir, args.evalue, query_type, args.min_word_size)
        
        # print()
        
        # Step 4: Format output
        # print("Step 4: Formatting output...")
        format_output_file(raw_results_full, args.output_file)
        
    finally:
        # Clean up temporary files - this will run even if an exception occurs
        # print("Step 5: Cleaning up temporary files...")
        temp_files_to_remove = [
            query_fasta_dest,  # query.fasta
            raw_results_full,  # blast_results.out
        ]
        
        for temp_file in temp_files_to_remove:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                    # print(f"✓ Removed temporary file: {os.path.basename(temp_file)}")
            except Exception as e:
                print(f"Warning: Could not remove {temp_file}: {e}")
        
        # print("✓ Basic cleanup completed")
        
        # Perform final cleanup of any remaining BLAST+ temporary files
        # This only removes blast_tmp_* files, NOT the target blast_db files
        cleanup_blast_temp_files(cache_dir)
    
    # print()
    
    # Display summary
    """
    print("=== Search Summary ===")
    print(f"Query sequences: {get_sequence_count(args.query_fasta)}")
    print(f"Genome sequences: {get_sequence_count(genome_fasta)}")
    """
    
    # Count results (subtract header line)
    try:
        with open(args.output_file, 'r') as f:
            result_count = sum(1 for line in f) - 1  # Subtract header
        print(f"Search results: {result_count}")
    except Exception:
        print("Search results: Unknown")
    
    """
    print()
    print("Search completed successfully!")
    """
    print(f"Blast results saved to: {args.output_file}")
    
    # Final cleanup check - this ensures we catch any temporary files that might have been missed
    # This only removes blast_tmp_* files, NOT the target blast_db files
    # print()
    cleanup_blast_temp_files(cache_dir)

if __name__ == "__main__":
    main()
