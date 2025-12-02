#!/usr/bin/env python3
"""
NCBI Datasets v2alpha Download Script (Python Version)
Downloads DNA and protein FASTA and GFF files for a list of genome accessions

This script provides the same functionality as ncbi-download.sh but uses the
coral/ncbi/download.py module for better code reuse and cross-platform support.
"""

import argparse
import os
import sys
import time
from pathlib import Path

from ncbi import download_and_extract_by_accession, FILE_TYPE_PATTERNS

def print_status(message: str, level: str = "INFO"):
    """Print a timestamped status message with color coding."""
    from datetime import datetime
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    # Color codes for different levels
    colors = {
        "INFO": "\033[0;34m",    # Blue
        "SUCCESS": "\033[0;32m", # Green
        "WARNING": "\033[1;33m", # Yellow
        "ERROR": "\033[0;31m",   # Red
        "NC": "\033[0m"          # No Color
    }
    
    color = colors.get(level, colors["INFO"])
    print(f"{color}[{timestamp}] {message}{colors['NC']}")

def show_usage():
    """Display usage information."""
    print("""Usage: python ncbi-download.py [OPTIONS] <accession_file|accession_id>

Downloads protein FASTA and GFF files from NCBI Datasets API v2alpha

Arguments:
  accession_file    File containing one genome accession per line (e.g., GCF_000001405.40)
  accession_id      A single genome accession ID (e.g., GCF_000001405.40)

Options:
  -o, --output DIR     Output directory (default: ncbi-downloads)
  -d, --delay SECONDS  Delay between requests (default: 1)
  -k, --keep-zip       Keep original zip files after extraction
  -h, --help           Show this help message

Valid annotation types (automatically included):
  GENOME_FASTA, PROT_FASTA, GENOME_GFF

Examples:
  python ncbi-download.py -o my_genomes accessions.txt
  python ncbi-download.py -o my_genomes GCF_000001405.40
""")

def process_accession(accession: str, output_dir: str, delay: float = 1.0):
    """
    Process a single accession by downloading and extracting all file types.
    
    Args:
        accession: The NCBI genome accession
        output_dir: The output directory for downloads
        delay: Delay between requests in seconds
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        print_status(f"Processing: {accession}", "INFO")
        
        # Use the coral/ncbi/download.py module
        cache_dir = download_and_extract_by_accession(accession, output_dir)
        
        if cache_dir and os.path.exists(cache_dir):
            print_status(f"Successfully processed: {accession}", "SUCCESS")
            return True
        else:
            print_status(f"Failed to process: {accession}", "ERROR")
            return False
            
    except Exception as e:
        print_status(f"Error processing {accession}: {str(e)}", "ERROR")
        return False
    finally:
        # Apply delay if specified
        if delay > 0:
            time.sleep(delay)

def process_accession_file(accession_file: str, output_dir: str, delay: float = 1.0):
    """
    Process multiple accessions from a file.
    
    Args:
        accession_file: Path to file containing accessions
        output_dir: The output directory for downloads
        delay: Delay between requests in seconds
        
    Returns:
        tuple: (success_count, failed_count, total_count)
    """
    success_count = 0
    failed_count = 0
    total_count = 0
    
    try:
        with open(accession_file, 'r') as f:
            accessions = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        
        total_count = len(accessions)
        print_status(f"Starting download of {total_count} genome(s)", "INFO")
        print_status(f"Output directory: {output_dir}", "INFO")
        print_status(f"Annotation types: {', '.join(FILE_TYPE_PATTERNS.keys())}", "INFO")
        print_status(f"Delay between requests: {delay}s", "INFO")
        
        for i, accession in enumerate(accessions, 1):
            print_status(f"[{i}/{total_count}] Processing: {accession}", "INFO")
            
            if process_accession(accession, output_dir, delay):
                success_count += 1
            else:
                failed_count += 1
                
    except FileNotFoundError:
        print_status(f"Error: Accession file '{accession_file}' not found", "ERROR")
        return 0, 0, 0
    except Exception as e:
        print_status(f"Error reading accession file: {str(e)}", "ERROR")
        return 0, 0, 0
    
    return success_count, failed_count, total_count

def process_single_accession(accession: str, output_dir: str, delay: float = 1.0):
    """
    Process a single accession.
    
    Args:
        accession: The NCBI genome accession
        output_dir: The output directory for downloads
        delay: Delay between requests in seconds
        
    Returns:
        tuple: (success_count, failed_count, total_count)
    """
    print_status(f"Starting download of 1 genome", "INFO")
    print_status(f"Output directory: {output_dir}", "INFO")
    print_status(f"Annotation types: {', '.join(FILE_TYPE_PATTERNS.keys())}", "INFO")
    print_status(f"Delay between requests: {delay}s", "INFO")
    
    success_count = 0
    failed_count = 0
    
    if process_accession(accession, output_dir, delay):
        success_count = 1
    else:
        failed_count = 1
    
    return success_count, failed_count, 1

def main():
    """Main function to handle command line arguments and execute downloads."""
    parser = argparse.ArgumentParser(
        description="Download NCBI genome datasets using the coral/ncbi/download.py module",
        add_help=False  # We'll handle help manually to match the bash script
    )
    
    parser.add_argument("-o", "--output", default="ncbi-downloads",
                       help="Output directory (default: ncbi-downloads)")
    parser.add_argument("-d", "--delay", type=float, default=1.0,
                       help="Delay between requests in seconds (default: 1)")
    parser.add_argument("-k", "--keep-zip", action="store_true",
                       help="Keep original zip files after extraction")
    parser.add_argument("-h", "--help", action="store_true",
                       help="Show this help message")
    
    # Parse known args first to check for help
    args, remaining = parser.parse_known_args()
    
    if args.help or len(remaining) == 0:
        show_usage()
        sys.exit(0)
    
    # Parse the remaining arguments
    if len(remaining) != 1:
        print_status("Error: Exactly one argument (accession file or ID) is required", "ERROR")
        show_usage()
        sys.exit(1)
    
    accession_arg = remaining[0]
    output_dir = args.output
    delay = args.delay
    
    # Create output directory
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        print_status(f"Error: Could not create output directory '{output_dir}': {e}", "ERROR")
        sys.exit(1)
    
    # Process accessions
    if os.path.isfile(accession_arg):
        # Process as accession file
        success_count, failed_count, total_count = process_accession_file(
            accession_arg, output_dir, delay
        )
    else:
        # Process as single accession
        success_count, failed_count, total_count = process_single_accession(
            accession_arg, output_dir, delay
        )
    
    # Print summary
    print()
    print_status("=== DOWNLOAD SUMMARY ===", "INFO")
    print_status(f"Successful downloads: {success_count}", "SUCCESS")
    print_status(f"Failed downloads: {failed_count}", "ERROR")
    print_status(f"Total processed: {total_count}", "INFO")
    
    # List what was downloaded
    if success_count > 0:
        print()
        print_status(f"Files extracted to: {output_dir}/ncbi_dataset/data", "INFO")
        print_status("Each genome is in its own subdirectory:", "INFO")
        
        data_dir = os.path.join(output_dir, "ncbi_dataset", "data")
        if os.path.exists(data_dir):
            for item in os.listdir(data_dir):
                item_path = os.path.join(data_dir, item)
                if os.path.isdir(item_path) and not item.startswith('.'):
                    print(f"  {item}")
    
    # Exit with appropriate code
    if failed_count == 0:
        print_status("All downloads completed successfully!", "SUCCESS")
        sys.exit(0)
    else:
        print_status("Some downloads failed. Check the log above for details.", "WARNING")
        sys.exit(1)

if __name__ == "__main__":
    main()
