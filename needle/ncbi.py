import os
import subprocess
import zipfile
import glob
import shutil
import requests
import xml.etree.ElementTree as ET
import csv


BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"

# Mapping of file types to their expected filename patterns
# This makes it easy to add new file types by just adding to this mapping
FILE_TYPE_PATTERNS = {
    "GENOME_FASTA": "*.fna",
    "PROT_FASTA": "protein.faa",
    "GENOME_GFF": "genomic.gff"
}

def download_and_extract_by_accession(accession: str, cache_dir: str) -> str:
    """
    Downloads and extracts all file types for a given accession to a cache directory.
    
    This function implements the same directory structure as scripts/ncbi-download.sh:
    cache_dir/
    └── ncbi_dataset/
        └── data/
            ├── dataset_catalog.json
            ├── assembly_data_report.jsonl
            └── {accession}/
                ├── {accession}*.fna (GENOME_FASTA)
                ├── genomic.gff (GENOME_GFF)
                └── protein.faa (PROT_FASTA)
    
    Args:
        accession: The NCBI genome accession (e.g., 'GCA_000507305.1')
        cache_dir: The cache directory that follows the ncbi-downloads structure
        
    Returns:
        Path to the extracted accession directory in the cache
        
    Raises:
        Exception: If download or extraction fails
    """
    # Create cache directory structure
    cache_data_dir = os.path.join(cache_dir, "ncbi_dataset", "data")
    accession_cache_dir = os.path.join(cache_data_dir, accession)
    
    # Check if files already exist in cache
    if os.path.exists(accession_cache_dir):
        # Only check if GENOME_FASTA pattern exists - this is the main file we need
        genome_pattern = FILE_TYPE_PATTERNS["GENOME_FASTA"]
        genome_files = glob.glob(os.path.join(accession_cache_dir, genome_pattern))
        if genome_files:
            return accession_cache_dir
    
    # Download all file types by concatenating the keys from the mapping
    annotation_types = ",".join(FILE_TYPE_PATTERNS.keys())
    url = f"{BASE_URL}/{accession}/download?include_annotation_type={annotation_types}"
    
    # Create temporary directories
    os.makedirs(cache_dir, exist_ok=True)
    temp_zip = os.path.join(cache_dir, f"{accession}.zip")
    temp_extract_dir = os.path.join(cache_dir, accession)
    
    try:
        # Download the zip file
        print(f"Downloading {accession} from NCBI...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(temp_zip, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        # Extract the zip file
        print(f"Extracting {accession}...")
        with zipfile.ZipFile(temp_zip, 'r') as zip_ref:
            zip_ref.extractall(temp_extract_dir)
        
        # Create the cache directory structure
        os.makedirs(cache_data_dir, exist_ok=True)
        
        # Move the accession directory to the flattened location
        if os.path.exists(accession_cache_dir):
            shutil.rmtree(accession_cache_dir)
        
        # Find the nested accession directory
        nested_accession_dir = os.path.join(temp_extract_dir, "ncbi_dataset", "data", accession)
        if os.path.exists(nested_accession_dir):
            shutil.move(nested_accession_dir, accession_cache_dir)
        else:
            raise FileNotFoundError(f"Expected directory structure not found for {accession}")
        
        # Copy metadata files if they exist
        metadata_files = ["assembly_data_report.jsonl", "dataset_catalog.json"]
        for metadata_file in metadata_files:
            src_path = os.path.join(temp_extract_dir, "ncbi_dataset", "data", metadata_file)
            dst_path = os.path.join(cache_data_dir, metadata_file)
            if os.path.exists(src_path):
                shutil.copy2(src_path, dst_path)
        
        print(f"Successfully cached {accession}")
        return accession_cache_dir
        
    except Exception as e:
        # Clean up on failure
        if os.path.exists(temp_zip):
            os.remove(temp_zip)
        if os.path.exists(temp_extract_dir):
            shutil.rmtree(temp_extract_dir, ignore_errors=True)
        if os.path.exists(accession_cache_dir):
            shutil.rmtree(accession_cache_dir, ignore_errors=True)
        raise e
    
    finally:
        # Clean up temporary files
        if os.path.exists(temp_zip):
            os.remove(temp_zip)
        if os.path.exists(temp_extract_dir):
            shutil.rmtree(temp_extract_dir, ignore_errors=True)

def download_and_extract_fasta(accession: str, output_dir: str, file_type: str, cache_dir: str) -> str:
    """
    Downloads and extracts FASTA file for a given accession and file type.
    Uses the shared cache to avoid re-downloading files.
    
    Args:
        accession: The NCBI genome accession
        output_dir: Directory to copy the final FASTA file to
        file_type: Type of FASTA file (e.g., 'GENOME_FASTA', 'PROT_FASTA')
        cache_dir: The cache directory that follows the ncbi-downloads structure
        
    Returns:
        Path to the extracted FASTA file in the output directory
    """
    # First, ensure the file is in the cache
    accession_cache_dir = download_and_extract_by_accession(accession, cache_dir)
    
    # Find the appropriate file using the mapping
    if file_type not in FILE_TYPE_PATTERNS:
        raise ValueError(f"Unknown file type: {file_type}. Valid types: {list(FILE_TYPE_PATTERNS.keys())}")
    
    pattern = FILE_TYPE_PATTERNS[file_type]
    files = glob.glob(os.path.join(accession_cache_dir, pattern))
    
    if not files:
        raise FileNotFoundError(f"No {file_type} file found for {accession} using pattern '{pattern}'")
    
    # Copy the file to the output directory
    os.makedirs(output_dir, exist_ok=True)
    source_file = files[0]
    
    # Determine the appropriate file extension based on the file type
    if file_type == "GENOME_FASTA":
        file_ext = ".fna"
    elif file_type == "PROT_FASTA":
        file_ext = ".faa"
    else:
        # For other file types, extract extension from the source file
        file_ext = os.path.splitext(source_file)[1]
    
    final_file_path = os.path.join(output_dir, f"{accession}{file_ext}")
    shutil.copy2(source_file, final_file_path)
    
    return final_file_path

def download_and_extract_prot_fasta(accession: str, output_dir: str, cache_dir: str) -> str:
    """
    Downloads and extracts protein FASTA file for a given accession.
    Uses the shared cache to avoid re-downloading files.
    
    Args:
        accession: The NCBI genome accession
        output_dir: Directory to copy the final protein FASTA file to
        cache_dir: The cache directory that follows the ncbi-downloads structure
        
    Returns:
        Path to the extracted protein FASTA file in the output directory
    """
    return download_and_extract_fasta(accession, output_dir, 'PROT_FASTA', cache_dir) 


def _parse_taxonomy_xml(root):
    """
    Parse taxonomy XML root and return (order, class_, family, genus, species) as strings.
    """
    order = class_ = family = genus = species = ''
    for taxon in root.findall('.//LineageEx/Taxon'):
        rank = taxon.findtext('Rank')
        name = taxon.findtext('ScientificName')
        if rank == 'order':
            order = name
        elif rank == 'class':
            class_ = name
        elif rank == 'family':
            family = name
        elif rank == 'genus':
            genus = name
        elif rank == 'species':
            species = name
    # If the current node is species, use its name
    for node in root.findall('.//Taxon'):
        rank = node.findtext('Rank')
        name = node.findtext('ScientificName')
        if rank == 'species':
            species = name
    return order, class_, family, genus, species

def parse_taxonomy_tsv(taxonomy_tsv):
    """
    Yield each row of a taxonomy TSV/CSV as a dict, robust to column order and delimiter.
    """
    with open(taxonomy_tsv, newline='') as f:
        sample = f.read(1024)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample)
        except csv.Error:
            dialect = csv.excel_tab  # fallback to tab
        reader = csv.DictReader(f, delimiter=dialect.delimiter)
        for row in reader:
            yield row

def fetch_and_append_taxonomy(genome_acc: str, taxonomy_tsv: str):
    """
    Fetch genome and taxonomy info and append to taxonomy_tsv.
    """

    fieldnames = ['Genome Accession', 'TaxID', 'Order', 'Class', 'Family', 'Genus', 'Species', 'Genome Name', 'Organism']
    # If file does not exist or is empty, create it with header
    if not os.path.exists(taxonomy_tsv) or os.path.getsize(taxonomy_tsv) == 0:
        with open(taxonomy_tsv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

    # Read all rows, filter out any with matching accession or taxid
    rows = []
    for row in parse_taxonomy_tsv(taxonomy_tsv):
        if row.get('Genome Accession', '') != genome_acc:
            rows.append(row)

    # Lookup UID for accession
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={genome_acc}&retmode=json"
    r = requests.get(url)
    r.raise_for_status()
    uid_list = r.json()['esearchresult']['idlist']
    if not uid_list:
        print(f"Warning: No UID found for accession {genome_acc}. Skipping taxonomy.", flush=True)
        return
    uid = uid_list[0]
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={uid}&retmode=json"
    r = requests.get(url)
    r.raise_for_status()
    doc = r.json()['result'][uid]
    species = doc.get('speciesname', '')
    taxid = doc.get('taxid', '')
    genome_name = doc.get('assemblyname', '')
    organism = doc.get('organism', '')

    # Fetch taxonomy info from NCBI taxonomy using taxid
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={taxid}&retmode=xml"
    r = requests.get(url)
    r.raise_for_status()
    root = ET.fromstring(r.text)
    order, class_, family, genus, species_from_tax = _parse_taxonomy_xml(root)
    # Prefer species from taxonomy, fallback to assembly summary if not present
    final_species = species_from_tax or species
    new_row = {
        'Order': order,
        'Class': class_,
        'Family': family,
        'Genus': genus,
        'Species': final_species,
        'Genome Accession': genome_acc,
        'Genome Name': genome_name,
        'Organism': organism,
        'TaxID': taxid
    }
    rows.append(new_row)
    # Write all rows back using DictWriter
    with open(taxonomy_tsv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row) 
