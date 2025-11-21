import os
import subprocess
import zipfile
import glob
import shutil
import requests

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