import os
from pathlib import Path


class DefaultPath(object):

    @staticmethod
    def ncbi_download_dir():
        return "ncbi-downloads"

    @staticmethod
    def ncbi_genome_dir(genome_accession):
        return os.path.join(DefaultPath.ncbi_download_dir(), "ncbi_dataset/data", genome_accession)

    @staticmethod
    def ncbi_genome_faa(genome_accession):
        return os.path.join(DefaultPath.ncbi_genome_dir(genome_accession), "protein.faa")

    @staticmethod
    def ncbi_genome_gff(genome_accession):
        return os.path.join(DefaultPath.ncbi_genome_dir(genome_accession), "genomic.gff")

    @staticmethod
    def ncbi_genome_fna(genome_accession):
        ncbi_download_dir = os.path.join(DefaultPath.ncbi_genome_dir(genome_accession))
        fna_file = None
        for pattern in ["*.fna", "*.fasta"]:
            matches = list(Path(ncbi_download_dir).glob(pattern))
            if matches:
                fna_file = str(matches[0])
        return fna_file

    @staticmethod
    def kegg_hmm_dir():
        return "kegg-downloads/profiles"

    @staticmethod
    def kegg_hmm(hmm_model):
        return os.path.join(DefaultPath.kegg_hmm_dir(), hmm_model+".hmm")
