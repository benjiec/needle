import os
import tempfile
import unittest
from needle.blast import Results


class TestParseBlastResults(unittest.TestCase):

    def test_parse_ncbi_header_and_extract_sequences_reverse(self):
        """Verify reverse-direction hit: target sequence normalized to 5'->3' and reverse-complemented for translation."""
        # Create synthetic FASTA and TSV with NCBI-style headers; reverse-direction on target (sstart > send)
        with tempfile.TemporaryDirectory() as tmpdir:
            query_fasta_path = os.path.join(tmpdir, "query.faa")
            target_fasta_path = os.path.join(tmpdir, "target.faa")
            results_tsv_path = os.path.join(tmpdir, "results.tsv")

            # >Q1
            # MNOPQRST
            with open(query_fasta_path, "w") as f:
                f.write(">Q1 synthetic query\n")
                f.write("MNOPQRST\n")

            # >S1 DNA: AAACCCGGGTTT
            with open(target_fasta_path, "w") as f:
                f.write(">S1 synthetic subject\n")
                f.write("AAACCCGGGTTT\n")

            header = "\t".join(Results.PRODUCER_HEADER)
            # qseqid sseqid evalue pident qstart qend sstart send sseq
            row = ["Q1", "S1", "1e-5", "99.9", "2", "5", "10", "7", "ABCD"]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(row) + "\n")

            res = Results(results_tsv_path, query_fasta_path, target_fasta_path)
            matches = res.matches()
            self.assertEqual(len(matches), 1)
            m = matches[0]
            self.assertEqual(m.query_accession, "Q1")
            self.assertEqual(m.target_accession, "S1")
            self.assertAlmostEqual(m.e_value, 1e-5, places=12)
            self.assertAlmostEqual(m.identity, 99.9, places=3)
            self.assertEqual(m.query_start, 2)
            self.assertEqual(m.query_end, 5)
            self.assertEqual(m.target_start, 10)
            self.assertEqual(m.target_end, 7)

            # Query positions 2..5 from MNOPQRST -> "NOPQ"
            self.assertEqual(m.query_sequence, "NOPQ")
            # Target 7..10 on AAACCCGGGTTT -> "GGGT"; reverse-complement -> "ACCC"
            self.assertEqual(m.target_sequence, "ACCC")
            # Ensure translated target (revcomp path) is a valid amino-acid sequence length <= query segment
            _ = m.target_sequence_translated()
            # Matched sequence preserved
            self.assertEqual(m.matched_sequence, "ABCD")

    def test_parse_ncbi_header_and_extract_sequences_forward(self):
        """Verify forward-direction hit: target sequence 5'->3' matches the query AA after translation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            query_fasta_path = os.path.join(tmpdir, "q.faa")
            target_fasta_path = os.path.join(tmpdir, "t.fna")
            results_tsv_path = os.path.join(tmpdir, "r.tsv")

            # Query protein containing MEF... to match positions 1..3
            with open(query_fasta_path, "w") as f:
                f.write(">Qfwd\n")
                f.write("MEFGHIKLMN\n")

            # Target DNA forward: ATG GAA TTT -> MEF at positions 5..13
            with open(target_fasta_path, "w") as f:
                f.write(">Tfwd\n")
                f.write("NNNNATGGAATTTNNNN\n")  # 1..4 N; 5..13 coding; 14..17 N

            header = "\t".join(Results.PRODUCER_HEADER)
            # Forward direction: sstart < send
            row = ["Qfwd", "Tfwd", "0", "100.0", "1", "3", "5", "13", "MEF"]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(row) + "\n")

            res = Results(results_tsv_path, query_fasta_path, target_fasta_path)
            ms = res.matches()
            self.assertEqual(len(ms), 1)
            m = ms[0]
            # Query substring 1..3 -> "MEF"
            self.assertEqual(m.query_sequence, "MEF")
            # Target DNA 5..13 5'->3'
            self.assertEqual(m.target_sequence, "ATGGAATTT")
            # Translated target equals matched_sequence equals query segment
            self.assertEqual(m.target_sequence_translated(), "MEF")
            self.assertEqual(m.matched_sequence, "MEF")

    def test_missing_query_fasta_only_target_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # Only target fasta provided
            target_fasta_path = os.path.join(tmpdir, "target.fna")
            results_tsv_path = os.path.join(tmpdir, "results.tsv")

            with open(target_fasta_path, "w") as f:
                f.write(">S2\n")
                f.write("AACCGGTTAACC\n")

            header = "\t".join(Results.PRODUCER_HEADER)
            # qseqid sseqid evalue pident qstart qend sstart send sseq
            row = ["Q2", "S2", "2e-6", "90.0", "3", "6", "5", "9", "CCGGT"]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(row) + "\n")

            res = Results(results_tsv_path, target_fasta_path=target_fasta_path)
            matches = res.matches()
            self.assertEqual(len(matches), 1)
            m = matches[0]
            self.assertIsNone(m.query_sequence)
            self.assertEqual(m.target_sequence, "GGTTA")  # positions 5..9 on AACCGGTTAACC

    def test_missing_target_fasta_only_query_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            query_fasta_path = os.path.join(tmpdir, "query.faa")
            results_tsv_path = os.path.join(tmpdir, "results.tsv")

            with open(query_fasta_path, "w") as f:
                f.write(">Q3\n")
                f.write("MNOPQRSTUVW\n")

            header = "\t".join(Results.PRODUCER_HEADER)
            row = ["Q3", "S3", "5e-4", "88.5", "4", "8", "2", "6", "ABCDE"]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(row) + "\n")

            res = Results(results_tsv_path, query_fasta_path=query_fasta_path)
            matches = res.matches()
            self.assertEqual(len(matches), 1)
            m = matches[0]
            self.assertEqual(m.query_sequence, "PQRST")  # positions 4..8 on MNOPQRSTUVW
            self.assertIsNone(m.target_sequence)


if __name__ == "__main__":
    unittest.main()
