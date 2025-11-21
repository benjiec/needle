import os
import tempfile
import unittest

from needle.blast import Results, group_matches, ProteinMatch


class TestBlastResults(unittest.TestCase):
    def test_parse_ncbi_header_and_extract_sequences(self):
        # Create synthetic FASTA and TSV with NCBI-style headers
        with tempfile.TemporaryDirectory() as tmpdir:
            query_fasta_path = os.path.join(tmpdir, "query.faa")
            target_fasta_path = os.path.join(tmpdir, "target.faa")
            results_tsv_path = os.path.join(tmpdir, "results.tsv")

            # >Q1
            # MNOPQRST
            with open(query_fasta_path, "w") as f:
                f.write(">Q1 synthetic query\n")
                f.write("MNOPQRST\n")

            # >S1
            # ABCDEFGHIJKLMNOPQRST
            with open(target_fasta_path, "w") as f:
                f.write(">S1 synthetic subject\n")
                f.write("ABCDEFGHIJKLMNOPQRST\n")

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
            # Target positions min(10,7)=7..10 from ABCDEFGHIJ... -> "GHIJ"
            self.assertEqual(m.target_sequence, "GHIJ")
            # Matched sequence preserved
            self.assertEqual(m.matched_sequence, "ABCD")

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

    def test_group_matches_and_proteinmatch_booleans(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            query_fasta_path = os.path.join(tmpdir, "query.faa")
            results_tsv_path = os.path.join(tmpdir, "results.tsv")

            # Query length = 30 (1..30)
            with open(query_fasta_path, "w") as f:
                f.write(">QX\n")
                f.write("A" * 30 + "\n")

            # Build TSV rows (NCBI headers)
            header = "\t".join(Results.PRODUCER_HEADER)
            rows = []
            # Complete coverage with small gaps in the middle (<=10): covers_start_to_end True, likely_complete True
            # On target Sx, intervals within 10_000 distance
            rows.append(["QX", "Sx", "1e-20", "90.0", "1", "5", "1000", "1020", "AAAAA"])
            rows.append(["QX", "Sx", "1e-20", "90.0", "10", "20", "1100", "1120", "AAAAA"])
            rows.append(["QX", "Sx", "1e-20", "90.0", "25", "30", "1200", "1220", "AAAAA"])

            # Only 5' end group on same target Sy
            rows.append(["QX", "Sy", "1e-10", "85.0", "1", "10", "5000", "5020", "BBBBB"])

            # Only 3' end group on same target Sz
            rows.append(["QX", "Sz", "1e-10", "85.0", "21", "30", "7000", "7020", "CCCCC"])

            # Overlapping results on target So (query overlap True)
            rows.append(["QX", "So", "1e-15", "95.0", "5", "15", "8000", "8020", "DDDDD"])
            rows.append(["QX", "So", "1e-15", "95.0", "12", "22", "8030", "8050", "EEEEE"])

            # Big intron split on target Si: two clusters separated by > max_intron_length (use 10_000 default)
            rows.append(["QX", "Si", "1e-12", "92.0", "1", "5", "100000", "100020", "FFFFF"])
            rows.append(["QX", "Si", "1e-12", "92.0", "6", "10", "120500", "120520", "GGGGG"])  # distance ~20500

            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")

            res = Results(results_tsv_path, query_fasta_path=query_fasta_path)
            pms = group_matches(res)

            # Index by target_id for assertions (note Si should produce two ProteinMatch instances)
            by_target = {}
            for pm in pms:
                by_target.setdefault(pm.target_id, []).append(pm)

            # Sx: complete with small gaps
            self.assertIn("Sx", by_target)
            self.assertEqual(len(by_target["Sx"]), 1)
            pm_sx = by_target["Sx"][0]
            self.assertTrue(pm_sx.covers_start_to_end)
            self.assertTrue(pm_sx.likely_complete)
            self.assertFalse(pm_sx.query_overlap)
            self.assertEqual(pm_sx.query_start, 1)
            self.assertEqual(pm_sx.query_end, 30)
            self.assertEqual(pm_sx.target_start, 1000)
            self.assertEqual(pm_sx.target_end, 1220)

            # Sy: only 5' end
            self.assertIn("Sy", by_target)
            self.assertEqual(len(by_target["Sy"]), 1)
            pm_sy = by_target["Sy"][0]
            self.assertFalse(pm_sy.covers_start_to_end)
            self.assertFalse(pm_sy.likely_complete)
            self.assertFalse(pm_sy.query_overlap)
            self.assertEqual(pm_sy.query_start, 1)
            self.assertEqual(pm_sy.query_end, 10)

            # Sz: only 3' end
            self.assertIn("Sz", by_target)
            self.assertEqual(len(by_target["Sz"]), 1)
            pm_sz = by_target["Sz"][0]
            self.assertFalse(pm_sz.covers_start_to_end)
            self.assertFalse(pm_sz.likely_complete)
            self.assertFalse(pm_sz.query_overlap)
            self.assertEqual(pm_sz.query_start, 21)
            self.assertEqual(pm_sz.query_end, 30)

            # So: overlapping query intervals
            self.assertIn("So", by_target)
            self.assertEqual(len(by_target["So"]), 1)
            pm_so = by_target["So"][0]
            self.assertTrue(pm_so.query_overlap)
            # Covers start-to-end? No
            self.assertFalse(pm_so.covers_start_to_end)
            self.assertFalse(pm_so.likely_complete)

            # Si: two clusters due to large intron
            self.assertIn("Si", by_target)
            self.assertEqual(len(by_target["Si"]), 2)
            # Ensure clusters' query ranges correspond
            ranges = sorted((pm.query_start, pm.query_end) for pm in by_target["Si"])
            self.assertEqual(ranges, [(1, 5), (6, 10)])

    def test_same_target_far_apart_split_by_max_intron_length(self):
        # Explicit test: same target ID, distance > max_intron_length creates separate groups
        with tempfile.TemporaryDirectory() as tmpdir:
            query_fasta_path = os.path.join(tmpdir, "q.faa")
            results_tsv_path = os.path.join(tmpdir, "r.tsv")
            with open(query_fasta_path, "w") as f:
                f.write(">Q\n")
                f.write("A" * 40 + "\n")
            header = "\t".join(Results.PRODUCER_HEADER)
            rows = [
                ["Q", "T", "1e-5", "90", "1", "10", "1000", "1020", "AAAAAAAAAA"],
                ["Q", "T", "1e-5", "90", "12", "20", "30000", "30020", "AAAAAAAAAA"],  # far apart
            ]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")
            res = Results(results_tsv_path, query_fasta_path=query_fasta_path)
            pms = group_matches(res, max_intron_length=10000)
            by_target = {}
            for pm in pms:
                by_target.setdefault(pm.target_id, []).append(pm)
            self.assertIn("T", by_target)
            self.assertEqual(len(by_target["T"]), 2)  # split due to large intron distance

    def test_same_target_close_overlapping_grouped_with_overlap_true(self):
        # Same target, overlapping query regions, close on target -> grouped with overlap=True
        with tempfile.TemporaryDirectory() as tmpdir:
            query_fasta_path = os.path.join(tmpdir, "q.faa")
            results_tsv_path = os.path.join(tmpdir, "r.tsv")
            with open(query_fasta_path, "w") as f:
                f.write(">Q\n")
                f.write("A" * 50 + "\n")
            header = "\t".join(Results.PRODUCER_HEADER)
            rows = [
                ["Q", "T2", "1e-12", "95", "5", "20", "1000", "1020", "X" * 10],
                ["Q", "T2", "1e-12", "95", "15", "30", "1030", "1050", "X" * 10],  # overlap in query; close on target
            ]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")
            res = Results(results_tsv_path, query_fasta_path=query_fasta_path)
            pms = group_matches(res, max_intron_length=10000)
            self.assertEqual(len(pms), 1)
            pm = pms[0]
            self.assertEqual(pm.target_id, "T2")
            self.assertTrue(pm.query_overlap)

    def test_transitive_clustering_chain_within_threshold(self):
        # Three matches: first and third farther than threshold apart,
        # but each consecutive pair within threshold -> grouped together via transitivity
        with tempfile.TemporaryDirectory() as tmpdir:
            query_fasta_path = os.path.join(tmpdir, "q.faa")
            results_tsv_path = os.path.join(tmpdir, "r.tsv")
            with open(query_fasta_path, "w") as f:
                f.write(">Q\n")
                f.write("A" * 60 + "\n")
            header = "\t".join(Results.PRODUCER_HEADER)
            rows = [
                ["Q", "T3", "1e-10", "90", "1", "10", "1000", "1020", "Z" * 10],    # end at 1020
                ["Q", "T3", "1e-10", "90", "11", "20", "1030", "1050", "Z" * 10],  # start 1030 (gap 10)
                ["Q", "T3", "1e-10", "90", "21", "30", "1065", "1085", "Z" * 10],  # start 1065 (gap 15 from previous)
            ]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")
            res = Results(results_tsv_path, query_fasta_path=query_fasta_path)
            # Threshold 20 -> pairs (1-2) gap 10, (2-3) gap 15 -> all grouped
            pms = group_matches(res, max_intron_length=20)
            self.assertEqual(len(pms), 1)
            pm = pms[0]
            self.assertEqual(pm.target_id, "T3")
            self.assertEqual(pm.query_start, 1)
            self.assertEqual(pm.query_end, 30)


if __name__ == "__main__":
    unittest.main()


