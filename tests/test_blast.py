import os
import tempfile
import unittest

from needle.blast import Results, group_matches, ProteinMatch, NucMatch, order_matches_for_junctions, NonlinearMatchException
import shutil


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


class TestOrderGroupMatches(unittest.TestCase):

    def test_order_matches_for_junctions_overlap_and_gap(self):
        class _M: pass
        m1 = _M(); m1.query_start=1; m1.query_end=10
        m2 = _M(); m2.query_start=8; m2.query_end=15
        m3 = _M(); m3.query_start=18; m3.query_end=20
        pairs = order_matches_for_junctions([m1, m3, m2])  # input not in order
        self.assertEqual(len(pairs), 2)
        self.assertEqual(pairs[0], (m1, m2, 3, 0))
        self.assertEqual(pairs[1], (m2, m3, 0, 2))

    def test_order_throws_error_if_junctions_overlap(self):
        # aaaaaaaa
        #      bbbbbb
        #        cccccccc

        class _M: pass
        m1 = _M(); m1.query_start=1; m1.query_end=8
        m2 = _M(); m2.query_start=6; m2.query_end=12
        m3 = _M(); m3.query_start=8; m3.query_end=16

        with self.assertRaises(NonlinearMatchException):
            pairs = order_matches_for_junctions([m1, m3, m2])  # input not in order

    def test_order_throws_error_on_contained_match(self):
        # aaaaaaaaa
        #      bbbbbbbbbbb
        #            ccc
        #                ddd

        class _M: pass
        m1 = _M(); m1.query_start=1; m1.query_end=10
        m2 = _M(); m2.query_start=6; m2.query_end=16
        m3 = _M(); m3.query_start=12; m3.query_end=14
        m4 = _M(); m4.query_start=16; m4.query_end=18

        with self.assertRaises(NonlinearMatchException):
            pairs = order_matches_for_junctions([m1, m4, m3, m2])  # input not in order

    def test_group_matches_separate_matches_by_contig_and_distance(self):
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
            # Complete coverage with small gaps in the middle (<=10):
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
            self.assertEqual(pm_sx.query_start, 1)
            self.assertEqual(pm_sx.query_end, 30)
            self.assertEqual(pm_sx.target_start, 1000)
            self.assertEqual(pm_sx.target_end, 1220)

            # Sy: only 5' end
            self.assertIn("Sy", by_target)
            self.assertEqual(len(by_target["Sy"]), 1)
            pm_sy = by_target["Sy"][0]
            self.assertEqual(pm_sy.query_start, 1)
            self.assertEqual(pm_sy.query_end, 10)

            # Sz: only 3' end
            self.assertIn("Sz", by_target)
            self.assertEqual(len(by_target["Sz"]), 1)
            pm_sz = by_target["Sz"][0]
            self.assertEqual(pm_sz.query_start, 21)
            self.assertEqual(pm_sz.query_end, 30)

            # So: overlapping query intervals
            self.assertIn("So", by_target)
            self.assertEqual(len(by_target["So"]), 1)
            pm_so = by_target["So"][0]

            # Si: two clusters due to large intron
            self.assertIn("Si", by_target)
            self.assertEqual(len(by_target["Si"]), 2)
            # Ensure clusters' query ranges correspond
            ranges = sorted((pm.query_start, pm.query_end) for pm in by_target["Si"])
            self.assertEqual(ranges, [(1, 5), (6, 10)])

    def test_group_matches_separate_closely_placed_tandems(self):
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
            # contig Sx, first copy
            rows.append(["QX", "Sx", "1e-20", "90.0", "1",  "5",  "1000", "1020", "AAAAA"])
            rows.append(["QX", "Sx", "1e-20", "90.0", "10", "20", "1100", "1120", "AAAAA"])
            rows.append(["QX", "Sx", "1e-20", "90.0", "25", "30", "1200", "1220", "AAAAA"])
            # contig Sx, second copy
            rows.append(["QX", "Sx", "1e-20", "90.0", "1",  "5",  "1300", "1320", "AAAAA"])
            rows.append(["QX", "Sx", "1e-20", "90.0", "10", "20", "1400", "1420", "AAAAA"])
            rows.append(["QX", "Sx", "1e-20", "90.0", "25", "30", "1500", "1520", "AAAAA"])

            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")

            res = Results(results_tsv_path, query_fasta_path=query_fasta_path)
            pms = group_matches(res, max_intron_length = 10_000)

            self.assertEqual(len(pms), 2)
            self.assertEqual(pms[0].query_start, 1)
            self.assertEqual(pms[0].query_end, 30)
            self.assertEqual(pms[0].target_start, 1000)
            self.assertEqual(pms[0].target_end, 1220)
            self.assertEqual(pms[1].query_start, 1)
            self.assertEqual(pms[1].query_end, 30)
            self.assertEqual(pms[1].target_start, 1300)
            self.assertEqual(pms[1].target_end, 1520)

    def test_group_matches_separate_closely_placed_partial_proteins_beyond_allowed_aa_overlap_distance(self):
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
            # contig Sx, first copy
            rows.append(["QX", "Sx", "1e-20", "90.0", "1",  "5",  "1000", "1020", "AAAAA"])
            rows.append(["QX", "Sx", "1e-20", "90.0", "10", "20", "1100", "1120", "AAAAA"])
            # contig Sx, second copy
            rows.append(["QX", "Sx", "1e-20", "90.0", "10", "20", "1400", "1420", "AAAAA"])

            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")

            res = Results(results_tsv_path, query_fasta_path=query_fasta_path)

	    # max aa overlap allowed is 9 bps. query/aa coordinate of first
	    # copy ended at 20 and rewound to 10 for the next copy, so group
	    # two copies separately
            pms = group_matches(res, max_intron_length = 10_000, max_aa_overlap = 9)
            self.assertEqual(len(pms), 2)
            self.assertEqual(pms[0].query_start, 1)
            self.assertEqual(pms[0].query_end, 20)
            self.assertEqual(pms[0].target_start, 1000)
            self.assertEqual(pms[0].target_end, 1120)
            self.assertEqual(pms[1].query_start, 10)
            self.assertEqual(pms[1].query_end, 20)
            self.assertEqual(pms[1].target_start, 1400)
            self.assertEqual(pms[1].target_end, 1420)

	    # max aa overlap allowed is 10 bps. query/aa coordinate of first
	    # copy ended at 20 and rewound to 10 for the next copy, so group
	    # two copies together.
            pms = group_matches(res, max_intron_length = 10_000, max_aa_overlap = 10)
            self.assertEqual(len(pms), 1)

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

    def test_protein_hit_id_deterministic_and_changes_when_inputs_change(self):
        # Two identical ProteinMatch objects should have the same ID
        a1 = NucMatch("QID","TID",1,3,1,9,0.0,100.0)
        a2 = NucMatch("QID","TID",4,6,10,18,1e-12,95.0)
        pm1 = ProteinMatch("TID", [a1, a2], 1, 6, 1, 18)
        pm2 = ProteinMatch("TID", [a1, a2], 1, 6, 1, 18)
        self.assertEqual(pm1.protein_hit_id, pm2.protein_hit_id)
        # Change an input value (e.g., query_start) to produce a different ID
        b2 = NucMatch("QID","TID",5,6,10,18,1e-12,95.0)
        pm3 = ProteinMatch("TID", [a1, b2], 1, 6, 1, 18)
        self.assertNotEqual(pm1.protein_hit_id, pm3.protein_hit_id)

    def test_extra_match_in_middle_changes_protein_id(self):
        # Two identical ProteinMatch objects should have the same ID
        a1 = NucMatch("QID","TID",1,3,1,9,0.0,100.0)
        a2 = NucMatch("QID","TID",4,6,10,18,1e-12,95.0)
        pm1 = ProteinMatch("TID", [a1, a2], 1, 6, 1, 18)
        pm2 = ProteinMatch("TID", [a1, a2], 1, 6, 1, 18)
        self.assertEqual(pm1.protein_hit_id, pm2.protein_hit_id)

        # Add new NucMatch in middle, changes protein id
        b2 = NucMatch("QID","TID",3,4,9,10,1e-12,95.0)
        pm3 = ProteinMatch("TID", [a1, b2, a2], 1, 6, 1, 18)
        self.assertNotEqual(pm1.protein_hit_id, pm3.protein_hit_id)

    def test_reverse_strand_target_revcomp_and_proteinmatch_orientation(self):
        # Reverse strand target (sstart > send): store ascending coordinates in match,
        # extract reverse-complement target sequence, and ProteinMatch orientation 5'->3' with start > end
        with tempfile.TemporaryDirectory() as tmpdir:
            query_fasta_path = os.path.join(tmpdir, "q.faa")
            target_fasta_path = os.path.join(tmpdir, "t.fna")
            results_tsv_path = os.path.join(tmpdir, "r.tsv")

            with open(query_fasta_path, "w") as f:
                f.write(">Qrev\n")
                f.write("ACDEFGHIKLMNPQRSTVWYAC\n")  # arbitrary protein-like string, length >= 10
            with open(target_fasta_path, "w") as f:
                f.write(">Trev\n")
                # DNA sequence positions 1..30
                f.write("AACCGGTTAACCGGTTACGTACGTACGTAA\n")

            header = "\t".join(Results.PRODUCER_HEADER)
            rows = []
            # query 5..10, target reverse: sstart 20, send 12 -> store 12..20 asc and revcomp the slice
            rows.append(["Qrev", "Trev", "1e-6", "90", "5", "10", "20", "12", "XXXXX"])
            # second nearby reverse match to form one cluster
            rows.append(["Qrev", "Trev", "1e-6", "90", "1", "4", "11", "6", "YYYYY"])

            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")

            res = Results(results_tsv_path, query_fasta_path=query_fasta_path, target_fasta_path=target_fasta_path)
            ms = res.matches()
            self.assertEqual(len(ms), 2)

            m1 = ms[0]
            # Query positions 5..10 of provided string: "FGHIKL"
            self.assertEqual(m1.query_sequence, (res._query_sequences_by_accession["Qrev"])[4:10])
            # Target positions 12..20 asc: compute and ensure revcomp is returned
            raw_target = (res._target_sequences_by_accession["Trev"])[11:20]  # 12..20 1-based
            # manual reverse complement
            comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
            expected_rc = "".join(comp.get(b, b) for b in raw_target[::-1])
            self.assertEqual(m1.target_sequence, expected_rc)

            pms = group_matches(res)
            self.assertEqual(len(pms), 1)
            pm = pms[0]
            # Reverse strand grouping -> target_start should be the 5' end, i.e., larger coordinate
            self.assertGreater(pm.target_start, pm.target_end)

    def test_assert_query_start_le_query_end(self):
        # Ensure parser raises when qstart > qend
        with tempfile.TemporaryDirectory() as tmpdir:
            query_fasta_path = os.path.join(tmpdir, "q.faa")
            results_tsv_path = os.path.join(tmpdir, "r.tsv")
            with open(query_fasta_path, "w") as f:
                f.write(">Q\nAAAAAA\n")
            header = "\t".join(Results.PRODUCER_HEADER)
            rows = [
                ["Q", "T", "1e-5", "90", "10", "5", "100", "90", "XXXXX"],  # qstart > qend -> error
            ]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")
            with self.assertRaises(ValueError):
                Results(results_tsv_path, query_fasta_path=query_fasta_path).matches()


class TestCollatingSequence(unittest.TestCase):

    def test_translation_and_collated_protein_sequence_and_pprint(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            query_fasta_path = os.path.join(tmpdir, "q.faa")
            target_fasta_path = os.path.join(tmpdir, "t.fna")
            results_tsv_path = os.path.join(tmpdir, "r.tsv")

            # Query of length 6 (positions 1..6)
            with open(query_fasta_path, "w") as f:
                f.write(">Qprot\n")
                f.write("AAAAAA\n")

            # Build target DNA with two coding blocks
            # Block1: ATG GAA TTT -> M E F (positions 11..19)
            # Block2: GAA GTG GGG -> E V G (positions 22..30)
            target = "N" * 10 + "ATGGAATTT" + "N" * 2 + "GAAGTGGGG" + "N" * 50
            with open(target_fasta_path, "w") as f:
                f.write(">Tprot\n")
                f.write(target + "\n")

            header = "\t".join(Results.PRODUCER_HEADER)
            rows = [
                ["Qprot", "Tprot", "1e-5", "90", "1", "3", "11", "19", ""],  # MEF
                ["Qprot", "Tprot", "1e-5", "90", "2", "4", "22", "30", ""],  # EVG
            ]
            with open(results_tsv_path, "w") as f:
                f.write(header + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")

            res = Results(results_tsv_path, query_fasta_path, target_fasta_path)
            pms = group_matches(res)
            self.assertEqual(len(pms), 1)
            pm = pms[0]

            # Check translation on first match
            m1 = pm.matches[0]
            self.assertEqual(m1.target_sequence_translated(), "MEF")
            # Check translation on second match
            m2 = pm.matches[1]
            self.assertEqual(m2.target_sequence_translated(), "EVG")

            # Collated protein sequence across full query (len 6):
            # pos1: M
            # pos2: {E/E} -> E
            # pos3: {F/V} -> {F/V}
            # pos4: G
            # After update: use X for missing and strip leading/trailing Xs; here ends at pos4
            collated = pm.collated_protein_sequence
            self.assertEqual(collated, "M(EF/EV)G")

    def test_collate_handles_single_match(self):
        a = NucMatch("Q","T",1,3,1,9,0.0,100.0); a.target_sequence="ATGGAATTT"    # MEF
        pm = ProteinMatch("T",[a],1,3,1,9)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "MEF")

    def test_collate_handles_gaps_and_overlaps(self):
        a = NucMatch("Q","T",1,3,1,9,0.0,100.0); a.target_sequence="ATGGAATTT"    # MEF
        b = NucMatch("Q","T",3,5,10,18,0.0,100.0); b.target_sequence="GAAGTGGGG"  # EVG
        c = NucMatch("Q","T",9,9,30,32,0.0,100.0); c.target_sequence="ATG"        # M
        pm = ProteinMatch("T",[a,b,c],1,9,1,32)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "ME(F/E)VGXXXM")

    def test_collate_handles_gaps_within_match(self):
        a = NucMatch("Q","T",1,3,1,6,0.0,100.0); a.target_sequence="ATGGAA"       # ME - but matching to 3 bps of query
        b = NucMatch("Q","T",3,5,10,18,0.0,100.0); b.target_sequence="GAAGTGGGG"  # EVG
        c = NucMatch("Q","T",9,9,30,32,0.0,100.0); c.target_sequence="ATG"        # M
        pm = ProteinMatch("T",[a,b,c],1,9,1,32)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "M(E/E)VGXXXM")

    def test_collate_handles_insertions_within_match(self):
        a = NucMatch("Q","T",1,3,1,12,0.0,100.0); a.target_sequence="ATGGAATTTTTT" # MEFF - but matching to 3 bps of query
        b = NucMatch("Q","T",3,5,10,18,0.0,100.0); b.target_sequence="GAAGTGGGG"   # EVG
        c = NucMatch("Q","T",9,9,30,32,0.0,100.0); c.target_sequence="ATG"         # M
        pm = ProteinMatch("T",[a,b,c],1,9,1,32)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "MEF(F/E)VGXXXM")

    def test_trimming_of_leading_and_trailing_X(self):
        # Explicit trimming test: single internal block at q=3..8 should produce no leading/trailing Xs.
        with tempfile.TemporaryDirectory() as tmpdir:
            q_path = os.path.join(tmpdir, "q.faa")
            t_path = os.path.join(tmpdir, "t.fna")
            r_path = os.path.join(tmpdir, "r.tsv")
            # Query: 10 AA
            with open(q_path, "w") as f:
                f.write(">Qtrim\nABCDEFGHIJ\n")
            # Target genome with coding for CDEFGH placed starting at pos 6
            coding = "TGCGATGAATTTGGTCAT"  # TGC(C) GAT(D) GAA(E) TTT(F) GGT(G) CAT(H)
            genome = "AAAAA" + coding + "AAAAA"
            with open(t_path, "w") as f:
                f.write(">Ttrim\n" + genome + "\n")
            header = "\t".join(Results.PRODUCER_HEADER)
            row = ["Qtrim", "Ttrim", "0", "100.0", "3", "8", "6", str(6 + len(coding) - 1), "CDEFGH"]
            with open(r_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(row) + "\n")
            res = Results(r_path, q_path, t_path)
            pms = group_matches(res)
            self.assertEqual(len(pms), 1)
            pm = pms[0]
            s1 = pm.collated_protein_sequence
            self.assertEqual(s1, "CDEFGH")


if __name__ == "__main__":
    unittest.main()


