import os
import tempfile
import unittest

from needle.blast import Results
from needle.match import group_matches, ProteinMatch, NucMatch, order_matches_for_junctions, NonlinearMatchException
from needle.match import extract_subsequence, extract_subsequence_strand_sensitive


class TestOrderGroupMatches(unittest.TestCase):

    @staticmethod
    def makeM(query_start, query_end, target_start, target_end):
        return NucMatch(
            query_accession=None, target_accession=None, e_value=0, identity=None,
            query_start=query_start, query_end=query_end, target_start=target_start, target_end=target_end)

    def test_order_matches_for_junctions_overlap_and_gap(self):
        m1 = self.makeM(1, 10, 1, 10)
        m2 = self.makeM(8, 15, 8, 15)
        m3 = self.makeM(18, 20, 18, 20)

        pairs = order_matches_for_junctions([m1, m3, m2])  # input not in order
        self.assertEqual(len(pairs), 2)
        self.assertEqual(pairs[0], (m1, m2, 3, 0))
        self.assertEqual(pairs[1], (m2, m3, 0, 2))

    def test_order_throws_error_if_junctions_overlap(self):
        # aaaaaaaa
        #      bbbbbb
        #        cccccccc

        m1 = self.makeM(1, 8, 1, 8)
        m2 = self.makeM(6, 12, 6, 12)
        m3 = self.makeM(8, 16, 8, 16)

        with self.assertRaises(NonlinearMatchException):
            pairs = order_matches_for_junctions([m1, m3, m2])  # input not in order

    def test_order_throws_error_on_contained_match(self):
        # aaaaaaaaa
        #      bbbbbbbbbbb
        #            ccc
        #                ddd

        m1 = self.makeM(1, 10, 1, 10)
        m2 = self.makeM(6, 16, 6, 16)
        m3 = self.makeM(12, 14, 12, 14)
        m4 = self.makeM(16, 18, 16, 18)

        with self.assertRaises(NonlinearMatchException):
            pairs = order_matches_for_junctions([m1, m4, m3, m2])  # input not in order

    def test_order_throws_error_if_query_coordinates_do_not_align_with_target_coordinates(self):

        m1 = self.makeM(1, 10, 1, 10)
        m2 = self.makeM(6, 16, 16, 10)

        with self.assertRaises(NonlinearMatchException):
            pairs = order_matches_for_junctions([m1, m2])

    def test_order_throws_error_if_query_order_is_not_same_as_target_order(self):

        m1 = self.makeM(1, 10, 12, 18)
        m2 = self.makeM(6, 16, 1, 10)

        with self.assertRaises(NonlinearMatchException):
            pairs = order_matches_for_junctions([m1, m2])

    def test_order_throws_error_if_query_overlap_is_too_large(self):

        m1 = self.makeM(6, 16, 1, 10)
        m2 = self.makeM(10, 20, 12, 22)

        pairs = order_matches_for_junctions([m1, m2], max_overlap_len=10)
        with self.assertRaises(NonlinearMatchException):
            pairs = order_matches_for_junctions([m1, m2], max_overlap_len=4)

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
            pms = group_matches(res.matches())

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
            pms = group_matches(res.matches(), max_intron_length = 10_000)

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
            pms = group_matches(res.matches(), max_intron_length = 10_000, max_aa_overlap = 9)
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
            pms = group_matches(res.matches(), max_intron_length = 10_000, max_aa_overlap = 10)
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
            pms = group_matches(res.matches(), max_intron_length=10000)
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
            pms = group_matches(res.matches(), max_intron_length=10000)
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
            pms = group_matches(res.matches(), max_intron_length=20)
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

            pms = group_matches(res.matches())
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
            pms = group_matches(res.matches())
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
            pms = group_matches(res.matches())
            self.assertEqual(len(pms), 1)
            pm = pms[0]
            s1 = pm.collated_protein_sequence
            self.assertEqual(s1, "CDEFGH")

if __name__ == "__main__":
    unittest.main()

class TestExtractSubsequence(unittest.TestCase):

    def test_forward_in_bounds(self):
        seq = "ABCDEFGH"  # 1..8
        # 2..5 => BCDE
        self.assertEqual(extract_subsequence(seq, 2, 5), "BCDE")
        # 1..8 => full
        self.assertEqual(extract_subsequence(seq, 1, 8), "ABCDEFGH")

    def test_reversed_coords_same_slice(self):
        seq = "ABCDEFGH"
        # reversed coordinates (5..2) should return same slice as (2..5)
        self.assertEqual(extract_subsequence(seq, 5, 2), "BCDE")

    def test_end_beyond_length_is_clamped(self):
        seq = "ABCDEFGH"  # len=8
        # 6..12 => 6..8 => FGH
        self.assertEqual(extract_subsequence(seq, 6, 12), "FGH")
        # reversed with beyond-length
        self.assertEqual(extract_subsequence(seq, 12, 6), "FGH")

    def test_start_beyond_length_returns_none(self):
        seq = "ABCDEFGH"  # len=8
        # left (min) > len(seq) -> None
        self.assertIsNone(extract_subsequence(seq, 9, 10))
        self.assertIsNone(extract_subsequence(seq, 10, 9))

    def test_invalid_coords_non_positive(self):
        seq = "ABCDEFGH"
        self.assertIsNone(extract_subsequence(seq, 0, 5))
        self.assertIsNone(extract_subsequence(seq, 5, 0))
        self.assertIsNone(extract_subsequence(seq, -1, 2))

    def test_none_full_sequence_returns_none(self):
        self.assertIsNone(extract_subsequence(None, 1, 3))

    def test_empty_sequence_returns_none(self):
        # len("") == 0, left=min(1,1)=1 > 0 -> None
        self.assertIsNone(extract_subsequence("", 1, 1))

    def test_reversed_coords_strand_sensitie(self):
        seq = "TAGTCAAA"
        # reversed coordinates (5..2) should return same slice as (2..5)
        self.assertEqual(extract_subsequence_strand_sensitive(seq, 5, 2), "GACT")


