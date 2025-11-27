import os
import tempfile
import unittest

from needle.blast import Results, group_matches, ProteinMatch
from needle.blast import order_matches_for_junctions, generate_transition_candidates, score_and_select_best_transition, stitch_cleaned_sequence, Candidate
import shutil
import types
import needle.blast as blast_mod
import tempfile as _tempfile


class TestBlastResults(unittest.TestCase):
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
            self.assertEqual(m.target_start, 7)
            self.assertEqual(m.target_end, 10)

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

    def test_flanking_sequences_forward_and_reverse(self):
        # Validate upstream/downstream flanks for both forward and reverse matches
        with tempfile.TemporaryDirectory() as tmpdir:
            q_path = os.path.join(tmpdir, "q.faa")
            t_path = os.path.join(tmpdir, "t.fna")
            r_path = os.path.join(tmpdir, "r.tsv")

            # Query protein length >= 6 so we can map 2 codons easily
            with open(q_path, "w") as f:
                f.write(">Q\n")
                f.write("AAAAAA\n")

            # Target genome: 1..40
            # We'll place a forward block at 11..16 (6bp) and reverse block at 31..36 (6bp)
            # Flank window = 4
            genome = "NNNNN" + "AAAACG" + "NNNNN" + "TTGCAA" + "NNNNN" + "GGGGG"  # just a sequence
            # Forward target [11..16] = "AAAACG"
            # Reverse target [31..36] = "GGGGG" -> but length 5; adjust to "GGGGNN" not ideal. Use a precise genome:
            genome = "AAAAAAAAAATGGAATTTCCCCCTTCCAAGGTTNNNN"
            # positions:
            #  1..10  = AAAAAAAAAA
            # 11..19  = TGGAATTT (9bp) -> MEF
            # 20..24  = CCCCC
            # 25..30  = TTCCAA
            # 31..35  = GGTTN
            # Use reverse block at 25..30 (TTCCAA) so RC = TTGGAA -> translates as L? We'll only test flanks as DNA
            with open(t_path, "w") as f:
                f.write(">T\n")
                f.write(genome + "\n")

            header = "\t".join(Results.PRODUCER_HEADER)
            rows = [
                # Forward match: 11..19, query 1..3 (MEF)
                ["Q", "T", "1e-5", "90", "1", "3", "11", "19", ""],
                # Reverse match: sstart > send (30..25), query 4..5
                ["Q", "T", "1e-5", "90", "4", "5", "30", "25", ""],
            ]
            with open(r_path, "w") as f:
                f.write(header + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")

            # Use flank window = 4
            res = Results(r_path, q_path, t_path, target_sequence_flanking_window=4)
            ms = res.matches()
            self.assertEqual(len(ms), 2)
            m_fwd = ms[0]
            m_rev = ms[1]

            # Forward expected: genome[7..23] because target 11..19 and flank 4 -> (11-4)=7, (19+4)=23
            expected_forward_concat = genome[6:23]  # 1-based -> 0-based slice
            self.assertEqual(
                (m_fwd.target_sequence_upstream or "") + (m_fwd.target_sequence or "") + (m_fwd.target_sequence_downstream or ""),
                expected_forward_concat
            )

            # Reverse expected: take genome[21..34] (since 25..30 +/-4) and reverse-complement it
            # (25-4)=21, (30+4)=34
            from Bio.Seq import Seq as _Seq
            expected_rev_concat = str(_Seq(genome[20:34]).reverse_complement())
            self.assertEqual(
                (m_rev.target_sequence_upstream or "") + (m_rev.target_sequence or "") + (m_rev.target_sequence_downstream or ""),
                expected_rev_concat
            )

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

    # pprint method removed; related tests omitted.
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
            self.assertEqual(collated, "ME{F/V}G")

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

    # New unit tests for junction optimization procedures
    def test_order_matches_for_junctions_overlap_and_gap(self):
        # Build synthetic matches with query positions only; no file IO
        class _M:
            pass
        m1 = _M(); m1.query_start=1; m1.query_end=10
        m2 = _M(); m2.query_start=8; m2.query_end=15
        m3 = _M(); m3.query_start=18; m3.query_end=20
        # We need NucMatch-like, but only .query_start/.query_end are used
        pairs = order_matches_for_junctions([m1, m2, m3])  # type: ignore
        # Adjacent pairs: (m1,m2) overlap 3 (10-8+1), (m2,m3) gap 2 (18-15-1)
        self.assertEqual(len(pairs), 2)
        self.assertEqual(pairs[0][2], 3)  # overlap_len
        self.assertEqual(pairs[0][3], 0)  # gap_len
        self.assertEqual(pairs[1][2], 0)
        self.assertEqual(pairs[1][3], 2)

    def test_generate_transition_candidates_overlap_and_gap(self):
        left = "ABCDEFG"
        right = "DEFGHIJ"
        # Use overlap_len=3 to enumerate 4 splits k=0..3
        cands = generate_transition_candidates(left, right, overlap_len=3, gap_len=0, flank_window_max_length=2)
        self.assertEqual(len(cands), 4)  # L+1
        # Validate stitched sequences for each k
        left_len = len(left)
        right_len = len(right)
        for k, cand in enumerate(cands):
            exp_left_prefix = left[0:left_len - k]
            exp_right_suffix = right[k:right_len]
            exp_stitched = exp_left_prefix + exp_right_suffix
            self.assertEqual(cand.stitched_local, exp_stitched)
            # window is last 2 of left_prefix + first 2 of right_suffix
            self.assertEqual(cand.window_seq, exp_left_prefix[-2:] + exp_right_suffix[:2])
        # Gap: single candidate with Xs
        cands_gap = generate_transition_candidates("AAA", "BBB", overlap_len=0, gap_len=2, flank_window_max_length=2)
        self.assertEqual(len(cands_gap), 1)
        self.assertIn("XX", cands_gap[0].stitched_local)

    def test_score_and_select_best_transition_fallback(self):
        # Ensure fallback deterministic heuristic works without hmmsearch
        c1 = Candidate(split_k=0, window_seq="AAAA", stitched_local="AAA")
        c2 = Candidate(split_k=1, window_seq="ZZZZ", stitched_local="AAAA")
        # c2 has longer stitched_local, should be picked
        best = score_and_select_best_transition([c1, c2], hmm_file_name="no_such_hmm.hmm")
        self.assertIs(best, c2)

    def test_stitch_cleaned_sequence_basic(self):
        # Two blocks, overlap 2 AA. Select a split to keep left fully (k=0)
        class _MM:
            pass
        left = _MM(); right = _MM()
        left.query_start=1; left.query_end=5
        right.query_start=4; right.query_end=8
        # Provide aa_by_match map
        aa_map = {id(left):"ABCDE", id(right):"DEFGH"}
        pairs = order_matches_for_junctions([left, right])  # type: ignore
        self.assertEqual(pairs[0][2], 2)
        # Choose split k=2 so overlap removed -> result ABCDE + FGH
        cand = Candidate(split_k=2, window_seq="", stitched_local="ABCDEFGH")
        stitched = stitch_cleaned_sequence([left, right], {0: cand}, aa_map)  # type: ignore
        self.assertEqual(stitched, "ABCDEFGH")

    def test_stitch_cleaned_sequence_multiple_blocks_mixed(self):
        # Three blocks: first-second overlap 2, second-third gap 3
        class _MM:
            pass
        a=_MM(); b=_MM(); c=_MM()
        a.query_start=1; a.query_end=5      # "ABCDE"
        b.query_start=4; b.query_end=9      # "DEFGHI"
        c.query_start=13; c.query_end=15    # "KLM" (gap of 3 AA from 10..12)
        aa_map = {id(a):"ABCDE", id(b):"DEFGHI", id(c):"KLM"}
        pairs = order_matches_for_junctions([a,b,c])  # type: ignore
        # pairs[0]: overlap 2; choose k=1 (keep A..F + GHI)
        cand0 = Candidate(split_k=1, window_seq="", stitched_local="ABCDEFGHI")
        # pairs[1]: gap 3; candidate should be left (DEFGHI) + XXX + right (KLM)
        cand1 = Candidate(split_k=None, window_seq="", stitched_local="DEFGHIXXXKLM")
        stitched = stitch_cleaned_sequence([a,b,c], {0:cand0, 1:cand1}, aa_map)  # type: ignore
        self.assertEqual(stitched, "ABCDEFGHIXXXKLM")

    def test_score_and_select_best_transition_with_mocked_hmmsearch_domtbl(self):
        # Mock tempfile.TemporaryDirectory to a known path and subprocess.run to write domtblout
        import tempfile as _tempfile
        import subprocess as _subprocess
        import shutil as _shutil
        tmp_root = _tempfile.mkdtemp()

        class _FakeTD:
            def __enter__(self_):
                return tmp_root
            def __exit__(self_, exc_type, exc, tb):
                _shutil.rmtree(tmp_root, ignore_errors=True)

        def _fake_run(cmd, check, stdout, stderr):
            # Write domtbl with two lines; make cand_1 have higher score at parts[13]
            domtbl_path = os.path.join(tmp_root, "out.domtbl")
            with open(domtbl_path, "w") as f:
                # columns: ensure index 13 exists; parts[0]=target name 'cand_i'
                def line(name, score):
                    parts = [name] + ["x"] * 12 + [str(score)] + ["x"] * 5
                    return " ".join(parts) + "\n"
                f.write(line("cand_0", 10.0))
                f.write(line("cand_1", 50.0))
            class _P: returncode = 0
            return _P()

        orig_td = _tempfile.TemporaryDirectory
        orig_run = _subprocess.run
        try:
            _tempfile.TemporaryDirectory = _FakeTD  # type: ignore
            _subprocess.run = _fake_run  # type: ignore
            c1 = Candidate(split_k=0, window_seq="AAAA", stitched_local="LEFT")
            c2 = Candidate(split_k=1, window_seq="BBBB", stitched_local="RIGHT")
            best = score_and_select_best_transition([c1, c2], hmm_file_name="ignored.hmm")
            self.assertIs(best, c2)  # cand_1 has higher score in mocked domtbl
        finally:
            _tempfile.TemporaryDirectory = orig_td  # type: ignore
            _subprocess.run = orig_run  # type: ignore

    def test_hmm_cleaned_target_protein_sequence_integration_with_mock_scoring(self):
        # Build three synthetic NucMatch blocks directly: overlap then gap
        from needle.blast import NucMatch, ProteinMatch
        a = NucMatch(query_accession="Q", target_accession="T",
                     query_start=1, query_end=3, target_start=1, target_end=9,
                     e_value=0.0, identity=100.0, on_reverse_strand=False)
        b = NucMatch(query_accession="Q", target_accession="T",
                     query_start=3, query_end=5, target_start=10, target_end=18,
                     e_value=0.0, identity=100.0, on_reverse_strand=False)
        c = NucMatch(query_accession="Q", target_accession="T",
                     query_start=9, query_end=9, target_start=30, target_end=32,
                     e_value=0.0, identity=100.0, on_reverse_strand=False)
        # Set DNA sequences so translations are known
        a.target_sequence = "ATGGAATTT"      # MEF
        b.target_sequence = "GAAGTGGGG"      # EVG
        c.target_sequence = "ATG"            # M
        pm = ProteinMatch(
            target_id="T",
            matches=[a, b, c],
            query_start=1,
            query_end=9,
            target_start=1,
            target_end=32,
            covers_start_to_end=False,
            likely_complete=False,
            query_overlap=True,
        )
        # Monkeypatch scoring to always choose split_k=1 on overlaps
        orig = blast_mod.score_and_select_best_transition
        def _fake_scorer(cands, hmm_file_name):
            for c in cands:
                if c.split_k == 1:
                    return c
            return cands[0]
        try:
            blast_mod.score_and_select_best_transition = _fake_scorer  # type: ignore
            cleaned = pm.hmm_cleaned_target_protein_sequence("dummy.hmm", flank_window_max_length=5)
        finally:
            blast_mod.score_and_select_best_transition = orig  # type: ignore
        # Assertions:
        # a) leading/trailing X removed: our sequence starts at pos1 and ends at pos9; trimmed by method -> no leading/trailing X
        self.assertFalse(cleaned.startswith("X"))
        self.assertFalse(cleaned.endswith("X"))
        # b) gap (q6..q8) is filled with X
        self.assertIn("XXX", cleaned)
        # c) best transition (k=1) used around overlap between first two blocks:
        # k=1 makes left_prefix 'ME' + right_suffix 'VG' => "MEVG" at the start
        self.assertTrue(cleaned.startswith("MEVG"))
        # Expected final stitched AA: "MEF"(1..3) and "EVG"(3..5) overlap by 1 AA at pos3.
        # Choosing k=1 yields "MEVG" for the first junction, and the 3-AA gap (6..8) is filled with "XXX",
        # followed by "M" at pos9 from the third block => "MEVGXXXM".
        self.assertEqual(cleaned, "MEVGXXXM")

    def test_group_matches_sets_hmm_cleaned_sequence_when_hmm_dir(self):
        # Use mocked scoring to avoid calling hmmsearch; Results should pass hmm_directory to grouping.
        with tempfile.TemporaryDirectory() as d:
            hmm_dir = os.path.join(d, "hmms")
            os.makedirs(hmm_dir, exist_ok=True)
            # Create dummy HMM file named after query accession
            hmm_path = os.path.join(hmm_dir, "Qhmm.hmm")
            with open(hmm_path, "w") as f:
                f.write("HMMER3/f")  # minimal content
            # Build two overlapping matches to force multiple candidates
            from needle.blast import NucMatch, ProteinMatch
            a = NucMatch(query_accession="Qhmm", target_accession="T",
                         query_start=1, query_end=3, target_start=1, target_end=9,
                         e_value=0.0, identity=100.0, on_reverse_strand=False)
            b = NucMatch(query_accession="Qhmm", target_accession="T",
                         query_start=3, query_end=5, target_start=10, target_end=18,
                         e_value=0.0, identity=100.0, on_reverse_strand=False)
            a.target_sequence = "ATGGAATTT"  # MEF
            b.target_sequence = "GAAGTGGGG"  # EVG
            # Create a Results with hmm_directory, but we won't use file IO; instead monkeypatch scoring.
            # We'll construct a fake results.tsv to satisfy parser with one cluster
            r_path = os.path.join(d, "r.tsv")
            header = "\t".join(Results.PRODUCER_HEADER)
            with open(r_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(["Qhmm","T","0","100","1","3","1","9",""]) + "\n")
                f.write("\t".join(["Qhmm","T","0","100","3","5","10","18",""]) + "\n")
            # Also minimal query/target to enable translation (not needed here since we set sequences)
            q_path = os.path.join(d, "q.faa")
            t_path = os.path.join(d, "t.fna")
            with open(q_path, "w") as f:
                f.write(">Qhmm\nMEFEVG\n")
            with open(t_path, "w") as f:
                # positions 1..9: ATGGAATTT -> MEF; 10..18: GAAGTGGGG -> EVG
                f.write(">T\nATGGAATTTGAAGTGGGG\n")
            # Monkeypatch the score function
            orig = blast_mod.score_and_select_best_transition
            def _fake_scorer(cands, hmm_file_name):
                for c in cands:
                    if c.split_k == 1:
                        return c
                return cands[0]
            try:
                blast_mod.score_and_select_best_transition = _fake_scorer  # type: ignore
                res = Results(r_path, q_path, t_path, hmm_directory=hmm_dir)
                pms = group_matches(res)
                self.assertEqual(len(pms), 1)
                pm = pms[0]
                self.assertEqual(pm.hmm_cleaned_protein_sequence, "MEVG")
            finally:
                blast_mod.score_and_select_best_transition = orig  # type: ignore

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
            s2 = pm.hmm_cleaned_target_protein_sequence("dummy.hmm", flank_window_max_length=10)
            self.assertEqual(s2, "CDEFGH")

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


if __name__ == "__main__":
    unittest.main()


