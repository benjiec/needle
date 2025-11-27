import os
import tempfile
import unittest

from needle.hits import (
    order_matches_for_junctions,
    generate_transition_candidates,
    score_and_select_best_transition,
    stitch_cleaned_sequence,
    Candidate,
    hmm_clean_protein,
    hmm_clean,
    adjust_target_coordinates,
)
import needle.hits as hits_mod

from needle.blast import Results, group_matches, ProteinMatch, NucMatch


class TestHits(unittest.TestCase):
    def test_order_matches_for_junctions_overlap_and_gap(self):
        class _M: pass
        m1 = _M(); m1.query_start=1; m1.query_end=10
        m2 = _M(); m2.query_start=8; m2.query_end=15
        m3 = _M(); m3.query_start=18; m3.query_end=20
        pairs = order_matches_for_junctions([m1, m2, m3])  # type: ignore
        self.assertEqual(len(pairs), 2)
        self.assertEqual(pairs[0][2], 3)
        self.assertEqual(pairs[0][3], 0)
        self.assertEqual(pairs[1][2], 0)
        self.assertEqual(pairs[1][3], 2)

    def test_generate_transition_candidates_overlap_and_gap(self):
        # Use explicit overlap example with conceptual overlap length 4
        left = "ABCDEFX"
        right = "yefghij"
        cands = generate_transition_candidates(left, right, overlap_len=4, gap_len=0, flank_window_max_length=2)
        self.assertEqual(len(cands), 5)
        # Explicitly assert stitched sequences for visual verification
        stitched = [c.stitched_local for c in cands]
        self.assertEqual(
            stitched,
            [
                "ABCyefghij",     # k=0: left trims 4, right trims 0
                "ABCDefghij",     # k=1: left trims 3, right trims 1
                "ABCDEfghij",     # k=2: left trims 2, right trims 2
                "ABCDEFghij",     # k=3: left trims 1, right trims 3
                "ABCDEFXhij",     # k=4: left trims 0, right trims 4
            ],
        )
        # Window seq still formed from tail(left_prefix,2) + head(right_suffix,2)
        windows = [c.window_seq for c in cands]
        self.assertEqual(
            windows,
            [
                "BCye",  # "BC" + "YE"
                "CDef",  # "CD" + "EF"
                "DEfg",  # "DE" + "FG"
                "EFgh",  # "EF" + "GH"
                "FXhi",  # "FX" + "HI"
            ],
        )
        cands_gap = generate_transition_candidates("AAA", "BBB", overlap_len=0, gap_len=2, flank_window_max_length=2)
        self.assertEqual(len(cands_gap), 1)
        self.assertIn("XX", cands_gap[0].stitched_local)

    def test_score_and_select_best_transition_fallback(self):
        c1 = Candidate(split_k=0, window_seq="AAAA", stitched_local="AAA")
        c2 = Candidate(split_k=1, window_seq="ZZZZ", stitched_local="AAAA")
        best = score_and_select_best_transition([c1, c2], hmm_file_name="no_such_hmm.hmm")
        self.assertIs(best, c2)

    def test_stitch_cleaned_sequence_basic(self):
        class _MM: pass
        left = _MM(); right = _MM()
        left.query_start=1; left.query_end=5
        right.query_start=4; right.query_end=8
        aa_map = {id(left):"ABCDE", id(right):"DEFGH"}
        pairs = order_matches_for_junctions([left, right])  # type: ignore
        self.assertEqual(pairs[0][2], 2)
        cand = Candidate(split_k=2, window_seq="", stitched_local="ABCDEFGH")
        stitched = stitch_cleaned_sequence([left, right], {0: cand}, aa_map)  # type: ignore
        self.assertEqual(stitched, "ABCDEFGH")

    def test_stitch_cleaned_sequence_multiple_blocks_mixed(self):
        class _MM: pass
        a=_MM(); b=_MM(); c=_MM()
        a.query_start=1; a.query_end=5
        b.query_start=4; b.query_end=9
        c.query_start=13; c.query_end=15
        aa_map = {id(a):"ABCDE", id(b):"DEFGHI", id(c):"KLM"}
        pairs = order_matches_for_junctions([a,b,c])  # type: ignore
        cand0 = Candidate(split_k=1, window_seq="", stitched_local="ABCDEFGHI")
        cand1 = Candidate(split_k=None, window_seq="", stitched_local="DEFGHIXXXKLM")
        stitched = stitch_cleaned_sequence([a,b,c], {0:cand0, 1:cand1}, aa_map)  # type: ignore
        self.assertEqual(stitched, "ABCDEFGHIXXXKLM")

    def test_score_and_select_best_transition_with_mocked_hmmsearch_domtbl(self):
        import tempfile as _tempfile
        import subprocess as _subprocess
        import shutil as _shutil
        tmp_root = _tempfile.mkdtemp()
        class _FakeTD:
            def __enter__(self_): return tmp_root
            def __exit__(self_, exc_type, exc, tb): _shutil.rmtree(tmp_root, ignore_errors=True)
        def _fake_run(cmd, check, stdout, stderr):
            domtbl_path = os.path.join(tmp_root, "out.domtbl")
            with open(domtbl_path, "w") as f:
                def line(name, score):
                    parts = [name] + ["x"] * 12 + [str(score)] + ["x"] * 5
                    return " ".join(parts) + "\n"
                f.write(line("cand_0", 10.0)); f.write(line("cand_1", 50.0))
            class _P: returncode = 0
            return _P()
        orig_td = _tempfile.TemporaryDirectory; orig_run = _subprocess.run
        try:
            _tempfile.TemporaryDirectory = _FakeTD  # type: ignore
            _subprocess.run = _fake_run  # type: ignore
            c1 = Candidate(split_k=0, window_seq="AAAA", stitched_local="LEFT")
            c2 = Candidate(split_k=1, window_seq="BBBB", stitched_local="RIGHT")
            best = score_and_select_best_transition([c1, c2], hmm_file_name="ignored.hmm")
            self.assertIs(best, c2)
        finally:
            _tempfile.TemporaryDirectory = orig_td  # type: ignore
            _subprocess.run = orig_run  # type: ignore

    def test_hmm_cleaned_protein_integration_with_mock_scoring(self):
        a = NucMatch("Q","T",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        b = NucMatch("Q","T",3,5,10,18,0.0,100.0,False); b.target_sequence="GAAGTGGGG"
        c = NucMatch("Q","T",9,9,30,32,0.0,100.0,False); c.target_sequence="ATG"
        pm = ProteinMatch("T",[a,b,c],1,9,1,32,False,False,True)
        orig = hits_mod.score_and_select_best_transition
        def _fake(cands, hmm): 
            for x in cands:
                if x.split_k == 1: return x
            return cands[0]
        try:
            hits_mod.score_and_select_best_transition = _fake  # type: ignore
            cleaned_pm = hmm_clean_protein(pm, "dummy.hmm", flank_window_max_length=5)
            cleaned = cleaned_pm.collated_protein_sequence
        finally:
            hits_mod.score_and_select_best_transition = orig  # type: ignore
        self.assertFalse(cleaned.startswith("X")); self.assertFalse(cleaned.endswith("X"))
        self.assertIn("XXX", cleaned)
        self.assertEqual(cleaned, "MEVGXXXM")

    def test_hmm_clean_protein_adjusts_overlap_coordinates(self):
        # Overlap: a(1..5), b(4..9) => overlap 2; choose k=1; c(13..15) should shift by 1
        # Use full-codon sequences to ensure AA coverage equals expected lengths
        a = NucMatch("Q","T",1,5,1,15,0.0,100.0,False); a.target_sequence="ATG"*5         # 'M'*5
        b = NucMatch("Q","T",4,9,16,33,0.0,100.0,False); b.target_sequence="GAA"*6        # 'E'*6
        c = NucMatch("Q","T",13,15,40,48,0.0,100.0,False); c.target_sequence="ATG"*3      # 'M'*3
        pm = ProteinMatch("T",[a,b,c],1,15,1,48,False,False,True)
        orig = hits_mod.score_and_select_best_transition
        def _fake(cands, hmm):
            for x in cands:
                if x.split_k == 1:
                    return x
            return cands[0]
        try:
            hits_mod.score_and_select_best_transition = _fake  # type: ignore
            cleaned_pm = hmm_clean_protein(pm, "dummy.hmm", flank_window_max_length=5)
        finally:
            hits_mod.score_and_select_best_transition = orig  # type: ignore
        nm = cleaned_pm.matches
        self.assertEqual(len(nm), 3)
        # a.end reduced by 1
        self.assertEqual(nm[0].query_start, 1)
        self.assertEqual(nm[0].query_end, 4)
        # b.start becomes a.end+1 = 5; b end remains 9 (no shift of downstream blocks)
        self.assertEqual(nm[1].query_start, 5)
        self.assertEqual(nm[1].query_end, 9)
        # c unchanged (no shifting of downstream matches)
        self.assertEqual(nm[2].query_start, 13)
        self.assertEqual(nm[2].query_end, 15)

    def test_hmm_clean_adjusts_overlap_coordinates_in_grouped(self):
        # Use hmm_clean on grouped matches and verify new coordinates
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            hmm_dir = os.path.join(d, "hmms"); os.makedirs(hmm_dir, exist_ok=True)
            with open(os.path.join(hmm_dir, "Qhmm.hmm"), "w") as f: f.write("HMMER3/f")
            # Build TSV with overlap: (1..5) and (4..9)
            r_path = os.path.join(d, "r.tsv"); q_path = os.path.join(d, "q.faa"); t_path = os.path.join(d, "t.fna")
            header = "\t".join(Results.PRODUCER_HEADER)
            with open(r_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(["Qhmm","T","0","100","1","5","1","15",""]) + "\n")
                f.write("\t".join(["Qhmm","T","0","100","4","9","16","33",""]) + "\n")
                f.write("\t".join(["Qhmm","T","0","100","13","15","40","48",""]) + "\n")
            with open(q_path, "w") as f: f.write(">Qhmm\n" + "A"*20 + "\n")
            # Provide target with full-codon blocks: positions 1..15 (ATG*5), 16..33 (GAA*6), 40..48 (ATG*3)
            with open(t_path, "w") as f:
                f.write(">T\n" + "ATG"*5 + "GAA"*6 + "N"*(39-33) + "ATG"*3 + "\n")
            orig = hits_mod.score_and_select_best_transition
            def _fake(cands, hmm):
                for x in cands:
                    if x.split_k == 1: return x
                return cands[0]
            try:
                hits_mod.score_and_select_best_transition = _fake  # type: ignore
                res = Results(r_path, q_path, t_path)
                pms = group_matches(res)
                cleaned = hmm_clean(pms, hmm_dir, flank_window_max_length=5)
                self.assertEqual(len(cleaned), 1)
                nm = cleaned[0].matches
                self.assertEqual((nm[0].query_start, nm[0].query_end), (1,4))
                self.assertEqual((nm[1].query_start, nm[1].query_end), (5,9))
                self.assertEqual((nm[2].query_start, nm[2].query_end), (13,15))
            finally:
                hits_mod.score_and_select_best_transition = orig  # type: ignore

    def test_hmm_clean_applies_to_grouped_matches(self):
        with tempfile.TemporaryDirectory() as d:
            hmm_dir = os.path.join(d, "hmms"); os.makedirs(hmm_dir, exist_ok=True)
            with open(os.path.join(hmm_dir, "Qhmm.hmm"), "w") as f: f.write("HMMER3/f")
            # Build simple dataset
            r_path = os.path.join(d, "r.tsv")
            header = "\t".join(Results.PRODUCER_HEADER)
            with open(r_path, "w") as f:
                f.write(header + "\n")
                f.write("\t".join(["Qhmm","T","0","100","1","3","1","9",""]) + "\n")
                f.write("\t".join(["Qhmm","T","0","100","3","5","10","18",""]) + "\n")
            q_path = os.path.join(d, "q.faa"); t_path = os.path.join(d, "t.fna")
            with open(q_path, "w") as f: f.write(">Qhmm\nMEFEVG\n")
            with open(t_path, "w") as f: f.write(">T\nATGGAATTTGAAGTGGGG\n")
            orig = hits_mod.score_and_select_best_transition
            def _fake(cands, hmm):
                for x in cands:
                    if x.split_k == 1: return x
                return cands[0]
            try:
                hits_mod.score_and_select_best_transition = _fake  # type: ignore
                res = Results(r_path, q_path, t_path)
                pms = group_matches(res)
                cleaned = hmm_clean(pms, hmm_dir, flank_window_max_length=5)
                self.assertEqual(len(cleaned), 1)
                self.assertEqual(cleaned[0].collated_protein_sequence, "MEVG")
            finally:
                hits_mod.score_and_select_best_transition = orig  # type: ignore

    def test_adjust_target_coordinates_gap_keeps_blocks(self):
        left = NucMatch("q","t",1,5,100,114,0.0,100.0,False); left.target_sequence="A"*15
        right = NucMatch("q","t",8,12,200,214,0.0,100.0,False); right.target_sequence="C"*15
        cand = Candidate(split_k=None, window_seq="", stitched_local="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        self.assertEqual((new_left.query_start, new_left.query_end), (1,5))
        self.assertEqual((new_right.query_start, new_right.query_end), (8,12))
        self.assertEqual(len(new_left.target_sequence or ""), 15)
        self.assertEqual(len(new_right.target_sequence or ""), 15)

    def test_adjust_target_coordinates_overlap_k0_no_change(self):
        left = NucMatch("q","t",1,5,100,114,0.0,100.0,False); left.target_sequence="A"*15
        right = NucMatch("q","t",5,9,200,214,0.0,100.0,False); right.target_sequence="C"*15
        cand = Candidate(split_k=0, window_seq="", stitched_local="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        self.assertEqual((new_left.query_start, new_left.query_end), (1,5))
        self.assertEqual((new_right.query_start, new_right.query_end), (5,9))
        self.assertEqual(len(new_left.target_sequence or ""), 15)
        self.assertEqual(len(new_right.target_sequence or ""), 15)

    def test_adjust_target_coordinates_overlap_k1_trims_and_adjacent(self):
        left = NucMatch("q","t",1,5,100,114,0.0,100.0,False); left.target_sequence="A"*15
        right = NucMatch("q","t",4,9,200,217,0.0,100.0,False); right.target_sequence="C"*18
        cand = Candidate(split_k=1, window_seq="", stitched_local="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        # left reduced by 1 AA at end: 1..4
        self.assertEqual((new_left.query_start, new_left.query_end), (1,4))
        # right unchanged
        self.assertEqual((new_right.query_start, new_right.query_end), (4,9))
        # DNA trimmed by 3 bases on left only
        self.assertEqual(len(new_left.target_sequence or ""), 12)
        self.assertEqual(len(new_right.target_sequence or ""), 18)

    def test_adjust_target_coordinates_overlap_trim_all_from_left(self):
        # left is length 1 AA, overlap 1, k=1 removes entire left
        left = NucMatch("q","t",1,1,100,102,0.0,100.0,False); left.target_sequence="A"*3
        right = NucMatch("q","t",1,4,200,211,0.0,100.0,False); right.target_sequence="C"*12
        cand = Candidate(split_k=1, window_seq="", stitched_local="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        # left becomes zero-length
        self.assertTrue(new_left.query_end < new_left.query_start)
        # right unchanged
        self.assertEqual(new_right.query_start, 1)
        self.assertEqual((new_right.query_end - new_right.query_start + 1), 4)
        self.assertEqual(len(new_right.target_sequence or ""), 12)


if __name__ == "__main__":
    unittest.main()


