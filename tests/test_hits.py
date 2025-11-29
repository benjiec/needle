import os
import tempfile
import unittest
import io

from needle.blast import order_matches_for_junctions

from needle.hits import (
    generate_transition_candidates,
    score_and_select_best_transition,
    stitch_cleaned_sequence,
    Candidate,
    hmm_clean_protein,
    hmm_clean,
    adjust_target_coordinates,
    write_protein_row,
    write_nucmatch_rows,
    write_fasta_record,
    export_protein_hits,
    get_hmmsearch_score_eval,
)
import needle.hits as hits_mod

from needle.blast import Results, group_matches, ProteinMatch, NucMatch


class TestHits(unittest.TestCase):
    def test_generate_transition_candidates_overlap(self):
        left = "ABCDEFX"
        right = "yefghij"

        cands = generate_transition_candidates(left, right, overlap_len=4, gap_len=0, overlap_flanking_len=2)
        self.assertEqual(len(cands), 5)

        self.assertEqual(
            [c.left_trimmed for c in cands],
            [
                4,     # k=0: left trims 4, right trims 0
                3,     # k=1: left trims 3, right trims 1
                2,     # k=2: left trims 2, right trims 2
                1,     # k=3: left trims 1, right trims 3
                0,     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.stitched for c in cands],
            [
                "ABCyefghij",     # k=0: left trims 4, right trims 0
                "ABCDefghij",     # k=1: left trims 3, right trims 1
                "ABCDEfghij",     # k=2: left trims 2, right trims 2
                "ABCDEFghij",     # k=3: left trims 1, right trims 3
                "ABCDEFXhij",     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.right_kept for c in cands],
            [
                "yefghij",     # k=0: left trims 4, right trims 0
                "efghij",     # k=1: left trims 3, right trims 1
                "fghij",     # k=2: left trims 2, right trims 2
                "ghij",     # k=3: left trims 1, right trims 3
                "hij",     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.window_seq for c in cands],
            [
                "BCyefghi",     # k=0: left trims 4, right trims 0
                "BCDefghi",     # k=1: left trims 3, right trims 1
                "BCDEfghi",     # k=2: left trims 2, right trims 2
                "BCDEFghi",     # k=3: left trims 1, right trims 3
                "BCDEFXhi",     # k=4: left trims 0, right trims 4
            ],
        )

    def test_generate_transition_candidates_overlap_same_as_full(self):
        left = "DEFX"
        right = "yefg"

        cands = generate_transition_candidates(left, right, overlap_len=4, gap_len=0, overlap_flanking_len=2)
        self.assertEqual(len(cands), 5)

        self.assertEqual(
            [c.left_trimmed for c in cands],
            [
                4,     # k=0: left trims 4, right trims 0
                3,     # k=1: left trims 3, right trims 1
                2,     # k=2: left trims 2, right trims 2
                1,     # k=3: left trims 1, right trims 3
                0,     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.stitched for c in cands],
            [
                "yefg",     # k=0: left trims 4, right trims 0
                "Defg",     # k=1: left trims 3, right trims 1
                "DEfg",     # k=2: left trims 2, right trims 2
                "DEFg",     # k=3: left trims 1, right trims 3
                "DEFX",     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.right_kept for c in cands],
            [
                "yefg",     # k=0: left trims 4, right trims 0
                "efg",     # k=1: left trims 3, right trims 1
                "fg",     # k=2: left trims 2, right trims 2
                "g",     # k=3: left trims 1, right trims 3
                "",     # k=4: left trims 0, right trims 4
            ],
        )
        self.assertEqual(
            [c.window_seq for c in cands],
            [
                "yefg",     # k=0: left trims 4, right trims 0
                "Defg",     # k=1: left trims 3, right trims 1
                "DEfg",     # k=2: left trims 2, right trims 2
                "DEFg",     # k=3: left trims 1, right trims 3
                "DEFX",     # k=4: left trims 0, right trims 4
            ],
        )
  
    def test_generate_transition_candidates_gap(self):
        cands_gap = generate_transition_candidates("AAA", "bbb", overlap_len=0, gap_len=2, overlap_flanking_len=2)
        self.assertEqual(len(cands_gap), 1)
        self.assertEqual(cands_gap[0].assigned_overlap_to_left, None)
        self.assertEqual(cands_gap[0].left_trimmed, 0)
        self.assertEqual(cands_gap[0].stitched, "AAAXXbbb")
        self.assertEqual(cands_gap[0].right_kept, "XXbbb")
        self.assertEqual(cands_gap[0].window_seq, "AAXXbb")

    def test_stitch_cleaned_sequence_basic(self):
        class _MM: pass
        left = _MM(); right = _MM()
        left.query_start=1; left.query_end=5
        right.query_start=4; right.query_end=8
        aa_map = {id(left):"ABCDE", id(right):"DEFGH"}
        pairs = order_matches_for_junctions([left, right])  # type: ignore
        self.assertEqual(pairs[0][2], 2)
        cand = Candidate(assigned_overlap_to_left=1, window_seq="", stitched="ABCDEFGH", left_trimmed=2, right_kept="DEFGH")
        stitched = stitch_cleaned_sequence([(left, right, None, None)], {0: cand}, aa_map)  # type: ignore
        self.assertEqual(stitched, "ABCDEFGH")

    def test_stitch_cleaned_sequence_multiple_blocks_mixed(self):
        class _MM: pass
        a=_MM(); b=_MM(); c=_MM()
        a.query_start=1; a.query_end=5
        b.query_start=4; b.query_end=9
        c.query_start=13; c.query_end=15
        aa_map = {id(a):"ABCDE", id(b):"DEFGHI", id(c):"KLM"}
        pairs = order_matches_for_junctions([a,b,c])  # type: ignore
        cand0 = Candidate(assigned_overlap_to_left=1, window_seq="", stitched="ABCDEFGHI", left_trimmed=1, right_kept="EFGHI")
        cand1 = Candidate(assigned_overlap_to_left=None, window_seq="", stitched="DEFGHIXXXKLM", left_trimmed=0, right_kept="XXXKLM")
        stitched = stitch_cleaned_sequence(
          [(a,b, None, None), (b,c, None, None)],
          {0:cand0, 1:cand1}, aa_map)
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
                    parts = [name, "x1", "x2", "x3", "x4", "x5", "3.2", str(score), "x8"]
                    return " ".join(parts) + "\n"
                f.write("# target name        accession   tlen query name           accession   qlen   E-value  score\n")
                f.write(line("cand_0", 10.0))
                f.write(line("cand_1", 50.0))
            class _P: returncode = 0
            return _P()

        orig_td = _tempfile.TemporaryDirectory; orig_run = _subprocess.run
        try:
            _tempfile.TemporaryDirectory = _FakeTD  # type: ignore
            _subprocess.run = _fake_run  # type: ignore
            c1 = Candidate(assigned_overlap_to_left=0, window_seq="AAAA", stitched="LEFT", left_trimmed=None, right_kept="")
            c2 = Candidate(assigned_overlap_to_left=1, window_seq="BBBB", stitched="RIGHT", left_trimmed=None, right_kept="")
            best = score_and_select_best_transition([c1, c2], hmm_file_name="ignored.hmm")
            self.assertIs(best, c2)
        finally:
            _tempfile.TemporaryDirectory = orig_td  # type: ignore
            _subprocess.run = orig_run  # type: ignore

    def test_hmm_cleaned_protein_integration_with_mock_scoring(self):
        a = NucMatch("Q","T",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"    # MEF
        b = NucMatch("Q","T",3,5,10,18,0.0,100.0,False); b.target_sequence="GAAGTGGGG"  # EVG
        c = NucMatch("Q","T",9,9,30,32,0.0,100.0,False); c.target_sequence="ATG"        # M
        pm = ProteinMatch("T",[a,b,c],1,9,1,32)

        orig = hits_mod.score_and_select_best_transition
        def _fake(cands, hmm): 
            for x in cands:
                if x.assigned_overlap_to_left == 1: return x
            return cands[0]
        try:
            hits_mod.score_and_select_best_transition = _fake
            cleaned_pm = hmm_clean_protein(pm, "dummy.hmm", overlap_flanking_len=5)
            cleaned = cleaned_pm.collated_protein_sequence
        finally:
            hits_mod.score_and_select_best_transition = orig

        self.assertEqual(cleaned, "MEFVGXXXM")

    def test_hmm_clean_protein_adjusts_overlap_coordinates(self):
        # Overlap: a(1..5), b(4..9) => overlap 2; choose k=1; c(13..15) should shift by 1

        a = NucMatch("Q","T",1,5,1,15,0.0,100.0,False); a.target_sequence="ATG"*5         # 'M'*5
        b = NucMatch("Q","T",4,9,16,33,0.0,100.0,False); b.target_sequence="GAA"*6        # 'E'*6
        c = NucMatch("Q","T",13,15,40,48,0.0,100.0,False); c.target_sequence="ATG"*3      # 'M'*3
        pm = ProteinMatch("T",[a,b,c],1,15,1,48)

        orig = hits_mod.score_and_select_best_transition
        def _fake(cands, hmm):
            for x in cands:
                if x.assigned_overlap_to_left == 1:
                    return x
            return cands[0]
        try:
            hits_mod.score_and_select_best_transition = _fake  # type: ignore
            cleaned_pm = hmm_clean_protein(pm, "dummy.hmm", overlap_flanking_len=5)
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

    def test_adjust_target_coordinates_gap_keeps_blocks(self):
        # query acc, target acc, query start, query end, target start, target end
        left = NucMatch("q","t",1,5,100,114,0.0,100.0,False); left.target_sequence="A"*15
        right = NucMatch("q","t",8,12,200,214,0.0,100.0,False); right.target_sequence="C"*15
        cand = Candidate(assigned_overlap_to_left=None, window_seq="", stitched="", left_trimmed=None, right_kept="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        self.assertEqual((new_left.query_start, new_left.query_end), (1,5))
        self.assertEqual((new_right.query_start, new_right.query_end), (8,12))
        self.assertEqual(new_left.target_sequence, "A"*15)
        self.assertEqual(new_right.target_sequence, "C"*15)

    def test_adjust_target_coordinates_overlap_k0_no_change(self):
        # query acc, target acc, query start, query end, target start, target end
        left = NucMatch("q","t",1,5,100,114,0.0,100.0,False); left.target_sequence="A"*15
        right = NucMatch("q","t",5,9,200,214,0.0,100.0,False); right.target_sequence="C"*15
        cand = Candidate(assigned_overlap_to_left=0, window_seq="", stitched="", left_trimmed=1, right_kept="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        self.assertEqual((new_left.query_start, new_left.query_end), (1,4))
        self.assertEqual((new_right.query_start, new_right.query_end), (5,9))
        self.assertEqual(new_left.target_sequence, "A"*12)
        self.assertEqual(new_right.target_sequence, "C"*15)

    def test_adjust_target_coordinates_overlap_k1_trims_and_adjacent(self):
        # query acc, target acc, query start, query end, target start, target end
        left = NucMatch("q","t",1,5,100,114,0.0,100.0,False); left.target_sequence="A"*15
        right = NucMatch("q","t",4,9,200,217,0.0,100.0,False); right.target_sequence="C"*18
        cand = Candidate(assigned_overlap_to_left=1, window_seq="", stitched="", left_trimmed=1, right_kept="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        self.assertEqual((new_left.query_start, new_left.query_end), (1,4))
        self.assertEqual((new_right.query_start, new_right.query_end), (5,9))
        self.assertEqual(new_left.target_sequence, "A"*12)
        self.assertEqual(new_right.target_sequence, "C"*15)

    def test_adjust_target_coordinates_overlap_trim_all_from_left(self):
        # query acc, target acc, query start, query end, target start, target end
        left = NucMatch("q","t",1,1,100,102,0.0,100.0,False); left.target_sequence="A"*3
        right = NucMatch("q","t",1,4,200,211,0.0,100.0,False); right.target_sequence="C"*12
        cand = Candidate(assigned_overlap_to_left=0, window_seq="", stitched="", left_trimmed=1, right_kept="")
        new_left, new_right = adjust_target_coordinates(left, right, cand)
        self.assertEqual((new_left.query_start, new_left.query_end), (1,0))
        self.assertEqual((new_right.query_start, new_right.query_end), (1,4))
        self.assertEqual(new_left.target_sequence, "")
        self.assertEqual(new_right.target_sequence, "C"*12)

    def test_write_protein_row(self):
        a = NucMatch("QX","TX",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        pm = ProteinMatch("TX",[a],1,3,1,9)
        pid = pm.protein_hit_id
        buf = io.StringIO()
        write_protein_row(buf, "GENOME1", pm, 1e-5, 42.0)
        line = buf.getvalue().strip()
        cols = line.split("\t")
        self.assertEqual(cols[0], pid)
        self.assertEqual(cols[1], "QX")
        self.assertEqual(cols[2], "TX")
        self.assertEqual(cols[3], "GENOME1")
        self.assertEqual(float(cols[4]), 1e-5)
        self.assertEqual(float(cols[5]), 42.0)
        self.assertTrue(cols[6].startswith("MEF"))

    def test_write_nucmatch_rows(self):
        a = NucMatch("QX","TX",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        b = NucMatch("QX","TX",4,6,10,18,0.0,90.0,False); b.target_sequence="GAAGTGGGG"
        pm = ProteinMatch("TX",[a,b],1,6,1,18)
        pid = pm.protein_hit_id
        buf = io.StringIO()
        write_nucmatch_rows(buf, pm)
        rows = [r for r in buf.getvalue().splitlines() if r]
        self.assertEqual(len(rows), 2)
        cols0 = rows[0].split("\t")
        self.assertEqual(cols0[0], pid)
        self.assertEqual(cols0[1], "TX")
        self.assertEqual(cols0[2], "1")
        self.assertEqual(cols0[3], "9")
        self.assertEqual(cols0[4], "1")
        self.assertEqual(cols0[5], "3")
        cols1 = rows[1].split("\t")
        self.assertEqual(cols1[0], pid)
        self.assertEqual(cols1[2], "10")
        self.assertEqual(cols1[3], "18")
        self.assertEqual(cols1[4], "4")
        self.assertEqual(cols1[5], "6")

    def test_write_fasta_record(self):
        a = NucMatch("QX","TX",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        pm = ProteinMatch("TX",[a],1,3,1,9)
        pid = pm.protein_hit_id
        buf = io.StringIO()
        write_fasta_record(buf, pm)
        s = buf.getvalue().splitlines()
        self.assertEqual(s[0], f">{pid}")
        self.assertTrue(s[1].startswith("MEF"))

    def test_export_protein_hits_filters_and_writes(self):
        # pm1: single block => eligible
        a = NucMatch("Q1","T1",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        pm1 = ProteinMatch("T1",[a],1,3,1,9, hmm_file="ignored.hmm")
        # pm2: overlapping blocks => not single sequence
        b1 = NucMatch("Q2","T2",1,3,1,9,0.0,100.0,False); b1.target_sequence="ATGGAATTT"
        b2 = NucMatch("Q2","T2",3,5,10,18,0.0,100.0,False); b2.target_sequence="GAAGTGGGG"
        pm2 = ProteinMatch("T2",[b1,b2],1,5,1,18, hmm_file="ignored.hmm")
        with tempfile.TemporaryDirectory() as d:
            p1 = os.path.join(d, "prot.tsv")
            p2 = os.path.join(d, "nuc.tsv")
            p3 = os.path.join(d, "prot_fastas")
            # mock hmm scoring
            orig = hits_mod.get_hmmsearch_score_eval
            try:
                hits_mod.get_hmmsearch_score_eval = lambda hmm, seq: (55.0, 1e-6)
                export_protein_hits("GENOMEZ", [pm1, pm2], p1, p2, p3)
            finally:
                hits_mod.get_hmmsearch_score_eval = orig
            with open(p1) as f:
                lines = [l.strip() for l in f if l.strip()]
            with open(p2) as f:
                lines2 = [l.strip() for l in f if l.strip()]
            faa_path = os.path.join(p3, "Q1.faa")
            with open(faa_path) as f:
                lines3 = [l.strip() for l in f if l.strip()]
            # headers + 1 line for pm1 only
            self.assertEqual(len(lines), 2)
            self.assertIn("GENOMEZ", lines[1])
            self.assertIn("1e-06", lines[1])
            # nuc rows: only pm1 yields rows => 1 line + header
            self.assertEqual(len(lines2), 2)
            # fasta: 2 lines (header+seq)
            self.assertEqual(len(lines3), 2)

    def test_export_protein_hits_appends_not_overwrites(self):
        a = NucMatch("QX","TX",1,3,1,9,0.0,100.0,False); a.target_sequence="ATGGAATTT"
        pm = ProteinMatch("TX",[a],1,3,1,9, hmm_file="ignored.hmm")
        with tempfile.TemporaryDirectory() as d:
            p1 = os.path.join(d, "prot.tsv")
            p2 = os.path.join(d, "nuc.tsv")
            p3 = os.path.join(d, "fastas")
            orig = hits_mod.get_hmmsearch_score_eval
            try:
                hits_mod.get_hmmsearch_score_eval = lambda hmm, seq: (10.0, 2e-5)
                export_protein_hits("GENOME1", [pm], p1, p2, p3)
                # capture initial counts
                with open(p1) as f:
                    lines_prot_1 = [l for l in f if l.strip()]
                with open(p2) as f:
                    lines_nuc_1 = [l for l in f if l.strip()]
                faa_path = os.path.join(p3, "QX.faa")
                with open(faa_path) as f:
                    lines_faa_1 = [l for l in f if l.strip()]
                # second export appends
                export_protein_hits("GENOME1", [pm], p1, p2, p3)
                with open(p1) as f:
                    lines_prot_2 = [l for l in f if l.strip()]
                with open(p2) as f:
                    lines_nuc_2 = [l for l in f if l.strip()]
                with open(faa_path) as f:
                    lines_faa_2 = [l for l in f if l.strip()]
            finally:
                hits_mod.get_hmmsearch_score_eval = orig
            # proteins.tsv: header + 1 row initially; after append expect +1 row
            self.assertEqual(len(lines_prot_1), 2)
            self.assertEqual(len(lines_prot_2), 3)
            # nuc.tsv: header + 1 row initially; after append expect +1 row
            self.assertEqual(len(lines_nuc_1), 2)
            self.assertEqual(len(lines_nuc_2), 3)
            # faa file: 2 lines initially; after append expect +2 lines
            self.assertEqual(len(lines_faa_1), 2)
            self.assertEqual(len(lines_faa_2), 4)

if __name__ == "__main__":
    unittest.main()
