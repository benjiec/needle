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
    compute_three_frame_translations,
    parse_hmmsearch_domtbl,
    hmmsearch_to_dna_coords,
    find_matches_at_locus
)
import needle.hits as hits_mod

from needle.blast import Results, group_matches, ProteinMatch, NucMatch


class TestHMMSearch(unittest.TestCase):

    def test_parse_hmmsearch_domtbl_returns_matches(self):

        with tempfile.TemporaryDirectory() as temp_dir:

            temp_file_path = os.path.join(temp_dir, 'domtbl.txt')
            with open(temp_file_path, 'w') as f:
                f.write("# target name  accession  tlen  query name  accession  qlen  E-value  score  bias  #  of  c-Evalue  i-Evalue  score  bias  from  to  from  to\n")
                f.write("cand_0 _ _ _ _ _ 0.1 100 _ _ _ _ _ _ _ 11 12 13 14\n")
                f.write("cand_0 _ _ _ _ _ 0.2 101 _ _ _ _ _ _ _ 21 22 23 24\n")
                f.close()

            matches = parse_hmmsearch_domtbl(temp_file_path)
            self.assertEqual(len(matches), 2)
            self.assertEqual(matches[0]["target_name"], "cand_0")
            self.assertEqual(matches[0]["evalue"], 0.1)
            self.assertEqual(matches[0]["score"], 100)
            self.assertEqual(matches[0]["hmm_from"], 11)
            self.assertEqual(matches[0]["hmm_to"], 12)
            self.assertEqual(matches[0]["target_from"], 13)
            self.assertEqual(matches[0]["target_to"], 14)
            self.assertEqual(matches[1]["target_name"], "cand_0")
            self.assertEqual(matches[1]["evalue"], 0.2)
            self.assertEqual(matches[1]["score"], 101)
            self.assertEqual(matches[1]["hmm_from"], 21)
            self.assertEqual(matches[1]["hmm_to"], 22)
            self.assertEqual(matches[1]["target_from"], 23)
            self.assertEqual(matches[1]["target_to"], 24)

    def test_parse_hmmsearch_domtbl_asserts_has_expected_headers(self):

        with tempfile.TemporaryDirectory() as temp_dir:

            temp_file_path = os.path.join(temp_dir, 'domtbl.txt')
            with open(temp_file_path, 'w') as f:
                # BAD unexpected header
                f.write("# target name  accession  tlen  query name  accession  qlen  E-value  score  from  to  from  to\n")
                f.write("cand_0 _ _ _ _ _ 0.1 100 11 12 13 14\n")
                f.write("cand_0 _ _ _ _ _ 0.2 101 21 22 23 24\n")
                f.close()

            with self.assertRaises(AssertionError):
                parse_hmmsearch_domtbl(temp_file_path)


class TestCleaningSequenceWithHMM(unittest.TestCase):

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


class TestRefiningHitsWithHMM(unittest.TestCase):

    def test_three_frame_translation_computes_per_frame_sequence_and_dna_coordinates(self):
        genomic = "ATGCGATGACTTCGTTATGCTT"

        # fwd
        self.assertEqual(
          compute_three_frame_translations(genomic, 2, 17),
          [
            (2, 16, "CDDFV"),
            (3, 17, "AMTSL"),
            (4, 15, "R*LR"),
          ]
        )

        # rev
        # TGCGATGACTTCGTT -> AACGAAGTCATCGCA
        self.assertEqual(
          compute_three_frame_translations(genomic, 16, 2),
          [
            (16, 2, "NEVIA"),
            (15, 4, "TKSS"),
            (14, 3, "RSHR")
          ]
        )

    def test_hmmsearch_to_dna_coords_converts_aa_coordinates_from_results_to_dna_coordinates_on_fwd_strand(self):

        orig = hits_mod.hmmsearch
        fake_matches = [
            dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=5, hmm_to=10, target_from=2, target_to=7),
            dict(target_name="cand_1", score=100, evalue=0.2, hmm_from=9, hmm_to=15, target_from=7, target_to=13),
            dict(target_name="cand_2", score=100, evalue=0.3, hmm_from=16, hmm_to=20, target_from=16, target_to=20)
        ]

        try:
            hits_mod.hmmsearch = lambda h,s: fake_matches

            translations = [
              (101, 199, "A"*33),
              (102, 197, "L"*32),
              (103, 198, "C"*32),
            ]

            matches = hmmsearch_to_dna_coords("hmmfile", translations)

            self.assertEqual(len(matches), 3)

            self.assertEqual(matches[0]["target_name"], "cand_0")
            self.assertEqual(matches[0]["hmm_from"], 5)
            self.assertEqual(matches[0]["hmm_to"], 10)
            # frame 0, target at 2 which is dna 4-6, so adjusted coordinate is 104
            self.assertEqual(matches[0]["target_from"], 104)
            # frame 0, target at 7 which is dna 19-21, so adjusted coordinate is 121
            self.assertEqual(matches[0]["target_to"], 121)
            self.assertEqual(matches[0]["evalue"], 0.1)
            self.assertEqual(matches[0]["matched_sequence"], "A"*6)

            self.assertEqual(matches[1]["target_name"], "cand_1")
            self.assertEqual(matches[1]["hmm_from"], 9)
            self.assertEqual(matches[1]["hmm_to"], 15)
            # frame 1, target at 7 which is dna 19-21, so adjusted coordinate from 102 is 120
            self.assertEqual(matches[1]["target_from"], 120)
            # frame 1, target at 13 which is dna 37-39, so adjusted coordinate from 102 is 140
            self.assertEqual(matches[1]["target_to"], 140)
            self.assertEqual(matches[1]["evalue"], 0.2)
            self.assertEqual(matches[1]["matched_sequence"], "L"*7)

            self.assertEqual(matches[2]["target_name"], "cand_2")
            self.assertEqual(matches[2]["hmm_from"], 16)
            self.assertEqual(matches[2]["hmm_to"], 20)
            # frame 2, target at 16 which is dna 46-48, so adjusted coordinate from 103 is 148
            self.assertEqual(matches[2]["target_from"], 148)
            # frame 2, target at 20 which is dna 58-60, so adjusted coordinate from 103 is 162
            self.assertEqual(matches[2]["target_to"], 162)
            self.assertEqual(matches[2]["evalue"], 0.3)
            self.assertEqual(matches[2]["matched_sequence"], "C"*5)

        finally:
            hits_mod.hmmsearch = orig

    def test_hmmsearch_to_dna_coords_converts_aa_coordinates_from_results_to_dna_coordinates_on_rev_strand(self):

        orig = hits_mod.hmmsearch
        fake_matches = [
            dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=5, hmm_to=10, target_from=2, target_to=7),
            dict(target_name="cand_1", score=100, evalue=0.2, hmm_from=9, hmm_to=15, target_from=7, target_to=13),
            dict(target_name="cand_2", score=100, evalue=0.3, hmm_from=16, hmm_to=20, target_from=16, target_to=20)
        ]

        try:
            hits_mod.hmmsearch = lambda h,s: fake_matches

            translations = [
              (199, 101, "A"*33),
              (198, 103, "L"*32),
              (197, 102, "C"*32),
            ]

            matches = hmmsearch_to_dna_coords("hmmfile", translations)

            self.assertEqual(len(matches), 3)

            self.assertEqual(matches[0]["target_name"], "cand_0")
            self.assertEqual(matches[0]["hmm_from"], 5)
            self.assertEqual(matches[0]["hmm_to"], 10)
            # frame 0, target at 2 which is dna 4-6, so adjusted coordinate from 199 decreasing is 196
            self.assertEqual(matches[0]["target_from"], 196)
            # frame 0, target at 7 which is dna 19-21, so adjusted coordinate from 199 decreasing is 179
            self.assertEqual(matches[0]["target_to"], 179)
            self.assertEqual(matches[0]["evalue"], 0.1)
            self.assertEqual(matches[0]["matched_sequence"], "A"*6)

            self.assertEqual(matches[1]["target_name"], "cand_1")
            self.assertEqual(matches[1]["hmm_from"], 9)
            self.assertEqual(matches[1]["hmm_to"], 15)
            # frame 1, target at 7 which is dna 19-21, so adjusted coordinate from 198 decreasing is 180
            self.assertEqual(matches[1]["target_from"], 180)
            # frame 1, target at 13 which is dna 37-39, so adjusted coordinate from 198 decreasing is 160
            self.assertEqual(matches[1]["target_to"], 160)
            self.assertEqual(matches[1]["evalue"], 0.2)
            self.assertEqual(matches[1]["matched_sequence"], "L"*7)

            self.assertEqual(matches[2]["target_name"], "cand_2")
            self.assertEqual(matches[2]["hmm_from"], 16)
            self.assertEqual(matches[2]["hmm_to"], 20)
            # frame 2, target at 16 which is dna 46-48, so adjusted coordinate from 197 decreasing is 152
            self.assertEqual(matches[2]["target_from"], 152)
            # frame 2, target at 20 which is dna 58-60, so adjusted coordinate from 197 decreasing is 138
            self.assertEqual(matches[2]["target_to"], 138)
            self.assertEqual(matches[2]["evalue"], 0.3)
            self.assertEqual(matches[2]["matched_sequence"], "C"*5)

        finally:
            hits_mod.hmmsearch = orig

    def test_hmmsearch_to_dna_coords_ignores_hmmsearch_result_not_on_fwd_strand(self):

        orig = hits_mod.hmmsearch
        fake_matches = [
            dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=5, hmm_to=10, target_from=2, target_to=5),
            dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=5, hmm_to=10, target_from=8, target_to=5),
            dict(target_name="cand_1", score=100, evalue=0.2, hmm_from=9, hmm_to=15, target_from=7, target_to=13),
            dict(target_name="cand_2", score=100, evalue=0.3, hmm_from=16, hmm_to=20, target_from=16, target_to=20)
        ]

        try:
            hits_mod.hmmsearch = lambda h,s: fake_matches

            translations = [
              (199, 101, "A"*33),
              (198, 103, "L"*32),
              (197, 102, "C"*32),
            ]

            matches = hmmsearch_to_dna_coords("hmmfile", translations)

            self.assertEqual(len(matches), 3)

            self.assertEqual(matches[0]["target_name"], "cand_0")
            self.assertEqual(matches[0]["hmm_from"], 5)
            self.assertEqual(matches[0]["hmm_to"], 10)
            self.assertEqual(matches[1]["target_name"], "cand_1")
            self.assertEqual(matches[1]["hmm_from"], 9)
            self.assertEqual(matches[1]["hmm_to"], 15)
            self.assertEqual(matches[2]["target_name"], "cand_2")
            self.assertEqual(matches[2]["hmm_from"], 16)
            self.assertEqual(matches[2]["hmm_to"], 20)

        finally:
            hits_mod.hmmsearch = orig

    def test_find_matches_at_locus_incrementally_search_for_more_matches(self):

        orig = hits_mod.hmmsearch_to_dna_coords

        searched = []
        def fake_hmmsearch_to_dna_coords(_, translations):
            first_match = [
                dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=5, hmm_to=10, target_from=10001, target_to=10018, matched_sequence="F"*6),
                dict(target_name="cand_1", score=100, evalue=0.2, hmm_from=9, hmm_to=15, target_from=11021, target_to=11041, matched_sequence="F"*7),
                dict(target_name="cand_2", score=100, evalue=0.3, hmm_from=16, hmm_to=20, target_from=11051, target_to=11065, matched_sequence="F"*5)
            ]
            second_match = [
                # add a new match
                dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=1, hmm_to=4, target_from=9001, target_to=9012, matched_sequence="F"*4),
                dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=5, hmm_to=10, target_from=10001, target_to=10018, matched_sequence="F"*6),
                dict(target_name="cand_1", score=100, evalue=0.2, hmm_from=9, hmm_to=15, target_from=11021, target_to=11041, matched_sequence="F"*7),
                dict(target_name="cand_2", score=100, evalue=0.3, hmm_from=16, hmm_to=20, target_from=11051, target_to=11065, matched_sequence="F"*5)
            ]

            searched.append((translations[0][0], translations[0][1]))

            if translations[0][0] == 10001: # initial
                return first_match
            elif translations[0][0] == 8001: # step is 2000, this is second search
                return second_match
            elif translations[0][0] == 6001: # step is 2000, this is third search, return same match
                return second_match
            raise Exception("Should not be searching anymore")

        try:
            hits_mod.hmmsearch_to_dna_coords = fake_hmmsearch_to_dna_coords
          
            old_matches = [
                NucMatch(query_accession="Q", target_accession="T", query_start=5, query_end=10, target_start=10001, target_end=10018, e_value=0.1, identity=None)
            ]

            new_matches = find_matches_at_locus(
                old_matches,
                "T"*15000,
                10001, 12001, "hmmfile", step=2000
            )

            self.assertNotEqual(new_matches, None)
            self.assertEqual(len(new_matches), 4)

            # does not search beyond sequence limit
            self.assertEqual(searched, [(10001, 12001), (8001, 14000), (6001, 15000)])

            self.assertEqual(new_matches[0].query_accession, "Q")
            self.assertEqual(new_matches[0].target_accession, "T")

            self.assertEqual(new_matches[0].query_start, 1)
            self.assertEqual(new_matches[0].query_end, 4)
            self.assertEqual(new_matches[0].target_start, 9001)
            self.assertEqual(new_matches[0].target_end, 9012)

            self.assertEqual(new_matches[1].query_start, 5)
            self.assertEqual(new_matches[1].query_end, 10)
            self.assertEqual(new_matches[1].target_start, 10001)
            self.assertEqual(new_matches[1].target_end, 10018)

            self.assertEqual(new_matches[2].query_start, 9)
            self.assertEqual(new_matches[2].query_end, 15)
            self.assertEqual(new_matches[2].target_start, 11021)
            self.assertEqual(new_matches[2].target_end, 11041)

            self.assertEqual(new_matches[3].query_start, 16)
            self.assertEqual(new_matches[3].query_end, 20)
            self.assertEqual(new_matches[3].target_start, 11051)
            self.assertEqual(new_matches[3].target_end, 11065)

        finally:
            hits_mod.hmmsearch_to_dna_coords = orig

    def test_find_matches_at_locus_stops_searching_if_found_nonlinear_overlap(self):

        orig = hits_mod.hmmsearch_to_dna_coords

        searched = []
        def fake_hmmsearch_to_dna_coords(_, translations):
            first_match = [
                dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=5, hmm_to=10, target_from=10018, target_to=10001, matched_sequence="F"*6),
                dict(target_name="cand_1", score=100, evalue=0.2, hmm_from=9, hmm_to=15, target_from=11041, target_to=11021, matched_sequence="F"*7),
                dict(target_name="cand_2", score=100, evalue=0.3, hmm_from=16, hmm_to=20, target_from=11065, target_to=11051, matched_sequence="F"*5)
            ]

            second_match = [
                dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=1, hmm_to=11, target_from=8550, target_to=8518, matched_sequence="F"*11),
                dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=5, hmm_to=10, target_from=10018, target_to=10001, matched_sequence="F"*6),
                dict(target_name="cand_1", score=100, evalue=0.2, hmm_from=9, hmm_to=15, target_from=11041, target_to=11021, matched_sequence="F"*7),
                dict(target_name="cand_2", score=100, evalue=0.3, hmm_from=16, hmm_to=20, target_from=11065, target_to=11051, matched_sequence="F"*5)
            ]

            searched.append((translations[0][0], translations[0][1]))

            if translations[0][0] == 12001: # initial
                return first_match
            elif translations[0][0] == 14001: # step is 2000, this is second search
                return first_match
            raise Exception("Should not be searching anymore")

        try:
            hits_mod.hmmsearch_to_dna_coords = fake_hmmsearch_to_dna_coords
          
            old_matches = [
                NucMatch(query_accession="Q", target_accession="T", query_start=5, query_end=10, target_start=10018, target_end=10001, e_value=0.1, identity=None)
            ]

            new_matches = find_matches_at_locus(
                old_matches,
                "A"*20000,
                12001, 10001, "hmmfile", step=2000
            )

            self.assertNotEqual(new_matches, None)
            self.assertEqual(len(new_matches), 3)
            self.assertEqual(searched, [(12001, 10001), (14001, 8002)])

        finally:
            hits_mod.hmmsearch_to_dna_coords = orig

    def test_find_matches_at_locus_incrementally_search_on_rev_strand_as_well(self):

        orig = hits_mod.hmmsearch_to_dna_coords

        searched = []
        def fake_hmmsearch_to_dna_coords(_, translations):
            first_match = [
                dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=5, hmm_to=10, target_from=10018, target_to=10001, matched_sequence="F"*6),
                dict(target_name="cand_1", score=100, evalue=0.2, hmm_from=9, hmm_to=15, target_from=11041, target_to=11021, matched_sequence="F"*7),
                dict(target_name="cand_2", score=100, evalue=0.3, hmm_from=16, hmm_to=20, target_from=11065, target_to=11051, matched_sequence="F"*5)
            ]

            searched.append((translations[0][0], translations[0][1]))

            if translations[0][0] == 12001: # initial
                return first_match
            elif translations[0][0] == 14001: # step is 2000, this is second search
                return first_match
            raise Exception("Should not be searching anymore")

        try:
            hits_mod.hmmsearch_to_dna_coords = fake_hmmsearch_to_dna_coords
          
            old_matches = [
                NucMatch(query_accession="Q", target_accession="T", query_start=5, query_end=10, target_start=10018, target_end=10001, e_value=0.1, identity=None)
            ]

            new_matches = find_matches_at_locus(
                old_matches,
                "A"*15000,
                12001, 10001, "hmmfile", step=2000
            )

            self.assertNotEqual(new_matches, None)
            self.assertEqual(len(new_matches), 3)
            self.assertEqual(searched, [(12001, 10001), (14001, 8002)])

        finally:
            hits_mod.hmmsearch_to_dna_coords = orig

    def test_find_matches_at_locus_ensures_hmm_matched_sequence_matches_translated_sequence(self):

        orig = hits_mod.hmmsearch_to_dna_coords
        expected_aa = "F"*6

        def fake_hmmsearch_to_dna_coords(_, translations):
            first_match = [
                dict(target_name="cand_0", score=100, evalue=0.1, hmm_from=5, hmm_to=10, target_from=10001, target_to=10018, matched_sequence=expected_aa)
            ]
            return first_match

        try:
            hits_mod.hmmsearch_to_dna_coords = fake_hmmsearch_to_dna_coords
          
            old_matches = [
                NucMatch(query_accession="Q", target_accession="T", query_start=5, query_end=10, target_start=10001, target_end=10018, e_value=0.1, identity=None)
            ]

            # target sequence on DNA matches what HMM says
            new_matches = find_matches_at_locus(
                old_matches,
                "C"*10000+"T"*18+"C"*20000,
                10001, 20000, "hmmfile", step=2000, force_extend=True
            )
            self.assertEqual(len(new_matches), 1)
            self.assertEqual(new_matches[0].target_sequence_translated(), expected_aa)

            # target sequence does not match what HMM says, for some odd reason
            with self.assertRaises(AssertionError):
                new_matches = find_matches_at_locus(
                    old_matches,
                    # not 18 T's
                    "C"*10000+"T"*14+"C"*20000,
                    10001, 20000, "hmmfile", step=2000, force_extend=True
                )

        finally:
            hits_mod.hmmsearch_to_dna_coords = orig


if __name__ == "__main__":
    unittest.main()
