import os
import tempfile
import unittest
import random

from needle.match import group_matches, ProteinHit, Match, order_matches_for_junctions, NonlinearMatchException
from needle.match import extract_subsequence, extract_subsequence_strand_sensitive


class TestOrderGroupMatches(unittest.TestCase):

    @staticmethod
    def makeM(query_start, query_end, target_start, target_end, query_accession=None, target_accession=None):
        return Match(
            query_accession=query_accession, target_accession=target_accession,
            e_value=0, identity=None,
            query_start=query_start, query_end=query_end, target_start=target_start, target_end=target_end
        )

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

        m1 = self.makeM( 5, 10, 1001, 1020, "Q", "S1")
        m2 = self.makeM(15, 20, 1101, 1120, "Q", "S1")
        m3 = self.makeM(25, 30, 1201, 1220, "Q", "S1")
        m4 = self.makeM(35, 40, 1301, 1320, "Q", "S2")
        m5 = self.makeM(45, 50, 1401, 1420, "Q", "S2")
        m6 = self.makeM(55, 60, 1501, 1520, "Q", "S2")
        m7 = self.makeM(65, 70, 11601, 11620, "Q", "S2")
        m8 = self.makeM(70, 80, 11701, 11720, "Q", "S2")
        m9 = self.makeM(75, 90, 11801, 11820, "Q", "S2")

        matches = [m1, m2, m3, m4, m5, m6, m7, m8, m9]
        random.shuffle(matches)
        pms = group_matches(matches)
        pms = sorted(pms, key=lambda p: (p.target_accession, p.query_start))

        self.assertEqual(len(pms), 3)
        self.assertEqual(pms[0].target_accession, "S1")
        self.assertEqual(pms[0].query_start, 5)
        self.assertEqual(pms[0].query_end, 30)
        self.assertEqual(pms[0].target_start, 1001)
        self.assertEqual(pms[0].target_end, 1220)
        self.assertEqual(pms[1].target_accession, "S2")
        self.assertEqual(pms[1].query_start, 35)
        self.assertEqual(pms[1].query_end, 60)
        self.assertEqual(pms[1].target_start, 1301)
        self.assertEqual(pms[1].target_end, 1520)
        self.assertEqual(pms[2].target_accession, "S2")
        self.assertEqual(pms[2].query_start, 65)
        self.assertEqual(pms[2].query_end, 90)
        self.assertEqual(pms[2].target_start, 11601)
        self.assertEqual(pms[2].target_end, 11820)

    def test_group_matches_separate_matches_on_reverse_strand(self):

        m1 = self.makeM(25, 30, 1020, 1001, "Q", "S1")
        m2 = self.makeM(15, 20, 1120, 1101, "Q", "S1")
        m3 = self.makeM( 5, 10, 1220, 1201, "Q", "S1")
        m4 = self.makeM(35, 40, 11301, 11320, "Q", "S1")
        m5 = self.makeM(45, 50, 11401, 11420, "Q", "S1")
        m6 = self.makeM(55, 60, 11501, 11520, "Q", "S1")

        matches = [m1, m2, m3, m4, m5, m6]
        random.shuffle(matches)
	# reduce max_overlap_len to really short, so that 25->30 then 15->20
	# would look like an overlap of the same protein if we didn't handle
	# reverse strand
        pms = group_matches(matches, max_overlap_len=4)
        pms = sorted(pms, key=lambda p: (p.target_accession, p.query_start))

        self.assertEqual(len(pms), 2)
        self.assertEqual(pms[0].target_accession, "S1")
        self.assertEqual(pms[0].query_start, 5)
        self.assertEqual(pms[0].query_end, 30)
        self.assertEqual(pms[0].target_start, 1220)
        self.assertEqual(pms[0].target_end, 1001)
        self.assertEqual(pms[1].target_accession, "S1")
        self.assertEqual(pms[1].query_start, 35)
        self.assertEqual(pms[1].query_end, 60)
        self.assertEqual(pms[1].target_start, 11301)
        self.assertEqual(pms[1].target_end, 11520)

    def test_group_matches_separate_matches_on_different_strands(self):
        m1 = self.makeM( 5, 18, 1001, 1020, "Q", "S1")
        m2 = self.makeM(15, 28, 1101, 1120, "Q", "S1")
        m3 = self.makeM(25, 38, 1201, 1220, "Q", "S1")

        matches = [m1, m2, m3]
        random.shuffle(matches)

        pms = group_matches(matches)
        pms = sorted(pms, key=lambda p: (p.target_accession, p.query_start))
        self.assertEqual(len(pms), 1)

        m3 = self.makeM(25, 38, 1220, 1201, "Q", "S1")
        matches = [m1, m2, m3]
        random.shuffle(matches)

        pms = group_matches(matches)
        pms = sorted(pms, key=lambda p: (p.target_accession, p.query_start))
        self.assertEqual(len(pms), 2)

        self.assertEqual(pms[0].target_accession, "S1")
        self.assertEqual(pms[0].query_start, 5)
        self.assertEqual(pms[0].query_end, 28)
        self.assertEqual(pms[0].target_start, 1001)
        self.assertEqual(pms[0].target_end, 1120)
        self.assertEqual(pms[1].target_accession, "S1")
        self.assertEqual(pms[1].query_start, 25)
        self.assertEqual(pms[1].query_end, 38)
        self.assertEqual(pms[1].target_start, 1220)
        self.assertEqual(pms[1].target_end, 1201)

    def test_group_matches_separate_closely_placed_tandems(self):

        m1 = self.makeM( 5, 18, 1001, 1020, "Q", "S1")
        m2 = self.makeM(15, 28, 1101, 1120, "Q", "S1")
        m3 = self.makeM(25, 38, 1201, 1220, "Q", "S1")
        m4 = self.makeM(34, 48, 1301, 1520, "Q", "S1")

        matches = [m1, m2, m3, m4]
        random.shuffle(matches)

	# use max_overlap_len to control what is a tandem vs normal overlapping
	# matches to query, here we are saying if overlap is more than
	# specified, it's not an overlap. so the first 3 have overlap of 4
	# between them, but the next one has larger overlap
        pms = group_matches(matches, max_overlap_len=4)
        pms = sorted(pms, key=lambda p: (p.target_accession, p.query_start))

        self.assertEqual(len(pms), 2)
        self.assertEqual(pms[0].target_accession, "S1")
        self.assertEqual(pms[0].query_start, 5)
        self.assertEqual(pms[0].query_end, 38)
        self.assertEqual(pms[0].target_start, 1001)
        self.assertEqual(pms[0].target_end, 1220)
        self.assertEqual(pms[1].target_accession, "S1")
        self.assertEqual(pms[1].query_start, 34)
        self.assertEqual(pms[1].query_end, 48)
        self.assertEqual(pms[1].target_start, 1301)
        self.assertEqual(pms[1].target_end, 1520)

    def test_group_matches_separate_closely_placed_tandems_on_reverse_strand(self):

        m1 = self.makeM(34, 48, 1101, 1020, "Q", "S1")
        m2 = self.makeM(25, 38, 1201, 1120, "Q", "S1")
        m3 = self.makeM(15, 28, 1301, 1220, "Q", "S1")
        m4 = self.makeM( 5, 18, 1601, 1520, "Q", "S1")

        matches = [m1, m2, m3, m4]
        random.shuffle(matches)

	# use max_overlap_len to control what is a tandem vs normal overlapping
	# matches to query, here we are saying if overlap is more than
	# specified, it's not an overlap. so the first 3 have overlap of 4
	# between them, but the next one has larger overlap
        pms = group_matches(matches, max_overlap_len=4)
        pms = sorted(pms, key=lambda p: (p.target_accession, p.query_start))

        self.assertEqual(len(pms), 2)
        self.assertEqual(pms[0].target_accession, "S1")
        self.assertEqual(pms[0].query_start, 5)
        self.assertEqual(pms[0].query_end, 38)
        self.assertEqual(pms[0].target_start, 1601)
        self.assertEqual(pms[0].target_end, 1120)
        self.assertEqual(pms[1].target_accession, "S1")
        self.assertEqual(pms[1].query_start, 34)
        self.assertEqual(pms[1].query_end, 48)
        self.assertEqual(pms[1].target_start, 1101)
        self.assertEqual(pms[1].target_end, 1020)

    def test_protein_hit_id_deterministic_and_changes_when_inputs_change(self):
        # Two identical ProteinHit objects should have the same ID
        a1 = Match("QID","TID",1,3,1,9,0.0,100.0)
        a2 = Match("QID","TID",4,6,10,18,1e-12,95.0)
        pm1 = ProteinHit([a1, a2], 1, 6, 1, 18)
        pm2 = ProteinHit([a1, a2], 1, 6, 1, 18)
        self.assertEqual(pm1.protein_hit_id, pm2.protein_hit_id)
        # Change an input value (e.g., query_start) to produce a different ID
        b2 = Match("QID","TID",5,6,10,18,1e-12,95.0)
        pm3 = ProteinHit([a1, b2], 1, 6, 1, 18)
        self.assertNotEqual(pm1.protein_hit_id, pm3.protein_hit_id)

    def test_extra_match_in_middle_changes_protein_id(self):
        # Two identical ProteinHit objects should have the same ID
        a1 = Match("QID","TID",1,3,1,9,0.0,100.0)
        a2 = Match("QID","TID",4,6,10,18,1e-12,95.0)
        pm1 = ProteinHit([a1, a2], 1, 6, 1, 18)
        pm2 = ProteinHit([a1, a2], 1, 6, 1, 18)
        self.assertEqual(pm1.protein_hit_id, pm2.protein_hit_id)

        # Add new Match in middle, changes protein id
        b2 = Match("QID","TID",3,4,9,10,1e-12,95.0)
        pm3 = ProteinHit([a1, b2, a2], 1, 6, 1, 18)
        self.assertNotEqual(pm1.protein_hit_id, pm3.protein_hit_id)

    def test_target_sequence_and_target_sequence_translated_on_fwd_strand(self):

        # fwd
        m1 = self.makeM(34, 48, 1001, 1020, "Q", "S1")
        m1.target_sequence = "ATGATGATG"
        self.assertEqual(m1.target_sequence_translated(), "MMM")

        # rev, but use target_sequence which is already reversed
        m2 = self.makeM(25, 38, 1201, 1120, "Q", "S1")
        m2.target_sequence = "ATGATGATG"
        self.assertEqual(m2.target_sequence_translated(), "MMM")

    def test_translation_and_collated_protein_sequence(self):
        m1 = self.makeM(1, 3, 11, 19, "Q", "S1")
        m1.target_sequence = "ATGGAATTT"
        m2 = self.makeM(2, 4, 22, 30, "Q", "S1")
        m2.target_sequence = "GAAGTGGGG"

        self.assertEqual(m1.target_sequence_translated(), "MEF")
        self.assertEqual(m2.target_sequence_translated(), "EVG")

        pm = ProteinHit([m1, m2], 1, 4, 11, 30)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "M(EF/EV)G")

    def test_collated_protein_sequence_does_not_include_leading_Xs(self):
        m1 = self.makeM(2, 4, 11, 19, "Q", "S1")
        m1.target_sequence = "ATGGAATTT"
        m2 = self.makeM(5, 7, 22, 30, "Q", "S1")
        m2.target_sequence = "GAAGTGGGG"

        self.assertEqual(m1.target_sequence_translated(), "MEF")
        self.assertEqual(m2.target_sequence_translated(), "EVG")

        pm = ProteinHit([m1, m2], 2, 7, 11, 30)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "MEFEVG")

    def test_collate_handles_single_match(self):
        a = Match("Q","T",1,3,1,9,0.0,100.0); a.target_sequence="ATGGAATTT"    # MEF
        pm = ProteinHit([a],1,3,1,9)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "MEF")

    def test_collate_handles_gaps_and_overlaps(self):
        a = Match("Q","T",1,3,1,9,0.0,100.0); a.target_sequence="ATGGAATTT"    # MEF
        b = Match("Q","T",3,5,10,18,0.0,100.0); b.target_sequence="GAAGTGGGG"  # EVG
        c = Match("Q","T",9,9,30,32,0.0,100.0); c.target_sequence="ATG"        # M
        pm = ProteinHit([a,b,c],1,9,1,32)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "ME(F/E)VGXXXM")

    def test_collate_handles_gaps_within_match(self):
        a = Match("Q","T",1,3,1,6,0.0,100.0); a.target_sequence="ATGGAA"       # ME - but matching to 3 bps of query
        b = Match("Q","T",3,5,10,18,0.0,100.0); b.target_sequence="GAAGTGGGG"  # EVG
        c = Match("Q","T",9,9,30,32,0.0,100.0); c.target_sequence="ATG"        # M
        pm = ProteinHit([a,b,c],1,9,1,32)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "M(E/E)VGXXXM")

    def test_collate_handles_insertions_within_match(self):
        a = Match("Q","T",1,3,1,12,0.0,100.0); a.target_sequence="ATGGAATTTTTT" # MEFF - but matching to 3 bps of query
        b = Match("Q","T",3,5,10,18,0.0,100.0); b.target_sequence="GAAGTGGGG"   # EVG
        c = Match("Q","T",9,9,30,32,0.0,100.0); c.target_sequence="ATG"         # M
        pm = ProteinHit([a,b,c],1,9,1,32)
        collated = pm.collated_protein_sequence
        self.assertEqual(collated, "MEF(F/E)VGXXXM")


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
