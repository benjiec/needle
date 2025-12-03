import os
import tempfile
import unittest
import io

from needle.io import (
    write_protein_row,
    write_nucmatch_rows,
    write_fasta_record,
    export_protein_hits,
)
import needle.io as io_mod

from needle.blast import Results, group_matches, ProteinMatch, NucMatch


class TestIO(unittest.TestCase):
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
            orig = io_mod.get_hmmsearch_score_eval
            try:
                io_mod.get_hmmsearch_score_eval = lambda hmm, seq: (55.0, 1e-6)
                export_protein_hits("GENOMEZ", [pm1, pm2], p1, p2, p3)
            finally:
                io_mod.get_hmmsearch_score_eval = orig
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
            orig = io_mod.get_hmmsearch_score_eval
            try:
                io_mod.get_hmmsearch_score_eval = lambda hmm, seq: (10.0, 2e-5)
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
                io_mod.get_hmmsearch_score_eval = orig
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
