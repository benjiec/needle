import os
import tempfile
import unittest

import importlib


class TestGFFParser(unittest.TestCase):
    def setUp(self):
        # Skip tests if gffutils is not available
        if importlib.util.find_spec("gffutils") is None:
            self.skipTest("gffutils not installed")

    def _write_gff(self, path: str, content: str):
        with open(path, "w") as f:
            f.write(content)

    def test_parse_simple_plus_and_minus(self):
        from needle.gff import parse_gff_to_hits

        gff_text = "\n".join(
            [
                "##gff-version 3",
                # Protein A on chr1, '+' strand
                # CDS1: 51 nt -> 17 aa (phase 0) => aa 1..17
                "chr1\t.\tCDS\t100\t150\t.\t+\t0\tprotein_id=protA",
                # CDS2: 61 nt -> (61-2)=59 -> 19 aa (phase 2) => aa 18..36
                "chr1\t.\tCDS\t200\t260\t.\t+\t2\tprotein_id=protA",
                # CDS3: 9 nt -> (9-1)=8 -> 2 aa (phase 1) => aa 37..38
                "chr1\t.\tCDS\t300\t308\t.\t+\t1\tprotein_id=protA",
                # Protein B on chr2, '-' strand
                # CDS order in transcript is descending by start for '-'
                # B1: 31 nt -> 10 aa (phase 0) => aa 1..10
                "chr2\t.\tCDS\t900\t930\t.\t-\t0\tprotein_id=protB",
                # B2: 11 nt -> (11-2)=9 -> 3 aa (phase 2) => aa 11..13
                "chr2\t.\tCDS\t850\t860\t.\t-\t2\tprotein_id=protB",
            ]
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            gff_path = os.path.join(tmpdir, "input.gff3")
            self._write_gff(gff_path, gff_text)

            hits = parse_gff_to_hits(gff_path)
            by_pid = {h._protein_hit_id: h for h in hits}

            # protA checks
            self.assertIn("protA", by_pid)
            pa = by_pid["protA"]
            self.assertEqual(pa.query_start, 1)
            self.assertEqual(pa.query_end, 38)
            self.assertEqual(pa.target_accession, "chr1")
            self.assertFalse(pa.on_reverse_strand)
            self.assertEqual(pa.target_start, 100)
            self.assertEqual(pa.target_end, 308)
            # matches aa spans
            self.assertEqual(len(pa.matches), 3)
            self.assertEqual((pa.matches[0].query_start, pa.matches[0].query_end), (1, 17))
            self.assertEqual((pa.matches[1].query_start, pa.matches[1].query_end), (18, 36))
            self.assertEqual((pa.matches[2].query_start, pa.matches[2].query_end), (37, 38))
            # target orientation within matches
            for m in pa.matches:
                self.assertLessEqual(m.target_start, m.target_end)

            # protB checks
            self.assertIn("protB", by_pid)
            pb = by_pid["protB"]
            self.assertEqual(pb.query_start, 1)
            self.assertEqual(pb.query_end, 13)
            self.assertEqual(pb.target_accession, "chr2")
            self.assertTrue(pb.on_reverse_strand)
            # reverse strand -> start > end (5' to 3' of gene)
            self.assertGreater(pb.target_start, pb.target_end)
            self.assertEqual(pb.target_start, 930)
            self.assertEqual(pb.target_end, 850)
            # matches aa spans and reversed target coords
            self.assertEqual(len(pb.matches), 2)
            self.assertEqual((pb.matches[0].query_start, pb.matches[0].query_end), (1, 10))
            self.assertEqual((pb.matches[1].query_start, pb.matches[1].query_end), (11, 13))
            for m in pb.matches:
                self.assertGreater(m.target_start, m.target_end)

    def test_ignores_cds_without_protein_id(self):
        from needle.gff import parse_gff_to_hits
        gff_text = "\n".join(
            [
                "##gff-version 3",
                "chr1\t.\tCDS\t10\t40\t.\t+\t0\tID=cds1",
                "chr1\t.\tCDS\t50\t80\t.\t+\t0\tprotein_id=with_id",
            ]
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            gff_path = os.path.join(tmpdir, "input.gff3")
            self._write_gff(gff_path, gff_text)
            hits = parse_gff_to_hits(gff_path)
            self.assertEqual(len(hits), 1)
            self.assertEqual(hits[0]._protein_hit_id, "with_id")


