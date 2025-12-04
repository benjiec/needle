import os
import tempfile
import unittest

from needle.hmm import parse_hmmsearch_domtbl

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
            self.assertEqual(matches[0]["query_from"], 11)
            self.assertEqual(matches[0]["query_to"], 12)
            self.assertEqual(matches[0]["target_from"], 13)
            self.assertEqual(matches[0]["target_to"], 14)
            self.assertEqual(matches[1]["target_name"], "cand_0")
            self.assertEqual(matches[1]["evalue"], 0.2)
            self.assertEqual(matches[1]["score"], 101)
            self.assertEqual(matches[1]["query_from"], 21)
            self.assertEqual(matches[1]["query_to"], 22)
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


if __name__ == "__main__":
    unittest.main()
