import gzip
import os
import tempfile
import unittest
from pathlib import Path

from needle.ortholog import (
	ModuleDefinitionCycleError,
	load_modules,
)


def _write_text(path: Path, content: str) -> None:
	path.write_text(content, encoding="utf-8")


def _write_gz_text(path: Path, content: str) -> None:
	with gzip.open(path, "wt", encoding="utf-8") as f:
		f.write(content)


class TestOrthologModule(unittest.TestCase):
	def setUp(self) -> None:
		self.tmpdir = tempfile.TemporaryDirectory()
		self.tmp = Path(self.tmpdir.name)
		return super().setUp()

	def tearDown(self) -> None:
		self.tmpdir.cleanup()
		return super().tearDown()

	def _paths(self) -> dict:
		return {
			"modules_tsv": self.tmp / "modules.tsv",
			"orthologs_tsv": self.tmp / "ko_list.tsv",
			"module_defs_tsv": self.tmp / "module_ko_list.tsv",
			"ortholog_fasta": self.tmp / "ko_all_consensus.fasta.gz",
		}

	def test_same_module_object_for_orthologs(self):
		paths = self._paths()
		_write_text(
			paths["modules_tsv"],
			"Module ID\tModule Name\nM00001\tModule One\n",
		)
		_write_text(
			paths["orthologs_tsv"],
			"Ortholog ID\tOrtholog Name\nK00001\tKO One\nK00002\tKO Two\n",
		)
		_write_text(
			paths["module_defs_tsv"],
			"Module ID\tModule Definition\nM00001\tK00001,K00002\n",
		)
		_write_gz_text(
			paths["ortholog_fasta"],
			">K00001 some extra\nMSEQVENCEONE\n>K00002|more\nMSEQUENCE_TWO\n",
		)

		result = load_modules(
			paths["modules_tsv"],
			paths["orthologs_tsv"],
			paths["module_defs_tsv"],
			paths["ortholog_fasta"],
			module_id="M00001",
		)
		self.assertIn("M00001", result)
		mod = result["M00001"]
		self.assertEqual(mod.id, "M00001")
		self.assertEqual(len(mod.orthologs), 2)
		for o in mod.orthologs:
			self.assertIs(o.module, mod)

	def test_recursive_module_definition_expansion(self):
		paths = self._paths()
		_write_text(
			paths["modules_tsv"],
			"Module ID\tModule Name\nM00001\tModule One\nM00002\tModule Two\n",
		)
		_write_text(
			paths["orthologs_tsv"],
			"Ortholog ID\tOrtholog Name\nK00001\tKO One\nK00002\tKO Two\n",
		)
		# M00001 references K00001 and module M00002; M00002 has K00002
		_write_text(
			paths["module_defs_tsv"],
			"Module ID\tModule Definition\nM00001\tK00001,M00002\nM00002\tK00002\n",
		)
		_write_gz_text(
			paths["ortholog_fasta"],
			">K00001 x\nAAA\n>K00002 y\nBBB\n",
		)

		# Load only M00001, ensure it recursively includes K00002
		part = load_modules(
			paths["modules_tsv"],
			paths["orthologs_tsv"],
			paths["module_defs_tsv"],
			paths["ortholog_fasta"],
			module_id="M00001",
		)
		self.assertEqual(set(part.keys()), {"M00001"})
		m1 = part["M00001"]
		self.assertEqual({o.id for o in m1.orthologs}, {"K00001", "K00002"})

		# Load all and ensure M00002 contains only K00002
		full = load_modules(
			paths["modules_tsv"],
			paths["orthologs_tsv"],
			paths["module_defs_tsv"],
			paths["ortholog_fasta"],
		)
		self.assertIn("M00001", full)
		self.assertIn("M00002", full)
		self.assertEqual({o.id for o in full["M00002"].orthologs}, {"K00002"})

	def test_cycle_detection(self):
		paths = self._paths()
		_write_text(
			paths["modules_tsv"],
			"Module ID\tModule Name\nM10000\tA\nM10001\tB\n",
		)
		_write_text(
			paths["orthologs_tsv"],
			"Ortholog ID\tOrtholog Name\nK99999\tZ\n",
		)
		# Cycle: M10000 -> M10001 -> M10000
		_write_text(
			paths["module_defs_tsv"],
			"Module ID\tModule Definition\nM10000\tM10001\nM10001\tM10000\n",
		)
		_write_gz_text(
			paths["ortholog_fasta"],
			">K99999\nZZZ\n",
		)
		with self.assertRaises(ModuleDefinitionCycleError):
			load_modules(
				paths["modules_tsv"],
				paths["orthologs_tsv"],
				paths["module_defs_tsv"],
				paths["ortholog_fasta"],
				module_id="M10000",
			)

	def test_create_consensus_aa_fasta(self):
		paths = self._paths()
		_write_text(
			paths["modules_tsv"],
			"Module ID\tModule Name\nM00003\tModule Three\n",
		)
		_write_text(
			paths["orthologs_tsv"],
			"Ortholog ID\tOrtholog Name\nK01000\tKO A\nK01001\tKO B\nK01002\tKO C\n",
		)
		_write_text(
			paths["module_defs_tsv"],
			"Module ID\tModule Definition\nM00003\tK01002,K01000,K01001\n",
		)
		# Provide sequences for K01000 and K01002; leave K01001 missing to ensure it's skipped
		_write_gz_text(
			paths["ortholog_fasta"],
			">K01000 desc\nMAMA\n>K01002 other\nTTTTGG\n",
		)
		result = load_modules(
			paths["modules_tsv"],
			paths["orthologs_tsv"],
			paths["module_defs_tsv"],
			paths["ortholog_fasta"],
			module_id="M00003",
		)
		m3 = result["M00003"]
		out_path = self.tmp / "out.faa"
		m3.create_consensus_aa_fasta(out_path)
		text = out_path.read_text(encoding="utf-8").strip().splitlines()
		# Expect sorted by KO id: K01000 then K01002; K01001 missing (no sequence)
		self.assertEqual(text[0], ">K01000")
		self.assertEqual(text[1], "MAMA")
		self.assertEqual(text[2], ">K01002")
		self.assertEqual(text[3], "TTTTGG")


if __name__ == "__main__":
	unittest.main()


