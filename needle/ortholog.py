from __future__ import annotations

import csv
import gzip
import io
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union


class ModuleDefinitionCycleError(RuntimeError):
	"""Raised when a cycle is detected while expanding module definitions."""


_KO_ID_RE = re.compile(r"K\d{5}")
_MODULE_ID_RE = re.compile(r"M\d{5}")


@dataclass
class Module:
	id: str
	name: Optional[str] = None
	orthologs: List["Ortholog"] = field(default_factory=list)

	def create_consensus_aa_fasta(self, path: Union[str, Path]) -> None:
		output_path = Path(path)
		records = []
		for ortholog in sorted(self.orthologs, key=lambda o: o.id):
			if ortholog.consensus_aa_sequence:
				records.append((ortholog.id, ortholog.consensus_aa_sequence))
		with output_path.open("w", encoding="utf-8") as f:
			for ortholog_id, seq in records:
				f.write(f">{ortholog_id}\n")
				# Wrap to 60 chars per FASTA line for readability
				for i in range(0, len(seq), 60):
					f.write(seq[i : i + 60] + "\n")


@dataclass
class Ortholog:
	id: str
	name: Optional[str]
	consensus_aa_sequence: Optional[str]
	module: Module


def _read_modules_tsv(path: Union[str, Path]) -> Dict[str, Optional[str]]:
	"""
	Returns mapping: Module ID -> Module Name (name may be None).
	Expected headers (from example): 'Module ID', 'Module Name'
	"""
	result: Dict[str, Optional[str]] = {}
	with Path(path).open("r", encoding="utf-8") as f:
		reader = csv.DictReader(f, delimiter="\t")
		if "Module ID" not in reader.fieldnames or "Module Name" not in reader.fieldnames:
			raise ValueError("modules.tsv must have headers: 'Module ID' and 'Module Name'")
		for row in reader:
			module_id = (row.get("Module ID") or "").strip()
			if not module_id:
				continue
			module_name = (row.get("Module Name") or "").strip() or None
			result[module_id] = module_name
	return result


def _read_orthologs_tsv(path: Union[str, Path]) -> Dict[str, Optional[str]]:
	"""
	Returns mapping: Ortholog ID -> Ortholog Name (name may be None).
	Expected headers (from example): 'Ortholog ID', 'Ortholog Name'
	"""
	result: Dict[str, Optional[str]] = {}
	with Path(path).open("r", encoding="utf-8") as f:
		reader = csv.DictReader(f, delimiter="\t")
		if "Ortholog ID" not in reader.fieldnames or "Ortholog Name" not in reader.fieldnames:
			raise ValueError("orthologs_tsv must have headers: 'Ortholog ID' and 'Ortholog Name'")
		for row in reader:
			ko_id = (row.get("Ortholog ID") or "").strip()
			if not ko_id:
				continue
			ko_name = (row.get("Ortholog Name") or "").strip() or None
			result[ko_id] = ko_name
	return result


def _read_module_defs_tsv(path: Union[str, Path]) -> Dict[str, List[str]]:
	"""
	Returns mapping: Module ID -> list of tokens (KO IDs or Module IDs).
	Expected headers (from example): 'Module ID', 'Module Definition'
	"""
	result: Dict[str, List[str]] = {}
	with Path(path).open("r", encoding="utf-8") as f:
		reader = csv.DictReader(f, delimiter="\t")
		if "Module ID" not in reader.fieldnames or "Module Definition" not in reader.fieldnames:
			raise ValueError("module_ko_list.tsv must have headers: 'Module ID' and 'Module Definition'")
		for row in reader:
			module_id = (row.get("Module ID") or "").strip()
			if not module_id:
				continue
			defn = (row.get("Module Definition") or "").strip()
			tokens = [t.strip() for t in defn.split(",") if t.strip()] if defn else []
			result[module_id] = tokens
	return result


def _open_fasta_maybe_gzip(path: Union[str, Path]) -> io.TextIOBase:
	p = Path(path)
	if p.suffix == ".gz":
		return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8")
	return p.open("r", encoding="utf-8")


def _read_ortholog_consensus_fasta(path: Union[str, Path]) -> Dict[str, str]:
	"""
	Reads a FASTA where the sequence ID contains the ortholog ID (KO) somewhere
	in the header. Additional content may follow the ID. Supports .gz files.
	Returns mapping: KO ID -> consensus amino acid sequence.
	"""
	sequences: Dict[str, str] = {}
	current_id: Optional[str] = None
	current_seq_parts: List[str] = []

	with _open_fasta_maybe_gzip(path) as f:
		for raw_line in f:
			line = raw_line.rstrip("\n\r")
			if not line:
				continue
			if line.startswith(">"):
				# Flush previous record
				if current_id is not None:
					if current_seq_parts:
						sequences.setdefault(current_id, "".join(current_seq_parts))
					current_seq_parts = []
				# Extract KO id from header using regex
				match = _KO_ID_RE.search(line)
				current_id = match.group(0) if match else None
			else:
				if current_id is not None:
					current_seq_parts.append(line.strip())
	# Flush last record
	if current_id is not None and current_seq_parts:
		sequences.setdefault(current_id, "".join(current_seq_parts))
	return sequences


def _expand_module_ko_set(
	module_id: str,
	module_defs: Dict[str, List[str]],
	memo: Dict[str, Set[str]],
	stack: Set[str],
) -> Set[str]:
	if module_id in memo:
		return memo[module_id]
	if module_id in stack:
		raise ModuleDefinitionCycleError(f"Cycle detected while expanding module {module_id}")
	stack.add(module_id)
	tokens = module_defs.get(module_id, [])
	result: Set[str] = set()
	for token in tokens:
		if _KO_ID_RE.fullmatch(token):
			result.add(token)
		elif _MODULE_ID_RE.fullmatch(token):
			child_set = _expand_module_ko_set(token, module_defs, memo, stack)
			if child_set:
				result.update(child_set)
		else:
			raise ValueError(f"Unrecognized token in module definition for {module_id}: {token!r}")
	stack.remove(module_id)
	memo[module_id] = result
	return result


def load_modules(
	modules_tsv: Union[str, Path],
	orthologs_tsv: Union[str, Path],
	module_defs_tsv: Union[str, Path],
	ortholog_fasta: Union[str, Path],
	module_id: Optional[str] = None,
) -> Dict[str, Module]:
	"""
	Load modules and orthologs from provided files.

	- modules_tsv: TSV with headers 'Module ID', 'Module Name'
	- orthologs_tsv: TSV with headers 'Ortholog ID', 'Ortholog Name'
	- module_defs_tsv: TSV with headers 'Module ID', 'Module Definition' (comma-delimited list of KO IDs or Module IDs)
	- ortholog_fasta: FASTA (optionally .gz) with sequence IDs containing KO IDs; extra header content allowed
	- module_id: if provided, only load this module (its KO set expands recursively)

	Returns mapping: Module ID -> Module instance.
	"""
	module_names = _read_modules_tsv(modules_tsv)
	ko_names = _read_orthologs_tsv(orthologs_tsv)
	module_defs = _read_module_defs_tsv(module_defs_tsv)
	ko_sequences = _read_ortholog_consensus_fasta(ortholog_fasta)

	memo: Dict[str, Set[str]] = {}
	all_module_ids: Set[str] = set(module_names.keys()) | set(module_defs.keys())

	if module_id is not None:
		target_ids = {module_id}
	else:
		target_ids = all_module_ids

	# Compute KO sets for target modules
	module_to_kos: Dict[str, Set[str]] = {}
	for mid in sorted(target_ids):
		kos = _expand_module_ko_set(mid, module_defs, memo, set())
		module_to_kos[mid] = kos

	# Build Module and Ortholog objects
	result: Dict[str, Module] = {}
	for mid in sorted(module_to_kos.keys()):
		module_obj = Module(id=mid, name=module_names.get(mid))
		result[mid] = module_obj

	for mid, ko_set in module_to_kos.items():
		module_obj = result[mid]
		for ko_id in sorted(ko_set):
			ortholog_obj = Ortholog(
				id=ko_id,
				name=ko_names.get(ko_id),
				consensus_aa_sequence=ko_sequences.get(ko_id),
				module=module_obj,
			)
			module_obj.orthologs.append(ortholog_obj)

	return result


