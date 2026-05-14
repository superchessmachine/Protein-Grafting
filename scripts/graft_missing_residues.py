#!/usr/bin/env python3
"""
Build a hybrid protein model by retaining experimental coordinates wherever
present and filling unresolved residue positions from a complete model.

The script matches residues by polymer sequence position. It is intended for
cases where an experimental mmCIF/PDB and a complete predicted/model structure
share the same sequence register, possibly with an N-terminal construct prefix
that should be removed from the output.
"""

import argparse
import json
import math
import os
from dataclasses import dataclass, replace
from typing import Dict, Iterable, List, Optional, Tuple

import gemmi
import numpy as np


STANDARD_AA = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


@dataclass
class AtomRec:
    name: str
    resname: str
    chain: str
    seq_id: Optional[int]
    resseq: Optional[int]
    icode: str
    x: float
    y: float
    z: float
    occ: float
    b: float
    element: str
    charge: str = ""
    record: str = "ATOM"


@dataclass
class StructureData:
    sequence: Dict[int, str]
    numbering: Dict[int, int]
    atoms_by_seq: Dict[int, List[AtomRec]]
    heterogens: List[AtomRec]


def _clean_token(value) -> str:
    if value is None:
        return ""
    if value is False:
        return ""
    value = str(value).strip()
    return "" if value in {".", "?"} else value


def _as_int(value) -> Optional[int]:
    value = _clean_token(value)
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _infer_element(atom_name: str, element: str = "") -> str:
    element = _clean_token(element).upper()
    if element:
        return element
    stripped = atom_name.strip()
    if not stripped:
        return ""
    if len(stripped) >= 2 and stripped[:2].upper() in {"ZN", "FE", "MG", "MN", "CA", "CL", "NA"}:
        return stripped[:2].upper()
    return stripped[0].upper()


def _get_loop(block: gemmi.cif.Block, tags: Iterable[str]):
    try:
        return block.find(list(tags))
    except Exception:
        return []


def _get_category(block: gemmi.cif.Block, prefix: str) -> Dict[str, List[str]]:
    try:
        return block.get_mmcif_category(prefix)
    except Exception:
        return {}


def _category_rows(category: Dict[str, List[str]]):
    if not category:
        return
    length = max((len(values) for values in category.values()), default=0)
    for index in range(length):
        yield {
            key: values[index] if index < len(values) else "?"
            for key, values in category.items()
        }


def _find_entity_for_chain(block: gemmi.cif.Block, chain_id: str) -> Optional[str]:
    for row in _category_rows(_get_category(block, "_atom_site.")):
        if row.get("group_PDB") == "ATOM" and _clean_token(row.get("label_asym_id")) == chain_id:
            return _clean_token(row.get("label_entity_id"))
    return None


def _load_cif(path: str, chain_id: str, keep_heterogens: List[str]) -> StructureData:
    block = gemmi.cif.read(path).sole_block()
    entity_id = _find_entity_for_chain(block, chain_id)
    if entity_id is None:
        raise SystemExit(f"{path}: could not find polymer chain {chain_id!r}.")

    sequence: Dict[int, str] = {}
    for row in _get_loop(
        block,
        [
            "_entity_poly_seq.entity_id",
            "_entity_poly_seq.mon_id",
            "_entity_poly_seq.num",
        ],
    ):
        row_entity, mon_id, num = list(row)
        if _clean_token(row_entity) == entity_id:
            seq_id = _as_int(num)
            if seq_id is not None:
                sequence[seq_id] = _clean_token(mon_id).upper()

    numbering: Dict[int, int] = {}
    for row in _get_loop(
        block,
        [
            "_pdbx_poly_seq_scheme.asym_id",
            "_pdbx_poly_seq_scheme.seq_id",
            "_pdbx_poly_seq_scheme.pdb_seq_num",
            "_pdbx_poly_seq_scheme.auth_seq_num",
        ],
    ):
        asym_id, seq_id_raw, pdb_num_raw, auth_num_raw = list(row)
        if _clean_token(asym_id) != chain_id:
            continue
        seq_id = _as_int(seq_id_raw)
        if seq_id is None:
            continue
        numbering[seq_id] = _as_int(pdb_num_raw) or _as_int(auth_num_raw) or seq_id

    atoms_by_seq: Dict[int, List[AtomRec]] = {}
    heterogens: List[AtomRec] = []
    best_alt_rank = {"": 0, "A": 1}
    seen_alt: Dict[Tuple[int, str], str] = {}

    keep_set = {name.upper() for name in keep_heterogens}

    for row in _category_rows(_get_category(block, "_atom_site.")):
        group = _clean_token(row.get("group_PDB"))
        comp_id = _clean_token(row.get("label_comp_id")).upper()
        label_asym = _clean_token(row.get("label_asym_id"))
        atom_name = _clean_token(row.get("auth_atom_id")) or _clean_token(row.get("label_atom_id"))
        altloc = _clean_token(row.get("label_alt_id"))

        if group == "ATOM" and label_asym == chain_id:
            seq_id = _as_int(row.get("label_seq_id"))
            if seq_id is None:
                continue
            if altloc not in best_alt_rank:
                continue
            key = (seq_id, atom_name)
            previous = seen_alt.get(key)
            if previous is not None and best_alt_rank[altloc] >= best_alt_rank[previous]:
                continue
            if previous is not None:
                atoms_by_seq[seq_id] = [
                    atom for atom in atoms_by_seq.get(seq_id, []) if atom.name != atom_name
                ]
            seen_alt[key] = altloc
            atoms_by_seq.setdefault(seq_id, []).append(
                AtomRec(
                    name=atom_name,
                    resname=comp_id,
                    chain=chain_id,
                    seq_id=seq_id,
                    resseq=_as_int(row.get("auth_seq_id")) or numbering.get(seq_id) or seq_id,
                    icode=_clean_token(row.get("pdbx_PDB_ins_code")),
                    x=float(row.get("Cartn_x")),
                    y=float(row.get("Cartn_y")),
                    z=float(row.get("Cartn_z")),
                    occ=float(_clean_token(row.get("occupancy")) or 1.0),
                    b=float(_clean_token(row.get("B_iso_or_equiv")) or 0.0),
                    element=_infer_element(atom_name, row.get("type_symbol", "")),
                    charge=_clean_token(row.get("pdbx_formal_charge")),
                    record="ATOM",
                )
            )
        elif group == "HETATM" and (not keep_set or comp_id in keep_set):
            heterogens.append(
                AtomRec(
                    name=atom_name,
                    resname=comp_id,
                    chain=_clean_token(row.get("auth_asym_id")) or label_asym or chain_id,
                    seq_id=None,
                    resseq=_as_int(row.get("auth_seq_id")),
                    icode=_clean_token(row.get("pdbx_PDB_ins_code")),
                    x=float(row.get("Cartn_x")),
                    y=float(row.get("Cartn_y")),
                    z=float(row.get("Cartn_z")),
                    occ=float(_clean_token(row.get("occupancy")) or 1.0),
                    b=float(_clean_token(row.get("B_iso_or_equiv")) or 0.0),
                    element=_infer_element(atom_name, row.get("type_symbol", "")),
                    charge=_clean_token(row.get("pdbx_formal_charge")),
                    record="HETATM",
                )
            )

    if not sequence:
        sequence = {seq_id: atoms[0].resname for seq_id, atoms in atoms_by_seq.items()}

    return StructureData(
        sequence=sequence,
        numbering=numbering,
        atoms_by_seq=atoms_by_seq,
        heterogens=heterogens,
    )


def _parse_pdb_atom_line(line: str) -> AtomRec:
    name = line[12:16].strip()
    return AtomRec(
        name=name,
        resname=line[17:20].strip().upper(),
        chain=(line[21:22].strip() or "A"),
        seq_id=_as_int(line[22:26]),
        resseq=_as_int(line[22:26]),
        icode=line[26:27].strip(),
        x=float(line[30:38]),
        y=float(line[38:46]),
        z=float(line[46:54]),
        occ=float(line[54:60].strip() or 1.0),
        b=float(line[60:66].strip() or 0.0),
        element=_infer_element(name, line[76:78].strip()),
        charge=line[78:80].strip(),
        record=line[0:6].strip(),
    )


def _load_pdb(path: str, chain_id: str, keep_heterogens: List[str]) -> StructureData:
    atoms_by_seq: Dict[int, List[AtomRec]] = {}
    heterogens: List[AtomRec] = []
    sequence: Dict[int, str] = {}
    keep_set = {name.upper() for name in keep_heterogens}

    with open(path) as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom = _parse_pdb_atom_line(line)
            if atom.record == "ATOM" and atom.chain == chain_id and atom.seq_id is not None:
                atoms_by_seq.setdefault(atom.seq_id, []).append(atom)
                sequence.setdefault(atom.seq_id, atom.resname)
            elif atom.record == "HETATM" and (not keep_set or atom.resname in keep_set):
                heterogens.append(atom)

    numbering = {seq_id: seq_id for seq_id in sequence}
    if not atoms_by_seq:
        raise SystemExit(f"{path}: could not find polymer chain {chain_id!r}.")
    return StructureData(sequence, numbering, atoms_by_seq, heterogens)


def load_structure(path: str, chain_id: str, keep_heterogens: List[str]) -> StructureData:
    ext = os.path.splitext(path.lower())[1]
    if ext in {".cif", ".mmcif"}:
        return _load_cif(path, chain_id, keep_heterogens)
    return _load_pdb(path, chain_id, keep_heterogens)


def _atom_lookup(atoms: List[AtomRec], atom_name: str) -> Optional[AtomRec]:
    atom_name = atom_name.strip().upper()
    for atom in atoms:
        if atom.name.strip().upper() == atom_name:
            return atom
    return None


def _point(atom: AtomRec) -> np.ndarray:
    return np.array([atom.x, atom.y, atom.z], dtype=float)


def _fit_transform(
    experimental: StructureData,
    model: StructureData,
    seq_ids: Iterable[int],
    atom_names: List[str],
) -> Tuple[np.ndarray, np.ndarray, float, int]:
    model_points = []
    experimental_points = []
    for seq_id in seq_ids:
        exp_atoms = experimental.atoms_by_seq.get(seq_id)
        model_atoms = model.atoms_by_seq.get(seq_id)
        if not exp_atoms or not model_atoms:
            continue
        for atom_name in atom_names:
            exp_atom = _atom_lookup(exp_atoms, atom_name)
            model_atom = _atom_lookup(model_atoms, atom_name)
            if exp_atom is not None and model_atom is not None:
                experimental_points.append(_point(exp_atom))
                model_points.append(_point(model_atom))

    if len(model_points) < 3:
        raise SystemExit("Need at least three shared atoms to align model coordinates.")

    p = np.asarray(model_points)
    q = np.asarray(experimental_points)
    p_centroid = p.mean(axis=0)
    q_centroid = q.mean(axis=0)
    h = (p - p_centroid).T @ (q - q_centroid)
    u, _, vt = np.linalg.svd(h)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T
    translation = q_centroid - p_centroid @ rotation.T
    transformed = p @ rotation.T + translation
    rmsd = math.sqrt(float(np.mean(np.sum((transformed - q) ** 2, axis=1))))
    return rotation, translation, rmsd, len(model_points)


def _identity_transform() -> Tuple[np.ndarray, np.ndarray]:
    return np.eye(3), np.zeros(3)


def _transform_atom(atom: AtomRec, rotation: np.ndarray, translation: np.ndarray) -> AtomRec:
    xyz = np.array([atom.x, atom.y, atom.z], dtype=float) @ rotation.T + translation
    return replace(atom, x=float(xyz[0]), y=float(xyz[1]), z=float(xyz[2]))


def _format_atom_line(serial: int, atom: AtomRec, chain: str, resseq: int, record: str) -> str:
    name = atom.name.strip()
    atom_name_fmt = name.rjust(4)
    resname = atom.resname.strip().upper()[:3]
    chain = (chain or "A")[:1]
    icode = (atom.icode or " ")[:1]
    element = _infer_element(name, atom.element).rjust(2)
    charge = (atom.charge or "").rjust(2)
    return (
        f"{record:<6s}{serial:5d} {atom_name_fmt} {resname:>3s} {chain}"
        f"{resseq:4d}{icode}   "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"{atom.occ:6.2f}{atom.b:6.2f}          {element}{charge}\n"
    )


def _distance(atom_a: Optional[AtomRec], atom_b: Optional[AtomRec]) -> Optional[float]:
    if atom_a is None or atom_b is None:
        return None
    return float(np.linalg.norm(_point(atom_a) - _point(atom_b)))


def _collapse_ranges(values: List[int]) -> List[Dict[str, int]]:
    if not values:
        return []
    values = sorted(values)
    ranges = []
    start = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append({"start": start, "end": previous, "count": previous - start + 1})
        start = previous = value
    ranges.append({"start": start, "end": previous, "count": previous - start + 1})
    return ranges


def _parse_atom_names(value: str) -> List[str]:
    names = [name.strip().upper() for name in value.split(",") if name.strip()]
    return names or ["CA"]


def _parse_ranges(values: List[str]) -> List[Tuple[int, int]]:
    ranges: List[Tuple[int, int]] = []
    for value in values:
        value = value.strip()
        if not value:
            continue
        if "-" in value:
            start_raw, end_raw = value.split("-", 1)
            start = int(start_raw)
            end = int(end_raw)
        else:
            start = end = int(value)
        if end < start:
            start, end = end, start
        ranges.append((start, end))
    return ranges


def _seq_in_ranges(seq_id: int, ranges: List[Tuple[int, int]]) -> Optional[Tuple[int, int]]:
    for start, end in ranges:
        if start <= seq_id <= end:
            return start, end
    return None


def build_hybrid(args: argparse.Namespace) -> Dict:
    keep_heterogens = [name.strip().upper() for name in args.keep_heterogen if name.strip()]
    experimental = load_structure(args.experimental, args.experimental_chain, keep_heterogens)
    model = load_structure(args.model, args.model_chain, [])
    model_ranges = _parse_ranges(args.model_range)

    seq_ids = [
        seq_id
        for seq_id in sorted(experimental.sequence)
        if seq_id > args.drop_prefix and experimental.sequence[seq_id] in STANDARD_AA
    ]
    if not seq_ids:
        raise SystemExit("No output residues remain after applying --drop-prefix.")

    missing_seq_ids = [seq_id for seq_id in seq_ids if seq_id not in experimental.atoms_by_seq]
    filled_seq_ids = [seq_id for seq_id in missing_seq_ids if seq_id in model.atoms_by_seq]
    unfilled_seq_ids = [seq_id for seq_id in missing_seq_ids if seq_id not in model.atoms_by_seq]

    if args.no_align:
        rotation, translation = _identity_transform()
        align_rmsd = None
        align_atom_count = 0
    else:
        align_seq_ids = [seq_id for seq_id in seq_ids if seq_id in experimental.atoms_by_seq and seq_id in model.atoms_by_seq]
        rotation, translation, align_rmsd, align_atom_count = _fit_transform(
            experimental,
            model,
            align_seq_ids,
            _parse_atom_names(args.align_atoms),
        )

    range_transforms: Dict[Tuple[int, int], Tuple[np.ndarray, np.ndarray]] = {}
    local_alignments = []
    for start, end in model_ranges:
        if args.local_align_flank <= 0:
            range_transforms[(start, end)] = (rotation, translation)
            local_alignments.append(
                {
                    "range": {"start": start, "end": end, "count": end - start + 1},
                    "enabled": False,
                }
            )
            continue

        anchor_seq_ids = [
            seq_id
            for seq_id in list(range(start - args.local_align_flank, start))
            + list(range(end + 1, end + args.local_align_flank + 1))
            if seq_id in experimental.atoms_by_seq and seq_id in model.atoms_by_seq
        ]
        try:
            local_rotation, local_translation, local_rmsd, local_atom_count = _fit_transform(
                experimental,
                model,
                anchor_seq_ids,
                _parse_atom_names(args.local_align_atoms),
            )
            range_transforms[(start, end)] = (local_rotation, local_translation)
            local_alignments.append(
                {
                    "range": {"start": start, "end": end, "count": end - start + 1},
                    "enabled": True,
                    "flank": args.local_align_flank,
                    "anchor_seq_ids": anchor_seq_ids,
                    "atom_count": local_atom_count,
                    "rmsd_angstrom": round(local_rmsd, 3),
                    "atoms": _parse_atom_names(args.local_align_atoms),
                }
            )
        except SystemExit as exc:
            range_transforms[(start, end)] = (rotation, translation)
            local_alignments.append(
                {
                    "range": {"start": start, "end": end, "count": end - start + 1},
                    "enabled": True,
                    "fallback": "global",
                    "reason": str(exc),
                }
            )

    out_lines = []
    out_lines.append("REMARK Hybrid model: experimental coordinates retained where present.\n")
    out_lines.append("REMARK Missing polymer residues filled from the supplied complete model.\n")
    out_lines.append(f"REMARK Dropped polymer sequence positions <= {args.drop_prefix}.\n")
    if align_rmsd is not None:
        out_lines.append(f"REMARK Model-to-experimental alignment RMSD {align_rmsd:.3f} A over {align_atom_count} atoms.\n")
    for start, end in model_ranges:
        out_lines.append(f"REMARK Forced model-coordinate range: sequence positions {start}-{end}.\n")

    output_atoms_by_seq: Dict[int, List[AtomRec]] = {}
    sources: Dict[int, str] = {}
    serial = 1

    for seq_id in seq_ids:
        forced_range = _seq_in_ranges(seq_id, model_ranges)
        if forced_range is not None and seq_id in model.atoms_by_seq:
            local_rotation, local_translation = range_transforms.get(forced_range, (rotation, translation))
            source_atoms = [
                _transform_atom(atom, local_rotation, local_translation)
                for atom in model.atoms_by_seq[seq_id]
            ]
            source = "model"
        elif seq_id in experimental.atoms_by_seq:
            source_atoms = experimental.atoms_by_seq[seq_id]
            source = "experimental"
        elif seq_id in model.atoms_by_seq:
            source_atoms = [
                _transform_atom(atom, rotation, translation)
                for atom in model.atoms_by_seq[seq_id]
            ]
            source = "model"
        else:
            continue

        if args.numbering == "compact":
            resseq = seq_ids.index(seq_id) + 1
        elif args.numbering == "label":
            resseq = seq_id
        else:
            resseq = experimental.numbering.get(seq_id, seq_id - args.drop_prefix)

        output_atoms_by_seq[seq_id] = source_atoms
        sources[seq_id] = source
        for atom in source_atoms:
            atom_out = replace(atom, resname=experimental.sequence.get(seq_id, atom.resname))
            out_lines.append(_format_atom_line(serial, atom_out, args.output_chain, resseq, "ATOM"))
            serial += 1

    out_lines.append("TER\n")

    kept_heterogens = []
    for atom in experimental.heterogens:
        if atom.resseq is None:
            continue
        chain = args.heterogen_chain or atom.chain or args.output_chain
        out_lines.append(_format_atom_line(serial, atom, chain, atom.resseq, "HETATM"))
        kept_heterogens.append({"resname": atom.resname, "chain": chain, "resseq": atom.resseq, "atom": atom.name})
        serial += 1

    out_lines.append("END\n")

    with open(args.out, "w") as handle:
        handle.writelines(out_lines)

    junctions = []
    for left, right in zip(seq_ids, seq_ids[1:]):
        if left not in output_atoms_by_seq or right not in output_atoms_by_seq:
            continue
        left_atoms = output_atoms_by_seq[left]
        right_atoms = output_atoms_by_seq[right]
        c_n = _distance(_atom_lookup(left_atoms, "C"), _atom_lookup(right_atoms, "N"))
        ca_ca = _distance(_atom_lookup(left_atoms, "CA"), _atom_lookup(right_atoms, "CA"))
        if c_n is None:
            continue
        if c_n < args.min_peptide_bond or c_n > args.max_peptide_bond:
            junctions.append(
                {
                    "left_seq_id": left,
                    "right_seq_id": right,
                    "left_number": experimental.numbering.get(left, left - args.drop_prefix),
                    "right_number": experimental.numbering.get(right, right - args.drop_prefix),
                    "left_source": sources.get(left),
                    "right_source": sources.get(right),
                    "status": "short" if c_n < args.min_peptide_bond else "long",
                    "c_n_distance_angstrom": round(c_n, 3),
                    "ca_ca_distance_angstrom": round(ca_ca, 3) if ca_ca is not None else None,
                }
            )

    report = {
        "experimental": args.experimental,
            "model": args.model,
            "output": args.out,
            "drop_prefix": args.drop_prefix,
            "output_residue_count": len(output_atoms_by_seq),
            "experimental_residue_count": sum(1 for seq_id in seq_ids if seq_id in experimental.atoms_by_seq),
            "filled_residue_count": len(filled_seq_ids),
            "overwritten_experimental_residue_count": sum(
                1 for seq_id in seq_ids if seq_id in experimental.atoms_by_seq and _seq_in_ranges(seq_id, model_ranges)
            ),
            "unfilled_residue_count": len(unfilled_seq_ids),
            "filled_ranges": _collapse_ranges(filled_seq_ids),
            "overwritten_experimental_ranges": _collapse_ranges(
                [seq_id for seq_id in seq_ids if seq_id in experimental.atoms_by_seq and _seq_in_ranges(seq_id, model_ranges)]
            ),
            "unfilled_ranges": _collapse_ranges(unfilled_seq_ids),
            "model_ranges": [
                {"start": start, "end": end, "count": end - start + 1}
                for start, end in model_ranges
            ],
            "retained_heterogens": kept_heterogens,
            "alignment": {
                "enabled": not args.no_align,
                "atom_count": align_atom_count,
                "rmsd_angstrom": round(align_rmsd, 3) if align_rmsd is not None else None,
                "atoms": _parse_atom_names(args.align_atoms),
            },
            "local_alignments": local_alignments,
            "junction_warnings": junctions,
            "min_peptide_bond_angstrom": args.min_peptide_bond,
            "max_peptide_bond_angstrom": args.max_peptide_bond,
        }

    if args.report:
        with open(args.report, "w") as handle:
            json.dump(report, handle, indent=2)
            handle.write("\n")

    return report


def main():
    parser = argparse.ArgumentParser(
        description="Retain experimental coordinates and graft missing residues from a complete model."
    )
    parser.add_argument("--experimental", required=True, help="Experimental mmCIF/PDB structure.")
    parser.add_argument("--model", required=True, help="Complete model mmCIF/PDB structure in the same sequence register.")
    parser.add_argument("--out", required=True, help="Output hybrid PDB.")
    parser.add_argument("--report", help="Optional JSON report path.")
    parser.add_argument("--experimental-chain", default="A", help="Experimental polymer chain ID (default: A).")
    parser.add_argument("--model-chain", default="A", help="Model polymer chain ID (default: A).")
    parser.add_argument("--output-chain", default="A", help="Output polymer chain ID (default: A).")
    parser.add_argument("--drop-prefix", type=int, default=0, help="Drop sequence positions <= this value.")
    parser.add_argument(
        "--numbering",
        choices=["scheme", "label", "compact"],
        default="scheme",
        help="Output residue numbering: experimental scheme, sequence labels, or compact 1..N.",
    )
    parser.add_argument(
        "--keep-heterogen",
        action="append",
        default=[],
        help="Retain experimental heterogen residue name; repeat for multiple names. If omitted, all heterogens are retained.",
    )
    parser.add_argument("--heterogen-chain", help="Override output chain ID for retained heterogens.")
    parser.add_argument(
        "--model-range",
        action="append",
        default=[],
        metavar="START-END",
        help="Force this sequence-position range to use model coordinates even where experimental atoms exist; repeat as needed.",
    )
    parser.add_argument(
        "--local-align-flank",
        type=int,
        default=0,
        help="For each --model-range, align model coordinates using this many flanking residues on each side.",
    )
    parser.add_argument(
        "--local-align-atoms",
        default="N,CA,C,O",
        help="Comma-separated atom names used for local model-range superposition (default: N,CA,C,O).",
    )
    parser.add_argument("--no-align", action="store_true", help="Do not superpose model coordinates onto experimental coordinates.")
    parser.add_argument(
        "--align-atoms",
        default="CA",
        help="Comma-separated atom names used for model-to-experimental superposition (default: CA).",
    )
    parser.add_argument(
        "--min-peptide-bond",
        type=float,
        default=1.1,
        help="C-N distance threshold for reporting short junction warnings in Angstrom (default: 1.1).",
    )
    parser.add_argument(
        "--max-peptide-bond",
        type=float,
        default=1.8,
        help="C-N distance threshold for reporting junction warnings in Angstrom (default: 1.8).",
    )
    args = parser.parse_args()

    report = build_hybrid(args)
    print(f"Wrote hybrid model: {args.out}")
    if args.report:
        print(f"Wrote graft report: {args.report}")
    print(f"Filled residues from model: {report['filled_residue_count']}")
    if report["junction_warnings"]:
        print(f"Junction warnings: {len(report['junction_warnings'])}")


if __name__ == "__main__":
    main()
