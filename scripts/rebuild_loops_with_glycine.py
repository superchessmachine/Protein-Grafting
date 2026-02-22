#!/usr/bin/env python3
"""
Rebuild loop gaps with glycine residues using PDBFixer.

This script scans each protein chain for missing internal residues, replaces the
templates for those residues with glycine, and lets PDBFixer build coordinates
to bridge the gaps. If the PDB file lacks sequence metadata, the script infers
missing loops by looking for residue-number gaps and large spatial jumps between
adjacent Cα atoms. By default, terminal gaps are skipped since they are not true
loops, but you can opt in to rebuilding them as well.

Example:
    python scripts/rebuild_loops_with_glycine.py \\
        --in data/Final_MARTINI_Fan.pdb \\
        --out data/Final_MARTINI_Fan_loops.pdb
"""

import argparse
import copy
import json
import math
from typing import Dict, Iterable, List, Tuple, Optional

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from openmm import unit


def _chain_label(chain, index: int) -> str:
    """
    Return a human-readable identifier for a chain.
    """
    cid = getattr(chain, "id", None)
    if cid is None:
        return str(index)
    cid = str(cid).strip()
    return cid if cid else str(index)


def _residue_label(res) -> str:
    """
    Build a residue label like GLY42 for summary messages.
    """
    if res is None:
        return ""
    rid = getattr(res, "id", "").strip()
    name = getattr(res, "name", "UNK").strip() or "UNK"
    return f"{name}{rid}" if rid else name


def _force_residue(entry, residue_name: str):
    """
    Return an object representing the same template but using the desired residue name.
    """
    if isinstance(entry, str):
        return residue_name.upper()

    new_entry = copy.copy(entry)
    for attr in ("name", "residueName", "resname"):
        if hasattr(new_entry, attr):
            try:
                setattr(new_entry, attr, residue_name.upper())
                return new_entry
            except Exception:
                pass
    return residue_name.upper()


def _pos_to_tuple_nm(pos) -> Tuple[float, float, float]:
    if hasattr(pos, "value_in_unit"):
        v = pos.value_in_unit(unit.nanometer)
        return float(v[0]), float(v[1]), float(v[2])
    return float(pos[0]), float(pos[1]), float(pos[2])


def _distance_angstrom(pos_a, pos_b) -> float:
    ax, ay, az = _pos_to_tuple_nm(pos_a)
    bx, by, bz = _pos_to_tuple_nm(pos_b)
    dx = ax - bx
    dy = ay - by
    dz = az - bz
    return math.sqrt(dx * dx + dy * dy + dz * dz) * 10.0


def _build_atom_position_map(fixer: PDBFixer) -> Dict:
    atoms = list(fixer.topology.atoms())
    positions = list(fixer.positions)
    return {atom: positions[i] for i, atom in enumerate(atoms)}


def _get_atom_position(atom_positions: Dict, residue, atom_name: str):
    atom_name = atom_name.strip().upper()
    for atom in residue.atoms():
        if atom.name.strip().upper() == atom_name:
            return atom_positions.get(atom)
    return None


def _parse_resseq(res_id: Optional[str]) -> Optional[int]:
    if not res_id:
        return None
    res_id = str(res_id).strip()
    if not res_id:
        return None

    digits = []
    for ch in res_id:
        if ch.isdigit():
            digits.append(ch)
        elif ch == "-" and not digits:
            digits.append(ch)
        elif digits:
            break
    if not digits:
        return None
    try:
        return int("".join(digits))
    except ValueError:
        return None


def _collect_chain_infos(
    chains: List,
    chain_residues: Dict[int, List],
    atom_positions: Dict,
) -> List[List[Dict]]:
    chain_infos: List[List[Dict]] = []
    for chain_index, chain in enumerate(chains):
        residues = chain_residues.get(chain_index, [])
        info_list = []
        for res in residues:
            info_list.append(
                {
                    "residue": res,
                    "seq": _parse_resseq(getattr(res, "id", "")),
                    "ca": _get_atom_position(atom_positions, res, "CA"),
                    "atoms": [
                        (
                            atom.name,
                            atom.element,
                            atom_positions.get(atom),
                        )
                        for atom in res.atoms()
                    ],
                }
            )
        chain_infos.append(info_list)
    return chain_infos


def _extend_missing_entries(
    missing: Dict[Tuple[int, int], List],
    additions: Dict[Tuple[int, int], int],
    residue_name: str,
):
    for key, count in additions.items():
        if count <= 0:
            continue
        missing.setdefault(key, [])
        missing[key].extend([residue_name.upper()] * count)


def _detect_numbering_gaps(
    chain_infos: List[List[Dict]],
    base_counts: Dict[Tuple[int, int], int],
) -> Tuple[Dict[Tuple[int, int], int], List[Dict]]:
    additions: Dict[Tuple[int, int], int] = {}
    logs: List[Dict] = []

    for chain_index, residues in enumerate(chain_infos):
        if len(residues) < 2:
            continue

        prev_seq = residues[0]["seq"]
        for idx in range(1, len(residues)):
            curr_seq = residues[idx]["seq"]
            if prev_seq is None or curr_seq is None:
                prev_seq = curr_seq
                continue

            gap = curr_seq - prev_seq - 1
            if gap > 0:
                key = (chain_index, idx)
                existing = base_counts.get(key, 0) + additions.get(key, 0)
                additions[key] = additions.get(key, 0) + gap
                logs.append(
                    {
                        "type": "numbering",
                        "chain_index": chain_index,
                        "insert_index": idx,
                        "before": residues[idx - 1]["residue"],
                        "after": residues[idx]["residue"],
                        "added": gap,
                        "total": existing + gap,
                        "details": {"prev_seq": prev_seq, "next_seq": curr_seq},
                    }
                )
            prev_seq = curr_seq

    return additions, logs


def _detect_distance_gaps(
    chain_infos: List[List[Dict]],
    base_counts: Dict[Tuple[int, int], int],
    threshold: float,
    ideal_spacing: float,
) -> Tuple[Dict[Tuple[int, int], int], List[Dict]]:
    additions: Dict[Tuple[int, int], int] = {}
    logs: List[Dict] = []

    for chain_index, residues in enumerate(chain_infos):
        if len(residues) < 2:
            continue

        for idx in range(len(residues) - 1):
            ca_a = residues[idx]["ca"]
            ca_b = residues[idx + 1]["ca"]
            if ca_a is None or ca_b is None:
                continue

            distance = _distance_angstrom(ca_a, ca_b)
            if distance <= threshold:
                continue

            est_missing = max(1, int(round(distance / ideal_spacing)) - 1)
            key = (chain_index, idx + 1)
            existing = base_counts.get(key, 0) + additions.get(key, 0)
            if est_missing <= existing:
                continue

            add = est_missing - existing
            additions[key] = additions.get(key, 0) + add
            logs.append(
                {
                    "type": "distance",
                    "chain_index": chain_index,
                    "insert_index": idx + 1,
                    "before": residues[idx]["residue"],
                    "after": residues[idx + 1]["residue"],
                    "added": add,
                    "total": existing + add,
                    "details": {
                        "distance": distance,
                        "threshold": threshold,
                        "estimated_total": est_missing,
                    },
                }
            )

    return additions, logs


def _format_detection_log(entry: Dict, chains: List) -> str:
    chain_label = _chain_label(chains[entry["chain_index"]], entry["chain_index"])
    before_label = _residue_label(entry.get("before"))
    after_label = _residue_label(entry.get("after"))
    added = entry.get("added", 0)
    total = entry.get("total", added)

    if entry.get("type") == "numbering":
        prev_seq = entry["details"].get("prev_seq")
        next_seq = entry["details"].get("next_seq")
        return (
            f"Detected numbering gap in chain {chain_label}: residues jump from {prev_seq} to {next_seq} "
            f"between {before_label} and {after_label}. Scheduling {total} GLY residue(s) (added {added})."
        )

    distance = entry["details"].get("distance")
    threshold = entry["details"].get("threshold")
    est_total = entry["details"].get("estimated_total", total)
    return (
        f"Detected geometric gap in chain {chain_label} between {before_label} and {after_label}: "
        f"Cα–Cα distance {distance:.2f} Å exceeds threshold {threshold:.2f} Å. "
        f"Estimating {est_total} missing residue(s), adding {added} GLY to reach {total}."
    )


def _summarize_insertions(
    chains: List,
    chain_residues: Dict[int, List],
    missing: Dict[Tuple[int, int], Iterable],
    residue_name: str,
) -> Tuple[Dict[Tuple[int, int], List], List[str], int]:
    """
    Convert the missing residue templates to glycine and build log strings.
    """
    glycine_missing: Dict[Tuple[int, int], List[str]] = {}
    summary_lines: List[str] = []
    total = 0

    for key, residues in missing.items():
        res_list = list(residues)

        if not isinstance(key, tuple) or len(key) != 2:
            glycine_missing[key] = [_force_residue(r, residue_name) for r in res_list]
            total += len(res_list)
            summary_lines.append(f"[WARN] Unrecognized missing-residue key {key}; forcing {residue_name.upper()}.")
            continue
        chain_index, insert_index = key
        if chain_index >= len(chains):
            glycine_missing[key] = [_force_residue(r, residue_name) for r in res_list]
            total += len(res_list)
            summary_lines.append(
                f"[WARN] Chain index {chain_index} outside available chains; inserted {len(res_list)} {residue_name.upper()} residue(s)."
            )
            continue
        glycine_missing[key] = [_force_residue(r, residue_name) for r in res_list]
        total += len(res_list)

        chain = chains[chain_index]
        chain_label = _chain_label(chain, chain_index)
        chain_res = chain_residues.get(chain_index, [])

        before = chain_res[insert_index - 1] if insert_index > 0 and insert_index <= len(chain_res) else None
        after = chain_res[insert_index] if insert_index < len(chain_res) else None

        before_label = _residue_label(before) if before is not None else "start"
        after_label = _residue_label(after) if after is not None else "end"

        summary_lines.append(
            f"Chain {chain_label}: inserting {len(res_list)} {residue_name.upper()} residue(s) between {before_label} and {after_label}"
        )

    return glycine_missing, summary_lines, total


def _filter_terminal_gaps(
    missing: Dict[Tuple[int, int], Iterable],
    chain_lengths: Dict[int, int],
) -> Tuple[Dict[Tuple[int, int], List], int]:
    """
    Drop terminal gaps unless explicitly requested.
    """
    filtered: Dict[Tuple[int, int], List] = {}
    skipped = 0

    for key, residues in missing.items():
        res_list = list(residues)
        if not isinstance(key, tuple) or len(key) != 2:
            filtered[key] = res_list
            continue

        chain_index, insert_index = key
        chain_len = chain_lengths.get(chain_index, 0)
        if insert_index == 0 or insert_index == chain_len:
            skipped += len(res_list)
            continue

        filtered[key] = res_list

    return filtered, skipped


def _normalize_missing_entries(missing: Dict[Tuple[int, int], Iterable], residue_name: str):
    for key in list(missing.keys()):
        residues = missing[key]
        missing[key] = [_force_residue(res, residue_name) for res in list(residues)]


def _build_gap_entries(
    missing: Dict[Tuple[int, int], List],
    chains: List,
    chain_residues: Dict[int, List],
    log_lookup: Dict[Tuple[int, int], List[Dict]],
    base_keys: Optional[set],
    residue_name: str,
) -> List[Dict]:
    entries: List[Dict] = []
    base_keys = base_keys or set()

    for key in sorted(missing.keys()):
        chain_index, insert_index = key
        seq_list = [_force_residue(res, residue_name) for res in missing[key]]
        missing[key] = seq_list

        chain_label = str(chain_index)
        chain_res = chain_residues.get(chain_index, [])
        chain_obj = chains[chain_index] if 0 <= chain_index < len(chains) else None
        if chain_obj is not None:
            chain_label = _chain_label(chain_obj, chain_index)

        before = chain_res[insert_index - 1] if insert_index > 0 and insert_index <= len(chain_res) else None
        after = chain_res[insert_index] if insert_index < len(chain_res) else None

        notes = []
        if log_lookup.get(key):
            notes.extend([_format_detection_log(log, chains) for log in log_lookup[key]])
        elif key in base_keys:
            notes.append("Listed by existing sequence metadata.")

        entry = {
            "chain_index": chain_index,
            "chain_id": chain_label,
            "insert_index": insert_index,
            "before_residue": _residue_label(before) or None,
            "after_residue": _residue_label(after) or None,
            "count": len(seq_list),
            "sequence": seq_list,
        }
        if notes:
            entry["notes"] = notes
        entries.append(entry)

    return entries


def _write_gap_report(path: str, entries: List[Dict], source_path: str, residue_name: str):
    payload = {
        "source_structure": source_path,
        "default_residue": residue_name.upper(),
        "instructions": (
            "Edit each gap's 'sequence' list (N->C order) with your preferred residue names. "
            "Save this JSON and feed it to scripts/rebuild_loops_from_config.py to rebuild with the custom sequence."
        ),
        "gaps": entries,
    }
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def rebuild_loops(args: argparse.Namespace):
    fixer = PDBFixer(filename=args.inp)

    if args.detect_only and not args.gap_report:
        print("Detect-only mode requested without --gap-report; no config file will be written.")

    if not args.keep_heterogens:
        fixer.removeHeterogens(keepWater=True)

    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    chain_residues = {idx: list(chain.residues()) for idx, chain in enumerate(chains)}
    chain_lengths = {idx: len(res_list) for idx, res_list in chain_residues.items()}

    missing = {key: list(val) for key, val in (getattr(fixer, "missingResidues", {}) or {}).items()}
    base_missing_keys = set(missing.keys())

    atom_positions = _build_atom_position_map(fixer)
    chain_infos = _collect_chain_infos(chains, chain_residues, atom_positions)

    detection_logs: List[Dict] = []
    existing_counts = {key: len(val) for key, val in missing.items()}

    numbering_additions, numbering_logs = _detect_numbering_gaps(chain_infos, existing_counts)
    if numbering_additions:
        _extend_missing_entries(missing, numbering_additions, args.residue_name)
        for key, count in numbering_additions.items():
            existing_counts[key] = existing_counts.get(key, 0) + count
    detection_logs.extend(numbering_logs)

    distance_additions, distance_logs = _detect_distance_gaps(
        chain_infos,
        existing_counts,
        args.distance_threshold,
        args.ideal_ca_spacing,
    )
    if distance_additions:
        _extend_missing_entries(missing, distance_additions, args.residue_name)
        for key, count in distance_additions.items():
            existing_counts[key] = existing_counts.get(key, 0) + count
    detection_logs.extend(distance_logs)

    log_lookup: Dict[Tuple[int, int], List[Dict]] = {}
    for log in detection_logs:
        key = (log["chain_index"], log["insert_index"])
        log_lookup.setdefault(key, []).append(log)

    if not args.include_termini and missing:
        missing, skipped = _filter_terminal_gaps(missing, chain_lengths)
        if skipped:
            print(f"Skipped {skipped} residue(s) at chain termini (use --include-termini to include them).")

    if detection_logs:
        for log in detection_logs:
            print(_format_detection_log(log, chains))

    _normalize_missing_entries(missing, args.residue_name)
    gap_entries = _build_gap_entries(missing, chains, chain_residues, log_lookup, base_missing_keys, args.residue_name)

    if args.gap_report:
        _write_gap_report(args.gap_report, gap_entries, args.inp, args.residue_name)
        print(f"Wrote gap config to {args.gap_report}")

    if args.detect_only:
        print("Detect-only mode requested; skipping coordinate rebuild.")
        return

    if missing:
        gly_missing, summaries, total = _summarize_insertions(chains, chain_residues, missing, args.residue_name)
        fixer.missingResidues = gly_missing

        if total == 0:
            print("No internal loop gaps required rebuilding.")
        else:
            if hasattr(fixer, "addMissingResidues"):
                fixer.addMissingResidues()
            else:
                print("This PDBFixer build inserts new residues during addMissingAtoms(); skipping explicit addMissingResidues step.")
            for line in summaries:
                print(line)
            print(f"Inserted {total} {args.residue_name.upper()} residue(s) across all chains.")
    else:
        print("No missing residues detected; skipping loop rebuild.")

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    if args.add_hydrogens:
        fixer.addMissingHydrogens(pH=args.ph)

    with open(args.out, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    print(f"Wrote: {args.out}")


def main():
    ap = argparse.ArgumentParser(description="Rebuild missing loops with glycine residues using PDBFixer.")
    ap.add_argument("--in", dest="inp", required=True, help="Input PDB file.")
    ap.add_argument("--out", required=True, help="Output PDB file with glycine bridges.")
    ap.add_argument("--add-hydrogens", action="store_true", help="Add hydrogens after rebuilding loops.")
    ap.add_argument("--ph", type=float, default=7.0, help="pH to use when adding hydrogens (default 7.0).")
    ap.add_argument("--keep-heterogens", action="store_true", help="Keep heterogens (default removes them).")
    ap.add_argument("--include-termini", action="store_true",
                    help="Also rebuild gaps at termini (default only rebuilds internal loops).")
    ap.add_argument("--distance-threshold", type=float, default=6.0,
                    help="Detect geometric gaps when adjacent Cα atoms are farther apart than this (Å). Default 6.0 Å.")
    ap.add_argument("--ideal-ca-spacing", type=float, default=3.8,
                    help="Ideal Cα-Cα spacing (Å) used to estimate how many residues to insert when only distance "
                         "information is available. Default 3.8 Å.")
    ap.add_argument("--residue-name", default="GLY",
                    help="Residue name to use when automatically filling gaps (default: GLY).")
    ap.add_argument("--gap-report",
                    help="Optional JSON path to store detected gaps for manual editing.")
    ap.add_argument("--detect-only", action="store_true",
                    help="Only detect/write the gap report without rebuilding coordinates.")
    args = ap.parse_args()

    if not args.residue_name:
        args.residue_name = "GLY"
    args.residue_name = args.residue_name.strip().upper() or "GLY"

    rebuild_loops(args)


if __name__ == "__main__":
    main()
