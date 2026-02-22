#!/usr/bin/env python3
"""
Apply a manually edited gap configuration to rebuild loops with custom sequences.

First run `rebuild_loops_with_glycine.py --gap-report gaps.json --detect-only` to generate
the configuration file. Edit each gap's `sequence` list to reflect the residues you want
to insert (N->C order), then run this script:

    python scripts/rebuild_loops_from_config.py \\
        --in data/Final_MARTINI_Fan.pdb \\
        --config gaps.json \\
        --out data/Final_MARTINI_Fan_custom_loops.pdb
"""

import argparse
import json
from typing import Dict, List, Tuple

from pdbfixer import PDBFixer
from openmm.app import PDBFile


def _chain_label(chain, index: int) -> str:
    cid = getattr(chain, "id", None)
    if cid is None:
        return str(index)
    cid = str(cid).strip()
    return cid if cid else str(index)


def _residue_label(res) -> str:
    if res is None:
        return ""
    rid = getattr(res, "id", "").strip()
    name = getattr(res, "name", "UNK").strip() or "UNK"
    return f"{name}{rid}" if rid else name


def _normalize_sequence(seq_list: List, count: int, fallback: str) -> List[str]:
    fallback = fallback.upper()
    normalized: List[str] = []
    for entry in seq_list:
        if isinstance(entry, str) and entry.strip():
            normalized.append(entry.strip().upper())
        else:
            normalized.append(fallback)

    if count is None or count <= 0:
        count = len(normalized) or 1

    if len(normalized) < count:
        normalized.extend([fallback] * (count - len(normalized)))
    elif len(normalized) > count:
        normalized = normalized[:count]
    return normalized


def _load_gap_config(path: str, fallback: str) -> Tuple[Dict[Tuple[int, int], List[str]], List[Dict]]:
    with open(path) as f:
        data = json.load(f)

    entries = data.get("gaps", [])
    missing: Dict[Tuple[int, int], List[str]] = {}
    for entry in entries:
        chain_index = entry.get("chain_index")
        insert_index = entry.get("insert_index")
        if chain_index is None or insert_index is None:
            continue

        try:
            chain_index = int(chain_index)
            insert_index = int(insert_index)
        except ValueError:
            continue

        count = entry.get("count")
        if isinstance(count, str) and count.isdigit():
            count = int(count)
        elif not isinstance(count, int):
            count = None

        seq_list = entry.get("sequence", [])
        normalized = _normalize_sequence(list(seq_list), count, fallback)
        if not normalized:
            continue

        missing[(chain_index, insert_index)] = normalized

    return missing, entries


def _summarize_config(missing: Dict[Tuple[int, int], List[str]], chains, chain_residues):
    for key in sorted(missing.keys()):
        chain_index, insert_index = key
        seq_list = missing[key]
        chain = chains[chain_index] if 0 <= chain_index < len(chains) else None
        chain_label = _chain_label(chain, chain_index) if chain is not None else str(chain_index)
        chain_res = chain_residues.get(chain_index, [])
        before = chain_res[insert_index - 1] if insert_index > 0 and insert_index <= len(chain_res) else None
        after = chain_res[insert_index] if insert_index < len(chain_res) else None
        before_label = _residue_label(before) if before is not None else "start"
        after_label = _residue_label(after) if after is not None else "end"
        print(
            f"Config: chain {chain_label}: inserting {len(seq_list)} residue(s) "
            f"between {before_label} and {after_label} (sequence: {', '.join(seq_list)})"
        )


def rebuild_from_config(args: argparse.Namespace):
    missing, entries = _load_gap_config(args.config, args.fallback_residue)
    if not missing:
        raise SystemExit("No usable gap entries were found in the config file.")

    fixer = PDBFixer(filename=args.inp)

    if not args.keep_heterogens:
        fixer.removeHeterogens(keepWater=True)

    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    chain_residues = {idx: list(chain.residues()) for idx, chain in enumerate(chains)}

    _summarize_config(missing, chains, chain_residues)
    fixer.missingResidues = missing

    if hasattr(fixer, "addMissingResidues"):
        fixer.addMissingResidues()
    else:
        print("This PDBFixer build inserts new residues during addMissingAtoms(); skipping explicit addMissingResidues step.")

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    if args.add_hydrogens:
        fixer.addMissingHydrogens(pH=args.ph)

    with open(args.out, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    print(f"Wrote: {args.out}")


def main():
    ap = argparse.ArgumentParser(description="Rebuild loops from a JSON gap configuration.")
    ap.add_argument("--in", dest="inp", required=True, help="Input PDB file.")
    ap.add_argument("--config", required=True, help="Gap configuration JSON (from rebuild_loops_with_glycine.py).")
    ap.add_argument("--out", required=True, help="Output rebuilt PDB file.")
    ap.add_argument("--fallback-residue", default="GLY",
                    help="Fallback residue name if an entry omits a sequence (default: GLY).")
    ap.add_argument("--add-hydrogens", action="store_true", help="Add hydrogens after rebuilding loops.")
    ap.add_argument("--ph", type=float, default=7.0, help="pH to use when adding hydrogens (default 7.0).")
    ap.add_argument("--keep-heterogens", action="store_true", help="Keep heterogens (default removes them).")
    args = ap.parse_args()

    args.fallback_residue = (args.fallback_residue or "GLY").strip().upper() or "GLY"
    rebuild_from_config(args)


if __name__ == "__main__":
    main()
