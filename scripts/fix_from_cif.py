#!/usr/bin/env python3
"""
Fix a structure from mmCIF (.cif/.mmcif) or PDB using PDBFixer.

Key points:
- For most PDBFixer builds, mmCIF is supported via PDBFixer(filename=...)
- Missing residues are inferred from mmCIF polymer sequence metadata when present.
- addMissingAtoms() is the method that actually builds what can be built.

Usage:
  python fix_from_cif.py --in input_structure.cif --out prepared_structure.pdb --out-format pdb --ph 7.4
  python fix_from_cif.py --in input_structure.cif --out prepared_structure.cif --out-format cif
"""

import argparse
import os
import inspect
from pdbfixer import PDBFixer
from openmm.app import PDBFile, PDBxFile

def _make_fixer(input_path: str) -> PDBFixer:
    """
    Create a PDBFixer object from an input structure, handling API differences.
    """
    # 1) Most versions: filename=... works for both .pdb and .cif/.mmcif
    try:
        return PDBFixer(filename=input_path)
    except Exception as e_filename:
        # 2) Some builds support passing a PDBxFile object, but the kw name varies.
        # We'll inspect the signature and try supported options.
        ext = os.path.splitext(input_path.lower())[1]
        if ext in [".cif", ".mmcif"]:
            cif = PDBxFile(input_path)

            sig = inspect.signature(PDBFixer.__init__)
            params = sig.parameters

            # Common possibilities across forks/builds
            for kw in ("pdbxfile", "pdbxFile", "pdbfile", "pdbFile"):
                if kw in params:
                    try:
                        return PDBFixer(**{kw: cif})
                    except Exception:
                        pass

        # If we got here, we truly can't load it with this installation.
        raise RuntimeError(
            f"Could not load '{input_path}' with this PDBFixer build.\n"
            f"First failure (filename=...): {repr(e_filename)}\n"
            f"Try updating:  conda install -c conda-forge pdbfixer openmm\n"
            f"Or provide a PDB instead."
        )

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input .cif/.mmcif or .pdb")
    ap.add_argument("--out", dest="out", required=True, help="Output file path")
    ap.add_argument("--out-format", choices=["pdb", "cif"], default="pdb",
                    help="Output format (default: pdb)")
    ap.add_argument("--ph", type=float, default=None,
                    help="If set, add hydrogens at this pH (e.g. 7.4)")
    ap.add_argument("--keep-heterogens", action="store_true",
                    help="Keep ligands/ions. Default removes heterogens except water (safer for rebuilding).")
    args = ap.parse_args()

    fixer = _make_fixer(args.inp)

    # Detect issues
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        fixer.replaceNonstandardResidues()

    # Heterogens: default remove (except water) to avoid clashes with rebuilt termini
    if not args.keep_heterogens:
        fixer.removeHeterogens(keepWater=True)

    # Missing residues/atoms
    fixer.findMissingResidues()
    fixer.findMissingAtoms()

    # NOTE: This is the actual "build" step in released PDBFixer
    fixer.addMissingAtoms()

    if args.ph is not None:
        fixer.addMissingHydrogens(pH=args.ph)

    # Write output
    if args.out_format == "pdb":
        with open(args.out, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    else:
        with open(args.out, "w") as f:
            PDBxFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    print(f"Wrote: {args.out}")

if __name__ == "__main__":
    main()
