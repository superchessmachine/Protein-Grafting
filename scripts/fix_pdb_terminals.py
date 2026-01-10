#!/usr/bin/env python3
"""
Rebuild missing N- and C-terminal residues using PDBFixer.
"""

import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input PDB")
    ap.add_argument("--out", dest="out", required=True, help="Output PDB")
    ap.add_argument("--ph", type=float, default=None, help="Add hydrogens at this pH")
    args = ap.parse_args()

    fixer = PDBFixer(filename=args.inp)

    # Identify problems
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()

    # Convert things like MSE → MET
    if fixer.nonstandardResidues:
        fixer.replaceNonstandardResidues()

    # Remove heterogens except water (recommended for rebuilding termini)
    fixer.removeHeterogens(keepWater=True)

    # Re-run detection after cleanup
    fixer.findMissingResidues()
    fixer.findMissingAtoms()

    # 🔑 This step rebuilds BOTH missing atoms AND missing residues
    fixer.addMissingAtoms()

    # Optional hydrogens
    if args.ph is not None:
        fixer.addMissingHydrogens(pH=args.ph)

    # Write result
    with open(args.out, "w") as f:
        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            f,
            keepIds=True
        )

    print(f"Fixed PDB written to: {args.out}")

if __name__ == "__main__":
    main()
