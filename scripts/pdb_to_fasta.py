#!/usr/bin/env python3
"""
Extract protein sequence(s) from a PDB and write FASTA.
Uses MDAnalysis and only standard residues.

Usage:
  python pdb_to_fasta_mda.py --in input.pdb --out output.fasta
"""

import argparse
from collections import defaultdict

import MDAnalysis as mda
from MDAnalysis.lib.util import convert_aa_code

def resname_to_oneletter(resname: str) -> str:
    """
    Convert 3-letter residue name -> 1-letter.
    Returns 'X' if unknown.
    """
    try:
        one = convert_aa_code(resname, "one")
        if len(one) == 1:
            return one
        return "X"
    except Exception:
        return "X"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input PDB path")
    ap.add_argument("--out", dest="out", required=True, help="Output FASTA path")
    ap.add_argument("--include-nonstd", action="store_true",
                    help="Include non-standard residues as 'X' rather than skipping them.")
    args = ap.parse_args()

    u = mda.Universe(args.inp)

    # Protein residues only (MDAnalysis knows standard protein residues)
    protein = u.select_atoms("protein")
    if protein.n_atoms == 0:
        raise SystemExit("No protein atoms found (selection 'protein' is empty).")

    # Build sequences per chain (segid if present, else chainID)
    # MDAnalysis stores chain IDs in residue.segid or residue.chainID depending on parser/version.
    seq_by_chain = defaultdict(list)

    for res in protein.residues:
        # Determine chain label
        chain = None
        # Try common attributes in a sane order
        if hasattr(res, "segid") and str(res.segid).strip():
            chain = str(res.segid).strip()
        elif hasattr(res, "chainID") and str(res.chainID).strip():
            chain = str(res.chainID).strip()
        else:
            chain = "A"  # fallback

        one = resname_to_oneletter(res.resname)

        if one == "X" and not args.include_nonstd:
            # Skip unknown/nonstandard residues by default
            continue

        seq_by_chain[chain].append(one)

    # Write FASTA
    with open(args.out, "w") as f:
        for chain in sorted(seq_by_chain.keys()):
            seq = "".join(seq_by_chain[chain])
            header = f">protein_chain_{chain}"
            f.write(header + "\n")
            # wrap at 80 chars
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

    print(f"Wrote FASTA to: {args.out}")
    for chain in sorted(seq_by_chain.keys()):
        print(f"Chain {chain}: {len(seq_by_chain[chain])} residues (after filtering)")

if __name__ == "__main__":
    main()
