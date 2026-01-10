#!/usr/bin/env python3
"""
Extract protein sequence(s) from a PDB and write FASTA using OpenMM.

Usage:
  python pdb_to_fasta_openmm.py --in input.pdb --out output.fasta
"""

import argparse
from collections import defaultdict
from openmm.app import PDBFile

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    # common variants you might want mapped:
    "MSE":"M",  # selenomethionine -> methionine
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input PDB path")
    ap.add_argument("--out", dest="out", required=True, help="Output FASTA path")
    ap.add_argument("--unknown-as-x", action="store_true",
                    help="Write unknown residues as X instead of skipping them.")
    args = ap.parse_args()

    pdb = PDBFile(args.inp)
    top = pdb.topology

    seq_by_chain = defaultdict(list)

    for chain in top.chains():
        chain_id = chain.id if chain.id not in (None, "") else "A"
        for res in chain.residues():
            name = res.name.strip().upper()
            one = AA3_TO_1.get(name, None)
            if one is None:
                if args.unknown_as_x:
                    one = "X"
                else:
                    continue
            seq_by_chain[chain_id].append(one)

    with open(args.out, "w") as f:
        for chain_id in sorted(seq_by_chain.keys()):
            seq = "".join(seq_by_chain[chain_id])
            f.write(f">protein_chain_{chain_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

    print(f"Wrote FASTA to: {args.out}")
    for chain_id in sorted(seq_by_chain.keys()):
        print(f"Chain {chain_id}: {len(seq_by_chain[chain_id])} residues")

if __name__ == "__main__":
    main()
