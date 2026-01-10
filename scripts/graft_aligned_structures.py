#!/usr/bin/env python3
"""
Graft experimental coordinates from an aligned template onto a full-length AF3/AlphaFold model,
producing a single clean chain PDB.

Assumptions:
- The experimental subset PDB (pdb1) and AF3/AlphaFold model (pdb2) are already aligned in PyMOL.
- The experimental FASTA (fasta1) is a prefix of the AF3 FASTA (fasta2) once you remove any inserted
  N-terminal residues that only exist in the experimental construct.

Behavior:
- Output sequence = (inserted experimental prefix) + fasta2
- The output contains a single chain 'A' renumbered sequentially starting at residue 1
- Coordinates/occ/B factors are overwritten by experimental values when atoms are available, otherwise AF3 values are kept
"""

import argparse
from dataclasses import dataclass
from typing import Dict, List, Tuple


@dataclass
class AtomRec:
    serial: int
    name: str          # atom name (4 chars stripped)
    altloc: str        # altLoc
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occ: float
    b: float
    element: str
    charge: str
    record: str        # "ATOM  " or "HETATM"


def parse_pdb_atoms(path: str, keep_hetatm: bool = False) -> List[AtomRec]:
    """
    Minimal PDB ATOM/HETATM parser. Keeps only one altloc per atom:
    prefers altLoc ' ' then 'A'. Drops other altLocs.
    """
    atoms: List[AtomRec] = []
    # For altloc filtering: per (chain, resseq, icode, atomname) keep best
    best_alt_rank = {" ": 0, "A": 1}
    seen_best: Dict[Tuple[str, int, str, str], str] = {}

    with open(path, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or (keep_hetatm and line.startswith("HETATM"))):
                continue

            record = line[0:6]
            serial = int(line[6:11])
            name = line[12:16].strip()
            altloc = line[16:17]
            resname = line[17:20].strip()
            chain = line[21:22].strip() or "A"
            resseq = int(line[22:26])
            icode = line[26:27]

            # coords/occ/b can be blank in some junk PDBs; assume valid here
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            occ_str = line[54:60].strip()
            b_str = line[60:66].strip()
            occ = float(occ_str) if occ_str else 1.00
            b = float(b_str) if b_str else 0.00

            element = line[76:78].strip()
            charge = line[78:80].strip()

            key = (chain, resseq, icode, name)
            if altloc not in best_alt_rank:
                # ignore weird altlocs by default
                continue

            if key in seen_best:
                prev = seen_best[key]
                # keep the better altloc (rank lower is better)
                if best_alt_rank[altloc] < best_alt_rank.get(prev, 999):
                    # replace: remove old atom record for that key
                    # (rare; but keep logic correct)
                    atoms = [a for a in atoms if not (a.chain == chain and a.resseq == resseq and a.icode == icode and a.name == name)]
                    seen_best[key] = altloc
                else:
                    continue
            else:
                seen_best[key] = altloc

            atoms.append(AtomRec(
                serial=serial, name=name, altloc=altloc, resname=resname, chain=chain,
                resseq=resseq, icode=icode, x=x, y=y, z=z, occ=occ, b=b,
                element=element, charge=charge, record=record.strip()
            ))

    return atoms


def residues_in_order(atoms: List[AtomRec], chain_id: str = "A") -> List[Tuple[int, str]]:
    """
    Return unique residues in file order for a given chain:
    list of (resseq, icode) pairs.
    """
    seen = set()
    order = []
    for a in atoms:
        if a.chain != chain_id:
            continue
        key = (a.resseq, a.icode)
        if key not in seen:
            seen.add(key)
            order.append(key)
    return order


def atoms_by_residue(atoms: List[AtomRec], chain_id: str = "A") -> Dict[Tuple[int, str], List[AtomRec]]:
    d: Dict[Tuple[int, str], List[AtomRec]] = {}
    for a in atoms:
        if a.chain != chain_id:
            continue
        key = (a.resseq, a.icode)
        d.setdefault(key, []).append(a)
    return d


def read_fasta_sequence(path: str) -> str:
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq).strip()


def format_pdb_atom_line(serial: int, atom: AtomRec, chain: str, resseq: int, icode: str,
                         resname: str, x: float, y: float, z: float, occ: float, b: float,
                         element: str, charge: str) -> str:
    """
    Write a standard ATOM line. We always output ATOM (not HETATM) for protein.
    """
    # Atom name formatting: right-justify in 4 columns as PDB expects
    atom_name = atom.name
    atom_name_fmt = atom_name.rjust(4)

    altloc = " "  # force blank altloc in output
    chain_fmt = chain if chain else "A"
    icode_fmt = icode if icode else " "

    element_fmt = (element or "").rjust(2)
    charge_fmt = (charge or "").rjust(2)

    return (
        f"ATOM  {serial:5d} {atom_name_fmt}{altloc}{resname:>3s} {chain_fmt}"
        f"{resseq:4d}{icode_fmt}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occ:6.2f}{b:6.2f}          "
        f"{element_fmt}{charge_fmt}\n"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb1", default="experimental_aligned_subset.pdb",
                    help="Experimental aligned PDB (subset + optional extra N-terminal residues)")
    ap.add_argument("--pdb2", default="full_length_af3_model.pdb",
                    help="Full-length AlphaFold3/AlphaFold PDB in the same frame as pdb1")
    ap.add_argument("--fasta1", default="experimental_subset.fasta", help="FASTA for pdb1")
    ap.add_argument("--fasta2", default="full_length_model.fasta", help="FASTA for pdb2")
    ap.add_argument("--out", default="grafted_model.pdb", help="Output PDB path")
    ap.add_argument("--extra-prefix-length", type=int, default=1,
                    help="Number of N-terminal residues that only exist in the experimental structure (default: 1).")
    args = ap.parse_args()

    seq1 = read_fasta_sequence(args.fasta1)
    seq2 = read_fasta_sequence(args.fasta2)

    if not seq1 or not seq2:
        raise SystemExit("Could not read sequences from fasta1/fasta2.")

    insert_len = args.extra_prefix_length
    if insert_len < 0:
        raise SystemExit("--extra-prefix-length must be >= 0.")
    if len(seq1) <= insert_len:
        raise SystemExit("fasta1 must contain more residues than the inserted prefix.")

    core1 = seq1[insert_len:]  # experimental sequence without inserted residues
    if not seq2.startswith(core1):
        raise SystemExit("Expected fasta1 (excluding the inserted prefix) to be a prefix of fasta2.")

    # Load atoms
    atoms1 = parse_pdb_atoms(args.pdb1, keep_hetatm=False)
    atoms2 = parse_pdb_atoms(args.pdb2, keep_hetatm=False)

    # Assume single chain A for both inputs; if your chain IDs differ, adjust here.
    chain1 = "A"
    chain2 = "A"

    res_order1 = residues_in_order(atoms1, chain1)
    res_order2 = residues_in_order(atoms2, chain2)

    # Basic sanity checks vs fasta lengths (not fatal if off due to missing coordinates)
    if len(res_order2) < len(seq2) - 5:
        print(f"[warn] pdb2 has fewer coordinate residues ({len(res_order2)}) than fasta2 length ({len(seq2)}).")
    if len(res_order1) < len(seq1) - 5:
        print(f"[warn] pdb1 has fewer coordinate residues ({len(res_order1)}) than fasta1 length ({len(seq1)}).")

    byres1 = atoms_by_residue(atoms1, chain1)
    byres2 = atoms_by_residue(atoms2, chain2)

    # Build mapping of sequence positions -> residue keys in each structure, using file order.
    # We'll use the residue order list positions as "sequence positions" for coordinate-bearing residues.
    # This script assumes pdb1 covers the N-terminus of pdb2 (after removing the inserted prefix),
    # so residues follow the same order once numbering is adjusted.
    n1 = len(res_order1)
    n2 = len(res_order2)

    if n1 < 1:
        raise SystemExit("pdb1 has no residues in chain A.")
    if n2 < 1:
        raise SystemExit("pdb2 has no residues in chain A.")
    if insert_len > n1:
        raise SystemExit("pdb1 does not contain enough residues for the specified inserted prefix.")

    # Determine overlap length in residues based on fasta1 (coordinate-bearing subset expected)
    overlap_out_end = len(seq1)  # residues 1..overlap_out_end are experimental-covered (inserted prefix included)

    # Create output atom records by iterating over pdb2 atoms as template, shifting resseq by insert_len,
    # then overwriting coords with pdb1 where applicable. Inserted residues are emitted first.
    out_lines: List[str] = []
    atom_serial = 1

    # ---- Write inserted prefix residues from pdb1 (if any) ----
    for idx in range(insert_len):
        ins_key = res_order1[idx]
        ins_atoms = byres1.get(ins_key, [])
        if not ins_atoms:
            raise SystemExit(f"Could not find atoms for inserted residue index {idx+1} in pdb1.")

        ins_resname = ins_atoms[0].resname
        out_resseq = idx + 1
        for a in ins_atoms:
            # element fallback: if blank, infer from atom name first letter
            elem = a.element.strip() if a.element else (a.name[0].upper() if a.name else "")
            out_lines.append(format_pdb_atom_line(
                serial=atom_serial,
                atom=a,
                chain="A",
                resseq=out_resseq,
                icode=" ",
                resname=ins_resname,
                x=a.x, y=a.y, z=a.z,
                occ=a.occ, b=a.b,
                element=elem,
                charge=a.charge
            ))
            atom_serial += 1

    # ---- Prepare quick lookup for experimental atoms for residues after the inserted prefix ----
    # Map output residue index -> dict(atom_name -> AtomRec from pdb1)
    exp_atom_lookup: Dict[int, Dict[str, AtomRec]] = {}

    # pdb1 residue index i (1-based in its own residue order)
    # output residue index = i  (inserted residues keep numbering aligned)
    for i1 in range(insert_len + 1, min(overlap_out_end, n1) + 1):
        key1 = res_order1[i1 - 1]
        atoms_res1 = byres1.get(key1, [])
        if atoms_res1:
            exp_atom_lookup[i1] = {a.name: a for a in atoms_res1}

    # ---- Now write pdb2 atoms, shifted by insert_len in residue numbering, with overrides ----
    # We iterate residue-by-residue in pdb2 order to keep output clean and grouped.
    for i2, key2 in enumerate(res_order2, start=1):  # i2 is 1-based position in pdb2 residue order
        out_res_index = i2 + insert_len  # shift to accommodate inserted prefix

        atoms_res2 = byres2.get(key2, [])
        if not atoms_res2:
            continue

        resname2 = atoms_res2[0].resname

        # If this output residue is within experimental-covered region (2..overlap_out_end),
        # overwrite coordinates for atoms that exist in exp lookup.
        exp_atoms_for_res = exp_atom_lookup.get(out_res_index)

        for a2 in atoms_res2:
            x, y, z = a2.x, a2.y, a2.z
            occ, b = a2.occ, a2.b

            if exp_atoms_for_res is not None:
                a1 = exp_atoms_for_res.get(a2.name)
                if a1 is not None:
                    x, y, z = a1.x, a1.y, a1.z
                    occ, b = a1.occ, a1.b  # take experimental occ/B too

            elem = a2.element.strip() if a2.element else (a2.name[0].upper() if a2.name else "")
            out_lines.append(format_pdb_atom_line(
                serial=atom_serial,
                atom=a2,
                chain="A",
                resseq=out_res_index,
                icode=" ",
                resname=resname2,
                x=x, y=y, z=z,
                occ=occ, b=b,
                element=elem,
                charge=a2.charge
            ))
            atom_serial += 1

    # Add TER and END
    out_lines.append("TER\n")
    out_lines.append("END\n")

    with open(args.out, "w") as f:
        f.writelines(out_lines)

    print(f"Wrote grafted structure: {args.out}")
    print(f"Output residues: 1..{len(seq2)+insert_len} (sequence = inserted prefix + fasta2)")
    print(f"Experimental override region in output: residues 1..{overlap_out_end}")
    print("Chain: A (single chain), atom serials/residue numbers renumbered cleanly.")


if __name__ == "__main__":
    main()
