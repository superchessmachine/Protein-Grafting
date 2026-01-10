# AF3 Structural Grafting Workflow

This repository collects the small utilities I use to stitch experimental coordinates onto an AlphaFold3/AlphaFold (AF3) model after aligning the structures in PyMOL. Each script is standalone so you can mix and match them inside any project, but keeping them together with a short playbook makes it easier to reproduce the full procedure.

## Repository Layout
- `scripts/` – all helper scripts (structure cleanup, FASTA extraction, grafting)
- `requirements.txt` – Python dependencies for the scripts

## Environment Setup
1. Use Python 3.9+ in a virtual environment or conda environment.
2. Install the required packages: `pip install -r requirements.txt`
3. Install PyMOL separately (not included here) for the manual structural alignment step.

## Workflow Outline
The logical order below mirrors the original "step1/step2/step3" folders while replacing the protein-specific files with generic placeholders.

1. **Prepare the experimental structure**
   - If you pull a structure from the PDB as an mmCIF file, convert and clean it with `scripts/fix_from_cif.py`. This fills in missing atoms, removes heterogens (unless you keep them), and optionally adds hydrogens.
   - If you already have a PDB, run `scripts/fix_pdb_terminals.py` to rebuild missing N- and C-terminal residues before anything else touches the file.
2. **Generate FASTA sequences**
   - Use `scripts/pdb_to_fasta.py` (MDAnalysis backend) or `scripts/pdb_to_fasta_pdbfixer.py` (OpenMM backend) to extract sequences from both the cleaned experimental structure and the AF3 model. The choice depends on the libraries you have installed.
3. **Align structures in PyMOL**
   - Load the experimental PDB and the AF3 model.
   - Trim/align so that the experimental construct exactly matches the region that should replace the AF3 coordinates. Save the aligned experimental subset as `experimental_aligned_subset.pdb` and the aligned AF3 model as `full_length_af3_model.pdb` (names can be anything; they just need to match the CLI arguments).
4. **Graft experimental coordinates onto the AF3 model**
   - Run `scripts/graft_aligned_structures.py` with the aligned PDBs plus the FASTA files you generated above. This script inserts any extra N-terminal residues from the experimental construct and overwrites the AF3 coordinates for the overlapping region using the experimental atoms. All residues are renumbered sequentially for a tidy output.

## Script Reference
Below are the most common invocations; adjust paths and options as needed. None of the examples contain identifying information—replace placeholders with your own filenames.

```bash
# 1) Clean an mmCIF file and export a PDB
python scripts/fix_from_cif.py --in input_structure.cif --out prepared_structure.pdb --out-format pdb --ph 7.4

# 2) Rebuild termini in an existing PDB
python scripts/fix_pdb_terminals.py --in prepared_structure.pdb --out prepared_structure_fixed.pdb --ph 7.0

# 3a) Extract FASTA using MDAnalysis
python scripts/pdb_to_fasta.py --in prepared_structure_fixed.pdb --out experimental_subset.fasta

# 3b) Extract FASTA using OpenMM (alternative backend)
python scripts/pdb_to_fasta_pdbfixer.py --in full_length_af3_model.pdb --out full_length_model.fasta --unknown-as-x

# 4) Graft the aligned structures
python scripts/graft_aligned_structures.py \
  --pdb1 experimental_aligned_subset.pdb \
  --pdb2 full_length_af3_model.pdb \
  --fasta1 experimental_subset.fasta \
  --fasta2 full_length_model.fasta \
  --out grafted_model.pdb \
  --extra-prefix-length 1
```

`--extra-prefix-length` reflects how many residues exist at the N-terminus of the experimental construct but not in the AF3 model (set it to 0 if there is no insertion). The script assumes that once those residues are removed, `fasta1` is a prefix of `fasta2` and that both PDBs are already in the same coordinate frame.

## Tips
- Treat the scripts as building blocks—you can skip or reorder steps if your structures already satisfy certain requirements.
- Keep raw downloads separate from this repo to avoid accidentally committing proprietary data.
- When in doubt, inspect intermediate files (especially the FASTA files) before grafting to ensure the sequences agree.
# Protein-Grafting
