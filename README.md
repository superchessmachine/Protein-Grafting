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
5. **Fill unresolved residues from a complete model (optional)**
   - If the experimental structure and complete model share the same sequence register, use `scripts/graft_missing_residues.py` to retain experimental coordinates wherever present and copy only unresolved residue positions from the complete model. This path can also drop an N-terminal construct prefix and retain selected experimental heterogens such as metal ions.
6. **Rebuild unresolved loops (optional)**
   - If the grafted model is missing internal loops, run `scripts/rebuild_loops_with_glycine.py` to detect gaps by sequence numbering and Cα spacing. Use `--gap-report loops.json --detect-only` to review the proposed insertions, edit the JSON to the residue identities you want, then call `scripts/rebuild_loops_from_config.py` to rebuild the loops (and optionally add hydrogens) with your edited sequence.

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

# 5) Fill unresolved residue positions from a complete same-register model
python scripts/graft_missing_residues.py \
  --experimental experimental_structure.cif \
  --model complete_model.cif \
  --out hybrid_model.pdb \
  --report hybrid_model_report.json \
  --drop-prefix 0 \
  --keep-heterogen ZN

# 6) Detect/rebuild missing loops with glycine and emit an editable config
python scripts/rebuild_loops_with_glycine.py \
  --in grafted_model.pdb \
  --out grafted_model_glycine_loops.pdb \
  --gap-report loops_config.json \
  --residue-name GLY \
  --detect-only

# 7) Rebuild loops from the edited config with custom sequences
python scripts/rebuild_loops_from_config.py \
  --in grafted_model_glycine_loops.pdb \
  --config loops_config.json \
  --out grafted_model_custom_loops.pdb \
  --add-hydrogens \
  --ph 7.4
```

`--extra-prefix-length` reflects how many residues exist at the N-terminus of the experimental construct but not in the AF3 model (set it to 0 if there is no insertion). The script assumes that once those residues are removed, `fasta1` is a prefix of `fasta2` and that both PDBs are already in the same coordinate frame.

`graft_missing_residues.py` is useful when you want the experimental structure to remain the source of truth. It matches by polymer sequence position, optionally superposes the complete model onto the resolved experimental residues, drops any sequence prefix you specify, retains experimental coordinates wherever they exist, and fills only absent polymer residues from the model. For numbered PDB inputs where unresolved experimental residues are completely absent from ATOM records, add `--use-model-sequence` so the complete model defines the output residue list. For difficult insertions, repeat `--model-range START-END` and set `--local-align-flank N` to use model coordinates for a bounded segment aligned to nearby experimental anchors. The JSON report lists filled ranges, overwritten experimental ranges, retained heterogens, alignment RMSD, local range fits, sequence source, and peptide-bond junction warnings so you can inspect any geometrically difficult insertions before using the hybrid model downstream.

## Loop Rebuilding Utilities
- `scripts/rebuild_loops_with_glycine.py` drives the automatic detection step. PDBFixer-supplied metadata is combined with residue numbering gaps and Cα–Cα distances to estimate how many residues should sit between two coordinates. Everything is filled with the residue you pass via `--residue-name` (GLY by default), and heterogens are removed unless `--keep-heterogens` is set. Add `--include-termini` to patch missing residues at the start or end of a chain. When you supply `--gap-report`, the script stores every detected gap (with the surrounding residues, chain index, insert position, and notes explaining which heuristic flagged it). Passing `--detect-only` stops after writing the JSON so you can review it before touching coordinates.
- `scripts/rebuild_loops_from_config.py` rebuilds the loops using an edited JSON report. For each `sequence` list in the config file, replace the placeholder glycine strings with the residue names you want in N→C order. The script replays the insertions with PDBFixer, optionally adds hydrogens (`--add-hydrogens --ph 7.0`), and retains heterogens when `--keep-heterogens` is supplied. Use `--fallback-residue` if any entry is left blank; those positions will default to that residue.

## Tips
- Treat the scripts as building blocks—you can skip or reorder steps if your structures already satisfy certain requirements.
- Keep raw downloads separate from this repo to avoid accidentally committing proprietary data.
- When in doubt, inspect intermediate files (especially the FASTA files) before grafting to ensure the sequences agree.
