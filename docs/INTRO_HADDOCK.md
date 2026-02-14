# HADDOCK Quickstart (Local Pipeline)

This is a concise, battle-tested checklist to get HADDOCK3 runs going in this repo.

## Environment
- Activate local venv:
  ```bash
  source .venv/bin/activate
  ```
- HADDOCK3 is invoked as `.venv/bin/haddock3` by the pipeline; no need to add it to PATH.

## Inputs & References
- Put complex PDBs in `data/` (e.g., `refold_XXXX_Epi.pdb`).
- Receptor chains are inferred from `pdbAG/<base>.pdb`, where `base` = complex name with `refold_` prefix removed and `_Epi` suffix removed. Remaining chains are ligand.
- Do **not** pre-create `results/run_<ID>`; HADDOCK must create/manage its run_dir.

## One-Command Pipeline (recommended)
```bash
python scripts/run_pipeline.py --complexes refold_8IUZ_Epi.pdb
```
What it does per complex:
1) Infer receptor vs ligand chains from `pdbAG/<base>.pdb`.
2) Split complex → `data/receptor_<base>.pdb`, `data/ligand_<base>.pdb`.
3) Merge to single chains → `data/receptor_single_<base>.pdb` (A), `data/ligand_single_<base>.pdb` (X), and merged complex.
4) Build intra-ligand lock restraints `restraints/intra_lig_lock_<base>.tbl` (keeps H/L orientation fixed).
5) Build interface AIRs `restraints/air_<base>.tbl` (receptor hotspots from B-factor > -20 to ligand interface residues within 6 Å).
6) Write config `docs/haddock3_config_<base>.cfg` (ncores=20, sampling=20; ambig_fname=air_*, unambig_fname=lock_*; molecules=merged singles).
7) Run HADDOCK3 with `--restart 0` via `.venv/bin/haddock3`.
8) Archive inputs/config/restraints into `results/run_<base>/` after the run finishes; clean intermediates from `data/`.

Logs:
- Pipeline: `docs/pipeline_logs/<base>.log`
- HADDOCK run log: `docs/pipeline_logs/haddock_<base>.log`
- Summary appended to `docs/LOG.md`

## Manual Run (if needed)
```bash
source .venv/bin/activate
haddock3 docs/haddock3_config_<base>.cfg
```
(Ensure `results/run_<base>` does not pre-exist.)

## Graph Conversion (flexref → interface-only graphs)
After a successful run:
```bash
source .venv/bin/activate
python scripts/flexref_to_graphs.py --runs results/run_<base>/2_flexref
# Outputs: Graphs/flexref/<run_id>.pth
```
Graphs contain only interface nodes/edges (ligand X ↔ receptor non-X within 10 Å).

## Troubleshooting
- If HADDOCK aborts complaining about existing run_dir, remove `results/run_<base>` and rerun, or rely on the pipeline (it cleans the run_dir before running).
- If haddock3 not found, ensure venv is activated (`.venv/bin/haddock3`).
- Missing reference PDB: ensure `pdbAG/<base>.pdb` exists (base = complex name stripped of `refold_` and `_Epi`).

Stay disciplined: run the pipeline, let HADDOCK create its run_dir, and keep `data/` clean.
