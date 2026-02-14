# HADDOCK_test — Guided Docking + Batch Pipeline

This repo is a scaffold plus a batch pipeline to prepare, run, and post-process guided HADDOCK3 docking jobs. It handles splitting/merging chains, building restraints, running HADDOCK, archiving per-run inputs, and cleaning intermediates. Graph conversion scripts are included for flexref outputs.

## Activation
Create + activate the local venv:
```bash
python3 -m venv .venv
source .venv/bin/activate
```

## Key Paths
- `data/` — input PDBs. **Per-run intermediates are cleaned after pipeline runs.** Example input included: `data/example_input.pdb` (antibody–antigen complex with epitope hotspots in B-factors).
- `restraints/` — AIRs, intra-ligand locks, other restraint files.
- `scripts/` — helpers: split/merge, restraints, pipeline, graph conversion.
- `results/` — per-run archives (`results/run_<ID>/`). Do **not** pre-create/seed these before HADDOCK runs; HADDOCK needs to create/manage its run_dir.
- `docs/` — configs/logs; `docs/pipeline_logs/` per-run logs.
- `Graphs/` — graph outputs (e.g., flexref interface graphs).
- `pdbAG/` — receptor reference PDBs used to infer receptor chains (strip `refold_` prefix and `_Epi` suffix from complex name).

## Prepare & run (full example)
Prereqs: activate `.venv` and have `haddock3` installed inside it.

1) Put complexes in `data/`. Example provided: `data/example_input.pdb`.
2) Ensure a receptor reference exists at `pdbAG/<base>.pdb` (base = filename stem). For the example, copy/link your receptor reference to `pdbAG/example_input.pdb`.
3) Hotspots: encode receptor hotspots in B-factors (> -20) of the reference PDB; pipeline builds AIRs to ligand residues within 6 Å.
4) Run the pipeline on the example:
```bash
python scripts/run_pipeline.py --complexes example_input.pdb
```
5) Inspect outputs in `results/run_example_input/` and logs in `docs/pipeline_logs/`.

Batch multiple inputs:
```bash
python scripts/run_pipeline.py --complexes example_input.pdb another.pdb third.pdb
```

## Batch Pipeline (recommended)
Pipeline steps per complex:
1) Infer receptor chains from `pdbAG/<base>.pdb` (base = complex stem without `refold_` and `_Epi`). Remaining chains = ligand.
2) Split complex → `data/receptor_<base>.pdb`, `data/ligand_<base>.pdb`.
3) Merge to single chains → `data/receptor_single_<base>.pdb` (A), `data/ligand_single_<base>.pdb` (X), and merged complex `data/complex_merged_<base>.pdb`.
4) Build intra-ligand lock restraints `restraints/intra_lig_lock_<base>.tbl` (keeps heavy/light orientation fixed).
5) Build interface AIRs `restraints/air_<base>.tbl` (receptor hotspots from B-factor > -20 to ligand interface residues within 6 Å).
6) Write config `docs/haddock3_config_<base>.cfg` (ncores=20, sampling=20; ambig_fname=air_*, unambig_fname=lock_*; molecules = merged singles).
7) Run HADDOCK3 (`.venv/bin/haddock3 --restart 0 ...`); archive splits/merges/restraints/config into `results/run_<base>/` **after** HADDOCK finishes. Do not pre-create `results/run_<base>` before the run.
8) Clean intermediates from `data/` (split/merged PDBs, merged complex).
Logs: `docs/pipeline_logs/<base>.log` and appended to `docs/LOG.md`. Run log: `docs/pipeline_logs/haddock_<base>.log`.

## Flow (plain text)
- Input PDBs in `data/` (e.g., `data/example_input.pdb`).
- Infer receptor from `pdbAG/<base>.pdb` (base = filename stem, strip `refold_`/`_Epi`). Remaining chains = ligand.
- Split complex → receptor/ligand; merge to single-chain receptor (A) + ligand (X).
- Build restraints: intra-ligand lock + AIRs from receptor hotspots (B-factor > -20) to ligand residues within 6 Å.
- Write `docs/haddock3_config_<base>.cfg`.
- Run HADDOCK3.
- Archive `results/run_<base>/` inputs + logs; clean `data/` intermediates.

### How to prepare PDB inputs
- Place complexes in `data/` (example: `data/example_input.pdb`). Filename stem = `<base>`.
- Place receptor reference at `pdbAG/<base>.pdb` (strip any `refold_` prefix and `_Epi` suffix). Chain IDs in `pdbAG/<base>.pdb` define receptor; remaining chains in the complex become ligand.
- Hotspots: set receptor residue B-factors > -20 for residues to attract ligand; pipeline builds AIRs to ligand residues within 6 Å.
- Ligand should be present in the complex PDB. If antibody H/L are separate chains, pipeline will merge to a single chain X and add intra-ligand lock restraints to keep orientation.

## Helper Scripts
- `scripts/prepare_from_complex.py` — single split/merge using receptor inferred from `pdbAG/<name>.pdb`.
- `scripts/merge_chains.py` — merge chains to one chain, renumber residues.
- `scripts/make_air.py` — build AIRs from residue lists.
- `scripts/merge_two_pdbs.py` — merge receptor_single + ligand_single into one PDB.
- `scripts/flexref_to_graphs.py` — convert flexref PDBs (`results/<ID>/2_flexref/*.pdb`) to **interface-only** PyG graphs (ligand X ↔ receptor non-X within 10 Å). Output: `Graphs/flexref/<ID>.pth`. **Dependencies:** `torch`, `torch_geometric`.
  - Ligand chain(s) default to `X`; override with `--lig-chains` for cases like ZDOCK outputs where ligand is chain `B` (or multiple chains): `python3 scripts/flexref_to_graphs.py --runs ZDOCK_rank1 --lig-chains B`.

## Notes
- HADDOCK prefers single-chain ligands; antibody H/L will drift apart unless you supply lock restraints. The pipeline builds intra-ligand lock restraints (`restraints/intra_lig_lock_<base>.tbl`) when using `run_pipeline.py`.
- `data/` is kept clean; per-run artifacts are archived under `results/run_<base>/`.
- Receptor chain inference depends on matching a `pdbAG/<base>.pdb` file (base = complex name without `refold_` and `_Epi`).
- Graphs store only interface nodes/edges (ligand X ↔ receptor non-X). Default ligand chain for graphing is `X`; set `--lig-chains` if your ligand is not `X`.
- 2026-02-13 (refinement finding): Docking succeeds only when Ab–Ag start with some (even artificial) physical contact; beta-factor hotspot guidance can still place the Ab, but if the initial pose is far with no contact, the docking run fails. Batch of 93 samples took ~5–6 hours (~2–4 minutes per sample).
