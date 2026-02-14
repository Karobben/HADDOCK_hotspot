#!/usr/bin/env python3
"""
Prepare receptor/ligand splits and merged single-chain PDBs from a complex PDB.
Receptor chains are inferred from a reference PDB in ../pdbAG matching the complex name
(after stripping prefix 'refold_' and suffix '_Epi'). Remaining chains are ligand.

Outputs (in data/):
  receptor.pdb / ligand.pdb           (split by chain sets)
  receptor_single.pdb / ligand_single.pdb (merged, renumbered to chain A/X)

Usage:
  python prepare_from_complex.py --complex 8IUZ_Epi.pdb

Notes:
- If reference PDB is missing, aborts with a clear message.
- Chain letters in the reference PDB define receptor chains.
- Merging uses scripts/merge_chains.py logic.
"""
import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
PDBAG = ROOT.parent / "pdbAG"
MERGER = ROOT / "scripts" / "merge_chains.py"


def strip_name(name: str) -> str:
    base = name
    if base.startswith("refold_"):
        base = base[len("refold_"):]
    if base.endswith("_Epi"):
        base = base[:-len("_Epi")]
    return base


def get_chain_set(pdb_path: Path):
    chains = set()
    with pdb_path.open() as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chains.add(line[21].strip())
    chains.discard("")
    return chains


def write_chains(in_pdb: Path, out_pdb: Path, keep: set):
    with in_pdb.open() as fin, out_pdb.open("w") as fout:
        for line in fin:
            if not line.startswith("ATOM"):
                continue
            ch = line[21].strip()
            if ch in keep:
                fout.write(line)


def run_merger(inp: Path, outp: Path, chain_id: str):
    os.system(f"{MERGER} --in {inp} --out {outp} --chain {chain_id} --start 1")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", required=True, help="Complex PDB filename in data/")
    args = ap.parse_args()

    complex_pdb = DATA / args.complex
    if not complex_pdb.exists():
        sys.exit(f"Complex PDB not found: {complex_pdb}")

    # find receptor reference
    stripped = strip_name(complex_pdb.stem)
    ref_pdb = PDBAG / f"{stripped}.pdb"
    if not ref_pdb.exists():
        sys.exit(f"Reference receptor PDB not found: {ref_pdb}")

    receptor_chains = get_chain_set(ref_pdb)
    if not receptor_chains:
        sys.exit("No chains found in receptor reference")

    # all chains in complex
    all_chains = get_chain_set(complex_pdb)
    ligand_chains = all_chains - receptor_chains
    if not ligand_chains:
        sys.exit("No ligand chains inferred; check reference PDB")

    # split
    rec_out = DATA / "receptor.pdb"
    lig_out = DATA / "ligand.pdb"
    write_chains(complex_pdb, rec_out, receptor_chains)
    write_chains(complex_pdb, lig_out, ligand_chains)

    # merge single-chain
    rec_single = DATA / "receptor_single.pdb"
    lig_single = DATA / "ligand_single.pdb"
    run_merger(rec_out, rec_single, "A")
    run_merger(lig_out, lig_single, "X")

    print(f"Receptor chains: {sorted(receptor_chains)} -> {rec_out}, {rec_single}")
    print(f"Ligand chains: {sorted(ligand_chains)} -> {lig_out}, {lig_single}")


if __name__ == "__main__":
    main()
