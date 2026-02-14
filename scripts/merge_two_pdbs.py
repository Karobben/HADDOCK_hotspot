#!/usr/bin/env python3
"""
Merge two PDBs (receptor + ligand) into a single PDB file.
Keeps records as-is; writes concatenated ATOM/HETATM lines.
Usage:
  python merge_two_pdbs.py --rec data/receptor_single.pdb --lig data/ligand_single.pdb --out data/complex_merged.pdb
"""
import argparse
from pathlib import Path

def merge(rec: Path, lig: Path, out: Path):
    lines = []
    for p in (rec, lig):
        with p.open() as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    lines.append(line)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(lines))
    print(f"Merged {rec.name} + {lig.name} -> {out}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rec", required=True)
    ap.add_argument("--lig", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    merge(Path(args.rec), Path(args.lig), Path(args.out))

if __name__ == "__main__":
    main()
