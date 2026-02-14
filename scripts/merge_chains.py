#!/usr/bin/env python3
"""Merge all chains in a PDB into a single chain and renumber residues.

Usage:
  python merge_chains.py --in data/ligand.pdb --out data/ligand_single.pdb --chain X --start 1

Defaults: chain=X, start=1.

Notes:
- ATOM/HETATM records only.
- Renumbers residues sequentially in the order they appear.
- Preserves insertion codes by discarding them (simple renumber).
"""
import argparse


def merge_chains(inp, outp, chain_id="X", start=1):
    new_lines = []
    resmap = {}  # (orig_chain, resid, icode) -> new_resid
    next_resid = start
    with open(inp) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            orig_chain = line[21]
            resid = line[22:26]
            icode = line[26]
            key = (orig_chain, resid, icode)
            if key not in resmap:
                resmap[key] = next_resid
                next_resid += 1
            new_res = resmap[key]
            # Build new line: set chain, set resid, clear icode
            newline = list(line.rstrip("\n"))
            # chain id
            while len(newline) < 22:
                newline.append(' ')
            newline[21] = chain_id
            # resid columns 23-26 (1-based) -> indices 22:26 in 0-based slicing
            res_str = f"{new_res:4d}"
            for i, ch in enumerate(res_str):
                if 22 + i < len(newline):
                    newline[22 + i] = ch
                else:
                    newline.append(ch)
            # insertion code at col 27 (index 26)
            if len(newline) < 27:
                newline.append(' ')
            else:
                newline[26] = ' '
            new_lines.append("".join(newline) + "\n")
    with open(outp, "w") as out:
        out.writelines(new_lines)
    print(f"Wrote {outp} with chain '{chain_id}' starting at resid {start}. Total residues: {next_resid-start}.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input PDB")
    ap.add_argument("--out", dest="outp", required=True, help="Output PDB")
    ap.add_argument("--chain", dest="chain_id", default="X", help="New chain ID (default X)")
    ap.add_argument("--start", dest="start", type=int, default=1, help="Starting residue number (default 1)")
    args = ap.parse_args()
    merge_chains(args.inp, args.outp, args.chain_id, args.start)


if __name__ == "__main__":
    main()
