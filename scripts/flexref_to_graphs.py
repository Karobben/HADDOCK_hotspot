#!/usr/bin/env python3
"""
Convert HADDOCK flexref PDBs (results/<ID>/2_flexref/*.pdb) into interface-only graph objects.
Outputs a single .pth per run: Graphs/flexref/<ID>.pth containing a list of PyG Data objects
(one per PDB), keeping only interface nodes/edges (ligand X ↔ receptor non-X within cutoff).
"""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import numpy as np
import torch
from scipy.spatial import cKDTree
from torch_geometric.data import Data

AA3_TO_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
    'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}
AA20 = np.array(list("ARNDCQEGHILKMFPSTWYV"), dtype="U1")
AA_TO_IDX = {aa: i for i, aa in enumerate(AA20.tolist())}
EYE20 = np.eye(20, dtype=np.float32)

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
OUTDIR = ROOT / "Graphs" / "flexref"


def get_ca_arrays(pdb_file: Path):
    res3_list = []
    chain_list = []
    resid_list = []
    coords_list = []
    with open(pdb_file, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            res3 = line[17:20].strip()
            chain = line[21].strip()
            try:
                resnum = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            res3_list.append(res3)
            chain_list.append(chain)
            resid_list.append(resnum)
            coords_list.append((x, y, z))
    if not coords_list:
        return (np.empty((0,), dtype="U3"), np.empty((0,), dtype="U1"), np.empty((0,), dtype=np.int32), np.empty((0, 3), dtype=np.float32))
    return (np.array(res3_list, dtype="U3"), np.array(chain_list, dtype="U1"), np.asarray(resid_list, dtype=np.int32), np.asarray(coords_list, dtype=np.float32))


def pdb_to_interface_graph(pdb_file: Path, chain_lig=("X",), cutoff: float = 10.0):
    res3, chains, resnums, coords = get_ca_arrays(pdb_file)
    chain_lig = np.array(chain_lig, dtype="U1")
    chain_rec = np.array([c for c in np.unique(chains) if c not in set(chain_lig.tolist())], dtype="U1")
    if chain_lig.size == 0:
        # if user passes empty ligand set, treat everything as receptor -> no interface
        return None
    if coords.shape[0] == 0:
        return None

    # cross-chain (lig<->rec) pairs only
    tree = cKDTree(coords)
    pairs = tree.query_pairs(cutoff)
    pairs = np.array(list(pairs), dtype=np.int64) if pairs else np.empty((0,2),dtype=np.int64)
    if pairs.size == 0:
        return None
    # filter to ligand-receptor pairs
    lig_mask = np.isin(chains, chain_lig)
    rec_mask = np.isin(chains, chain_rec)
    keep_lr = (lig_mask[pairs[:,0]] & rec_mask[pairs[:,1]]) | (rec_mask[pairs[:,0]] & lig_mask[pairs[:,1]])
    pairs = pairs[keep_lr]
    if pairs.size == 0:
        return None

    nodes_inter = np.unique(pairs.reshape(-1))
    if nodes_inter.size == 0:
        return None

    # build edges among interface nodes within cutoff (re-query on subset)
    sub_idx = nodes_inter
    sub_coords = coords[sub_idx]
    sub_tree = cKDTree(sub_coords)
    sub_pairs = sub_tree.query_pairs(cutoff)
    sub_pairs = np.array(list(sub_pairs), dtype=np.int64) if sub_pairs else np.empty((0,2),dtype=np.int64)
    # directed edges + self-loops
    if sub_pairs.size > 0:
        e0, e1 = sub_pairs[:,0], sub_pairs[:,1]
        edge_index = np.concatenate([np.stack([e0,e1],axis=0), np.stack([e1,e0],axis=0)], axis=1)
    else:
        edge_index = np.empty((2,0), dtype=np.int64)
    n_sub = sub_coords.shape[0]
    self_edges = np.arange(n_sub, dtype=np.int64)
    edge_index = np.concatenate([edge_index, np.stack([self_edges,self_edges], axis=0)], axis=1)

    # map residue/chain/aa to subset
    aa1 = np.array([AA3_TO_1.get(r, "X") for r in res3], dtype="U1")
    idx = np.fromiter((AA_TO_IDX.get(a, -1) for a in aa1), dtype=np.int64, count=len(aa1))
    x_all = np.zeros((len(aa1), 20), dtype=np.float32)
    ok = idx >= 0
    x_all[ok] = EYE20[idx[ok]]

    resi_all = np.array([f"{a}{c}{r}" for a,c,r in zip(aa1.tolist(), chains.tolist(), resnums.tolist())], dtype=object)

    x = x_all[sub_idx]
    pos = coords[sub_idx]
    resi = resi_all[sub_idx]

    # edge attributes
    src, dst = edge_index
    edge_dist = np.linalg.norm(pos[src] - pos[dst], axis=1).astype(np.float32)
    sub_chains = chains[sub_idx]
    same_chain = (sub_chains[src] == sub_chains[dst])
    edge_attr = np.zeros((edge_index.shape[1], 2), dtype=np.float32)
    edge_attr[same_chain, 0] = 1.0
    edge_attr[~same_chain, 1] = 1.0
    edge_attr *= edge_dist[:, None]

    return {
        "x": x,
        "edge_index": edge_index,
        "edge_attr": edge_attr,
        "resi": resi,
        "pos": pos,
        "key": os.path.basename(pdb_file)
    }


def dic_to_pyg(d):
    return Data(
        x=torch.from_numpy(d["x"]).float(),
        edge_index=torch.from_numpy(d["edge_index"]).long(),
        edge_attr=torch.from_numpy(d["edge_attr"]).float(),
        resi=d["resi"],
        pos=torch.from_numpy(d["pos"]).float(),
        key=d["key"],
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="*", help="Run dirs containing 2_flexref (default: all results/*/2_flexref)")
    ap.add_argument("--lig-chains", nargs="*", default=["X"], help="Ligand chain IDs (default: X)")
    args = ap.parse_args()

    if args.runs:
        run_dirs = [Path(r) for r in args.runs]
    else:
        run_dirs = [p for p in RESULTS.glob("run_*/2_flexref") if p.is_dir()]

    OUTDIR.mkdir(parents=True, exist_ok=True)

    for rdir in run_dirs:
        # allow passing either the 2_flexref dir or a run dir containing it
        if rdir.name == "2_flexref" and rdir.is_dir():
            flex_dir = rdir
            base_dir = rdir.parent
        else:
            flex_dir = rdir / "2_flexref"
            base_dir = rdir
        if not flex_dir.is_dir():
            # fallback: if rdir itself contains pdbs, treat it as the flex dir
            pdbs_direct = list(rdir.glob("*.pdb"))
            if pdbs_direct:
                flex_dir = rdir
                base_dir = rdir
            else:
                continue
        run_id = base_dir.name or flex_dir.name or flex_dir.stem
        out_path = OUTDIR / f"{run_id}.pth"

        pdbs = sorted(p for p in flex_dir.glob("*.pdb") if p.is_file())
        if not pdbs:
            continue
        graphs = []
        for pdb in pdbs:
            g = pdb_to_interface_graph(pdb, chain_lig=args.lig_chains)
            if g is not None:
                graphs.append(dic_to_pyg(g))
        if not graphs:
            continue
        torch.save(graphs, out_path)
        print(f"saved {len(graphs)} interface graphs -> {out_path}")


if __name__ == "__main__":
    main()
