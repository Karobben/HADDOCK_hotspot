#!/usr/bin/env python3
"""
Batch pipeline to prepare HADDOCK inputs from complex PDBs.

For each complex PDB in data/:
- Infer receptor chains from ../pdbAG/<base>.pdb where base = complex name stripped of 'refold_' prefix and '_Epi' suffix.
- Split complex into receptor/ligand PDBs (base-specific filenames).
- Merge receptor chains to single chain A; ligand chains to single chain X (base-specific filenames).
- Build intra-ligand lock restraints (H/L interface from split ligand, CA-CA within 8 Å) to keep H/L orientation fixed.
- Build AIR restraints linking receptor hotspots (CA B-factor > -20 on chain A of the complex) to ligand interface residues (CA within 6 Å to receptor in merged files).
- Write a per-complex HADDOCK3 config pointing to base-specific files.
- Log actions to docs/pipeline_logs/<base>.log and append summary to docs/LOG.md.

Usage:
  python scripts/run_pipeline.py --complexes 8IUZ_Epi.pdb another.pdb
"""
import argparse
from pathlib import Path
import datetime
import math
import sys
import subprocess

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
PDBAG = ROOT / "pdbAG"
REST = ROOT / "restraints"
DOCS = ROOT / "docs"
LOG = DOCS / "LOG.md"
PIPELOG_DIR = DOCS / "pipeline_logs"
MERGER = ROOT / "scripts" / "merge_chains.py"


def strip_name(name: str) -> str:
    base = name
    if base.startswith("refold_"):
        base = base[len("refold_"):]
    if base.endswith("_Epi"):
        base = base[:-len("_Epi")]
    return base


def chain_set(pdb: Path):
    chains = set()
    with pdb.open() as f:
        for line in f:
            if line.startswith("ATOM"):
                ch = line[21].strip()
                if ch:
                    chains.add(ch)
    return chains


def write_chains(src: Path, dst: Path, keep: set):
    with src.open() as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith("ATOM") and line[21].strip() in keep:
                fout.write(line)


def run_merger(inp: Path, outp: Path, chain_id: str):
    cmd = ["python", str(MERGER), "--in", str(inp), "--out", str(outp), "--chain", chain_id, "--start", "1"]
    subprocess.run(cmd, check=True)


def load_ca(pdb: Path, chain_filter=None):
    ca = []
    with pdb.open() as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom = line[12:16].strip()
            if atom != "CA":
                continue
            chain = line[21].strip()
            if chain_filter and chain not in chain_filter:
                continue
            resid = int(line[22:26])
            b = float(line[60:66])
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            ca.append((chain, resid, b, (x,y,z)))
    return ca


def closest_pairs(ca1, ca2, cutoff):
    pairs = []
    for _, r1, _, c1 in ca1:
        for _, r2, _, c2 in ca2:
            d = math.sqrt(sum((c1[i]-c2[i])**2 for i in range(3)))
            if d <= cutoff:
                pairs.append((d, r1, r2))
    return pairs


def build_intra_lock(lig_pdb: Path, out_tbl: Path, cutoff=8.0, disable=False):
    if disable:
        with out_tbl.open("w") as f:
            f.write("! Intra-ligand lock disabled for unguided run\n")
        return 0
    caH = load_ca(lig_pdb, chain_filter={"H"})
    caL = load_ca(lig_pdb, chain_filter={"L"})
    # For single-chain ligands (no H or no L), write an empty restraint file so HADDOCK doesn't fail
    if not caH or not caL:
        with out_tbl.open("w") as f:
            f.write("! No H/L pair found; single-chain ligand, no intra-lig lock applied\n")
        return 0
    pairs = closest_pairs(caH, caL, cutoff)
    pairs = sorted(pairs)[:12]
    with out_tbl.open("w") as f:
        f.write("! Intra-ligand restraints to lock H/L interface\n")
        for d, rH, rL in pairs:
            f.write(f"assign (segid H and resid {rH} and name CA)\n"
                    f"       (segid L and resid {rL} and name CA)\n"
                    f"       {d:.2f} 0.5 0.5\n")
    return len(pairs)


def build_air(complex_pdb: Path, rec_single: Path, lig_single: Path, out_tbl: Path, bthr=-20.0, cutoff=6.0, disable=False):
    ca_complex_A = [r for r in load_ca(complex_pdb, chain_filter={"A"}) if r[2] > bthr]
    ca_rec = load_ca(rec_single)
    ca_lig = load_ca(lig_single)

    if disable:
        # Unguided mode: write an empty AIR file so HADDOCK runs without any ambiguous restraints
        with out_tbl.open("w") as f:
            f.write("! AIR disabled for unguided full-surface run\n")
        return 0, set(), set()

    # interface residues by CA distance cutoff
    iface_pairs = closest_pairs(ca_rec, ca_lig, cutoff)
    rec_if = {p[1] for p in iface_pairs}
    lig_if = {p[2] for p in iface_pairs}
    rec_hot = {r[1] for r in ca_complex_A}
    # fallback: if no hotspots survive B-factor filter, use receptor interface residues as hotspots
    if not rec_hot:
        rec_hot = rec_if
    if not rec_hot or not lig_if:
        return 0, rec_if, lig_if
    with out_tbl.open("w") as f:
        f.write("! AIRs between receptor hotspots (A) and ligand interface (X)\n")
        for r in sorted(rec_hot):
            for l in sorted(lig_if):
                f.write(f"assign (segid A and resid {r} and name CA)\n"
                        f"       (segid X and resid {l} and name CA)\n"
                        f"       2.0 2.0 0.0\n")
    return len(rec_hot)*len(lig_if), rec_if, lig_if


def write_config(base: str, air_tbl: Path, lock_tbl: Path, rec_single: Path, lig_single: Path, run_dir: Path, rigid_only: bool = False, sampling: int = 20, full_surface: bool = False):
    cfg = DOCS / f"haddock3_config_{base}.cfg"
    if rigid_only:
        # rigid + minimal flex, no restraints; for full-surface we also omit ambig/unambig
        cfg.write_text(f"""# Auto-generated rigid+flex config for {base}
run_dir = "{run_dir}"
mode = "local"
ncores = 20

molecules = [
    "{rec_single}",
    "{lig_single}"
    ]

[topoaa]

[rigidbody]
sampling = {sampling}
ranair = true
mol_fix_origin_1 = false
mol_fix_origin_2 = false

[flexref]
# keep it simple for speed; no ambig_fname needed here
mol_fix_origin_1 = false
mol_fix_origin_2 = false

[caprieval]
""")
    else:
        if full_surface:
            # full-surface: no restraints referenced, keep ranair true
            cfg.write_text(f"""# Auto-generated full-surface config for {base}
run_dir = "{run_dir}"
mode = "local"
ncores = 20

molecules = [
    "{rec_single}",
    "{lig_single}"
    ]

[topoaa]

[rigidbody]
sampling = {sampling}
ranair = true
mol_fix_origin_1 = false
mol_fix_origin_2 = false

[flexref]
# no ambig/unambig for full-surface
mol_fix_origin_1 = false
mol_fix_origin_2 = false

[emref]
mol_fix_origin_1 = false
mol_fix_origin_2 = false

[caprieval]
""")
        else:
            # masked (guided) run with restraints
            cfg.write_text(f"""# Auto-generated config for {base}
run_dir = "{run_dir}"
mode = "local"
ncores = 20

molecules = [
    "{rec_single}",
    "{lig_single}"
    ]

[topoaa]

[rigidbody]
ambig_fname = "{air_tbl}"
unambig_fname = "{lock_tbl}"
sampling = {sampling}
mol_fix_origin_1 = false
mol_fix_origin_2 = false

[flexref]
ambig_fname = "{air_tbl}"
unambig_fname = "{lock_tbl}"
mol_fix_origin_1 = false
mol_fix_origin_2 = false

[emref]
ambig_fname = "{air_tbl}"
unambig_fname = "{lock_tbl}"
mol_fix_origin_1 = false
mol_fix_origin_2 = false
""")
    return cfg

def run_haddock(cfg: Path, base: str):
    # run haddock3 for this config; capture output to per-base log
    log_file = PIPELOG_DIR / f"haddock_{base}.log"
    haddock_exec = ROOT / ".venv" / "bin" / "haddock3"
    cmd = [str(haddock_exec),  str(cfg)]
    with log_file.open("w") as lf:
        res = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    return res.returncode, log_file


def append_log(base: str, lines: list):
    PIPELOG_DIR.mkdir(parents=True, exist_ok=True)
    plog = PIPELOG_DIR / f"{base}.log"
    ts = datetime.datetime.now().isoformat(timespec='seconds')
    with plog.open("a") as f:
        f.write(f"[{ts}]\n" + "\n".join(lines) + "\n\n")
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.touch()
    with LOG.open("a") as f:
        f.write(f"- {base}: " + "; ".join(lines) + "\n")


def process_one(complex_name: str, full_surface: bool = False, rigid_only: bool = False, sampling: int = 20):
    complex_pdb = DATA / complex_name
    if not complex_pdb.exists():
        raise FileNotFoundError(f"Complex PDB not found: {complex_pdb}")
    base_full = complex_pdb.stem
    base_ref = strip_name(base_full)
    ref_pdb = PDBAG / f"{base_ref}.pdb"
    if not ref_pdb.exists():
        raise FileNotFoundError(f"Reference receptor PDB not found: {ref_pdb}")

    variant_suffix = "_fullsurf" if full_surface else ""
    if rigid_only:
        variant_suffix += "_rigid"
    base_variant = f"{base_full}{variant_suffix}"

    rec_chains = chain_set(ref_pdb)
    all_chains = chain_set(complex_pdb)
    lig_chains = all_chains - rec_chains
    if not lig_chains:
        raise RuntimeError(f"No ligand chains inferred for {complex_name}")

    rec_split = DATA / f"receptor_{base_full}.pdb"
    lig_split = DATA / f"ligand_{base_full}.pdb"
    write_chains(complex_pdb, rec_split, rec_chains)
    write_chains(complex_pdb, lig_split, lig_chains)

    rec_single = DATA / f"receptor_single_{base_full}.pdb"
    lig_single = DATA / f"ligand_single_{base_full}.pdb"
    run_merger(rec_split, rec_single, "A")
    run_merger(lig_split, lig_single, "X")

    merged_complex = DATA / f"complex_merged_{base_full}.pdb"
    try:
        from merge_two_pdbs import merge as merge_two
    except ImportError:
        merge_two = None
    if merge_two:
        try:
            merge_two(rec_single, lig_single, merged_complex)
        except Exception:
            pass

    # restraint files
    lock_tbl = REST / f"intra_lig_lock_{base_full}.tbl"
    lock_count = build_intra_lock(lig_split, lock_tbl, disable=full_surface)

    air_tbl = REST / (f"air_full_{base_full}.tbl" if full_surface else f"air_{base_full}.tbl")
    air_count, rec_if, lig_if = build_air(complex_pdb, rec_single, lig_single, air_tbl, disable=full_surface)

    run_dir = ROOT / (f"results_full/run_{base_full}" if full_surface else f"results/run_{base_full}")
    cfg = write_config(base_variant, air_tbl, lock_tbl, rec_single, lig_single, run_dir, rigid_only=rigid_only, sampling=sampling, full_surface=full_surface)

    # ensure run_dir is clean before running (HADDOCK will create it)
    if run_dir.exists():
        import shutil
        shutil.rmtree(run_dir, ignore_errors=True)

    rc, run_log = run_haddock(cfg, base_variant)

    # archive inputs alongside results after run (do not pre-create before haddock)
    if not run_dir.exists():
        run_dir.mkdir(parents=True, exist_ok=True)
    for src in [rec_split, lig_split, rec_single, lig_single, merged_complex, lock_tbl, air_tbl, cfg]:
        try:
            if src.exists():
                (run_dir / src.name).write_bytes(src.read_bytes())
        except Exception:
            pass

    # clean intermediates from data
    for f in [rec_split, lig_split, rec_single, lig_single, merged_complex]:
        try:
            if f.exists():
                f.unlink()
        except Exception:
            pass

    lines = [
        f"receptor chains {sorted(rec_chains)}; ligand chains {sorted(lig_chains)}",
        f"split -> {rec_split.name}, {lig_split.name}; merged -> {rec_single.name}, {lig_single.name}",
        f"intra lock pairs: {lock_count} -> {lock_tbl.name}",
        f"AIR pairs: {air_count} (rec_if {sorted(rec_if)}, lig_if {sorted(lig_if)}) -> {air_tbl.name}",
        f"config: {cfg.name}",
        f"run haddock3 rc={rc}, log={run_log.name}"
    ]
    append_log(base_variant, lines)
    return base_variant


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complexes", nargs="+", required=True, help="Complex PDB filenames in data/")
    ap.add_argument("--full-surface", action="store_true", help="Disable AIR mask (full surface docking)")
    ap.add_argument("--rigid-only", action="store_true", help="Rigid-body only (skip flexref/emref)")
    ap.add_argument("--sampling", type=int, default=20, help="Sampling for rigidbody")
    args = ap.parse_args()
    for comp in args.complexes:
        try:
            process_one(comp, full_surface=args.full_surface, rigid_only=args.rigid_only, sampling=args.sampling)
        except Exception as e:
            sys.stderr.write(f"[ERROR] {comp}: {e}\n")
            continue


if __name__ == "__main__":
    main()
