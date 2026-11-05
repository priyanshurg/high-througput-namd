#!/usr/bin/env python3
"""
High-throughput ion-binding-site pipeline (CHARMM36m + VMD/psfgen + NAMD).

Per PDB/ENT file:
  1. Clean: first model only, altloc A (else first altloc), drop waters / ligands /
     hydrogens, KEEP ions (NA, K, NH4). Skip the entry if any non-canonical residue.
  2. Build PSF with psfgen, solvate (padding 20 A)                  -> build.tcl
  3. Minimise 5000 steps (everything free)                           -> min.namd
  4. Freeze protein residues with NO atom within 7 A of an ion       -> select_fixed.tcl
  5. 500 ps NVT equilibration, dcdfreq 5 ps (100 frames)             -> eq.namd

Usage example:
  python prep_run.py --input-dir ../PDB-Database/nh4_protein_only --out-dir runs_nh4 \
      --ff-dir /path/to/toppar \
      --jobs 4 --cores 4 --prep-only

Drop the --prep-only flag to run the full minimisation + equilibration. The script will
skip any entries that have already been processed (unless --force is specified).
Requires: Python 3.8+, Biopython, VMD, NAMD, CHARMM36m force field files (top_all36_prot.rtf, par_all36m_pro
"""
import argparse, math, os, shlex, shutil, subprocess, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent

CANON = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
         "MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
# HIS -> HSD (delta-protonated) by default; change here if you want HSE
HIS_MAP = {"HIS":"HSD","HID":"HSD","HSD":"HSD","HIE":"HSE","HSE":"HSE","HIP":"HSP","HSP":"HSP"}
# PDB ion resname -> (CHARMM resname, CHARMM atom name)   <-- VERIFY the NH4 entry against your topology
ION_MAP = {"NA":("SOD","SOD"), "SOD":("SOD","SOD"),
           "K":("POT","POT"),  "POT":("POT","POT"),
           "NH4":("NH4","NZ")}
WATER = {"HOH","WAT","H2O","DOD","TIP","TIP3","SOL"}
BACKBONE = {"N","CA","C"}


# ----------------------------------------------------------------- PDB cleaning
def parse(path):
    atoms = []
    with open(path) as fh:
        for ln in fh:
            rec = ln[:6].strip()
            if rec == "ENDMDL":
                break                                   # first model only
            if rec not in ("ATOM", "HETATM"):
                continue
            ln = ln.rstrip("\n").ljust(80)
            try:
                atoms.append(dict(
                    rec=rec, name=ln[12:16].strip(), alt=ln[16].strip(), res=ln[17:20].strip(),
                    chain=ln[21], resid=int(ln[22:26]), icode=ln[26].strip(),
                    xyz=(float(ln[30:38]), float(ln[38:46]), float(ln[46:54])),
                    elem=ln[76:78].strip().upper()))
            except ValueError:
                continue
    return atoms

def is_h(a):
    if a["elem"]:
        return a["elem"] in ("H", "D")
    return a["name"].lstrip("0123456789")[:1] in ("H", "D")

def pick_altlocs(atoms):
    alts = {}
    for a in atoms:
        if a["alt"]:
            alts.setdefault((a["chain"], a["resid"], a["icode"]), set()).add(a["alt"])
    choice = {k: ("A" if "A" in v else sorted(v)[0]) for k, v in alts.items()}
    return [a for a in atoms
            if not a["alt"] or a["alt"] == choice[(a["chain"], a["resid"], a["icode"])]]

def clean(path):
    """-> (protein_residues, ions, skip_reason)"""
    atoms = pick_altlocs(parse(path))
    order, groups = [], {}
    for a in atoms:
        k = (a["chain"], a["resid"], a["icode"], a["res"])
        if k not in groups:
            groups[k] = []; order.append(k)
        groups[k].append(a)
    prot, ions = [], []
    for k in order:
        chain, resid, icode, rn = k
        at = groups[k]
        if rn in WATER:
            continue
        if rn in CANON or rn in HIS_MAP:
            prot.append((chain, resid, icode, HIS_MAP.get(rn, rn), [a for a in at if not is_h(a)]))
        elif at[0]["rec"] == "HETATM" and rn in ION_MAP:
            ions.append((rn, [a for a in at if not is_h(a)]))
        elif at[0]["rec"] == "ATOM":
            return None, None, f"non-standard residue {rn} {chain}{resid} in ATOM records"
        elif BACKBONE <= {a["name"] for a in at}:
            return None, None, f"modified amino acid {rn} {chain}{resid}"
        # any other HETATM (ligand, glycerol, SO4, other metals ...) is simply dropped
    if not prot:
        return None, None, "no protein residues"
    if not ions:
        return None, None, "no NA/K/NH4 ion found"
    return prot, ions, None


# ------------------------------------------------------- segments / psfgen inputs
def _xyz(at, name):
    for a in at:
        if a["name"] == name:
            return a["xyz"]
    return None

def build_segments(prot):
    segs, cur, prev = [], [], None
    for r in prot:
        chain, resid, icode, rn, at = r
        if prev is not None:
            if chain != prev[0]:
                brk = True
            else:
                c, n = _xyz(prev[4], "C"), _xyz(at, "N")
                brk = math.dist(c, n) > 2.0 if (c and n) else abs(resid - prev[1]) > 1
            if brk:
                segs.append(cur); cur = []
        cur.append(r); prev = r
    if cur:
        segs.append(cur)
    return segs

def pdbline(i, name, res, resid, xyz, seg):
    nm = (" " + name.ljust(3)) if len(name) < 4 else name
    return (f"ATOM  {i % 100000:5d} {nm} {res:>3s}  {resid:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00      {seg:<4s}")

def write_inputs(prot, ions, wd):
    segs = build_segments(prot)
    manifest, rmap, sg = [], ["segid\tnew_resid\tchain\torig_resid\ticode\tresname"], {}
    for si, seg in enumerate(segs, 1):
        segid = f"P{si:03d}"
        lines, n = [], 0
        for ri, (chain, resid, icode, rn, at) in enumerate(seg, 1):
            last = ri == len(seg)
            rmap.append(f"{segid}\t{ri}\t{chain.strip() or '-'}\t{resid}\t{icode or '-'}\t{rn}")
            for a in at:
                nm = a["name"]
                if last and nm == "O":     nm = "OT1"
                elif last and nm == "OXT": nm = "OT2"
                elif nm == "OXT":          continue
                n += 1
                lines.append(pdbline(n, nm, rn, ri, a["xyz"], segid))
                if rn == "CYS" and nm == "SG":
                    sg[(segid, ri)] = a["xyz"]
        (wd / f"{segid}.pdb").write_text("\n".join(lines) + "\nEND\n")
        first = {"GLY": "GLYP", "PRO": "PROP"}.get(seg[0][3], "NTER")
        manifest.append(f"protein {segid} {segid}.pdb {first}")
    for k, (rn, at) in enumerate(ions, 1):
        segid = f"I{k:03d}"
        crn, can = ION_MAP[rn]
        a = next((x for x in at if x["name"] == can or x["name"] == "N"), at[0])
        (wd / f"{segid}.pdb").write_text(pdbline(1, can, crn, 1, a["xyz"], segid) + "\nEND\n")
        manifest.append(f"ion {segid} {segid}.pdb {crn}")
    (wd / "segments.txt").write_text("\n".join(manifest) + "\n")
    crn0, can0 = ION_MAP[ions[0][0]]            # salt cation = same type as the first site ion
    (wd / "salt.txt").write_text(f"{crn0} {can0}\n")
    (wd / "resmap.tsv").write_text("\n".join(rmap) + "\n")
    keys = list(sg)
    ss = [f"{a[0]} {a[1]} {b[0]} {b[1]}" for i, a in enumerate(keys) for b in keys[i+1:]
          if math.dist(sg[a], sg[b]) < 2.5]
    if ss:
        (wd / "disu.txt").write_text("\n".join(ss) + "\n")
    return len(segs), len(prot), len(ions)


# ------------------------------------------------------------------ run one entry
def sh(cmd, log, cwd, env):
    with open(log, "w") as f:
        subprocess.run(cmd, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT,
                       stdin=subprocess.DEVNULL, check=True)

def process(pdb, a):
    eid = Path(pdb).stem
    wd = Path(a.out_dir).resolve() / eid
    if (wd / "eq.dcd").exists() and not a.force:
        return eid, "already_done", ""
    prot, ions, why = clean(pdb)
    if why:
        return eid, "skipped", why
    wd.mkdir(parents=True, exist_ok=True)
    nseg, nres, nion = write_inputs(prot, ions, wd)
    for f in ("build.tcl", "select_fixed.tcl", "common.namd", "min.namd", "eq.namd"):
        shutil.copy(HERE / f, wd / f)
    env = os.environ.copy()
    env.update(FF_DIR=str(Path(a.ff_dir).resolve()), EXTRA_TOP=a.extra_top, EXTRA_PAR=a.extra_par,
               PADDING=str(a.padding), POCKET_CUTOFF=str(a.pocket_cutoff), SALT_CONC=str(a.salt_conc))
    namd = shlex.split(a.namd) + [f"+p{a.cores}", "+devices", str(a.gpuid)]
    vmd = shlex.split(a.vmd) + ["-dispdev", "text", "-e"]
    try:
        sh(vmd + ["build.tcl"], wd / "build.log", wd, env)
        if not (wd / "solvated.psf").exists():
            return eid, "failed", "psfgen/solvate did not produce solvated.psf (see build.log)"
        if a.prep_only:
            return eid, "prepped", f"{nseg} segs, {nres} res, {nion} ions"
        sh(namd + ["min.namd"], wd / "min.log", wd, env)
        sh(vmd + ["select_fixed.tcl"], wd / "select.log", wd, env)
        sh(namd + ["eq.namd"], wd / "eq.log", wd, env)
    except subprocess.CalledProcessError as e:
        return eid, "failed", f"{e.cmd[-1]} exited {e.returncode} (check the .log in {wd})"
    return eid, "ok", f"{nseg} segs, {nres} res, {nion} ions"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--ff-dir", required=True, help="dir with top_all36_prot.rtf, par_all36m_prot.prm, toppar_water_ions_namd.str")
    p.add_argument("--extra-top", default="", help="colon-separated extra RTF/STR files (e.g. NH4 topology)")
    p.add_argument("--extra-par", default="", help="colon-separated extra PRM/STR files (e.g. NH4 parameters)")
    p.add_argument("--padding", type=float, default=20.0)
    p.add_argument("--salt-conc", type=float, default=0.15, help="M; 0 = neutralise only")
    p.add_argument("--pocket-cutoff", type=float, default=7.0)
    p.add_argument("--jobs", type=int, default=1, help="structures processed in parallel")
    p.add_argument("--cores", type=int, default=4, help="NAMD cores per structure")
    p.add_argument("--vmd", default="vmd")
    p.add_argument("--namd", default="namd3", help='e.g. "namd2" or "namd3" or "charmrun ++local namd2"')
    p.add_argument("--prep-only", action="store_true", help="stop after build/solvate (inspect first!)")
    p.add_argument("--force", action="store_true")
    p.add_argument("--gpu-id", default="0", help='specify GPU ID to use for NAMD (e.g. "0" or "0,1" or "0,1,2,3"). Can only specify one GPU ID if using gpu-resident mode.')
    a = p.parse_args()

    files = sorted(set(Path(a.input_dir).glob("*.ent")) | set(Path(a.input_dir).glob("*.pdb")))
    print(f"{len(files)} structures found")
    Path(a.out_dir).mkdir(parents=True, exist_ok=True)
    rows = []
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(process, str(f), a) for f in files]
        for fu in as_completed(futs):
            r = fu.result(); rows.append(r)
            print("\t".join(r), flush=True)
    with open(Path(a.out_dir) / "summary.tsv", "w") as f:
        f.write("entry\tstatus\tnote\n")
        for r in sorted(rows):
            f.write("\t".join(r) + "\n")

if __name__ == "__main__":
    main()