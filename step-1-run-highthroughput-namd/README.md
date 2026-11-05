# Step 1: High-Throughput NAMD Pipeline

This directory contains the core pipeline scripts, VMD/psfgen templates, NAMD configurations, and CHARMM36m force field files for preparing and simulating protein ion-binding pockets in high throughput.

## Files Overview

| File | Type | Description |
|:---|:---|:---|
| [`prep_run.py`](prep_run.py) | Python Script | Parallel batch driver that orchestrates PDB cleaning, segmentation, VMD/psfgen, and NAMD execution. |
| [`build.tcl`](build.tcl) | VMD / TCL | Headless VMD script using `psfgen` and `solvate` to build the PSF, add TIP3P water padding, autoionize, and write cell dimensions. |
| [`select_fixed.tcl`](select_fixed.tcl) | VMD / TCL | Analyzes minimized coordinates, selects residues within the pocket cutoff (default 7.0 Å), freezes distant protein (`beta=1`), and outputs `fixed.pdb` and `pocket_residues.txt`. |
| [`common.namd`](common.namd) | NAMD Config | Common simulation parameters: CHARMM36m parameter loading, PME electrostatics, 10–12 Å cutoff, periodic boundary conditions. |
| [`min.namd`](min.namd) | NAMD Config | 5,000 steps unconstrained conjugate-gradient energy minimization. |
| [`eq.namd`](eq.namd) | NAMD Config | 500 ps NVT equilibration at 310 K with Langevin dynamics, fixed bulk atoms, and trajectory output (`eq.dcd`). |
| [`toppar/`](toppar/) | Directory | Core CHARMM36m topology RTF and parameter PRM files. |
| [`inputs-toppar/`](inputs-toppar/) | Directory | Extended parameter streams for various lipids, ions, cofactors, and ligands. |

## Quick Start

```bash
# Dry-run: clean, build PSF, solvate and ionize only (inspect before running MD)
python prep_run.py \
    --input-dir ../PDB-Database/nh4_protein_only \
    --out-dir ../step-1-out/runs-nh4 \
    --ff-dir ./toppar \
    --jobs 4 \
    --cores 4 \
    --prep-only

# Full run: minimization, pocket selection, and NVT equilibration
python prep_run.py \
    --input-dir ../PDB-Database/nh4_protein_only \
    --out-dir ../step-1-out/runs-nh4 \
    --ff-dir ./toppar \
    --jobs 4 \
    --cores 4
```

See the [Root README](../README.md) for full documentation, CLI options, and pipeline details.
