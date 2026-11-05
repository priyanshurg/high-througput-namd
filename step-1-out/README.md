# Step 1: Simulation Outputs

This directory is the designated output location for the high-throughput molecular dynamics pipeline executed by [`prep_run.py`](../step-1-run-highthroughput-namd/prep_run.py).

> **Note**: Simulation outputs (trajectories, PSF/PDB files, binary coordinates, restart files) are excluded from git version control due to their large size.

## Expected Directory Structure

When running `prep_run.py`, an output directory (such as `runs-nh4/` or `runs-pot/`) is created with subdirectories for each processed PDB entry:

```text
step-1-out/
├── README.md
└── runs-nh4/
    ├── summary.tsv                      # Overall status manifest for all runs
    ├── pdb1cpd/
    │   ├── P001.pdb, P002.pdb ...       # Individual protein segment coordinate files
    │   ├── I001.pdb ...                 # Bound ion coordinate files
    │   ├── segments.txt                 # Manifest of segments loaded into psfgen
    │   ├── salt.txt                     # Salt cation specification
    │   ├── disu.txt                     # Detected disulfide bridge connections (if any)
    │   ├── resmap.tsv                   # Mapping: segid & new_resid <-> original chain & resid
    │   ├── prot_ions.psf / .pdb         # Initial complex with guessed hydrogens
    │   ├── solvated.psf / .pdb          # System after TIP3P water box addition
    │   ├── autoion.psf / .pdb           # System after autoionization
    │   ├── ionized.psf / .pdb           # Final neutralized, solvated simulation system
    │   ├── cell.tcl                     # Periodic boundary cell dimensions & origin
    │   ├── charge.txt                   # Net charge verification
    │   ├── build.log                    # PSF building and solvation log
    │   ├── min.namd / min.log           # 5,000-step unconstrained minimization setup & log
    │   ├── min.coor / .vel / .xsc       # Minimized coordinates, velocities, and cell vectors
    │   ├── select_fixed.tcl / select.log# Pocket identification script and log
    │   ├── fixed.pdb                    # Restraint reference: beta=0 (free pocket/water/ion), beta=1 (frozen)
    │   ├── pocket_residues.txt          # Residues within the pocket cutoff (e.g. 7 Å of ion)
    │   ├── eq.namd / eq.log             # 500 ps NVT pocket-restrained equilibration setup & log
    │   ├── eq.dcd                       # Trajectory coordinates (1 ps / frame)
    │   └── eq.restart.*                 # Trajectory checkpoint restart files
    └── pdb1d7i/
        └── ...
```

## Status Summary File (`summary.tsv`)

A run summary is generated at the root of each run directory:

```tsv
entry       status    note
pdb1cpd     ok        1 segs, 289 res, 1 ions
pdb1d7h     skipped   non-standard residue MSE A105 in ATOM records
pdb1d7i     ok        1 segs, 312 res, 2 ions
```

### Possible Status Values:
- **`ok`**: Structure fully built, minimized, and equilibrated with `eq.dcd` generated.
- **`prepped`**: Structure prepared through solvation (`--prep-only` mode).
- **`already_done`**: Skipped because an existing `eq.dcd` was found (use `--force` to re-run).
- **`skipped`**: Filtered out during initial cleaning (e.g., non-canonical backbone residues, missing ions, or no protein records).
- **`failed`**: An error occurred during VMD/psfgen or NAMD execution (inspect the corresponding `.log` file in the folder).
