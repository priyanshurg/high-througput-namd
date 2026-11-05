# Step 0: PDB Download Scripts

This directory contains automated acquisition scripts for querying the **RCSB PDB Search API v2** and downloading experimental protein-only structures with non-covalently bound ions.

## Scripts Overview

| Script | Ion Target | Non-Polymer Comp ID | Destination Directory |
|:---|:---|:---|:---|
| [`download-nh4.py`](download-nh4.py) | Ammonium ($\text{NH}_4^+$) | `NH4` | `../PDB-Database/nh4_protein_only` |
| [`download-pot.py`](download-pot.py) | Potassium ($\text{K}^+$) | `K` | `../PDB-Database/pot_protein_only` |
| [`download-sod.py`](download-sod.py) | Sodium ($\text{Na}^+$) | `NA` | `../PDB-Database/sod_protein_only` |

## How It Works

Each script performs the following:
1. **API Query Payload**:
   - Matches polymer entity: `entity_poly.rcsb_entity_polymer_type == "Protein"`
   - Matches non-polymer ligand: `rcsb_nonpolymer_instance_annotation.comp_id` to the respective ion ID
   - Filters out covalently bound cases: `rcsb_nonpolymer_instance_annotation.type == "HAS_NO_COVALENT_LINKAGE"`
2. **Pagination**:
   - Queries the API in batches of 100 entries (`paginate: { start, rows }`) until all matching structures are found.
3. **Automated Download**:
   - Uses Biopython (`Bio.PDB.PDBList`) to download all coordinate files in standard PDB format (`.ent` files) into `../PDB-Database/<ion>_protein_only/`.

## Prerequisites

```bash
pip install requests biopython
```

## Usage

Run any script directly from this directory:

```bash
# Download all NH4+ bound protein structures
python download-nh4.py

# Download all K+ bound protein structures
python download-pot.py

# Download all Na+ bound protein structures
python download-sod.py
```
