import requests
import os
from Bio.PDB import PDBList

def download_sod_proteins(download_dir="../PDB-Database/sod_protein_only"):
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Combined Query: sod ligand AND Protein only
    query_payload = {
                        "query": {
                            "type": "group",
                            "nodes": [
                            {
                                "type": "group",
                                "nodes": [
                                {
                                    "type": "terminal",
                                    "service": "text",
                                    "parameters": {
                                    "attribute": "entity_poly.rcsb_entity_polymer_type",
                                    "operator": "exact_match",
                                    "value": "Protein"
                                    }
                                },
                                {
                                    "type": "group",
                                    "logical_operator": "and",
                                    "nodes": [
                                    {
                                        "type": "terminal",
                                        "service": "text",
                                        "parameters": {
                                        "attribute": "rcsb_nonpolymer_instance_annotation.comp_id",
                                        "operator": "exact_match",
                                        "value": "NA"
                                        }
                                    },
                                    {
                                        "type": "terminal",
                                        "service": "text",
                                        "parameters": {
                                        "attribute": "rcsb_nonpolymer_instance_annotation.type",
                                        "operator": "exact_match",
                                        "value": "HAS_NO_COVALENT_LINKAGE",
                                        "negation": False
                                        }
                                    }
                                    ],
                                    "label": "nested-attribute"
                                }
                                ],
                                "logical_operator": "and"
                            }
                            ],
                            "logical_operator": "and",
                            "label": "text"
                        },
                        "return_type": "entry",
                        "request_options": {
                            "paginate": {
                            "start": 0,
                            "rows": 100
                            },
                            "results_content_type": [
                            "experimental"
                            ],
                            "sort": [
                            {
                                "sort_by": "score",
                                "direction": "desc"
                            }
                            ],
                            "scoring_strategy": "combined"
                        }
                        }
    all_pdb_ids = []
    start = 0
    rows = 100

    print("Searching for Protein-only structures with sod+...")

    # Pagination loop to fetch all results
    while True:
        query_payload["request_options"]["paginate"]["start"] = start
        response = requests.post(search_url, json=query_payload)
        
        if response.status_code == 204:  # No more results
            break
        
        response.raise_for_status()
        data = response.json()
        
        batch_ids = [hit['identifier'] for hit in data.get('result_set', [])]
        all_pdb_ids.extend(batch_ids)
        
        if len(batch_ids) < rows:
            break
        start += rows

    if not all_pdb_ids:
        print("No matches found.")
        return

    print(f"Total structures found: {len(all_pdb_ids)}")

    # Download Step
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    pdbl = PDBList()
    pdbl.download_pdb_files(all_pdb_ids, pdir=download_dir, file_format="pdb")

if __name__ == "__main__":
    download_sod_proteins()