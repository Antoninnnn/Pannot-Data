# main.py
import time
from get_id_list import get_swissprot_ids
from uniprot_utils import get_protein_info

if __name__ == "__main__":
    # Step 1: get SwissProt IDs
    protein_ids = get_swissprot_ids(limit=100, max_results=5)
    print("Retrieved IDs:", protein_ids)

    # Step 2: loop and fetch info
    for pid in protein_ids:
        info = get_protein_info(pid)
        print(f"\nProtein {pid}:")
        print("  ID:", info.get("primaryAccession"))
        print("  Entry Name:", info.get("uniProtkbId"))
        print("  Protein Name:", info.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"))
        time.sleep(0.2)  # polite delay