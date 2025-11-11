#!/usr/bin/env python3
import sys
import os

# Add src/Composer directory to sys.path so we can import get_id_list
composer_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src', 'Composer')
sys.path.append(composer_path)

from src.Composer.get_id_list import get_swissprot_ids

if __name__ == "__main__":
    protein_ids = get_swissprot_ids(limit=10, max_results=1000)
    print("Retrieved", len(protein_ids), "SwissProt IDs")
    print("List of the first 10 ids: ", protein_ids[:10])  # print first 10 IDs
    # Save to a text file
    with open("swissprot_ids.txt", "w") as f:
        for pid in protein_ids:
            f.write(pid + "\n")
    print("SwissProt IDs saved to swissprot_ids.txt")