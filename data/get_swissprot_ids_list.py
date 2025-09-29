#!/usr/bin/env python3
import sys
import os

# Add src/Composer directory to sys.path so we can import get_id_list
composer_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src', 'Composer')
sys.path.append(composer_path)

import get_id_list

if __name__ == "__main__":
    # UniProtKB Swiss-Prot has ~570k+ reviewed entries (2025), so set a high max_results
    # Adjust limit for API page size (default 500 works well)
    protein_ids = get_id_list.get_swissprot_ids(max_results=100, limit=500)

    print(f"Retrieved {len(protein_ids)} Swiss-Prot accession IDs")

    # Save to file
    output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "swissprot_ids.txt")
    with open(output_file, "w") as f:
        for pid in protein_ids:
            f.write(pid + "\n")

    print(f"Swiss-Prot IDs saved to {output_file}")