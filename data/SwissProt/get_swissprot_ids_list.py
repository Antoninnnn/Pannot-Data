import sys
import os
# 添加上上级目录的src/Composer到系统路径
composer_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'Composer')
sys.path.append(composer_path)
import get_id_list

if __name__ == "__main__":
    protein_ids = get_id_list.get_swissprot_ids(limit=200, max_results=300)
    print("Retrieved", len(protein_ids), "SwissProt IDs")
    print("List of the first 10 ids: ", protein_ids[:10])  # print first 10 IDs
    # Save to a text file
    with open("swissprot_ids.txt", "w") as f:
        for pid in protein_ids:
            f.write(pid + "\n")
    print("SwissProt IDs saved to swissprot_ids.txt")
