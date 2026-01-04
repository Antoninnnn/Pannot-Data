import sys
import os
import json

# Add src directory to sys.path
project_root = os.path.abspath(os.getcwd())
src_path = os.path.join(project_root, "src")
sys.path.append(src_path)

# r/w paths
data_path = os.path.join(project_root, "data")
id_list_path = os.path.join(data_path, "enzyme_id_list.txt")

from src.Composer.uniprot_utils import get_protein_info, extract_uniprot_info, extract_uniprot_comments, extract_uniprot_cross_references, extract_uniprot_features


def test1():
    with open(id_list_path, "r", encoding="utf-8") as f:

        for line in f:
            protein_id = line.strip()
            print(json.dumps(extract_uniprot_features(protein_id), indent = 4))

    print("done")

def extract_uniprot_features_test(protein_id):
    feature_classification = {
        "MoleculeProcessing": ["Initiator methionine", "Signal", "Transit peptide",
                               "Propeptide", "Chain", "Peptide"],
        "Regions": ["Topological domain", "Transmembrane", "Intramembrane",
                    "Domain", "Repeat", "Calcium binding",
                    "Zinc finger", "DNA binding", "Nucleotide binding",
                    "Region", "Coiled coil", "Motif",
                    "Compositional bias"],
        "Sites": ["Active site", "Metal binding", "Binding site", "Site"],
        "AminoAcidModifications": ["Non-standard residue", "Modified residue", "Lipidation",
                                   "Glycosylation", "Disulfide bond", "Cross-link"],
        "NaturalVariations": ["Alternative sequence", "Natural variant"],
        "ExperimentalInfo": ["Mutagenesis", "Sequence uncertainty", "Sequence conflict",
                             "Non-adjacent residues", "Non-terminal residue"],
        "SecondaryStructure": ["Helix", "Turn", "Beta strand"]
    }
    data = get_protein_info(protein_id)

    sample_count = 0
    location_count = 0
    alternative_sequence_count = 0
    ligand_count = 0
    ligand_part_count = 0
    description_count = 0
    for ftr in data.get("features", []):
        if ftr["type"] == "Metal binding":
            type = ftr["type"]
            sample_count += 1
            if ftr.get("location"):
                location_count += 1
            if ftr.get("alternativeSequence"):
                print(f"yes {type}")
                alternative_sequence_count += 1
            if not ftr.get("alternativeSequence"):
                print(f"no {type}")
            if ftr.get("ligand"):
                ligand_count += 1
            if ftr.get("ligandPart"):
                ligand_part_count += 1
            if ftr.get("description"):
                description_count += 1

    statistics_regions = {
        "sample_count": sample_count,
        "location_count": location_count,
        "alternative_sequence_count": alternative_sequence_count,
        "ligand_count": ligand_count,
        "ligand_part_count": ligand_part_count,
        "description_count": description_count,
    }

    return statistics_regions



def test2():
    with open(id_list_path, "r", encoding="utf-8") as f:
        for line in f:
            protein_id = line.strip()
            sample = extract_uniprot_cross_references(protein_id)
            print(json.dumps(sample, indent = 4))


if __name__ == "__main__":
    test2()
