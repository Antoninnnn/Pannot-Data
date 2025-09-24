# uniprot_utils.py
import requests

def get_protein_info(protein_id):
    """
    Fetch detailed protein information for a given SwissProt ID.

    Parameters:
        protein_id (str): UniProt accession ID (e.g., "P69905")

    Returns:
        dict: JSON response with protein information
    """
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
    response = requests.get(url)
    response.raise_for_status()
    return response.json()


def get_protein_fields(protein_id, fields):
    """

    Some useful fields you can request:

    accession

    id (entry name)

    gene_names

    organism_name

    length

    sequence

    protein_name

    ec

    go_f (GO molecular function)

    go_p (GO biological process)

    go_c (GO cellular component)

    keywords

    xref_pdb (cross references to PDB)

    xref_ensembl (cross references to Ensembl)
    """


    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}"
    params = {
        "fields": fields,
        "format": "json"
    }
    response = requests.get(url, params=params)
    response.raise_for_status()
    return response.json()

def list_uniprot_attributes(protein_id):
    # 获取蛋白质信息
    data = get_protein_info(protein_id)
    # 打印蛋白质ID的顶级键信息
    print(f"Top-level keys for {protein_id}:")
    # 遍历并打印所有键
    for key in data.keys():
        print(" -", key)

    # 返回所有键的列表
    return list(data.keys())


def extract_uniprot_fields(protein_id):
    data = get_protein_info(protein_id)
    
    info = {
        "Accession": data.get("primaryAccession"),
        "EntryName": data.get("uniProtkbId"),
        "Organism": data.get("organism", {}).get("scientificName"),
        "GeneSymbols": [g["geneName"]["value"] for g in data.get("genes", []) if "geneName" in g],
        "ProteinName": data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
        "Sequence": data.get("sequence", {}).get("value"),
        "Length": data.get("sequence", {}).get("length"),
        "MolecularWeight": data.get("sequence", {}).get("molWeight"),
        "ECNumbers": [ec["value"] for ec in data.get("proteinDescription", {}).get("recommendedName", {}).get("ecNumbers", [])],
        "Keywords": [kw["name"] for kw in data.get("keywords", [])],
        "CrossReferences": [ref["id"] for ref in data.get("dbReferences", [])],
    }
    return info

if __name__ == "__main__":
    protein_id = "P69905"
    fields = extract_uniprot_fields(protein_id)
    for k, v in fields.items():
        print(f"{k}: {v}")

    protein_id = "P69905"  # Hemoglobin alpha
    attributes = list_uniprot_attributes(protein_id)
    print("\nAttributes:", attributes)
