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
    """
    Safely extract common UniProt fields for a given protein.
    Returns a dict with None or [] if fields are missing.
    """
    data = get_protein_info(protein_id)

    # --- helpers ---
    def safe_get(d, keys, default=None):
        """Safely walk nested dict with a list of keys."""
        for k in keys:
            if isinstance(d, dict) and k in d:
                d = d[k]
            else:
                return default
        return d

    def safe_list(lst, key):
        """Extract a list of values for given key from list of dicts."""
        return [item.get(key) for item in lst if key in item]

    # --- extraction ---
    info = {
        "Accession": data.get("primaryAccession"),
        "EntryName": data.get("uniProtkbId"),
        "Organism": safe_get(data, ["organism", "scientificName"]),
        "GeneSymbols": [g["geneName"]["value"] for g in data.get("genes", []) if "geneName" in g],
        "ProteinName": safe_get(data, ["proteinDescription", "recommendedName", "fullName", "value"]),
        "Sequence": safe_get(data, ["sequence", "value"]),
        "Length": safe_get(data, ["sequence", "length"]),
        "MolecularWeight": safe_get(data, ["sequence", "molWeight"]),
        "ECNumbers": [ec["value"] for ec in safe_get(data, ["proteinDescription", "recommendedName", "ecNumbers"], [])],
        "Keywords": safe_list(data.get("keywords", []), "name"),
        "CrossReferences": safe_list(data.get("dbReferences", []), "id"),
        "Publications": [
            safe_get(ref, ["citation", "title"])
            for ref in data.get("references", [])
            if safe_get(ref, ["citation", "title"])
        ],
        "Functions": [
            safe_get(c, ["texts", 0, "value"])
            for c in data.get("comments", [])
            if c.get("commentType") == "FUNCTION" and safe_get(c, ["texts", 0, "value"])
        ],
    }

    return info

# Example usage
if __name__ == "__main__":
    protein_id = "P69905"  # Hemoglobin alpha
    fields = extract_uniprot_fields(protein_id)
    for k, v in fields.items():
        print(f"{k}: {v}")

    protein_id = "P69905"  # Hemoglobin alpha
    attributes = list_uniprot_attributes(protein_id)
    print("\nAttributes:", attributes)
