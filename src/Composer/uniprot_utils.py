import requests

def get_protein_info(protein_id: str) -> dict:
    """
    Fetch the full UniProtKB record for a protein accession as JSON.

    Parameters
    ----------
    protein_id : str
        UniProtKB accession ID (e.g., "P01308").

    Returns
    -------
    dict
        The UniProtKB entry in JSON format.
    """

    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}"
    response = requests.get(url, params={"format": "json"}, timeout=30)
    response.raise_for_status()
    return response.json()


def get_protein_fields(protein_id: str, fields: str) -> dict:
    """
    Retrieve specific UniProtKB fields for a given protein accession.

    Parameters
    ----------
    protein_id : str
        UniProtKB accession ID (e.g., "P69905").
    fields : str
        Comma-separated list of UniProtKB fields to retrieve.
        See: https://www.uniprot.org/help/return_fields
        Example: "accession,id,gene_names,organism_name,length,sequence,ec,go_f,go_p,go_c,keywords"

    Returns
    -------
    dict
        The JSON response containing the requested fields for the protein.
        If the accession does not exist or has no results, an empty dictionary is returned.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}"
    params = {
        "fields": fields,
        "format": "json"
    }
    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()
    return response.json()


def get_protein_key_value(protein_id: str, key_path: str):
    """
    Retrieve the value corresponding to a nested key path in a UniProtKB JSON entry.
    Stops at lists: will return the entire list if the key resolves to one.

    Parameters
    ----------
    protein_id : str
        UniProt accession ID (e.g., "P01308").
    key_path : str
        Dot-separated path of keys (e.g., "entryAudit" or "organism.scientificName").

    Returns
    -------
    object
        The value corresponding to the given key path.
    """
    data = get_protein_info(protein_id)
    keys = key_path.split(".")
    cur = data
    for k in keys:
        cur = cur[k]
    return cur