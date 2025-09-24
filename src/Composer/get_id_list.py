import requests

def get_swissprot_ids(limit=500, max_results=5000):
    """
    Fetch SwissProt protein accession IDs from UniProtKB.
    
    Parameters:
        limit (int): number of results per request (max 500).
        max_results (int): maximum number of results to fetch in total.
    
    Returns:
        ids (list[str]): list of UniProt accession IDs.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    query = "reviewed:true"  # SwissProt = reviewed entries
    fields = "accession"
    
    ids = []
    cursor = None
    
    while len(ids) < max_results:
        params = {
            "query": query,
            "fields": fields,
            "format": "json",
            "size": limit,
        }
        if cursor:
            params["cursor"] = cursor
        
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        
        for result in data.get("results", []):
            ids.append(result["primaryAccession"])
        
        # Stop if no more results
        if "nextCursor" not in data:
            break
        
        cursor = data["nextCursor"]
    
    return ids


# Example usage:
if __name__ == "__main__":
    protein_ids = get_swissprot_ids(limit=200, max_results=300)
    print("Retrieved", len(protein_ids), "SwissProt IDs")
    print("List of the first 10 ids: ", protein_ids[:10])  # print first 10 IDs