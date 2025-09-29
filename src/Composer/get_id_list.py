import requests

def get_swissprot_ids(max_results: int = 5000, limit: int = 500) -> list[str]:
    """
    Retrieve UniProtKB Swiss-Prot accession IDs from the UniProtKB API.

    Parameters
    ----------
    max_results : int
        Maximum number of accession IDs to fetch in total.
    limit : int
        Maximum number of results to fetch per request (page size).

    Returns
    -------
    list[str]
        List of retrieved accession IDs.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": "reviewed:true",   # Swiss-Prot entries
        "fields": "accession",      # Only accession IDs
        "format": "json",
        "size": limit
    }

    accessions = []
    cursor = None

    while len(accessions) < max_results:
        if cursor:
            params["cursor"] = cursor

        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()

        # Extract accessions
        for entry in data.get("results", []):
            if "primaryAccession" in entry:
                accessions.append(entry["primaryAccession"])
                if len(accessions) >= max_results:
                    break

        # Pagination: stop if no nextCursor
        cursor = data.get("nextCursor")
        if not cursor:
            break

    return accessions


if __name__ == "__main__":
    ids = get_swissprot_ids(max_results=20, limit=10)
    print(f"Retrieved {len(ids)} Swiss-Prot accession IDs:")
    for i, acc in enumerate(ids, 1):
        print(f"{i}: {acc}")