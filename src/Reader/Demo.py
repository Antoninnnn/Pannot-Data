import gzip
import io
import json
import time
from pathlib import Path

import pandas as pd
import requests

BASE = "https://rest.uniprot.org"

# 1) TSV → DataFrame (cursor pagination)
def uniprotkb_tsv_to_df(query, fields, size=500, reviewed=True, taxon=None, max_pages=None):
    """
    Return a pandas DataFrame for a UniProtKB query using TSV format + cursor pagination.
    """
    params = {
        "query": query,
        "fields": ",".join(fields),
        "format": "tsv",
        "size": size,
        "compressed": "true",   # server will gzip the response
    }
    if reviewed is not None:
        params["query"] = f"({params['query']}) AND (reviewed:{str(reviewed).lower()})"
    if taxon:
        params["query"] = f"({params['query']}) AND (organism_id:{taxon})"

    url = f"{BASE}/uniprotkb/search"
    all_chunks = []
    next_cursor = None
    page = 0

    while True:
        if next_cursor:
            params["cursor"] = next_cursor
        r = requests.get(url, params=params, timeout=60)
        r.raise_for_status()

        # decompress gzip -> TSV
        buf = io.BytesIO(r.content)
        with gzip.GzipFile(fileobj=buf) as gz:
            tsv = gz.read().decode("utf-8", errors="replace")

        # first page includes header, later pages too; let pandas handle it
        chunk = pd.read_csv(io.StringIO(tsv), sep="\t")
        # When no results remain, UniProt still returns headers only
        if chunk.shape[0] == 0:
            break
        all_chunks.append(chunk)

        # parse next cursor from Link header if present
        # Example header: Link: <...&cursor=xxxx&size=500>; rel="next"
        link = r.headers.get("Link", "")
        next_cursor = None
        for part in link.split(","):
            if 'rel="next"' in part:
                # cursor is in the URL; extract after 'cursor=' and before '&' or '>'
                s = part.find("cursor=")
                if s != -1:
                    s += len("cursor=")
                    e = part.find("&", s)
                    if e == -1:
                        e = part.find(">", s)
                    next_cursor = part[s:e]
        page += 1
        if max_pages and page >= max_pages:
            break
        # be a good citizen
        time.sleep(0.2)

    return pd.concat(all_chunks, ignore_index=True) if all_chunks else pd.DataFrame()

# Example: build your “instruction-tuning dataframe”
fields = [
    "accession", "id", "organism_name", "gene_primary",
    "protein_name", "ec", "go_f", "go_p", "go_c", "keywordid",
    "xref_pfam", "xref_pdb", "length", "sequence"
]
df = uniprotkb_tsv_to_df(
    query="reviewed:true AND existence:1",  # Swiss-Prot, protein level evidence
    fields=fields,
    size=500,
    taxon=9606,   # Human as example
    max_pages=2   # demo cap; remove to get all
)

# print(df.head())
# 2) FASTA for a query (stream to file)
def download_fasta(query, out_path, reviewed=True, taxon=None):
    params = {"query": query, "format": "fasta", "compressed": "true", "size": 500}
    if reviewed is not None:
        params["query"] = f"({params['query']}) AND (reviewed:{str(reviewed).lower()})"
    if taxon:
        params["query"] = f"({params['query']}) AND (organism_id:{taxon})"

    url = f"{BASE}/uniprotkb/stream"
    r = requests.get(url, params=params, timeout=120)
    r.raise_for_status()
    Path(out_path).write_bytes(r.content)  # this is gzipped FASTA
    return out_path

# download_fasta("enzyme", "uniprot_enzyme.fa.gz", reviewed=True)

# 3) Full JSON → JSONL archive (rich records for later mining)
def jsonl_dump_by_accessions(accessions, out_jsonl):
    url_tpl = f"{BASE}/uniprotkb/{{acc}}.json"
    with open(out_jsonl, "w", encoding="utf-8") as w:
        for acc in accessions:
            r = requests.get(url_tpl.format(acc=acc), timeout=60)
            if r.status_code == 200:
                w.write(json.dumps(r.json(), ensure_ascii=False) + "\n")
            time.sleep(0.1)
    return out_jsonl

# # Example: dump JSON for first 50 from the dataframe above
# if not df.empty:
#     jsonl_dump_by_accessions(df["Entry"].head(50).tolist(), "uniprot_records.jsonl")