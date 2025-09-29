#!/usr/bin/env python3
import os
import sys
import json
from typing import List, Dict

# --- Make src importable (so we can do "from Composer import ...") ---
HERE = os.path.dirname(os.path.abspath(__file__))          # /Pannot-Data/data
PROJECT_ROOT = os.path.dirname(HERE)                       # /Pannot-Data
SRC_PATH = os.path.join(PROJECT_ROOT, "src")
COMPOSER_PATH = os.path.join(SRC_PATH, "composer")
if SRC_PATH not in sys.path:
    sys.path.append(SRC_PATH)
if COMPOSER_PATH not in sys.path:
    sys.path.append(COMPOSER_PATH)

# --- Import from Composer package ---
from Composer import template_utils

TEMPLATES_PATH = os.path.join(HERE, "Template", "GO_templates_demo.json")
IDS_PATH = os.path.join(HERE, "swissprot_ids.txt")
OUTPUT_DIR = os.path.join(HERE, "Raw_Sample")
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "GO_samples")  # JSONL file (no extension per spec)

def load_templates(path: str) -> List[Dict[str, str]]:
    """Load the list of GO templates from JSON."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, list):
        raise ValueError("GO_templates_demo.json must contain a JSON list.")
    return data

def load_accession_ids(path: str) -> List[str]:
    """Load accession IDs from a txt file (one per line)."""
    ids = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            acc = line.strip()
            if acc:
                ids.append(acc)
    return ids

def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    templates = load_templates(TEMPLATES_PATH)
    accessions = load_accession_ids(IDS_PATH)

    wrote = 0
    with open(OUTPUT_PATH, "w", encoding="utf-8") as out_f:
        for t_idx, tpl in enumerate(templates, start=1):
            for acc in accessions:
                try:
                    filled = template_utils.fill_template_with_uniprot_fields(tpl, acc)
                    record = {
                        "accession": acc,
                        "template_index": t_idx,
                        "instruction": filled.get("instruction", ""),
                        "input": filled.get("input", ""),
                        "output": filled.get("output", "")
                    }
                    out_f.write(json.dumps(record, ensure_ascii=False) + "\n")
                    wrote += 1
                except Exception as e:
                    err_record = {
                        "accession": acc,
                        "template_index": t_idx,
                        "error": f"{type(e).__name__}: {e}"
                    }
                    out_f.write(json.dumps(err_record, ensure_ascii=False) + "\n")

    print(f"Generated {wrote} samples to {OUTPUT_PATH}")
    print("Format: JSON Lines (one JSON object per line).")

if __name__ == "__main__":
    main()