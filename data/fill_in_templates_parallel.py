#!/usr/bin/env python3
"""
fill_in_templates_parallel.py

Run preprocessing for ONE shard file (one process handles one txt file).
Designed for SLURM job arrays / distributed runs.

Example:
  python data/fill_in_templates_parallel.py --id_list_path data/train_ids/train_id_7.txt

Outputs (organized per shard):
  data/Sample/train_id_7/EC Numbers/ec_numbers_samples.tsv
  data/Sample/train_id_7/Cofactors/cofactors_samples.tsv
  ...
"""

import sys
import os
import argparse

# -------------------------
# Path setup
# -------------------------
project_root = os.path.abspath(os.getcwd())
src_path = os.path.join(project_root, "src")
sys.path.append(src_path)

data_path = os.path.join(project_root, "data")
maps_path = os.path.join(project_root, "maps")

ec_map_path = os.path.join(maps_path, "ec_mapping.json")

template_path = os.path.join(data_path, "Template_pt")

# Templates
ec_template_a_path = os.path.join(template_path, "EC Numbers", "ec_numbers_templates_a.json")
ec_template_b_path = os.path.join(template_path, "EC Numbers", "ec_numbers_templates_b.json")
ec_template_c_path = os.path.join(template_path, "EC Numbers", "ec_numbers_templates_c.json")
ec_template_d_path = os.path.join(template_path, "EC Numbers", "ec_numbers_templates_d.json")
ec_template_e_path = os.path.join(template_path, "EC Numbers", "ec_numbers_templates_e.json")
ec_template_f_path = os.path.join(template_path, "EC Numbers", "ec_numbers_templates_f.json")
ec_template_g_path = os.path.join(template_path, "EC Numbers", "ec_numbers_templates_g.json")
ec_template_h_path = os.path.join(template_path, "EC Numbers", "ec_numbers_templates_h.json")
ec_template_i_path = os.path.join(template_path, "EC Numbers", "ec_numbers_templates_i.json")
ec_template_j_path = os.path.join(template_path, "EC Numbers", "ec_numbers_templates_j.json")

ca_template_a_path = os.path.join(template_path, "Catalytic Activity", "catalytic_activity_templates_a.json")
ca_template_b_path = os.path.join(template_path, "Catalytic Activity", "catalytic_activity_templates_b.json")

cf_template_a_path = os.path.join(template_path, "Cofactors", "cofactors_templates_a.json")
cf_template_b_path = os.path.join(template_path, "Cofactors", "cofactors_templates_b.json")
cf_template_c_path = os.path.join(template_path, "Cofactors", "cofactors_templates_c.json")
cf_template_d_path = os.path.join(template_path, "Cofactors", "cofactors_templates_d.json")

sl_template_a_path = os.path.join(template_path, "Subcellular Location", "subcellular_location_templates_a.json")
sl_template_b_path = os.path.join(template_path, "Subcellular Location", "subcellular_location_templates_b.json")
sl_template_c_path = os.path.join(template_path, "Subcellular Location", "subcellular_location_templates_c.json")
sl_template_d_path = os.path.join(template_path, "Subcellular Location", "subcellular_location_templates_d.json")

pw_template_a_path = os.path.join(template_path, "Pathway", "pathway_templates_a.json")
pw_template_b_path = os.path.join(template_path, "Pathway", "pathway_templates_b.json")
pw_template_c_path = os.path.join(template_path, "Pathway", "pathway_templates_c.json")
pw_template_d_path = os.path.join(template_path, "Pathway", "pathway_templates_d.json")

go_template_a_path = os.path.join(template_path, "Gene Ontology", "gene_ontology_templates_a.json")
go_template_b_path = os.path.join(template_path, "Gene Ontology", "gene_ontology_templates_b.json")
go_template_c_path = os.path.join(template_path, "Gene Ontology", "gene_ontology_templates_c.json")

rg_template_a_path = os.path.join(template_path, "Regions", "regions_templates_a.json")
rg_template_b_path = os.path.join(template_path, "Regions", "regions_templates_b.json")

st_template_a_path = os.path.join(template_path, "Sites", "sites_templates_a.json")
st_template_b_path = os.path.join(template_path, "Sites", "sites_templates_b.json")
st_template_c_path = os.path.join(template_path, "Sites", "sites_templates_c.json")
st_template_d_path = os.path.join(template_path, "Sites", "sites_templates_d.json")


# -------------------------
# Imports from your project
# -------------------------
from src.Composer.uniprot_utils import (
    extract_uniprot_info,
    extract_uniprot_comments,
    extract_uniprot_features,
    extract_uniprot_cross_references,
)
from src.Composer.template_utils import (
    create_samples_ec,
    create_samples_ca,
    create_samples_cf,
    create_sample_cf_residue_level,
    create_samples_sl,
    create_samples_pw,
    create_samples_go,
    create_samples_rg,
    create_samples_st,
)


# -------------------------
# Output helpers
# -------------------------
def make_output_paths(sample_root: str) -> dict:
    """
    Task-specific TSV paths under sample_root.
    """
    return {
        "ec_sample_path": os.path.join(sample_root, "EC Numbers", "ec_numbers_samples.tsv"),
        "ca_sample_a_path": os.path.join(sample_root, "Catalytic Activity", "catalytic_activity_samples_a.tsv"),
        "ca_sample_b_path": os.path.join(sample_root, "Catalytic Activity", "catalytic_activity_samples_b.tsv"),
        "cf_sample_path": os.path.join(sample_root, "Cofactors", "cofactors_samples.tsv"),
        "cf_sample_c_path": os.path.join(sample_root, "Cofactors", "cofactors_samples_residue_level.tsv"),
        "sl_sample_path": os.path.join(sample_root, "Subcellular Location", "subcellular_location_samples.tsv"),
        "pw_sample_path": os.path.join(sample_root, "Pathway", "pathway_samples.tsv"),
        "go_sample_a_path": os.path.join(sample_root, "Gene Ontology", "gene_ontology_samples_a.tsv"),
        "go_sample_b_path": os.path.join(sample_root, "Gene Ontology", "gene_ontology_samples_b.tsv"),
        "go_sample_c_path": os.path.join(sample_root, "Gene Ontology", "gene_ontology_samples_c.tsv"),
        "rg_sample_path": os.path.join(sample_root, "Regions", "regions_samples.tsv"),
        "st_sample_path": os.path.join(sample_root, "Sites", "sites_samples.tsv"),
    }


def shard_name_from_idfile(id_list_path: str) -> str:
    # train_id_7.txt -> train_id_7
    return os.path.splitext(os.path.basename(id_list_path))[0]


def ensure_parent_dir(path: str) -> None:
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)


# -------------------------
# Main
# -------------------------
def main():
    parser = argparse.ArgumentParser(description="Fill templates for ONE train_id shard file (parallel-safe).")
    parser.add_argument(
        "--id_list_path",
        type=str,
        required=True,
        help="Path to a single txt file containing protein IDs (one per line).",
    )
    parser.add_argument(
        "--out_root",
        type=str,
        default=os.path.join(data_path, "Sample"),
        help="Base output directory. Shard subfolder will be created under this directory.",
    )
    parser.add_argument(
        "--log_every",
        type=int,
        default=50,
        help="Print progress every N proteins (default: 50).",
    )
    args = parser.parse_args()

    if not os.path.exists(args.id_list_path):
        raise FileNotFoundError(f"ID list file not found: {args.id_list_path}")

    shard = shard_name_from_idfile(args.id_list_path)
    shard_out_root = os.path.join(args.out_root, shard)
    os.makedirs(shard_out_root, exist_ok=True)

    out = make_output_paths(shard_out_root)
    for p in out.values():
        ensure_parent_dir(p)

    print(f"=== Shard: {shard} ===")
    print(f"ID file: {args.id_list_path}")
    print(f"Output root: {shard_out_root}")

    processed = 0
    skipped_empty = 0

    with open(args.id_list_path, "r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            protein_id = line.strip()
            if not protein_id:
                skipped_empty += 1
                continue

            # Fetch UniProt info
            meta_info = extract_uniprot_info(protein_id)
            comments_info = extract_uniprot_comments(protein_id)
            cross_references_info = extract_uniprot_cross_references(protein_id)
            features_info = extract_uniprot_features(protein_id)

            # Comments block
            create_samples_ec(
                meta_info, ec_map_path,
                ec_template_a_path, ec_template_b_path, ec_template_c_path, ec_template_d_path,
                ec_template_e_path, ec_template_f_path, ec_template_g_path, ec_template_h_path,
                ec_template_i_path, ec_template_j_path,
                out["ec_sample_path"]
            )
            create_samples_ca(
                comments_info, meta_info,
                ca_template_a_path, ca_template_b_path,
                out["ca_sample_a_path"], out["ca_sample_b_path"]
            )
            create_samples_cf(
                comments_info, meta_info,
                cf_template_a_path, cf_template_b_path,
                out["cf_sample_path"]
            )
            create_sample_cf_residue_level(
                comments_info, features_info, meta_info,
                cf_template_c_path, cf_template_d_path,
                out["cf_sample_c_path"]
            )
            create_samples_sl(
                comments_info, meta_info,
                sl_template_a_path, sl_template_b_path, sl_template_c_path, sl_template_d_path,
                out["sl_sample_path"]
            )
            create_samples_pw(
                comments_info, meta_info,
                pw_template_a_path, pw_template_b_path, pw_template_c_path, pw_template_d_path,
                out["pw_sample_path"]
            )

            # Cross-references block
            create_samples_go(
                cross_references_info, meta_info,
                go_template_a_path, go_template_b_path, go_template_c_path,
                out["go_sample_a_path"], out["go_sample_b_path"], out["go_sample_c_path"]
            )

            # Features block
            create_samples_rg(
                features_info, meta_info,
                rg_template_a_path, rg_template_b_path,
                out["rg_sample_path"]
            )
            create_samples_st(
                features_info, meta_info,
                st_template_a_path, st_template_b_path, st_template_c_path, st_template_d_path,
                out["st_sample_path"]
            )

            processed += 1
            if args.log_every > 0 and processed % args.log_every == 0:
                print(f"[{shard}] processed {processed} proteins")

    print(f"=== Done shard: {shard} ===")
    print(f"Processed: {processed}")
    if skipped_empty:
        print(f"Skipped empty lines: {skipped_empty}")
    print(f"Outputs written under: {shard_out_root}")


if __name__ == "__main__":
    main()
