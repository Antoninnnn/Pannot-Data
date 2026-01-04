import sys
import os

# Add src directory to sys.path
project_root = os.path.abspath(os.getcwd())
src_path = os.path.join(project_root, "src")
sys.path.append(src_path)

# r/w paths
data_path = os.path.join(project_root, "data")
maps_path = os.path.join(project_root, "maps")
id_list_path = os.path.join(data_path, "enzyme_id_list_test.txt")

ec_map_path = os.path.join(maps_path, "ec_mapping.json")

template_path = os.path.join(data_path, "Template_pt")
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

sample_path = os.path.join(data_path, "Sample")
ec_sample_path = os.path.join(sample_path, "EC Numbers", "ec_numbers_samples.jsonl")
ca_sample_a_path = os.path.join(sample_path, "Catalytic Activity", "catalytic_activity_samples_a.jsonl")
ca_sample_b_path = os.path.join(sample_path, "Catalytic Activity", "catalytic_activity_samples_b.jsonl")
cf_sample_path = os.path.join(sample_path, "Cofactors", "cofactors_samples.jsonl")
cf_sample_c_path = os.path.join(sample_path, "Cofactors", "cofactors_samples_residue_level.jsonl")
sl_sample_path = os.path.join(sample_path, "Subcellular Location", "subcellular_location_samples.jsonl")
pw_sample_path = os.path.join(sample_path, "Pathway", "pathway_samples.jsonl")
go_sample_a_path = os.path.join(sample_path, "Gene Ontology", "gene_ontology_samples_a.jsonl")
go_sample_b_path = os.path.join(sample_path, "Gene Ontology", "gene_ontology_samples_b.jsonl")
go_sample_c_path = os.path.join(sample_path, "Gene Ontology", "gene_ontology_samples_c.jsonl")
rg_sample_path = os.path.join(sample_path, "Regions", "regions_samples.jsonl")
st_sample_path = os.path.join(sample_path, "Sites", "sites_samples.jsonl")

from src.Composer.uniprot_utils import (extract_uniprot_info, extract_uniprot_comments, extract_uniprot_features,
                                        extract_uniprot_cross_references)
from src.Composer.template_utils import (create_samples_ec, create_samples_ca, create_samples_cf,
                                         create_sample_cf_residue_level, create_samples_sl, create_samples_pw,
                                         create_samples_go, create_samples_rg, create_samples_st)

def main():
    with open(id_list_path, "r", encoding="utf-8") as f:
        for i, line in enumerate(f):
            protein_id = line.strip()
            meta_info = extract_uniprot_info(protein_id)
            comments_info = extract_uniprot_comments(protein_id)
            cross_references_info = extract_uniprot_cross_references(protein_id)
            features_info = extract_uniprot_features(protein_id)

            # Sample creation with information in the "comments" block
            create_samples_ec(meta_info, ec_map_path, ec_template_a_path, ec_template_b_path, ec_template_c_path, ec_template_d_path, ec_template_e_path, ec_template_f_path, ec_template_g_path, ec_template_h_path, ec_template_i_path, ec_template_j_path, ec_sample_path)
            create_samples_ca(comments_info, meta_info, ca_template_a_path, ca_template_b_path, ca_sample_a_path, ca_sample_b_path)
            create_samples_cf(comments_info, meta_info, cf_template_a_path, cf_template_b_path, cf_sample_path)
            create_sample_cf_residue_level(comments_info, features_info, meta_info, cf_template_c_path, cf_template_d_path, cf_sample_c_path)
            create_samples_sl(comments_info, meta_info, sl_template_a_path, sl_template_b_path, sl_template_c_path, sl_template_d_path, sl_sample_path)
            create_samples_pw(comments_info, meta_info, pw_template_a_path, pw_template_b_path, pw_template_c_path, pw_template_d_path, pw_sample_path)

            # Sample creation with information in the "cross-references" block
            create_samples_go(cross_references_info, meta_info, go_template_a_path, go_template_b_path, go_template_c_path, go_sample_a_path, go_sample_b_path, go_sample_c_path)

            # Sample creation with information in the "features" block
            create_samples_rg(features_info, meta_info, rg_template_a_path, rg_template_b_path, rg_sample_path)
            create_samples_st(features_info, meta_info, st_template_a_path, st_template_b_path, st_template_c_path, st_template_d_path, st_sample_path)

            if not (i % 50) and i != 0:
                print("Filled with 50 proteins")


if __name__ == "__main__":
    main()