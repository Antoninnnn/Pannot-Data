import sys
import os
import re
import json
import random
import warnings


import csv

DEFAULT_FIELDS = [
    "task",
    "instructions",
    "sequence",
    "text_label",
    "StartLoc",
    "EndLoc",
    "accession_id",
    "sequence_length",
]

def append_row_to_table(row, table_path, fmt="tsv", fieldnames=DEFAULT_FIELDS):
    os.makedirs(os.path.dirname(table_path), exist_ok=True)
    delimiter = "\t" if fmt.lower() == "tsv" else ","
    write_header = (not os.path.exists(table_path)) or (os.path.getsize(table_path) == 0)

    with open(table_path, "a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter)
        if write_header:
            writer.writeheader()
        writer.writerow(row)

# Add src directory to sys.path
project_root = os.path.dirname(os.path.abspath(os.getcwd()))
src_path = os.path.join(project_root, "src")
sys.path.append(src_path)

from src.Composer.uniprot_utils import natural_join

PATTERN = re.compile(r"\{([^{}]+)\}")

def load_random_entry(json_path):
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, list) or not data:
        raise ValueError("JSON file must be a non-empty list.")
    return random.choice(data)


def append_sample_to_jsonl(sample, jsonl_path):
    os.makedirs(os.path.dirname(jsonl_path), exist_ok=True)
    with open(jsonl_path, "a", encoding="utf-8") as f:
        f.write(json.dumps(sample, ensure_ascii=False) + "\n")


def build_sample(sample_info, meta_info, template, ec_task=False,
                 ec_number=None, ec_first_digit=None, ec_second_digit=None,
                 ec_third_digit=None, ec_map=None):

    def replace_placeholders(text, mapping):
        try:
            return PATTERN.sub(lambda m: mapping[m.group(1)], text)
        except KeyError as e:
            raise ValueError(f"Missing field for placeholder: {e.args[0]}")

    def replace_placeholders_ec(text, ec, first_digit, second_digit, third_digit, map_ec):
        def repl(m):
            placeholder = m.group(1)
            if placeholder == "ECNumbers":
                return ec
            elif placeholder == "1":
                return map_ec[first_digit]["name"]
            elif placeholder == "2":
                return map_ec[first_digit]["subclasses"][second_digit]["name"]
            elif placeholder == "3":
                return map_ec[first_digit]["subclasses"][second_digit]["subsubclasses"][third_digit]
            else:
                raise ValueError(f"Unknown placeholder: {placeholder}")
        return PATTERN.sub(repl, text)

    sample_input = replace_placeholders(template["input"], meta_info)

    if ec_task:
        sample_output = replace_placeholders_ec(
            template["output"], ec_number, ec_first_digit, ec_second_digit, ec_third_digit, ec_map
        )
        start_loc, end_loc = None, None
    else:
        sample_output = replace_placeholders(template["output"], sample_info) if sample_info else template["output"]

        # Pull span if present in sample_info (Regions/Sites/BindingSite/etc.)
        if isinstance(sample_info, dict):
            start_loc = sample_info.get("StartLocation", None)
            end_loc = sample_info.get("EndLocation", None)
        else:
            start_loc, end_loc = None, None

    sample_meta = {
        "accession_id": meta_info.get("Accession"),
        # "entry_name": meta_info.get("EntryName"),
        "sequence_length": meta_info.get("Length"),
        # "molecular_weight": meta_info.get("MolecularWeight")
    }

    sample = {
        "task": template["task"],
        "instructions": template["instruction"],   # if you want field name "instructions"
        "sequence": sample_input,                  # if template["input"] is "{Sequence}"
        "text_label": sample_output,
        "StartLoc": start_loc,
        "EndLoc": end_loc,
        "meta_data": sample_meta
    }

    return sample

def build_sample_csv(sample_info, meta_info, template, ec_task=False,
                 ec_number=None, ec_first_digit=None, ec_second_digit=None,
                 ec_third_digit=None, ec_map=None):
    """
    Returns a flat dict (one row) suitable for CSV/TSV writing.
    Columns:
      task, instructions, sequence, text_label, StartLoc, EndLoc, accession_id, sequence_length
    """

    def replace_placeholders(text, mapping):
        try:
            return PATTERN.sub(lambda m: str(mapping[m.group(1)]), text)
        except KeyError as e:
            raise ValueError(f"Missing field for placeholder: {e.args[0]}")

    def replace_placeholders_ec(text, ec, first_digit, second_digit, third_digit, map_ec):
        def repl(m):
            placeholder = m.group(1)
            if placeholder == "ECNumbers":
                return ec
            elif placeholder == "1":
                return map_ec[first_digit]["name"]
            elif placeholder == "2":
                return map_ec[first_digit]["subclasses"][second_digit]["name"]
            elif placeholder == "3":
                return map_ec[first_digit]["subclasses"][second_digit]["subsubclasses"][third_digit]
            else:
                raise ValueError(f"Unknown placeholder: {placeholder}")
        return PATTERN.sub(repl, text)

    # Template input is typically "{Sequence}" and should be resolvable from meta_info
    sequence_text = replace_placeholders(template["input"], meta_info)

    if ec_task:
        text_label = replace_placeholders_ec(
            template["output"], ec_number, ec_first_digit, ec_second_digit, ec_third_digit, ec_map
        )
        start_loc, end_loc = None, None
    else:
        # Fill normal placeholders (Content/Description/StartLocation/EndLocation/etc.)
        text_label = replace_placeholders(template["output"], sample_info) if sample_info else template["output"]

        if isinstance(sample_info, dict):
            start_loc = sample_info.get("StartLocation")
            end_loc = sample_info.get("EndLocation")
            # Optional: normalize single-site records to have EndLoc = StartLoc
            if start_loc is not None and end_loc is None:
                end_loc = start_loc
        else:
            start_loc, end_loc = None, None

    # Flatten meta into top-level columns for TSV/CSV
    row = {
        "task": template.get("task"),
        "instructions": template.get("instruction"),
        "sequence": sequence_text,
        "text_label": text_label,
        "StartLoc": start_loc,
        "EndLoc": end_loc,
        "accession_id": meta_info.get("Accession"),
        "sequence_length": meta_info.get("Length"),
    }

    return row


def gen_and_append(sample_info, meta_info, template_path, sample_path, ec_task=False, ec_number=None, ec_first_digit=None, ec_second_digit=None, ec_third_digit=None, ec_map=None):
    template = load_random_entry(template_path)
    sample = build_sample_csv(sample_info, meta_info, template, ec_task, ec_number, ec_first_digit, ec_second_digit, ec_third_digit, ec_map)
    append_row_to_table(sample, sample_path)


def create_samples_ec(meta_info, ec_map_path, template_a_path, template_b_path, template_c_path, template_d_path, template_e_path, template_f_path, template_g_path, template_h_path, template_i_path, template_j_path, sample_path):
    ec_numbers = meta_info["ECNumbers"]
    with open(ec_map_path, "r", encoding="utf-8") as f:
        ec_map = json.load(f)

    for ec_number in ec_numbers:
        ec_split = ec_number.split(".")
        first_digit = ec_split[0]
        second_digit = ec_split[1]
        third_digit = ec_split[2]
        first_three_digits = ".".join(ec_split[:3])

        if first_three_digits in ["1.1.1", "1.1.2", "1.1.3", "1.1.4", "1.1.5", "1.1.7",
                                  "1.1.9", "1.1.98", "1.1.99", "1.2.1", "1.2.2", "1.2.3",
                                  "1.2.4", "1.2.5", "1.2.7", "1.2.98", "1.2.99", "1.3.1",
                                  "1.3.2", "1.3.3", "1.3.4", "1.3.5", "1.3.7", "1.3.8",
                                  "1.3.98", "1.3.99", "1.4.1", "1.4.2", "1.4.3", "1.4.4",
                                  "1.4.5", "1.4.7", "1.4.9", "1.4.98", "1.4.99", "1.5.1",
                                  "1.5.3", "1.5.4", "1.5.5", "1.5.7", "1.5.8", "1.5.98",
                                  "1.5.99", "1.6.1", "1.6.2", "1.6.3", "1.6.4", "1.6.5",
                                  "1.6.6", "1.6.7", "1.6.8", "1.6.99", "1.7.1", "1.7.2",
                                  "1.7.3", "1.7.5", "1.7.6", "1.7.7", "1.7.99", "1.8.1",
                                  "1.8.2", "1.8.3", "1.8.4", "1.8.5", "1.8.6", "1.8.7",
                                  "1.8.98", "1.8.99", "1.9.3", "1.9.6", "1.9.98", "1.9.99",
                                  "1.10.1", "1.10.2", "1.10.3", "1.10.5", "1.10.9", "1.10.99",
                                  "1.11.2", "1.12.1", "1.12.2", "1.12.5", "1.12.7", "1.12.98",
                                  "1.12.99", "1.13.11", "1.13.12", "1.14.11", "1.14.12", "1.14.13",
                                  "1.14.14", "1.14.15", "1.14.16", "1.14.17", "1.14.18", "1.14.19",
                                  "1.14.20", "1.14.21", "1.16.1", "1.16.2", "1.16.3", "1.16.5",
                                  "1.16.8", "1.16.9", "1.16.98", "1.16.99", "1.17.1", "1.17.2",
                                  "1.17.3", "1.17.4", "1.17.5", "1.17.7", "1.17.8", "1.17.9",
                                  "1.17.98", "1.17.99", "1.18.1", "1.18.3", "1.18.6", "1.18.96",
                                  "1.18.99", "1.19.1", "1.19.6", "1.20.1", "1.20.2", "1.20.4",
                                  "1.20.9", "1.20.98", "1.20.99", "1.21.1", "1.21.3", "1.21.4",
                                  "1.21.98", "1.21.99", "1.22.1", "1.23.1", "1.23.5"]:
            template_path = template_a_path
        elif first_three_digits in ["1.11.1", "2.1.1", "2.1.2", "2.1.3", "2.1.4", "2.1.5",
                                    "2.2.1", "2.6.1", "2.6.2", "2.6.3", "2.7.1", "2.7.2",
                                    "2.7.3", "2.7.4", "2.7.5", "2.7.6", "2.7.7", "2.7.8",
                                    "2.7.9", "2.7.10", "2.7.11", "2.7.12", "2.7.13", "2.7.14",
                                    "2.7.99", "2.8.1", "2.8.2", "2.8.3", "2.8.5", "2.9.1",
                                    "2.10.1", "3.1.1", "3.1.2", "3.1.3", "3.1.4", "3.1.5",
                                    "3.1.6", "3.1.7", "3.1.8", "3.1.11", "3.1.12", "3.1.13",
                                    "3.1.14", "3.1.15", "3.1.16", "3.1.21", "3.1.22", "3.1.25",
                                    "3.1.26", "3.1.27", "3.1.30", "3.1.31", "3.3.1", "3.3.2",
                                    "3.4.11", "3.4.13", "3.4.14", "3.4.15", "3.4.16", "3.4.17",
                                    "3.4.18", "3.4.19", "3.4.21", "3.4.22", "3.4.23", "3.4.24",
                                    "3.4.25", "3.4.26", "3.4.99", "3.13.2", "6.1.1", "6.1.2",
                                    "6.2.1", "6.3.1", "6.3.2", "6.3.3", "6.3.4", "6.3.5"]:
            template_path = template_b_path
        elif first_three_digits in ["2.3.1", "2.4.99", "4.2.2", "4.2.3", "5.1.1", "5.1.2",
                                    "5.1.3", "5.1.99", "5.3.1", "5.3.2", "5.3.3", "5.3.4",
                                    "5.4.1", "5.4.3", "5.4.4", "5.4.99"]:
            template_path = template_c_path
        elif first_three_digits in ["2.3.2", "2.4.1", "2.4.2", "3.2.1", "3.2.2", "4.1.1",
                                    "4.1.2", "4.1.3", "4.1.99", "4.2.1", "4.2.99", "4.3.1",
                                    "4.3.2", "4.3.3", "4.3.99", "5.3.99", "5.4.2", "5.5.1"]:
            template_path = template_d_path
        elif first_three_digits in ["2.3.3"]:
            template_path = template_e_path
        elif first_three_digits in ["2.6.99", "2.8.4", "3.6.3", "3.6.4", "3.6.5", "3.13.1",
                                    "6.6.1"]:
            template_path = template_f_path
        elif first_three_digits in ["3.5.1", "3.5.2", "3.5.3", "3.5.4", "3.5.5", "3.5.99",
                                    "3.6.1", "3.6.2", "3.7.1", "3.8.1", "3.8.2"]:
            template_path = template_g_path
        elif first_three_digits in ["1.13.99", "1.14.99", "7.1.1", "7.1.2", "7.1.3", "7.2.1",
                                    "7.2.2", "7.2.4", "7.3.2", "7.4.2", "7.6.2"]:
            template_path = template_h_path
        elif first_three_digits in ["1.15.1", "2.5.1", "3.9.1", "3.10.1", "3.11.1", "3.12.1",
                                    "6.4.1", "6.5.1"]:
            template_path = template_i_path
        elif first_three_digits in ["1.97.1", "4.4.1", "4.5.1", "4.6.1", "4.7.1", "4.99.1",
                                    "5.2.1", "5.99.1"]:
            template_path = template_j_path
        else:
            warnings.warn(f"Unknown EC number {ec_number}")
            continue

        gen_and_append(None, meta_info, template_path, sample_path, True, ec_number, first_digit, second_digit, third_digit, ec_map)


def create_samples_ca(comments_info, meta_info, template_a_path, template_b_path, sample_a_path, sample_b_path):
    ca_infos = comments_info["CatalyticActivity"]
    for ca_info in ca_infos:
        if ca_info["Reaction"]:
            gen_and_append(ca_info, meta_info, template_a_path, sample_a_path)
        else:
            raise ValueError(f'Invalid data on catalytic activity: {ca_info}')

        if ca_info["Substrates"] and ca_info["Products"]:
            gen_and_append(ca_info, meta_info, template_b_path, sample_b_path)


def create_samples_cf(comments_info, meta_info, template_a_path, template_b_path, sample_path):
    cf_infos = comments_info["Cofactors"]
    for cf_info in cf_infos:
        if cf_info["Cofactors"] and cf_info["Notes"]:
            cf_info["Notes"] = cf_info["Notes"][0].lower() + cf_info["Notes"][1:]
            gen_and_append(cf_info, meta_info, template_a_path, sample_path)
        elif cf_info["Cofactors"] and not cf_info["Notes"]:
            gen_and_append(cf_info, meta_info, template_b_path, sample_path)
        else:
            raise ValueError(f'Invalid data on cofactors: {cf_info}')


def create_sample_cf_residue_level(comments_info, feature_info, meta_info, template_c_path, template_d_path, sample_path):
    cofactors_at_different_sites = comments_info["Cofactors"]
    binding_sites = feature_info["BindingSite"]
    sample_infos = []
    for cofactors_at_one_site in cofactors_at_different_sites:
        for binding_site in binding_sites:
            if binding_site["Ligand"] in cofactors_at_one_site["ListOfCofactors"]:
                record = {
                    "Cofactors": cofactors_at_one_site["Cofactors"],
                    "StartLocation": binding_site["StartLocation"],
                    "EndLocation": binding_site["EndLocation"]
                }
                if record not in sample_infos:
                    sample_infos.append(record)

    for sample_info in sample_infos:
        if sample_info["StartLocation"] == sample_info["EndLocation"]:
            gen_and_append(sample_info, meta_info, template_c_path, sample_path)
        else:
            gen_and_append(sample_info, meta_info, template_d_path, sample_path)



def create_samples_sl(comments_info, meta_info, template_a_path, template_b_path, template_c_path, template_d_path, sample_path):
    sl_infos = comments_info["SubcellularLocation"]
    for sl_info in sl_infos:
        if sl_info["Name"] == "Canonical" and sl_info["SubcellularLocations"] and sl_info["Notes"]:
            sl_info["Notes"] = sl_info["Notes"][0].lower() + sl_info["Notes"][1:]
            gen_and_append(sl_info, meta_info, template_a_path, sample_path)
        elif sl_info["Name"] == "Canonical" and sl_info["SubcellularLocations"] and not sl_info["Notes"]:
            gen_and_append(sl_info, meta_info, template_b_path, sample_path)
        elif sl_info["Name"] != "Canonical" and sl_info["SubcellularLocations"] and sl_info["Notes"]:
            sl_info["Notes"] = sl_info["Notes"][0].lower() + sl_info["Notes"][1:]
            gen_and_append(sl_info, meta_info, template_c_path, sample_path)
        elif sl_info["Name"] != "Canonical" and sl_info["SubcellularLocations"] and not sl_info["Notes"]:
            gen_and_append(sl_info, meta_info, template_d_path, sample_path)
        else:
            raise ValueError(f'Invalid data on subcellular location: {sl_info}')


def create_samples_pw(comments_info, meta_info, template_a_path, template_b_path, template_c_path, template_d_path, sample_path):
    pw_infos = comments_info["Pathway"]
    for pw_info in pw_infos:
        pw_info["SuperPathway"] = pw_info["SuperPathway"][0].lower() + pw_info["SuperPathway"][1:]

        if pw_info["SuperPathway"] and pw_info["Pathway"] and pw_info["SubPathway"] and pw_info["Step"]:
            gen_and_append(pw_info, meta_info, template_a_path, sample_path)
        elif pw_info["SuperPathway"] and pw_info["Pathway"] and not pw_info["SubPathway"] and not pw_info["Step"]:
            gen_and_append(pw_info, meta_info, template_b_path, sample_path)
        elif pw_info["SuperPathway"] and not pw_info["Pathway"] and not pw_info["SubPathway"] and not pw_info["Step"]:
            gen_and_append(pw_info, meta_info, template_c_path, sample_path)
        elif pw_info["SuperPathway"] and pw_info["Pathway"] and pw_info["SubPathway"] and not pw_info["Step"]:
            gen_and_append(pw_info, meta_info, template_d_path, sample_path)
        else:
            raise ValueError(f'Invalid data on pathway: {pw_info}')


def create_samples_go(cross_references_info, meta_info, template_a_path, template_b_path, template_c_path, sample_a_path, sample_b_path, sample_c_path):
    cross_references_info["GoBP"] = natural_join(cross_references_info["GoBP"], "and")
    cross_references_info["GoMF"] = natural_join(cross_references_info["GoMF"], "and")
    cross_references_info["GoCC"] = natural_join(cross_references_info["GoCC"], "and")
    if cross_references_info["GoBP"]:
        gen_and_append(cross_references_info, meta_info, template_a_path, sample_a_path)
    if cross_references_info["GoMF"]:
        gen_and_append(cross_references_info, meta_info, template_b_path, sample_b_path)
    if cross_references_info["GoCC"]:
        gen_and_append(cross_references_info, meta_info, template_c_path, sample_c_path)


def create_samples_rg(features_info, meta_info, template_a_path, template_b_path, sample_path):
    rg_infos = features_info["Regions"]
    for rg_info in rg_infos:
        if rg_info["Description"]:
            gen_and_append(rg_info, meta_info, template_a_path, sample_path)
        else:
            gen_and_append(rg_info, meta_info, template_b_path, sample_path)


def create_samples_st(features_info, meta_info, template_a_path, template_b_path, template_c_path, template_d_path, sample_path):
    st_infos = features_info["Sites"]
    for st_info in st_infos:
        if st_info["StartLocation"] == st_info["EndLocation"]:
            if st_info["Description"]:
                gen_and_append(st_info, meta_info, template_a_path, sample_path)
            else:
                gen_and_append(st_info, meta_info, template_b_path, sample_path)
        else:
            if st_info["Description"]:
                gen_and_append(st_info, meta_info, template_c_path, sample_path)
            else:
                gen_and_append(st_info, meta_info, template_d_path, sample_path)