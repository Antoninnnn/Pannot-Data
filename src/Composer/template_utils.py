import re
import random
import warnings
import json
import csv

from .uniprot_preprocessing import (
    format_cofactor, format_subcellular_location, format_pathway, format_ec_numbers,
    split_and_format_go, format_regions, format_sites
)
from .constants.ec_template_groups import EC_TEMPLATE_GROUPS


_PLACEHOLDER_PATTERN = re.compile(r"\{([^{}]+)\}")
_META_FIELDS = (
    "Accession", "EntryName", "Organism", "TaxonomyId",
    "Existence", "Sequence", "Length", "MolecularWeight",
    "IsEnzyme",
)

def load_templates(template_path):
    """Load a template pool from a JSON file."""
    with open(template_path, "r", encoding="utf-8") as f:
        templates = json.load(f)

    if not isinstance(templates, list):
        raise ValueError(f"Template file {template_path} must contain a list of templates.")

    return templates


def write_samples(f, samples, fmt):
    """Write samples to jsonl/csv/tsv."""
    if not samples:
        return 0

    fmt = (fmt or "").lower()

    if fmt == "jsonl":
        text = "".join(json.dumps(s, ensure_ascii=False) + "\n" for s in samples)
        f.write(text)
        return len(samples)

    elif fmt in ("csv", "tsv"):
        # Use key order from the first sample as CSV/TSV columns
        fieldnames = list(samples[0].keys())

        try:
            need_header = (f.tell() == 0)
        except Exception:
            need_header = False

        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter=("\t" if fmt == "tsv" else ","),
        )

        if need_header:
            writer.writeheader()

        def _csv_safe(v):
            if v is None:
                return ""
            if isinstance(v, (dict, list)):
                return json.dumps(v, ensure_ascii=False)
            return v

        for s in samples:
            writer.writerow({k: _csv_safe(s.get(k)) for k in fieldnames})

        return len(samples)

    else:
        raise ValueError(f"Unsupported fmt={fmt!r}, only 'csv', 'tsv', 'jsonl' are allowed.")


def default_schema(info, template, sample_input, sample_output):
    """Build a default instruction sample dictionary."""
    sample_meta = {
        "accession": info["Accession"],
        "entry_name": info["EntryName"],
        "organism": info["Organism"],
        "taxonomy_id": info["TaxonomyId"],
        "existence": info["Existence"],
        "length": info["Length"],
        "molecular_weight": info["MolecularWeight"],
        "is_enzyme": info["IsEnzyme"],
    }

    return {
        "index": template["index"],
        "task": template["task"],
        "instruction": template["instruction"],
        "input": sample_input,
        "output": sample_output,
        "meta": sample_meta,
    }


def my_schema(info, template, sample_input, sample_output):
    """
    Build a customized instruction sample dictionary.
    Modify this schema function to match your own dataset format and requirements.
    """
    return {
        "task": template["task"],
        "instructions": template["instruction"],
        "sequence": info["Sequence"],
        "text_label": sample_output,
        "StartLoc": info.get("StartLocation"),
        "EndLoc": info.get("EndLocation"),
        "accession_id": info["Accession"],
        "sequence_length": info["Length"],
    }


def build_sample(info, template, schema=my_schema):
    """Build a single instruction sample by filling placeholders in a template."""

    def replace_placeholder(values, text):
        try:
            return _PLACEHOLDER_PATTERN.sub(lambda m: values[m.group(1)], text)
        except KeyError as e:
            raise ValueError(f"Missing field for placeholder: {e.args[0]}.")

    sample_input = replace_placeholder(info, template["input"])
    sample_output = replace_placeholder(info, template["output"])

    return schema(info, template, sample_input, sample_output)


def build_meta_info(meta):
    """Extract shared UniProt metadata fields for instruction samples."""
    return {key: meta[key] for key in _META_FIELDS}


def build_catalytic_activity_samples(catalytic_reactions, meta, templates):
    """Build one instruction sample per catalytic reaction for a UniProt entry."""
    samples = []

    for reaction in catalytic_reactions:
        info = {
            **build_meta_info(meta),
            "Reaction": reaction,
        }

        template = random.choice(templates)
        samples.append(build_sample(info, template))

    return samples


def build_cofactor_samples(formatted_cofactor_records, meta, templates_a, templates_b):
    """Build one instruction sample per cofactor record for a UniProt entry."""
    samples = []

    for record in formatted_cofactor_records:
        info = {
            **build_meta_info(meta),
            "Cofactors": record["Cofactors"],
            "Notes": record["Notes"],
        }

        pool = templates_a if record["Notes"] else templates_b
        samples.append(build_sample(info, random.choice(pool)))

    return samples


def build_subcellular_location_samples(
    formatted_subcellular_location_records, meta,
    templates_a, templates_b, templates_c, templates_d,
):
    """Build one instruction sample per subcellular-location record for a UniProt entry."""
    samples = []

    for record in formatted_subcellular_location_records:
        name = record["Name"]
        notes = record["Notes"]

        info = {
            **build_meta_info(meta),
            "Name": name,
            "SubcellularLocations": record["SubcellularLocations"],
            "Notes": notes,
        }

        if not name:
            pool = templates_a if notes else templates_b
        else:
            pool = templates_c if notes else templates_d

        samples.append(build_sample(info, random.choice(pool)))

    return samples


def build_pathway_samples(
    formatted_pathway_records, meta,
    templates_a, templates_b, templates_c, templates_d,
):
    """Build one instruction sample per pathway record for a UniProt entry."""
    samples = []

    for record in formatted_pathway_records:
        sp = record["SuperPathway"]
        pw = record["Pathway"]
        sub = record["SubPathway"]
        step = record["Step"]

        info = {
            **build_meta_info(meta),
            "SuperPathway": sp,
            "Pathway": pw,
            "SubPathway": sub,
            "Step": step,
        }

        if sp and pw and sub and step:
            pool = templates_a
        elif sp and pw and not sub and not step:
            pool = templates_b
        elif sp and not pw and not sub and not step:
            pool = templates_c
        elif sp and pw and sub and not step:
            pool = templates_d
        else:
            warnings.warn("Unexpected format of pathway record.")
            continue

        samples.append(build_sample(info, random.choice(pool)))

    return samples


def build_go_samples(formatted_go_records, meta, templates_bp, templates_mf, templates_cc):
    """Build up to three GO instruction samples (BP/MF/CC) for a UniProt entry."""
    bp = formatted_go_records["BP"]
    mf = formatted_go_records["MF"]
    cc = formatted_go_records["CC"]

    info = {**build_meta_info(meta), "BP": bp, "MF": mf, "CC": cc}

    return {
        "BP": build_sample(info, random.choice(templates_bp)) if bp else None,
        "MF": build_sample(info, random.choice(templates_mf)) if mf else None,
        "CC": build_sample(info, random.choice(templates_cc)) if cc else None,
    }


def build_ec_samples(
    formatted_ec_numbers, meta, templates_a, templates_b, templates_c,
    templates_d, templates_e, templates_f, templates_g,
    templates_h, templates_i, templates_j, templates_k,
):
    """Build one instruction sample per EC number record for a UniProt entry."""
    samples = []

    for record in formatted_ec_numbers:
        first_three = record["FirstThreeDigits"]

        info = {
            **build_meta_info(meta),
            "ECNumber": record["ECNumber"],
            "ClassInfo": record["ClassInfo"],
            "SubclassInfo": record["SubclassInfo"],
            "SubSubclassInfo": record["SubSubclassInfo"],
            "FirstThreeDigits": first_three,
        }

        if first_three in EC_TEMPLATE_GROUPS["A"]:
            pool = templates_a
        elif first_three in EC_TEMPLATE_GROUPS["B"]:
            pool = templates_b
        elif first_three in EC_TEMPLATE_GROUPS["C"]:
            pool = templates_c
        elif first_three in EC_TEMPLATE_GROUPS["D"]:
            pool = templates_d
        elif first_three in EC_TEMPLATE_GROUPS["E"]:
            pool = templates_e
        elif first_three in EC_TEMPLATE_GROUPS["F"]:
            pool = templates_f
        elif first_three in EC_TEMPLATE_GROUPS["G"]:
            pool = templates_g
        elif first_three in EC_TEMPLATE_GROUPS["H"]:
            pool = templates_h
        elif first_three in EC_TEMPLATE_GROUPS["I"]:
            pool = templates_i
        elif first_three in EC_TEMPLATE_GROUPS["J"]:
            pool = templates_j
        else:
            pool = templates_k

        samples.append(build_sample(info, random.choice(pool)))

    return samples


def build_regions_samples(formatted_regions_records, meta, templates_a, templates_b):
    """Build one instruction sample per regions record for a UniProt entry."""
    samples = []

    for record in formatted_regions_records:
        info = {
            **build_meta_info(meta),
            "StartLocation": record["StartLocation"],
            "EndLocation": record["EndLocation"],
            "Content": record["Content"],
            "Description": record["Description"],
        }

        pool = templates_a if record["Description"] else templates_b

        samples.append(build_sample(info, random.choice(pool)))

    return samples


def build_sites_samples(
        formatted_sites_record, meta, templates_a,
        templates_b, templates_c, templates_d
):
    """Build one instruction sample per sites record for a UniProt entry."""
    samples = []

    for record in formatted_sites_record:
        info = {
            **build_meta_info(meta),
            "StartLocation": record["StartLocation"],
            "EndLocation": record["EndLocation"],
            "Content": record["Content"],
            "Description": record["Description"],
            # "Ligand": record["Ligand"],
            # "LigandPart": record["LigandPart"],
        }

        same = record["StartLocation"] == record["EndLocation"]
        pool = (templates_a if record["Description"] else templates_b) if same else \
            (templates_c if record["Description"] else templates_d)

        samples.append(build_sample(info, random.choice(pool)))

    return samples


def gen_catalytic_activity_samples(comments, meta, templates):
    """Generate catalytic activity instruction samples."""
    catalytic_reactions = comments["CatalyticActivity"]
    return build_catalytic_activity_samples(catalytic_reactions, meta, templates)


def gen_cofactor_samples(comments, meta, templates_a, templates_b):
    """Generate cofactor instruction samples."""
    formatted_cofactors = format_cofactor(comments["Cofactor"])
    return build_cofactor_samples(formatted_cofactors, meta, templates_a, templates_b)


def gen_subcellular_location_samples(comments, meta, templates_a, templates_b, templates_c, templates_d):
    """Generate subcellular location instruction samples."""
    formatted_subcellular_locations = format_subcellular_location(comments["SubcellularLocation"])
    return build_subcellular_location_samples(formatted_subcellular_locations, meta, templates_a, templates_b, templates_c, templates_d)


def gen_pathway_samples(comments, meta, templates_a, templates_b, templates_c, templates_d):
    """Generate pathway instruction samples."""
    formatted_pathways = format_pathway(comments["Pathway"])
    return build_pathway_samples(formatted_pathways, meta, templates_a, templates_b, templates_c, templates_d)


def gen_go_samples(cross_references, meta, templates_bp, templates_mf, templates_cc):
    """Generate gene ontology (GO) instruction samples."""
    formatted_gene_ontologies = split_and_format_go(cross_references["GeneOntology"])
    return build_go_samples(formatted_gene_ontologies, meta, templates_bp, templates_mf, templates_cc)


def gen_ec_samples(
        meta, templates_a, templates_b, templates_c,
        templates_d, templates_e, templates_f, templates_g,
        templates_h, templates_i, templates_j, templates_k,
):
    """Generate EC number instruction samples."""
    formatted_ec_numbers = format_ec_numbers(meta["ECNumber"])
    ec_samples = build_ec_samples(
        formatted_ec_numbers, meta, templates_a, templates_b, templates_c,
        templates_d, templates_e, templates_f, templates_g,
        templates_h, templates_i, templates_j, templates_k,
    )
    return ec_samples


def gen_regions_samples(features, meta, templates_a, templates_b):
    """Generate regions instruction samples."""
    formatted_regions = format_regions(features["Regions"])
    return build_regions_samples(formatted_regions, meta, templates_a, templates_b)


def gen_sites_samples(
        features, meta, templates_a,
        templates_b, templates_c, templates_d
):
    """Generate sites instruction samples."""
    formatted_sites = format_sites(features["Sites"])
    return build_sites_samples(
        formatted_sites, meta, templates_a,
        templates_b, templates_c, templates_d
    )