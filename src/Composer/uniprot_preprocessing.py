from .constants.ec_classification import EC_CLASSIFICATION_MAP
from .constants.feature_description import FEATURE_CONTENT_BY_CATEGORY


def natural_join(lst, conj="and", default=None):
    """Join a list of strings into a natural-language phrase using commas and a conjunction."""
    if not lst:
        return default

    # Normalize & filter
    items = [
        str(x).strip()
        for x in lst
        if x is not None and str(x).strip()
    ]

    if not items:
        return default

    if len(items) == 1:
        return items[0]
    elif len(items) == 2:
        return f" {conj} ".join(items)
    else:
        return ", ".join(items[:-1]) + f", {conj} " + items[-1]


def safe_space_join(lst, default=None):
    """Safely join a list of strings using a single space."""
    if not lst:
        return default

    # Normalize & filter
    items = [
        str(x).strip()
        for x in lst
        if x is not None and str(x).strip()
    ]

    if not items:
        return default

    return " ".join(items)


def lowercase_first_note(notes):
    """Lowercase the first character of the first non-empty note."""
    if not notes:
        return notes

    processed = []
    first_done = False

    for note in notes:
        if note and not first_done:
            note = note.strip()
            if note:
                note = note[0].lower() + note[1:]
                first_done = True
        processed.append(note)

    return processed


def format_cofactor(cofactors):
    """Format cofactor records into template-ready styles."""
    formatted = []

    for record in cofactors:
        cofactors_str = natural_join(record["Cofactors"])
        notes_str = safe_space_join(lowercase_first_note(record["Notes"]))

        formatted.append({
            "Cofactors": cofactors_str,
            "Notes": notes_str,
        })

    return formatted


def format_subcellular_location(subcellular_locations):
    """Format subcellular location records into template-ready styles."""
    formatted = []

    for record in subcellular_locations:
        loc_items = []

        for subloc in record["SubcellularLocations"]:
            location = subloc.get("Location")
            orientation = subloc.get("Orientation")
            topology = subloc.get("Topology")

            if not location:
                continue

            if topology and orientation:
                loc_items.append(f"{location} ({topology}; {orientation})")
            elif orientation:
                loc_items.append(f"{location} ({orientation})")
            elif topology:
                loc_items.append(f"{location} ({topology})")
            else:
                loc_items.append(location)

        subcellular_locations_str = natural_join(loc_items)
        notes_str = safe_space_join(lowercase_first_note(record["Notes"]))

        formatted.append({
            "Name": record["Name"],
            "SubcellularLocations": subcellular_locations_str,
            "Notes": notes_str,
        })

    return formatted


def format_pathway(pathways):
    """Format pathway records into template-ready styles."""
    formatted = []

    for record in pathways:
        for pwy in record:
            parts = pwy.split(": ")
            if len(parts) == 2:
                pathway_str, step = parts
            elif len(parts) == 1:
                pathway_str = parts[0]
                step = None
            else:
                raise ValueError(f"Invalid format of pathway (too many ': '): {pwy!r}")

            levels = pathway_str.split("; ")
            if len(levels) == 3:
                super_pathway, pathway, sub_pathway = levels
            elif len(levels) == 2:
                super_pathway, pathway = levels
                sub_pathway = None
            elif len(levels) == 1:
                super_pathway = levels[0]
                pathway = None
                sub_pathway = None
            else:
                raise ValueError(f"Invalid format of pathway hierarchy: {pwy!r}")

            # Flatten the nested pathway list structure
            formatted.append({
                "SuperPathway": super_pathway,
                "Pathway": pathway,
                "SubPathway": sub_pathway,
                "Step": step,
            })

    return formatted


def split_go(gene_ontologies):
    """Split Gene Ontology (GO) records into BP, MF, and CC categories."""
    bp, mf, cc = [], [], []

    for record in gene_ontologies:
        go_term = record["Term"]
        go_id = record["ID"]

        if len(go_term) < 3:
            raise ValueError(f"Invalid GO term format: {go_term!r}")

        category = go_term[0]
        term = go_term[2:]

        if category == "P":
            bp.append({"ID": go_id, "Term": term})
        elif category == "F":
            mf.append({"ID": go_id, "Term": term})
        elif category == "C":
            cc.append({"ID": go_id, "Term": term})
        else:
            raise ValueError(f"Invalid GO term category: {go_term!r}")

    split_record = {
        "BP": bp,
        "MF": mf,
        "CC": cc,
    }

    return split_record


def format_go(split_record):
    """Format gene ontology (GO) records into template-ready styles based on the split records."""

    def format_item(item):
        return f"{item['Term']} ({item['ID']})"

    bp_items = [format_item(item) for item in split_record["BP"]]
    mf_items = [format_item(item) for item in split_record["MF"]]
    cc_items = [format_item(item) for item in split_record["CC"]]

    formatted_bp = natural_join(bp_items)
    formatted_mf = natural_join(mf_items)
    formatted_cc = natural_join(cc_items)

    formatted = {
        "BP": formatted_bp,
        "MF": formatted_mf,
        "CC": formatted_cc,
    }

    return formatted


def split_and_format_go(gene_ontologies):
    """Split GO annotations into BP/MF/CC and format each category as a template-ready string."""
    split_record = split_go(gene_ontologies)
    formatted_record = format_go(split_record)
    return formatted_record


def format_ec_numbers(ec_numbers):
    """Format EC numbers into template-ready dictionaries."""
    formatted = []

    for ec in ec_numbers:
        parts = str(ec).strip().split(".")
        if len(parts) != 4:
            raise ValueError(f"Invalid EC number: {ec!r}")

        a, b, c, _ = parts  # 4th field exists, but we don't need it for these features
        prefix_abc = ".".join(parts[:3])

        class_info = EC_CLASSIFICATION_MAP.get(a)
        subclass_info = (
            class_info.get("subclasses", {}).get(b)
            if class_info else None
        )
        subsubclass_label = (
            subclass_info.get("subsubclasses", {}).get(c)
            if subclass_info else None
        )

        formatted.append(
            {
                "ECNumber": ec,
                "ClassInfo": class_info["name"] if class_info else None,
                "SubclassInfo": subclass_info["name"] if subclass_info else None,
                "SubSubclassInfo": subsubclass_label,
                "FirstThreeDigits": prefix_abc,
            }
        )

    return formatted


def _format_location(value, modifier):
    """Format a location with optional modifier."""
    if modifier == "UNKNOWN":
        return "unknown"

    if modifier != "EXACT":
        return f"{value} ({modifier.lower()})"

    return str(value)


def _build_basic_feature_record(record, category):
    """Build shared feature fields for a given category (e.g., 'Regions', 'Sites')."""
    return {
        "StartLocation": _format_location(record["Start"], record["StartModifier"]),
        "EndLocation": _format_location(record["End"], record["EndModifier"]),
        "Description": record["Description"],
        "Content": FEATURE_CONTENT_BY_CATEGORY[category][record["Type"]],
    }


def format_regions(regions):
    """Format region annotations into template-ready dictionaries."""
    formatted = []

    for record in regions:
        formatted.append(_build_basic_feature_record(record, "Regions"))

    return formatted


def _format_ligand_record(name, note):
    """Format ligand/ligand_part into a single display string, or None."""
    return (
        f"{name} ({note})" if name and note
        else name if name
        else note if note
        else None
    )


def format_sites(sites):
    """Format site annotations into template-ready dictionaries."""
    formatted = []

    for record in sites:
        basic_record = _build_basic_feature_record(record, "Sites")

        ligand_record = _format_ligand_record(
            record.get("Ligand"),
            record.get("LigandNote"),
        )

        basic_record["Ligand"] = ligand_record

        # Ligand part (reserved)
        # basic_record["LigandPart"] = _format_ligand_record(
        #     record.get("LigandPart"),
        #     record.get("LigandPartNote"),
        # )

        formatted.append(basic_record)

    return formatted