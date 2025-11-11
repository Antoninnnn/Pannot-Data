import requests
import warnings

def get_protein_info(protein_id):
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
    response = requests.get(url)
    response.raise_for_status()
    return response.json()


def safe_get(d, keys, default=None):
    for k in keys:
        if isinstance(d, dict) and k in d:
            d = d[k]
        else:
            return default
    return d


def safe_list(lst, key):
    return [item.get(key) for item in lst if key in item]


def natural_join(lst, conj="and"):
    if not lst:
        return ""
    else:
        if len(lst) == 1:
            return lst[0]
        elif len(lst) == 2:
            return f" {conj} ".join(lst)
        elif len(lst) >= 3:
            return ", ".join(lst[:-1]) + f", {conj} " + lst[-1]


def warn_general(protein_id, field, field_idx, field_type, message):
    warnings.warn(f"[General Warning] Protein: {protein_id} | Field: {field}[{field_idx}] | Field Type: {field_type} | {message}")


def warn_unexpected_key_missing(protein_id, field, field_idx, field_type, path):
    warnings.warn(f"[Unexpected Key Missing] Protein: {protein_id} | Field: {field}[{field_idx}] | Field Type: {field_type} | Missing Path: {path}")


def extract_uniprot_info(protein_id):
    data = get_protein_info(protein_id)

    info = {
        "Accession": data.get("primaryAccession"),
        "EntryName": data.get("uniProtkbId"),
        "Organism": safe_get(data, ["organism", "scientificName"]),
        "GeneSymbols": [g["geneName"]["value"] for g in data.get("genes", []) if "geneName" in g],
        "ProteinName": safe_get(data, ["proteinDescription", "recommendedName", "fullName", "value"]),
        "Sequence": safe_get(data, ["sequence", "value"]),
        "Length": safe_get(data, ["sequence", "length"]),
        "MolecularWeight": safe_get(data, ["sequence", "molWeight"]),
        "ECNumbers": [ec["value"] for ec in safe_get(data, ["proteinDescription", "recommendedName", "ecNumbers"], [])],
        "Keywords": safe_list(data.get("keywords", []), "name"),
        "CrossReferences": safe_list(data.get("dbReferences", []), "id"),
        "Publications": [safe_get(ref, ["citation", "title"]) for ref in data.get("references", []) if safe_get(ref, ["citation", "title"])],
        "Functions": [safe_get(c, ["texts", 0, "value"]) for c in data.get("comments", []) if c.get("commentType") == "FUNCTION" and safe_get(c, ["texts", 0, "value"])]
    }

    return info


def extract_uniprot_comments(protein_id):
    data = get_protein_info(protein_id)

    field = "comments"

    catalytic_activity_collection = []
    cofactors_collection = []
    subcellular_location_collection = []
    pathway_collection = []
    for i, cmt in enumerate(data.get("comments", [])):
        # Catalytic Activity
        # Ref: UniProt Help / Catalytic activity: https://www.uniprot.org/help/catalytic_activity
        # Examples: "O36015", "P17050", and "A0A0S3QTD0"
        if cmt["commentType"] == "CATALYTIC ACTIVITY":
            # Get the catalytic reaction of the protein
            rxn = cmt.get("reaction", {})
            if rxn:
                reaction = rxn["name"]
                left_and_right = reaction.split(" = ")
                if len(left_and_right) == 2:
                    substrates = natural_join(left_and_right[0].split(" + "), "and")
                    products = natural_join(left_and_right[1].split(" + "), "and")
                elif len(left_and_right) == 1:
                    substrates = None
                    products = None
                else:
                    warn_general(protein_id, field, i, cmt["commentType"], f'["reaction"] has an unexpected format: {reaction}')
                    continue
            else:
                warn_unexpected_key_missing(protein_id, field, i, cmt["commentType"], '["reaction"]')
                continue

            record = {
                "Reaction": reaction,
                "Substrates": substrates,
                "Products": products,
            }

            catalytic_activity_collection.append(record)


        # Cofactors
        # Ref: UniProt Help / Cofactor: https://www.uniprot.org/help/cofactor
        # Examples: "P00918", "G3XD94", "Q90240", "P38289", "P9WP55", "B3QM53", "O15304", "P26662", "A9CEQ7", and "Q9F1R6"
        elif cmt["commentType"] == "COFACTOR":
            # Get the cofactors of the protein
            cofs = cmt.get("cofactors", [])
            if cofs:
                list_of_cofactors = [cof["name"] for cof in cofs]
                cofactors = natural_join(list_of_cofactors, "and")
            else:
                list_of_cofactors = None
                cofactors = None

            # Get the notes on the cofactors of the protein
            txts = cmt.get("note", {}).get("texts", [])
            if txts:
                notes = [txt["value"] for txt in txts]
                if len(notes) == 1:
                    notes = notes[0]
                elif len(notes) > 1:
                    warn_general(protein_id, field, i, cmt["commentType"], f'["note"]["texts"] has a length of {str(len(notes))} (expected 1)')
                    continue
            else:
                notes = None

            record = {
                "Cofactors": cofactors,
                "ListOfCofactors": list_of_cofactors,
                "Notes": notes
            }

            cofactors_collection.append(record)


        # Subcellular Location
        # Ref: UniProt Help / Subcellular location: https://www.uniprot.org/help/subcellular_location
        # Examples: "Q9ULX6", "Q02630", "P38182", "O43687"
        elif cmt["commentType"] == "SUBCELLULAR LOCATION":
            # Get the name of the molecule (alternative product / isoform) if present
            name = cmt.get("molecule", "Canonical")

            # Get the subcellular locations of the protein
            sublocs = cmt.get("subcellularLocations", [])
            if sublocs:
                subcellular_locations_local = []
                for j, subloc in enumerate(sublocs):
                    location = subloc["location"]["value"] if "location" in subloc else None
                    topology = subloc["topology"]["value"] if "topology" in subloc else None
                    orientation = subloc["orientation"]["value"] if "orientation" in subloc else None

                    if location and not orientation and not topology:
                        subcellular_location_local = location
                    elif location and topology and not orientation:
                        subcellular_location_local = f"{location} ({topology})"
                    elif location and topology and orientation:
                        subcellular_location_local = f"{location} ({topology}; {orientation})"
                    elif location and orientation and not topology:
                        subcellular_location_local = f"{location} ({orientation})"
                    else:
                        warn_unexpected_key_missing(protein_id, field, i, cmt["commentType"], f'["subcellularLocations"][{str(j)}]["location"]')
                        continue

                    subcellular_locations_local.append(subcellular_location_local)

                subcellular_locations = natural_join(subcellular_locations_local, "and")
            else:
                warn_unexpected_key_missing(protein_id, field, i, cmt["commentType"], '["subcellularLocations"]') # Only "P07884", "Q8N884"
                continue

            # Get the notes on the subcellular locations of the protein
            txts = cmt.get("note", {}).get("texts", [])
            if txts:
                notes = [txt["value"] for txt in txts]
                if len(notes) == 1:
                    notes = notes[0]
                elif len(notes) > 1:
                    warn_general(protein_id, field, i, cmt["commentType"], f'["note"]["texts"] has a length of {str(len(notes))} (expected 1)')
                    continue
            else:
                notes = None

            subcellular_location = {
                "Name": name,
                "SubcellularLocations": subcellular_locations,
                "Notes": notes
            }

            subcellular_location_collection.append(subcellular_location)


        # Pathway
        # Ref: UniProt Help / Pathway: https://www.uniprot.org/help/pathway
        # Examples: "P77731", "P20580", "Q47316", and "P00562"
        elif cmt["commentType"] == "PATHWAY":
            # Get the pathways the protein involved in
            txts = cmt.get("texts", [])
            if txts:
                pwys = [txt["value"] for txt in txts]
                if len(pwys) == 1:
                    pwys = pwys[0]

                    pathway_step = pwys.split(": ")
                    hierarchical_pathway = pathway_step[0]
                    if len(pathway_step) == 2:
                        step = pathway_step[1]
                        superpathway_pathway_subpathway = hierarchical_pathway.split("; ")
                        if len(superpathway_pathway_subpathway) == 3:
                            superpathway = superpathway_pathway_subpathway[0]
                            pathway = superpathway_pathway_subpathway[1]
                            subpathway = superpathway_pathway_subpathway[2]
                        else:
                            warn_general(protein_id, field, i, cmt["commentType"], f'["texts"][0]["value"] has an unexpected format when step exists: {hierarchical_pathway}')
                            continue

                    elif len(pathway_step) == 1:
                        step = None
                        superpathway_pathway_subpathway = hierarchical_pathway.split("; ")
                        if len(superpathway_pathway_subpathway) == 2:   # Commonest
                            superpathway = superpathway_pathway_subpathway[0]
                            pathway = superpathway_pathway_subpathway[1]
                            subpathway = None
                        elif len(superpathway_pathway_subpathway) == 1: # Rarely appear
                            superpathway = superpathway_pathway_subpathway[0]
                            pathway = None
                            subpathway = None
                        elif len(superpathway_pathway_subpathway) == 3: # Only "Q94JQ4"
                            superpathway = superpathway_pathway_subpathway[0]
                            pathway = None
                            subpathway = None
                        else:
                            warn_general(protein_id, field, i, cmt["commentType"], f'["texts"][0]["value"] has an unexpected format when step does not exist: {hierarchical_pathway}')
                            continue

                    else:
                        warn_general(protein_id, field, i, cmt["commentType"], f'["texts"][0]["value"] has an unexpected format: {pwys}')
                        continue

                    record = {
                        "SuperPathway": superpathway,
                        "Pathway": pathway,
                        "SubPathway": subpathway,
                        "Step": step,
                    }

                    pathway_collection.append(record)

                else:
                    warn_general(protein_id, field, i, cmt["commentType"], f'["texts"] has a length of {str(len(pwys))} (expected 1)')
                    continue
            else:
                warn_unexpected_key_missing(protein_id, field, i, cmt["commentType"], '["texts"]')
                continue


    info = {
        "CatalyticActivity": catalytic_activity_collection,
        "Cofactors": cofactors_collection,
        "SubcellularLocation": subcellular_location_collection,
        "Pathway": pathway_collection,
    }

    return info


def extract_uniprot_cross_references(protein_id):
    data = get_protein_info(protein_id)

    field = "uniProtKBCrossReferences"

    go_bps = []
    go_mfs = []
    go_ccs = []
    for i, cref in enumerate(data.get("uniProtKBCrossReferences", [])):
        # Gene Ontology (GO)
        # Ref: UniProt Help / Gene Ontology (GO): https://www.uniprot.org/help/gene_ontology
        # Examples: "P04637", "P00533", "P68871", "P69905", "P60709", "P07437", "P04406", "P99999", "P01308", "P00722", "P0A8V2", and "P02144"
        if cref["database"] == "GO":
            # Get the GO ID of the protein
            go_id = cref.get("id")
            if not go_id:
                warn_unexpected_key_missing(protein_id, field, i, cref["database"], '["id"]')
                continue

            # Get the GO term of the protein
            props = cref.get("properties", [])
            if len(props) == 2:
                value = []
                for prop in props:
                    if prop["key"] == "GoTerm":
                        value.append(prop["value"])

                if len(value) != 1:
                    warn_general(protein_id, field, i, cref["database"], f'["properties"] contains {str(len(value))} attributes with ["key"] "GoTerm" (expected 1)')
                    continue

                value = value[0]
                go_term = value[2:]
                go_sample = f"{go_term} ({go_id})"
                if value[0] == "P":
                    go_bps.append(go_sample)
                elif value[0] == "F":
                    go_mfs.append(go_sample)
                elif value[0] == "C":
                    go_ccs.append(go_sample)
                else:
                    warn_general(protein_id, field, i, cref["database"], f'["value"] annotates an invalid GO type: {value}')
                    continue

            else:
                warn_general(protein_id, field, i, cref["database"], f'Unexpected number of ["properties"]: {len(props)}')
                continue

    info = {
        "GoBP": go_bps,
        "GoMF": go_mfs,
        "GoCC": go_ccs
    }

    return info


# Sequence Annotation
# Ref: UniProt Help / Sequence annotation (Features): https://www.uniprot.org/help/sequence_annotation
# Examples: "P04637", "P38398", "Q14596"
def extract_uniprot_features(protein_id):
    feature_classification = {
        "MoleculeProcessing": ["Initiator methionine", "Signal", "Transit peptide",
                               "Propeptide", "Chain", "Peptide"],
        "Regions": ["Topological domain", "Transmembrane", "Intramembrane",
                    "Domain", "Repeat", "Calcium binding",
                    "Zinc finger", "DNA binding", "Nucleotide binding",
                    "Region", "Coiled coil", "Motif",
                    "Compositional bias"],
        "Sites": ["Active site", "Metal binding", "Site"], # "Metal binding" not found
        "AminoAcidModifications": ["Non-standard residue", "Modified residue", "Lipidation",
                                   "Glycosylation", "Disulfide bond", "Cross-link"],
        "NaturalVariations": ["Alternative sequence", "Natural variant"],
        "ExperimentalInfo": ["Mutagenesis", "Sequence uncertainty", "Sequence conflict",
                             "Non-adjacent residues", "Non-terminal residue"],
        "SecondaryStructure": ["Helix", "Turn", "Beta strand"]
    }
    map_to_the_content = {
        # Molecule processing
        "Initiator methionine": "the initiator methionine cleavage",
        "Signal": "a signal sequence that targets proteins to the secretory pathway or periplasmic space",
        "Transit peptide": "a transit peptide for organelle targeting",
        "Propeptide": "a propeptide cleaved during protein maturation or activation",
        "Chain": "a mature polypeptide chain",
        "Peptide": "an active peptide in the mature protein",

        # Regions
        "Topological domain": "non-membrane regions of membrane-spanning proteins",
        "Transmembrane": "a membrane-spanning region",
        "Intramembrane": "a region located within the membrane without crossing it",
        "Domain": "a modular protein domain",
        "Repeat": "repeated sequence motifs or domains",
        "Calcium binding": "regions that coordinate calcium ions",
        "Zinc finger": "zinc finger motifs within the protein",
        "DNA binding": "a DNA-binding domain",
        "Nucleotide binding": "nucleotide-binding regions such as ATP or GTP sites",
        "Region": "a region of interest in the protein sequence",
        "Coiled coil": "coiled-coil regions within the protein",
        "Motif": "short sequence motifs of biological interest",
        "Compositional bias": "a region of compositional bias in the protein",

        # Sites
        "Active site": "amino acid residues directly involved in enzymatic activity",
        "Binding site": "a binding site for a chemical group such as a coenzyme or prosthetic group",
        "Metal binding": "metal-binding sites that coordinate metal ions",
        "Site": "a specific amino acid site of interest within the sequence",

        # Amino acid modifications
        "Non-standard residue": "non-standard amino acids such as selenocysteine or pyrrolysine",
        "Modified residue": "post-translationally modified residues excluding lipids, glycans, and cross-links",
        "Lipidation": "covalently attached lipid groups",
        "Glycosylation": "covalently attached glycan groups",
        "Disulfide bond": "cysteine residues forming disulfide bonds",
        "Cross-link": "covalently linked residues between proteins",

        # Natural variations
        "Alternative sequence": "amino acid variants producing alternate protein isoforms",
        "Natural variant": "natural amino acid variants of the protein",

        # Experimental info
        "Mutagenesis": "experimentally altered sites created by mutagenesis",
        "Sequence uncertainty": "regions of sequence uncertainty",
        "Sequence conflict": "sequence discrepancies of unknown origin",
        "Non-adjacent residues": "non-consecutive residues within the sequence",
        "Non-terminal residue": "an internal residue not located at the protein terminus",

        # Secondary structure
        "Helix": "helical regions within the experimentally determined structure",
        "Turn": "turn regions within the experimentally determined structure",
        "Beta strand": "beta-strand regions within the experimentally determined structure"
    }
    data = get_protein_info(protein_id)
    field = "features"

    # "When a feature is known to extend beyond the position that is given in this section,
    # the endpoint specification will be preceded by '<' (less than)
    # for features which continue to the N-terminal direction
    # or by '>' (greater than) for features which continue to the C-terminal direction.
    # Example: 'P62756'
    # Unknown endpoints are denoted by a question mark '?'.
    # Example: 'P78586'
    # Uncertain endpoints are denoted by a question mark '?' before the position, e.g. '?42'.
    # Example: 'Q3ZC31'"
    def extract_location(ft):
        return {
            "StartLocation": str(ft["location"]["start"]["value"]),
            "EndLocation": str(ft["location"]["end"]["value"])
        }

    def extract_description(ft):
        return ft.get("description", "")

    def extract_ligand_and_ligand_part(ft):
        ligand = ft.get("ligand", {}).get("name", "")
        ligand_note = ft.get("ligand", {}).get("note", "")
        ligand_part = ft.get("ligandPart", {}).get("name", "")
        if ligand:
            if ligand_part and ligand_note:
                ligand = f"{ligand} ({ligand_part}) ({ligand_note})"
            elif ligand_part and not ligand_note:
                ligand = f"{ligand} ({ligand_part})"
            elif not ligand_part and ligand_note:
                ligand = f"{ligand} ({ligand_note})"
        return ligand

    def extract_alternative_sequence(ft):
        return {
            "OriginalSequence": ft.get("alternativeSequences", {}).get("originalSequence", ""),
            "AlternativeSequences": natural_join(ft.get("alternativeSequences", {}).get("alternativeSequences", ""), "or")
        }

    molecule_processing_collection = []
    regions_collection = []
    sites_collection = []
    binding_sites_collection = []
    amino_acid_modifications_collection = []
    natural_variations_collection = []
    experimental_info_collection = []

    for i, ftr in enumerate(data.get("features", [])):
        field_type = ftr["type"]

        # Molecule processing, Regions, Sites (exclude Binding Site) & Amino acid modifications: Only description
        if field_type in feature_classification["Regions"] or field_type in feature_classification["Sites"]\
                or field_type in feature_classification["MoleculeProcessing"]\
                or field_type in feature_classification["AminoAcidModifications"]:
            record = {
                "FeatureType": field_type,
                "Content": map_to_the_content[field_type],
                "StartLocation": extract_location(ftr)["StartLocation"],
                "EndLocation": extract_location(ftr)["EndLocation"],
                "Description": extract_description(ftr),
            }

        # Binding site: Description & Ligand
        elif field_type == "Binding site":
            record = {
                "FeatureType": field_type,
                "Content": map_to_the_content[field_type],
                "StartLocation": extract_location(ftr)["StartLocation"],
                "EndLocation": extract_location(ftr)["EndLocation"],
                "Description": extract_description(ftr),
                "Ligand": extract_ligand_and_ligand_part(ftr)
            }

        # Natural variations & Experimental info: Description & Alternative Sequence
        elif field_type in feature_classification["NaturalVariations"] or field_type in feature_classification["ExperimentalInfo"]:
            record = {
                "FeatureType": field_type,
                "Content": map_to_the_content[field_type],
                "StartLocation": extract_location(ftr)["StartLocation"],
                "EndLocation": extract_location(ftr)["EndLocation"],
                "Description": extract_description(ftr),
                "OriginalSequence": extract_alternative_sequence(ftr)["OriginalSequence"],
                "AlternativeSequences": extract_alternative_sequence(ftr)["AlternativeSequences"]
            }

        # Secondary structure:
        elif field_type in feature_classification["SecondaryStructure"]:
            continue

        else:
            warn_general(protein_id, field, i, "", f"Invalid feature type: {field_type}")
            continue

        if field_type in feature_classification["MoleculeProcessing"]:
            molecule_processing_collection.append(record)
        elif field_type in feature_classification["Regions"]:
            regions_collection.append(record)
        elif field_type in feature_classification["Sites"]:
            sites_collection.append(record)
        elif field_type == "Binding site":
            binding_sites_collection.append(record)
        elif field_type in feature_classification["AminoAcidModifications"]:
            amino_acid_modifications_collection.append(record)
        elif field_type in feature_classification["NaturalVariations"]:
            natural_variations_collection.append(record)
        elif field_type in feature_classification["ExperimentalInfo"]:
            experimental_info_collection.append(record)

    info = {
        "MoleculeProcessing": molecule_processing_collection,
        "Regions": regions_collection,
        "Sites": sites_collection,
        "BindingSite": binding_sites_collection,
        "AminoAcidModifications": amino_acid_modifications_collection,
        "NaturalVariations": natural_variations_collection,
        "ExperimentalInfo": experimental_info_collection
    }

    return info