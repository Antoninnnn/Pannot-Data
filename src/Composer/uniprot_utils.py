import requests
import warnings

from .constants.feature_classification import FEATURE_CLASSIFICATION_MAP


def get_protein_info(protein_id):
    """Retrieve the full UniProtKB JSON record for a given protein ID."""
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
    response = requests.get(url)
    response.raise_for_status()
    return response.json()


def get_proteins_info(protein_ids):
    """Retrieve full UniProtKB JSON records for a list of protein accession IDs."""
    query = " OR ".join(f"accession:{pid}" for pid in protein_ids)
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "json",
        "size": len(protein_ids),
    }
    r = requests.get(url, params=params)
    r.raise_for_status()
    return r.json()["results"]


def safe_get(d, keys, default=None):
    """Safely retrieve a nested value from a UniProtKB JSON dictionary."""
    for k in keys:
        if isinstance(d, dict) and k in d:
            d = d[k]
        else:
            return default
    return d


def safe_list(lst, key):
    """Extract a list of values associated with a given key from a list of dictionaries."""
    return [item.get(key) for item in lst if key in item]


def extract_uniprot_meta(data):
    """Extract core metadata from a UniProtKB entry given a protein JSON record."""
    ec_numbers = [ec["value"] for ec in safe_get(data, ["proteinDescription", "recommendedName", "ecNumbers"], [])]

    info = {
        "Accession": data.get("primaryAccession"),
        "EntryName": data.get("uniProtkbId"),
        "Organism": safe_get(data, ["organism", "scientificName"]),
        "TaxonomyId": safe_get(data, ["organism", "taxonId"]),
        "Existence": data.get("proteinExistence"),
        "Sequence": safe_get(data, ["sequence", "value"]),
        "Length": safe_get(data, ["sequence", "length"]),
        "MolecularWeight": safe_get(data, ["sequence", "molWeight"]),
        "ECNumber": ec_numbers,
        "IsEnzyme": bool(ec_numbers),
    }

    return info


def extract_uniprot_comments(data):
    """Extract structured functional annotations from UniProt comments blocks."""

    def extract_comments_note(comment):
        """Extract free-text notes from a UniProt comment block."""
        return safe_list(safe_get(comment, ["note", "texts"], []), "value")

    catalytic_reactions = []
    cofactors = []
    subcellular_locations = []
    pathways = []

    for i, cmt in enumerate(data.get("comments", [])):
        # Catalytic Activity
        # Ref: UniProt Help / Catalytic activity: https://www.uniprot.org/help/catalytic_activity
        if cmt["commentType"] == "CATALYTIC ACTIVITY":
            # Get the catalytic reaction of the protein
            rxn = safe_get(cmt, ["reaction", "name"])

            if rxn:
                catalytic_reactions.append(rxn)
            else:
                warnings.warn(f"Protein {data.get('primaryAccession')}: empty CATALYTIC ACTIVITY comment in comment #{i}")

        # Cofactor
        # Ref: UniProt Help / Cofactor: https://www.uniprot.org/help/cofactor
        elif cmt["commentType"] == "COFACTOR":
            # Get the cofactors of the protein
            cofs = safe_list(cmt.get("cofactors", []), "name")

            # Get the notes on the cofactors of the protein
            cofs_notes = extract_comments_note(cmt)

            if cofs or cofs_notes:
                cofactors.append({
                    "Cofactors": cofs,      # list[str]: cofactor names
                    "Notes": cofs_notes,    # list[str]: notes on the cofactors
                })
            else:
                warnings.warn(f"Protein {data.get('primaryAccession')}: empty COFACTOR comment in comment #{i}")

        # Subcellular Location
        # Ref: UniProt Help / Subcellular location: https://www.uniprot.org/help/subcellular_location
        elif cmt.get("commentType") == "SUBCELLULAR LOCATION":
            # Get the name of the molecule (alternative product / isoform) if present
            name = cmt.get("molecule")

            # Get the subcellular locations of the protein
            sublocs = []
            for entry in cmt.get("subcellularLocations", []):
                location = safe_get(entry, ["location", "value"])
                orientation = safe_get(entry, ["orientation", "value"])
                topology = safe_get(entry, ["topology", "value"])

                if location or orientation or topology:
                    sublocs.append({
                        "Location": location,
                        "Orientation": orientation,
                        "Topology": topology,
                    })

            # Get the notes on the subcellular locations of the protein
            sublocs_notes = extract_comments_note(cmt)

            if sublocs or sublocs_notes or name:
                subcellular_locations.append({
                    "Name": name,                       # str | None: alternative product / isoform name
                    "SubcellularLocations": sublocs,    # list[dict]: each entry describes one subcellular location
                    "Notes": sublocs_notes,             # list[str]: notes on subcellular location
                })
            else:
                warnings.warn(f"Protein {data.get('primaryAccession')}: empty SUBCELLULAR LOCATION comment in comment #{i}")

        # Pathway
        # Ref: UniProt Help / Pathway: https://www.uniprot.org/help/pathway
        elif cmt["commentType"] == "PATHWAY":
            # Get the pathways the protein involved in
            pwys = safe_list(cmt.get("texts", []), "value")

            if pwys:
                pathways.append(pwys)
            else:
                warnings.warn(f"Protein {data.get('primaryAccession')}: empty PATHWAY comment in comment #{i}")

    info = {
        "CatalyticActivity": catalytic_reactions,
        "Cofactor": cofactors,
        "SubcellularLocation": subcellular_locations,
        "Pathway": pathways,
    }

    return info


def extract_uniprot_cross_references(data):
    """Extract structured functional annotations from UniProt cross-referenced databases."""
    gene_ontologies = []

    for idx, cref in enumerate(data.get("uniProtKBCrossReferences", [])):
        # Gene Ontology (GO)
        # Ref: UniProt Help / Gene Ontology (GO): https://www.uniprot.org/help/gene_ontology
        if cref["database"] == "GO":
            # Get the GO ID of the protein
            go_id = cref.get("id")

            # Get the GO term of the protein
            go_terms = [
                entry["value"]
                for entry in cref.get("properties", [])
                if entry.get("key") == "GoTerm"
            ]

            if len(go_terms) != 1:
                warnings.warn(
                    f"Protein {data.get('primaryAccession')}: invalid GO annotation "
                    f"in cross-references #{idx}"
                )
                continue

            go_term = go_terms[0]

            if go_id and go_term:
                gene_ontologies.append({
                    "ID": go_id,          # str: Gene Ontology identifier
                    "Term": go_term,      # str: Gene Ontology term (with namespace prefix)
                })
            else:
                warnings.warn(
                    f"Protein {data.get('primaryAccession')}: incomplete GO annotation "
                    f"in cross-reference #{idx}"
                )

    info = {
        "GeneOntology": gene_ontologies
    }

    return info


# Sequence Annotation
# Ref: UniProt Help / Sequence annotation (Features): https://www.uniprot.org/help/sequence_annotation
def extract_uniprot_features(data):
    """Extract structured functional annotations from UniProt features blocks."""
    regions = []
    sites = []

    # --- Reserved categories (kept intentionally for future use) ---
    # molecule_processing = []
    # amino_acid_modifications = []
    # natural_variations = []
    # experimental_info = []
    # secondary_structure = []

    for i, ftr in enumerate(data.get("features", [])):
        feature_type = ftr["type"]

        basic_record = {
            "Type": feature_type,
            "Start": safe_get(ftr, ["location", "start", "value"]),
            "StartModifier": safe_get(ftr, ["location", "start", "modifier"]),
            "End": safe_get(ftr, ["location", "end", "value"]),
            "EndModifier": safe_get(ftr, ["location", "end", "modifier"]),
            "Description": ftr.get("description"),
        }

        # Regions
        if feature_type in FEATURE_CLASSIFICATION_MAP["Regions"]:
            regions.append(basic_record)

        # Sites
        elif feature_type in FEATURE_CLASSIFICATION_MAP["Sites"]:
            sites.append(
                {
                    **basic_record,
                    "Ligand": safe_get(ftr, ["ligand", "name"]),
                    "LigandNote": safe_get(ftr, ["ligand", "note"]),
                    "LigandPart": safe_get(ftr, ["ligandPart", "name"]),
                    "LigandPartNote": safe_get(ftr, ["ligandPart", "note"]),
                }
            )

        # Molecule processing (reserved)
        # elif feature_type in FEATURE_CLASSIFICATION_MAP["MoleculeProcessing"]:
        #     molecule_processing.append(basic_record)

        # Amino acid modifications (reserved)
        # elif feature_type in FEATURE_CLASSIFICATION_MAP["AminoAcidModifications"]:
        #     amino_acid_modifications.append(basic_record)

        # Natural variations (reserved)
        # elif feature_type in FEATURE_CLASSIFICATION_MAP["NaturalVariations"]:
        #     natural_variations.append(
        #         {
        #             **basic_record,
        #             "OriginalSequence": safe_get(
        #                 ftr, ["alternativeSequence", "originalSequence"]
        #             ),
        #             "AlternativeSequences": safe_get(
        #                 ftr, ["alternativeSequence", "alternativeSequences"]
        #             ),  # list[str] | None
        #         }
        #     )

        # Experimental info (reserved)
        # elif feature_type in FEATURE_CLASSIFICATION_MAP["ExperimentalInfo"]:
        #     experimental_info.append(
        #         {
        #             **basic_record,
        #             "OriginalSequence": safe_get(
        #                 ftr, ["alternativeSequence", "originalSequence"]
        #             ),
        #             "AlternativeSequences": safe_get(
        #                 ftr, ["alternativeSequence", "alternativeSequences"]
        #             ),
        #         }
        #     )

        # Secondary structure (reserved)
        # elif feature_type in FEATURE_CLASSIFICATION_MAP["SecondaryStructure"]:
        #     secondary_structure.append(basic_record)

    info = {
        "Regions": regions,
        "Sites": sites,
        # Reserved for future:
        # "MoleculeProcessing": molecule_processing,
        # "AminoAcidModifications": amino_acid_modifications,
        # "NaturalVariations": natural_variations,
        # "ExperimentalInfo": experimental_info,
        # "SecondaryStructure": secondary_structure,
    }

    return info