# Pannot-Data: A Multi-task, High-resolution Instruction Dataset for Text-based Protein Understanding

## Versioning
The current release is **v1.0**

## Release History

- **v1.0** (Jan 2026)
  - Initial public release of the Pannot-Data dataset construction pipeline
  - Instruction samples constructed via a hybrid human–LLM workflow, where data are retrieved from UniProtKB/Swiss-Prot, processed into a template-ready format, and instantiated through curated instruction templates
  - Coverage of eight core task families: catalytic activity, cofactor, pathway, subcellular location, gene ontology, EC numbers, regions, and sites
  - Unified CLI and runner for reproducible dataset generation

## Environment Setup

### Clone the Repository
```bash
git clone --branch develop https://github.com/Antoninnnn/Pannot-Data.git
cd Pannot-Data
```

### Prerequisites
- Python ≥ 3.10
- Conda (recommended)

### Create Environment
```bash
conda create -n pannot-data python=3.10
conda activate pannot-data
```
If Conda is not installed, please install it before proceeding.

### Install Dependencies for Dataset Construction
```bash
pip install -e .
```
This installs the dependencies required for dataset construction. The current pipeline only relies on the `requests` package and does not require any machine learning frameworks.

### Verify Installation
```bash
python -c "import Composer; print('Environment setup successful')"
```
If this command runs without errors and prints the message above, the environment setup is successful. Otherwise, an import error indicates an incomplete installation.

## Dataset Construction

### Overview

Pannot-Data is constructed by retrieving curated protein entries from UniProtKB/Swiss-Prot in JSON format via the UniProtKB REST API, and transforming structured annotations into instruction samples through a hybrid human–LLM co-authored instruction template framework. The current pipeline focuses on several task families derived from specific UniProtKB JSON blocks:

- **Comments-derived tasks** (from the UniProtKB JSON `comments` block):
  - Catalytic activity
  - Cofactor
  - Pathway
  - Subcellular location

- **Gene Ontology (GO)** (from the UniProtKB JSON `uniProtKBCrossReferences` block):
  - GO biological process (BP), molecular function (MF), cellular component (CC)

- **Residue-/region-level annotations** (from the UniProtKB JSON `features` block):
  - Regions (e.g., domains/regions)
  - Sites (residue-level site annotations)

- **EC numbers** (from the UniProtKB JSON `proteinDescription` → `recommendedName` → `ecNumbers`):
  - Enzyme Commission (EC) annotations associated with the recommended protein name

Since `comments`, `uniProtKBCrossReferences`, and `features` are list-based blocks with different schemas, the pipeline groups extraction and sample generation by these blocks (plus EC fields) to minimize redundant traversals and enable independent generators for each task family.

### Sample Generation

Instruction samples are generated via a unified runner that supports multiple task families in a single invocation. Users specify the target tasks, input protein IDs, and output configuration through command-line arguments.

A minimal example is shown below:
```bash
bash scripts/run.sh \
  --tasks comments cross_references ec features \
  --ids_rel demo_ids.txt \
  --out_subdir demo_samples \
  --template_subdir Template \
  --batch_size 100 \
  --output_format jsonl
```

This command retrieves UniProtKB/Swiss-Prot entries for the specified protein IDs and generates instruction samples for all supported tasks using a demo ID list containing 500 accession IDs.

**Key arguments:**

- `--tasks`: One or more task families to generate. Supported values include:
  - `comments`: Catalytic activity, cofactor, pathway, and subcellular location
  - `cross_references`: Gene Ontology (GO) annotations
  - `ec`: Enzyme Commission (EC) number annotations
  - `features`: Region- and site-level annotations

- `--ids_rel`: Path to a text file containing UniProt accession IDs (one per line), specified relative to the project data directory.

- `--out_subdir`: Name of the output subdirectory under the dataset output root (default: empty, outputs are written directly to `data/Sample`).

- `--template_subdir`: Name of the template directory under the project data directory (default: `Template`).

- `--batch_size`: Number of proteins processed per API request batch (default: `100`).

- `--output_format`: Output format for generated samples. Supported values: `jsonl`, `tsv`, `csv` (default: `jsonl`).

**Note:** The `--tasks` and `--ids_rel` arguments are required and have no default values; at least one task must be specified for the program to run.

You can use the following command to generate instruction samples across all supported tasks for Pannot pretraining:
```bash
bash scripts/run.sh \
  --tasks comments cross_references ec features \
  --ids_rel splits/train_condensed_ids.txt \
  --out_subdir train_condensed_samples \
  --template_subdir Template_pt \
  --batch_size 100 \
  --output_format tsv
```
This command generates the training samples using `splits/train_condensed_ids.txt`; to generate the full set of evaluation samples, you can replace `--ids_rel` with `splits/test_id_ids.txt`, `splits/test_ood_family_ids.txt`, `splits/val_id_ids.txt`, and `splits/val_ood_family_ids.txt`, and correspondingly set `--out_subdir` to `test_id_samples`, `test_ood_family_samples`, `val_id_samples`, and `val_ood_family_samples`, respectively.



