# Pannot-Data: A Multi-task, High-resolution Instruction Dataset for Text-based Protein Understanding

## Versioning
The current release is **v1.0**

## Release History

- **v1.0** (Jan 2026)
  - 1213

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

Pannot-Data is constructed by retrieving curated protein entries from UniProtKB/Swiss-Prot in JSON format via the UniProtKB REST API, and transforming structured annotations into multi-task, high-resolution instruction samples. The current pipeline focuses on several task families derived from specific UniProtKB JSON blocks:

- **Comments-derived tasks** (from the `comments` block; a list of typed comment objects):
  - Catalytic activity
  - Cofactor
  - Pathway
  - Subcellular location

- **Gene Ontology (GO)** (from the `uniProtKBCrossReferences` block; a list of cross-reference records):
  - GO biological process (BP), molecular function (MF), cellular component (CC)

- **Residue-/region-level annotations** (from the `features` block; a list of sequence feature records):
  - Regions (e.g., domains/regions)
  - Sites (residue-level site annotations)

- **EC numbers** (from `proteinDescription` → `recommendedName` → `ecNumbers`):
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
  --batch_size 100 \
  --output_format tsv
```

This command retrieves UniProtKB/Swiss-Prot entries for the specified protein IDs and generates instruction samples for the selected tasks.

```md
**Key arguments:**

- `--tasks`: One or more task families to generate. Supported values include:
  - `comments`: Catalytic activity, cofactor, pathway, and subcellular location
  - `cross_references`: Gene Ontology (GO) annotations
  - `ec`: Enzyme Commission (EC) number annotations
  - `features`: Region- and site-level annotations

- `--ids_rel`: Path to a text file containing UniProt accession IDs (one per line), specified relative to the project data directory.

- `--out_subdir`: Name of the output subdirectory under the dataset output root.

- `--batch_size`: Number of proteins processed per API request batch.

- `--output_format`: Output format for generated samples (e.g., `tsv` or `csv`).
```
