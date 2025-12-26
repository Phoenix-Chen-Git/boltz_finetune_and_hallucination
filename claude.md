# Boltz-2 Finetuning Project Summary

This document summarizes the current state of the environment, data, and the planned workflow for finetuning Boltz-2 on the steroid hormone sensor dataset.

## 1. Environment State
- **Micromamba Environment**: `boltz` (Python 3.11)
- **Installation**: `boltz` package installed in editable mode with `[cuda]` extras.
- **Dependencies**: `openpyxl`, `redis`, `mmseqs2` (partially configured).
- **GPU**: NVIDIA A800-SXM4-80GB.

## 2. Model Assets & Cache (`~/.boltz`)
All required inference assets are downloaded and verified:
- `boltz2_conf.ckpt`: Structure prediction weights (~3GB).
- `boltz2_aff.ckpt`: Affinity prediction weights (~2.5GB).
- `mols/`: Extracted molecule database.
- `ccd.pkl`: Chemical Component Dictionary.

## 3. Data Analysis: Steroid Hormone Sensor Screening
- **File**: `finetune/Steroid hormone sensor screening_to Baker lab.xlsx`
- **Content**:
    - Hundreds of protein sequences (variants/mutants).
    - Normalized affinity values for: Estradiol, Progesterone, Corticosterone, Testosterone.
- **Goal**: Finetune Boltz-2 to improve its affinity prediction accuracy for these specific sensors and ligands.

## 4. Current Progress
- Successfully ran structure inference on `examples/prot_no_msa.yaml`.
- Verified that the model can load weights and process sequences without MSA.
- Data extraction logic for the Excel file has been tested.

## 5. Planned Finetuning Workflow
To finetune the model, we need to transform the Excel data into the Boltz training format:

1. **Structure Generation**: Generate 3D structures for the sequences in the Excel file (since they are variants, we can use Boltz-2 inference or AlphaFold).
2. **MSA Generation**: Create Multiple Sequence Alignments (`.a3m` files) for the unique sequences.
3. **Data Preprocessing**:
    - Convert sequences and structures to `.npz` format.
    - Create a `manifest.json` mapping sequences, structures, MSAs, and target affinity values.
4. **Training Configuration**:
    - Customize `scripts/train/configs/full.yaml` to focus on the affinity loss.
    - Set the `boltz2_aff.ckpt` as the starting checkpoint.
5. **Execution**: Run `scripts/train/train.py` using the processed dataset.

## 6. Precise Next Steps: Data Preparation Script

We will implement `scripts/prepare_finetune_dataset.py` to automate the transformation of the Excel data.

### Step 6.1: Ligand Mapping
We will use the following SMILES strings for the four hormones:
- **Estradiol**: `CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O`
- **Progesterone**: `CC(=O)C1CCC2C1(CCC3C2CCC4=CC(=O)CCC34C)C`
- **Corticosterone**: `CC12CCC3C(C1CCC2(C(=O)CO)O)CCC4=CC(=O)CCC34C`
- **Testosterone**: `CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C`

### Step 6.2: Sequence Extraction & MSA Generation
1.  **Unique Sequences**: Extract all unique sequences from the four sheets.
2.  **Local MSA**: Run `mmseqs2` against a sequence database (e.g., UniRef100) to generate `.a3m` files for each unique sequence.
    - Files will be saved as `<sequence_hash>.a3m`.
3.  **Caching**: Use a local cache to avoid re-generating MSAs for identical sequences across different hormones.

### Step 6.3: Structure Initialization (Complex Generation)
Since we need 3D complexes for affinity prediction:
1.  **Reference Prediction**: Run `boltz predict` for one representative variant per hormone to get a high-confidence binding pose.
2.  **Variant Modeling**: For other variants, we will either:
    - **A**: Run Boltz inference for all (most accurate but slow).
    - **B**: Use the reference pose and swap the sequence (fast, assuming mutations don't shift the ligand significantly).

### Step 6.4: Training Manifest Assembly
Create a `manifest.json` where each entry contains:
- `id`: Unique identifier (e.g., `Estradiol_Variant_1`).
- `protein_sequence`: From Excel.
- `ligand_smiles`: From our mapping.
- `structure_path`: Path to the `.npz` file.
- `msa_path`: Path to the `.a3m` file.
- `target_affinity`: The `Norm aff` value from Excel (mapped to `log10(IC50)` if necessary).

## 7. Configuration for Finetuning
We will create a new config `scripts/train/configs/affinity_finetune.yaml` that:
- Loads `boltz2_aff.ckpt`.
- Sets the `affinity_prediction` flag to `True`.
- Points to the custom dataset directory.
- Minimizes the structure/diffusion weights to focus learning on the affinity head.

