# 🦠 FluRisk: Influenza Virus Risk Assessment Tool

[![Python](https://img.shields.io/badge/Python-3.6+-blue.svg)](https://python.org/) [![Diamond](https://img.shields.io/badge/Diamond-2.1.11-orange.svg)](https://github.com/bbuchfink/diamond)

`FluRisk` (`flupre`) is a comprehensive risk assessment tool for influenza A viruses. The tool identifies amino acid molecular markers associated with five phenotypic traits and calculates risk scores for early warning and risk assessment. These phenotypes include human adaptation, transmissibility, receptor binding ability, virulence, and drug resistance. Additionally, the tool provides automatic sequence annotation, antigenic relationship assessment between viruses and vaccines, and prediction of host specificity, pathogenicity, and receptor binding preferences.

## 📋 Overview

The `FluRisk` pipeline consists of six primary modules:

| Module                     | Function   | Description                                                  |
| -------------------------- | ---------- | ------------------------------------------------------------ |
| 🏷️ **Annotation**           | `anno`     | Sequence annotation, antigenic relationship assessment, and reassortment risk evaluation |
| 🔍 **Extraction**           | `extract`  | Molecular marker identification for five phenotypic traits   |
| 🎯 **Host Prediction**      | `predh`    | Host classification (avian/human) based on adaptation markers |
| ⚡ **Virulence Prediction** | `predv`    | Virulence level prediction (pathogenic/non-pathogenic)       |
| 🔗 **Binding Prediction**   | `predr`    | Receptor binding preference prediction (α2,3-linked sialic acid/α2,6-linked sialic acid) |
| 📊 **Risk Assessment**      | `risk2tab` | Comprehensive risk scoring for all phenotypic traits         |

**Input Format**: Individual FASTA files containing all sequences of a single influenza virus strain, or directories where each file represents sequences from separate influenza virus strains.

------

## ⚙️ Prerequisites

### System Dependencies

#### 💎 Diamond Installation

Diamond is a high-performance sequence alignment tool required for protein sequence comparison. Install using one of the following methods:

**Conda Installation (Recommended):**

```bash
conda install -c bioconda diamond=2.1.11
```

**Direct Binary Download:** Visit the [Diamond official website](https://github.com/bbuchfink/diamond) to download the appropriate version for your system. This tool has been validated with diamond-2.1.11.

------

## 🔧 Installation

### Environment Setup

```bash
conda create -n FluRisk-env python=3.6
conda activate FluRisk-env

# Install required dependencies
conda install -c bioconda diamond=2.1.11
conda install -c conda-forge git-lfs

# Initialize Git LFS
git lfs install

# Clone the repository
git clone https://github.com/viralInformatics/FluRisk
cd FluRisk

pip install .
```

------

## 🚀 Usage

The `FluRisk` tool provides six subcommands for comprehensive influenza virus risk assessment.

> **💡 Tip**: For a complete risk assessment, execute modules sequentially: `anno` → `extract` → `predh` → `predv` → `predr` → `risk2tab`

### 1. 🏷️ Sequence Annotation (`anno`)

Performs sequence annotation through DIAMOND/BLAST alignment against influenza virus host databases, evaluates antigenic relationships between viruses and vaccines.

**Syntax:**

```bash
flupre anno -i <input> -o <output_directory> [options]
```

**Parameters:**

- `-i, --input` (required): Input FASTA file or directory containing FASTA files
- `-u, --updated_directory`: Directory for standardized FASTA files (default: `standardized_fasta/`)
- `-o, --output_directory`: Output directory for annotation results (default: `result/`)

**Example:**

```bash
cd test/
flupre anno -i test_files/ -o annotation_results/ -u standardized_sequences/
```

**Output Files:**

- ```
  standardized_sequences/
  ```

  : Standardized protein sequences with normalized sequence IDs

  - For nucleotide input: `.trans2protein.fasta.stdName.annotation.fa`
  - For protein input: `.stdName.annotation.fa`

- ```
  annotation_results/
  ```

  : Contains four annotation result files

  - `antigenInfo.blast.csv`: Top 5 most similar HA proteins and their antigenic/genetic relationships
  - `antigenInfo.csv`: Predicted antigenic and genetic relationships between virus strains and vaccine strains (for H1, H3, H5, H7, H9)
  - `for_anno.csv`: Sequence annotation results file

### 2. 🔍 Marker Extraction (`extract`)

Identifies amino acid molecular markers associated with five influenza virus phenotypic traits.

**Syntax:**

```bash
flupre extract -i <input> -a <annotation_path> -o <output_directory> [options]
```

**Parameters:**

- `-i, --input` (required): Directory containing standardized FASTA sequences (from `anno` step)
- `-a, --anno_path` (required): Directory containing annotation result files (from `anno` step)
- `-o, --output_directory`: Output directory for marker extraction results (default: current directory). Five folders will be created in the current directory: `adaptation/`, `virulence/`, `binding/`, `transmissibility/`, and `resistance/`. Each folder contains marker identification result files for the corresponding phenotype.
- `-p, --prefix`: Optional prefix for output filenames

**Example:**

```bash
flupre extract -i standardized_sequences/ -a annotation_results/ -p strain1
```

**Output Structure:**

```
adaptation/       # Mammalian Adaptability markers
binding/          # Human Receptor Binding Ability markers
transmissibility/ # Mammalian Transmissibility markers
virulence/        # Mammalian Virulence markers
resistance/       # Drug Resistance markers
```

Each directory contains CSV files with identified molecular markers for the corresponding phenotype (one file per strain).

### 3. 🎯 Host Prediction (`predh`)

Predicts influenza virus host classification (avian/human) based on mammalian adaptation-related molecular markers.

**Syntax:**

```bash
flupre predh -i <input> [options]
```

**Parameters:**

- `-i, --input` (required): Input CSV file or directory containing adaptation marker data. Specify the directory containing adaptability molecular markers generated in the previous step.
- `-t, --threshold`: Prediction probability threshold (default: 0.5)
- `-o, --output_directory`: Output directory (default: `host_predictions/`)
- `-p, --prefix`: Output filename prefix

**Example:**

```bash
flupre predh -i adaptation/ -t 0.5 -o host_predictions/ -p experiment1
```

**Output Files:**

- `xxx_prediction.csv`: Host prediction results for all input strains
- `combined_markers_matrix.csv`: Feature matrix used for model prediction

### 4. ⚡ Virulence Prediction (`predv`)

Predicts influenza virus virulence level (pathogenic/non-pathogenic) based on mammalian virulence-related molecular markers.

**Syntax:**

```bash
flupre predv -i <input> [options]
```

**Parameters:**

- `-i, --input` (required): Input CSV file or directory containing virulence marker data. Specify the directory containing virulence molecular markers generated in the previous step.
- `-t, --threshold`: Prediction probability threshold (default: 0.5)
- `-o, --output_directory`: Output directory (default: `virulence_predictions/`)
- `-p, --prefix`: Output filename prefix

**Example:**

```bash
flupre predv -i virulence/ -t 0.5 -o virulence_predictions/ -p experiment1
```

**Output Files:**

- `xxx_prediction.csv`: Virulence prediction results for all input strains
- `combined_markers_matrix.csv`: Feature matrix used for model prediction

### 5. 🔗 Binding Prediction (`predr`)

Predicts influenza virus receptor binding preference (α2,3/α2,6-linked sialic acid) based on human receptor binding-related molecular markers.

**Syntax:**

```bash
flupre predr -i <input> [options]
```

**Parameters:**

- `-i, --input` (required): Input CSV file or directory containing binding marker data. Specify the directory containing human receptor binding-related molecular markers generated in the previous step.
- `-t, --threshold`: Prediction probability threshold (default: 0.5)
- `-o, --output_directory`: Output directory (default: `binding_predictions/`)
- `-p, --prefix`: Output filename prefix

**Example:**

```bash
flupre predr -i binding/ -t 0.5 -o binding_predictions/ -p experiment1
```

**Output Files:**

- `xxx_prediction.csv`: Receptor binding preference prediction results for all input strains
- `combined_markers_matrix.csv`: Feature matrix used for model prediction

### 6. 🕸️ Risk Assessment (`risk2tab`)

Calculates comprehensive risk scores for five phenotypic traits based on identified molecular markers and prediction results.

**Syntax:**

```bash
risk2tab -i <input> -a <annotation_path> [options]
```

**Parameters:**

- `-i, --input` (required): Directory containing phenotype marker results
- `-a, --anno_path` (required): Directory containing annotation result files (default: `result/`)
- `-hp, --host_predictions`: Directory containing host prediction results (default: `host_predictions/`)
- `-vp, --virulence_predictions`: Directory containing virulence prediction results (default: `virulence_predictions/`)
- `-bp, --binding_predictions`: Directory containing binding prediction results (default: `binding_predictions/`)
- `-o, --output_directory`: Output directory for risk assessment results

**Example:**

```bash
# Basic risk assessment
risk2tab -i ./ -a annotation_results/ -o radar_file/

# Complete risk assessment with all predictions
risk2tab -i ./ -a annotation_results/ -hp host_predictions/ -vp virulence_predictions/ -bp binding_predictions/ -o radar_file/
```

**Output Files:**

- `[strain]_risk_values.csv`: Individual strain risk scores for four phenotypes and six drugs
- `summary_risk_values.csv`: Comprehensive risk assessment summary for all strains

**Risk Score Table Structure:**

| Strain  | Mammalian Transmissibility | Human Receptor Binding Ability | Mammalian Adaptability | Mammalian Virulence | zanamivir | oseltamivir | ...  |
| ------- | -------------------------- | ------------------------------ | ---------------------- | ------------------- | --------- | ----------- | ---- |
| sample1 | 0.78                       | 0.62                           | 0.91                   | 0.45                | 0.12      | 0.34        | ...  |

------

## 📊 Output Formats

### Annotation Output

- Standardized protein sequences with normalized sequence identifiers
- Antigenic relationship assessment files (CSV format)
- Reassortment risk evaluation results
- Comprehensive sequence annotation tables

### Marker Extraction Output

- Phenotype-specific molecular marker files organized by trait type
- Detailed marker annotation with protein positions and functional information
- Standardized CSV format for downstream analysis

### Prediction Output

- Host classification results with confidence scores
- Virulence level predictions with probability measures
- Receptor binding preference classifications
- Feature matrices used for machine learning predictions

### Risk Assessment Output

- Individual strain risk profiles across all phenotypic traits
- Comprehensive risk summary tables for comparative analysis
- Drug resistance risk scores for six antiviral compounds

------

## 🤖 Model Information

The prediction modules employ ensemble machine learning models trained on curated influenza sequence datasets. The tool integrates multiple phenotypic prediction models to provide comprehensive risk assessment capabilities for influenza surveillance and pandemic preparedness.

**Phenotype Descriptions:**

- **Mammalian Transmissibility**: Ability to transmit between mammalian hosts
- **Human Receptor Binding Ability**: Preference for human-type (α2,6-linked sialic acid) receptors
- **Mammalian Adaptability**: Adaptation to mammalian host environments
- **Mammalian Virulence**: Pathogenicity in mammalian hosts
- **Drug Resistance**: Resistance to antiviral medications

------

## 📦 Dependencies

All Python dependencies are specified in the `requirements.txt` file and are automatically installed during the package installation process. Key dependencies include:

- BioPython for sequence processing
- Pandas and NumPy for data manipulation
- Scikit-learn for machine learning components
- Diamond for sequence alignment

------

## 📝 Notes

> ⚠️ **Important**: It is essential to use the same versions of the software and models that were used during the initial training and validation for consistency and accuracy.

> 💡 **Note**: The Probability column in prediction results represents the probability of the assigned class, not the positive class probability.

For further assistance or to report issues, please visit the [GitHub Issues](https://github.com/viralInformatics/FluRisk/issues) section of the project repository.

------

<div align="center"></div>
