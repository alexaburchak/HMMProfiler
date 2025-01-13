# DNA and AA Analysis Pipeline

## Overview
This pipeline processes DNA and amino acid (AA) sequences from FASTQ files. It performs quality filtering, translation, HMMER searches, and clustering analysis. The pipeline supports separate workflows for DNA and AA data.

## Prerequisites
Ensure the following dependencies are installed:
- **Bash** (compatible shell environment)
- **seqkit** (for sequence manipulation)
- **pigz** (for parallel gzip compression)
- **hmmsearch** and **nhmmer** (from HMMER suite)
- **R** with necessary scripts:
  - `alligator_translate.R`
  - `alligator_trim.R`
  - `no_loop_LD.R`

Organize required files:
- Place HMM models in the appropriate directories:
  - `~/Desktop/currentHMM/models/AA/`
  - `~/Desktop/currentHMM/models/DNA/`

## Directory Structure
The pipeline creates the following directories within the specified output directory:
- `filtered_reads/DNA`: Filtered DNA sequences in FASTA format
- `filtered_reads/AA`: Translated AA sequences in FASTA format
- `counts`: Count files generated during AA translation
- `hmmerout`: HMMER output files
- `trimmed_reads`: Trimmed sequence files
- `trim_tables`: Tables generated after trimming
- `lv_tables`: Tables with Levenshtein distances and clustering information

## How to Use

1. **Run the Script**:
   ```bash
   ./pipeline.sh

2. **Input Parameters**
- Provide the path to the directory containing FASTQ files.
- Select the pipeline type (AA for amino acids or DNA for DNA sequences).
- For clustering:
  - Specify query size (default: 100).
  - Choose a sampling method (top, bottom, or random; default: random).
  - Set the number of clusters (default: 5).
  - Enter the column names for VH and VL.

3. **Outputs**
- The results are saved in:
  - `~/Desktop/AA_pipeline_out` (for the AA pipeline)
  - `~/Desktop/DNA_pipeline_out` (for the DNA pipeline)

## Workflow

### Step 1: Prepare Input Files
- Unzips `.fastq.gz` files in the specified input directory.

### Step 2: Quality Filtering
- Filters sequences using `seqkit` with a minimum quality threshold of `20`.
- Converts filtered FASTQ files to FASTA format.

### Step 3: Translation (AA Pipeline)
- Translates DNA sequences into amino acid sequences in all reading frames using `alligator_translate.R`.

### Step 4: HMMER Searches
- Runs `hmmsearch` (AA pipeline) or `nhmmer` (DNA pipeline) with provided HMM models.
- Saves HMMER output in the `hmmerout` directory.

### Step 5: Trimming
- Trims sequences based on HMMER results using `alligator_trim.R`.

### Step 6: Clustering
- Calculates Levenshtein distances and clusters sequences using `no_loop_LD.R`.

## Environment Variables

- **`QUERY_SIZE`**: Size of query sequences for clustering (default: `100`).
- **`METHOD`**: Sampling method for queries (`top`, `bottom`, or `random`; default: `random`).
- **`NUM_CLUST`**: Number of clusters (default: `5`).

## Notes

- Ensure the necessary R scripts are accessible at the defined paths.
- Replace placeholder HMM file references with appropriate file paths for the DNA pipeline.
- Verify that all tools are in the system PATH or use full paths to the executables.

