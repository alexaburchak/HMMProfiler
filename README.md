# DNA and AA Analysis Pipeline

## Overview
This pipeline processes DNA and amino acid (AA) sequences from FASTQ files. It performs quality filtering, translation, HMMER searches, and sequence trimming. The pipeline supports separate workflows for DNA and AA data.

## Install System Dependencies
This pipeline requires several tools and dependencies to be installed before execution. You can use one of the following installation methods: 

### Option 1: Using APT (Debian/Ubuntu)
```bash
sudo apt update && sudo apt install -y \
    jq \
    pigz \
    seqkit \
    hmmer \
    nhmmer \
    r-base
```

### Option 2: Using Conda (Preferred)
```bash
# Create conda environment
conda create -n pipeline_env
conda acitvate pipeline_env

# Install dependencies 
conda install -c bioconda jq pigz seqkit  
conda install -c conda-forge r-base
conda install -c bioconda hmmer
```

### Option 3: Using Homebrew (MacOS/Linux)
```bash
brew install jq pigz seqkit hmmer
brew install --cask r
```

## How to Use

1. **Configure the Environment**
Make sure that the following directory paths are correctly set up in your configuration (configs.json):

- fastq_dir: Directory containing FASTQ files.
- output_base_dir: Directory where output files will be stored.
- pipeline_mode: Either aa or dna.
- hmm_file: Path to the HMM profile file.
- min_quality: Minimum quality score for filtering.
- best_start_position: Set to auto or a specific frame.
- gzip_compression_level: Compression level for output files.

2. **Run the Script**:
   ```bash
   ./pipe_with_config.sh

3. **Outputs**
- The results are saved in:
  - `~/Desktop/AA_pipeline_out` (for the AA pipeline)
  - `~/Desktop/DNA_pipeline_out` (for the DNA pipeline)

## Workflow

### Step 1: Prepare Input Files
- Unzips `.fastq.gz` files in the specified input directory.

### Step 2: Determine translation frame (optional)
- If best_start_position is set to 'auto', a subset of sequences from each file will be processed to determine the optimal translation frame. 

### Step 3: Quality Filtering
- Filters sequences using `seqkit` with a minimum quality threshold of `20`.
- Converts filtered FASTQ files to FASTA format.

### Step 4: Translation (AA Pipeline)
- Translates DNA sequences into amino acid sequences in all reading frames using `alligator_translate.R`.

### Step 5: HMMER Searches
- Runs `hmmsearch` (AA pipeline) or `nhmmer` (DNA pipeline) with provided HMM models.
- Saves HMMER output in the `hmmerout` directory.

### Step 6: Trimming
- Trims sequences based on HMMER results using `alligator_trim.R`.

## Notes
- Ensure the necessary R scripts are accessible at the defined paths.


