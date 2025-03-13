# Pipeline for Identifying Functional Regions in Antibodies from NGS Data

## Overview
This pipeline processes next-generation sequencing (FASTQ) data to identify and quantify functional regions in antibodies. It first translates nucleotide sequences in all six reading frames using SeqKit, then runs HMMER to search for matches against a given model. Finally, it merges and trims the matched sequences, counts unique sequence combinations, and outputs the count results as a CSV file.

## Install System Dependencies
This pipeline requires several tools and dependencies to be installed before execution. You can use one of the following installation methods: 

### Option 1: Using Conda (Preferred)
```bash
# Create conda environment
conda create -n ngs_env python==3.11.7
conda acitvate ngs_env

# Install dependencies 
conda install -c bioconda seqkit -y
conda install -c conda-forge libcxx=16 -y
conda install -c bioconda hmmer -y
```

### Option 2: Using Homebrew (MacOS/Linux)
```bash
brew install seqkit hmmer
```

## How to Use

1. **Configure Necessary Parameters**
Make sure that the following directory paths are correctly set up in your configuration file (pipe_configs.json):

- input_pairs: 
  - fastq_path: Path to FASTQ file.
  - model_path: Path to the HMM profile file.
- counts_outpath: Path to write all unique sequence combinations found in matches output and their frequency.
- min_quality: Minimum quality score for FASTQ filtering.

2. **Run the Script**:
Assuming your config file is named pipe_configs.json, from the command line you can run: 
```bash
node main.js -c pipe_configs.json
```

3. **Outputs**
- The pipeline outputs one csv file that is saved to the assigned output path:
  1. `counts_outpath`: All unique sequence combinations found and their frequency.  
    - {model_name}_seq: Trimmed sequences for each model searched. Each model will have its own CSV column. 
    - count: Count of occurences of each combination of sequences in matches_outpath. 

## Workflow

### Step 1: Set Parameters
- Reads config file and extract parameters. 
- Defines path for final output file.

### Step 2: Quality Filtering and Translation (fullSeqKit)
- Filters sequences using `seqkit` with a minimum quality threshold of set in your config file.
- Translates DNA sequences into amino acid sequences in all reading frames.
- Outputs translated reads to FASTA format.

### Step 3: Run HMMER (runHMMSearch)
- Runs `hmmsearch` with provided HMM model.
- Saves HMMER output to temporary domain table file (.tbl). 

### Step 4: Extract Best HMMER Hits (extractBestHMMHits)
- Parses the .domtblout file to extract the highest-scoring hit per target.
- Generates a BED file containing alignment coordinates for the best matches. 

### Step 5: Trim Sequences Based on HMMER Hits (trimSeqs)
- Uses seqkit subseq to extract trimmed sequences based on coordinates from the BED file.
- Outputs trimmed sequences in FASTA format.

### Step 6: Map Trimmed Sequences to Targets (mapFastaSeqs)
- Reads the trimmed FASTA file and maps sequences to their corresponding target names and models.
- Organizes data in a structured map format.

### Step 7: Count Unique Sequence Combinations (countSeqs) and Write Results to CSV (writeCSV)
- Aggregates unique sequence combinations based on mapped sequences.
- Counts occurrences across all targets and models.
- Saves final results as a CSV file to `counts_outpath`.

## Notes
- The name of your model should contain some information about the type of structure being identified (CDR, VH/VL sequence, etc.) for successful sequence pairing in countSeqs(). This information will be extracted and used to name the sequence columns in the final counts_outpath. 
