# Pipeline for Identifying and Quantifying VH-VL Pairs from Next-Generation Sequencing Data

## Overview
This pipeline processes next-generation sequencing (FASTQ) data to identify and quantify VH-VL antibody sequence pairs. It first translates nucleotide sequences in all six reading frames using SeqKit, then runs HMMER to search for matches against a given model. Finally, it merges and trims the matched sequences, counts unique VH-VL pairs, and outputs the results as CSV files.

## Install System Dependencies
This pipeline requires several tools and dependencies to be installed before execution. You can use one of the following installation methods: 

### Option 1: Using Conda (Preferred)
```bash
# Create conda environment
conda create -n pipeline_env
conda acitvate pipeline_env

# Install dependencies 
conda install -c bioconda jq seqkit  
conda install -c bioconda hmmer
```

### Option 2: Using Homebrew (MacOS/Linux)
```bash
brew install jq seqkit hmmer
```

## How to Use

1. **Configure Necessary Parameters**
Make sure that the following directory paths are correctly set up in your configuration file (pipe_configs.json):

- input_pairs: 
  - fastq_path: Path to FASTQ file.
  - model_path: Path to the HMM profile file.
- matches_outpath: Path to write all matches found by hmmsearch and their sequence information.
- counts_outpath: Path to write all unique VH/VL pairs found in matches output and their frequency.
- min_quality: Minimum quality score for FASTQ filtering.

2. **Run the Script**:
Assuming your config file is named pipe_configs.json, from the command line you can run: 
```bash
node main.js -c pipe_configs.json
```

3. **Outputs**
- The pipeline outputs two csv files that are saved to the assigned output paths:
  1. `matches_outpath`: All matches found by hmmsearch and their sequence information. 
    - target_name: Name of read (from FASTQ file).
    - score: Bit score calculated by hmmsearch.
    - e_value: E-value calculated by hmmsearch.
    - ali_from: Start position of match.
    - ali_to: End position of match.
    - sequence: Full sequence.
    - FASTQ_filename: Name of FASTQ file where match was found. 
    - model_name: HMM used to run hmmsearch. 
    - trimmed_seq: Sequence trimmed to match using ali_from and ali_to. 
    - seq_len: Length of trimmed sequence. 
    - start_pos: Start position for translation (can be 1, 2, 3, -1, -2, or -3). 


  2. `counts_outpath`: All unique VH/VL pairs found in ngs_hmm_matches and their frequency.  
    - vh_sequence: Trimmed heavy sequence.
    - vl_sequence: Trimmed light sequence. 
    - count: Count of occurences of each VH/VL pair in matches_outpath. 

## Workflow

### Step 1: Set Parameters
- Reads config file and extract parameters. 
- Defines paths for final output files.

### Step 2: Quality Filtering and Translation (fullSeqKit)
- Filters sequences using `seqkit` with a minimum quality threshold of set in your config file.
- Translates DNA sequences into amino acid sequences in all reading frames.
- Converts translated reads to FASTA format.

### Step 3: Run HMMER (runHMMSearch)
- Runs `hmmsearch` with provided HMM model.
- Saves HMMER output to temporary .tbl file. 

### Step 4: Parse HMMER Output (parseFullHMMOutput)
- Extracts relevant columns from hmmer output. 
- Keeps the highest scoring match per target. 

### Step 5: Merge HMMER Output with FASTA sequences (mergeData) and Write to CSV (writeCSV)
- Trims sequences based on ali_from and ali_to coordinates from hmmsearch. 
- Generates final output file `matches_outpath`.

### Step 6: Count Unique VH/VL Pairs (countSeqs) and Write to CSV (writeCSV)
- Reads `matches_outpath` line-by-line to identify VH/VL pairs by target_name, FASTQ_filename and model_name.
- Stores pairs in an object and counts frequency of each pair. 
- Generates final output file `counts_outpath`

## NOTES
- The name of your model should contain either 'VH' or 'VL' for successful sequence pairing in countSeqs(). 
