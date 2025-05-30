# HMMProfiler: A Duo of Pipelines for Identifying Functional Regions and Matching Sequences from NGS Data

## Overview
This repository contains two main pipelines for processing next-generation sequencing (FASTQ) data to identify and quantify functional regions of interest. The first pipeline (counts_pipeline.js) processes raw sequencing data, translates nucleotide sequences, and counts unique sequence combinations. The second pipeline (matches_pipeline.js) takes the output from the first pipeline and matches query sequences against identified functional regions to find the closest matches. These pipelines provide an end-to-end solution for processing NGS data, identifying functional regions, and matching them against query sequences. While tailored for antibody repertoire analysis, the core logic is adaptable to other structured sequence domains provided appropriate profile HMMs. 

This repository includes the current Node.js-based implementation of the pipelines. A legacy version written in Bash/R is included in `bash_pipeline/` for reference only.

## Getting Started
1. **Clone the Repository:**
```bash
git clone https://github.com/alexaburchak/HMMProfiler.git
cd HMMProfiler
```

2. **Install Node.js Dependencies:**
Make sure you have Node.js installed (can be downloaded [here](https://nodejs.org/en/download)), then run:
```bash
npm install
```
This will install all required packages listed in `package.json`.

3. **Linter and Config Files:**
- This project uses a linter configured via `biome.json` to ensure consistent code style.
- Configuration files (`count_parameters.json`, `matches_parameters.json`) are provided as templates under `test_configs`—please edit them to suit your data and environment.

## Install System Dependencies
Both pipelines require several tools and dependencies to be installed before execution. You can use one of the following installation methods: 

### Option 1: Using Conda (Preferred)
```bash
# Create conda environment
conda create -n ngs_env 
conda activate ngs_env

# Install dependencies 
conda install -c bioconda seqkit=2.9.0 -y
conda install -c conda-forge libcxx=16 -y
conda install -c bioconda hmmer=3.4 -y
```

### Option 2: Using Homebrew (MacOS)
```bash
brew install seqkit hmmer
```

## Building a Custom Hidden Markov Model with HMMER
Some example HMM profiles are provided under `models/AA`, but you can also generate your own HMM using [HMMER](http://eddylab.org/software/hmmer/Userguide.pdf):

1. **Prepare your aligned FASTA file(`.sto`, `.fa`, `.fasta`):**
- I used the [MSA](https://www.bioconductor.org/packages/release/bioc/html/msa.html) package in R (V1.40.0), but you can also use a multiple sequence alignment tool like [Clustal Omega](http://www.clustal.org/omega/). 

2. **Build the HMM profile:**
```bash
hmmbuild your_model.hmm your_alignment.fa
```
The resulting `.hmm` file can be referenced in your configuration file under `model_path`.


## Pipeline 1: Identifying and Quantifying Functional Regions `counts_pipeline.js`

### How to use 

1. **Prepare Configuration File:** `count_parameters.json`
  - `input_pairs`: Array of objects containing: 
    - `fastq_path`: Path to FASTQ file.
    - `model_path`: Path to the HMM profile.
  - `counts_outpath`: Path to write all unique sequence combinations found in matches output and their frequency.
  - `min_quality`: Minimum quality score for FASTQ filtering.
  - `hmm_coverage`: Minimum required model coverage for hmmsearch hits.

2. **Run the Script**
You can test the counts pipeline by running:
```bash
cd HMMProfiler

# You will need to change "USERNAME" to your own username in counts_parameters.json
node counts_pipeline.js -c test_configs/count_parameters.json
```

3. **Output:** `counts_outpath`
  - `{model_name}_seq`: Trimmed sequences for each model searched. Each model will have its own CSV column. 
  - `Count`: Count of occurrences of each combination of sequences. 
  - `Total_Count`: Count of occurrences of all possible sequence combinations (sum of `Count` column). 
  - `Frequency`: Frequency of each sequence combination relative to all detected combinations.

### Workflow

**1. Load Configuration** - Reads the config file and defines output paths.

**2. Quality Filtering and Translation** - Filters sequences, translates them into amino acid sequences, and saves results in FASTA format.

**3. Run HMMER** - Searches sequences against the provided HMM model using `hmmsearch`.

**4. Extract Best HMMER Hits** - Generates a BED file containing alignment coordinates for the highest-scoring match per sequence. 

**5. Trim Sequences Based on HMMER Hits** - Extracts trimmed sequences based on alignment coordinates from the BED file.

**6. Map Trimmed Sequences to Targets** - Maps trimmed sequences to their corresponding target names and models.

**7. Count Unique Sequence Combinations** - Aggregates and counts unique sequence occurrences.

**8. Save Counts to CSV** - Outputs a count and frequency for each unique sequence combination.

## Pipeline 2: Matching Query Sequences `matches_pipeline.js`

### How to use

1. **Prepare Configuration File:** `matches_parameters.json`
  - `queryEntries`: Array of query objects containing:
    - `name`: Query name.
    - `sequences`: Array of amino acid strings. 
  - `libraries`: Array of objects containing: 
    - `name`: Name of library.
    - `model_paths`: Array of paths to HMM profiles.
    - `counts_path`: Path to the output CSV from `counts_pipeline.js`
  - `output_path`: Path to write CSV output of all detected matches. 
  - `max_LD`: Maximum Levenshtein distance for matching sequences. 

2. **Run the Script**
You can test the matching pipeline by running:
```bash
cd HMMProfiler

# You will need to change "USERNAME" to your own username in matches_parameters.json
node matches_pipeline.js -c test_configs/matches_parameters.json
```

3. **Output:** `output_path`
  - `Query`: Query name from input. 
  - `{modelName}_Match`: Sequence from `counts_path` matched to the query sequence.
  - `{modelName}_LD`: Distance between the query sequence and closest match.
  - `Library`: Name of the counts library being matched. 
  - Additional columns from `counts_pipeline.js` output. 

### Workflow

**1. Load Configuration** - Reads input parameters from match_parameters.json.

**2. Validate Headers** - Confirm that `counts_path` column names match the provided model names.

**3. Run HMMER on Query Sequences** - Searches for query sequences against the specified HMM model.

**4. Extract Best Hits** - Identifies the best-matching functional region.

**5. Trim Sequences** - Extracts relevant portions of query sequences.

**6. Find Closest Matches** - Compares query sequences to identified functional regions using Levenshtein distance.

**7. Save Matches to CSV** - Outputs the closest matches for each query sequence.

## Notes
- Ensure that model names contain information about the type of structure (e.g., CDR, VH/VL) for accurate sequence pairing. This information will be extracted and used to name the sequence columns in `counts_outpath`. 

- Temporary files generated during processing are automatically cleaned up.