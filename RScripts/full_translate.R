# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

# Check if package is installed (install it if not)
install_if_needed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install packages 
install_if_needed("dplyr")
install_if_needed("stringr")

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("Biostrings")
}

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(Biostrings)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
fasta_file <- args[1]  
counts_fwd <- args[2]
counts_rev <- args[3]
fasta_output_dir <- args[4]  
best_start_pos_file <- args[5]

# Extract base name of input file
base_name <- tools::file_path_sans_ext(basename(fasta_file))

# Extract start position from csv
best_start_pos <- read.table(best_start_pos_file, header = FALSE, sep = ",")
frame <- as.integer(best_start_pos[1, 1])
cat("Translating from position:", frame)

# Read the DNA sequences and compute reverse complement
DNA_seqs <- readDNAStringSet(fasta_file)
rev_DNA_seqs <- Biostrings::reverseComplement(DNA_seqs)

# Function to save counts to csv 
generate_counts <- function(sequence, output_filepath){
  # Convert DNAStringSet to df
  sequence_df <- data.frame(
    target_name = names(sequence),
    Sequence = as.character(sequence)
  )
  
  # Generate count table 
  counts <- sequence_df %>%
    group_by(Sequence) %>%
    summarise(
      count = n(),
      target_name = head(target_name, 1),
      .groups = 'drop'
    )
  
  # Save counts as CSV
  write.csv(counts, file = output_filepath, row.names = FALSE)
}

# Function to translate sequences
translate_frames <- function(fwd_seqs, rev_seqs, frame, base_name, 
                             output_dir, counts_fwd, counts_rev) {
  options(warn = -1)

  # Subset and translate sequences
  fsubseqs <- subseq(fwd_seqs, start = frame)
  ftranslated <- Biostrings::translate(fsubseqs, if.fuzzy.codon = "X")
  rsubseqs <- subseq(rev_seqs, start = frame)
  rtranslated <- Biostrings::translate(rsubseqs, if.fuzzy.codon = "X")
  
  # Create output filenames
  foutput_file <- file.path(output_dir, paste0(base_name, "_fwd_AA.fasta"))
  routput_file <- file.path(output_dir, paste0(base_name, "_rev_AA.fasta"))
  
  # Write translated sequences to files
  writeXStringSet(ftranslated, filepath = foutput_file, format = "fasta")
  writeXStringSet(rtranslated, filepath = routput_file, format = "fasta")

  # Generate counts for forward and reverse sequences
  generate_counts(ftranslated, counts_fwd)
  generate_counts(rtranslated, counts_rev)
  
  # Restore warning settings to default
  options(warn = 0)
}

# Process forward and reverse strands
translate_frames(DNA_seqs, rev_DNA_seqs, frame, base_name, fasta_output_dir, counts_fwd, counts_rev)
message("Processed and saved FASTA files for: ", fasta_file)