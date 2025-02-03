suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(stringr)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
fasta_file <- args[1]         
fasta_output_dir <- args[2]   
numeric_prefix <- args[3]     # Numeric prefix for output files

# Read DNA sequences and compute reverse complement
DNA_seqs <- readDNAStringSet(fasta_file)
rev_DNA_seqs <- Biostrings::reverseComplement(DNA_seqs)

# Function to translate sequences in all reading frames
translate_frames <- function(fwd_seqs, rev_seqs, prefix, output_dir) {
  set.seed(123)  
  options(warn = -1)
  
  for (frame in 1:3) {
    # Subset sequences for the current frame
    fsubseqs <- subseq(fwd_seqs, start = frame)
    rsubseqs <- subseq(rev_seqs, start = frame)
    
    # Sample 1000 sequences 
    fsubseqs_sample <- fsubseqs[sample(length(fsubseqs), 1000)]
    rsubseqs_sample <- rsubseqs[sample(length(rsubseqs), 1000)]
    
    # Translate sequences
    ftranslated <- Biostrings::translate(fsubseqs_sample, if.fuzzy.codon = "X")
    rtranslated <- Biostrings::translate(rsubseqs_sample, if.fuzzy.codon = "X")
    
    # Create output filenames using the numeric prefix
    foutput_file <- file.path(output_dir, paste0(prefix, "_fwd", frame, ".fasta"))
    routput_file <- file.path(output_dir, paste0(prefix, "_rev", frame, ".fasta"))
    
    # Write translated sequences to files
    writeXStringSet(ftranslated, filepath = foutput_file, format = "fasta")
    writeXStringSet(rtranslated, filepath = routput_file, format = "fasta")
  }
  # Restore warning settings to default
  options(warn = 0)
}

# Process forward and reverse strands
translate_frames(DNA_seqs, rev_DNA_seqs, numeric_prefix, fasta_output_dir)