library(Biostrings)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
fasta_file <- args[1]  
count_output_dir <- args[2]
fasta_output_dir <- args[3]  

# Extract base name of input file
base_name <- tools::file_path_sans_ext(basename(fasta_file))

# Extract unique file naming (first occurrence of bc10XX)
bc10XX <- str_extract(base_name, "bc10\\w+")

# Read the DNA sequences and compute reverse complement
DNA_seqs <- readDNAStringSet(fasta_file)
rev_DNA_seqs <- Biostrings::reverseComplement(DNA_seqs)

# Function to save counts to csv 
generate_counts <- function(sequence, output_dir, prefix){
  # Convert DNAStringSet to df
  sequence_df <- data.frame(
    target_name = names(sequence),
    Sequence = as.character(sequence)
  )
  
  # Generate count table 
  counts <- sequence_df %>%
    group_by(target_name) %>%
    summarise(
      count = n(),
      Sequence = first(Sequence),
      .groups = 'drop'
    )
  
  # Save counts as CSV
  output_file <- file.path(output_dir, paste0(prefix, "_counts.csv"))
  write.csv(counts, file = output_file, row.names = FALSE)
}

# Function to translate sequences
translate_frames <- function(fwd_seqs, rev_seqs, bc10XX, 
                             output_dir, count_output_dir) {
  for (frame in 1:3) {
    # Subset and translate sequences
    fsubseqs <- subseq(fwd_seqs, start = frame)
    ftranslated <- Biostrings::translate(fsubseqs, if.fuzzy.codon = "X")
    rsubseqs <- subseq(rev_seqs, start = frame)
    rtranslated <- Biostrings::translate(rsubseqs, if.fuzzy.codon = "X")
    
    # Create output filenames
    foutput_file <- file.path(output_dir, paste0(bc10XX, "_fwd_AA_frame", frame, ".fasta"))
    routput_file <- file.path(output_dir, paste0(bc10XX, "_rev_AA_frame", frame, ".fasta"))
    
    # Write translated sequences to files
    writeXStringSet(ftranslated, filepath = foutput_file, format = "fasta")
    writeXStringSet(rtranslated, filepath = routput_file, format = "fasta")
    message("Saved: ", foutput_file, " and ", routput_file)

    # Generate counts for forward and reverse sequences
    generate_counts(ftranslated, count_output_dir, paste0(bc10XX, "_fwd_frame", frame))
    generate_counts(rtranslated, count_output_dir, paste0(bc10XX, "_rev_frame", frame))
  }
}

# Process forward and reverse strands
translate_frames(DNA_seqs, rev_DNA_seqs, bc10XX, fasta_output_dir, count_output_dir)
message("Processed and saved 6 FASTA files for: ", fasta_file)