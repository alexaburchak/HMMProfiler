# Description: Uses msa to perform multiple sequence alignment on queries obtained from
# a csv file. Writes aligned sequences to fasta files to generate hmms 

# Rscript msa.R <input.csv>

# Load necessary packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

BiocManager::install("Biostrings")
library(Biostrings)

BiocManager::install("msa")
library(msa)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the user provided the correct number of arguments (1 or 2)
if (length(args) < 1 || length(args) > 2) {
  stop("Please provide one mandatory argument: <path/to/input.csv>")
}

# Assign arguments to variables
csv_path <- args[1]

# Prompt user to provide chain and sequence type 
cat("Please provide the chain type (Enter 'heavy', 'light', 'kappa' or 'lambda'):")
chain <- tolower(trimws(readLines("stdin", n = 1)))

cat("Do you want to translate sequences before alignment? (Enter 'Y' or 'N'): ")
type <- tolower(trimws(readLines("stdin", n = 1)))

# Validate the chain input
if (!(chain %in% c("heavy", "light", "kappa", "lambda"))) {
  stop("Invalid chain type. Please enter 'heavy', 'light', 'kappa' or 'lambda'.")
}

# Validate the type input
if (!(type %in% c("dna", "aa"))) {
  stop("Invalid response. Please enter 'Y' or 'N'.")
}

# Assign columns based on chain 
if (chain == "heavy") {
  v_call_col <- "v_call_heavy"
  d_call_col <- "d_call_heavy"
  j_call_col <- "j_call_heavy"
} else {
  v_call_col <- "v_call_light"
  d_call_col <- "d_call_light"
  j_call_col <- "j_call_light"
}

# Print column names 
cat("V Call Column:", v_call_col, "\n")
cat("D Call Column:", d_call_col, "\n")
cat("J Call Column:", j_call_col, "\n")


# Function to process sequences from CSV files
align_sequences <- function(csv_path, chain, type, 
                              v_call_col, d_call_col, j_call_col) {
  
  # Read the CSV file
  sequence_csv <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  if (chain == "lambda") {
    # Filter lambda sequences for IGL
    filtered_csv <- subset(sequence_csv,
                           grepl("IGL", sequence_csv[[v_call_col]]) |
                             grepl("IGL", sequence_csv[[d_call_col]]) |
                             grepl("IGL", sequence_csv[[j_call_col]]))
    
    # Convert lambda sequences to DNAStringSet object
    DNA_sequence <- DNAStringSet(filtered_csv$sequence_light)
    names(DNA_sequence) <- filtered_csv$sequence_id_light
    
  } else if (chain == "kappa") {
    # Filter kappa sequences for IGK
    filtered_csv <- subset(sequence_csv,
                           grepl("IGK", sequence_csv[[v_call_col]]) |
                             grepl("IGK", sequence_csv[[d_call_col]]) |
                             grepl("IGK", sequence_csv[[j_call_col]]))
    
    # Convert kappa sequences to DNAStringSet object
    DNA_sequence <- DNAStringSet(filtered_csv$sequence_light)
    names(DNA_sequence) <- filtered_csv$sequence_id_light
    
  } else if (chain == "heavy") {
    # Filter heavy sequences for IGHV3
    filtered_csv <- subset(sequence_csv,
                           grepl("IGHV3", sequence_csv[[v_call_col]]) |
                             grepl("IGHV3", sequence_csv[[d_call_col]]) |
                             grepl("IGHV3", sequence_csv[[j_call_col]]))
    
    # Convert heavy sequences to DNAStringSet object
    DNA_sequence <- DNAStringSet(filtered_csv$sequence_heavy)
    names(DNA_sequence) <- filtered_csv$sequence_id_heavy
    
  } else {
    # Convert all light sequences to DNAStringSet object
    DNA_sequence <- DNAStringSet(sequence_csv$sequence_light)
    names(DNA_sequence) <- sequence_csv$sequence_id_light
  }
  
  if (type == "y") {
    # Translate DNA seq to AA 
    translated_seq <- Biostrings::translate(DNA_sequence, if.fuzzy.codon = "X")

    # Perform MSA on translated amino acid sequences
    AA_alignment <- msa(translated_seq)
    
    # Convert alignment to AAStringSet
    AA_aligned_set <- as(AA_alignment, "AAStringSet")
    
    # Save aligned sequences to FASTA file
    writeXStringSet(AA_aligned_set, filepath = file.path("~/Desktop/", paste0("aligned_", chain, "_AA.fasta")))
    
  } else if (type == "n") {
    # Perform MSA on DNA sequences
    DNA_alignment <- msa(DNA_sequence)
    
    # Convert alignment to DNAStringSet
    DNA_aligned_set <- as(DNA_alignment, "DNAStringSet")
    
    # Save aligned sequences to FASTA file
    writeXStringSet(DNA_aligned_set, filepath = file.path("~/Desktop/", paste0("aligned_", chain, "_DNA.fasta")))
  }
  print("Output saved to desktop!")
}

# Example usages
# align_sequences(
#   csv_path = "~/Downloads/SRR12483521_1_Paired_All.csv",
#   chain = "heavy",
#   type = "n",
#   v_call_col = "v_call_heavy",
#   d_call_col = "d_call_heavy",
#   j_call_col = "j_call_heavy"
# )
# 
# align_sequences(
#   csv_path = "~/Downloads/full_SRR12875356_1_Paired_All.csv",
#   chain = "kappa",
#   type = "n",
#   v_call_col = "v_call_light",
#   d_call_col = "d_call_light",
#   j_call_col = "j_call_light"
# )


