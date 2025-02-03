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
hmmer_tbl <- args[1]  
count_file <- args[2]
csv_output <- args[3] 
fasta_output <- args[4]
# need to add optional second hmmer file

# read_domtblout function (adapted from rhmmer package)
read_domtblout <- function(file){
  # Column types for domtblout format
  col_types <- readr::cols(
    target_name         = readr::col_character(),
    domain_accession    = readr::col_character(),
    domain_len          = readr::col_integer(),
    query_name          = readr::col_character(),
    query_accession     = readr::col_character(),
    qlen                = readr::col_integer(),
    e_value     = readr::col_double(),
    score      = readr::col_double(),
    bias       = readr::col_double(),
    domain_N            = readr::col_integer(),
    domain_of           = readr::col_integer(),
    domain_cevalue      = readr::col_double(),
    domain_ievalue      = readr::col_double(),
    domain_score        = readr::col_double(),
    domain_bias         = readr::col_double(),
    hmm_from            = readr::col_integer(),
    hmm_to              = readr::col_integer(),
    ali_from            = readr::col_integer(),
    ali_to              = readr::col_integer(),
    env_from            = readr::col_integer(),
    env_to              = readr::col_integer(),
    acc                 = readr::col_double()
  )
  
  N <- length(col_types$cols)
  
  # the line delimiter should always be just "\n", even on Windows
  lines <- readr::read_lines(file, lazy=FALSE, progress=FALSE)
  
  table <- sub(
    pattern = sprintf("(%s).*", paste0(rep('\\S+', N), collapse=" +")),
    replacement = '\\1',
    x=lines,
    perl = TRUE
  ) %>%
    gsub(pattern="  *", replacement="\t") %>%
    paste0(collapse="\n") %>%
    readr::read_tsv(
      col_names=names(col_types$cols),
      comment='#',
      na='-',
      col_types = col_types,
      lazy=FALSE,
      progress=FALSE
    )
  
  descriptions <- lines[!grepl("^#", lines, perl=TRUE)] %>%
    sub(
      pattern = sprintf("%s *(.*)", paste0(rep('\\S+', N), collapse=" +")),
      replacement = '\\1',
      perl = TRUE
    )
  
  table$description <- descriptions[!grepl(" *#", descriptions, perl=TRUE)]
  
  table
}

# Function to trim sequences 
process_hmmer <- function(hmmer_tbl,count_file, sequence) {
  # Read hmmerouts 
  hmmer_tbl <- read_domtblout(hmmer_tbl)
  counts <- read.table(count_file, header = TRUE, sep = ",")
  
  # Trim sequences based on HMMER alignments
  hmmer_tbl <- hmmer_tbl %>%
    group_by(target_name) %>%
    slice_max(domain_score, with_ties = FALSE) %>% # keep highest scoring domain 
    ungroup() %>%
    dplyr::select(target_name, domain_len, domain_score, e_value, ali_from, ali_to) %>%
    inner_join(counts, by = "target_name") %>% # add sequence info
    mutate(
      trimmed_seq = substr(Sequence, ali_from, ali_to), # trim AA sequence 
      trimmed_seq_len = nchar(trimmed_seq) 
    ) 
    
  return(hmmer_tbl)
}

# Function to trim sequences (heavy + light)
process_duo_hmmer <- function(heavy_hmmer_tbl, light_hmmer_tbl, sequence) {
  # Read hmmerouts 
  heavy_hmmer_tbl <- read_domtblout(heavy_hmmer_tbl)
  light_hmmer_tbl <- read_domtblout(light_hmmer_tbl)
  sequence <- readAAStringSet(sequence)
  
  # Convert AAStringSet to df
  sequence_df <- data.frame(
    target_name = names(sequence),
    Sequence = as.character(sequence)
  )
  
  # Trim sequences based on HMMER alignments
  heavy_hmmer_tbl <- heavy_hmmer_tbl %>%
    group_by(target_name) %>%
    slice_max(domain_score, with_ties = FALSE) %>% # keep highest scoring domain 
    ungroup() %>%
    dplyr::select(target_name, domain_len, domain_score, ali_from, ali_to) %>%
    inner_join(sequence_df, by = "target_name") %>% # add sequence info
    mutate(
      trimmed_seq = substr(Sequence, ali_from, ali_to), # trim AA sequence 
      trimmed_seq_len = nchar(trimmed_seq) 
    ) %>%
    rename_with(~ paste0("H_", .), -target_name)
  
  light_hmmer_tbl <- light_hmmer_tbl %>%
    group_by(target_name) %>%
    slice_max(domain_score, with_ties = FALSE) %>% # keep highest scoring domain 
    ungroup() %>%
    dplyr::select(target_name, domain_len, domain_score, ali_from, ali_to) %>%
    inner_join(sequence_df, by = "target_name") %>% # add sequence info
    mutate(
      trimmed_seq = substr(Sequence, ali_from, ali_to), # trim AA sequence 
      trimmed_seq_len = nchar(trimmed_seq) 
    ) %>%
    rename_with(~ paste0("L_", .), -target_name)
  
  # Merge heavy and light data
  combined_tbl <- heavy_hmmer_tbl %>%
    inner_join(light_hmmer_tbl, by = "target_name")
  
  return(combined_tbl)
}

# Write table of results to csv 
table <- process_hmmer(hmmer_tbl, count_file, sequence)
write.csv(table, file = csv_output)

# Write trimmed sequences to fasta file 
trimmed_seqs <- AAStringSet(table$trimmed_seq)
names(trimmed_seqs) <- table$target_name
writeXStringSet(trimmed_seqs, filepath = fasta_output)