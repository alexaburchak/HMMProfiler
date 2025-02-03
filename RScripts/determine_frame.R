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
install_if_needed("readr")
install_if_needed("tidyr")

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# Get command line arguments 
args <- commandArgs(trailingOnly = TRUE)
hmmer_directory <- args[1]

# Get all tbl files in the directory
tbl_files <- list.files(hmmer_directory, pattern = "\\.tbl$", full.names = TRUE)

# Extract numerical prefix and frame from file names
file_info <- data.frame(
  file_path = tbl_files,
  stringsAsFactors = FALSE
) %>%
  mutate(
    file_base = basename(file_path),  # Extract only the file name
    group = as.integer(sub("_.*", "", file_base)),  # Extract the group
    frame = sub(".*_(fwd\\d|rev\\d)\\.tbl", "\\1", file_base)  # Extract the frame
  )


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

# Read and combine all files
all_frames <- file_info %>%
  rowwise() %>%
  mutate(data = list(read_domtblout(file_path))) %>%
  unnest(cols = c(data)) %>%
  ungroup()

# Filter for low e-values and count occurrences per group and frame
low_evalue_counts <- all_frames %>%
  filter(e_value <= 1e-5) %>%
  count(group, frame, name = "low_evalue_count") %>%
  arrange(group, desc(low_evalue_count))

# Determine the best frame for each group
best_frames <- low_evalue_counts %>%
  group_by(group) %>%
  slice_max(order_by = low_evalue_count, n = 1) %>%
  ungroup() %>%
  separate(frame, into = c("direction", "start_pos"), sep = "(?<=\\D)(?=\\d)") %>%
  mutate(start_pos = as.numeric(start_pos))

# Determine the best start position for translation 
best_start_pos <- best_frames %>%
  count(start_pos) %>%
  arrange(desc(n)) %>%
  slice(1) %>%
  pull(start_pos)

# Save best frame to csv file 
write.table(data.frame(best_start_pos), 
            file = "~/Desktop/best_start_pos.csv", 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = ",")
