args <- commandArgs(trailingOnly = TRUE)
heavy_hmmer_tbl <- args[1]  
light_hmmer_tbl <- args[2] 
sequence <- args[3]  
output <- args[4]

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
process_hmmer <- function(heavy_hmmer_tbl, light_hmmer_tbl, sequence) {
  # Read hmmerouts 
  heavy_hmmer_tbl <- read_domtblout(heavy_hmmer_tbl)
  light_hmmer_tbl <- read_domtblout(light_hmmer_tbl)
  
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
      H_trimmed_seq = substr(Sequence, ali_from, ali_to), # trim AA sequence 
      H_trimmed_seq_len = nchar(trimmed_seq) 
    ) 
  
  light_hmmer_tbl <- light_hmmer_tbl %>%
    group_by(target_name) %>%
    slice_max(domain_score, with_ties = FALSE) %>% # keep highest scoring domain 
    ungroup() %>%
    dplyr::select(target_name, domain_len, domain_score, ali_from, ali_to) %>%
    inner_join(sequence_df, by = "target_name") %>% # add sequence info
    mutate(
      L_trimmed_seq = substr(Sequence, ali_from, ali_to), # trim AA sequence 
      L_trimmed_seq_len = nchar(trimmed_seq) 
    ) 
  
  # Merge heavy and light data
  combined_tbl <- heavy_hmmer_tbl %>%
    inner_join(light_hmmer_tbl, by = "target_name")
  
  return(combined_tbl)
}

table <- process_hmmer(heavy_hmmer_tbl, light_hmmer_tbl, sequence)
write.csv(table, file = output)
