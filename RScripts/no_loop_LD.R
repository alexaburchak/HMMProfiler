# Load required packages
library(stringdist)
library(dplyr)
library(tidyr)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_csv <- args[1]
vh_col <- args[2]
vl_col <- args[3]
output_csv <- args[4]

# Read the input table
table <- read.csv(input_csv, stringsAsFactors = FALSE)

# Get environment variables
query_size <- as.numeric(Sys.getenv("QUERY_SIZE", unset = 100))  # Default: 100
method <- Sys.getenv("METHOD", unset = "random")  # Default: random
num_clust <- as.numeric(Sys.getenv("NUM_CLUST", unset = 5))  # Default: 5
max_dist <- as.numeric(Sys.getenv("MAX_DIST", unset = 5))  # Default: 5


# Test variables 
table <- read.csv("~/Desktop/pipeline_out/trim_tables/P18157_1012_S12_L001.csv", stringsAsFactors = FALSE)
vh_col <- "fwd_trimmed_seq"
vl_col <- "rev_trimmed_seq"
output_csv <- "~/Desktop/lv_test.csv"
query_size <- 100
method <- "random"
num_clust <- 5
max_dist <- 5

# Function to create query set by sampling HMMER results
generate_query <- function(table, sample_size, method) {
  set.seed(123)
  
  if (method == "top") {
    query_data <- table %>% arrange(desc(count)) %>% head(sample_size)
  } else if (method == "bottom") {
    query_data <- table %>% arrange(count) %>% head(sample_size)
  } else if (method == "random") {
    query_data <- table %>% sample_n(sample_size)
  } else {
    stop("Invalid method. Choose 'top', 'bottom', or 'random'.")
  }
  
  return(query_data)
}

# Calculate Levenshtein distances
calculate_levenshtein_all <- function(table, vh_col, vl_col, 
                                      query_size = 10, 
                                      method = "random",
                                      cluster = "heavy",
                                      num_clust = 5,
                                      max_dist = 5) {  # Add max_dist as an argument
  # Generate query data
  query_data <- generate_query(table, query_size, method)
  
  # Compute pairwise Levenshtein distances
  cat("Calculating Levenshtein distances for all queries...\n")
  vh_distances <- stringdist::stringdistmatrix(query_data[[vh_col]], table[[vh_col]], method = "lv")
  vl_distances <- stringdist::stringdistmatrix(query_data[[vl_col]], table[[vl_col]], method = "lv")
  
  # Add row and column names
  rownames(vh_distances) <- query_data$target_name  
  colnames(vh_distances) <- table$target_name
  
  rownames(vl_distances) <- query_data$target_name
  colnames(vl_distances) <- table$target_name
  
  # Clustering options
  if(cluster == "heavy"){
    # Filter rows and columns where all distances are <= max_dist
    filtered_indices <- apply(vh_distances, 2, function(column) all(column <= 20))
    query_data <- query_data[,filtered_indices]
    lev_distances <- vh_distances[filtered_indices, filtered_indices]
    
    if (nrow(query_data) == 0) {
      stop("No rows satisfy the max distance condition. Try increasing max_dist.")
    }
    
  }else if(cluster == "light"){
    # Filter rows and columns where all distances are <= max_dist
    filtered_indices <- apply(vl_distances, 1, function(row) all(row <= max_dist))
    query_data <- query_data[filtered_indices, ]
    lev_distances <- vl_distances[filtered_indices, filtered_indices]
    
    if (nrow(query_data) == 0) {
      stop("No rows satisfy the max distance condition. Try increasing max_dist.")
    }
    
  }else{
    cat("Clustering combined heavy and light distances...")
    # Cluster combined vh and vl 
    combined_distances <- sqrt(vh_distances^2 + vl_distances^2)
    
    # Filter rows and columns where all distances are <= max_dist
    filtered_indices <- apply(combined_distances, 1, function(row) all(row <= max_dist))
    query_data <- query_data[filtered_indices, ]
    lev_distances <- combined_distances[filtered_indices, filtered_indices]
    
    if (nrow(query_data) == 0) {
      stop("No rows satisfy the max distance condition. Try increasing max_dist.")
    }
  }
  
  # Hierarchical clustering
  hc <- hclust(dist(lev_distances))  
  clusters <- cutree(hc, k = num_clust) 
  query_data$cluster <- clusters
  
  # Plot dendrogram
  cat("Generating dendrogram plot...\n")
  png("dendrogram.png")  # Save as PNG file
  plot(hc, main = "Hierarchical Clustering of Filtered Query Sequences")
  dev.off()
  
  return(list(
    lev_distances = lev_distances,
    query_data = query_data,
    dendrogram = hc
  ))
}

# Run calculation
results <- calculate_levenshtein_all(
  table = table, 
  vh_col = vh_col, 
  vl_col = vl_col, 
  query_size = 100, 
  method = "random",
  cluster = "heavy",
  num_clust = 5,
  max_dist = 10
)

# Extract outputs
lev_distances <- results$lev_distances
query_data <- results$query_data

# Write outputs
#write.csv(lev_distances, "~/Desktop/pipeline_out/distance_matrix.csv", row.names = TRUE)
write.csv(query_data, "~/Desktop/pipeline_out/query_data.csv", row.names = FALSE)

cat("Outputs saved to pipeline_out:\n")
#cat("- Levenshtein distance matrix: distance_matrix.csv\n")
cat("- Query data with clusters: query_data.csv\n")
cat("- Dendrogram: dendrogram.png\n")
