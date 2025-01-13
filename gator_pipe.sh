#!/bin/bash

# Helper Functions
create_directories() {
    local base_output_dir=$1
    mkdir -p "$base_output_dir/filtered_reads/DNA" \
             "$base_output_dir/filtered_reads/AA" \
             "$base_output_dir/counts" \
             "$base_output_dir/hmmerout" \
             "$base_output_dir/trimmed_reads" \
             "$base_output_dir/trim_tables" \
             "$base_output_dir/lv_tables"
}

unzip_fastq_files() {
    local fastq_dir=$1
    echo "Unzipping FASTQ files..."
    for fastq_file in "$fastq_dir"/*.fastq.gz; do
        if [ -f "$fastq_file" ]; then
            echo "Unzipping $fastq_file..."
            pigz -d "$fastq_file"
        fi
    done
}


# Pipeline for AA analysis 
aa_pipeline() {
    local fastq_dir=$1
    local output_dir=~/Desktop/AA_pipeline_out
    local filter_AA_reads="$output_dir/filtered_reads/AA"
    create_directories "$output_dir"
    unzip_fastq_files "$fastq_dir"

    # Translate in all reading frames
    for file in "$fastq_dir"/*.fastq; do
        if [ -f "$file" ]; then
            # Define filenames
            local base_name=$(basename "${file}" .fastq)
            local filtered_fastq="filtered_${base_name}.fastq"
            local DNA_reads="$output_dir/filtered_reads/DNA/${base_name}.fasta"

            # Quality filtering
            echo "Filtering $file..."
            seqkit seq --min-qual 20 "$file" -o "$filtered_fastq"

            # Convert FASTQ to FASTA
            echo "Converting $filtered_fastq to FASTA..."
            seqkit fq2fa "$filtered_fastq" -o "$DNA_reads"

            # Translate DNA to AA sequences
            echo "Translating $DNA_reads..."
            Rscript ../thesisRCode/alligator_translate.R "$DNA_reads" "$output_dir/counts" "$filter_AA_reads"

            # Clean intermediate files
            pigz -6 "$file" "$DNA_reads"
            rm "$filtered_fastq" 

        fi
    done

    # HMMER processing
    echo "Running HMMER searches..."
    for file in "$filter_AA_reads"/*.fasta; do
        local base_name=$(basename "$file" .fasta)
        local heavy_hmmerfile="$output_dir/hmmerout/${base_name}_Hhmmer.tbl"
        local kappa_hmmerfile="$output_dir/hmmerout/${base_name}_Khmmer.tbl"

        # Run HMMER 
        echo "Running hmmsearch for $file..." 
        hmmsearch --domtblout "$heavy_hmmerfile" ~/Desktop/currentHMM/models/AA/VH-prot-mod.hmm "$file"
        hmmsearch --domtblout "$kappa_hmmerfile" ~/Desktop/currentHMM/models/AA/kappa_AA.hmm "$file"

        # Trim reads 
        echo "Trimming $file..." 
        Rscript ../thesisRCode/alligator_trim.R "$heavy_hmmerfile" "$kappa_hmmerfile" "$file" 

        # Zip intermediate files
        pigz -6 "$file" "$heavy_hmmerfile" "$kappa_hmmerfile" 
    done
}


# Pipeline for DNA Analysis 
dna_pipeline() {
    local fastq_dir=$1
    local output_dir=~/Desktop/DNA_pipeline_out
    create_directories "$output_dir"
    unzip_fastq_files "$fastq_dir"

    # Filter sequences and generate count tables 
    for file in "$fastq_dir"/*.fastq; do
        if [ -f "$file" ]; then
            # Define filenames
            local base_name=$(basename "${file}" .fastq)
            local filtered_fastq="filtered_${base_name}.fastq"
            local DNA_reads="$output_dir/filtered_reads/DNA/${base_name}.fasta"

            # Quality filtering
            echo "Filtering $file..."
            seqkit seq --min-qual 20 "$file" -o "$filtered_fastq"

            # Convert FASTQ to FASTA
            echo "Converting $filtered_fastq to FASTA..."
            seqkit fq2fa "$filtered_fastq" -o "$DNA_reads"

            # Clean intermediate files
            pigz -6 "$file" 
            rm "$filtered_fastq" 

        fi
    done

    # HMMER processing
    echo "Running HMMER searches..."
    for file in "$output_dir/filtered_reads/DNA"/*.fasta; do
        local base_name=$(basename "$file" .fasta)
        local heavy_hmmerfile="$output_dir/hmmerout/${base_name}_Hhmmer.tbl"
        local kappa_hmmerfile="$output_dir/hmmerout/${base_name}_Khmmer.tbl"

        # Run HMMER 
        echo "Running hmmsearch for $file..." 
        nhmmer --tblout "$heavy_hmmerfile" ~/Desktop/currentHMM/models/DNA/ NEED HMM FILE "$file"
        nhmmer --tblout "$kappa_hmmerfile" ~/Desktop/currentHMM/models/DNA/ NEED HMM FILE "$file"

        # Trim reads 
        echo "Trimming $file..." 
        Rscript ../thesisRCode/alligator_trim.R "$heavy_hmmerfile" "$kappa_hmmerfile" "$file" 

        # Zip intermediate files
        pigz -6 "$file" "$heavy_hmmerfile" "$kappa_hmmerfile" 
    done
}

# Main Script Logic
read -p "Enter the pathway to the directory containing FASTQ files: " fastq_dir
read -p "Which pipeline would you like to run? (AA/DNA): " pipeline_choice

if [[ $pipeline_choice == "AA" ]]; then
    aa_pipeline "$fastq_dir"
elif [[ $pipeline_choice == "DNA" ]]; then
    dna_pipeline "$fastq_dir"
else
    echo "Invalid choice. Please enter AA or DNA."
    exit 1
fi

# Clustering 
read -p "Enter query size (default 100): " QUERY_SIZE
QUERY_SIZE=${QUERY_SIZE:-100}  # Default to 100 if no input is provided

read -p "Enter sampling method ('top', 'bottom', or 'random'; default 'random'): " METHOD
METHOD=${METHOD:-random}  # Default to 'random' if no input is provided

read -p "Enter number of desired clusters (default 5): " NUM_CLUST
NUM_CLUST=${NUM_CLUST:-5}  # Default to 5 if no input is provided

read -p "Enter the name of the VH column: " vh_col
read -p "Enter the name of the VL column: " vl_col

# Export the environment variables
export QUERY_SIZE
export METHOD
export NUM_CLUST

# Calculate Levenshtein distances + basic clustering of queries 
tables="$output_dir/trim_tables"
lv_tables="$output_dir/lv_tables"
for table in "$tables"/*.csv; do 
    if [ -f "$table" ]; then
        echo "Calculating Levenshtein distances for $table..."
        Rscript ~/Desktop/thesisRCode/no_loop_LD.R "$table" "$vh_col" "$vl_col" "$lv_tables"
    fi
done

echo "Pipeline completed!"