Working parameters: 
- fastq_dir: Directory containing input FASTQ files.
- output_base_dir: Base directory for output files.
- pipeline_mode: Set to "aa" for Amino Acid pipeline or "dna" for DNA pipeline.
- gzip_compression_level: Compression level for pigz (1-9)
- hmm_file: Default HMM file
- min_quality: Minimum quality score for filtering
- best_start_position: "auto" (determine dynamically) or a specific frame (e.g., 1, 2, or 3)

Clustering steps for when this works: 
- query_size: Default query size
- sampling_method: "top", "bottom", "random"
- num_clusters: Number of desired clusters
- vh_column_name: Column name for VH in CSV files
- vl_column_name: Column name for VL in CSV files