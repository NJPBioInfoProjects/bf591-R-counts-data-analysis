#Imports
library(tidyverse)
library(DESeq2)

#' Load a tsv located at specific location `filename` into a tibble
#'
#'
#' @param filename (str): the path to a specific file (ie 'file/path/to/file.tsv')
#'
#' @return tibble: a (g x 1+m) tibble with a 'gene' column followed by
#' sample names as column names.
#'
#' @note Column 'gene' should be first and the only column to contain strings.
#' Data in sample_name columns CANNOT be strings
#'
#' @example `verse_counts <- read_data('verse_counts.tsv')`

read_data <- function(filename){
  data <- read_delim(filename)
  
  return(data)
}


#' Filter out genes with zero variance
#'
#'
#' @param verse_counts tibble: a (g x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns with sample names as column names.
#'
#' @return tibble: a (n x 1+m) tibble with a 'gene' column followed by m columns
#' of raw counts with genes that have zero variance across samples removed
#'
#' @note (g >= n)
#'
#' @example `filtered_counts <- filter_zero_var_genes(verse_counts)`

filter_zero_var_genes <- function(verse_counts) {
  filtered_cts <- verse_counts %>%
    filter(!apply(select(., -gene), 1, function(x) length(unique(x)) == 1))
  
  return(filtered_cts)    
}


#' Extract time point information from sample name
#'
#'
#' @param str string: sample name from count data.
#'
#' @return string: string character representing sample time point
#'
#' @example `timepoint_from_sample("vAd_1")`
#' output:`"Ad"`

# timepoint_from_sample which extracts the ages of the subjects (P0, P4, P7, and Ad) from the sample names
timepoint_from_sample <- function(x) {
  timepoint <- substr(x, 2, 3)
  
  
  return(timepoint)
}


#' Grab sample replicate number from sample name
#'
#'
#' @param str  string: sample name from count data.
#'
#' @return string: string character represent sample replicate number
#'
#' @example `sample_replicate("vAd_1")`
#' output: `"1"`

# sample_replicate which extracts the sample replicate number (“1” or “2”)
sample_replicate <- function(x) {
  replicate <- substr(x, 5, 6)
  
  return(replicate)
}


#' Generate sample-level metadata from sample names.
#'
#' Will include columns named "sample", "timepoint", and "replicate" that store
#' sample names, sample time points, and sample replicate, respectively.
#'
#'
#' @param sample_names vector: character vector of length (_S_) consisting of sample
#' names from count data.
#'
#' @return tibble: a (_S_ x 3) tibble with column names "sample",
#' "timepoint", and "replicate". "sample"holds sample_names; "timepoint"
#' stores sample time points; and "replicate" stores sample replicate
#'
#' @note _S_ < m
#'
#' @example `meta <- meta_info_from_labels(colnames(count_data)[colnames(count_data)!='gene'])`

meta_info_from_labels <- function(sample_names) {
  
  meta_info <- tibble(
    sample = sample_names,
    timepoint = timepoint_from_sample(sample_names),
    replicate = sample_replicate(sample_names)
  )  
  
  return(meta_info)
}


#' Calculate total read counts for each sample in a count data.
#'
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @return tibble of read totals from each sample. A tibble can be `(1 x _S_)` 
#' with sample names as columns names OR `(_S_ x 2)` with columns ("sample", "value")
#'
#' @examples `get_library_size(count_data)`

get_library_size <- function(count_data) {
  size <- count_data %>%
    select(-gene) %>%
    summarize(across(everything(), sum)) %>%
    pivot_longer(cols = everything(), names_to = 'sample', values_to = 'value')
  
  return(size)
}


#' Normalize raw count data to counts per million WITH pseudocounts using the
#' following formula:
#'     count / (sample_library_size/10^6)
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m columns of cpm normalized read counts
#'
#' @examples
#' `normalize_by_cpm(count_data)`

normalize_by_cpm <- function(count_data) {
  library_sizes <- colSums(count_data[-1])  # Exclude gene column for sums
  
  # Apply CPM normalization formula to each column
  df_normalized <- count_data %>%
    mutate(across(-gene, ~ . / library_sizes[as.character(cur_column())] * 1e6))
  
  return(df_normalized)
}

#' Normalize raw count matrix using DESeq2
#'
#' @param count_data tibble: a (n x 1+m) tibble with a 'gene' column followed
#' by m raw counts columns of read counts

#' @param meta_data tibble: sample-level information tibble corresponding to the
#' count matrix columns
#'
#' @return tibble: DESeq2 normalized count matrix
#' @export
#'
#' @examples
#' `deseq_normalize(count_data, meta_data)`

deseq_normalize <- function(count_data, meta_data) {
  # Prepare the count matrix (remove 'gene' column and set row names)
  count_matrix <- as.matrix(count_data %>% select(-gene))
  rownames(count_matrix) <- count_data$gene
  
  # Ensure meta_data aligns with the count matrix columns
  meta_data <- meta_data %>%
    filter(sample %in% colnames(count_matrix)) %>%
    arrange(match(sample, colnames(count_matrix)))
  
  # Convert meta_data to a DataFrame for DESeq2 (to avoid tibble row name warning)
  meta_data_df <- as.data.frame(meta_data)
  rownames(meta_data_df) <- meta_data_df$sample  # Set sample names as row names
  
  # Create DESeqDataSet with the trivial design formula
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = meta_data_df,
    design = ~ 1
  )
  
  # Run DESeq2 normalization (this does not fit the model, just normalizes)
  dds <- estimateSizeFactors(dds)
  
  # Extract the normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Convert to a tibble
  normalized_counts_tibble <- as_tibble(normalized_counts, rownames = "gene")
  
  return(normalized_counts_tibble)
}


#' Perform and plot PCA using processed data.
#'
#' PCA is performed over genes, and samples should be colored by time point.
#' Both `y` and `x` axis should have percent of explained variance included.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param meta tibble: sample-level meta information (_S_ x 3)
#' @param title string: title for plot
#'
#' @return ggplot: scatter plot showing each sample in the first two PCs.
#'
#' @examples
#' `plot_pca(data, meta, "Raw Count PCA")`

plot_pca <- function(data, meta_data, title = "PCA Plot") {
  data <- as_tibble(data)
  
  # Convert to matrix
  count_matrix <- data %>%
    select(-gene) %>%   
    t() %>%             
    as.data.frame()      
  
  # Ensure the row names of the count matrix are the sample names
  rownames(count_matrix) <- colnames(data)[-1]
  
  # Perform PCA on the transposed matrix)
  pca_result <- prcomp(count_matrix, scale. = TRUE)
  
  # Calculate the percentage of variance explained by each PC
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
  
  # Extract the first two PCs and convert to a data frame
  pca_data <- as.data.frame(pca_result$x[, 1:2])
  colnames(pca_data) <- c("PC1", "PC2")
  
  # Add sample names to the PCA data and join with metadata
  pca_data <- pca_data %>%
    mutate(sample = rownames(pca_data)) %>%
    inner_join(meta_data, by = "sample")
  
  # Create the scatter plot
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = timepoint)) +  # Use 'condition' from metadata
    geom_point() +
    labs(
      title = title,
      x = paste0("PC1: ", round(variance_explained[1], 1), "% variance"),
      y = paste0("PC2: ", round(variance_explained[2], 1), "% variance")
    )
  
  return(p)
}


#' Plot gene count distributions for each sample using boxplots.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale the `y` axis to log10 values.
#' Default is FALSE, and y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: boxplot show gene count distributions for each sample
#'
#' @example `plot_sample_distributions(data, scale_y_axis=TRUE, title='Raw Count Distributions')`

plot_sample_distributions <- function(data, scale_y_axis=FALSE, title="") {
    # Transform the data to a long format for ggplot2
    long_data <- data %>%
      pivot_longer(
        cols = -gene,  # Exclude 'gene' column from pivoting
        names_to = "sample", 
        values_to = "count"
      )
    
    # Create the base ggplot object
    p <- ggplot(long_data, aes(x = sample, y = count, color = sample)) +
      geom_boxplot() +
      #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = title,
        x = "sample",
        y = "counts"
      )
    
    # Apply log10 scaling to the y-axis if requested
    if (scale_y_axis) {
      p <- p + scale_y_log10() +
        labs(title = title, y = "counts")
    }
  
  # Example usage: Plot with log10 scale on the y-axis
}


#' Plot relationship between mean read counts and variability over all genes.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale to y-axis to log10 values. Default
#' is false, and the y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: A scatter plot where the x-axis is the rank of gene ordered by mean
#' count over all samples, and the y-axis is the observed variance of the
#' given gene. Each dot should have their transparency increased. The scatter
#' plot should also be accompanied by a line representing the average mean and
#' variance values.
#'
#' @example `plot_variance_vs_mean(data, scale_y_axis=TRUE, title='variance vs mean (raw counts)')`

plot_variance_vs_mean <- function(data, scale_y_axis = FALSE, title = "") {
  #calculate mean and variance
  var <- data %>% select(-gene) %>% apply(1, var)
  mean <- data %>% select(-gene) %>% apply(1, mean)
  
  var_mean <- tibble(
    variance = var,
    mean = mean
  )
  
  rownames(var_mean) <- filtered$gene
  
  p <- var_mean %>% ggplot() +
    geom_point(aes(x=rank(mean), y=log10(var))) +
    geom_smooth(aes(x = rank(mean), y = log10(var))) +
    ggtitle(title)
  
  return(p)
}

