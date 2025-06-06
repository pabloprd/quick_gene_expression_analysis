#Summmary Stats and Differential Gene expression:
#By: Paul Parodi
#Last updated: 6/6/2025



#' Pairwise differential gene expression with Limma
#'
#' This function runs pairwise differential expression using the Limma package between the two groups you are comparing in a column
#'
#' @param normalized_counts  Dataframe of your normalized counts with sample names at the top and genes as the index.
#' @param metadata_table Metadata table associated with your dataframe
#' @param testing_column Column you are conducting the pairwise comparison in
#' @param column1 One of the categorical group names you are comparing 'Male'
#' @param column2 Another of the categorical group names you are comparing 'Female'
#'
#'
#' @return A data matrix containing the differential gene expression between two groups. Returns measurements logFC, Average expression, test statistic, P. Value, adjusted P. Value, and B statistic
#'
#' @examples
#' # Limma_dge_prep(normalized_counts_dataframe, metadata_table, 'Sex", 'M', 'F')
#'
#' @export
Limma_dge_prep <- function(normalized_counts, metadata_table, testing_column, column1, column2){
  #May have to edit this later to not remove the last three columns as I'll just insert the table with those columns taken out

  #Creating matrices and groups
  matrix_of_counts <- as.matrix(normalized_counts) # Make the matrix of just the counts. change variable name to be better.
  sample_names <- metadata_table[[testing_column]] # The sample names of the metadata column you are looking at
  group <- factor(sample_names) # Grouping the factors of the sample_names
  design <- model.matrix(~0 + group) # Creating a matrix of the new groups
  colnames(design) <- levels(group) # Make the different levels of the group the column names of the new thing.

  #First linear model
  fit <- lmFit(matrix_of_counts, design)

  #Contrast matrix
  contrast_name = paste0(column1, 'v', column2) # Make the title
  contrast_formula <- paste0(column1, "-" ,column2) # Putting in the title

  contrast.matrix <- limma::makeContrasts(contrasts = setNames(contrast_formula, contrast_name),
                                          levels = design) #Creating the contrast matrix that compares the groups

  #More refined linear model
  fit2 <- contrasts.fit(fit, contrast.matrix) #making the linear model based off this contrast matrix
  fit2 <- eBayes(fit2) #running Bayesian correction on the model
  Results <- topTable(fit2, adjust.method = 'fdr', number = Inf) #Outputting the files to results.

  return(Results)
}


#' Calculating summary stats and pairwise gene expression
#'
#' This function calculates the summary stats and whether the two groups you are comparing have statistically significant differences in gene expression
#'
#' @param normalized_counts  Dataframe of your normalized counts with sample names at the top and genes as the index.
#' @param metadata_table Metadata table associated with your dataframe
#' @param testing_column Column you are conducting the pairwise comparison in
#' @param column1 One of the categorical group names you are comparing ex: 'Male'
#' @param column2 Another of the categorical group names you are comparing ex: 'Female'
#' @param sample_char The character that the sample name starts with
#' @param sample_col_name The column name in the metadata that has the sample names. It automatically assigns it as "Sample Name"
#'
#' @return A data matrix containing the summary stats and differential gene expression between two groups. Returns measurements means, samples, differences of each sample,  logFC, Average expression, test statistic, P. Value, adjusted P. Value, and B statistic
#'
#' @examples
#' # sumstats_and_dge(normalized_counts_dataframe, metadata_table, 'Sex", 'M', 'F', 'GS', "Sample Name)
#'
#' @export
sumstats_and_dge <- function(normalized_counts, metadata_table, testing_column, column1, column2, sample_char = NULL, sample_col_name = "Sample Name")
{
  #Core DGE analysis (unchanged)
  dge_limma_table <- Limma_dge_prep(normalized_counts, metadata_table, testing_column, column1, column2)

  #Unified summary stats calculation
  get_group_stats <- function(group){
    tbl <- normalized_counts %>%
      gene_tables(metadata_table, sample_col_name, testing_column, group, sample_char)

    tbl_stats <- t(apply(tbl, 1, \(row) {
      numeric_row <- as.numeric(row)
      c(mean = mean(numeric_row, na.rm = TRUE),
        total_samples = sum(!is.na(numeric_row)),
        Max_difference = abs(max(numeric_row) - min(numeric_row)))
    }))

    stats_df <- data.frame(tbl_stats)
    colnames(stats_df) <- paste0(group, c('_mean', '_total_samples', '_Max_difference'))
    stats_df
  }
  #Process both groups
  ss1 <- get_group_stats(column1)
  ss2 <- get_group_stats(column2)[rownames(ss1),] #organizing the row names

  #Final Integration
  cbind(ss1, ss2, dge_limma_table[rownames(ss1), ])

}

