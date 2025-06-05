#         TESTING DOCUMENT
#SUMSTATS AND DGE


#' @export
sumstats_and_dge <- function(normalized_counts_table, metadata_table,
                             testing_column, column1, column2, sample_char, sample_col_name = "Sample Name"){
  dge_limma_table <- Limma_dge_prep(normalized_counts_table,  metadata_table,
                                    testing_column, column1, column2)
  gene_table1 <- gene_tables(normalized_counts_table, metadata_table, sample_col_name,testing_column, column1, sample_char)
  gene_table2 <- gene_tables(normalized_counts_table,metadata_table,sample_col_name, testing_column,  column2, sample_char)
  ss1 <- Sumstatsselect(gene_table1, column1)
  ss2 <- Sumstatsselect(gene_table2, column2)
  #Resorting the tables so they bind correctly.
  ss2_reordered <- ss2[rownames(ss1), ]
  final_table <- cbind(ss1, ss2)
  #Resorting again
  final_table <- final_table[rownames(dge_limma_table), ]
  final_table <- cbind(final_table, dge_limma_table)

}


#SUMSTATSSELECT
#Function that gets the sum stats and also pares it down to three sum stats.
#' @export
Sumstatsselect <- function(table, column_name){
  sum_stats_table <- t(apply(table, 1,Sum_stats))
  sum_stats_table <- data.frame(sum_stats_table)
  sum_stats_table <- sum_stats_table %>% dplyr::select('mean', 'total_samples', 'Max_difference')
  names(sum_stats_table) <- paste0(column_name, c('_mean', '_total_samples', '_Max_difference'))
  return(sum_stats_table)
}


#' @export
Sum_stats <- function(row){
  numeric_row <- as.numeric(row)
  c(
    mean = mean(numeric_row, na.rm = TRUE),
    sd = sd(numeric_row, na.rm = TRUE),
    median = median(numeric_row, na.rm = TRUE),
    Max = max(numeric_row),
    Min = min(numeric_row),
    total_samples = sum(!is.na(numeric_row)),
    Max_difference = abs(max(numeric_row) - min(numeric_row))
  )
}


#' @export
Limma_dge_prep <- function(normalized_counts, metadata_table, testing_column, column1, column2){
  #May have to edit this later to not remove the last three columns as I'll just insert the table with those columns taken out
  exprs <- as.matrix(normalized_counts) #Make the matrix of just the counts.
  sample_names <- metadata_table[[testing_column]] #The sample names of the metadata column you are looking at
  group <- factor(sample_names) #Grouping the factors of the sample_names
  design <- model.matrix(~0 + group) #creating a matrix of the new groups
  colnames(design) <- levels(group) #Make the different levels of the group the column names of the new thing.
  fit <- lmFit(exprs, design) #Create a linear model
  contrast_name = paste0(column1, 'v', column2) #Make the title
  contrast_formula <- paste0(column1, "-" ,column2) #putting in the title
  contrast.matrix <- limma::makeContrasts(contrasts = setNames(contrast_formula, contrast_name),
                                          levels = design) #Creating the contrast matrix that compares the groups
  fit2 <- contrasts.fit(fit, contrast.matrix) #making the linear model based off this contrast matrix
  fit2 <- eBayes(fit2) #running bayesian correction on the model
  Results <- topTable(fit2, adjust.method = 'fdr', number = Inf) #Outputting the files to results.
  return(Results)
}
