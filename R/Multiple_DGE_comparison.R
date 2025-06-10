#FUNCTIONS THAT DO MULTIPLE PAIRWISE COMPARISONS ON GENE EXPRESSION
#Created by: Paul Parodi
#Last updated: 5/30/2025


#' Multiple differential gene expression with Limma
#'
#' This function runs multiple differential gene expression using the Limma package between the multiple groups you are comparing in a column
#'
#' @param normalized_counts  Dataframe of your normalized counts with sample names at the top and genes as the index.
#' @param metadata_table Metadata table associated with your dataframe
#' @param testing_column Column you are conducting the pairwise comparison in
#'
#' @return a matrix of the comparisons
#'
#' @examples
#' # Limma_mdge_prep(normalized_counts_dataframe, metadata_table, 'Income')
#'
#' @export
Limma_mdge_prep <- function(normalized_counts, metadata_table, testing_column){
  exprs <- as.matrix(normalized_counts) #Creates matrix of normalized counts
  sample_names <- metadata_table[[testing_column]] #Gets sample names of metadata
  group <- factor(sample_names) #converting a vector of values into a categorical variable with a fixed set of possile values called levels
  design <- model.matrix(~0 + group) #creates a design matrix for the linear modeling where each column corresponds to one level of the group factor
  colnames(design) <- levels(group) #Making the column names of design be named after each level of the group factor
  groups <- colnames(design) #stores the names of the columns in the design matrix (which are group names like control and treatment) into a variable called groups
  contrast_pairs <- combn(groups, 2, simplify = FALSE) #Finds all possible pairs of group names
  contrast_strings <- sapply(contrast_pairs, function(pair) paste(pair, collapse = "-")) #for each pair it joins the two group with "_"
  contrast.matrix <- do.call(makeContrasts, c(lapply(contrast_strings,
                                                     function(c) parse(text = c)[[1]]), list(levels = design))) #Outputs a matrix where each tells limma which group comparisons are wanted to make in the analysis
  fit <- lmFit(exprs, design) #Making linear model
  fit2 <- contrasts.fit(fit, contrast.matrix) #
  fit2 <- eBayes(fit2) #A new fitted model that contains results for my specified contrasts, rather than just the original group coefficients.
  return(fit2)
}


#' Multiple differential gene expression with Limma outputting the dataframe
#'
#' This function makes a dataframe out of the multiple differential gene expression matrix.
#'
#' @param normalized_counts  Dataframe of your normalized counts with sample names at the top and genes as the index.
#' @param metadata_table Metadata table associated with your dataframe
#' @param testing_column Column you are conducting the pairwise comparison in
#'
#' @return A dataframe containing the differential gene expression between all groups. Returns measurements Average expression, test statistic, P. Value, adjusted P. Value, and F statistic
#'
#' @examples
#' # Limma_mdge_results(normalized_counts_dataframe, metadata_table, 'Income')
#'
#' @export
Limma_mdge_results <- function(normalized_counts, metadata_table, testing_column){
  linear_model <- Limma_mdge_prep(normalized_counts, metadata_table, testing_column) #here is the linear model
  Results <- topTable(linear_model, adjust.method = 'fdr', number = Inf) #Outputting the files to results.
  return(Results)
}


#' Multiple differential gene expression venn diagram
#'
#' This function makes a venn diagram out of the multiple differential gene expression matrix.
#'
#' @param normalized_counts  Dataframe of your normalized counts with sample names at the top and genes as the index.
#' @param metadata_table Metadata table associated with your dataframe
#' @param testing_column Column you are conducting the pairwise comparison in
#'
#' @return a venn diagram of the matrices results
#'
#' @examples
#' # Limma_mdge_vd_results(normalized_counts_dataframe, metadata_table, 'Income')
#'
#' @export
Limma_mdge_vd_results <- function(normalized_counts, metadata_table, testing_column){
  linear_model <- Limma_mdge_prep(normalized_counts, metadata_table, testing_column) #Linear model
  result3 <- decideTests(linear_model, adjust = "BH", p.value = 0.05, lfc = log2(2)) #Venn diagram model
  return(result3)
}
