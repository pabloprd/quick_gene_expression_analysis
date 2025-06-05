#Outputting the linear models for differential gene expression
#' @export
Limma_mdge_prep <- function(normalized_counts, metadata_table, testing_column){
  exprs <- as.matrix(normalized_counts)
  sample_names <- metadata_table[[testing_column]]
  group <- factor(sample_names)
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  groups <- colnames(design)
  contrast_pairs <- combn(groups, 2, simplify = FALSE)
  contrast_strings <- sapply(contrast_pairs, function(pair) paste(pair, collapse = "-"))
  contrast.matrix <- do.call(makeContrasts, c(lapply(contrast_strings,
                                                     function(c) parse(text = c)[[1]]), list(levels = design)))
  fit <- lmFit(exprs, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  return(fit2)
}

#Outputting the results table for DGE
#' @export
Limma_mdge_results <- function(normalized_counts, metadata_table, testing_column){
  linear_model <- Limma_mdge_prep(normalized_counts, metadata_table, testing_column)
  Results <- topTable(linear_model, adjust.method = 'fdr', number = Inf) #Outputting the files to results.
}

#Outputting results that can be plotted into a Venn diagram, up regulated, down regulated etc.
#' @export
Limma_mdge_vd_results <- function(normalized_counts, metadata_table, testing_column){
  linear_model <- Limma_mdge_prep(normalized_counts, metadata_table, testing_column)
  result3 <- decideTests(linear_model, adjust = "BH", p.value = 0.05, lfc = log2(2))
}
