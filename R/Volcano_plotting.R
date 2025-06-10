#Volcano plotting 
#By: Paul Parodi
#Last updated: 6/6/2025

#' Marking regulation
#'
#' Makes a new column marking upregulated or downregulated depending on p.value and log fold change
#'
#' @param table dataframe of normalized counts
#' @param col1 Fold change
#' @param col2 adjusted p. value
#' @param col1sig significance you care about for adjusted p value
#' @param col2sig significance you care about for fold change
#'
#'
#' @return a new dataframe which has a column that says upregulated or downregulated
#'
#' @examples
#' regulation_level(dataframe, fold_change_col, adj_pval_col, 2.0, 0.05)
#'
#' @export
regulation_level <- function(table, col1, col2, col1sig, col2sig){
  table$regulation <- "NO"
  table$regulation[table[[col1]] > col1sig & table[[col2]] < col2sig] <- "UP"
  table$regulation[table[[col1]] < -col1sig & table[[col2]] < col2sig] <- "DOWN"
  return(table)
}






#' Volcano plotting gene expression
#'
#' Volcano plot that shows gene expression that is upregulated or downregulated between groups
#'
#' @param table dataframe of normalized counts
#' @param fc_col Fold change
#' @param adj_pval_col adjusted p. value
#' @param fc_cutoff fold change cutoff
#' @param pvalue p value
#' @param gene_amount amount of genes you want to be labeled
#' @param title Title of the plot you are making
#'
#' @return Volcano plot showing clearly which genes are upregulated or downregulated based on your parameters
#'
#' @examples
#' regulation_level(dataframe, fold_change_col, adj_pval_col, 2.0, 0.05, 30, 'African Ancestry Vs European Ancestry')
#'
#' @export
volcano_regulation_expression <- function(table, fc_col, adj_pval_col, fc_cutoff, pvalue, gene_amount, title){
  #Creating a regulation column for the table
  regulation_table <- regulation_level(table, fc_col, adj_pval_col, fc_cutoff, pvalue)

  #Turning the indices into a column called gene symbol
  regulation_table <- rownames_to_column(regulation_table, var = 'gene_symbol')


  #New line fixing factor order
  regulation_table$regulation <- factor(
    regulation_table$regulation,
    levels = c("DOWN", "NO", "UP")  # Explicit order
  )

  print(table(regulation_table$regulation))

  #Selecting the top amount
  regulation_table$delabel <- ifelse(regulation_table$gene_symbol %in% head(regulation_table[order(regulation_table[[adj_pval_col]]), "gene_symbol"], gene_amount), regulation_table$gene_symbol, NA)

  labeled_regulation_plot <- ggplot(data = regulation_table , aes(x = logFC, y = -log10(adj.P.Val), col = regulation, label = delabel)) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(pvalue), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_colour_manual(values = c("blue", "gray", "red"), # to set the colors of our variable
                       labels = c("Downregulated", "Not significant", "Upregulated"),   drop = FALSE)+ # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    geom_text_repel(max.overlaps = Inf)+
    ggtitle(title)

  return(labeled_regulation_plot)
}
