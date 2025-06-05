#Makes a new column marking upregulated or downregulated depending on certain values
#' @export
regulation_level <- function(table, col1, col2, col1sig, col2sig){
  table$regulation <- "NO"
  table$regulation[table[[col1]] > col1sig & table[[col2]] < col2sig] <- "UP"
  table$regulation[table[[col1]] < -col1sig & table[[col2]] < col2sig] <- "DOWN"
  return(table)
}

#Plotting the gene expression
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

  #Selecting the top 20 most
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
