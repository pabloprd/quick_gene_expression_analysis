#PLOTTING FUNCTIONS
#By: Paul Parodi
#6/6/2025

#' Single gene boxplot
#'
#' This function boxplots a column of your choice and a gene of your choice and allows you to compare them.
#'
#' @param counts Dataframe of normalized counts
#' @param metadata_table Metadata table associated with your dataframe
#' @param column Metadata table column you are interested in age, sex, etc.
#' @param gene Gene you want to visualize the expression of ex: 'CD4'
#' @param sample_names column that the sample names are under in the metadata.
#'
#'
#' @return a nice boxplot comparing groups
#'
#' @examples
#' # Single_gene_boxplot(normalized_counts, metadata_tablee , 'Ancestry', 'CCL4', "Sample Name")

#' @export
Single_gene_boxplot <- function(counts, metadata_table, column, gene, sample_names = "Sample Name"){
  merged_gene_dataframe <-counts %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = paste0(sample_names), values_to = "Expression") %>%
    inner_join(metadata_table, by = paste0(sample_names)) %>%
    filter(Gene == gene)

  merged_gene_dataframe[[column]] <- as.factor(merged_gene_dataframe[[column]])

  ggplot(merged_gene_dataframe, aes(x = .data[[column]], y = Expression, fill = .data[[column]])) +
    geom_boxplot(color = "black", size = 0.7) +
    scale_fill_viridis_d(option = "viridis", alpha = 0.8) +
    geom_jitter(width = 0.2, color = "gray20", alpha = 0.7) +
    labs(title = paste("Normalized Expression of", gene), x = column , y = "Expression") +

    #coord_cartesian(xlim = c(0, 5), ylim = c(-10, 6))+ #You can take this out potentially and just have it vary between groups

    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      legend.position = "right"
    )
}


#' Plotting group ontology
#'
#' Function to produce ontology plots for a gene list that you specify
#'
#' @param gene_list gene list
#' @param level level of ontology you want to look at 1,2,3,4
#'
#'
#' @return a plot of gene counts per ontology group
#'
#' @examples
#' group_ontology_plot(list_of_genes, 3)


#' @export
group_ontology_plot <- function(gene_list, level = 1){
  ggo <- groupGO(gene     = gene_list,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = level,
                 readable = TRUE)
  return(barplot(ggo))
}



#The enrichment of the genes and stuff
#Plots gene ratios below a certain p-value that you specify and shows the function they are associated with
#Gene ratio: statistic commonly used in GO or pathway enrichment analyses to describe the
#Proportion of your input genes that are associated with a specific GO term or pathway
# Gene ratio = Number of input genes annotated to a term/total number of input genes

#' @export
enrichment_plot <- function(gene_list, pval = 0.05, qval = 0.01){
  ego <- enrichGO(gene         = gene_list,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",       # BP: Biological Process
                  pAdjustMethod= "BH",
                  pvalueCutoff = pval,
                  qvalueCutoff = qval,
                  readable     = TRUE)

  return(dotplot(ego))
}



