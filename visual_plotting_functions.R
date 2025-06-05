#Plotting functions


#Boxplot plotting function

#' @export
Single_gene_boxplot <- function(counts, metadata_table, column, gene, column_name){
  merged_gene_dataframe <-counts %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = paste0(column_name), values_to = "Expression") %>%
    inner_join(metadata_table, by = paste0(column_name)) %>%
    filter(Gene == gene)

  merged_gene_dataframe[[column]] <- as.factor(merged_gene_dataframe[[column]])

  ggplot(merged_gene_dataframe, aes(x = .data[[column]], y = Expression, fill = .data[[column]])) +
    geom_boxplot(color = "black", size = 0.7) +
    scale_fill_viridis_d(option = "viridis", alpha = 0.8) +
    geom_jitter(width = 0.2, color = "gray20", alpha = 0.7) +
    labs(title = paste("Expression of", gene), x = column , y = "Expression") +

    #coord_cartesian(xlim = c(0, 5), ylim = c(-10, 6))+ #You can take this out potentially and just have it vary between groups

    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      legend.position = "right"
    )
}

#Heat mapping function
#' @export
heatmapping_function <-function(comparison_table, column_name_TO, top_gene_nums = 0, bottom_gene_nums = 0){
  if (top_gene_nums != 0) {
    top <- top_selection(comparison_table, column_name_TO, top_gene_nums)
    if(bottom_gene_nums!= 0){
      bottom <- bottom_selection(comparison_table, column_name_TO, bottom_gene_nums)
      combo <- rbind(top, bottom)
    }
    else{
      combo <- top
    }
  }
  else{
    bottom <- bottom_selection(comparison_table, column_name_TO, bottom_gene_nums)
    combo <- bottom
  }
  combo <- combo[ , !(names(combo) %in% c("P.Value", "adj.P.Val"))]
  combo <- arrange(combo, column_name_TO)
  data <- as.matrix(combo)
  heat <- pheatmap((data), legend = TRUE, cluster_rows = FALSE)
  return(heat)
}


#' @export
heatmapping_data_export <-function(comparison_table, column_name_TO, top_gene_nums = 0, bottom_gene_nums = 0){
  if (top_gene_nums != 0) {
    top <- top_selection(comparison_table, column_name_TO, top_gene_nums)
    if(bottom_gene_nums!= 0){
      bottom <- bottom_selection(comparison_table, column_name_TO, bottom_gene_nums)
      combo <- rbind(top, bottom)
    }
    else{
      combo <- top
    }
  }
  else{
    bottom <- bottom_selection(comparison_table, column_name_TO, bottom_gene_nums)
    combo <- bottom
  }
  combo <- combo[ , !(names(combo) %in% c("P.Value", "adj.P.Val"))]
  #combo <- combo[ , "logFC", drop = FALSE]
  combo <- arrange(combo, column_name_TO)
  data <- as.matrix(combo)
  return(data)
}





#What the group doing ontology wise plotting
#' @export
group_ontology_plot <- function(gene_list, level){
  ggo <- groupGO(gene     = gene_list,
                 OrgDb    = org.Hs.eg.db,
                 ont      = "BP",
                 level    = level,
                 readable = TRUE)
  return(barplot(ggo))
}

#The enrichment of the genes and stuff
#' @export
enrichment_plot <- function(gene_list, pval, qval){
  ego <- enrichGO(gene         = gene_list,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",       # BP: Biological Process
                  pAdjustMethod= "BH",
                  pvalueCutoff = pval,
                  qvalueCutoff = qval,
                  readable     = TRUE)

  return(dotplot(ego))
}
