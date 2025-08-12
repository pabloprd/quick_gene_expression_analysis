# GENE ONTOLOGY ANNOTATIONS
#By: Paul Parodi
#Last updated: 6/5/2025



#' Cluster profiler Annotations
#'
#' Makes an annotation dataframe with every annotation associated with genes in the panel.
#'
#' @param normalized_counts data frame of normalized counts
#' @param universe_symbols list of gene eids that you are comparing against (leave this blank)
#' @param pvalcutoff p value cutoff you want to implement (set to 0.05)
#' @param qvalcutoff q value cutoff you want to implement (set tp 0.01)
#' @param OrgDB genome wide annotation database (automatically set as org.Hs.eg.db which is thhe homo sapiens database)
#' @param ont Ontology. Can either be: BP (Biological processes), MF (Molecular Function), CC (Cellular Components), ALL
#'
#'
#' @return a data matrix of annotated genes using the clusterprofiler(GO) database
#'
#' @examples
#' # GO_annotations(normalized_counts_df)
#'
#' @export
GO_annotations <- function(normalized_counts, universe_symbols = NULL,
                           pvalcutoff = 0.05, qvalcutoff = 0.01,
                           OrgDb = org.Hs.eg.db, ont = "BP") {
  # Convert gene symbols to ENTREZ IDs
  gene_symbols <- rownames(normalized_counts)
  gene_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
  universe_ids <- if (!is.null(universe_symbols)) {
    bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
  } else NULL

  ego <- enrichGO(
    gene          = gene_ids,
    universe      = universe_ids,
    OrgDb         = OrgDb,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = pvalcutoff,
    qvalueCutoff  = qvalcutoff,
    readable      = TRUE
  )
  as.data.frame(ego)
}



#' Returns the description of function of the gene annotation
#'
#'
#' returns the description of a function of the gene you are looking for, grabs the one with the highest pvalue related to the GO database.
#'
#' @param gene gene of interest that you are looking for
#' @param annotations annotation data frame
#' @param comparison_table table of comparisons you are looking at.
#'
#' @return returns first annotation associated with the gene
#'
#' @examples
#' # GO_annotations(normalized_counts_df)
#'
#' @export
GO_pvalue_finder <- function(gene, annotations, comparison_table){
  sorted_annotations <- annotations[order(annotations$p.adjust), ]
  for(i in 1:nrow(sorted_annotations)){
    row <- sorted_annotations[i, ]
    if (gene %in% unlist(strsplit(row$geneID, "/"))){
      return(row$Description)
    }
  }
  return(NA)
}



#' Assigning annotations to comparison table
#'
#' Returns the comparison table you input with the top annotation for each gene using the GO database
#'
#' @param comparison_table gene of interest that you are looking for
#' @param normalized_counts annotation data frame
#' @param annotations table of comparisons you are looking at.
#'
#' @return returns comparison table with added annotations
#'
#' @examples
#' # GO_annotations(Sex_comparison_table, normalized_counts, annotations)
#'
#' @export
GO_definitions_to_table <- function(comparison_table, normalized_counts, annotations){
  #annotations <- GO_annotations(comparison_table, normalized_counts)
  #list_of_genes <- normalized_counts$Name #this is the issue
  list_of_genes <- rownames(normalized_counts)
  comparison_table$Description <- NA
  for(gene in list_of_genes){
    description <- GO_pvalue_finder(gene, annotations, comparison_table)
    row_num <- which(rownames(comparison_table) == gene)
    comparison_table$Description[row_num] = description

  }
  return(comparison_table)
}


#' Annotation selection for certain genes
#'
#' Returns the annotation selection for certain genes. Enter in the gene, the annotation data frame you created and the number of annotations you want returned to you
#' associated with that gene. Be aware that there won't always be enough annotations.
#'
#' @param gene Gene you are interested in.
#' @param annotation_df The annotation data frame you are using
#' @param num  Number of annotations you want returned (may be less total annotations than requested)
#'
#' @return returns the annotation selection for certain genes
#'
#' @examples
#' # GO_anno_selection('CCL4', annotation_df, 5)
#'
#' @export
Go_anno_selection <- function(gene, annotation_df, num){
  sorted_annotations <- annotation_df[order(annotation_df$pvalue), ]
  table <- sorted_annotations %>% filter(grepl(gene, geneID))
  selected_table <- table %>% slice_min(order_by = p.adjust, n = num)
  return(selected_table)
}




#' List of entrez ids
#'
#' Make a list of entrez ids associated with the genes you submitted
#'
#' @param normalized_counts dataframe of counts, can be any dataframe that has indices as genes
#' @param OrgDb Database you are using to get the entrez_ids from: Autofills to (org.Hs.eg.db)
#'
#' @return returns a list of entrez_ids
#'
#' @examples
#' entrez_id_list(counts_dataframe)
#'
#' @export
entrez_id_list <- function(normalized_counts, OrgDb = org.Hs.eg.db) {
  genes <- rownames(normalized_counts)
  #Getting the entrez_ids
  entrez_ids <- bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID',  OrgDb = org.Hs.eg.db)
  return(entrez_ids$ENTREZID)
}



