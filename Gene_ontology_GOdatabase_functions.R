#GO Annotation stuff


#Gene Ontology for just the genes you are looking at. You can change the pvalue and q value to much higher ones
#If you need it.
#' @export
GO_annotations_panel <- function(raw_counts, table, normalized_counts, pvalcutoff = 0.05, qvalcutoff = 0.01){
  list_of_genes <- raw_counts$Name
  genes <- rownames(normalized_counts)

  #Getting the entrez_ids
  all_entrez_ids <- bitr(list_of_genes, fromType = 'SYMBOL', toType = 'ENTREZID',  OrgDb = org.Hs.eg.db)
  entrez_ids <- bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID',  OrgDb = org.Hs.eg.db)

  #Universe vectors
  valid_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")  #Getting all the valid entrezid that are in the GO database
  universe_vec <- as.character(all_entrez_ids$ENTREZID)
  universe_vec <- universe_vec[universe_vec %in% valid_entrez]
  universe_vec <- unique(universe_vec)                     #Get the unique and valid entrez ids for the panel

  #Calling the GO database
  ego <- enrichGO(gene          = as.character(entrez_ids$ENTREZID),
                  universe      = universe_vec, # Background gene list (optional, defaults to all genes in the annotation database)
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP", # Options: "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component), "ALL"
                  pAdjustMethod = "BH", # Benjamini-Hochberg for multiple testing correction
                  pvalueCutoff  = pvalcutoff,
                  qvalueCutoff  = qvalcutoff,
                  readable = TRUE) # More stringent cutoff for adjusted p-value
  return(data.frame(ego))
}


#Gene Ontology for the whole genome. Don't have to change the pvalues but you can.
#' @export
GO_annotations_global <- function(raw_counts, table, normalized_counts, pvalcutoff = 0.05, qvalcutoff = 0.01){
  list_of_genes <- raw_counts$Name
  genes <- rownames(normalized_counts)
  entrez_ids <- bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID',  OrgDb = org.Hs.eg.db)

  #Calling the GO database
  ego <- enrichGO(gene          = as.character(entrez_ids$ENTREZID),
                  universe      = , # Background gene list (optional, defaults to all genes in the annotation database)
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP", # Options: "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component), "ALL"
                  pAdjustMethod = "BH", # Benjamini-Hochberg for multiple testing correction
                  pvalueCutoff  = pvalcutoff,
                  qvalueCutoff  = qvalcutoff,
                  readable = TRUE) # More stringent cutoff for adjusted p-value
  return(data.frame(ego))
}

#GO_PVALUE_FINDER: RETURNS THE DESCRIPTION OF FUNCTION OF THE GENE ANNOTATION
#WITH THE LOWEST ADJUSTED P-VALUE (HIGHEST LIKELIHOOD)
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


#GO_DEFINITIONS_TO_TABLE: TAKES IN TABLES, ANNOTATIONS, AND NORMALIZED COUNTS TO OUTPUT A TABLE
#WITH THE TOP ANNOTATION ALREADY ATTACHED
#' @export
GO_definitions_to_table <- function(comparison_table, normalized_counts, annotations){
  #annotations <- GO_annotations(comparison_table, normalized_counts)
  #list_of_genes <- normalized_counts$Name #this is the issue
  list_of_genes <- rownames(normalized_counts)
  comparison_table$Description <- NA
  for(gene in list_of_genes){
    print(gene)
    description <- GO_pvalue_finder(gene, annotations, comparison_table)
    row_num <- which(rownames(comparison_table) == gene)
    comparison_table$Description[row_num] = description

  }
  return(comparison_table)
}

#' @export
Go_anno_selection <- function(gene, annotation_df, num){
  sorted_annotations <- annotation_df[order(annotation_df$pvalue), ]
  table <- sorted_annotations %>% filter(grepl(gene, geneID))
  selected_table <- table %>% slice_min(order_by = p.adjust, n = num)
  return(selected_table)
}

#Go annotations lists
#' @export
entrez_id_list <- function(raw_counts, normalized_counts){
  list_of_genes <- raw_counts$Name
  genes <- rownames(normalized_counts)
  #Getting the entrez_ids
  all_entrez_ids <- bitr(list_of_genes, fromType = 'SYMBOL', toType = 'ENTREZID',  OrgDb = org.Hs.eg.db)
  entrez_ids <- bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID',  OrgDb = org.Hs.eg.db)
  return(entrez_ids$ENTREZID)
}



