#GENE ONTOLOGY GPROFILER FUNCTIONS
#By: Paul Parodi
#Last updated: 6/6/2025




#' gprofiler Gene ontology
#'
#' Makes an annotation table using gprofiler.
#'
#' @param eids list of gene eids
#'
#' @return a data matrix of annotated genes using the gprofiler database
#'
#' @examples
#' # gprofiler(list_of_eids)
#'
#' @export
gprofiler <- function(eids){
  gprofiletable <- gost(
    query = eids,
    organism = "hsapiens",         # For human
    correction_method = "bonferroni",     # Multiple testing correction
    sources = c("GO:BP"),# Databases to query
    evcodes = TRUE
  )
  return(gprofiletable)
}

#Vestigial function as of now.
gprofiler_barplot <- function(gplot ,num1, num2){
  top_terms <- gplot$result[num1:num2, ]
  gbarplot <- ggplot(top_terms, aes(x = reorder(term_name, --log10(p_value)), y = -log10(p_value), fill = source)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "Term", y = "-log10(p-value)")
  return(gbarplot)
}


#' Mapping the gene ids to gprofiler terms
#'
#' Turns your dataframe into an annotated dataframe with gprofiler gene terms.
#'
#' @param Table Your dataframe of normalized counts, or any dataframe that has genes as the index and the gene ids
#' @param id_col Name of the column in df with comma-seperated gene IDs (as string)
#' @param id_to_name_map Named vector (names are gene IDs, values are gene names) or data frame with
#' @param new_col The new column name that you want for the gene names
#'
#' @return a barplot of the counts in relation to each gene ontology term.
#'
#' @examples
#' # gene_gprofile_mapping(dataframe, geneID, id_to_name_mapping_list, 'gene_names')
#'
#' @export
gene_gprofile_mapping <- function(Table, id_col, id_to_name_map, new_col = "gene_names"){
  if (is.data.frame(id_to_name_map)) {
    id_to_name_map <- setNames(as.character(id_to_name_map[[2]]), as.character(id_to_name_map[[1]])
    )
  }
  map_ids <- function(id_string) {
    ids <- unlist(strsplit(id_string, ","))
    names <- id_to_name_map[ids]
    names[is.na(names)] <- ids[is.na(names)]  # fallback to ID if name not found
    paste(names, collapse = "/")
  }

  Table[[new_col]] <- vapply(Table[[id_col]], map_ids, character(1))
  return(Table)
}

#' Mapping genes to entrez ids
#'
#' Turns your dataframe into an annotated dataframe with gprofiler gene terms.
#'
#' @param Table Your dataframe of normalized counts, or any dataframe that has genes as the index and the gene ids
#' @param id_col Name of the column in df with comma-seperated gene IDs (as string)
#' @param id_to_name_map Named vector (names are gene IDs, values are gene names) or data frame with
#' @param new_col The new column name that you want for the gene names
#'
#' @return a dataframe of gprofiler annotated terms.
#'
#' @examples
#' # gene_entrez_gprofile_mapping(dataframe, column_of_ids, normalized counts)
#'
#' @export
gene_entrez_gprofile_mapping <- function(Table, id_col, norm_counts){
  genes <- rownames(norm_counts)
  E_ids <- entrez_id_list(norm_counts, norm_counts)
  id_to_name <- setNames(genes, E_ids)
  Final_table <- gene_gprofile_mapping(Table, id_col, id_to_name)
  return(Final_table)
}


