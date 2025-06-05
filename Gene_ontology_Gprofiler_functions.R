#GENE ONTOLOGY GPROFILER FUNCTIONS
#By: Paul Parodi
#Last updated: 6/4/2025


#Returning the annotations table for gprofiler2
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


#Making a bar plot based off those gprofiler2 annotations
#' @export
gprofiler_barplot <- function(gplot ,num1, num2){
  top_terms <- gplot$result[num1:num2, ]
  gbarplot <- ggplot(top_terms, aes(x = reorder(term_name, --log10(p_value)), y = -log10(p_value), fill = source)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "Term", y = "-log10(p-value)")
  return(gbarplot)
}


#Mapping the gene ids to the names in the table that you are interested in
gene_gprofile_mapping <- function(Table, id_col, id_to_name_map, new_col = "gene_names"){
  #id_col: Name of the column in df with comma-seperated gene IDs (as string)
  #id_to_name_map: Named vector (names are gene IDs, values are gene names) or data frame with
  #id and name columns
  #If id_to_name map is a dataframe, conver to named vector
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

#Mapping the entrez ids to the genes and then pasting them on the table.
gene_entrez_gprofile_mapping <- function(Table, id_col, norm_counts){
  genes <- rownames(norm_counts)
  E_ids <- entrez_id_list(norm_counts, norm_counts)
  id_to_name <- setNames(genes, E_ids)
  Final_table <- gene_gprofile_mapping(Table, id_col, id_to_name)
  return(Final_table)
}


