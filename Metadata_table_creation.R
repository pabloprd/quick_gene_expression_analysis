# Metadata Table functions
#By Paul Parodi
#Last updated: 6/9/2025


#' Combines two categorical columns into another metadata column
#'
#' This function combines two categorical columns into a new metadata column for pairwise comparisons
#'
#' @param metadata_table  metadata table with all the feature columns
#' @param column1 Column 1 that you wish to combine with another column
#' @param column2 Column 2 that you wish to combine with another column
#' @param new_title New title of the column you are making
#'
#'
#' @return a dataframe of the metadata with an extra column of combined categorical values you created
#'
#' @examples
#' #new_metadata_table(metadata, 'Ancestry', 'Sex', 'Sex_and_Ancestry')
#'
#' @export
new_metadata_table <- function(metadata_table, column1, column2, new_title){
  metadata_table[[new_title]] <- paste(metadata_table[[column1]], metadata_table[[column2]], sep = "")
  return(metadata_table)
}



#Redundant code that I may use later
gene_tables <- function(normcounts, table, sample_col_name, column_name, column_value, sample_char){
  specific_table <- table[table[[column_name]] == column_value,]
  samples <- specific_table[[sample_col_name]]
  #gene_list_table <-  normcounts %>% dplyr::select('Name', intersect(names(.), samples))
  gene_list_table <- normcounts[, intersect(colnames(normcounts), samples), drop = FALSE]
  return(gene_list_table)
}
