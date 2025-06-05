#Function for creating new metadata table with new columns that helps hone down the data
#Takes in metadata table, 1st column you're interested in combining, 2nd column you're interested in combining
#' @export
new_metadata_table <- function(metadata_table, column1, column2, new_title){
  metadata_table[[new_title]] <- paste(metadata_table[[column1]], metadata_table[[column2]], sep = "")
  return(metadata_table)
}



#Function to create tables specifically associated to one demographic
#GENE TABLES
#' @export
gene_tables <- function(normcounts, table, sample_col_name, column_name, column_value, sample_char){
  specific_table <- table[table[[column_name]] == column_value,]
  samples <- specific_table[[sample_col_name]]
  #gene_list_table <-  normcounts %>% dplyr::select('Name', intersect(names(.), samples))
  gene_list_table <- normcounts[, intersect(colnames(normcounts), samples), drop = FALSE]
  return(gene_list_table)
}
