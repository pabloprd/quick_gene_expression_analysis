#Selecting the top values of a table ordered by the column name you enter and returns the number you specified
#' @export
top_selection <- function(table, column_name, number){
  table1 <- table %>% slice_max(order_by = .data[[column_name]], n = number)
  return(table1)
}


#Selecting the bottom values of a table ordered by the column name you enter and returns the number you specified
#' @export
bottom_selection <- function(table, column_name, number){
  table1 <- table %>% slice_min(order_by = .data[[column_name]], n = number)
  return(table1)
}
