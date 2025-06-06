# SELECTION OF TOP AND BOTTOM VALUES OF A DATAFRAME
# By: Paul Parodi
# Last updated: 6/6/2025

#' Selecting top values
#'
#' Selecting the top values of a table ordered by the column name you enter and returns the number you specified
#'
#' @param table  A data frame
#' @param column_name The column name you are ordering the data frame by
#' @param number The  number of values you want the function to return from the top values.
#'
#' @return An organized, miniature dataframe of top values from the dataframe you selected
#'
#' @examples
#' # top_selection(dataframe, 'P.value', 5)
#'
#' @export
top_selection <- function(table, column_name, number){
  table1 <- table %>% slice_max(order_by = .data[[column_name]], n = number)
  return(table1)
}


#' Selecting bottom values
#'
#' Selecting the bottom values of a table ordered by the column name you enter and returns the number you specified
#'
#' @param table  A data frame
#' @param column_name The column name you are ordering the data frame by
#' @param number The  number of values you want the function to return from the bottom values.
#'
#' @return An organized, miniature dataframe of bottom values from the dataframe you selected
#'
#' @examples
#' # bottom_selection(dataframe, 'P.value', 5)
#'

#' @export
bottom_selection <- function(table, column_name, number){
  table1 <- table %>% slice_min(order_by = .data[[column_name]], n = number)
  return(table1)
}
