#### Script containing general use functions

#### Manipulate Tibbles ####
## Turn tibble into a matrix with rownames from a column
tbl_to_matrix <- function(x, row){
  mat <- as.matrix(select(x, -row))
  rownames(mat) <- pull(x, !!row)
  return(mat)
}

## Extract a list of variables by group from a data frame
tbl_var_to_list <- function(tbl, var, group=NULL){
  if (is.null(group)){
    group <- group_vars(tbl)[1]
  }
  items <- pull(tbl, !!var)
  labs <- pull(tbl, !!group)
  sapply(unique(labs), function(x){items[labs == x]}, simplify = FALSE)
}

#### Set NA type values ####
#set all NA values to a generic value
set_na <- function(x, val){
  x[is.na(x)] <- val
  return(x)
}

#set all inf values to a generic value
set_inf <- function(x, val){
  x[is.infinite(x)] <- val
  return(x)
}

#set all NaN values to a generic value
set_nan <- function(x, val){
  x[is.nan(x)] <- val
  return(x)
}
