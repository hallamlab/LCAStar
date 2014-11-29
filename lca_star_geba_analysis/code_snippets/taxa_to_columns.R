# function to split taxonomy into multiple columns by depth
taxa_to_columns <- function(in_vec, cols=4) {
  # split by;
  split_list <- strsplit(as.character(in_vec), ";")
  # apply matrix
  funny <- function(l,cols) {
    l[1:cols]
  }
  mat <- sapply(split_list, funny, cols=cols)
  mat <- t(mat)
  # add header
  header <- paste("level", 1:ncol(mat), sep="_")
  colnames(mat) <- header
  return(mat)
}



