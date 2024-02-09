rbind_df_list <- function(lst, new_col = "dataset") {
  lst_with_newcol <- mapply(x=names(lst), y=lst, FUN = function(x, y){
    y[,new_col] <- x
    y
  }, SIMPLIFY = FALSE)
  Reduce(rbind, lst_with_newcol)
}