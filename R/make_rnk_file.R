#' Make .rnk file
#' Input: txt file
#' Output: writes out .rnk file
#' 
#' make .rnk file (ie. for GSEA preranked)
#' @param file filename of some txt file
#' @param column column number of orig file you want to rank (>1)
#' @export
#' 

make_rnk_file = function(file, column){
  
  data = read.delim(file)
  data = data[, c(1, column)]
  name = sub(".txt", "", file)
  write.table(data, paste0(name, "_", colnames(data)[2],".rnk"), quote = F, sep = "\t", row.names = F, col.names = F)
}