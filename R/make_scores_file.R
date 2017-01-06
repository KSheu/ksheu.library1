#' Makes scores file for IGV
#'
#' Writes out minimized scores file 
#' @param file Score file (output of PCA_from_file)
#' @param ncomps number of PCs wanted
#' 
#' @export
#'
make_scores_file = function(file, ncomps=5){
  data = read.delim(file, header = T)
  data = data[, c(0:ncomps+1)]
  colnames(data)[1] <- "Sample" 
  name=sub(".txt","",file)
  write.table(data, paste(name, "_PC1-5.txt", sep = '') ,sep = "\t", row.names =F, quote = F)
}