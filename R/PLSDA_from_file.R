#'Partial Least Squares from file
#'
#' Writes out loadings
#'
#' @param file file for X matirx
#' @param sample.names Vector of sample names in X matrix
#' @param response.values Vector of response values in same order matching sample.names
#' @param comp number of components to compute
#' @param scale default=T
#' 
#' @export
#'


PLSDA_from_file = function(file, sample.names, response.values, comps = 2, scale = F, ind.names = F){
  require(mixOmics)
  require(dummies)
  data = read.table(file, sep='\t',header=T,stringsAsFactors=FALSE, quote = "")
  data = data[rowSums((data[, -1] == 0)) < ncol(data[-1]), ] #remove genes with no variance
  rownames(data) = data[,1]
  t.data = data.frame(t(data[,-1])) 
  t.data$group = (response.values[match(rownames(t.data), sample.names)])
  pls.fit = plsda(X = t.data[,-ncol(t.data)], Y = (t.data[,ncol(t.data)]), scale = scale, ncomp = comps)
  plotIndiv(pls.fit, legend = T,ind.names = ind.names)
  write.table(as.data.frame(pls.fit$loadings$X), paste0(gsub(".txt", "", file), "_PLSDA_Xloadings.txt"), sep = "\t", row.names = T, quote = F)
}