#' Partial Least Squares and predict second dataset
#'
#' Builds PLS model from training dataset and Writes out loadings.
#' Predicts classifications of test dataset and returns dataframe of predictions
#'
#' @param file file for X matirx
#' @param sample.names Vector of sample names in X matrix
#' @param response.values Vector of response values in same order matching sample.names
#' @param comps number of components to compute
#' @param scale default=T
#' @param ind.names Labels the samples, default =F
#' @param file2 file for test data matrix
#' @param sample.names2 Vector of sample names in 2nd dataset, if needed
#' @param response.values Vector of response values in same order matching sample.names2, if available
#' 
#' @export
#'


PLSDA_from_file_and_predict_second_dataset = function(file, file2, sample.names, response.values, sample.names2=NULL, response.values2=NULL, comps = 3, scale = F, ind.names = F){
  require(mixOmics)
  data = read.table(file, sep='\t',header=T,stringsAsFactors=FALSE, quote = "")
  data = data[rowSums((data[, -1] == 0)) < ncol(data[-1]), ] #remove genes with no variance
  data2 = read.table(file2, sep='\t',header=T,stringsAsFactors=FALSE, quote = "")
  
  common.genes = intersect_all(data[,1], data2[,1])
  # data = data[data[,1] %in% common.genes, ]
  # data2 = data2[data2[,1] %in% common.genes, ]
  data=data1[match(common.genes,data[,1]),]
  data2=data2[match(common.genes,data2[,1]),]
  
  rownames(data) = data[,1]
  t.data = data.frame(t(data[,-1])) 
  t.data$group = (response.values[match(rownames(t.data), sample.names)])
  pls.fit = plsda(X = t.data[,-ncol(t.data)], Y = (t.data[,ncol(t.data)]), scale = scale, ncomp = comps)
  plotIndiv(pls.fit, legend = T,ind.names = ind.names)
  write.table(as.data.frame(pls.fit$loadings$X), paste0(gsub(".txt", "", file), "_PLSDA_Xloadings.txt"), sep = "\t", row.names = T, quote = F)

  rownames(data2) = data2[,1]
  t.data2 = data.frame(t(data2[,-1])) 
  
  test.predict <- predict(pls.fit, t.data2, method = "max.dist")
  prediction <- as.data.frame(test.predict$class$max.dist[, comps])
  colnames(prediction)[1] <-"prediction"
  if (!is.null(sample.names2)) {
    prediction$actual = (response.values2[match(rownames(prediction), sample.names2)])
    print(table(prediction$prediction, droplevels(prediction$actual)))
  }
  return(prediction)
  
}



