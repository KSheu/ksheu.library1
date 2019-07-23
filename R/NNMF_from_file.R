#' Nonnegative matrix factorization
#' 
#' Input: matrix A (ex. rows as genes, cols as samples)
#' Output: writes to file decomposed matrices W (loadings) and H (scores) that minimize the loss between the original and decomposed matrices
#' 
#' Returns nnmf object
#' 
#' @param file file for A matrix (.txt or .RDS file)
#' @param rank rank for the decompostiton
#' @param method 'scd' (seq. coord-wise descent) or 'lee' (Lee's algo)
#' @param loss 'mse' (mean squared error) or 'mkl' (mean KL divergence) 
#' for sc data probably better to use KL divergence? (Frobenius norm is standard)
#' @param impute default=F Simple NNMF or impute 0 values by first setting them all to NA
#' @param verbose 0-2  
#' 
#' @export
#' 


NNMF_from_file <- function(file, rank = 5, method = "scd", loss = "mkl", impute = F, verbose = 2){
  require(NNLM)
  
  if(grepl(".txt$", file)==TRUE){
    matrix = read.delim(file)
    rownames(matrix) = matrix[,1]
    matrix = as.matrix(matrix[,-1])
  }else{
    matrix <-readRDS(file)
    matrix = as.matrix(matrix)
  }
  
  if(impute==T){
    matrix[matrix==0]<-NA
  }
  
  system.time(decomp.nmf <- nnmf(matrix, k = rank, method = method, loss = loss, verbose = verbose, check.k = F));
  print(decomp.nmf)
  W <-decomp.nmf$W
  H <-decomp.nmf$H
  decomp.nmf.hat = W%*%H
  
  W.save = cbind(gene=rownames(W), W)
  H.save = cbind(samples = colnames(H), t(H))
  
  #save data
  if(grepl(".txt$", file)==TRUE){
    name=sub(".txt","",file)
  }else{
    name=sub(".RDS","",file)
  }
  
  savename=paste(name,"_NNMF_Wloadings.txt",sep='');
  write.table(W.save,savename,sep='\t',row.names=FALSE,quote=FALSE);
  savename=paste(name,"_NNMF_Hscores.txt",sep='');
  write.table(H.save,savename,sep='\t',row.names=FALSE,quote=FALSE);
  
  return(decomp.nmf)
  #with decomp.nmf can predict on new samples
  # newH <- predict(decomp.nmf, matrix2, which = 'H');
  # str(newH)
 
}
