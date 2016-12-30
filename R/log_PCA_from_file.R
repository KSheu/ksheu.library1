#' PCA from file after taking log2
#' 
#' Reads file with samples in columns and variables in rows. Takes log2(x+1) and does PCA. Writes to file scores, loadings, eigenvalues.
#' 
#' @param file Filepath/filename of data matrix with no row numbering
#' @param center default=T
#' @param scale default=F
#' 
#' @importFrom stats prcomp screeplot
#' @importFrom utils read.delim read.table write.table
#' 
#' @export
#' 

log_PCA_from_file=function(file,center=TRUE,scale=FALSE) {
  
  data=read.delim(file, header = T)
  
  #log2(x+1)
  data = cbind("genename" = data[,1], log(data[,-1]+1, 2))#log2
  data$genename = as.character(data$genename)
  
  #remove rows that are all 0
  data= data[rowSums((data[,-1]==0))<ncol(data[-1]),]
  
  #t.data = t(data) #if genenames inrownames
  t.data=t(data[,-1]) ##subtract off the gene name
  pca<-prcomp(t.data,scale=scale,center=center);
  pca_scores=pca$x
  pca_scores=cbind("Score"=rownames(pca_scores),pca_scores)
  pca_loadings=pca$rotation
  pca_loadings=cbind("Loading"=data[,1],pca_loadings)#if genenames not in rownames
  #pca_loadings=cbind("Loading"=rownames(pca_loadings),pca_loadings)#if genenames in rownames
  pca_evalues=pca$sdev
  
  #save data
  name=sub(".txt","",file)
  savename=paste(name,"_prcomp_scores.txt",sep='');
  write.table(pca_scores,savename,sep='\t',row.names=FALSE,quote=FALSE);
  savename=paste(name,"_prcomp_loadings.txt",sep='');
  write.table(pca_loadings,savename,sep='\t',row.names=FALSE,quote=FALSE);
  savename=paste(name,"_prcomp_sdev.txt",sep='');
  write.table(pca_evalues,savename,sep='\t',row.names=FALSE,quote=FALSE);
  print(summary(pca))
  screeplot(pca)
  
}
