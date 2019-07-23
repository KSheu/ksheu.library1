#' Do Independent components analysis
#' The input matrix X (input as genes in rows, samples in columns) is a linear combination of non-Gaussian independent components.
#' As opposed to PCA, which assumes generally Gaussian distributed data.
#' Assume X = SA, where S are the independent components.
#' We want to find another matrix W , where XW = S, to unmix 
#' 
#' Input: matrix X (ex. rows as genes, cols as samples)
#' Output: plots only currently
#'  
#' @param data file for the input matrix X (.txt or .RDS file)
#' @param data.labels labels for samples with two columns: samp.names and type
#' @param comps number of components to be extracted, defualt = 5
#'  
#' @export
#' 


# data.labels = read.delim("F:/project_bioinfo/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_brain_patient_matched_regions.txt", 
#                   stringsAsFactors = F, nrows = 2)
# data = read.delim("F://project_bioinfo/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_brain_patient_matched_regions_format_log2_top75.txt")
# comps = 10
ICA_from_file = function(data, data.labels, comps = 5){
  
  require(fastICA)
  
  data = read.delim(data)
  ica<-fastICA((data[,-1]), n.comp= comps, alg.typ = "parallel", fun = "logcosh", 
               alpha = 1,method = "C", row.norm = FALSE, maxit = 200,tol = 0.0001, verbose = TRUE)
  
  pairs(ica$S, col=rainbow(3)[data[,1]])
  
  plot(ica$S[,1], ica$S[,1], col=rainbow(3)[data[,1]], xlab="Comp 1", ylab="Comp 2")
  
  X = ica$X #Pre-processed data
  
  A = t(ica$A) #ica loadings
  rownames(A) = colnames(data)[-1]
  A = as.data.frame(A)
  
  
  S = ica$S #ICA components
  rownames(S) = data[,1]
  
  
  K = ica$K #PCA loadings
  pca.scores = ica$X %*% ica$K
  
  W = ica$W # estimated un-mixing matrix
  
  plot(ica$X %*% ica$K, main = "PCA components")
  plot(A, main = "ICA components")
  
  A$type = data.labels[,2][match(rownames(A), gsub("-",".", data.labels[,1]))]
  ggplot(A, aes(A$V1, A$V3, color = type))+geom_point(size=4)+theme_bw(base_size = 16)+xlab("Component 1")+ylab("Component 3")
  ggplot(A, aes(A$V2, A$V4, color = type))+geom_point(size=4)+theme_bw(base_size = 16)+xlab("Component 2")+ylab("Component 4")
  
}
