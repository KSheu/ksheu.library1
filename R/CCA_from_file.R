#' Do Canonical correlation analysis
#' 
#' Input: matrix A (ex. rows as genes, cols as samples)
#' Input: matrix B (ex. rows as genes, cols as samples)
#' Output: writes to file scores and loadings of decomposed matrices that maximize the correlation between A and B, as well as p-values for each component
#' 
#' WARNING - script will internally convert all colnames to alphanumeric only.
#' 
#' Returns cca object, number of significant canonical dimensions (can be smaller than #of vars in smaller set)
#' 
#' @param file1 file for A matrix (.txt file)
#' @param file2 file for B matrix
#' @param rank rank for the decompostiton
#' @param pval_method default = "Wilks". Method for testing signficance of components (Wilks, Hotelling,Pillai,Roy)
#' @param standardize default = T
#' 
#' @export
#' 

# based on https://stats.idre.ucla.edu/r/dae/canonical-correlation-analysis/
# install.packages("CCA")
# install.packages("CCP")

CCA_from_file = function(file1, file2, rank, pval_method){
  require(ggplot2)
  require(GGally)
  require(CCA)
  require(CCP) #for testing significance of canonical vars
  
  #file1 = ("F:/SKCM_expression/GDSC_Cell_line_RMA_proc_basalExp_melanoma.txt")
  #file2 = ("F:/dependency_ceres/avana_ceres_gene_effects_MELANOMA.txt")
  
  #melanoma drug response/RNAi CCA on methylation and/or expression data
  data1 = read.delim(file1)
  data2 = read.delim(file2)
  
  f1 <- basename(file1)
  f2 <- basename(file2)
  f1 <- sub(f1, pattern = ".txt", replacement = "")
  f2 <- sub(f2, pattern = ".txt", replacement = "")
  
  rownames(data1) = data1[,1]
  rownames(data2) = data2[,1]
  
  colnames(data1) = toupper(gsub("[^[:alnum:]]","",colnames(data1)))
  colnames(data2) = toupper(gsub("[^[:alnum:]]","",colnames(data2)))
  #colnames(data2) = gsub("_SKIN", "", colnames(data2))
  
  t.data1 = data.frame(t(data1[,-1]))
  t.data2 = data.frame(t(data2[,-1]))
  
  
  samples = intersect_all(rownames(t.data1), rownames(t.data2))
  t.data1 = t.data1[(rownames(t.data1) %in% samples), ]
  t.data2 = t.data2[(rownames(t.data2) %in% samples), ]
  
  t.data1 = na.omit(t.data1)
  t.data2 = na.omit(t.data2)
  
  t.data1 = t.data1[order(rownames(t.data1)), ]
  t.data2 = t.data2[order(rownames(t.data2)), ]
  
  
  #cc1 <- cc(t.data1[,c(1:10)], t.data2[,c(1:10)])
  cc1 <- cc(t.data1, t.data2)
  
  # display the canonical correlations
  cc1$cor
  # raw canonical coefficients
  cc1[3:4]
  
  # compute canonical loadings -- this runs within cc() function call
  cc2 <- cc1$scores
  
  # display canonical loadings
  cc2[3:6]
  
  
  # tests of canonical dimensions
  rho <- cc1$cor
  ## Define number of observations, number of variables in first set, and number of variables in the second set.
  n <- dim(t.data1)[1]
  #p <- length(t.data1[,c(1:10)])
  p <- length(t.data1)
  #q <- length(t.data2[,c(1:10)])
  q <- length(t.data2)
  
  
  ## Calculate p-values using the F-approximations of different test statistics (Wilks, Hotelling,Pillai,Roy):
  pvals <- p.asym(rho, n, p, q, tstat = "Wilks")
  p.asym(rho, n, p, q, tstat = "Wilks")
  
  
  #compute the standardized canonical coefficients for easier comparison among vars if not same sd
  #can interpret in same manner as standardized regression coefficients
  # standardized psych canonical coefficients diagonal matrix of file1 sd's
  s1 <- diag(sqrt(diag(cov(t.data1))))
  s1 %*% cc1$xcoef
  # standardized acad canonical coefficients diagonal matrix of file2 sd's
  s2 <- diag(sqrt(diag(cov(t.data2))))
  s2 %*% cc1$ycoef
  
  #Writing tables out
  colnames(cc1$scores$corr.X.xscores) <- paste0("X", c(1:ncol(cc1$scores$corr.X.xscores)))
  write.table(x = cc1$scores$corr.X.xscores, file = paste(f1, f2, "corr.X.xscores.txt", sep = "_"), quote = F, row.names = T, sep = '\t') 
  
  colnames(cc1$scores$corr.X.yscores) <- paste0("X", c(1:ncol(cc1$scores$corr.X.yscores)))
  write.table(x= cc1$scores$corr.X.yscores, file = paste(f1, f2, "corr.X.yscores.txt", sep = "_"), quote = F, row.names = T, sep = '\t')
  
  colnames(cc1$scores$corr.Y.xscores) <- paste0("X", c(1:ncol(cc1$scores$corr.Y.xscores)))
  write.table(x=cc1$scores$corr.Y.xscores, file = paste(f1, f2, "corr.Y.xscores.txt", sep = "_"), quote = F, row.names = T, sep = '\t')
  
  colnames(cc1$scores$corr.Y.yscores) <- paste0("X", c(1:ncol(cc1$scores$corr.Y.yscores)))
  write.table(x=cc1$scores$corr.Y.yscores, file = paste(f1, f2, "corr.Y.yscores.txt", sep = "_"), quote = F, row.names = T, sep = '\t')
  
  colnames(cc1$scores$xscores) <- paste0("X", c(1:ncol(cc1$scores$xscores)))
  write.table(x=cc1$scores$xscores, file = paste(f1, f2, "xscores.txt", sep = "_"), quote = F, row.names = T, sep = '\t')
  
  colnames(cc1$scores$yscores) <- paste0("X", c(1:ncol(cc1$scores$yscores)))
  write.table(x=cc1$scores$yscores, file = paste(f1,f2, "yscores.txt", sep = "_"), quote = F, row.names = T, sep = '\t')
  
  #write.table(unlist(cc1), file = paste(file1, "_CORRELATIONS.txt"), quote = F, row.names = F) 
  #write.table(unlist(cc2), file = paste(file1, "_LOADINGS.txt"), quote = F, row.names = F)
  #write.table(unlist(pvals), file = paste(file1, "_PVALS.txt"), quote = F, row.names = F)
  
  plt.cc(cc1, var.label = T, int = 5)
  #plt.indiv(cc1)
  #plt.var(cc1)

}
