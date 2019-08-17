#' Implements tensor factorization methods, Tucker and PARAFAC/canonical polyadic decomp
#' 
#' Tucker:  If there is no truncationin one of the modes, then this is the same as the MPCA,mpca.  
#' If there is no truncation in all the modes, then this is the same as the HOSVD,hosvd. 
#' 
#' Input: Multidimensional array that gets converted to as.tensor
#' Output: Decomposed tensor object
#' 
#' @param matrix matrix that gets coverted to multidimensional array that get converted through as.tensor(arr)
#' @param dims dimensions to fold the input matrix
#' @param method "tucker" or "cp"
#' @param z_labs labels for the extra dimensions created, rows and colnames are the x and y axis labs
#' @param rank default = 3, number (if CP) or vector of numbers (if tucker) that is <= the ranks of the tensor
#' @param plot plot heamtmaps of eigenvectors
#' 
#' 
#' @export
#' 

# INPUT
# matrix = read.delim("F:/atac_melanoma/atac_melanoma_hg38/ltg_melanoma_ATAC_normalizedcounts.induced_log2.txt")
# matrix = read.delim("F:/atac_melanoma/atac_melanoma_hg38/ltg_melanoma_ATAC_normalizedcounts.txt")
# matrix = matrix[, -22] #rm 399_1wk
# matrix = matrix[, c("peak","X262_d0","X262_2hr","X262_LTg","X262_TNF",
#                 "X308_d0", "X308_2hr", "X308_LTg","X308_TNF",
#                 "X399_d0","X399_2hr","X399_LTg", "X399_TNF",
#                 "X3998MEL_d0","X3998MEL_2hr", "X3998MEL_LTg","X3998MEL_TNF",
#                 "X257A2_d0" ,"X257A2_2hr" , "X257A2_LTg","X257A2_TNF",
#                 "X370_d0","X370_2hr","X370_LTg","X370_TNF",
#                 "X381_d0","X381_2hr","X381_LTg","X381_TNF",
#                 "X410_d0","X410_2hr","X410_LTg","X410_TNF")]
# rownames(matrix) = matrix$peak
# 
# dims = c(nrow(matrix), 4, 8)
# method = "tucker"
# xlabs = rownames(matrix)
# ylabs = c('0hr','2hr','ltg','tnf')
# zlabs = c("X262","X308","X399","X3998MEL","X257A2","X370","X381","X410")
# ranks = c(10,4,8)
# plot = F
# savename="ltg_melanoma_ATAC_normalizedcounts"
# 
# tensor_decomp(matrix, dims, method, xlabs = xlabs, ylabs =ylabs, zlabs = zlabs, 
#               ranks = ranks, savename="ltg_melanoma_ATAC_normalizedcounts")

tensor_decomp = function(matrix, dims, method, xlabs = rownames(matrix), ylabs =colnames(matrix)[-1], zlabs = NULL, ranks = A@modes, savename="TEST", plot=F){
  
  require(rTensor);require(pheatmap)
  
  rownames(matrix) = matrix[,1] #genenames into rownames
  matrix = as.matrix(matrix[,-1])
  
  # convert to tensor
  arr = k_fold(mat = matrix, m = 1, modes = dims)
  arr
  
  # A = as.tensor(arr, drop = FALSE)
  A = arr
  A@modes
  A@num_modes
  dim(A)
  A@data[,,8]
  
  if(method=="tucker"){
    ## Tucker decomposition
    print("Running Tucker Decomp")
    tucker_decomp <- rTensor::tucker(A, ranks = ranks)
    str(tucker_decomp)
    
    
    A_approx <- ttl(tucker_decomp$Z, tucker_decomp$U, 1:3)
    fnorm(A - A_approx) / fnorm(A)
    tucker_decomp$norm_percent
    
    plot(tucker_decomp$all_resids)
    
    U.1 = data.frame(tucker_decomp$U[[1]]);rownames(U.1) = (xlabs)
    U.2 = data.frame(tucker_decomp$U[[2]]);rownames(U.2) = ylabs
    U.3 = data.frame(tucker_decomp$U[[3]]);rownames(U.3) = zlabs
    
    if (plot==T){
      pheatmap(U.1, cluster_cols = F, cluster_rows = T, scale="row", clustering_method = "ward.D2")
      pheatmap(U.2, cluster_cols = F, cluster_rows = F, scale="none")
      pheatmap(U.3, cluster_cols = F, cluster_rows = F, scale="none")
    }
    
    write.table(cbind(variable = rownames(U.1), U.1), paste0(savename,"_U.1_",method,".txt"), quote=F,row.names = F, sep="\t")
    write.table(cbind(variable = rownames(U.2), U.2), paste0(savename,"_U.2_",method,".txt"), quote=F,row.names = F, sep="\t")
    write.table(cbind(variable = rownames(U.3), U.3), paste0(savename,"_U.3_",method,".txt"), quote=F,row.names = F, sep="\t")
  }
  
  
  if(method=="cp"){
    ## cp decompostion
    print("Running Canonical Polyadic Decomp")
    cp_decomp <- cp(A, num_components = rank, max_iter = 50)
    str(cp_decomp)
    
    cp_decomp$conv #did it converge
    plot(cpD$all_resids)
    
    A_approx <- ttl(cp_decomp$Z, cp_decomp$U, 1:3)
    fnorm(A - A_approx) / fnorm(A)
    cp_decomp$norm_percent
  }
  
  
} 
