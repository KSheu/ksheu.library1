#' Implements tensor factorization methods, Tucker and PARAFAC/canonical polyadic decomp
#' 
#' Tucker:  If there is no truncationin one of the modes, then this is the same as the MPCA,mpca.  
#' If there is no truncation in all the modes, then this is the same as the HOSVD,hosvd. 
#' 
#' Input: Multidimensional array that gets converted to as.tensor
#' Output: Decomposed tensor object
#' 
#' @param arr multidimensional array that get converted through as.tensor(arr)
#' @param method "tucker" or "cp"
#' @param rank default = 3, number (if CP) or vector of numbers (if tucker) that is <= the ranks of the tensor
#' 
#' 
#' 
#' @export
#' 



tensor_decomp = function(arr, method, rank = tnsr@modes){
  
  require(rTensor)
  
  
  # convert to tensor
  A = as.tensor(arr, drop = FALSE)
  A@modes
  A@num_modes
  dim(A)
  A@data
  
  if(method=="tucker"){
    ## Tucker decomposition
    tucker_decomp <- rTensor::tucker(A, ranks = rank)
    str(tucker_decomp)
    
    
    A_approx <- ttl(tucker_decomp$Z, tucker_decomp$U, 1:3)
    fnorm(A - A_approx) / fnorm(A)
    tucker_decomp$norm_percent
    
    U.1 = data.frame(tucker_decomp$U[[1]]);rownames(U.1) = make.unique(rownames(data))
    U.2 = data.frame(tucker_decomp$U[[2]]);rownames(U.2) = unique(info$type)[-1]
    U.3 = data.frame(tucker_decomp$U[[3]]);rownames(U.3) = c(0,1,3,8)
    
  }
  
  
  if(method=="cp"){
    ## cp decompostion
    cp_decomp <- cp(A, num_components = rank, max_iter = 50)
    str(cp_decomp)
    
    cp_decomp$conv #did it converge
    plot(cpD$all_resids)
    
    A_approx <- ttl(tucker_decomp$Z, tucker_decomp$U, 1:3)
    fnorm(A - A_approx) / fnorm(A)
    cp_decomp$norm_percent
  }
  
  
}