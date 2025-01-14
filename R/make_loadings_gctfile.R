#'Format PCA loadings for IGV
#'
#'Outputs gct file in format for IGV
#'
#'@param file loadings file (output of PCA_from_file)
#'@param ncomps number of components
#'
#'
#'@export

make_loadings_gctfile = function(file, ncomps = 5){
  loadings = read.delim(file, header =T)
  loadings = loadings[, c(0:ncomps+1)]
  loadings = cbind ("Loading" = loadings[,1], loadings)
  name=sub(".txt","",file)
  out_file = paste0(name, "_PC1-5.gct")
  
  cat(paste0("#1.2\n",dim(loadings)[1],"\t",ncomps,"\n"), file=out_file)
  write.table(loadings, out_file,sep = "\t", row.names =F, quote = F, append = T)
}
