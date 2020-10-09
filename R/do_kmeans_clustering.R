#' Do k-means clustering on data matrix, genes in rownames, samples in columnnames
#' Input: Numeric only data matrix with gene names in rownames
#' Output: prints clustered heatmap
#' 
#' 
#' @param data_matrix matrix or dataframe to cluster, numeric only
#' @param k_clusters cluster number, default = 5
#' @param colseps vector of numbers for column separation gaps
#' @param annotation Annotation data frame with genenames in rownames
#'  
#' @export
#' 


do_kmeans_clustering = function(data_matrix, k_clusters = 5, colseps = NULL, annotation_row=NULL,annotation_col=NULL,
                                cluster_rows=F, cluster_cols=T,
                                show_rownames = T, show_colnames = F,
                                breaks=c(-4,seq(-1.5,1.5,length=100),4)){
  
  require(pheatmap)
  
  set.seed(1)
  scaled = t(scale(t(data_matrix)))
  
  kmcluster <- kmeans(scaled,iter.max = 1000, centers = k_clusters, nstart = 1) #do kmeans
  # kmcluster$cluster
  mat <- cbind(data_matrix, cluster = kmcluster$cluster)
  mat <- mat[order(mat$cluster), ] #put genes in order of cluster
  # mat <- mat[which(mat$cluster==5), ]
  #make row_seps
  count <- 0
  for(i in 1:k_clusters){
    count[i] <- length(kmcluster$cluster[kmcluster$cluster == i])
  }
  rowseps <- cumsum(count)
  
  p = pheatmap(mat[, -ncol(mat)], cluster_rows=cluster_rows,cluster_cols=cluster_cols, scale = "row", 
               clustering_method = "ward.D2",
               gaps_col = colseps,
               gaps_row =rowseps,
               annotation_row = (annotation_row),
               annotation_col = (annotation_col),
               # annotation_colors = mycolors,
               colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu"))[2:11])(103), 
               show_rownames = show_rownames, show_colnames = show_colnames,
               breaks=breaks, border_color=NA)
  print(p)
  
}
