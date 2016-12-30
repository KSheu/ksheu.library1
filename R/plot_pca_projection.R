#' Plot PCA projection and data points from original PCA
#' 
#' Plots projection from rotated.scores file (output of intersect_do_PCA_and_project_second_dataset) on top of original data
#' 
#' @param file scores file of original data
#' @param rotated.file File containing scores matrix of projected data
#' @param info.name Vector of sample names in original PCA 
#' @param info.type Vector of sample types in original PCA in the same order as names
#' @param info.name2 Vector of sample names of projected data
#' @param info.type2 Vector of sample types of projected data in the same order
#' @param title Title of the plot
#' @param labels default=T
#' @param PCx,PCy PCs to display
#' 
# @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
#' 
#' @export
#'

plot_pca_projection = function(file, rotated.file, info.name, info.type, info.name2, info.type2, title = "Projection", labels = F, PCx="PC1", PCy="PC2"){
  require(ggplot2)
  pc.scores = read.table(file, header = TRUE, row.names = 1)
  pc.scores.reduced = pc.scores[, 1:5]
  pc.scores.reduced$type = info.type[match(rownames(pc.scores.reduced), info.name)]
  
  projected_data = read.table(rotated.file,header = T, row.names = 1)
  projected_data.reduced = projected_data[,1:5]
  projected_data.reduced$type = info.type2[match(rownames(projected_data.reduced), info.name2)]
  
  pcx.y <- ggplot(pc.scores.reduced, aes_string(x=PCx,y=PCy)) +geom_point(size = I(2), aes("Type",color = factor(type))) +
    theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
          legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
          axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
    labs(title = title)+
    theme_bw()+
    if(labels==TRUE){geom_text(mapping = aes(label = rownames(pc.scores.reduced)), check_overlap = TRUE, size = 3)}
  
  pcx.y <- pcx.y + ggplot(projected_data.reduced, aes_string(x=PCx,y=PCy)) +geom_point(size = I(2), aes("Type",color = factor(type))) +
    theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
          legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
          axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
    labs(title = title)+
    theme_bw()+
    if(labels==TRUE){geom_text(mapping = aes(label = rownames(projected_data.reduced)), check_overlap = TRUE, size = 3)}
  
  pcx.y
  
}