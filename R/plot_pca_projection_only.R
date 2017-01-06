#' Plot PCA projection only
#' 
#' Plots projection from rotated.scores file (output of intersect_do_PCA_and_project_second_dataset)
#' 
#' @param rotated.file2 File containing scores matrix
#' @param info.name Vector of sample names
#' @param info.type Vector of sample types in the same order
#' @param title Title of the plot
#' @param labels default=T
#' @param PCx,PCy PCs to display
#' 
# @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
#' @export
#'
#'
plot_pca_projection_only= function(rotated.file2, info.name, info.type, title = "Projection", labels = TRUE, PCx="PC1", PCy="PC2"){
  require(ggplot2)
  projected_data = read.delim(rotated.file2)
  projected_data$type = info.type[match(projected_data$Sample, info.name)]

  pcx.y <- ggplot(projected_data, aes_string(x=PCx,y=PCy)) +geom_point(size = I(2), aes(color = factor(type))) +
    theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
          legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
          axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
    labs(title = title)+
    theme_bw()+
    if(labels==TRUE){geom_text(data = projected_data, mapping = aes(label = Sample), check_overlap = TRUE, size = 3)}
  pcx.y
  
}