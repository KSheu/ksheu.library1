#' Plot PCA
#' 
#' Plots PCA from scores file (output of PCA_from_file)
#' 
#' @param file File containing scores matrix
#' @param info.name Vector of sample names
#' @param info.type Vector of sample types in the same order
#' @param title Title of the plot
#' @param labels default=T
#' @param PCx,PCy PCs to display
#' 
# @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
#' 
#' @export
#' 
plot_pca = function(file, info.name, info.type, title = "", labels = TRUE, PCx="PC1", PCy="PC2"){  
  #Input: PCA scores file to be ploted
  ##process pca output and adds groupings
  require(ggplot2)
  table <- read.table(file, header = TRUE)
  table$type = info.type[match(table$Score, info.name)]
  
  
  pcx.y <- ggplot(table, aes_string(x=PCx,y=PCy)) +geom_point(size = I(2), aes("Type",color = factor(type))) +
    theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
          legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
          axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
    labs(title = title)+
    theme_bw()+
    if(labels==TRUE){geom_text(data = table, mapping = aes(label = Score), check_overlap = TRUE, size = 3)}
  pcx.y
  
}  