#' Plot PCA projection and data points from original PCA with all colored in
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
#' @param ellipse Construct confidence region based on groups in info.type, default = T
#' @param conf default = 0.95 
#' @param density plot x-y density plots
#' @param fliph,flipv flip plot hoirzontally or vertically
#' 
# @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
#' 
#' @export
#'


plot_pca_projection_all = function(file, rotated.file, info.name, info.type, info.name2, info.type2, title = "Projection", 
                                   labels = F, PCx="PC1", PCy="PC2", ellipse = F, conf = 0.95, density = F,
                                   fliph = F, flipv = F){
  require(ggplot2)
  pc.scores = read.table(file, header = TRUE, row.names = 1)
  pc.scores.reduced = pc.scores
  pc.scores.reduced$type = info.type[match(rownames(pc.scores.reduced), info.name)]
  
  projected_data = read.table(rotated.file,header = T, row.names = 1)
  projected_data.reduced = projected_data
  projected_data.reduced$type = info.type2[match(rownames(projected_data.reduced), info.name2)]
  
  combined.data = rbind(pc.scores.reduced, projected_data.reduced)
  colnames(combined.data)[ncol(combined.data)] <-"type" 
  combined.data = na.omit(combined.data) #don't plot if no annot
  if (fliph==T){combined.data[,PCx] = combined.data[,PCx]*-1}
  if (flipv==T){combined.data[,PCy] = combined.data[,PCy]*-1}
  
  
  pcx.y <- ggplot(combined.data, aes_string(x=PCx,y=PCy)) +geom_point(size = I(2), aes(color = factor(type))) +
    theme(legend.position="right",plot.title=element_text(size=30),legend.text=element_text(size=22),
          legend.title=element_text(size=20),axis.title=element_text(size=30),legend.background = element_rect(),
          axis.text.x = element_text(margin = margin(b=-2)),axis.text.y = element_text(margin = margin(l=-14)))+
    guides(color=guide_legend(title="Type"))+
    labs(title = title)+
    theme_bw()+
    if(labels==TRUE){geom_text(mapping = aes(label = rownames(combined.data)), check_overlap = TRUE, size = 3)}

  
  if(ellipse==TRUE){
    plot(combined.data[,c(PCx, PCy)], main=title)
    ord = ordiellipse(combined.data[,c(PCx, PCy)],combined.data$type, kind = "sd", conf = conf) 
    
    cov_ellipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
    {
      theta <- (0:npoints) * 2 * pi/npoints
      Circle <- cbind(cos(theta), sin(theta))
      t(center + scale * t(Circle %*% chol(cov)))
    }
    
    df_ell <- data.frame(matrix(ncol = 0, nrow = 0))
    for(g in (droplevels(combined.data$type))){
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(combined.data[combined.data$type==g,],
                                                       cov_ellipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                    ,type=g))
    }
    
    pcx.y2 = pcx.y + geom_path(data=df_ell, aes(x=df_ell[,PCx], y=df_ell[,PCy], colour = type), size=1, linetype=1)
    pcx.y2
  }else{
    pcx.y
  }
  # if(density==TRUE){
  #   
  #   # Marginal density plot of x (top panel) and y (right panel)
  #   xplot <- ggdensity(combined.data, PCx, fill = "type")+ clean_theme()
  #   yplot <- ggdensity(combined.data, PCy, fill = "type")+ rotate()+ clean_theme()
  #   # Arranging the plot
  #   (ggarrange(xplot, NULL, pcx.y, yplot,
  #                   ncol = 2, nrow = 2,  align = "hv",
  #                   widths = c(2, 1), heights = c(1, 2),
  #                   common.legend = TRUE))
  # }
  # else{
  #   print(pcx.y)
  # }
  

}
