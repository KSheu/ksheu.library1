#' Plot PCA with save option
#'
#' @param scores.file PCA score file (scores file, output of PCA_from_file)
#' @param info.name vector of sample names
#' @param info.type vector of sample types in same order as names
#' @param PCx,PCy PCs to display
#' @param aes_fill column header to color samples by
#' @param savename prefix to save file by
#' @param savePlot default=F
#' @param width,height size of the plot
#' @param title title of the plot
#' @param label column name to label samples by on PCA plot, default = F
#' @param savetype ".pdf" or ".png" or ".tiff", etc.
#' @param colors ie. colpalette <- c("Cluster 1"="#4F81BD","Cluster 2"="#C0504D","Cluster 3"="#8064A2","Cluster 4"="#9BBB59")
#' @param factor.levels ie. levels <- c("Undifferentiated", "Neural crest-like", "Transitory", "Melanocytic"))
#' 
#' @export


# a.scores=PCA score data.frame with column annotations
# aes_fill= column header to color samples by
# savename= prefix to save file by
# label= column name to label samples by on PCA plot, false if no label
# Function

plot_pca_and_save=function(scores.file, info.name, info.type, PCx="PC1",PCy="PC2", factor.levels = levels(info.type),aes_fill=NULL,savename,savePlot=F,width=3,height=2,title = "",label=F,savetype=".pdf",w=6,h=3,aspect.ratio=0.8,legend=F,legendname="default",colors=colpalette[1:4],psize=2.5,plotLabelname="none") {
  #plot graph
  require(ggplot2)
  savename <- paste(savename,aes_fill,sep="_")
  offset=1.1
  xmax=max(a.scores[PCx])*offset
  xmin=min(a.scores[PCx])*offset
  ymax=max(a.scores[PCy])*offset
  ymin=min(a.scores[PCy])*offset
  
  if(legendname=="default") legendname=aes_fill
  
  a.scores = read.table(scores.file, header = T)
  a.scores$type = info.type[match(a.scores$Score, info.name)]
  a.scores$type<-factor(a.scores$type, levels=factor.levels)
  
  p=ggplot(a.scores,aes_string(x=PCx,y=PCy))+geom_point(size=psize,pch=21,colour="black",aes_string(fill=aes_fill))+theme_bw()+labs(title = title)+theme(aspect.ratio=aspect.ratio,panel.border=element_rect(colour="black",size=1),legend.position="none",axis.title = element_text(size=10),axis.text = element_text(size=10))
  if(is.null(colpalette)) p=p+scale_fill_manual(legendname,values=colors) else p=p+scale_fill_manual(legendname,values=colors)
  if(label!=F) p=p+geom_text(aes_string(label=label),size=1,hjust=-0.3) 
  if(legend==T) p=p+theme(legend.key =element_blank(), legend.position="right", legend.text=element_text(size=11),panel.border=element_rect(colour="black",size=1),legend.key.size=unit(4,"mm"))+ guides(fill = guide_legend(override.aes = list(size=psize*1.5)))
  
  if(plotLabelname!="none") p=p+xlim(xmin,xmax)+ylim(ymin,ymax)+geom_text(data=data.frame("X"=xmin,"Y"=ymax,"lab"=plotLabelname),size=3,aes_string("X","Y",label="lab"),fontface=2,hjust=0,vjust=1)
  if(plotLabelname=="none") p <- p+xlim(xmin,xmax)+ylim(ymin,ymax)
  
  if(savePlot==T) ggsave(paste0(savename,"_",PCx,"_vs_",PCy,savetype),dpi=300,plot=p,width=w,height=h) 
  p
}