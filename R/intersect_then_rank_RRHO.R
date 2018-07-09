#RANK-RANK HYPERGEOMETRIC FUNCTIONS
##################################################
#BETTER VERSION THIS IS HOW IT WAS DONE IN THE PAPER
#INPUT: Takes two metric files (they are just names in column 1 and numbers in column 2, which will be used to do the ranking)
#This version takes the common names(i.e. genes) between the two files, and then performs ranking
#Also if gene name appears more than once throws out everything but first instance
#OUTPUT: Spits out a file in the format accepted by the "simple" version of the RRHO algorithm located at http://systems.crump.ucla.edu/rankrank/rankranksimple.php


create_RRHO_file_intersect_then_rank<-function(metric_file1,metric_file2,outfile.prefix=paste0("rankrank.",metric_file1,"_vs_",metric_file2),metric.col.file1=2,metric.col.file2=2,rank1="desc",rank2="desc") {
  
  rankcol1= metric.col.file1-1
  rankcol2= metric.col.file2-1
  
  metric1 <- read.table(metric_file1, row.names=NULL,check.names=F,header=F,stringsAsFactors=F,sep="\t")
  metric2 <- read.table(metric_file2, row.names=NULL,check.names=F,header=F,stringsAsFactors=F,sep="\t")
  metric1[,1] <- toupper(metric1[,1])
  metric2[,1] <- toupper(metric2[,1])
  metric1 <- metric1[!duplicated(metric1[,1]),]
  metric2 <- metric2[!duplicated(metric2[,1]),]
  rownames(metric1)<-make.names(metric1[,1])
  rownames(metric2)<-make.names(metric2[,1])
  metric1 <- metric1[,-1,drop=FALSE]
  metric2 <- metric2[,-1,drop=FALSE]
  
  metric1 <- metric1[,rankcol1,drop=FALSE]
  metric2 <- metric2[,rankcol2,drop=FALSE]
  rownames(metric1) <-toupper(rownames(metric1))
  rownames(metric2)<-toupper(rownames(metric2))
  
  metric1.names = rownames(metric1)
  metric2.names = rownames(metric2)
  
  common_names<-intersect(metric1.names,metric2.names)
  print(paste0("File1 genes=",length(metric1.names)," File2 genes=",length(metric2.names)," Intersect genes=",length(common_names)))
  metric1.positions<- match(common_names,metric1.names)
  metric2.positions<- match(common_names,metric2.names)
  
  #matched 
  mm1 <- metric1[metric1.positions,,drop=FALSE]
  mm2 <- metric2[metric2.positions,,drop=FALSE] 
  
  #matched and ordered in descending order(i.e row 1 is highest value)
  if(rank1=="desc"){
    mmo1 <- mm1[order(-mm1[,1]),,drop=FALSE]
  } else {
    mmo1 <- mm1[order(mm1[,1]),,drop=FALSE]
  }
  if(rank2=="desc"){
    mmo2 <- mm2[order(-mm2[,1]),,drop=FALSE]
  } else {
    mmo2 <- mm2[order(mm2[,1]),,drop=FALSE]
  }
  
  
  #matched,ordered,and ranked
  mmr1 <- transform(mmo1, rank = 1:nrow(mmo1))
  mmr2 <- transform(mmo2, rank = 1:nrow(mmo2))
  
  #now use ranks from set 1 as reference
  ord<- match(rownames(mmr1),rownames(mmr2)) 
  mmr2<-mmr2[ord,,drop=FALSE]
  
  
  both <- data.frame(rownames(mmr1),rownames(mmr2),mmr1[,2],mmr2[,2],mmr1[,1],mmr2[,1]) 
  colnames(both) <- c("Unigene","Gene_Symbol",paste0(metric_file1,".rank"),paste0(metric_file2,".rank"),paste0(metric_file1,".metric"),paste0(metric_file2,".metric"))
  write.table(both, file = paste0(outfile.prefix, ".txt"), col.name=T,row.names=F, sep="\t", quote=F)
}
