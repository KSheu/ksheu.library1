##GSEA-squared (last mod. 07/2019 ksheu)
library(ksheu.library1); library(dplyr); library(ggplot2);library(pheatmap);library(ggpubr);library(Matching);library(RColorBrewer)
source("C:/Users/msheu/Desktop/Katherine/signed-ks-test.R") #from Github ks.test2

#fig1e PLOT----
dfm = read.delim("F:/gsea/small_cell_rnk_GSEA/lungprostate_PCA_V1_c5.GseaPreranked.1515644251568/gsea_report_for_na_neg_1515644251568.xls")
dfm2 = read.delim("F:/gsea/small_cell_rnk_GSEA/lungprostate_PCA_V1_c5.GseaPreranked.1515644251568/gsea_report_for_na_pos_1515644251568.xls")
dfm = rbind(dfm, dfm2)
dfm = dfm[order(dfm$NES),]
dfm$rnk = seq(1:nrow(dfm))

#figure out important keywords
words = data.frame(table(unlist(strsplit(tolower(dfm$NAME), "_"))))
words = words[order(words$Freq, decreasing = T),]
words$Var1 = toupper(words$Var1)
str(words)
words.upper = words[which(words$Freq>5&words$Freq<500),]

set.seed(1);kspvals = data.frame(word = as.character(), pval = as.numeric(), ES = as.numeric())
for (i in seq(1:nrow(words.upper))){
  tryCatch({
    print(i)
    print(words.upper$Var1[i])
    test<-as.numeric(dfm$rnk[(grepl(words.upper$Var1[i], dfm$NAME))] ) 
    background<-as.numeric(dfm$rnk[!(grepl(words.upper$Var1[i], dfm$NAME))] ) 
    # pval = data.frame(word = words.upper$Var1[i], pval= ks.boot(test, background, n=100)$ks.boot$p.value)
    pval = data.frame(word = words.upper$Var1[i], pval= ks.test.2(test, background)$p, ES= ks.test.2(test, background)$ES)
    kspvals = rbind(kspvals, pval)
  },error=function(e){})
}

kspvals$freq = words.upper$Freq[match(kspvals$word, words.upper$Var1)]
kspvals = kspvals[order(kspvals$pval),]
# kspvals$rnk = seq(1:nrow(kspvals))
# kspvals$type = ifelse(grepl("NEURON|NEUROTRANSMITTER|SYNAP|VOLTAGE|CEREBRAL|CORTEX", kspvals$word), "neuro",
#                       ifelse(grepl("CELL_CYCLE|MITOSIS|MITOTIC|DNA_REP|CHROMOSOME_SEG|SPINDLE", kspvals$word), "cycle",
#                              ifelse(grepl("INFLAM|IMMUNE|IMMUNITY|INTERLUEKIN|LEUKOCYTE", kspvals$word), "immune",
#                                     ifelse(grepl("ADHE", kspvals$word), "adhesion",
#                                            ifelse(grepl("SPLIC", kspvals$word),"splicing","x_other")))))
# ggplot(kspvals, aes(type, rnk))+geom_point(aes(color = type), position = "jitter")+
#   theme(text = element_text(size=20)) +theme_classic()
# write.table(kspvals, "D:/gsea/small_cell_rnk_GSEA/lungprostate_PCA_V1_c5.GseaPreranked.1515644251568/word_freq_ks.boot.txt", quote =F, sep = "\t", row.names = F)
write.table(kspvals, "F:/gsea/small_cell_rnk_GSEA/lungprostate_PCA_V1_c5.GseaPreranked.1515644251568/word_freq_ks.test.2.txt", quote =F, sep = "\t", row.names = F)


dfm$type = ifelse(grepl("NEURO|NEUROTRANSMITTER|SYNAP|VOLTAGE|AXON|CEREBRAL|CORTEX", dfm$NAME), "neuro",
                  ifelse(grepl("CELL_CYCLE|MITOTIC|DNA_REPLICATION|CHROMOSOME_SEGREGATION|SPINDLE|CELL_DIVISION", dfm$NAME), "cycle",
                         ifelse(grepl("INFLAM|IMMUNE|IMMUNITY|INTERLUEKIN|LEUKOCYTE", dfm$NAME), "immune",
                                ifelse(grepl("ADHESION|ADHERENS", dfm$NAME), "adhesion",
                                       ifelse(grepl("SPLIC", dfm$NAME),"splicing","x_other")))))
dfm$logFDR = log10(dfm$FDR.q.val)
ggplot(dfm, aes(type, rnk))+geom_point(aes(color = type), position = "jitter")+
  theme_classic(base_size = 11)+ theme(text = element_text(size=20), axis.text.y = element_text(angle = -30, vjust = 1, hjust=1.2),legend.position = "none") +scale_x_discrete(limits=c("adhesion", "immune", "splicing", "cycle", "neuro"))+
  geom_abline(slope = 0, intercept =  4591,size = 1, linetype="dotted")+
  geom_abline(slope = 0, intercept =  488,size = 1,linetype="dotted")+ylab("rank")

plot(abs(dfm$NES), dfm$FDR.q.val)
dfm$sig = ifelse((dfm$FDR.q.val<0.05),1,0)
table(dfm$type, dfm$sig)


frame = data.frame(pval = as.numeric())
test<-as.numeric(dfm$rnk[(dfm$type=="neuro")] )
background<-as.numeric(dfm$rnk[!(dfm$type=="neuro")] )
frame = rbind(frame, -log10(ks.test.2(test, background)$p))
test<-as.numeric(dfm$rnk[(dfm$type=="cycle")] )
background<-as.numeric(dfm$rnk[!(dfm$type=="cycle")])
frame = rbind(frame, -log10(ks.test.2(test, background)$p))
test<-as.numeric(dfm$rnk[(dfm$type=="splicing")] )
background<-as.numeric(dfm$rnk[!(dfm$type=="splicing")])
frame = rbind(frame, -log10(ks.test.2(test, background)$p))
test<-as.numeric(dfm$rnk[(dfm$type=="immune")] )
background<-as.numeric(dfm$rnk[!(dfm$type=="immune")] )
frame = rbind(frame, -log10(ks.test.2(test, background)$p))
test<-as.numeric(dfm$rnk[(dfm$type=="adhesion")] )
background<-as.numeric(dfm$rnk[!(dfm$type=="adhesion")] )
frame = rbind(frame, -log10(ks.test.2(test, background)$p))
