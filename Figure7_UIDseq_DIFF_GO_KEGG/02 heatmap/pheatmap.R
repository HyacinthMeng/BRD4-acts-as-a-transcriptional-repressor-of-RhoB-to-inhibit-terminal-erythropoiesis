###https://cran.r-project.org/web/packages/pheatmap/index.html
  #It contains the inroduction of pheatmap#
####----- BY MY   2021 rewrite @ZJU-----###

###install and library pacakages################################################
#install.packages("pheatmap")



library(pheatmap)
#setwd("") 

##load file 
# rt=read.table("../Deseq2_UID_diffmRNAExp.txt",sep="\t",header=T,check.names=F,row.names = 1)

rt=read.table("selected_exp.txt",sep="\t",header=T,check.names=F,row.names = 1)

#rt=log2(rt+1)
View(rt)

Type=c(rep("L1",4),rep("R1",4))    #Group type
names(Type)=colnames(rt)
Type=as.data.frame(Type)



pdf(file="heatmap_UID.pdf",width = 7.0,height = 10)
#?pheatmap
pheatmap(rt, annotation=Type, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cluster_cols =F,
         scale="row",
         show_rownames=T,
         border_color=NA,
         fontsize_row=3,
         fontsize_col=5)

dev.off()

################################################################################
