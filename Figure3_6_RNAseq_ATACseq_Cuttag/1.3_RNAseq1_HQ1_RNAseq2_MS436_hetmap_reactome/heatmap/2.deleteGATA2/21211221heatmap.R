
## combine cuttag _ RNAseq information

x=read.table("RNA1_diffExp.txt",header=TRUE,sep="\t",quote = "")
head(x)

y=read.csv("gene_selected.txt")
head(y)

library(tidyverse)
by <- join_by(x$id==y$ID)
combine1 =right_join(x,y,by)
View(combine)
dim(combine)
write.table(combine,'RNA1_selected.xls',sep = '\t', quote = FALSE, row.names = FALSE)




#x=read.table("RNA2_diffExp.txt",header=TRUE,sep="\t",quote = "")
x=read.table("RNA2_allDiffdiffExp.txt",header=TRUE,sep="\t",quote = "")
head(x)

by <- join_by(x$id == y$ID)
combine2 =right_join(x,y,by)
View(combine)
dim(combine)
write.table(combine,'RNA2_selected.xls',sep = '\t', quote = FALSE, row.names = FALSE)




# Draw
library(readxl)
library(pheatmap)

rt1 = read.table("RNA1.txt",header=TRUE,sep="\t",quote = "",row.names = 1)
rt1=as.matrix(rt1)

Type=c(rep("RNA1_Con",4),rep("RNA1_Exp",4))    #Group type
names(Type)=colnames(rt1)
Type=as.data.frame(Type)



## At end of plotting, reset to previous settings:
?pheatmap
p1 = pheatmap(rt1, annotation_col=Type, 
         #annotation_colors = c("navy", "firebrick3"),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cluster_cols =F,
         #cluster_rows  =F,
         scale="row",
         border_color="white",
         fontsize_row=5,
         fontsize_col=5)
p1



rt2 = read.table("RNA2.txt",header=TRUE,sep="\t",quote = "",row.names = 1)
rt2=as.matrix(rt2)
#rt=log2(rt+1)
#rt[rt>15]=15
#View(rt2)
library(pheatmap)

#Type=c(rep("Cuttag_Con",3),rep("Cuttag_Exp",3))
Type2=c(rep("RNA2_Con",3),rep("RNA2_Exp",3))
names(Type2)=colnames(rt2)
Type2=as.data.frame(Type2)

p2 = pheatmap(rt2, annotation_col =Type2, 
            color = colorRampPalette(c("#1D1AEC", "white","#E3CC22"))(50),
            #cluster_cols =F,
            #cluster_rows  =F,
            border_color="white",
            scale="row",
            fontsize_row=5,
            fontsize_col=5)
p2
p1 + p2

pdf(file =  "1heatmap.pdf" ,
    width=7, height=6)

p3=p1 + p2
p3
dev.off()

install.packages("cowplot")
install.packages("patchwork")
install.packages("gridExtra")
install.packages("aplot")
dev.off()






BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
?`ComplexHeatmap-package`





library(devtools)
install_github("jokergoo/ComplexHeatmap")
