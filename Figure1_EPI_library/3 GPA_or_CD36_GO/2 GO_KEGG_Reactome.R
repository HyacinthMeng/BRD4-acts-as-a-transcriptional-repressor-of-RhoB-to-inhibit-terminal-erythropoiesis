### Introduction ：https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
####----- BY MY   2021 rewrite @ZJU-----###


###install and library pacakages################################################
#install.packages(c("colorspace"，"stringi","DOSE","clusterProfiler","pathview")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ReactomePA")
# BiocManager::install("DO.db",force = TRUE) #注意安装相应依赖包
# BiocManager::install("reactome.db")

###load pacakages########################################################################
library(DO.db)
library(ReactomePA)
library("ggplot2")
library("clusterProfiler")
library("org.Hs.eg.db")





rt=read.table("CD36_GPA_0.5_ID.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]

# pvalueFilter=0.05          
# qvalueFilter=0.05           


geneFC=rt$logFC
gene=rt$entrezID
names(geneFC)=gene

### GO ------------------------------------------------------------------------
GO <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =1, 
               qvalueCutoff = 1,
               ont="all",
               readable =T)
write.table(GO,file="GPA_GO.txt",sep="\t",quote=F,row.names = F)




# barplot
pdf(file="GPA_GO_bar.pdf",width = 10,height = 10)
barplot(GO, drop = TRUE, 
        showCategory =10,
        split="ONTOLOGY")+
        facet_grid(ONTOLOGY~., scale='free')
dev.off()

# dotplot
pdf(file="GPA_GO_bubble.pdf",width = 10,height = 10)
dotplot(GO,showCategory = 10,split="ONTOLOGY",
               orderBy = "GeneRatio") + 
              facet_grid(ONTOLOGY~., scale='free')
dev.off()





### KEGG  ----------------------------------------------------------------------
KEGG <- enrichKEGG(gene = gene, organism = "hsa", 
                 pvalueCutoff =1, qvalueCutoff =1)
write.table(KEGG,file="CD36_GPA_KEGG.txt",sep="\t",quote=F,row.names = F)

# barplot
pdf(file="CD36_GPA_KEGG_barplot.pdf",width = 10,height = 7)
barplot(KEGG, drop = TRUE, showCategory = 20)
dev.off()

# dotplot
pdf(file="CD36_GPA——KEGG_bubble.pdf",width = 10,height = 7)
dotplot(KEGG, showCategory = 20, orderBy = "GeneRatio")
dev.off()




### reactome Pathway  ---------------------------------------------------------
reactome <- enrichPathway(gene = gene,
                          pvalueCutoff =0.05, 
                          readable =T)
write.table(reactome,file="GPA_GO_reactome.txt",sep="\t",quote=F,row.names = F)

# dotplot
pdf(file="GPA_reactome_GO_bubble.pdf",width = 10,height = 15)
dotplot(reactome, showCategory=22)
dev.off()


# barplot
pdf(file="GPA_reactome_GO_bar.pdf",width = 10,height = 15)
barplot(reactome, showCategory=22)
dev.off()
###############################################################################

