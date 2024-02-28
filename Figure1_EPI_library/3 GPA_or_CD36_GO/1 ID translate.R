###https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
   ##It contains the inroduction of org.Hs.eg.db#
   
####----- BY MY   2021 rewrite @ZJU-----###

###   install and library pacakages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")

library("org.Hs.eg.db")
getwd()
setwd("")

# translate ID
rt=read.table("GPA_CD36_Diff_0.5.txt",sep="\t",check.names=F,header=T)
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="CD36_GPA_0.5_ID.txt",sep="\t",quote=F,row.names=F)

################################################################################
