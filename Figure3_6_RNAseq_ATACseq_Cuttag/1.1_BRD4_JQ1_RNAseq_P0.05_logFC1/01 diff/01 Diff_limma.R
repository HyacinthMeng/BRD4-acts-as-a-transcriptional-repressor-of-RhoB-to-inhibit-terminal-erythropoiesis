###https://bioconductor.org/packages/release/bioc/html/limma.html
#It contains the inroduction of limmma#
####----- BY MY   2021 rewrite @ZJU-----###

###install and library pacakages################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)

###set cutoff value############################################################
logFoldChange=1
adjustP=0.05
#setwd("")


###Load your file---genematrix#################################################
rt=read.table("symbol.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
rt=rt[rowMeans(rt)>0,]  #filter the low expression of gene#
#data=normalizeBetweenArrays(data)
rt=log2(rt+1)
###Diff analysis###############################################################
modType=c(rep("normal",4),rep("tumor",4))#change the group number #
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)

###output result###############################################################
#All result#
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F,row.names=T)

#All diff result#
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & P.Value < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F,row.names=T)

#All UP result#
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & P.Value < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F,row.names=T)

#All Down result#
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & P.Value < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F,row.names=T)

#Outout expression  of diff gene
hmExp=rt[row.names(diffSig),]
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)

################################################################################
