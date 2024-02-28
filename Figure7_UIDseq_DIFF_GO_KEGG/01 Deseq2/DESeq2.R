###########################################################################
# BY MY 20210809 
##########################################################################


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")



library(DESeq2)
library(limma)


# set cutoff
foldChange=1
padj=0.05



rt=read.table("UIDseq_symbol_L1_R1.txt",sep="\t",header=T,check.names=F)

# get tidy data
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

# avereps
data=avereps(data)
dim(data)
# 39623     8

# filter
data=data[rowMeans(data)>1,] 
# round
data=round(data,0)

# Group infomation
group_list=c(rep("Luc",4),rep("RHOB",4))
group_list = factor(group_list)


colData <- data.frame(row.names=colnames(data), 
                      group=group_list)


# Deseq
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = colData,
                              design = ~group)

dds2 <- DESeq(dds)
rld <- rlogTransformation(dds2)  
data_new=assay(rld)

#write.csv(data_new,file="Deseq2_gene_Norlization_log2.csv",quote=FALSE)

resultsNames(dds2)
res <-  results(dds2, contrast=c("group","RHOB","Luc"))   # Luc: control

resOrdered <- res[order(res$padj),]

resOrdered=as.data.frame(resOrdered)
resOrdered=na.omit(resOrdered)

# All result
write.csv(resOrdered,file="Deseq2_UID_Diff.csv",quote = FALSE)


# diffSig
diffSig = resOrdered[(resOrdered$padj < padj & (resOrdered$log2FoldChang>foldChange | resOrdered$log2FoldChange<(-foldChange))),]
write.table( diffSig, file="Deseq2_UID_diffSig.xls",sep="\t",quote=F,row.names=T)

# UP
diffUp = resOrdered[(resOrdered$padj < padj & (resOrdered$log2FoldChang>foldChange)),]
write.table( diffUp, file="Deseq2_UID_up.xls",sep="\t",quote=F,row.names=T)

# Down
diffDown = resOrdered[(resOrdered$padj < padj & (resOrdered$log2FoldChange<(-foldChange))),]
write.table( diffDown, file="Deseq2_UID_down.xls",sep="\t",quote=F,row.names=T)

diffExp=rbind(id=colnames(data_new),data_new[rownames(diffSig),])
write.table(diffExp,file="Deseq2_UID_diffmRNAExp.txt",sep="\t",quote=F,col.names=T) 


###-----------------------------------------------------------------------------
