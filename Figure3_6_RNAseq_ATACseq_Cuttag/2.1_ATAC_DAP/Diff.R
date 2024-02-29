
################################################################################
##== DAP Analysis
################################################################################

setwd("/home/zju/desktop/ZM/ATAC-ZM/2_2_1212_P0.1")
#BiocManager::install("DiffBind")
BiocManager::install("edgeR")
library(DiffBind)
library(edgeR)

# read file
library(DiffBind)
dbObj <- dba(sampleSheet="sampleID.csv")
View(dbObj)

dbObj <- dba.count(dbObj,bUseSummarizeOverlaps=TRUE))
??dba.plotPCA
dba.plotPCA(dbObj,label=DBA_ID)
plot(dbObj)



##   DAP
# Establishing a contrast 
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)

#  summary of results
dba.show(dbObj, bContrasts=T)
#  overlapping peaks identified by the two different tools (DESeq2 and edgeR)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)

# save result
# EdgeR
out <- as.data.frame(comp1.edgeR)
write.table(out, file="D_edgeR.txt", sep="\t", quote=F, col.names = NA)
# DESeq2
out <- as.data.frame(comp1.deseq)
write.table(out, file="D_deseq2.txt", sep="\t", quote=F, col.names = NA)

# Create bed files for each keeping only significant peaks (p  0.05 or other value) 
# EdgeR
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$p.value < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed, file="Exp_vs_Con_edgeR_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)

# DESeq2
out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$p.value < 0.05), 
                  c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="Exp_Con_deseq2_sig.bed", sep="\t", quote=F, row.names=F, col.names=F)





################################################################################

                             