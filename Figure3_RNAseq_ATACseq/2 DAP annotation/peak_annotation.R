# 表观 peak 注释
# 针对单个文件


# installl packages
#BiocManager::install("ChIPseeker")
#BiocManager::install("enrichplot")
#BiocManager::install("ChIPpeakAnno")
#BiocManager::install("EnsDb.Hsapiens.v79")
#BiocManager::install("ReactomePA")
#BiocManager::install("clusterProfiler")
#install.packages("devtools")
#install.packages("ggupset")
#install.packages("ggimage")

#devtools::install_github("GuangchuangYu/ggtree")
#devtools::install_github("YuLab-SMU/ChIPseeker")
 
# load packages
library("GenomicFeatures")
library("enrichplot")
library("ChIPpeakAnno")
library(ChIPseeker)
library(ggupset)
library(ggplot2)


#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

# 读取文件
# 读取文件为一个list 对象
#peaklist = c()
#for (id in list.files(pattern = ".bed")){
#  peak <- readPeakFile(id)
#  peaklist[id]= list(peak)
#}
#class(peaklist)

  
##== 1 读取数据
peak <- readPeakFile("deseq2.bed")


##== 2 创建物种 txdb
##这个可以读gtf或者gff3文件
txdb <-makeTxDbFromGFF(file="./Homo_sapiens.GRCh38.105.gtf",format="gtf")

 #annoData <- toGRanges(txdb, feature="gene")
#annoData[1:2]

##== 3 #查看peak在全基因组的位置: 峰值覆盖图
covplot(peak)
covplot = covplot(peak)
ggsave("1.covplot.png",plot=covplot,height= 10,width=14)
ggsave("1.covplot.pdf",plot=covplot,height= 10,width=14)

##== 4 peaks结合TSS 区域的情况
#与 TSS 区域结合的 ChIP 峰谱
#首先，为了计算与TSS区域结合的ChIP峰的轮廓，我们应该准备TSS区域，
#即TSS位点的侧翼序列。然后对齐映射到这些区域的峰，并生成tagMatrix。

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

# 与 TSS 区域的热图
tagHeatmap(tagMatrix, xlim =c(-3000, 3000), color="red")
Heatmap = tagHeatmap(tagMatrix, xlim =c(-3000, 3000), color="red")
ggsave("2.Heatmap.png",plot=Heatmap,height= 10,width=14)
ggsave("2.Heatmap.pdf",plot=Heatmap,height= 10,width=14)
# 或者以下代码一步即可
# peakHeatmap(peak, TxDb=txdb, upstream=3000, downstream=3000)

#与 TSS 区域结合的 ChIP 峰的平均分布：Average Profile of ChIP peaks binding to TSS region
TSSAvg_resample = plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            conf=0.95,resample = 100,
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
ggsave("3.TSSAvg_resample.png",plot=TSSAvg_resample,height= 10,width=14)
ggsave("3.TSSAvg_resample.pdf",plot=TSSAvg_resample,height= 10,width=14)

#或者使用以下命令将生成与上所示相同的图
# plotAvgProf2(peak, TxDb=txdb, upstream=3000, downstream=3000,
 #            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")




##== 4：chip 与其它区域关系： Profile of ChIP peaks binding to different regions

# Profile of ChIP peaks binding to body regions：可选择多种
body = plotPeakProf2(peak = peak, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb, weightCol = "V5",ignore_strand = F)
ggsave("4.body.png",plot=body,height= 10,width=14)
ggsave("4.body.pdf",plot=body,height= 10,width=14)
# get the profile ChIP peaks binding to gene body regions with no 
# flank extension or flank extension decided by actual length.

## The first method using getBioRegion(), getTagMatrix() and plotPeakProf() to plot in three steps.
genebody <- getBioRegion(TxDb = txdb,
                         by = "gene")
                         # type = "body")

matrix_no_flankextension <- getTagMatrix(peak,windows = genebody)
                                         #nbin = 800)

# 1. 创建画布
png(filename = "PeakProf.png",units = "px",bg = "white",res = 72)             

plotPeakProf(matrix_no_flankextension,conf = 0.95)
dev.off()

pdf(filename = "PeakProf.pdf", width = 5,height = 5)
plotPeakProf(matrix_no_flankextension,conf = 0.95)
dev.off()
# ......




##== 5 peaks注释
# 注意不同txdb （UCSC  ensemble等）
??annotatePeak

peakAnno <- ChIPseeker::annotatePeak(
  peak,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = "org.Hs.eg.db")

write.table(
  as.data.frame(peakAnno),
  "peak.annotation.csv",
  sep=",",
  row.names = F)
  #quote = "")

# 可视化基因组注释
png(filename = "AnnoPie.png", width = 480,height = 480,units = "px",bg = "white",res = 72)             
plotAnnoPie(peakAnno)
dev.off()

pdf("annoPie.pdf",width = 7,height = 5)              # 分辨率
plotAnnoPie(peakAnno)
dev.off()


bar = plotAnnoBar(peakAnno)
ggsave("6.bar_Anno.png",plot=bar,height= 10,width=14)
ggsave("6.bar_Anno.pdf",plot=bar,height= 10,width=14)

##== 一些注释会重叠
venn = vennpie(peakAnno)


##==  upsetplot（重叠）
upsetplot(peakAnno)


##== 组合图
png(filename = "upsetplot.png", width = 480,height = 480,units = "px",bg = "white",res = 72)             
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

pdf("upsetplot.pdf",width = 9,height = 5)              # 分辨率
upsetplot(peakAnno, vennpie=TRUE)
dev.off()


##== Visualize distribution of TF-binding loci relative to TSS

DistToTSS = plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

ggsave("8.DistToTSS_Anno.png",plot=DistToTSS,height= 10,width=14)
ggsave("8.DistToTSS_Anno.pdf",plot=DistToTSS,height= 10,width=14)

##== 6 Functional enrichment analysis功能富集分析
#BiocManager::install("topGO")
#BiocManager::install("ReactomePA")

library(DOSE)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(ReactomePA)
##没有相关的db包，则需要hub后进行保存进行调用，eg：拟南芥，但是拟南芥是有orgdb的包的
#require(AnnotationHub)
#hub <- AnnotationHub()
#query(hub, "Arabidopsis thaliana")
#Ara <- hub[["AH95951"]]
#keytypes(Ara)
#columns(Ara)
#length(keys(Ara))[1]
#head(keys(Ara))
#saveDb("Ara", file）
#as.GRanges(peakAnno) %>% head(3)



??enrichPathway

MY= as.data.frame(peakAnno)$ENTREZID
View(MY)
#reactome
pathway1 <- enrichPathway(as.data.frame(peakAnno)$ENTREZID,
                          pvalueCutoff = 2) #函数参数是ENTREZID
View(pathway1)
write.table(pathway1,'Down_reactome_pathway.xls',sep = '\t', quote = FALSE, row.names = FALSE)

#####
png(filename = "d_reactome.png", width = 480,height = 550,units = "px",bg = "white",res = 72)             

dotplot(pathway1,x = "GeneRatio",color = "p.adjust",showCategory = 10,size = NULL,
        split = NULL,font.size = 12,title = "",orderBy = "x",label_format = 30,
)
dev.off()


####
pdf("d-reactome.pdf",width = 9,height = 10)              # 分辨率
dotplot(pathway1,x = "GeneRatio",color = "p.adjust",showCategory = 10,size = NULL,
        split = NULL,font.size = 12,title = "",orderBy = "x",label_format = 30,
)
dev.off()


dotplot(pathway1)
write.table(pathway1,'d_reactome_pathway.xls',sep = '\t', quote = FALSE, row.names = FALSE)



#gene <- seq2gene(peak, tssRegion = c(-1000, 1000), 
#                 flankDistance = 3000, TxDb=txdb)
#pathway2 <- enrichPathway(gene)
#head(pathway2, 2)
#dotplot(pathway2)
       
##读取已经加载的org.db包
#library(AnnotationDbi)
#gly=loadDb("D:/R-4.1.3/library/AnnotationDbi/extdata/tair.orgdb")
#keytypes(tair)
       
##基因list文件读取
gene_list <- read.table('peak.annotation.csv',header=T,sep=",")
head(gene_list)
dim(gene_list)

 ##转换geneid,有的物种直接用这个就能转换了，有的需要提前转换好
??bitr
gene.df<-bitr(gene_list$geneId,fromType='ENSEMBL',
                     toType=c("ENTREZID", "SYMBOL", "REFSEQ"),
                     OrgDb='org.Hs.eg.db')
head(gene.df)
dim(gene.df)
??enrichGO
GO = enrichGO(gene_list$ENTREZID, OrgDb = "org.Hs.eg.db", 
               keyType="ENTREZID", pvalueCutoff=1,
               qvalueCutoff=1,ont='all')
head(GO)
write.table(GO,'GO.xls',sep = '\t', quote = FALSE, row.names = FALSE)
       
##Dotplot visualization
plot1=dotplot(GO, split="ONTOLOGY",showCategory = 20,label_format=100) + facet_grid(ONTOLOGY~., scale="free")
       
ggsave("9.d_GO_dotplot.png",plot=plot1,height= 16,width=14)
ggsave("9.d_GO_dotplot.pdf",plot=plot1,height= 16,width=14)

plot2= barplot(GO, split="ONTOLOGY",showCategory = 20,label_format=100) + facet_grid(ONTOLOGY~., scale="free")
ggsave("10.d_GO_barplot.png",plot=plot2,height= 16,width=14)
ggsave("10.d_GO_barplot.pdf",plot=plot2,height= 16,width=14)
       

#  KEGG analysis
??enrichKEGG
#kegg富集
R.utils::setOption("clusterProfiler.download.method",'auto')

KEGG = enrichKEGG(gene_list$ENTREZID, organism = "hsa",
                      pvalueCutoff=1,keyType = 'kegg',
                  )
head(KEGG)
write.table(KEGG,'d_KEGG.xls',sep = '\t', 
            quote = FALSE, row.names = FALSE)

plot3= dotplot(KEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")
                      
ggsave("11.d_KEGG_plot.png",plot=plot3,height= 10,width=14)
ggsave("11.d_KEGG_plot.pdf",plot=plot3,height= 10,width=14)                      

plot4= barplot(KEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")
ggsave("12.d-KEGG_bar.png",plot=plot4,height= 10,width=14)
ggsave("12.d-KEGG_bar.png",plot=plot4,height= 10,width=14)
                      








################################################################################
