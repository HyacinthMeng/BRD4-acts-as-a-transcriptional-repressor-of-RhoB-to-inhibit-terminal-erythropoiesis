#===========================================================
#  1. 读入两个差异分析表（分别是转录组和翻译组）并合并；  
##==========================================================
RNA =read.csv("diff1.csv")
Ribo=read.csv("edgeR_logFC_peak.annotation.csv")

dim(RNA)
dim(Ribo)

##合并两个表格。选择表头为"id"的列进行合并；
combine= merge(RNA,Ribo,
               by.x="SYMBOL",
               by.y="SYMBOL",
               suffixes = c(".RNA",".Ribo") ,
               all.x=F,
               all.y=F)

dim(combine)

#保存合并后的数据；
write.table(combine,"merge.txt",sep="\t",row.names=FALSE)

#====================
#  2. 数据准备 
#====================

dat = read.table("merge.txt",sep="\t",header = TRUE)
View(dat)
# 查看表头;
colnames(dat) 
dim(dat)
#取RNA的差异倍数，将作为X轴;
#取Ribo的差异倍数，将作为Y轴;
x=dat$logFC 
y=dat$log2FC

#计算两个组学差异倍数的相关性，并取2位小数;
cor = round(cor(x,y),2) 
cor  
#准备作为图形的标题;
main = paste("correlation=",cor,sep="")  
main

#数据筛选：
## RNA
#RNAdown:转录组显著下调的行;
RNAdown<-intersect(which(dat$logFC <= -1),
                   which(dat$P.Value  <= 0.05)) 
View(RNAdown)
#write.table(RNAdown,"RNAdown.txt",sep="\t",row.names=FALSE)
#RNAup：转录组显著上调的行;
RNAup<-intersect(which(dat$logFC >= 1), 
                 which(dat$P.Value  <=0.05)) 
RNAup
# RNA显著差异的行；
RNA_diff = intersect(which(abs(dat$logFC) >= 1),
                     which(dat$P.Value <=0.05)) 
RNA_diff
# RNA不显著差异的行；##
RNA_nodiff = union(which(abs(dat$logFC) < 1),
                   which(dat$P.Value > 0.05)) 


## cuttag
#Ribo显著下调的行 ；
Ribdown<-which(dat$log2FC <= -0.4)

#Ribdown<-intersect(which(dat$log2FC <= -0.4))
                   # which(dat$p.value  <= 0.05)) 
#Ribo显著上调的行；
Ribup<-which(dat$log2FC >= 0.4)
#Ribup<-intersect(which(dat$Fold >= 0.4))
                 # which(dat$p.value <=0.05)) 
# Ribo显著差异的行；

Ribo_diff = which(abs(dat$log2FC) >= 0.4 )
                      


# Ribo不显著差异的行；##
Ribo_nodiff = which(abs(dat$log2FC) < 0.4)
# Ribo_nodiff = union(which(abs(dat$log2FC) < 0.4),
#                    which(dat$p.value > 0.05))  

# 两个组学都下调；
samedown<-intersect(RNAdown,Ribdown)

samedn=dat[samedown,]
View(samedn)
write.table(samedn,"samedn.txt",sep="\t",row.names=FALSE)

# 两个组学都上调；
sameup<-intersect(RNAup,Ribup) 

sameup2=dat[sameup,]
View(sameup)
write.table(sameup2,"sameup.txt",sep="\t",row.names=FALSE)

# RNA下调，Ribo上调；
diff1<-intersect(RNAdown,Ribup) 

RNAdn_Cutup=dat[diff1,]

write.table(RNAdn_Cutup,"RNAdn_Cutup.txt",sep="\t",row.names=FALSE)

# RNA上调，Ribo下调；
diff2<-intersect(RNAup,Ribdown) 

RNAup_Cutdn=dat[diff2,]

write.table(RNAup_Cutdn,"RNAup_Cutdn.txt",sep="\t",row.names=FALSE)


#差异方向相同的基因数量；
homo=length(union(sameup,samedown)) 
homo #160

#差异方向相反的基因数量；
opp=length(union(diff1,diff2)) 
opp #106

# 仅仅RNA有差异的基因数量；
RNAonly=length(intersect(Ribo_nodiff,RNA_diff))
# 仅仅Ribo有差异的基因数量；
Riboonly=length(intersect(Ribo_diff,RNA_nodiff)) 

#生成图例所用标签；
RNA_change=paste("Transcription(",RNAonly,")",sep="")
Ribo_change=paste("H3K27me3(",Riboonly,")",sep="")
homo_change=paste("Homodirection(",homo,")",sep="")
opp_change=paste("Opposite(",opp,")",sep="")

#给图例准备颜色
id<-c(RNA_change,Ribo_change,homo_change,opp_change)
co<-c("#D8834E","#A2D1A6","#C63C3B","#4A4C9D")  


# opp
# 右下桔：D8834E；左下绿：A2D1A6 右上红：C63C3B 左上蓝：4A4C9D
#=====================
#   3. 图形绘制
#=====================
# 绘制所有点散点图，注意可以修改X轴和y轴范围；
# 可利用type = "n",记得改坐标轴范围大小
pdf("cuttag_rna1——1104.pdf", width=6, height=5)

plot(x,y,pch=16,cex = 0.1,
     xlim=c(-4,3.0),ylim=c(-2.5,2.4),
     cex.lab=1.2,
     #main=main ,
     col="white",
     xlab="RNAseq expression",
     ylab="cuttag expression")



#利用次级函数绘制不同象限的图，并加入图例；

# RNA diff
points(x[RNA_diff],y[RNA_diff],pch=16,cex=0.3,col="#4A4C9D")

# cuttag Diff
points(x[Ribo_diff],y[Ribo_diff],pch=16,cex=0.3,col="#A2D1A6")

# same Down and UP
## 右下桔：D8834E；左下绿：A2D1A6 右上红：C63C3B 左上蓝：4A4C9D

# 两个组学都下调；
# samedown<-intersect(RNAdown,Ribdown)
# 两个组学都上调；
# sameup<-intersect(RNAup,Ribup) 

# RNA下调，Ribo上调；--左上
# diff1<-intersect(RNAdown,Ribup) 

# RNA上调，Ribo下调；---右下
# diff2<-intersect(RNAup,Ribdown)


# 右上红：C63C3B
points(x[c(sameup)],y[c(sameup)],pch=16,cex=1,col="#C63C3B")

#左下绿：A2D1A6
points(x[c(samedown)],y[c(samedown)],pch=16,cex=1,col="#A2D1A6")


#左侧上蓝：D8834E
points(x[c(diff1)],y[c(diff1)],pch=16,cex=1,col="#4A4C9D")

# 右下桔：D8834E
points(x[c(diff2)],y[c(diff2)],pch=16,cex=1,col="#D8834E")


# 添加文字标签
RHOB :2.101075371,-0.512692746

?text
text(x = c(2.101075371), y = c(-0.512692746), 
     labels = c("RHOB"), cex = c( 0.5 , 0.5 ),adj = -0.5)

points(c(2.101075371),
       c(-0.512692746),pch=17,cex=1.2,col="black")

dev.off()
# text(sameup,gfg_data$y,labels = gfg_data$lab,pos = 4)


#
points(x[c(diff1)],y[c(diff1)],pch=16,cex=0.3,col="#D8834E")

legend("topright",legend=id,col=co,pch=16,bty = "n",cex = 0.8,text.font = 2)

#dev.off()
