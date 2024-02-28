#===========================================================
#  1. redad data   
##==========================================================
RNA =read.csv("diff2.csv")
ATAC=read.csv("peak_combine_peak.annotation.csv")

dim(RNA)
dim(ATAC)


combine= merge(RNA,ATAC,
               by.x="SYMBOL",
               by.y="SYMBOL",
               suffixes = c(".RNA",".ATAC") ,
               all.x=F,
               all.y=F)

dim(combine) #895

# merge data ；
write.table(combine,"merge.txt",sep="\t",row.names=FALSE)

#====================
#  2. prepare data
#====================

dat = read.table("merge.txt",sep="\t",header = TRUE)
View(dat)

colnames(dat) 
dim(dat)i

#  x and y axis
x=dat$logFC  # RNAseq
y=dat$Fold   # ATACseq


# filter  subgroup data
## RNA
#RNAdown:
RNAdown<-intersect(which(dat$logFC <= -1),
                   which(dat$P.Value  <= 0.05)) 
RNAdown

#RNAup：
RNAup<-intersect(which(dat$logFC >= 1), 
                 which(dat$P.Value  <=0.05)) 
RNAup

# RNA diff
RNA_diff = intersect(which(abs(dat$logFC) >= 1),
                     which(dat$P.Value <=0.05)) 
RNA_diff

# RNA_nodiff
RNA_nodiff = union(which(abs(dat$logFC) < 1),
                   which(dat$P.Value > 0.05)) 


## ATACseq
# ATACdown
ATACdown<-which(dat$Fold <= -0.6)

#ATACdown<-intersect(which(dat$log2FC <= -0.4))
                   # which(dat$p.value  <= 0.05)) 

# ATACup；
ATACup<-which(dat$Fold >= 0.6)
#ATACup<-intersect(which(dat$Fold >= 0.4))
                 # which(dat$p.value <=0.05)) 

# AATAC_diff
ATAC_diff = which(abs(dat$Fold) >= 0.6 )
View(ATAC_diff)                     


# ATAC_nodiff
ATAC_nodiff = which(abs(dat$Fold) < 0.6)
# ATAC_nodiff = union(which(abs(dat$Fold) < 0.4),
#                    which(dat$p.value > 0.05))  

# samedown
samedown<-intersect(RNAdown,ATACdown)

# sameup
sameup<-intersect(RNAup,ATACup) 

# RNA Down，ATACUP；
diff1<-intersect(RNAdown,ATACup) 

# RNA UP，ATAC Down；
diff2<-intersect(RNAup,ATACdown) 

#差异方向相同的基因数量；
homo=length(union(sameup,samedown)) 
homo #250

#差异方向相反的基因数量；
opp=length(union(diff1,diff2)) 
opp #212


RNAonly=length(intersect(ATAC_nodiff,RNA_diff))

ATAConly=length(intersect(ATAC_diff,RNA_nodiff)) 




#=====================
#   3. Draw figure
#=====================

plot(x,y,pch=16,cex = 0.1,
     xlim=c(-5,3.9),ylim=c(-6.5,6.5),
     cex.lab=1.2,
     #main=main ,
     col="white",
     xlab="RNAseq expression",
     ylab="ATACseq expression")



# ：C63C3B
points(x[c(sameup)],y[c(sameup)],pch=16,cex=0.5,col="#C63C3B")

#：A2D1A6
points(x[c(samedown)],y[c(samedown)],pch=16,cex=0.5,col="#A2D1A6")


#：D8834E
points(x[c(diff1)],y[c(diff1)],pch=16,cex=0.5,col="#4A4C9D")

# ：D8834E
points(x[c(diff2)],y[c(diff2)],pch=16,cex=0.5,col="#D8834E")


# add text
# RHOB :3.896509447,0.701979223


text(x = c(3.896509447), y = c(0.7019792236), 
     labels = c("RHOB"), cex = c( 0.5 , 0.5 ),adj = 1.2)

points(c(3.896509447),
       c(0.7019792236),pch=17,cex=0.5,col="black")

points(x[c(diff1)],y[c(diff1)],pch=16,cex=0.3,col="#D8834E")


###############################################################################
