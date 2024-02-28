#-----        http://mfuzz.sysbiolab.eu/      --------#
#BY MY  2021-11-20
###------Step1 install packages------###

# BiocManager::install("limma")

library(Mfuzz)
library(limma)


###------Step2 get mean matrix ------###

rt <-read.table("CB2_Diff.csv",header =T,sep=",",check.names = FALSE,row.names = 1)
rt2=t(rt)
rt3 = as.data.frame(rt2)




#------mean by group------#
out=aggregate(rt3[,2:ncol(rt3)],by=list(rt3$Type),mean) 
exp=t(out)

colnames(exp) =exp[1,]
exp_f=exp[2:nrow(exp),]
View(exp_f)        



###------Step3 ------###
exp_f <-read.table("mfuzz_exp_log2.txt",header =T,sep="\t",check.names = FALSE,row.names = 1)
View(exp_f)
exp_f<-as.matrix(exp_f)


eset <- new("ExpressionSet",exprs = exp_f) 
eset <-filter.NA(eset,thres =0.25)
eset <- filter.std(eset,min.std=0)

eset <- standardise(eset)

###------Step4 mfuster------###

set.seed(123)
c <- 4    # set mfuster num
m <- mestimate(eset) 
mf <- mfuzz(eset, c = c, m = m)

###------Step5 view result------###
str(mf)
summary(mf)


###------Step6 visuallize------###

library(RColorBrewer)

color.2 <- colorRampPalette(rev(c("#756bb1","#bcbddc", "#efedf5")))(50)


mfuzz.plot(eset,mf,mfrow=c(1,4),
           new.window= FALSE,
           time.labels= colnames(eset) ,
           colo = color.2)

pdf("mfuzz1.pdf")
mfuzz.plot(eset,mf,mfrow=c(3,3),
           new.window= FALSE,
           time.labels=colnames(exp_f)  
           )
dev.off()





library(RColorBrewer)

color <- colorRampPalette(rev(c("#E02401","#0F52BA", "#3E7C17")))(1000)
mfuzz.plot2(eset,mf,
            mfrow=c(2,3),
            new.window= FALSE,
            #time.labels=colnames(exp_f),
            #centre=TRUE ,
            colo =color )



# More fancy choice of colors
#BiocManager::install("marray")
library("marray")
mfuzzColorBar()
mfuzz.plot2(eset,mf,mfrow=c(2,3),
            colo="fancy",
            #  ax.col="red",
             bg = "#fdd49e",
            #col.axis="red",col.lab="white",
              #col.main="green",
            col.sub="blue",col="blue",cex.main=1.3,cex.lab=1.1)


pdf("mfuzz0.pdf")
mfuzz.plot2(eset, mf=mf,mfrow=c(3,3),centre=TRUE,x11=F,centre.lwd=0.2,
            time.labels=colnames(exp_f))
dev.off()
 
 
###------Step7 output result------###

dir.create(path="mfuzz",recursive = TRUE)


for(i in 1:4){
  potname<-names(mf$mfuster[unname(mf$mfuster)==i])
  write.csv(mf[[4]][potname,i],paste0("mfuzz","/mfuzz_",i,".csv"))
}

# cluster matrix
red_mfuster <- mf$mfuster
red_exp_mfuster <- cbind(exp_f[names(red_mfuster), ], red_mfuster)

write.table(red_exp_mfuster , 'all_exp_mfuster .txt', sep = '\t', col.names = NA, quote = FALSE)

################################################################################

