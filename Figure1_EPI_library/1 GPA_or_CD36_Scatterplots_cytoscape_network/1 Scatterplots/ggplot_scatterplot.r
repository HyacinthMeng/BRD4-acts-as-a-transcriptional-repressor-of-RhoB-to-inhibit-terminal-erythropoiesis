# install.packages("ggrepel")
library(ggrepel)
library(ggplot2)
library(tidyverse)
library(readxl)
library(ggsci)




getwd()
setwd("")
# dat = read.csv("ZM_EPI_CHM_LOGFCC.csv",sep=",",header = TRUE)

dat = read_excel("./final screening libarary result.xlsx")

colnames(dat) 
dim(dat)

View(dat)

#   You will get the sampe results dependeds on data features
# filter data
# GPA_CD36_Diff_0.5 = dat %>%
#    filter(abs(LOGFC_GPA) > 0.5  | LOGFC_CD36 < -0.5)

# filter data
GPA_CD36_Diff_0.5 = dat %>%
    filter(abs(LOGFC_GPA) > 0.5  | abs(LOGFC_CD36) > 0.5)


#GPA_CD36_Diff_0.5 = dat %>%
#    filter(LOGFC_GPA > 0.5  | LOGFC_CD36 < -0.5)

write.table(GPA_CD36_Diff_0.5,"GPA_CD36_Diff_0.5.txt",sep="\t",quote = FALSE,row.names = FALSE)


cat(sort(unique(dat$Target_Group)))

cat(sort(unique(GPA_CD36_Diff_0.5$Target_Group)))

GPA_CD36_Diff_0.5$Target_Group <- factor(GPA_CD36_Diff_0.5$Target_Group, 
		 levels = c("Bromodomain proteins (BRD)",
                    "Histone acetyltransferase and histone deacetylases",
                    "Histone methyltransferase and histone demethylase",
                    "Topoisomerase(TOP1, 2A,2B)",
                    "Others"))

p = ggplot(dat,aes(x=LOGFC_CD36,y=LOGFC_GPA )) +
    geom_point(color="grey")+
    #geom_point(color="grey")+#FF0000
    geom_point(aes(x=LOGFC_CD36,y=LOGFC_GPA,
                color= Target_Group, size = 0.8),GPA_CD36_Diff_0.5)+
    #geom_point(aes(x=LOGFC_CD36,y=LOGFC_GPA, color="red"),CD36_diff)+
    
	# add line
    geom_vline(xintercept=0) +
    geom_hline(yintercept =0)+
	
	# add text
    #geom_text(data=GPA_diff,mapping = aes(label=ID)) +
    geom_text_repel(data = GPA_CD36_Diff_0.5,
              mapping = aes(label=ID),
              vjust = "inward", hjust = "inward",
              angle = 0, size=2,
              max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
    # geom_text_repel(data = GPA_CD36_Diff_0.5, 
    #            mapping = aes(label=MOLENAME), vjust = "inward", hjust = "inward",
    #          angle = 0, size=2.3,
    #          max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
    theme_bw()+
    theme(panel.grid =element_blank() )+
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.position = c(0.45, 1),
            legend.justification = c(1,1.2))+
    scale_size_continuous(guide = FALSE)


# ggsci--- choose color 
p2<-p+scale_color_npg()+scale_fill_npg()
p2
ploty(p2)
ggsave("Epi_Group_f2.pdf",width = 22, height = 18, units = c("cm"))











##########################################################################
 
## 差异
View(GPA_CD36_insect)
ggplot(dat,aes(x=LOGFC_CD36,y=LOGFC_GPA )) +
    geom_point(color="#00b0eb")+
    #geom_point(color="grey")+#FF0000
    geom_point(aes(x=LOGFC_CD36,y=LOGFC_GPA,color="#0d1fc0", size = 1.5),GPA_diff_0.5_up)+
    geom_point(aes(x=LOGFC_CD36,y=LOGFC_GPA,color="#e20612", size = 1.5),GPA_CD36_diff_0.5_dn)+
    #geom_point(aes(x=LOGFC_CD36,y=LOGFC_GPA, color="red"),CD36_diff)+
    # 辅助线
    geom_vline(xintercept=c(-0.5,0.5),lty=4,lwd=0.8) +
    geom_hline(yintercept =c(-0.5,0.5),lty=4,lwd=0.8)+
    #geom_text(data=GPA_diff,mapping = aes(label=ID)) +
   # geom_text(data = GPA_diff_0.5_up,
   #           mapping = aes(label=ID),
   #           vjust = "inward", hjust = "inward",
   #           angle = 0, size=2)+
    geom_text_repel(data = GPA_diff_0.5_up,
              mapping = aes(label=ID),
              vjust = "inward", hjust = "inward",
              angle = 0, size=2)+
    geom_text_repel(data = GPA_CD36_diff_0.5_dn,
              mapping = aes(label=ID),
              vjust = "inward", hjust = "inward",
              angle = 0, size=2)+
    theme_bw()+
    theme(panel.grid =element_blank() )+
    theme(panel.background = element_rect(fill = "white")) +
    theme(legend.position = "none")

# "#e20612","#ffd401",
ggsave("小分子_red_blue.pdf",width = 22, height = 20, units = c("cm"))




install.packages("plotly")
install.packages("ggplotly")
library("plotly")
ggplotly(P1)
