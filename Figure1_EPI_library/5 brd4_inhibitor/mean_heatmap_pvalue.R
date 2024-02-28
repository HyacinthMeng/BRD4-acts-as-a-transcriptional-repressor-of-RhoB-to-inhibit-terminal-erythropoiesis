
# pheatmap & Diff 
# BY  MY 2023

# load packages
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(rstatix)


rt1 <- read.csv("brd4_inhibitor_used.csv",header = T,row.names = 1)

#rt1 = read_csv("brd4_inhibitor.csv")
#View(rt1)


annotation_row = data.frame(
  NongDu  = factor(c("con","con","con",rep ( rep(c("LOW", "Middle", "High"),each = 3),5) )),
  DrugClass = factor(rep(c("CTR","JQ", "ARV","A1874","Dbet57","MS436"), c(3, 9, 9, 9, 9, 9))))

row.names(annotation_row)=rownames(rt1)

# colors
ann_colors = list(
  NongDu = c( con="#FFF5F0",LOW="#FEE0D2", Middle="#FC9272",High="#CB181D"  ),
  DrugClass = c(CTR= "#FFF5F0",JQ = "#66C2A5" ,ARV= "#8DA0CB", A1874= "#E78AC3" ,Dbet57="#A6D854",MS436= "#FFD92F")
)


pheatmap(rt1, 
         #annotation_col=annotation_col, 
         annotation_row = annotation_row, 
         annotation_colors = ann_colors,
         color = colorRampPalette(c("navy", "white", "#CB181D"))(100),
         cluster_cols =F,
         cluster_rows  =F,
         scale="column",
          gaps_row = c(3,12,21,30,39,48),
         gaps_col = 1,
         border_color=NA,
         #border_color="white",
         # fontsize_row=5,
         fontsize_col=10,
         labels_row = NULL,
         labels_col = NULL,
         show_rownames = FALSE, 
         show_colnames =TRUE)




# get mean according to group -----------------------------------------------------

rt1 <- read.csv("brd4_inhibitor.csv",header = T)

rt_mean = rt1 %>%
  group_by(Group,Dose) %>% 
  summarise(mean_CD36 = mean(CD36),
            mean_GPA = mean(GPA))

# write.table(rt_mean,"mean.csv",sep = ",",row.names = FALSE)






##=== mean heatmap
rt1 <- read.csv("mean.csv",header = T,row.names = 1)

annotation_row = data.frame(
  NongDu  = factor(c("con",rep ( rep(c("LOW", "Middle", "High"),each = 1),5) )),
  DrugClass = factor(rep(c("CTR","MS436","JQ", "ARV","A1874","Dbet57"), c(1, 3, 3, 3, 3, 3))))

row.names(annotation_row)=rownames(rt1)


ann_colors = list(
  NongDu = c( con="#FFF5F0",LOW="#FEE0D2", Middle="#FC9272",High="#CB181D"  ),
  DrugClass = c(CTR= "#FFF5F0",JQ = "#66C2A5" ,ARV= "#8DA0CB", A1874= "#E78AC3" ,Dbet57="#A6D854",MS436= "#FFD92F")
)



pdf("mean_heatmap.pdf",width = 6,height = 10)
pheatmap(rt1, 
         #annotation_col=annotation_col, 
         annotation_row = annotation_row, 
         annotation_colors = ann_colors,
         color = colorRampPalette(c("navy", "white", "#CB181D"))(100),
         cluster_cols =F,
         cluster_rows  =F,
         scale="column",
         #scale="row",
         gaps_row = c(1,4,7,10,13,16),
         gaps_col = 1,
         border_color=NA,
         #border_color="white",
         # fontsize_row=5,
         fontsize_col=10,
         labels_row = NULL,
         labels_col = NULL,
         show_rownames = FALSE, 
         show_colnames =TRUE )

dev.off()



# add p value -----------------------------------------------------------------
rt2 <- read.csv("brd4_inhibitor.csv",header = T)

rt2 = rt2 %>% 
  mutate(G = paste(Group,"Dose"))



# t test：
stat.test_CD36 = rt2 %>%
  t_test(
    CD36 ~  G ,
    p.adjust.method = "holm"
  )

#View(stat.test_CD36)
write.table(stat.test_CD36,"p_CD36.csv",sep = ",")


CD36 = cbind(stat.test_CD36$p,stat.test_CD36$group2)
#View(CD36[1:45,])

pmt <- stat.test_CD36$p
pmt

# mark pvalue level
if (!is.null(pmt)){
  ssmt <- pmt<= 0.01
  pmt[ssmt] <-'**'
  smt <- pmt >0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt]<- ''
} else {
  pmt <- F
}

pmt
# keep the same order with the Drug dose
pmt1 = c("" , "", "*",  "*",    "*",  "*" , "**", "" ,  "*" , "*",  "" ,  "*"  ,"**", "" ,  "*" , "**"  )
pmt2 = c("" ,"",   "" ,  "" ,  "" ,  ""  , ""  , "*",  "*" , "*",  "",  "" ,  "",   "" ,  ""  , "**"  )
pmt3 = cbind(pmt1,pmt2)

pdf("mean_pvale_heatmap.pdf",width = 10,height = 10)
pheatmap(rt1, 
         #annotation_col=annotation_col, 
         annotation_row = annotation_row, 
         annotation_colors = ann_colors,
         color = colorRampPalette(c("navy", "white", "#CB181D"))(100),
         cluster_cols =F,
         cluster_rows  =F,
         scale="column",
         #scale="row",
         gaps_row = c(1,4,7,10,13,16),
         #gaps_col = 1,
         border_color=NA,
         #border_color="white",
         # fontsize_row=5,
         fontsize_col=10,
         labels_row = NULL,
         labels_col = NULL,
         show_rownames = FALSE, 
         show_colnames =TRUE,
         display_numbers = pmt3,
         fontsize_number = 20, number_color = "Black")

dev.off()


# #----------------------------------------------------------------------------






