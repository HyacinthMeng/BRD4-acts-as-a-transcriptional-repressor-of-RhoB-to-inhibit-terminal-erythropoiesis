# choose color style ：https://colorbrewer2.org/#type=sequential&scheme=YlOrBr&n=9
# BY MY 2023

library(RColorBrewer)
library(ggplot2)


pdata <- read.table("GO_selected.txt",header = T,sep="\t")
View(pdata)


# # attention: the order of level ----------------------------------------------
# factor
level = pdata[,3]  

# from A to Z
pdata$Description <- factor(pdata$Description,levels = level)

# from Z to A
pdata$Description <- factor(pdata$Description,levels = rev(level)) 


# Draw GO  ---------------------------------------------------------------------
ggplot(pdata) +
  aes(x = Description, y =-log10(p.adjust), 
      fill = ONTOLOGY,size=Count) +
  
  ## point
  geom_point(shape=21,
             color="black") +
  # scale size
  scale_size("Count",range=c(2,7))+     

  ## rect  
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=0 ,
                xmax=10.5, 
                ymin=-Inf, 
                ymax=Inf),
            fill='#ffc7c3',alpha = .12)+  
  
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=10.5 ,
                xmax=20.5, 
                ymin=-Inf, 
                ymax=Inf),
            fill='#d1e5c0',alpha = .08)+
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=20.5 ,
                xmax=30.5, 
                ymin=-Inf, 
                ymax=Inf),
            fill='#aad1e6',alpha = .08)+
  
  # add p0.05 rect   (log10(0.05)  == -1.30103)  
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=0 ,
                xmax=30.5, 
                ymin=-Inf, 
                ymax=1.30103),
            fill='#cbd5e8',alpha = .1)+
  
  geom_point(shape=21,
             color="black") +
  
  
  # theme settings
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.text.x = element_text(angle = 45,hjust=1,vjust=1),    
    
    legend.text = element_text(size = 8, face = "bold"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"), 
    plot.margin = margin(1,0.5,0.5,4, unit = "cm")) +
  
  theme(panel.grid.minor = element_line(size = 1,color = "#fafafa"))+
  xlab("GO enrichment analysis")+
  ylab("-log10(p.adjust)")+
  coord_flip() +  # 坐标旋转 
  
  # range setting ------
  # scale color
  scale_color_gradient(low = "#FF9EA3",high="#E82338", 
                       limits=c(2,40) 
                       )+
  # P value                                      )+
  scale_y_continuous(breaks=seq(0, 35, 5),limits = c(0, 35))  
  
 
   
ggsave("GO_pvalue1.30103.pdf", width =8 , height = 6.5)



 



#===============================================================================
