# GO  barplot & dotplot
# BY MY 2023


library(ggplot2)
library(RColorBrewer)


display.brewer.all()
color <- brewer.pal(3,"Dark2")
colorl <- rep(color,each=10)

# load data
pdata <- read.table("GO_selected.txt",header = T,sep="\t")
View(pdata)

# factor and order
level = pdata[,3]
pdata$Description <- factor(pdata$Description,levels = level)



# GO barplot 1--------------------------------------------------------------

p1 = ggplot(pdata) +
    aes(x = Description, y = Count, fill = ONTOLOGY) +
    geom_bar(stat = "identity",colour="black") +
    scale_fill_manual(values =color)+
    theme(
      axis.title=element_text(size=15,face="plain",color="black"),
      axis.text = element_text(size=12,face="plain",color="black"),
      axis.text.x = element_text(angle = 90,colour = colorl,hjust=1,vjust=0.6),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 8, face = "bold"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      #legend.direction = "horizontal",
      legend.position = c(0.8,0.9),
      legend.background = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = "black"),
      plot.background = element_blank())+
  
    # changge direction
    coord_flip()


# ggsci color style
library("ggsci")
p2_npg <- p1 + scale_fill_npg()
p2_npg

ggsave("GO_barplot_ggsci.pdf", width =8 , height = 5 )







# dotplot ---------------------------------------------------------------------


p3 = ggplot(pdata) +
     aes(x = Description, y =-log10(pvalue), fill = ONTOLOGY,size=Count) +
     geom_point(shape=21,color="black") +
  
     #scale size
     scale_size("Count",range=c(2,5))+
     scale_fill_manual(values =color)+
    
     theme(
      axis.title=element_text(size=15,face="plain",color="black"),
      axis.text = element_text(size=12,face="plain",color="black"),
      axis.text.x = element_text(angle = 90,hjust=1,vjust=0.6),
      axis.title.y = element_blank(),
      axis.text.y = element_text(colour = colorl),
      #legend.title = element_blank(),
      legend.text = element_text(size = 8, face = "bold"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      #legend.direction = "horizontal",
      #legend.position = c(0.5,0.9),
      legend.background = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = "black"),
      plot.background = element_blank())+
  
    xlab("Counts")+
    coord_flip()

library("ggsci")
p2_npg <- p3 + scale_fill_npg()
p2_npg

ggsave("GO_dotplot_ggsci.pdf", width =8 , height = 5 )

#-------------------------------------------------------------------------------


