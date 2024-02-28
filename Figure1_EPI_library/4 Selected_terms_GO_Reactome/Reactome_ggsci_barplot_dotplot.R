library(tidyverse)
library(ggplot2)

Reactome = read.table("Reactome_selected.txt",header=TRUE,sep="\t",quote = "") 

View(Reactome)


# count Enrichment Factor
Reactome <- separate(data=Reactome, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") 
Reactome <- separate(data=Reactome, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
Reactome <- mutate(Reactome, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
Reactome$p.adjust


# pvalue ---------------------------------- ---------------

p1 <- ggplot(Reactome,aes(Count,reorder(Description,Count))) +  # reorder(Description,Count)
      
        geom_point(aes(size=Count,
                       color=-1*log10(pvalue))) +
       
        # scale
        scale_size("Count",range=c(1,8))+
        
        scale_color_gradient(low = "#FFC1C4",high="#C81B2E",  
                           limits=c(0,40),     # setting range
        )+ 
  
        # labs & theme
        labs(color="-log10(pvalue)",size="Count", 
             shape="Group",x="Count",y="Reactome",title="Reactome") +
  
        theme_bw()+
        theme(panel.grid=element_blank())

p1

ggsave("Reactome_pvalue.pdf", width =8 , height = 5)

#-------------------------------------------------------------------------------