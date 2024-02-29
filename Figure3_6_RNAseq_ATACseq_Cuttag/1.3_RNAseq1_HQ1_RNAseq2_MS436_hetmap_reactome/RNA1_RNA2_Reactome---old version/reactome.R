library(tidyverse)
library(ggplot2)

GO = read.table("merged_reactome_UP_selected.txt",header=TRUE,sep="\t",quote = "")  #读取第1部分enrichGO分析输出的文件eG。
head(GO)
View(GO)


#计算Enrichment Factor
GO <- separate(data=GO, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") 
GO <- separate(data=GO, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
GO <- mutate(GO, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 


# Draw fig.
p <- ggplot(GO,aes(enrichment_factor,Description)) +  
  geom_point(aes(size=Count,color=pvalue,shape=Group)) +  
  #scale_color_gradient(low="red",high = "green") +  
  scale_color_gradient(low = '#d90424', high = '#374a89')+
  labs(color="pvalue",size="Count", shape="Group",x="Enrichment Factor",y="GO term",title="GO enrichment") +  
  theme_bw()
p


# pvalue
p <- ggplot(GO,aes(Group,Description)) +  
  geom_point(aes(size=Count,color=pvalue)) +  
  
  # 比例大小--缩放
  scale_size("Count",range=c(6,12))+
  #scale_color_gradient(low="red",high = "green") +  
  # scale_color_gradient(low = '#d90424', high = '#374a89')+
  scale_color_gradient(low = 'red',high="#d90424")+
  labs(color="pvalue",size="Count", shape="Group",x="GROUP",y="Reactome Pathway",title="Reactome enrichment") +  
  theme_bw()+
  theme(panel.grid=element_blank())
p


# p.adjust
p <- ggplot(GO,aes(Group,Description)) +  
  geom_point(aes(size=Count,color=p.adjust)) +  
  
  # 比例大小--缩放
  scale_size("Count",range=c(6,12))+
  #scale_color_gradient(low="red",high = "green") +  
  #scale_color_gradient(low = '#d90424', high = '#374a89')+
  scale_color_gradient(low = 'red',high="#d90424")+
  labs(color="p.adjust",size="Count", shape="Group",x="GROUP",y="Reactome Pathway",title="Reactome enrichment") +  
  theme_bw()+
  theme(panel.grid=element_blank())
p









#===============================================================================
