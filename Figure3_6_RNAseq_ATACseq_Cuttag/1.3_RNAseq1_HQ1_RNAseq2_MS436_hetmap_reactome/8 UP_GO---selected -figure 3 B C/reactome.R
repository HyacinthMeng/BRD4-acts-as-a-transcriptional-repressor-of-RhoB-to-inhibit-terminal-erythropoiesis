library(tidyverse)
library(ggplot2)

GO = read.table("RNA1_RNA2_UP_GO_selected.txt",header=TRUE,sep="\t",quote = "")  #读取第1部分enrichGO分析输出的文件eG。

View(GO)


#计算Enrichment Factor
GO <- separate(data=GO, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") 
GO <- separate(data=GO, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
GO <- mutate(GO, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 



# pvalue
p1 <- ggplot(GO,aes(Group,Description)) +  
  geom_point(aes(size=Count,color=pvalue)) +  
  
  # 比例大小--缩放
  scale_size("Count",range=c(6,10))+
  #scale_color_gradient(low="red",high = "green") +  
  # scale_color_gradient(low = '#d90424', high = '#374a89')+
  scale_color_gradient(low = 'red',high="#d90424")+
  labs(color="pvalue",size="Count", shape="Group",x="GROUP",y="GO:BP",title="GO:BP") +  
  theme_bw()+
  theme(panel.grid=element_blank())
  
p1
?ggsave
ggsave("RNA1_RNA2_UP_GO_combine_pvalue.pdf", width = 6, height = 5)

# p.adjust
p <- ggplot(GO,aes(Group,Description)) +  
  geom_point(aes(size=Count,color=p.adjust)) +  
  
  # 比例大小--缩放
  scale_size("Count",range=c(6,10))+
  #scale_color_gradient(low="red",high = "green") +  
  #scale_color_gradient(low = '#d90424', high = '#374a89')+
  scale_color_gradient(low = 'red',high="#d90424")+
  labs(color="p.adjust",size="Count", shape="Group",x="GROUP",y="GO:BP",title="GO:BP") +  
  theme_bw()+
  theme(panel.grid=element_blank())
p


ggsave("RNA1_RNA2_UP_GO_combine_padjust.pdf", width = 6, height = 5)






#===============================================================================
