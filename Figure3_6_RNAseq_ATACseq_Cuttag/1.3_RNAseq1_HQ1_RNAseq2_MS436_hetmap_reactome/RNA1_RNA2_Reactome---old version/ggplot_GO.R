
# GO 结果可视化
# BY Hyacninth

# Enrichment Factor = GeneRatio/BgRatio**
# GeneRatio：基因比，分子是富集到此GO term上的基因数，而分母是所有得输入基因数。
# BgRatio：背景比，分子是此GO term得基因数，分母则是所有被GO注释的基因数。


library(tidyverse)
GO = read.table("GO.txt",header=TRUE,sep="\t",quote = "")  #读取第1部分enrichGO分析输出的文件eG。
head(GO)

#计算Enrichment Factor
GO <- separate(data=GO, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") 
GO <- separate(data=GO, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
GO <- mutate(GO, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 


# BP\MF\CC各取排名前10的term

GOBP <- GO %>%  
  filter(ONTOLOGY=="BP") %>%  
  filter(row_number() >= 1,row_number() <= 10)

GOCC <- GO %>%  
  filter(ONTOLOGY=="CC") %>%  
  filter(row_number() >= 1,row_number() <= 10)

GOMF <- GO %>%  
  filter(ONTOLOGY=="MF") %>%  
  filter(row_number() >= 1,row_number() <= 10)


# 合并成新的表格
GO10 <- rbind(GOBP,GOMF,GOCC)


library(ggplot2)
p <- ggplot(GO10,aes(enrichment_factor,Description)) +  
  geom_point(aes(size=Count,color=qvalue,shape=ONTOLOGY)) +  
  #scale_color_gradient(low="red",high = "green") +  
  scale_color_gradient(low = '#d90424', high = '#374a89')+
  labs(color="pvalue",size="Count", shape="Ontology",x="Enrichment Factor",y="GO term",title="GO enrichment") +  
  theme_bw()
p



#===============================================================================
#