install.packages("tidyverse")
library(tidyverse)

rt = read.table("alldiffgene.txt",header = T)

rt2 = read.table("GSE53983_All_hs_countData.txt",sep="\t",header = T)


rt3=left_join(rt,rt2 ,by ="Gene")

View(rt3) #5853

write.table(rt3,"mfuzz_gene.txt",sep = "\t")



rowMeans(df)

rt4=rt3%>% 
  mutate(proerythroblast=rowMeans(rt3[,2:4]) ) %>% 
  mutate(early_basophilic=rowMeans(rt3[,5:7]) ) %>% 
  mutate(late_basophilic=rowMeans(rt3[,8:10]) ) %>% 
  mutate(polychromatic=rowMeans(rt3[,11:13]) ) %>% 
  mutate(orthochromatic=rowMeans(rt3[,11:13]) )

View(rt4)

rt5= rt4 %>% 
  select(Gene,proerythroblast,early_basophilic,late_basophilic,polychromatic,orthochromatic)
View(rt5)


rt6=as.matrix(rt5)
rownames(rt6)=rt6[,1]
exp=rt6[,2:ncol(rt6)]
dimnames=list(rownames(exp),colnames(exp))
rt6=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

rt6=log2(rt6+1)

View(rt6)

write.table(rt6,"mfuzz_exp_log2.txt",sep = "\t")

################################################################################

