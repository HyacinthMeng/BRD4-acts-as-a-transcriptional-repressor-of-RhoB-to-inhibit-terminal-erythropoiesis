library(tidyverse)
library(ggplot2)
library(scales)

GO = read.table("RNA1_RNA2_UP_Reactome_combine_21.txt",header=TRUE,sep="\t",quote = "")  #读取第1部分enrichGO分析输出的文件eG。

View(GO)


#计算Enrichment Factor
GO <- separate(data=GO, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") 
GO <- separate(data=GO, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
GO <- mutate(GO, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
GO$p.adjust
View(GO)

# PART1.1 ——pvalue 气泡图   demo 测试------------------------------------------------------


p1 <- ggplot(GO,aes(Count,reorder(Description,Count))) +    #reorder(Description,Count)
  geom_point(aes(size=Count,
                 color=-1*log10(pvalue))) +
  # 比例大小--缩放
  scale_size("Count",range=c(3,8))+
  # 自定义间隔--------------------------------------------------
scale_color_gradient(low = "#FF9EA3",high="#E82338",  #  "#FFE0E1"  #low = 'blue',high="#d90424"
                     limits=c(3,70),  # 自定义范围
                     #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                     #                     labels = c("3","10","20","30","40")  # 自定义标签
)+ 
  #labs(color=expression(-log[10](Qvalue))
  labs(color="-log10(pvalue)",size="Count", 
       shape="Group",x="Count",y="GO:BP",title="GO:BP") +
  theme_bw()+
  theme(panel.grid=element_blank())

p1
# library("ggsci")
# p2_npg <- p1 + scale_fill_npg()
# p2_npg

ggsave("Reactome_combine_pvalue——1107_bubble.pdf", width =8 , height = 5)






# PART1.1 pvalue  气泡图 -----------------------------------------------------
p1 <- ggplot(GO,aes(Group,reorder(Description,Count))) +    #reorder(Description,Count)
  geom_point(aes(size=Count,
                 color=-1*log10(pvalue))) +
  # 比例大小--缩放
  scale_size("Count",range=c(5,12))+
  scale_color_gradient(low = "#FF9EA3",high="#E82338",  #  "#FFE0E1"  #low = 'blue',high="#d90424"
                       limits=c(30,70),  # 自定义范围
                       #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                       #                     labels = c("3","10","20","30","40")  # 自定义标签
  )+ 
  # 自定义间隔--------------------------------------------------
  # scale_color_gradient2(low = muted("blue"),mid = "white",high=muted("#C81B2E"),  #  "#FFE0E1"  #low = 'blue',high="#d90424"
  #                    midpoint = 1.30103,  # log10(0.05)
  #                    limits=c(30,70),  # 自定义范围---根据数据交互结果
  #                    #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                      #                     labels = c("3","10","20","30","40")  # 自定义标签
#  )+ 
  #labs(color=expression(-log[10](Qvalue))
  labs(color="-log10(pvalue)",size="Count", 
       shape="Group",x="Count",y="Reactome",title="Reactome") +
  theme_bw()+
  theme(panel.grid=element_blank())

p1
# library("ggsci")
# p2_npg <- p1 + scale_fill_npg()
# p2_npg

ggsave("RNA1_2_Reactome_combine_pvalue——1107_bubble——final2.pdf", width =6 , height = 5)




# PART1.2——Padjust_bubble ---------------------------------------------------

p1 <- ggplot(GO,aes(Group,reorder(Description,Count))) +    #reorder(Description,Count)
  geom_point(aes(size=Count,
                 color=-log10(p.adjust))) +
  # 比例大小--缩放
  scale_size("Count",range=c(5,12))+
  
  scale_color_gradient(low = "#FF9EA3",high="#E82338",  #  "#FFE0E1"  #low = 'blue',high="#d90424"
                       limits=c(30,70),  # 自定义范围
                       #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                       #                     labels = c("3","10","20","30","40")  # 自定义标签
  )+ 
  # 自定义间隔--------------------------------------------------
  

 #  scale_color_gradient2(low = muted("blue"),mid = "white",high=muted("#C81B2E"),  #  "#FFE0E1"  #low = 'blue',high="#d90424"
 #                      midpoint = 1.30103,  # log10(0.05)
 #                     limits=c(30,70),  # 自定义范围---根据数据交互结果
                     #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                     #                     labels = c("3","10","20","30","40")  # 自定义标签
#  )+ 
  #labs(color=expression(-log[10](Qvalue))
  labs(color="-log10(p.adjust)",size="Count", 
       shape="Group",x="Count",y="Reactome",title="Reactome") +
  theme_bw()+
  theme(panel.grid=element_blank())

p1
# library("ggsci")
# p2_npg <- p1 + scale_fill_npg()
# p2_npg
log10(0.05)
ggsave("RNA1_2_Reactome_combine_p.adjust——1107_bubble_final2.pdf", width =7, height = 5)












p1 <- ggplot(GO,aes(Group,reorder(Description,Count))) +    #reorder(Description,Count)
  geom_point(aes(size=Count,
                 color=-1*log10(pvalue))) +
  # 比例大小--缩放
  scale_size("Count",range=c(5,12))+
  scale_color_gradient(low = "#FF9EA3",high="#E82338",  #  "#FFE0E1"  #low = 'blue',high="#d90424"
                       limits=c(30,70),  # 自定义范围
                       #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                       #                     labels = c("3","10","20","30","40")  # 自定义标签
  )+ 
  # 自定义间隔--------------------------------------------------
# scale_color_gradient2(low = muted("blue"),mid = "white",high=muted("#C81B2E"),  #  "#FFE0E1"  #low = 'blue',high="#d90424"
#                    midpoint = 1.30103,  # log10(0.05)
#                    limits=c(30,70),  # 自定义范围---根据数据交互结果
#                    #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
#                     labels = c("3","10","20","30","40")  # 自定义标签
#  )+ 
#labs(color=expression(-log[10](Qvalue))
labs(color="-log10(pvalue)",size="Count", 
     shape="Group",x="Count",y="Reactome",title="Reactome") +
  theme_bw()+
  theme(panel.grid=element_blank())

p1
# library("ggsci")
# p2_npg <- p1 + scale_fill_npg()
# p2_npg

ggsave("RNA1_2_Reactome_combine_pvalue——1107_bubble——final2.pdf", width =6 , height = 5)




# PART1.2——Padjust_bubble ---------------------------------------------------

p1 <- ggplot(GO,aes(Group,reorder(Description,Count))) +    #reorder(Description,Count)
  geom_point(aes(size=Count,
                 color=-log10(p.adjust))) +
  # 比例大小--缩放
  scale_size("Count",range=c(5,12))+
  
  scale_color_gradient(low = "#FF9EA3",high="#E82338",  #  "#FFE0E1"  #low = 'blue',high="#d90424"
                       limits=c(30,70),  # 自定义范围
                       #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                       #                     labels = c("3","10","20","30","40")  # 自定义标签
  )+ 
  # 自定义间隔--------------------------------------------------


#  scale_color_gradient2(low = muted("blue"),mid = "white",high=muted("#C81B2E"),  #  "#FFE0E1"  #low = 'blue',high="#d90424"
#                      midpoint = 1.30103,  # log10(0.05)
#                     limits=c(30,70),  # 自定义范围---根据数据交互结果
#                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
#                     labels = c("3","10","20","30","40")  # 自定义标签
#  )+ 
#labs(color=expression(-log[10](Qvalue))
labs(color="-log10(p.adjust)",size="Count", 
     shape="Group",x="Count",y="Reactome",title="Reactome") +
  theme_bw()+
  theme(panel.grid=element_blank())

p1
# library("ggsci")
# p2_npg <- p1 + scale_fill_npg()
# p2_npg
log10(0.05)
ggsave("RNA1_2_Reactome_combine_p.adjust——1107_bubble_final2.pdf", width =7, height = 5)













# PART2。1 pvalue_柱状图 -------------------------------------------------------
library(scales)
p1 <- ggplot(GO,aes(Count,reorder(Description,Count),
                    
                    fill=-log10(pvalue))) +    #reorder(Description,Count)
   
 geom_bar(aes(size=Count,
                 color=-1*log10(pvalue),fill=-log10(pvalue))),
stat = "identity") +
  # 比例大小--缩放
  scale_size("Count",range=c(3,8))+
  # 自定义间隔--------------------------------------------------
  scale_color_gradient2(low = muted("blue"),mid = "white",high=muted("#C81B2E"),  #  "#FFE0E1"  #low = 'blue',high="#d90424"
                      midpoint = 1.30103,  # log10(0.05)
                      limits=c(0,40),  # 自定义范围---根据数据交互结果
                      #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                      #                     labels = c("3","10","20","30","40")  # 自定义标签
)+ 
  #labs(color=expression(-log[10](Qvalue))
  labs(color="-log10(pvalue)",size="Count", 
       shape="Group",x="Count",y="Reactome",title="Reactome") +
  theme_bw()+
  theme(panel.grid=element_blank())

p1
# library("ggsci")
# p2_npg <- p1 + scale_fill_npg()
# p2_npg

ggsave("Reactome_combine_pvalue——1104_barplot——final.pdf", width =8 , height = 5)









ggplot(GO,aes(Count,reorder(Description,Count,fill=-log10(pvalue)))) +    #reorder(Description,Count)
  geom_bar(#aes(#size=Count,
                # color=-log10(p.adjust)),
           #stat = "identity"
           ) +
  # 比例大小--缩放
  scale_size("Count",range=c(3,8))+
  # 自定义间隔--------------------------------------------------
  scale_color_gradient2(low = muted("blue"),mid = "white",high=muted("#C81B2E"),  #  "#FFE0E1"  #low = 'blue',high="#d90424"
                      midpoint = 1.30103,  # log10(0.05)
                      limits=c(0,40),  # 自定义范围---根据数据交互结果
                      #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                      #                     labels = c("3","10","20","30","40")  # 自定义标签
)+ 
  #labs(color=expression(-log[10](Qvalue))
  labs(color="-log10(p.adjust)",size="Count", 
       shape="Group",x="Count",y="Reactome",title="Reactome") +
  theme_bw()+
  theme(panel.grid=element_blank())

p1

ggsave("Reactome_combine_p.adjust——1104_bubble_final.pdf", width =8 , height = 5)







