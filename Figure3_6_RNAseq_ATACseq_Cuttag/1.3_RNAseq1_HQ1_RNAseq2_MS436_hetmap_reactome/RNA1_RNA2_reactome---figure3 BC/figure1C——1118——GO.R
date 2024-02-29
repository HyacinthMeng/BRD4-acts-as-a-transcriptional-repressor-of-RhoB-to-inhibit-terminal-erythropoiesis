library(ggplot2)
library(RColorBrewer)


pdata <- read.table("GO_selected.txt",header = T,sep="\t")
View(pdata)

# 创建一个factor,通过数据自身 赋值与level )
level = pdata[,3]
pdata$Description <- factor(pdata$Description,levels = level)






# 20231119---理解 图层设计理念的应用-----------
# 球的颜色代表 GO term 分类 -----Go term----------------------------------------------------
# 添加p =0.05 的阴影区
ggplot(pdata) +
  aes(x = Description, y =-log10(p.adjust), 
      fill = ONTOLOGY,size=Count) +
  geom_point(shape=21,
             color="black",
             #   color="black"
  ) +
  #scale_fill_hue() +
  # 比例大小--缩放
  scale_size("Count",range=c(2,7))+     # 定义尺寸相对大小
  #scale_fill_manual(values =color)+   # 自定义颜色
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=0 ,
                xmax=10.5, 
                ymin=-Inf, 
                ymax=Inf),
            fill='#ffc7c3',alpha = .12)+   #"#e7c6c3""   #d1e5c0  #FFB6C1 粉色
  
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
  
  # 添加p 值阴影区域
  # geom_rect(data=pdata, inherit.aes=TRUE,
 #           aes(xmin=0 ,
  #              xmax=30.5, 
  #              ymin=-Inf, 
   #             ymax=1.30103),
  #          fill='#cbd5e8',alpha = .1)+
  
  
  
  geom_point(shape=21,
             color="black",
             #   color="black"
  ) +
  
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.text.x = element_text(angle = 45,hjust=1,vjust=1),    #横坐标字体倾斜设置
   # axis.title.y = element_blank(), 
  #  axis.text.y = element_text(colour = colorl),
    
    #legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "bold"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"), # 调整legend和图的边距 unit = "pt"
    #legend.direction = "horizontal",
    #legend.position = c(0.5,0.9), 
    # legend.background = element_blank(),
    
    plot.margin = margin(1,0.5,0.5,4, unit = "cm")
    
    #panel.background = element_rect(fill = "transparent",colour = "black"),
    #  plot.background = element_blank()
  ) +
  # theme_bw()+  #去除灰色背景
  # theme_classic()+   # 去除网格
   theme_bw() + theme(panel.grid=element_blank())+ #去除网格线
  
  
  theme(panel.grid.minor = element_line(size = 1,color = "#fafafa"))+ #添加网格线
  xlab("GO enrichment analysis")+
  ylab("-log10(p.adjust)")+
  coord_flip() +  # 坐标旋转 
  # scale_x_reverse()    #先对横坐标的level 进行数据反转 pdata$Description <- factor(pdata$Description,levels = rev(level))   #
  # scale_y_reverse()
  scale_color_gradient(low = "#FF9EA3",high="#E82338",  #  "#FFE0E1"  #low = 'blue',high="#d90424" low = "#FF9EA3",high="#E82338"
                       limits=c(2,40)  # 自定义范围   limits=c(1,12)
                       #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                       #                     labels = c("3","10","20","30","40")  # 自定义标签
  )+
    scale_y_continuous(breaks=seq(0, 35, 5),limits = c(0, 35))  # 设置横坐标  GO
  
  
ggsave("有p 阴影区域——GO_combine_pvalue——1121_bubble_pvalue下线=1.30103_width_7.pdf", width =8 , height = 6.5)


ggsave("无p 阴影区域——GO_combine_pvalue——1121_bubble_pvalue下线=1.30103_width_7.pdf", width =8 , height = 6.5)










pdata$Description <- factor(pdata$Description,levels = rev(level))   

# 正序
ggplot(pdata) +
  aes(x = Description, y =-log10(p.adjust), 
      fill = ONTOLOGY,size=Count) +
  geom_point(shape=21,
             color="black",
             #   color="black"
  ) +
  #scale_fill_hue() +
  # 比例大小--缩放
  scale_size("Count",range=c(2,7))+     # 定义尺寸相对大小
  #scale_fill_manual(values =color)+   # 自定义颜色
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=0 ,
                xmax=10.5, 
                ymin=-Inf, 
                ymax=Inf),
            fill='#aad1e6',alpha = .08)+   #"#e7c6c3""   #d1e5c0  #FFB6C1 粉色
  
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
            fill='#ffc7c3',alpha = .15)+
  
   # 添加p 值阴影区域
   geom_rect(data=pdata, inherit.aes=TRUE,
             aes(xmin=0 ,
                xmax=30.5, 
                ymin=-Inf, 
                ymax=1.30103),
                fill='#cbd5e8',alpha = .1)+
  
  
  
geom_point(shape=21,
           color="black",
           #   color="black"
) +
  
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.text.x = element_text(angle = 45,hjust=1,vjust=1),    #横坐标字体倾斜设置
    # axis.title.y = element_blank(), 
    #  axis.text.y = element_text(colour = colorl),
    
    #legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "bold"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"), # 调整legend和图的边距 unit = "pt"
    #legend.direction = "horizontal",
    #legend.position = c(0.5,0.9), 
    # legend.background = element_blank(),
    
    plot.margin = margin(1,0.5,0.5,4, unit = "cm")
    
    #panel.background = element_rect(fill = "transparent",colour = "black"),
    #  plot.background = element_blank()
  ) +
  # theme_bw()+  #去除灰色背景
  # theme_classic()+   # 去除网格
  theme_bw() + theme(panel.grid=element_blank())+ #去除网格线
  
  
  theme(panel.grid.minor = element_line(size = 1,color = "#fafafa"))+ #添加网格线
  xlab("GO enrichment analysis")+
  ylab("-log10(p.adjust)")+
  coord_flip() +  # 坐标旋转 
  # scale_x_reverse()    #先对横坐标的level 进行数据反转 pdata$Description <- factor(pdata$Description,levels = rev(level))   #
  # scale_y_reverse()
  scale_color_gradient(low = "#FF9EA3",high="#E82338",  #  "#FFE0E1"  #low = 'blue',high="#d90424" low = "#FF9EA3",high="#E82338"
                       limits=c(2,40)  # 自定义范围   limits=c(1,12)
                       #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                       #                     labels = c("3","10","20","30","40")  # 自定义标签
  )+
  scale_y_continuous(breaks=seq(0, 35, 5),limits = c(0, 35))  # 设置横坐标  GO


ggsave("排序——有p 阴影区域——GO_combine_pvalue——1121_bubble_pvalue下线=1.30103_width_7.pdf", width =8 , height = 6.5)


ggsave("排序——无p 阴影区域——GO_combine_pvalue——1121_bubble_pvalue下线=1.30103_width_7.pdf", width =8 , height = 6.5)




















#######################################################################

# 未添加p =0.05 的阴影区
ggplot(pdata) +
  aes(x = Description, y =-log10(p.adjust), 
      fill = ONTOLOGY,size=Count) +
  geom_point(shape=21,
             color="black",
             #   color="black"
  ) +
  #scale_fill_hue() +
  # 比例大小--缩放
  scale_size("Count",range=c(2,7))+     # 定义尺寸相对大小
  #scale_fill_manual(values =color)+   # 自定义颜色
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=0 ,
                xmax=10.5, 
                ymin=-Inf, 
                ymax=Inf),
            fill='#aad1e6',alpha = .08)+   # #ffc7c3 粉色
  
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=10.5 ,
                xmax=20.5, 
                ymin=-Inf, 
                ymax=Inf),
            fill='#d1e5c0',alpha = .08)+
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=20.5 ,
                xmax=31, 
                ymin=-Inf, 
                ymax=Inf),
            fill='#ffc7c3',alpha = .08)+   #aad1e6 蓝色
  # 添加p 值阴影区
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=0 ,
                xmax=30.5, 
                ymin=-Inf, 
                ymax=1.30103),
            fill='#cbd5e8',alpha = .10)+
  
  
  geom_point(shape=21,
             color="black",
             #   color="black"
  ) +
  
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.text.x = element_text(angle = 45,hjust=1,vjust=1),    #横坐标字体倾斜设置
    # axis.title.y = element_blank(), 
    #  axis.text.y = element_text(colour = colorl),
    
    #legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "bold"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"), # 调整legend和图的边距 unit = "pt"
    #legend.direction = "horizontal",
    #legend.position = c(0.5,0.9), 
    # legend.background = element_blank(),
    
    plot.margin = margin(1,0.5,0.5,4, unit = "cm")
    
    #panel.background = element_rect(fill = "transparent",colour = "black"),
    #  plot.background = element_blank()
  ) +
  # theme_bw()+  #去除灰色背景
  # theme_classic()+   # 去除网格
  # theme_bw() + theme(panel.grid=element_blank())+ #去除网格线
  
  
  theme(panel.grid.minor = element_line(size = 1,color = "#fafafa"))+ #添加网格线
  xlab("GO enrichment analysis")+
  ylab("-log10(p.adjust)")+
  coord_flip() +  # 坐标旋转 
  # scale_x_reverse()    #先对横坐标的level 进行数据反转 pdata$Description <- factor(pdata$Description,levels = rev(level))   #
  # scale_y_reverse()
  scale_color_gradient(low = "#FF9EA3",high="#E82338",  #  "#FFE0E1"  #low = 'blue',high="#d90424" low = "#FF9EA3",high="#E82338"
                       limits=c(2,40)  # 自定义范围   limits=c(1,12)
                       #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                       #                     labels = c("3","10","20","30","40")  # 自定义标签
  )+
  scale_y_continuous(breaks=seq(0, 35, 5),limits = c(0, 35))  # 设置横坐标  GO


ggsave("排序——有阴影区——GO_combine_pvalue——1121_bubble_pvalue下线=1.30103_width_7.pdf", width =10, height = 5)

  
  
   
?coord_flip





pdata$Description <- factor(pdata$Description,levels = rev(level))   #




# 球的颜色代表 p值-----Go term----------------------------------------------------
# 圆圈的颜色是 p.adjust
ggplot(pdata) +
  aes(x = Description ,y=Count ,size=Count) +

  geom_point (aes(size=Count,color= -log10(p.adjust))
             #   color="black"
  ) +
  #scale_fill_hue() +
  # 比例大小--缩放
  scale_size("Count",range=c(2,7))+     # 定义尺寸相对大小
  # scale_fill_manual(values =color)+   # 自定义颜色
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=0 ,
                xmax=10.5, 
                ymin=-Inf, 
                ymax=Inf),
            fill='#ffc7c3',alpha = .08)+   #"#e7c6c3""   #d1e5c0  #FFB6C1 粉色
  
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
  
  geom_point (aes(size=Count,color= -log10(p.adjust))
              #   color="black"
  ) +
  
  
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.text.x = element_text(angle = 45,hjust=1,vjust=1),    #横坐标字体倾斜设置
    # axis.title.y = element_blank(), 
    #  axis.text.y = element_text(colour = colorl),
    
    #legend.title = element_blank(),
    legend.text = element_text(size = 8, face = "bold"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"), # 调整legend和图的边距 unit = "pt"
    #legend.direction = "horizontal",
    #legend.position = c(0.5,0.9), 
    # legend.background = element_blank(),
    
    plot.margin = margin(1,0.5,0.5,4, unit = "cm")
    
    #panel.background = element_rect(fill = "transparent",colour = "black"),
    #  plot.background = element_blank()
  ) +
  # theme_bw()+  #去除灰色背景
  # theme_classic()+   # 去除网格
  # theme_bw() + theme(panel.grid=element_blank())+ #去除网格线
  
  
  theme(panel.grid.minor = element_line(size = 1,color = "#fafafa"))+ #添加网格线
  xlab("GO enrichment analysis")+
  ylab("Count")+
  coord_flip() +  # 坐标旋转 
# scale_x_reverse()    #先对横坐标的level 进行数据反转 pdata$Description <- factor(pdata$Description,levels = rev(level))   #
# scale_y_reverse()
  # 定义颜色
  scale_color_gradient(low = "#FF9EA3",high="#E82338",  #  "#FFE0E1"  #low = 'blue',high="#d90424" low = "#FF9EA3",high="#E82338"
                       limits=c(2,40)  # 自定义范围   limits=c(1,12)
                       #                     breaks = c("3","10","20","30","40"), # 自定义间隔    # 
                       #                     labels = c("3","10","20","30","40")  # 自定义标签
  )+
  scale_y_continuous(breaks=seq(0, 35, 10),limits = c(0, 35))  # 设置横坐标  GO


ggsave("Reactome_combine_padjust——1121_bubble_pvalue下线=4.pdf", width =6 , height = 5)








#===============================================================================
