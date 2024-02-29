#载入相关的R包；
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(dplyr)
#读入Excel数据；
dt <- read_csv("KEGG.csv")

dt1 = read.table("RNA1_UP_reactome.txt",header=TRUE,sep="\t",quote = "")

#预览数据；
head(dt1)
View(dt1)
class(dt1)
#提取用于作图的列；
#df <- dt[,c(1,2,3,7)]
#对InTerm_InList列进行拆分；
#df <- separate(df, InTerm_InList, sep = "/", into = c("Count", "InTerm"))

#转换数据类型；
dt1$Count <- as.numeric(dt1$Count)


#根据Pvalue选择Top10的pathway进行绘制；

fig_data <- dt1 %>% top_n(-2,pvalue) %>% arrange(Count)
fig_data


#转成因子，防止重新排序；
fig_data$Description <- factor(fig_data$Description,
                               levels=fig_data$Description,ordered=TRUE)
#预览数据；
fig_data







# 2. 绘制柱状图

#建立数据与图形（点）的映射关系，即确定点的(x,y)坐标，绘制散点图；
p1<-ggplot(fig_data, aes(Count,Description))+
  geom_point(aes(size=Count,color=pvalue))
#隐藏坐标轴的标题；
p2<-p1+labs(y="")
p2

#自定义颜色渐变；
p3<-p2+scale_colour_gradient(low="red",high="yellow")
p3

#设置x轴范围，避免点的溢出绘图区；
p4<-p3+scale_x_continuous(limits = c(5.5, 12.5),
                          breaks = c(6,8, 10, 12),
                          label = c("6","8", "10", "12"))
p4


#3. 自定义主题



#应用自带主题；
p4+theme_light()

#自定义图表主题，对图表做精细调整；
top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2
#隐藏坐标轴，并对字体样式、颜色、刻度长度等进行限定；
mytheme<-theme_light()+
  theme(axis.text=element_text(family = "sans",colour ="gray20",size = 10),
        panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.6,colour = "gray20"),
        axis.ticks.length = unit(1.2,units = "mm"),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
p4+mytheme









#
