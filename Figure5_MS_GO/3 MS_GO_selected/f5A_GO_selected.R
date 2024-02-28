


pdata = read.table("MS_selected_GO_selected.txt",header=TRUE,sep="\t",quote = "")  #读取第1部分enrichGO分析输出的文件eG。



level = pdata[,3]
# pdata$Description <- factor(pdata$Description,levels = level)
pdata$Description <- factor(pdata$Description,levels = rev(level))   

# View(pdata)


ggplot(pdata) +
  aes(x = Description ,y=-log10(p.adjust) ,size=Count) +
  
  geom_point (aes(size=Count,color= -log10(p.adjust)),
                 color="white"
  ) +
 
  scale_size("Count",range=c(2,7))+     
  scale_fill_manual(values ="white")+  
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=11.5 ,
                xmax=16, 
                ymin=-Inf, 
                ymax=Inf),
            fill='#fbb4ae',alpha = .08)+  
  
  # geom_rect(data=pdata, inherit.aes=TRUE,
  #          aes(xmin=0 ,
  #              xmax=11.5, 
  #              ymin=-Inf, 
  #              ymax=Inf),
  #          fill='#b3cde3',alpha = .12)+     
  
  geom_rect(data=pdata, inherit.aes=TRUE,
            aes(xmin=0 ,
                xmax=11.5, 
                ymin=-Inf, 
               ymax=Inf),
            fill='#ccebc5',alpha = .05)+   
  
  geom_point (aes(size=Count,color= -log10(p.adjust)) ) +
  
  
  scale_color_gradient(low = "#fcae91",high="#e41a1c",  
                       limits=c(3,40)  )+
  scale_y_continuous(breaks=seq(0, 35, 5),limits = c(0, 35)) +
  
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black"),
    axis.text.x = element_text(angle = 45,hjust=1,vjust=1),    
    
    legend.text = element_text(size = 8, face = "bold"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"), 
    plot.margin = margin(1,0.5,0.5,4, unit = "cm")) +
  
  theme_bw() + theme(panel.grid=element_blank())+
  
  xlab("GO enrichment analysis")+
  ylab("-log10(p.adjust)")+
  
  coord_flip()   

ggsave("MS_GO_combine_padjust.pdf", width =7.5 , height = 5)

##------------------------------------------------------------------------------

