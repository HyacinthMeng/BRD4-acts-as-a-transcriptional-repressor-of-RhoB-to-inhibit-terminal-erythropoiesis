###https://cran.r-project.org/web/packages/ggplot2/index.html
   #It contains the inroduction of pheatmap#
####----- BY MY   2021 rewrite @ZJU-----###

library(ggplot2)
###load file####################################################################
diff_stat <- read.table("limma.txt",,row.names=1,header=T,sep="\t")

#根据差异水平赋值颜色，例如 |log2FC| >=1 & FDR p-value < 0.05
diff_stat[which(diff_stat$adj.P.Val < 0.05 & diff_stat$logFC > 1),'diff'] <- 'up'
diff_stat[which(diff_stat$adj.P.Val< 0.05 & diff_stat$logFC < -1),'diff'] <- 'dowm'
diff_stat[!(diff_stat$diff %in% c('up', 'dowm')),'diff'] <- 'no'

p1 <- ggplot(diff_stat, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = diff), size = 1.0) +
  scale_colour_manual(limits = c('up', 'dowm', 'no'), values = c('red', 'blue', 'grey'), labels = c('UP', 'DOWN', 'NOT')) +
  labs(x = 'log2 Fold Change', y = '-log10 adj.P.Val')+
  
  
  geom_vline(xintercept = c(-1, 1), color = 'gray', linetype = 2, size = 0.5)+
  geom_hline(yintercept = -log10(0.05), color = 'gray', linetype = 2, size = 0.5)+
  
  theme(panel.grid =element_blank())+ 
  theme(title=element_text(size=12),
                 axis.title.x=element_text(size=12,),
                 axis.title.y=element_text(size=12,)) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.2, 0.9))

p1
#ggsave('volcano_plot7.pdf', p1, width = 6, height = 4)
#ggsave('volcano_plot1.png', p1, width = 3.5, height = 4.5)



#################################################################################
