###https://cran.r-project.org/web/packages/ggplot2/index.html
   #It contains the inroduction of pheatmap#
####----- BY MY   2021 rewrite @ZJU-----###

library(ggplot2)
###load file####################################################################

diff_stat <- read.table("Deseq2_UID_Diff.csv",
                        row.names=1,header=T,sep=",")


# |log2FC| >=1 & FDR p-value < 0.05
diff_stat[which(diff_stat$padj < 0.05 & diff_stat$log2FoldChange > 1),'diff'] <- 'up'
diff_stat[which(diff_stat$padj< 0.05 & diff_stat$log2FoldChange < -1),'diff'] <- 'dowm'
diff_stat[!(diff_stat$diff %in% c('up', 'dowm')),'diff'] <- 'no'

p1 <- ggplot(diff_stat, aes(x = log2FoldChange, y = -log10(padj))) +
  
        # point
        geom_point(aes(color = diff), size = 1.0) +
        scale_colour_manual(limits = c('up', 'dowm', 'no'), 
                            values = c('red', 'blue', 'grey'), 
                            labels = c('UP', 'DOWN', 'NOT')) +
        labs(x = 'log2 Fold Change', y = '-log10 padj') +
        
        # cutoff line
        geom_vline(xintercept = c(-1, 1), color = 'gray', 
                   linetype = 2, size = 0.5)+
        geom_hline(yintercept = -log10(0.05), color = 'gray', 
                   linetype = 2, size = 0.5) +
        
        
        # theme style
        theme_bw() +theme(panel.grid =element_blank()) +
        theme(title=element_text(size=12),
                       axis.title.x=element_text(size=12,),
                       axis.title.y=element_text(size=12,)) +
      
        theme(legend.title = element_blank(), 
              legend.key = element_rect(fill = 'transparent'), 
              legend.background = element_rect(fill = 'transparent'), 
              legend.position = c(0.2, 0.9))

p1


#ggsave('volcano_plot7.pdf', p1, width = 6, height = 4)


##-----------------------------------------------------------------------------
