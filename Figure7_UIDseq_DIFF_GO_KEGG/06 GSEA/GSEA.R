###############################################################################
# BY MY 
###############################################################################


rm(list = ls())


# BiocManager::install("biomaRt")

library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(GseaVis)

genelist_input <- fread(file="R1_ID_ALL_gene.txt", header = T, sep='\t', 
                        data.table = F)

genename <- as.character(genelist_input[,1]) 
genename

# ID -->ENTREZID
gene_map <- select(org.Hs.eg.db, 
                   keys=genename, 
                   keytype="SYMBOL", 
                   columns=c("ENTREZID"))

colnames(gene_map)[1]<-"Gene"

write.csv(as.data.frame(gene_map),"geneID2.csv",row.names =F)



Inner_data<-inner_join(gene_map,genelist_input,by = "Gene")
Inner_data<-Inner_data[,-1]
Inner_data<-na.omit(Inner_data)
Inner_data$log2FC<-sort(Inner_data$log2FC,decreasing = T)
Inner_data

geneList = Inner_data[,2]
names(geneList) = as.character(Inner_data[,1])
geneList


# GSEA analyse ----------------------------------------------------------------

#GSEA——GO
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", 
                      ont="all",  minGSSize = 10, maxGSSize = 1000, 
                      pvalueCutoff=1)

#GSEA——KEGG
KEGG_gseresult <- gseKEGG(geneList,  minGSSize = 10, 
                          maxGSSize = 1000, pvalueCutoff=1)

#GSEA——Reactome
Go_Reactomeresult <- gsePathway(geneList,  minGSSize = 10, 
                                maxGSSize = 1000, pvalueCutoff=1)

# save
write.table (Go_gseresult, file ="Go_gseresult3.csv", 
             sep =",", row.names =TRUE)

write.table (KEGG_gseresult, file ="KEGG_gseresult3.csv", 
             sep =",", row.names =TRUE)

write.table (Go_Reactomeresult, file ="Go_Reactomeresult3.csv", 
             sep =",", row.names =TRUE)






# Draw figure   ----------------------------------------------------------------

# ridgeplot
ridgeplot(Go_gseresult,10)


# gseaplot：
gseaplot(Go_Reactomeresult,1,pvalue_table = TRUE) 


#gseaplot2
library(ggsci)
library(cowplot)

gseaplot2(Go_gseresult, 681, pvalue_table = TRUE)


# Draw single -----------------------------------------------------------------
ID= c("GO:0030036")
gseaplot2(Go_gseresult, ID, 
          title = "Actin cytoskeleton organization",
          #pvalue_table = TRUE,
          color="#4daf4a")


ID= c("GO:0030036")
gseaplot2(Go_gseresult, ID, 
          title = "Actin cytoskeleton organization",
          #pvalue_table = TRUE,
          color="#4daf4a",
          subplots = 1:2, 
          rel_heights = c(2.5, 1))
ggsave("Actin cytoskeleton organization_2.pdf",height = 4,width = 6)




#------- subGroup  and color -------------------------------------------------------

# set： Group  and color-----------------------
ID= c("GO:0030261","GO:0006325","GO:0006338","GO:0006334",   
     "GO:0034314","GO:0045010","GO:0032956","GO:0030036",
     "GO:0030218")

pvalue_table = TRUE

Selected_GSEA_Color9 = c( "#B71729" ,"#F23545","#FFInner_dataE", # 一组
                          "#0066A5" ,"#3E8CCE","#B1CFF7","#9187D4","#D9E9FF",# 另一组
                          "#FFD6D8")
gseaplot2(Go_gseresult, ID,color=Selected_GSEA_Color9)


gseaplot2(Go_gseresult, ID,color=Selected_GSEA_Color9,
          base_size = 10,
          rel_heights = c(2.5, 0.5, 0.5),
          pvalue_table = TRUE)
  

ID = "GO:0030261"
gseaplot2(Go_gseresult,"GO:0030261",
          title = "",
         
          base_size = 14,
          rel_heights = c(1, 0.2, 0.4),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line")

library(GseaVis)
gseaNb(Go_gseresult,"GO:0030261",
         # addGene = mygene,
          addPval = T,
          pvalX = 0.75,pvalY = 0.8,
          pCol = 'black',
          pHjust = 0)


ggsave(paste("i",".pdf"),device = cairo_pdf,width = 10, height = 6)





##==========================================

Go_gseresult@result$ID

terms <- gsea_res@result$ID[1:3]

gseaNb(object = Go_gseresult,
       #geneSetID = c('GO:0006325','GO:0006335','GO:0006334','GO:0030261'),
       geneSetID = c("chromatin organization","chromatin remodeling"),
       subPlot = 2,
       termWidth = 35)
       #legend.position = c(0.8,0.8),
       #addPval = T,
       #pvalX = 0.05,pvalY = 0.1)

gseaNb(object = Go_gseresult,
       #geneSetID = c('GO:0006325','GO:0006335','GO:0006334','GO:0030261'),
       geneSetID = c("chromatin organization","chromatin remodeling",
                     "nucleosome assembly","chromosome condensation"))




# msigdbr  analysis -----------------------------------------------------------
library(msigdbr)

#  GO--C5
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
gsea_res <- GSEA(geneList, 
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 seed = 123)

write.table (gsea_res, file ="msigdbr_C5_Go_gseresult.csv", sep =",", 
             row.names =TRUE)



# C2——KEGG
KEGG_df = msigdbr(species = "Homo sapiens",category = "C2",
                  subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, entrez_gene)
head(KEGG_df)

KGG_gsea_res <- GSEA(geneList, 
                 TERM2GENE = KEGG_df,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 seed = 123)

write.table (KGG_gsea_res, file ="msigdbr_C2_KEGG_gseresult.csv", sep =",", 
             row.names =TRUE)


# C2——Reactome
REACTOME_df = msigdbr(species = "Homo sapiens",category = "C2",
                  subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, entrez_gene)
head(REACTOME_df)

REACTOME_gsea_res <- GSEA(geneList, 
                     TERM2GENE = REACTOME_df,
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     seed = 123)

write.table (REACTOME_gsea_res, file ="msigdbr_C2_REACTOME_gseresult.csv", 
             sep =",", row.names =TRUE)



# draw
gseaNb(object = gsea_res,
       geneSetID = "GOBP_CYTOPLASMIC_TRANSLATION",
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       addPval = T,
       pvalX = 0.1,pvalY = 0.1)

   

#-------------------------------------------------------------------------------    