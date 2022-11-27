library(clusterProfiler)
library(ReactomePA)

GO_cluster_Mm <- compareCluster(geneCluster = Kinase_DUSP_LPS_DE_list_rev_Mm,
                               fun = "enrichGO",
                               #keyType = "SYMBOL",
                               #universe      = as.character(unique(flatten(Kinase_DUSP_LPS_DE_list))),
                               OrgDb         = "org.Mm.eg.db",
                               ont           = "BP",
                               #pool = TRUE,
                               minGSSize = 10,
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 1,
                               qvalueCutoff  = 1,
                               readable      = TRUE)

KEGG_cluster_Mm <- compareCluster(geneCluster = Kinase_DUSP_LPS_DE_list_rev_Mm,
                             fun = "enrichKEGG",
                             #keyType = "SYMBOL",
                             #universe      = as.character(unique(flatten(Kinase_DUSP_LPS_DE_list))),
                             organism         = "mmu",
                             minGSSize = 10,
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 1,
                             qvalueCutoff  = 1)

MKEGG_cluster_Mm <- compareCluster(geneCluster = Kinase_DUSP_LPS_DE_list_rev_Mm,
                               fun = "enrichMKEGG",
                               #keyType = "SYMBOL",
                               #universe      = as.character(unique(flatten(Kinase_DUSP_LPS_DE_list))),
                               organism         = "mmu",
                               minGSSize = 10,
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 1,
                               qvalueCutoff  = 1)

RA_cluster_Mm <- compareCluster(geneCluster = Kinase_DUSP_LPS_DE_list_rev_Mm,
                                fun = "enrichPathway",
                                #keyType = "SYMBOL",
                                #universe      = as.character(unique(flatten(Kinase_DUSP_LPS_DE_list))),
                                organism         = "mouse",
                                minGSSize = 10,
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)



dotplot(GO_cluster_Mm,showCategory=30)
dotplot(KEGG_cluster_Mm,showCategory=30)
dotplot(MKEGG_cluster_Mm,showCategory=30)
dotplot(RA_cluster_Mm,showCategory=30)

GO_cluster_Hs <- compareCluster(geneCluster = Kinase_DUSP_LPS_DE_list_rev_Hs,
                                fun = "enrichGO",
                                #keyType = "SYMBOL",
                                #universe      = as.character(unique(flatten(Kinase_DUSP_LPS_DE_list))),
                                OrgDb         = "org.Hs.eg.db",
                                ont           = "BP",
                                #pool = TRUE,
                                minGSSize = 10,
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)

KEGG_cluster_Hs <- compareCluster(geneCluster = Kinase_DUSP_LPS_DE_list_rev_Hs,
                                  fun = "enrichKEGG",
                                  #keyType = "SYMBOL",
                                  #universe      = as.character(unique(flatten(Kinase_DUSP_LPS_DE_list))),
                                  organism         = "hsa",
                                  minGSSize = 10,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1)

MKEGG_cluster_Hs <- compareCluster(geneCluster = Kinase_DUSP_LPS_DE_list_rev_Hs,
                                   fun = "enrichMKEGG",
                                   #keyType = "SYMBOL",
                                   #universe      = as.character(unique(flatten(Kinase_DUSP_LPS_DE_list))),
                                   organism         = "hsa",
                                   minGSSize = 10,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 1,
                                   qvalueCutoff  = 1)

RA_cluster_Hs <- compareCluster(geneCluster = Kinase_DUSP_LPS_DE_list_rev_Hs,
                                fun = "enrichPathway",
                                #keyType = "SYMBOL",
                                #universe      = as.character(unique(flatten(Kinase_DUSP_LPS_DE_list))),
                                organism         = "human",
                                minGSSize = 10,
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 1,
                                qvalueCutoff  = 1,
                                readable      = TRUE)



dotplot(GO_cluster_Hs,showCategory=30)
dotplot(KEGG_cluster_Hs,showCategory=30)
dotplot(MKEGG_cluster_Hs,showCategory=30)
dotplot(RA_cluster_Hs,showCategory=30)

RA_df_Hs <- as.data.frame(RA_cluster_Hs)
RA_df_Hs$Organism <- "Human"
RA_df_Mm <- as.data.frame(RA_cluster_Mm)
RA_df_Mm$Organism <- "Mouse"

RA_df <- rbind(RA_df_Hs,RA_df_Mm)

eval_col <- function(x){
  eval(parse(text=as.character(x)))
}

RA_df$GeneRatio <- sapply(RA_df$GeneRatio,eval_col)
RA_df$BgRatio <- sapply(RA_df$BgRatio,eval_col)

write_csv(RA_df,"RA_df.csv")

