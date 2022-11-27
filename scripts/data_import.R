library(tidyverse)
library(readxl)

#read from modified xlsx
Kinase_DUSP_expression_LPS_logFC <- read_excel("Supplementary_tables_S1_to_S12_pure_data.xlsx", sheet = "S10_KinDUSP_expression_LPS_stim")

Kinase_DUSP_LPS_DE <- read_excel("Supplementary_tables_S1_to_S12_pure_data.xlsx", sheet = "S11_All_proteins_regulated_LPS")

Kinase_DUSP_expression_LPS_logFC_list <- split(t(Kinase_DUSP_expression_LPS_logFC[,-c(1:2)]), seq(ncol(Kinase_DUSP_expression_LPS_logFC[,-c(1:2)])))
Kinase_DUSP_expression_LPS_logFC_list <- setNames(Kinase_DUSP_expression_LPS_logFC_list, colnames(Kinase_DUSP_expression_LPS_logFC[,-c(1:2)])) 

Kinase_DUSP_expression_LPS_logFC_list <- lapply(Kinase_DUSP_expression_LPS_logFC_list, function(x){names(x) <- Kinase_DUSP_expression_LPS_logFC$Homolog; return(x)})

Kinase_DUSP_expression_LPS_logFC_list <- lapply(Kinase_DUSP_expression_LPS_logFC_list, function(x){sort(x,decreasing = T)})
Kinase_DUSP_expression_LPS_logFC_list <- lapply(Kinase_DUSP_expression_LPS_logFC_list, na.omit)

Kinase_DUSP_LPS_DE_list <- list()
Kinase_DUSP_LPS_DE_list$Upregulated_in_Monocytes <-  Kinase_DUSP_LPS_DE$Upregulated_in_Monocytes
Kinase_DUSP_LPS_DE_list$Downregulated_in_Monocytes <-  Kinase_DUSP_LPS_DE$Downregulated_in_Monocytes
Kinase_DUSP_LPS_DE_list$Upregulated_in_Dendritic_cells <- Kinase_DUSP_LPS_DE$Upregulated_in_Dendritic_cells
Kinase_DUSP_LPS_DE_list$Downregulated_in_Dendritic_cells <- Kinase_DUSP_LPS_DE$Downregulated_in_Dendritic_cells


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#Kinase_DUSP_LPS_DE_list <- split(t(Kinase_DUSP_LPS_DE), seq(ncol(Kinase_DUSP_LPS_DE)))
Kinase_DUSP_LPS_DE_list <- lapply(Kinase_DUSP_LPS_DE_list, na.omit)
Kinase_DUSP_LPS_DE_list <- lapply(Kinase_DUSP_LPS_DE_list, function(x) sapply(x,tolower))
Kinase_DUSP_LPS_DE_list <- lapply(Kinase_DUSP_LPS_DE_list, function(x) sapply(x,firstup))
Kinase_DUSP_LPS_DE_list <- lapply(Kinase_DUSP_LPS_DE_list, function(x) AnnotationDbi::mapIds(org.Mm.eg.db,keys=x,column = "ENTREZID",keytype = "SYMBOL",multiVals = "first"))
Kinase_DUSP_LPS_DE_list <- lapply(Kinase_DUSP_LPS_DE_list, na.omit)

Kinase_DUSP_LPS_DE_rev <- read_csv("DUSP_manuscript_differentialy_expressed_rev.csv")
Kinase_DUSP_LPS_DE_list_rev_Hs <- list()
Kinase_DUSP_LPS_DE_list_rev_Hs$Upregulated_in_human_Dendritic_cells <-  Kinase_DUSP_LPS_DE_rev$Upregulated_in_human_Dendritic_cells
Kinase_DUSP_LPS_DE_list_rev_Hs$Downregulated_in_human_Dendritic_cells <-  Kinase_DUSP_LPS_DE_rev$Downregulated_in_human_Dendritic_cells
Kinase_DUSP_LPS_DE_list_rev_Hs$Upregulated_in_human_Monocytes <- Kinase_DUSP_LPS_DE_rev$Upregulated_in_human_Monocytes
Kinase_DUSP_LPS_DE_list_rev_Hs$Downregulated_in_human_Monocytes <- Kinase_DUSP_LPS_DE_rev$Downregulated_in_human_Monocytes

Kinase_DUSP_LPS_DE_list_rev_Mm <- list()
Kinase_DUSP_LPS_DE_list_rev_Mm$Upregulated_in_murine_Dendritic_cells <-  Kinase_DUSP_LPS_DE_rev$Upregulated_in_murine_Dendritic_cells
Kinase_DUSP_LPS_DE_list_rev_Mm$Downregulated_in_murine_Dendritic_cells <-  Kinase_DUSP_LPS_DE_rev$Downregulated_in_murine_Dendritic_cells

Kinase_DUSP_LPS_DE_list_rev_Hs <- lapply(Kinase_DUSP_LPS_DE_list_rev_Hs, na.omit)
Kinase_DUSP_LPS_DE_list_rev_Hs <- lapply(Kinase_DUSP_LPS_DE_list_rev_Hs, function(x) AnnotationDbi::mapIds(org.Hs.eg.db,keys=x,column = "ENTREZID",keytype = "SYMBOL",multiVals = "first"))

Kinase_DUSP_LPS_DE_list_rev_Mm <- lapply(Kinase_DUSP_LPS_DE_list_rev_Mm, na.omit)
Kinase_DUSP_LPS_DE_list_rev_Mm <- lapply(Kinase_DUSP_LPS_DE_list_rev_Mm, function(x) AnnotationDbi::mapIds(org.Mm.eg.db,keys=x,column = "ENTREZID",keytype = "SYMBOL",multiVals = "first"))
