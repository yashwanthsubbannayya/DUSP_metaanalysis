
RA_df_curated <- read_csv("RA_df_curated.csv")

RA_df_curated <- RA_df_curated %>%
  filter(p.adjust<= 0.0075)

RA_df_curated$Description

RA_df_curated$Cluster <- factor(RA_df_curated$Cluster,levels = c("Upregulated_in_human_Monocytes","Downregulated_in_human_Monocytes","Upregulated_in_human_Dendritic_cells","Downregulated_in_human_Dendritic_cells","Upregulated_in_murine_Dendritic_cells","Downregulated_in_murine_Dendritic_cells"))

#RA_df_curated <- RA_df_curated %>%
#  group_by(Cluster)%>%
#  arrange(Cluster,p.adjust,Description)

RA_df_curated[,-c(2,4:6,8:11)] %>% 
  spread(Description,p.adjust) -> RA_matrix_curated_proto


RA_matrix_curated <- as.matrix(RA_matrix_curated_proto[,-1])
row.names(RA_matrix_curated) <- RA_matrix_curated_proto$Cluster

RA_matrix_curated[is.na(RA_matrix_curated)] <- 0

RA_matrix_curated_clust <- hclust(dist(t(RA_matrix_curated),method="manhattan"))

RA_matrix_curated_clust$order

RA_df_curated$Description <- factor(RA_df_curated$Description,levels = rev(colnames(RA_matrix_curated)[RA_matrix_curated_clust$order]))


plot <- ggplot(RA_df_curated,aes(y=Description,x=Cluster))+
  geom_point(aes(fill=log(p.adjust),size=GeneRatio),pch=21,colour="black")+
  scale_fill_gradient(low = "#D32726FF", high = "#040494FF", guide = "colourbar")+
  scale_size(range=c(1,10))+
  #scale_x_discrete()+
  #scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, " " , " "),width = 35))+
  scale_x_discrete(breaks=levels(RA_df_curated$Cluster),labels=c("Up in human\n monocytes","Down in human\n monocytes","Up in human\n dendritic cells","Down in human\n dendritic cells","Up in murine\n dendritic cells","Down in murine\n dendritic cells"))+
  labs(x = "DUSP and Kinase mediated\n regulation upon LPS stimulation", y= "Reactome Pathway")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text.y = element_text(size=8),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        panel.grid.major = element_line(colour = 'grey20', size=0.1, linetype='solid'),
        panel.grid.minor = element_line(colour = 'grey20', size=0.1, linetype='solid'),
        panel.border = element_blank(),
        panel.background = element_blank())+
  facet_wrap(.~Organism,scales = "free_x",ncol = 2)

print(plot)
ggsave("Hs_Mm_DUSP_KIN_RA.svg", plot, width = 9, height = 12)
ggsave("Hs_Mm_DUSP_KIN_RA.pdf", plot, width = 9, height = 12)

pdf("Hs_Mm_DUSP_KIN_RA_2.pdf", width = 9, height = 12)
print(plot)
dev.off()
