# library -----------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(msigdbr)
library(GSEABase)
library(patchwork)
library(ggrepel)

# prepare the dataset with all the annoration needed ---------------------- 
# locate the files 
# file <- dir("out/") %>%
#   str_subset(pattern = "res_") %>%
#   str_subset(pattern = ".txt") %>%
#   str_subset(pattern = "shr",negate = T)
test_plot <- read_csv("../../data/Results080324.csv") |> 
  # dplyr::select(gene = `Gene names`,log2FoldChange) |> 
  # mutate(log2FoldChange = as.numeric(log2FoldChange))
  dplyr::select(gene = `Gene name`,statistic,mean_CTR,mean_DSS,padj = `padj (FDR)`,p.value) |> 
  mutate(log2FC = mean_DSS - mean_CTR) |> 
  # make it negative in this case to simulate the FC
  mutate(neg_statistic = -statistic)
  # mutate(log2FoldChange = as.numeric(log2FoldChange))

# plot volcano
# add the info of the genename
plot_volcano <- test_plot %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FC)>1&!is.na(gene)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

# save in a table to genes only up and only down
plot_volcano

plot_volcano %>%
  filter(col == 1) %>%
  arrange(desc(log2FC)) %>%
  mutate(direction = case_when(log2FC>0~"UP",
                               T~"DOWN")) %>%
  write_tsv("../../out/table/table_significant_prot.tsv")

plot_volcano %>%
  ggplot(aes(x=log2FC,y=-log(padj),label = gene))+
  # geom_point()
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FC,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FC,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  # ggrepel::geom_text_repel(
  #   data = plot_volcano %>% group_by(conditionVsCX) %>% arrange(padj) %>% dplyr::slice(1:10),
  #   aes(label = symbol),segment.alpha=0.4) +
  # facet_wrap(~conditionVsCX)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none") +
  geom_text_repel(data = plot_volcano_label,size = 2,box.padding = 0.5,segment.alpha = 0.1,max.overlaps = 20,force_pull = 100)
ggsave("../../out/plot/vulcano_plot_pseudobulk.pdf",width = 10,height = 10)

plot_volcano %>%
  ggplot(aes(x=log2FC,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FC,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FC,y=-log(padj),col=factor(col)),alpha=0.5)+
  # geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  # ggrepel::geom_text_repel(
  #   data = plot_volcano %>% group_by(conditionVsCX) %>% arrange(padj) %>% dplyr::slice(1:10),
  #   aes(label = symbol),segment.alpha=0.4) +
  # facet_wrap(~conditionVsCX)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")

# read_csv("../../data/Organoids_Absinta_Results.csv") |>
#   dplyr::select(gene = `Gene names`,statistic,log2FoldChange,CTRL_Rep1,CTRL_Rep2,CTRL_Rep3,Treated_Rep1,Treated_Rep2,Treated_Rep3,mean_CTRL,mean_Treated,p.value,padj) |> 
#   mutate(neg_statistic = -statistic) |> 
#   mutate(log2FoldChange = as.numeric(log2FoldChange)) |> 
#   filter(gene %in% c("STAT3","GNB2","HRAS","GNB4","CSK","RAP1B","PRKACB","CRKL","GNB1","AKT1","AKT2","AKT3")) |> 
#   pivot_longer(names_to = "sample",values_to = "exp",CTRL_Rep1:Treated_Rep3) |> 
#   separate(sample,into = c("treat","rep"),remove = F) |>
#   mutate(exp_log = log1p(exp)) |> 
#   ggplot(aes(x=treat,y=exp_log))+facet_wrap(~gene,scales = "free")+
#   geom_boxplot()+
#   geom_point(position = position_jitter(width = 0.1))+
#   theme_bw()+
#   theme(strip.background = element_blank())

test_plot |> 
  ggplot(aes(x=log2FC,y=neg_statistic))+geom_point()

# use the statistic colum as the FC is too sparse
test_df <- read_csv("../../data/Results080324.csv") |>
  # dplyr::select(gene = `Gene names`,log2FoldChange) |> 
  # mutate(log2FoldChange = as.numeric(log2FoldChange))
  dplyr::select(gene = `Gene name`,statistic,mean_CTR,mean_DSS,padj = `padj (FDR)`,p.value) |> 
  mutate(log2FC = mean_DSS - mean_CTR)

test <- test_df |>
  split(f = test_df$gene) |> 
  lapply(function(x){
    gene <- x$gene |>
      str_split(pattern = ";") |> 
      unlist()
    data.frame(gene = gene,log2FC = x$log2FC)
  }) |> 
  bind_rows() |> 
  dplyr::filter(!is.na(log2FC))
  

# load the results 
results <- list(prot_DSS_vs_CTRL = test)

# GSEA -------------------------------------------------------------------- 
# use the FC dataset to create the ranked list of genes 
# Symbol or Entrez? 
# x <- results$res_MutvsWT_shr
list_ranks <- lapply(results, function(x){
  
  x <- dplyr::filter(x,!is.na(gene)) %>%
    # remove mito and ribo genes
    # filter(str_detect(gene,pattern = "^MT-",negate = T)) %>% 
    # filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate = T)) %>%
    # average logFC in case of duplicated genenames
    group_by(gene) %>%
    summarise(logFC = mean(log2FC))
  
  ranks <- setNames(x$logFC, x$gene)
  ranks
}) 
glimpse(list_ranks)

# score all the signatures in MsigDB from C2 category ---------------------
# library("msigdbr")
msigdbr_collections() %>%
  print(n=30)
#
reactome_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
# kegg_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
# get all the C2 terms
# cgp_gene_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
# C2_gene_sets = msigdbr(species = "Mus musculus", category = "C2")
# h_gene_sets <- msigdbr(species = "Mus musculus", category = "H")
head(reactome_gene_sets)

# format in order to be accepted by GSEA
pathways <- split(x = reactome_gene_sets$gene_symbol, f = reactome_gene_sets$gs_name)
# head(pathways)

# RUN GSEA ----------------------------------------------------------------
df_tables_GSEA_all <- lapply(list_ranks, function(x){
  fgsea(pathways, x, minSize=10, maxSize=500)  
}) %>%
  bind_rows(.id = "dataset") %>% 
  # the ladingEdge columns has to be re-arranged in order to save the file as a table (originally is a list) 
  mutate(leadingEdge = unlist(lapply(.$leadingEdge, function(x){
    paste0(x,collapse = "|")
  }))) %>%
  arrange(padj,-abs(NES)) 

dim(df_tables_GSEA_all)

head(df_tables_GSEA_all,n=20) 

# save the whole table
df_tables_GSEA_all %>%
  write_tsv("../../out/table/prot_DSS_vs_CTRL_all_KEGG.tsv")

# explore the correlation between each subset focussing on the ES
# library(GGally)
# test_NES <- df_tables_GSEA_all %>%
#   dplyr::select(dataset,ES,pathway) %>%
#   pivot_wider(names_from = dataset,values_from = ES) %>%
#   column_to_rownames("pathway")
# 
# ggpairs(test_NES)

# COLLAPSE REDUNDANT ------------------------------------------------------
# collapsing the similar pathways 

# split the dataset per type
list_tables_GSEA_all <- split(df_tables_GSEA_all,f = df_tables_GSEA_all$dataset)

names(list_tables_GSEA_all)

list_collapsedPathways <- lapply(names(list_tables_GSEA_all),function(x){
  collapsePathways(list_tables_GSEA_all[[x]], pathways, list_ranks[[x]])
}) %>%
  setNames(names(list_tables_GSEA_all))

# collapsedPathways_MG <- collapsePathways(df_tables_GSEA_all %>%
#                                            filter(dataset %in% "list_df_comparisons_MG_old_young"), pathways, list_ranks$list_df_comparisons_MG_old_young)
# glimpse(collapsedPathways_MG)
str(list_collapsedPathways)

list_mainPathways <- pmap(list(list_tables_GSEA_all,list_collapsedPathways),function(x,y){
  x %>%
    dplyr::filter(pathway %in% y$mainPathways) %>%
    arrange(padj,-abs(NES)) %>%
    pull(pathway) 
})

str(list_mainPathways)

# save list of non redundant terms
# chackt the order of the names is the same
sum(!names(list_tables_GSEA_all) == names(list_mainPathways))

# filter only the non redundant fro each comparison
df_tables_GSEA_all_non_redundant <- 
  pmap(list(list_tables_GSEA_all,list_mainPathways),function(x,y){
    x %>%
      dplyr::filter(pathway %in% y)
  }) %>%
  bind_rows()

# save the table
df_tables_GSEA_all_non_redundant %>%
  write_tsv(file = "../../out/table/prot_DSS_vs_CTRL_nonredundant_KEGG.tsv")

test <- df_tables_GSEA_all_non_redundant %>%
  group_by(dataset) %>%
  top_n(wt = padj*(-1),n = 5)

# test plot to show the main terms in each dataset
library(ggrepel)
df_tables_GSEA_all_non_redundant %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = (-log10(0.05)),linetype="dashed",col="gray")
ggsave("../../out/plot/GSEA_res_prot_DSS_vs_CTRL_nonredundant_KEGG.pdf",width = 14,height = 10)

# library(ggrepel)
df_tables_GSEA_all %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = (-log10(0.05)),linetype="dashed",col="gray")
ggsave("../../out/plot/GSEA_res_prot_DSS_vs_CTRL_all_KEGG.pdf",width = 15,height = 10)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# score all the signatures in MsigDB from C2 category ---------------------
# library("msigdbr")
msigdbr_collections() %>%
  print(n=30)
#
reactome_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")
# kegg_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
# get all the C2 terms
# cgp_gene_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
# C2_gene_sets = msigdbr(species = "Mus musculus", category = "C2")
# h_gene_sets <- msigdbr(species = "Mus musculus", category = "H")
head(reactome_gene_sets)

# format in order to be accepted by GSEA
pathways <- split(x = reactome_gene_sets$gene_symbol, f = reactome_gene_sets$gs_name)
# head(pathways)

# RUN GSEA ----------------------------------------------------------------
df_tables_GSEA_all <- lapply(list_ranks, function(x){
  fgsea(pathways, x, minSize=10, maxSize=500)  
}) %>%
  bind_rows(.id = "dataset") %>% 
  # the ladingEdge columns has to be re-arranged in order to save the file as a table (originally is a list) 
  mutate(leadingEdge = unlist(lapply(.$leadingEdge, function(x){
    paste0(x,collapse = "|")
  }))) %>%
  arrange(padj,-abs(NES)) 

dim(df_tables_GSEA_all)

head(df_tables_GSEA_all,n=20) 

# save the whole table
df_tables_GSEA_all %>%
  write_tsv("../../out/table/prot_DSS_vs_CTRL_all_REACTOME.tsv")

# explore the correlation between each subset focussing on the ES
# library(GGally)
# test_NES <- df_tables_GSEA_all %>%
#   dplyr::select(dataset,ES,pathway) %>%
#   pivot_wider(names_from = dataset,values_from = ES) %>%
#   column_to_rownames("pathway")
# 
# ggpairs(test_NES)

# COLLAPSE REDUNDANT ------------------------------------------------------
# collapsing the similar pathways 

# split the dataset per type
list_tables_GSEA_all <- split(df_tables_GSEA_all,f = df_tables_GSEA_all$dataset)

names(list_tables_GSEA_all)

list_collapsedPathways <- lapply(names(list_tables_GSEA_all),function(x){
  collapsePathways(list_tables_GSEA_all[[x]], pathways, list_ranks[[x]])
}) %>%
  setNames(names(list_tables_GSEA_all))

# collapsedPathways_MG <- collapsePathways(df_tables_GSEA_all %>%
#                                            filter(dataset %in% "list_df_comparisons_MG_old_young"), pathways, list_ranks$list_df_comparisons_MG_old_young)
# glimpse(collapsedPathways_MG)
str(list_collapsedPathways)

list_mainPathways <- pmap(list(list_tables_GSEA_all,list_collapsedPathways),function(x,y){
  x %>%
    dplyr::filter(pathway %in% y$mainPathways) %>%
    arrange(padj,-abs(NES)) %>%
    pull(pathway) 
})

str(list_mainPathways)

# save list of non redundant terms
# chackt the order of the names is the same
sum(!names(list_tables_GSEA_all) == names(list_mainPathways))

# filter only the non redundant fro each comparison
df_tables_GSEA_all_non_redundant <- 
  pmap(list(list_tables_GSEA_all,list_mainPathways),function(x,y){
    x %>%
      dplyr::filter(pathway %in% y)
  }) %>%
  bind_rows()

# save the table
df_tables_GSEA_all_non_redundant %>%
  write_tsv(file = "../../out/table/prot_DSS_vs_CTRL_nonredundant_REACTOME.tsv")

test <- df_tables_GSEA_all_non_redundant %>%
  group_by(dataset) %>%
  top_n(wt = padj*(-1),n = 5)

# test plot to show the main terms in each dataset
library(ggrepel)
df_tables_GSEA_all_non_redundant %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = (-log10(0.05)),linetype="dashed",col="gray")
ggsave("../../out/plot/GSEA_res_prot_DSS_vs_CTRL_nonredundant_REACTOME.pdf",width = 14,height = 10)

# library(ggrepel)
df_tables_GSEA_all %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = (-log10(0.05)),linetype="dashed",col="gray")
ggsave("../../out/plot/GSEA_res_prot_DSS_vs_CTRL_all_REACTOME.pdf",width = 15,height = 10)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# score all the signatures in MsigDB from C2 category ---------------------
# library("msigdbr")
msigdbr_collections() %>%
  print(n=30)
#
reactome_gene_sets <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP")
# kegg_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
# get all the C2 terms
# cgp_gene_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
# C2_gene_sets = msigdbr(species = "Mus musculus", category = "C2")
# h_gene_sets <- msigdbr(species = "Mus musculus", category = "H")
head(reactome_gene_sets)

# format in order to be accepted by GSEA
pathways <- split(x = reactome_gene_sets$gene_symbol, f = reactome_gene_sets$gs_name)
# head(pathways)

# RUN GSEA ----------------------------------------------------------------
df_tables_GSEA_all <- lapply(list_ranks, function(x){
  fgsea(pathways, x, minSize=10, maxSize=500)  
}) %>%
  bind_rows(.id = "dataset") %>% 
  # the ladingEdge columns has to be re-arranged in order to save the file as a table (originally is a list) 
  mutate(leadingEdge = unlist(lapply(.$leadingEdge, function(x){
    paste0(x,collapse = "|")
  }))) %>%
  arrange(padj,-abs(NES)) 

dim(df_tables_GSEA_all)

head(df_tables_GSEA_all,n=20) 

# save the whole table
df_tables_GSEA_all %>%
  write_tsv("../../out/table/prot_DSS_vs_CTRL_all_GOBP.tsv")

# explore the correlation between each subset focussing on the ES
# library(GGally)
# test_NES <- df_tables_GSEA_all %>%
#   dplyr::select(dataset,ES,pathway) %>%
#   pivot_wider(names_from = dataset,values_from = ES) %>%
#   column_to_rownames("pathway")
# 
# ggpairs(test_NES)

# COLLAPSE REDUNDANT ------------------------------------------------------
# collapsing the similar pathways 

# split the dataset per type
list_tables_GSEA_all <- split(df_tables_GSEA_all,f = df_tables_GSEA_all$dataset)

names(list_tables_GSEA_all)

list_collapsedPathways <- lapply(names(list_tables_GSEA_all),function(x){
  collapsePathways(list_tables_GSEA_all[[x]], pathways, list_ranks[[x]])
}) %>%
  setNames(names(list_tables_GSEA_all))

# collapsedPathways_MG <- collapsePathways(df_tables_GSEA_all %>%
#                                            filter(dataset %in% "list_df_comparisons_MG_old_young"), pathways, list_ranks$list_df_comparisons_MG_old_young)
# glimpse(collapsedPathways_MG)
str(list_collapsedPathways)

list_mainPathways <- pmap(list(list_tables_GSEA_all,list_collapsedPathways),function(x,y){
  x %>%
    dplyr::filter(pathway %in% y$mainPathways) %>%
    arrange(padj,-abs(NES)) %>%
    pull(pathway) 
})

str(list_mainPathways)

# save list of non redundant terms
# chackt the order of the names is the same
sum(!names(list_tables_GSEA_all) == names(list_mainPathways))

# filter only the non redundant fro each comparison
df_tables_GSEA_all_non_redundant <- 
  pmap(list(list_tables_GSEA_all,list_mainPathways),function(x,y){
    x %>%
      dplyr::filter(pathway %in% y)
  }) %>%
  bind_rows()

# save the table
df_tables_GSEA_all_non_redundant %>%
  write_tsv(file = "../../out/table/prot_DSS_vs_CTRL_nonredundant_GOBP.tsv")

test <- df_tables_GSEA_all_non_redundant %>%
  group_by(dataset) %>%
  top_n(wt = padj*(-1),n = 5)

# test plot to show the main terms in each dataset
library(ggrepel)
df_tables_GSEA_all_non_redundant %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = (-log10(0.05)),linetype="dashed",col="gray")
ggsave("../../out/plot/GSEA_res_prot_DSS_vs_CTRL_nonredundant_GOBP.pdf",width = 14,height = 10)

# library(ggrepel)
df_tables_GSEA_all %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = (-log10(0.05)),linetype="dashed",col="gray")
ggsave("../../out/plot/GSEA_res_prot_DSS_vs_CTRL_all_GOBP.pdf",width = 15,height = 10)
