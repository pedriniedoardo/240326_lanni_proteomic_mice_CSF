# AIM ---------------------------------------------------------------------
# run unbiased GSEA on the full proteomic dataset REACTOME annotation

# library -----------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(msigdbr)
library(GSEABase)
library(patchwork)
library(ggrepel)

# read in the data --------------------------------------------------------
# read in the DE table result
# remember to correct for the estimate of the log2FC in the dataset.
df <- read_csv("../../data/Results080324.csv") |> 
  # dplyr::select(gene = `Gene names`,log2FoldChange) |> 
  # mutate(log2FoldChange = as.numeric(log2FoldChange))
  dplyr::select(gene = `Gene name`,statistic,mean_CTR,mean_DSS,padj = `padj (FDR)`,p.value) |> 
  mutate(log2FC = mean_DSS - mean_CTR) |> 
  # make it negative in this case to simulate the FC
  mutate(neg_statistic = -statistic)

# wrangling ---------------------------------------------------------------
# clean the gene names and pull the ranking metric
test <- df |>
  split(f = df$gene) |> 
  lapply(function(x){
    gene <- x$gene |>
      str_split(pattern = ";") |> 
      unlist()
    data.frame(gene = gene,log2FC = x$log2FC)
  }) |> 
  bind_rows() |> 
  dplyr::filter(!is.na(log2FC))


# load the object in a list, split it by comparison
results <- list(prot_DSS_vs_CTRL = test)

# build the rank ----------------------------------------------------------
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


# pick the signature ------------------------------------------------------
# library("msigdbr")
msigdbr_collections() %>%
  print(n=30)

#
reactome_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")
head(reactome_gene_sets)

# format in order to be accepted by GSEA
pathways <- split(x = reactome_gene_sets$gene_symbol, f = reactome_gene_sets$gs_name)
# head(pathways)

# Run GSEA ----------------------------------------------------------------
df_tables_GSEA_all <- lapply(list_ranks, function(x){
  fgsea(pathways, x, minSize=10, maxSize=500)  
}) %>%
  bind_rows(.id = "dataset") %>% 
  # the ladingEdge columns has to be re-arranged in order to save the file as a table (originally is a list) 
  mutate(leadingEdge = unlist(lapply(.$leadingEdge, function(x){
    paste0(x,collapse = "|")
  }))) %>%
  arrange(padj,-abs(NES)) 

# check the result
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

# Collapse Redundant Terms ------------------------------------------------
# split the dataset per type
list_tables_GSEA_all <- split(df_tables_GSEA_all,f = df_tables_GSEA_all$dataset)

# check the comparisions elementes
names(list_tables_GSEA_all)

# collapse the pathways
list_collapsedPathways <- lapply(names(list_tables_GSEA_all),function(x){
  collapsePathways(list_tables_GSEA_all[[x]], pathways, list_ranks[[x]])
}) %>%
  setNames(names(list_tables_GSEA_all))

# check the result
str(list_collapsedPathways)

# pull the main pathways
list_mainPathways <- pmap(list(list_tables_GSEA_all,list_collapsedPathways),function(x,y){
  x %>%
    dplyr::filter(pathway %in% y$mainPathways) %>%
    arrange(padj,-abs(NES)) %>%
    pull(pathway) 
})

str(list_mainPathways)

# save list of non redundant terms
# check the order of the names is the same
sum(!names(list_tables_GSEA_all) == names(list_mainPathways))

# filter only the non redundant for each comparison
df_tables_GSEA_all_non_redundant <- 
  pmap(list(list_tables_GSEA_all,list_mainPathways),function(x,y){
    x %>%
      dplyr::filter(pathway %in% y)
  }) %>%
  bind_rows()

# save the table
df_tables_GSEA_all_non_redundant %>%
  write_tsv(file = "../../out/table/prot_DSS_vs_CTRL_nonredundant_REACTOME.tsv")

# plot GSEA results -------------------------------------------------------
# full GSEA table result
# library(ggrepel)
df_tables_GSEA_all %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "REACTOME_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = (-log10(0.05)),linetype="dashed",col="gray")
ggsave("../../out/plot/GSEA_res_prot_DSS_vs_CTRL_all_REACTOME.pdf",width = 15,height = 10)

# non redundant terms
# library(ggrepel)
df_tables_GSEA_all_non_redundant %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "REACTOME_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = (-log10(0.05)),linetype="dashed",col="gray")
ggsave("../../out/plot/GSEA_res_prot_DSS_vs_CTRL_nonredundant_REACTOME.pdf",width = 14,height = 10)
