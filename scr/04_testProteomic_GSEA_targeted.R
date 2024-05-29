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

# build the signature -----------------------------------------------------
# format in order to be accepted by GSEA
pathways <- readRDS("../../out/object/set_240528.rds")
head(pathways)

# how many genes are in the dataset
lapply(pathways,function(x){
  sum(x %in% names(list_ranks$prot_DSS_vs_CTRL))
})

# Run GSEA ----------------------------------------------------------------
df_tables_GSEA_all <- lapply(list_ranks, function(x){
  fgsea(pathways, x, minSize=1, maxSize=500)  
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
  write_tsv("../../out/table/prot_DSS_vs_CTRL_set240528.tsv")

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
ggsave("../../out/plot/GSEA_res_prot_DSS_vs_CTRL_set240528.pdf",width = 15,height = 10)

# in this case drop the significance threshold
df_tables_GSEA_all %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "REACTOME_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)
ggsave("../../out/plot/GSEA_res_prot_DSS_vs_CTRL_set240528_02.pdf",width = 15,height = 10)

# plot the profile of the GSEA
shortlist_pathway <- unique(df_tables_GSEA_all$pathway)

pdf(paste0("../../out/plot/plotProfiles_res_prot_DSS_vs_CTRL_set240528.pdf"),width = 6,height = 4) 
map(shortlist_pathway,function(name_sig){
  test <- plotEnrichment(pathways[[name_sig]], list_ranks$prot_DSS_vs_CTRL) + labs(title = name_sig) 
  test
})
dev.off()
