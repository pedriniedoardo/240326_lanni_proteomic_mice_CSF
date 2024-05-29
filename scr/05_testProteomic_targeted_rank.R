# AIM ---------------------------------------------------------------------
# plot the ranks of the tailored signatures

# library -----------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(msigdbr)
library(GSEABase)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)

# pull the DE results -----------------------------------------------------
# read in the DE table result
# remember to correct for the estimate of the log2FC in the dataset.
df <- read_csv("../../data/Results080324.csv") |> 
  # dplyr::select(gene = `Gene names`,log2FoldChange) |> 
  # mutate(log2FoldChange = as.numeric(log2FoldChange))
  dplyr::select(gene = `Gene name`,statistic,mean_CTR,mean_DSS,padj = `padj (FDR)`,p.value) |> 
  mutate(log2FC = mean_DSS - mean_CTR) |> 
  # make it negative in this case to simulate the FC
  mutate(neg_statistic = -statistic)

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

# pull the pathways reference ---------------------------------------------
# format in order to be accepted by GSEA
pathways <- readRDS("../../out/object/set_240528.rds")
head(pathways)

# how many genes are in the dataset
lapply(pathways,function(x){
  sum(x %in% names(list_ranks$prot_DSS_vs_CTRL))
})

# read in GSEA results ----------------------------------------------------
res_GSEA <- read_tsv("../../out/table/prot_DSS_vs_CTRL_set240528.tsv")

# subset the leading edges
df_leading_edges <- res_GSEA %>%
  pull(leadingEdge) %>%
  str_split(pattern = "\\|") %>%
  setNames(res_GSEA$pathway)

# plot heatmap with rank of genes -----------------------------------------
# define a shortlist of patways to focus on
shortlist_pathway <- res_GSEA$pathway

lapply(shortlist_pathway,function(path){

  # get the top and the bottom genes of the signatures in the rank
  df_rank <- results$prot_DSS_vs_CTRL %>%
    arrange(desc(log2FC)) %>%
    mutate(rank = nrow(.):1) %>%
    mutate(rank2 = 1:nrow(.))
  
  # filter the ranks of the gene in the signature
  id_UP <- df_rank %>%
    filter(gene %in% pathways[[path]]) %>%
    mutate(color = case_when(gene %in% df_leading_edges[[path]]~"red",T~"black"))
  # filter(color != "black")
  
  head(id_UP)
  
  # library(ComplexHeatmap)
  m <- matrix(df_rank$rank,ncol = 1)
  ha_up <- rowAnnotation(foo = anno_mark(at = id_UP$rank2, labels = id_UP$gene,
                                         labels_gp = gpar(col = id_UP$color,
                                                          fontsize = 10),
                                         link_width = unit(20, "mm"),
                                         extend = unit(100, "mm")))
  
  hm_up <- Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha_up)
  
  pdf(paste0("../../out/plot/rank_",path,".pdf"),width = 3,height = 15)
  draw(hm_up,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(100, 2, 2, 2), "mm"))
  dev.off()
})
