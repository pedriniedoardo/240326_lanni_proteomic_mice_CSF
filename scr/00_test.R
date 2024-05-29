df_long <- read_csv("../../data/Results080324.csv") %>%
  dplyr::select(gene = `Gene name`,CTR_10,CTR_3,CTR_5,CTR_6,CTR_7,CTR_8,CTR_9,DSS_10,DSS_11,DSS_4,DSS_5,DSS_6,DSS_7,DSS_8,DSS_9) %>%
  pivot_longer(names_to = "sample",values_to = "exp",-gene) %>%
  separate(sample,into = c("treat","rep"),remove = F,sep = "_")

# try the test
test <- df_long %>%
  filter(gene %in% c("SAA2"))

test %>%
  ggplot(aes(x=treat,y=exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1))+
  facet_wrap(~gene)+
  theme_bw()+
  theme(strip.background = element_blank())

lm(test,formula = exp~treat) %>%
  summary()

test %>%
  group_by(treat) %>%
  summarise(avg = mean(exp)) %>%
  pivot_wider(names_from = treat,values_from = avg) %>%
  mutate(logFC = DSS-CTR) %>%
  as.data.frame()

df_long2 <- read_csv("../../data/Results080324.csv") %>%
  dplyr::select(gene = `Gene name`,CTR = mean_CTR,DSS = mean_DSS) %>%
  pivot_longer(names_to = "treat",values_to = "avg_ref",-gene)

df_long %>%
  group_by(gene,treat) %>%
  summarise(avg = mean(exp)) %>%
  left_join(df_long2,by=c("gene","treat")) %>%
  # test if they are different
  mutate(delta = avg-avg_ref) %>%
  arrange(delta)


df_test <- read_csv("../../data/Results080324.csv") %>%
  dplyr::select(gene = `Gene name`,mean_CTR,mean_DSS,log2FoldChange_ref = log2FoldChange) %>%
  # try to replicate
  mutate(test = log2(mean_DSS/mean_CTR)) %>%
  mutate(delta = log2FoldChange_ref - test)

df_test %>%
  ggplot(aes(x=delta))+geom_histogram()


library(enrichR)
?enrichr()

# -------------------------------------------------------------------------
x <- plot_list_DOWN$DSS_vs_CTRL
y <- "DSS_vs_CTRL"
list_plot_DOWN <- pmap(list(plot_list_DOWN,names(plot_list_DOWN)), function(x,y){
  
  x %>%
    # dplyr::filter(annotation == "GO_Molecular_Function_2023") %>%
    # dplyr::filter(annotation == "GO_Cellular_Component_2023") %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 40)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(paste(y,"DOWN"))
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

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
reactome_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")
head(reactome_gene_sets)

# format in order to be accepted by GSEA
pathways <- split(x = reactome_gene_sets$gene_symbol, f = reactome_gene_sets$gs_name)

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

df_tables_GSEA_all <- read.table("../../out/table/prot_DSS_vs_CTRL_all_REACTOME.tsv",header = T)

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

pdf(paste0("../../out/plot/plotProfiles_res_prot_DSS_vs_CTRL_REACTOME_INNATE_IMMUNE_SYSTEM.pdf"),width = 6,height = 4) 
plotEnrichment(pathways[["REACTOME_INNATE_IMMUNE_SYSTEM"]], list_ranks$prot_DSS_vs_CTRL) + labs(title = "REACTOME_INNATE_IMMUNE_SYSTEM")
dev.off()

plotEnrichment(pathways[["REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION"]], list_ranks$prot_DSS_vs_CTRL) + labs(title = "REACTOME_ASPARAGINE_N_LINKED_GLYCOSYLATION")

plotEnrichment(pathways[["REACTOME_NEUTROPHIL_DEGRANULATION"]], list_ranks$prot_DSS_vs_CTRL) + labs(title = "REACTOME_NEUTROPHIL_DEGRANULATION")
plotEnrichment(pathways[["REACTOME_TRANSPORT_TO_THE_GOLGI_AND_SUBSEQUENT_MODIFICATION"]], list_ranks$prot_DSS_vs_CTRL) + labs(title = "REACTOME_TRANSPORT_TO_THE_GOLGI_AND_SUBSEQUENT_MODIFICATION")

plotEnrichment(pathways[["REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2"]], list_ranks$prot_DSS_vs_CTRL) + labs(title = "REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2")


plotEnrichment(pathways[["REACTOME_MEMBRANE_TRAFFICKING"]], list_ranks$prot_DSS_vs_CTRL) + labs(title = "REACTOME_MEMBRANE_TRAFFICKING")

# subset the leading edges
df_leading_edges <- df_tables_GSEA_all %>%
  pull(leadingEdge) %>%
  str_split(pattern = "\\|") %>%
  setNames(df_tables_GSEA_all$pathway)

# plot heatmap with rank of genes -----------------------------------------
# define a shortlist of patways to focus on
shortlist_pathway <- c("REACTOME_INNATE_IMMUNE_SYSTEM","REACTOME_MEMBRANE_TRAFFICKING")

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
  
  pdf(paste0("../../out/plot/rank_",path,".pdf"),width = 3,height = 20)
  draw(hm_up,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(100, 2, 2, 2), "mm"))
  dev.off()
})
