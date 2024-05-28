# libraries ---------------------------------------------------------------
library(tidyverse)
library(enrichR)
library(scales)
library(patchwork)

# run enrichr with the list of genes in the module
# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
# filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Cell"))
dbs %>%
  filter(str_detect(libraryName,pattern = "Human"))

dbs %>%
  filter(str_detect(libraryName,pattern = "MSigDB"))
dbs %>%
  filter(str_detect(libraryName,pattern = "SigDB"))
dbs %>%
  filter(str_detect(libraryName,pattern = "KEGG"))
dbs %>%
  filter(str_detect(libraryName,pattern = "Reactome"))
dbs %>%
  filter(str_detect(libraryName,pattern = "GO"))


dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2022","WikiPathway_2023_Human","GO_Biological_Process_2023")

# query -------------------------------------------------------------------
test_plot <- read_csv("../../data/Results080324.csv") |> 
  # dplyr::select(gene = `Gene names`,log2FoldChange) |> 
  # mutate(log2FoldChange = as.numeric(log2FoldChange))
  dplyr::select(gene = `Gene name`,statistic,mean_CTR,mean_DSS,padj = `padj (FDR)`,p.value) |> 
  mutate(log2FC = mean_DSS - mean_CTR) |> 
  # make it negative in this case to simulate the FC
  mutate(neg_statistic = -statistic)

list_genes_UP <- test_plot %>%
  mutate(test = "DSS_vs_CTRL") %>%
  filter(padj<0.05 & log2FC > 1&!is.na(gene)) %>%
  split(f = .$test) %>%
  map(function(x){
    x %>%
      pull(gene)
  })

list_genes_DOWN <- test_plot %>%
  mutate(test = "DSS_vs_CTRL") %>%
  filter(padj<0.05 & log2FC < -1&!is.na(gene)) %>%
  split(f = .$test) %>%
  map(function(x){
    x %>%
      pull(gene)
  })

# define the background
background <- test_plot$gene

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_enrichr_UP <- lapply(list_genes_UP,function(x){
  genes <- x
  out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  # out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>%
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_enrichr_UP %>%
  write_tsv("../../out/table/enrichR_DE_pseudobulk_UP_background.tsv")

list_enrichr_DOWN <- lapply(list_genes_DOWN,function(x){
  genes <- x
  out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  # out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>%
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_enrichr_DOWN %>%
  write_tsv("../../out/table/enrichR_DE_pseudobulk_DOWN_background.tsv")

# list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")
plot_list_UP <- list_enrichr_UP %>%
  split(f = .$comparison)

plot_list_DOWN <- list_enrichr_DOWN %>%
  split(f = .$comparison)

# library(scales)
list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_UP,nrow = 1)
ggsave("../../out/plot/enrichR_DE_pseudobulk_UP_background.pdf",width = 6,height = 12,limitsize = FALSE)

list_plot_DOWN <- pmap(list(plot_list_DOWN,names(plot_list_DOWN)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_DOWN,nrow = 1)
ggsave("../../out/plot/enrichR_DE_pseudobulk_DOWN_background.pdf",width = 6,height = 12,limitsize = FALSE)

