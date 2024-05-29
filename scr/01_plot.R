# AIM ---------------------------------------------------------------------
# generate some summary plots

# library -----------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(msigdbr)
library(GSEABase)
library(patchwork)
library(ggrepel)

# read in the dataset -----------------------------------------------------
# read in the DE table result
# remember to correct for the estimate of the log2FC in the dataset.
df <- read_csv("../../data/Results080324.csv") |> 
  # dplyr::select(gene = `Gene names`,log2FoldChange) |> 
  # mutate(log2FoldChange = as.numeric(log2FoldChange))
  dplyr::select(gene = `Gene name`,statistic,mean_CTR,mean_DSS,padj = `padj (FDR)`,p.value) |> 
  mutate(log2FC = mean_DSS - mean_CTR) |> 
  # make it negative in this case to simulate the FC
  mutate(neg_statistic = -statistic)

# volcano plot ------------------------------------------------------------
# add the to label the gene name
plot_volcano <- df %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FC)>1&!is.na(gene)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

# save in a table the shortlist version of the genes
df_sig <- plot_volcano %>%
  filter(col == 1) %>%
  arrange(desc(log2FC)) %>%
  mutate(direction = case_when(log2FC>0~"UP",
                               T~"DOWN"))
df_sig %>% 
  write_tsv("../../out/table/table_significant_prot.tsv")

# generate the volcano plot with labels for the features
plot_volcano %>%
  ggplot(aes(x=log2FC,y=-log(padj),label = gene))+
  # geom_point()
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FC,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FC,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none") +
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
  # label all the significant features
  geom_text_repel(data = df_sig,size = 2,box.padding = 0.5,segment.alpha = 0.1,max.overlaps = 20,force_pull = 100)
ggsave("../../out/plot/vulcano_plot_wLabels.pdf",width = 10,height = 10)

# plot a version of the volcano without labels
plot_volcano %>%
  ggplot(aes(x=log2FC,y=-log(padj),label = gene))+
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FC,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FC,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none") +
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/plot/vulcano_plot_woLabels.pdf",width = 10,height = 10)
