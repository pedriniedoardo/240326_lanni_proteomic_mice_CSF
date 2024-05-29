# AIM ---------------------------------------------------------------------
# generate UpSet from the tailored list of annotation

# libraries ---------------------------------------------------------------
library(tidyverse)
library(UpSetR) 

# read in the data --------------------------------------------------------
df <- read_tsv("../../data/Lanni_customSet.txt")

# wrangling ---------------------------------------------------------------
# divided the dataset into list of UP and DOWN genes and generate a list of genes
list_DE_up <- df %>%
  filter(direction == "UP") %>%
  split(f = .$set) %>%
  lapply(function(x){
    x %>%
      pull(gene) %>%
      unique()
  })

list_DE_down <- df %>%
  filter(direction == "DOWN") %>%
  split(f = .$set) %>%
  lapply(function(x){
    x %>%
      pull(gene) %>%
      unique()
  })

# plot --------------------------------------------------------------------
pdf("../../out/plot/upset_DSS_vs_CTRL_UP.pdf",width = 7,height = 6,onefile=FALSE)
upset(fromList(list_DE_up), order.by = "freq",nintersects = NA,nsets = 10)
dev.off()

pdf("../../out/plot/upset_DSS_vs_CTRL_DOWN.pdf",width = 7,height = 6,onefile=FALSE)
upset(fromList(list_DE_down), order.by = "freq",nintersects = NA,nsets = 10)
dev.off()

# extract intersections ---------------------------------------------------
# for UP list
df1_up <- lapply(list_DE_up,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

head(df1_up)

df2_up <- data.frame(gene=unique(unlist(list_DE_up)))
head(df2_up)

# now loop through each individual gene and pick the list of all the intersections they belong to
df_int_up <- lapply(df2_up$gene,function(x){
  # pull the name of the intersections
  intersection <- df1_up %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

head(df_int_up,n=20)

df_int_up %>%
  write_tsv("../../out/table/df_int_up.tsv")

# confirm the data and the list are congruent
df_int_up %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

# for DOWN list
df1_down <- lapply(list_DE_down,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

head(df1_down)

df2_down <- data.frame(gene=unique(unlist(list_DE_down)))
head(df2_down)

# now loop through each individual gene and pick the list of all the intersections they belong to
df_int_down <- lapply(df2_down$gene,function(x){
  # pull the name of the intersections
  intersection <- df1_down %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

head(df_int_down,n=20)

df_int_down %>%
  write_tsv("../../out/table/df_int_down.tsv")

# confirm the data and the list are congruent
df_int_down %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))
