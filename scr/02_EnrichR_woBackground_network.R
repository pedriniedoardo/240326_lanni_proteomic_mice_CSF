# libraries ---------------------------------------------------------------
library(epiR)
library(tidyverse)

# try to run the analysis from the enrichmenta table available ------------
# load the table of enrichment analsyis
df <- read_tsv("../../out/table/enrichR_DE_UP_woBackground.tsv")

# pull all the Hits genes
list_test2 <- df %>%
  filter(annotation == "Reactome_2022") %>%
  split(f=.$Term) %>%
  lapply(function(x){
    x %>%
      pull(Genes) %>%
      str_split(pattern = ";") %>%
      unlist() %>%
      unique()
  })

# pull all the unique genes
ref_gene2 <- data.frame(gene = unlist(list_test2) %>% unique())

# build the adiecency columns
list_reactome_gene2t2 <- lapply(list_test2,function(x){
  
  df_test <- x %>%
    data.frame(gene = .) %>%
    mutate(test = 1)
  
  ref_gene2 %>%
    left_join(df_test,by = c("gene")) %>%
    mutate(test = case_when(is.na(test)~0,
                            T~test))
})

# build the dataset with the pairing of the terms
df_pairs2 <- crossing(source = names(list_reactome_gene2t2),target = names(list_reactome_gene2t2)) %>%
  # filter out the equal
  filter(source != target) %>%
  # filter out the same combination swapped
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(names_to = "def",values_to = "term",-id) %>%
  group_by(id) %>%
  arrange(id,term) %>%
  mutate(test = paste0(term,collapse = "|")) %>%
  ungroup() %>%
  group_by(test) %>%
  summarise() %>%
  separate(test,into = c("source","target"),sep = "\\|")

# for each pair of terms mesure the kappa similarity score
df_similarity2 <- df_pairs2 %>%
  # slice(1:10) %>%
  pmap(function(source,target){
    mat_test2 <- table(list_reactome_gene2t2[[source]]$test,
                       list_reactome_gene2t2[[target]]$test)
    data.frame(source,target) %>%
      bind_cols(epi.kappa(mat_test2)$kappa)
  }) %>%
  bind_rows()

# filter only the edges with a similarity score > 0.3
df_similarity2_filter <- df_similarity2 %>%
  filter(est>0.3)

# use this as a edge table.
df_similarity2_filter %>%
  write_tsv(file = "../../out/table/enrichR_DE_UP_woBackground_edge_REACTOME.tsv")

# from the edge table build the network.
# export the node table and add the missing information
# unique terms
df_node <- read_csv("../../out/table/enrichR_DE_UP_woBackground_edge_REACTOME.tsv default node.csv")

id_select <- lapply(df_node,function(x){
  sum(is.na(x)) < 1
}) %>%
  unlist()

df_node_final <- df_node[,id_select] %>%
  left_join(df,by = c("name"="Term")) %>%
  dplyr::rename(glayCluster = `__glayCluster`)

df_node_final %>%
  write_csv("../../out/table/enrichR_DE_UP_woBackground_edge_REACTOME.tsv default node_full.csv")
