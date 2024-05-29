# AIM ---------------------------------------------------------------------
# build the signature object to be used by GSEA

# libraries ---------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(msigdbr)
library(GSEABase)
library(homologene)

# update homologene database for geneID conversion ------------------------
# updateHomologene(destfile = "../../data/Homologene_240528")
homologeneData2_240528 <- read_tsv("../../data/Homologene_240528") %>%
  as.data.frame()

# read in the files -------------------------------------------------------
# read in the gmt siganture files
file_gmt <- dir("../../data/signatures/set_240528",full.names = T) %>%
  str_subset(pattern = ".gmt")

# build the list for loading the files
list_name_gmt <- data.frame(file_gmt = file_gmt) %>%
  mutate(term = str_remove_all(file_gmt,pattern = "../../data/signatures/set_240528/") %>% str_sub(start = 1,end = -16)) %>%
  mutate(source = str_extract(file_gmt,pattern="Hs|Mm")) %>%
  # sort out the human from the mouse signatures
  split(f = .$source)

# wrangling ---------------------------------------------------------------
# build the signature file from the Hs dataset
pathways_gmt_Hs <- lapply(list_name_gmt$Hs$file_gmt, function(file){
  # read in the gene names (human)
  human_gene <- getGmt(file) %>%
    geneIds() %>%
    .[[1]]
  
  # build the conversion table to mouse
  df_homologene <- homologene(human_gene, inTax = 9606, outTax = 10090,db = homologeneData2_240528)
  df_homologene <- df_homologene %>%
    dplyr::rename("mouse_gene"="10090","human_gene"="9606","mouse_geneID"="10090_ID","human_geneID"="9606_ID")
  
  test_df <- data.frame(human_gene = human_gene) %>%
    left_join(df_homologene,by = "human_gene")
  
  return(test_df)
}) %>%
  setNames(list_name_gmt$Hs$term)

# build the signature file from the Mm dataset
pathways_gmt_Mm <- lapply(list_name_gmt$Mm$file_gmt, function(file){
  # read in the gene names (mouse)
  mouse_gene <- getGmt(file) %>%
    geneIds() %>%
    .[[1]]
  
  # build the conversion table to mouse
  df_homologene <- homologene(mouse_gene, inTax = 10090, outTax = 9606,db = homologeneData2_240528)
  # rename the columns
  df_homologene <- df_homologene %>%
    dplyr::rename("mouse_gene"="10090","human_gene"="9606","mouse_geneID"="10090_ID","human_geneID"="9606_ID")
  
  test_df <- data.frame(mouse_gene = mouse_gene) %>%
    left_join(df_homologene,by = "mouse_gene")
  
  return(test_df)
}) %>%
  setNames(list_name_gmt$Mm$term)

# build the unique pathway file
pathway <- c(
  # list from Hs
  lapply(pathways_gmt_Hs,function(geneset){
    # filter the non NA human gene symbol
    gene <- geneset %>%
      filter(!is.na(human_gene)) %>%
      pull(human_gene)
    return(gene)
    }),
  # list from Mm
  lapply(pathways_gmt_Mm,function(geneset){
    # filter the non NA human gene symbol
    gene <- geneset %>%
      filter(!is.na(human_gene)) %>%
      pull(human_gene)
    return(gene)
    })
)

# save the object ---------------------------------------------------------
saveRDS(pathway,"../../out/object/set_240528.rds")

