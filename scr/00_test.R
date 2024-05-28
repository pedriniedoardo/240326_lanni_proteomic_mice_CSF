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
