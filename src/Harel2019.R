library(dplyr)
library(tidyr)

pg <- readr::read_csv("data/Harel2019_S1B.csv", skip = 1)
colnames(pg) <- gsub("T: ", "", colnames(pg))
colnames(pg) <- gsub(" ", "_", colnames(pg))
experiment_cols <- colnames(pg)[1:116]


highDGAT <- pg %>% 
  filter(Gene_names == "DGAT1") %>% 
  select(all_of(experiment_cols)) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "patient",
    values_to = "log2ratio"
  ) %>% 
  filter(log2ratio >= quantile(log2ratio, prob=1-25/100, na.rm = T)) 

lowDGAT <- pg %>% 
  filter(Gene_names == "DGAT1") %>% 
  select(all_of(experiment_cols)) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "patient",
    values_to = "log2ratio"
  ) %>% 
  filter(log2ratio <= quantile(log2ratio, prob=1-75/100, na.rm = T))

genes <- readr::read_csv("data/Harel2019_S2_sigGenes.csv", col_names = F)[,1, drop=T]

sig <- pg %>% 
  #filter(Gene_names %in% genes) %>% 
  select(c(lowDGAT$patient, highDGAT$patient, "Gene_names")) %>% 
  mutate(p = apply(., 1, function(i){
    high <- as.numeric(i[names(i) %in% highDGAT$patient])
    low <- as.numeric(i[names(i) %in% lowDGAT$patient])
    
    res <- tryCatch(
      t.test(high, low)$p.value, error=function(x) NA 
    )
    return(res)
  }),
  adjp = p.adjust(p, method = "BH"),
  sig = adjp < 0.05
  ) %>% 
  filter(adjp < 0.05)
