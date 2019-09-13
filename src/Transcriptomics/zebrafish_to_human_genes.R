if (!requireNamespace("biomaRt", quietly = TRUE)){
  install.packages("biomaRt")
}
library(biomaRt)

if (!requireNamespace("readr", quietly = TRUE)){
  install.packages("readr")
}
library(readr)

if (!requireNamespace("dplyr", quietly = TRUE)){
  install.packages("dplyr")
}
library(dplyr)


# Import transcriptomics data
df <- readr::read_csv("data/Transcriptomics/GFP_DGAT1_sigresults_SO_DW Copy.csv") %>% 
  rename("ensembl_gene_id" = 1)

# Perform database search with BioMart:
# match zebrafish ensembl gene ids to zebrafish gene names and human gene names
ensembl_drerio <- useMart("ensembl", dataset = "drerio_gene_ensembl" )

G_list <- getBM(filters = "ensembl_gene_id", 
                attributes= c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_associated_gene_name"),
                values= df$ensembl,
                mart= ensembl_drerio)


# Merge the original dataframe and the dataframe with the human gene names
final <- merge(df, 
               G_list,
               by = "ensembl_gene_id") %>% 
  select(1, 8:9, 2:7)

# Export
write_csv(final, "results/Transcriptomics/GFP_DGAT_sigresults_SO_humanid.csv")
