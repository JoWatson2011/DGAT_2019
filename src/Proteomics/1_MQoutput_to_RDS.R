library(dplyr)
library(data.table)

# Read experiment names
experiments <- (fread(input = "data/Proteomics/summary_ST.txt", select = "Experiment") %>%
                  unique %>%
                  filter(grepl("STY", Experiment)) %>% 
                  mutate(Ligand = substr(Experiment, 5, 5),
                         Timepoint = substr(Experiment, 7, 9),
                         rep = substr(Experiment, 11, 12)) %>% 
                  group_by(Ligand) %>% 
                  arrange(Timepoint, .by_group=T))$Experiment

# Use experiment names to get column names for phospho (STY)Sites.txt
experiments_cols <- c(sapply(experiments, function(e){
  ML = paste("Ratio M/L normalized",e)
  HL = paste("Ratio H/L normalized",e)
  return(c(ML, HL))
}, USE.NAMES=F))

saveRDS(experiments_cols,"data/Proteomics/experiments_cols.RDS")

sty <- fread(input = "data/Proteomics/Phospho (STY)Sites.txt",
             select = c("id", "Protein", "Protein names",
                        "Gene names",	"Amino acid", "Position",
                        "Sequence window", "Reverse", "Potential contaminant",
                        "Localization prob", "Ratio M/L", "Ratio M/L normalized",
                        "Ratio H/L", "Ratio H/L normalized", "Ratio H/M",
                        "Ratio H/M normalized","Mod. peptide IDs", experiments_cols
             ))


saveRDS(sty, "data/Proteomics/sty.RDS")


# Read experiment names
experiments <- (fread(input = "data/Proteomics/summary_ST.txt", select = "Experiment") %>%
                  filter(Experiment != "") %>% 
                  unique())$Experiment

# Use experiment names to get column names for modificationSpecificPeptides.txt
experiments_cols <- c(sapply(experiments, function(e){
  E = paste("Experiment",e)
  ML = paste("Ratio M/L normalized",e) 
  HL = paste("Ratio H/L normalized",e)
  return(c(E, ML, HL))
}, USE.NAMES=F)) %>% sort()

