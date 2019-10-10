###
# 20190813 JW 
###


library(dplyr)
library(gplots)
library(data.table)
library(readr)
library(RColorBrewer)
library(ggplot2)


# Read experiment names
experiments <- (fread(input = "data/Proteomics/summary_ST.txt", select = "Experiment") %>%
                  unique %>% 
                  filter(grepl("PRO", Experiment)) %>% 
  mutate(Ligand = substr(Experiment, 5, 5),
         Timepoint = substr(Experiment, 7, 9),
         rep = substr(Experiment, 11, 12)) %>% 
  group_by(Ligand) %>% 
  arrange(Timepoint, .by_group=T))$Experiment

# Use experiment names to get column names for proteinGroups_ST.txt
experiments_cols <- c(sapply(experiments, function(e){
  ML = paste("Ratio M/L normalized",e) 
  HL = paste("Ratio H/L normalized",e)
  return(c(ML, HL))
}, USE.NAMES=F))


# Read proteinGroups_ST.txt, specifying only the columns of interest
proteinGroups <- fread(input = "data/Proteomics/proteinGroups_ST.txt",
             select = c("id", "Protein names",
                        "Gene names","Reverse", "Potential contaminant",
                        "Ratio M/L", "Ratio M/L normalized",
                        "Ratio H/L", "Ratio H/L normalized", "Ratio H/M",
                        "Ratio H/M normalized", experiments_cols
             ))

# Summary of variables proteinGroups_ST.txt will be filtered on.
# After running the script this information will be stored in the
# variable named report.
report <- vector()
report["Total Proteins Identified"] <- nrow(proteinGroups)
report["Potential Contaminants"] <- sum(proteinGroups$`Potential contaminant` == "+")
report["Reverse"] <- sum(proteinGroups$Reverse == "+")
report["No Quantitative Data (incomplete cases)"] <- proteinGroups %>% select(experiments_cols) %>%
  filter(rowSums(is.na(.)) == ncol(.)) %>% nrow()


# Filter proteinGroups.txt to remove contaminants, reverse.
proteinGroups_flt <- proteinGroups %>%
  filter(`Potential contaminant` != "+", `Reverse` != "+")

# Add info to report
report["Sites remaining following filtering"] <- nrow(proteinGroups_flt)
report["Incomplete cases in filtered dataset"] <- proteinGroups_flt %>% select(experiments_cols) %>%
  filter(rowSums(is.na(.)) == ncol(.)) %>% nrow()

# Calculate median of combined/total SILAC ratios.
# Stored in variable named medians and added to report
medians <- apply(proteinGroups[,11:16], 2, median, na.rm=T)
report["Median normalised SILAC ratios (H/L, M/l, H/M)"] <- paste(medians[c(2,4,6)], collapse = " ")

# Filter rows with no SILAC ratios
proteinGroups_flt_cc <- proteinGroups_flt[rowSums(is.na(proteinGroups_flt[experiments_cols])) < length(experiments_cols),] %>% 
  select(-"Ratio M/L", -"Ratio M/L normalized",
         -"Ratio H/L", -"Ratio H/L normalized", 
         -"Ratio H/M", -"Ratio H/M normalized")


# Define function to calculate median SILAC ratio between replicates.
silac_median <- function(filtered_proteinGroups, ML_or_HL, experiment_names){
  ratios <- filtered_proteinGroups %>% select("id",
                                    grep(paste("Ratio",ML_or_HL),
                                         colnames(filtered_proteinGroups)))
  
  colnames(ratios)[colnames(ratios) %in% experiment_names] <- 
    sapply(colnames(ratios)[colnames(ratios) %in% experiment_names],
           function(x) substr(x,22,max(nchar(x))))
  
  x <- seq(2, ncol(ratios), 2)
  
  for(i in x){
    name <- substr(colnames(ratios)[i],1,nchar(colnames(ratios)[i])-3)
    ratios[,name] <- apply(ratios[,i:(i+1)], 1, median, na.rm=F)
  }
  
  return(ratios)
}  

ML <- silac_median(proteinGroups_flt_cc, "M/L", experiments_cols)
HL <- silac_median(proteinGroups_flt_cc, "H/L", experiments_cols)

colnames(ML) <- c("id", paste0("PRO_30m", c("_01", "_02", "_03", "_med")))
colnames(HL) <- c("id", paste0("PRO_4h", c("_01", "_02", "_03", "_med")))

FINAL <- merge(ML, HL,by="id") %>%
  merge(.,proteinGroups_flt_cc[,1:7], by="id") %>%
  select("id", "Protein names", "Gene names",
         grep("_med", colnames(.)))

FINAL_3R <- FINAL[rowSums(is.na(FINAL[,grep("_med", colnames(FINAL))])) < length(grep("_med", colnames(FINAL))),]

saveRDS(FINAL_3R, "results/Proteomics/proteinFinal_3R.rds")
readr::write_csv(FINAL_3R, "results/Proteomics/proteinFinal_3R.csv")

### to calculate if one rep / 3 is NA

ML_2R <- proteinGroups_flt_cc %>% select("id",
                               grep("Ratio M/L",
                                    colnames(.))) %>% 
  mutate(keep = apply(., 1, function(x) sum(is.na(x)) < 2)) %>% 
  filter(keep == T) %>% 
  select(1:(ncol(.) - 1))

colnames(ML_2R)[colnames(ML_2R) %in% experiments_cols] <- 
  sapply(colnames(ML_2R)[colnames(ML_2R) %in% experiments_cols],
         function(x) substr(x,22,max(nchar(x))))

x <- seq(2, ncol(ML_2R), 2)

for(i in x){
  name <- substr(colnames(ML_2R)[i],1,nchar(colnames(ML_2R)[i])-3)
  ML_2R[,name] <- apply(ML_2R[,i:(i+1)], 1, median, na.rm=T)
}

HL_2R <- proteinGroups_flt_cc %>% select("id",
                               grep("Ratio H/L",
                                    colnames(.))) %>% 
  mutate(keep = apply(., 1, function(x) sum(is.na(x)) < 2)) %>% 
  filter(keep == T) %>% 
  select(1:(ncol(.) - 1))

colnames(HL_2R)[colnames(HL_2R) %in% experiments_cols] <- 
  sapply(colnames(HL_2R)[colnames(HL_2R) %in% experiments_cols],
         function(x) substr(x,22,max(nchar(x))))

x <- seq(2, ncol(HL_2R), 2)

for(i in x){
  name <- substr(colnames(HL_2R)[i],1,nchar(colnames(HL_2R)[i])-3)
  HL_2R[,name] <- apply(HL_2R[,i:(i+1)], 1, median, na.rm=F)
}

colnames(ML_2R) <- c("id", paste0("PRO_30m", c("_01", "_02", "_03", "_med")))
colnames(HL_2R) <- c("id", paste0("PRO_4h", c("_01", "_02", "_03", "_med")))

FINAL_2R <- merge(ML_2R, HL_2R,by="id") %>%
  merge(.,proteinGroups_flt_cc[,1:7], by="id") %>%
  select("id", "Protein names",
         "Gene names", grep("_med", colnames(.)))

saveRDS(FINAL_2R, "results/Proteomics/proFinal_2R.rds")
readr::write_csv(FINAL_2R, "results/Proteomics/proFinal_2R.csv")



#### FIGURES
# Pie chart to visualise propoprtion of missing values across all experiments.
pie(table(is.na(FINAL_x[,8:23])), labels = c("Not NA", "NA"), main = "Missing Ratios in Filtered Data")

#Heatmap correlation
cors <- proteinGroups %>%
  filter(id %in% FINAL_3R$id) %>% 
  select(grep("PRO", colnames(.))) %>% 
  cor(use = "pairwise.complete.obs")

tiff("results/Proteomics/figs/proCor_3R.tif")
heatmap.2(x = cors,
          col = RColorBrewer::brewer.pal(9, "Blues"),
          trace = "none",
          tracecol = NULL,
          key.title = " "
)
dev.off()

cors <- proteinGroups %>%
  filter(id %in% FINAL_2R$id) %>% 
  select(grep("PRO", colnames(.))) %>% 
  cor(use = "pairwise.complete.obs")

tiff("results/Proteomics/figs/proCor_2R.tif")
heatmap.2(x = cors,
          col = RColorBrewer::brewer.pal(9, "Blues"),
          trace = "none",
          tracecol = NULL,
          key.title = " "
)
dev.off()