###
# 20190813 JW 
###

library(dplyr)
library(data.table)
library(readr)

source("src/Proteomics/silac_median.R")
sty <- readRDS("data/Proteomics//sty.RDS")
experiments_cols <- readRDS("data/Proteomics/experiments_cols.RDS")

# Summary of variables Phospho (STY)Sites.txt will be filtered on.
# After running the script this information will be stored in the
# variable named report.
report <- vector()
report["Total Sites Identified"] <- nrow(sty)
report["Potential Contaminants"] <- sum(sty$`Potential contaminant` == "+")
report["Reverse"] <- sum(sty$Reverse == "+")
report["Localisation probability < 0.75"] <- sum(sty$`Localization prob` < 0.75)
report["No detected phosphorylation (incomplete cases)"] <-
  sty %>%
  select(experiments_cols) %>%
  filter(rowSums(is.na(.)) == ncol(.)) %>% nrow()


# Filter Phospho (STY)Sites.txt to remove contaminants, reverse, and low localisation
# probability sites. Localisation probability threshold can be changed here.
sty_flt <- sty %>%
  filter(`Potential contaminant` != "+", `Reverse` != "+", `Localization prob` > 0.75)

# Add info to report
report["Sites remaining following filtering"] <- nrow(sty_flt)
report["Incomplete cases in filtered dataset"] <- sty_flt %>% select(experiments_cols) %>%
  filter(rowSums(is.na(.)) == ncol(.)) %>% nrow()

# Calculate median of combined/total SILAC ratios.
# Stored in variable named medians and added to report
medians <- apply(sty[,11:16], 2, median, na.rm=T)
report["Median normalised SILAC ratios (H/L, M/l, H/M)"] <- paste(medians[c(2,4,6)], collapse = " ")

# Filter rows with no SILAC ratios
sty_flt_cc <- sty_flt[rowSums(is.na(sty_flt[experiments_cols])) < length(experiments_cols),] %>% 
  select(-"Ratio M/L", -"Ratio M/L normalized",
        -"Ratio H/L", -"Ratio H/L normalized", 
        -"Ratio H/M", -"Ratio H/M normalized")

saveRDS(sty_flt_cc, "data/Proteomics/sty_flt_cc.RDS")

ML <- silac_median(sty_flt_cc, "M/L", experiments_cols)
HL <- silac_median(sty_flt_cc, "H/L", experiments_cols)


# Combine data to make final output table


colnames(ML) <- c("id", paste0("STY_30m", c("_01", "_02", "_03", "_med")))
colnames(HL) <- c("id", paste0("STY_4h", c("_01", "_02", "_03", "_med")))

FINAL <- merge(ML, HL,by="id") %>% merge(.,sty_flt_cc[,1:7], by="id") %>% select("id", "Protein", "Protein names",
                                                                                 "Gene names", "Amino acid", "Position",
                                                                                 "Sequence window", grep("_med", colnames(.)))
FINAL_3r <- FINAL[rowSums(is.na(FINAL[,grep("_med", colnames(FINAL))])) < length(grep("_med", colnames(FINAL))),]

saveRDS(FINAL_3r, "results/Proteomics/styFinal_3R.rds")

### to calculate if one rep / 3 is NA

ML_2R <- sty_flt_cc %>% select("id",
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

HL_2R <- sty_flt_cc %>% select("id",
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

colnames(ML_2R) <- c("id", paste0("STY_30m", c("_01", "_02", "_03", "_med")))
colnames(HL_2R) <- c("id", paste0("STY_4h", c("_01", "_02", "_03", "_med")))

FINAL_2R <- merge(ML_2R, HL_2R,by="id") %>% merge(.,sty_flt_cc[,1:7], by="id") %>% select("id", "Protein", "Protein names",
                                                                                 "Gene names", "Amino acid", "Position",
                                                                                 "Sequence window", grep("_med", colnames(.)))
FINAL_2R <- FINAL_2R[rowSums(is.na(FINAL_2R[,grep("_med", colnames(FINAL_2R))])) < length(grep("_med", colnames(FINAL_2R))),]

saveRDS(FINAL_2R, "results/Proteomics/styFinal_2R.rds")


#### FIGURES
# Correlation heatmaps
cors <- cor(x = as.matrix(sty_flt_cc[,grep("STY", colnames(sty_flt_cc))]),
            use="pairwise.complete.obs") 


tiff("results/Proteomics/figs/cor_2R.tif")
gplots::heatmap.2(x = cors,
                  col = RColorBrewer::brewer.pal(9, "Blues"),
                  trace = "none",
                  tracecol = NULL,
                  key.title = " ")
dev.off()

# Pie charts to visualise propoprtion of missing values across all
# and phosphorylated residues

png("results/Proteomics/figs/Missing_3R.png")
pie(table(is.na(FINAL_3r[,grep("_med", colnames(FINAL_3r))])),
    labels = c(paste0("NA\n "),
               paste0("Not NA\n ")),
    main = "Missing Ratios in Filtered Data: 3 reps")
dev.off()

png("results/Proteomics/figs/Residues_3R.png")
pie(table(FINAL_3r$`Amino acid`), main= "Phosphorylated residues: 3 reps")
dev.off()

png("results/Proteomics/figs/Missing_2R.png")
pie(table(is.na(FINAL_2R[,grep("_med", colnames(FINAL_2R))])),
    labels = c(paste0("NA\n "),
               paste0("Not NA\n ")),
    main = "Missing Ratios in Filtered Data: 2 reps")
dev.off()

png("results/Proteomics/figs/Residues_2R.png")
pie(table(FINAL_2R$`Amino acid`), main= "Phosphorylated residues: 2 reps")
dev.off()

