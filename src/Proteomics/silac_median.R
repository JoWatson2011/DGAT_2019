silac_median <- function(filtered_data, ML_or_HL, experiment_names){
  ratios <- filtered_data %>% select("id",
                                    grep(paste("Ratio",ML_or_HL),
                                         colnames(filtered_data)))
  
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
