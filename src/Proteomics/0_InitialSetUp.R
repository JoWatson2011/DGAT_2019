#######
# JW 190814
# Install all required packages for pipeline
#######

if (!requireNamespace("dplyr", quietly = TRUE)){
  install.packages("dplyr")
}

if (!requireNamespace("data.table", quietly = TRUE)){
  install.packages("data.table")
}

if (!requireNamespace("ggplot2", quietly = TRUE)){
  install.packages("ggplot2")
}

if (!requireNamespace("gplots", quietly = TRUE)){
  install.packages("gplots")
}

if (!requireNamespace("readr", quietly = TRUE)){
  install.packages("readr")
}

if (!requireNamespace("RColorBrewer", quietly = TRUE)){
  install.packages("RColorBrewer")
}