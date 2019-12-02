library(readr)
library(dplyr)
library(ggplot2)

amp <- read_csv("data/ CraigsAmpData.csv") %>% 
  select(chr = `Zebrafish chromosome`,
         start = `Start position`,
         g = `Amplification_G-score`, 
         q = `Deletion_q-value`
) %>% 
  mutate(chr = sub("zv", "", chr)) %>% 
  mutate(q = ifelse(grepl("[:alpha:]", q), NA, as.numeric(q)))

chr18to20 <- amp %>% 
  filter(chr %in% c(18, 19, 20)) %>% 
  mutate(start = ifelse(chr == 18, yes = paste0("18.", start), 
                        no = ifelse(chr == 19, yes = paste0("19.", start), 
                               no = ifelse(chr == 20, paste0("20.", start), NA)
                               )
                        )
         ) #%>% 
  mutate(start = as.numeric(start))


cols <- c("#A33141", "#3B5998")

lapply(cols, function(col){
  ggplot(chr18to20, aes(x = start, y = g, group = 1)) + 
    geom_line(color = col) +
    geom_line(aes(x = 0:(nrow(chr18to20)-1), y= 0), alpha = 0) +
    geom_ribbon(aes(ymax = g), ymin = 0, color = col, fill = col) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          text = element_blank()
    )
  ggsave(paste0("results/CraigAmp/", col, "_1.png"))
  
  ggplot(chr18to20, aes(x = start, y = g, group = 1)) + 
    geom_line(color = col) + 
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          text = element_blank()
    )
  ggsave(paste0("results/CraigAmp/", col, "_2.png"))
  
  
  ggplot(chr18to20, aes(x = start, y = g, group = 1)) + 
    geom_line(color = col) +
    scale_y_continuous(name = "Amplification G-score") +
    scale_x_discrete(name = "Chromosome", 
                     labels = NULL
    ) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "grey"),
          axis.ticks.x = element_blank())
  ggsave(paste0("results/CraigAmp/", col, "_3.png"))
  
  ggplot(chr18to20, aes(x = start, y = g, group = 1)) + 
    geom_line(color = col) +
    geom_line(aes(x = 0:(nrow(chr18to20)-1), y= 0), alpha = 0) +
    geom_ribbon(aes(ymax = g), ymin = 0, color = col, fill = col) +
    scale_y_continuous(name = "Amplification G-score") +
    scale_x_discrete(name = "Chromosome", 
                     labels = NULL
    ) +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "grey"),
          axis.ticks.x = element_blank())
  ggsave(paste0("results/CraigAmp/", col, "_4.png"))
})
