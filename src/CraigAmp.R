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
         )

chr23to25 <- amp %>% 
  filter(chr %in% c(23, 24, 25)) %>% 
  mutate(start = ifelse(chr == 23, yes = paste0("18.", start), 
                        no = ifelse(chr == 24, yes = paste0("19.", start), 
                                    no = ifelse(chr == 25, paste0("20.", start), NA)
                        )
  )
  )

ggplot(chr18to20, aes(x = start, y = g, group = 1)) + 
    geom_line(color = "#d81d1d") + 
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
ggsave(paste0("results/CraigAmp/chr18_20_", "#d81d1d", "_2.png"))

ggplot(chr23to25, aes(x = start, y = g, group = 1)) + 
  geom_line(color = "#d81d1d") + 
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
ggsave(paste0("results/CraigAmp/chr23_25_", "#d81d1d", "_2.png"))
  