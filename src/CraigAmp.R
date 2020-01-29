library(readr)
library(dplyr)
library(ggplot2)

amp <- read_csv("data/ CraigsAmpData.csv") %>% 
  select(chr = `Zebrafish chromosome`,
         start = `Start position`,
         g = `Amplification_G-score`
) %>% 
  mutate(chr = sub("zv", "", chr)) 


chr1to3 <- amp %>%
  filter(chr %in% c(1, 2, 3)) %>%
  mutate(start = ifelse(chr == 1, yes = paste0("18.", start),
                        no = ifelse(chr == 2, yes = paste0("19.", start),
                                    no = ifelse(chr == 3, paste0("20.", start), NA)
                        )
  )
  )

chr15to17 <- amp %>%
  filter(chr %in% c(15, 16, 17)) %>%
  mutate(start = ifelse(chr == 15, yes = paste0("18.", start),
                        no = ifelse(chr == 16, yes = paste0("19.", start),
                                    no = ifelse(chr == 17, paste0("20.", start), NA)
                        )
  )
  )


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

ggplot(chr1to3, aes(x = start, y = g, group = 1)) + 
  geom_line(color = "#d81d1d") + 
  scale_y_continuous(limits = c(0, 1550)) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title = element_blank()
  )
ggsave(paste0("results/CraigAmp/chr1_3.png"))

ggplot(chr15to17, aes(x = start, y = g, group = 1)) + 
  geom_line(color = "#d81d1d") + 
  scale_y_continuous(limits = c(0, 1550)) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title = element_blank()
  )
ggsave(paste0("results/CraigAmp/chr15_17.png"))



ggplot(chr18to20, aes(x = start, y = g, group = 1)) + 
    geom_line(color = "#d81d1d") + 
  scale_y_continuous(limits = c(0, 1550)) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title = element_blank()
  )
ggsave(paste0("results/CraigAmp/chr18_20.png"))

ggplot(chr23to25, aes(x = start, y = g, group = 1)) + 
  geom_line(color = "#d81d1d") + 
  scale_y_continuous(limits = c(0, 1550)) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title = element_blank()
  )
ggsave(paste0("results/CraigAmp/chr23_24.png"))

  