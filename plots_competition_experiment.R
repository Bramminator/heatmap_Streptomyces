.libPaths('~/R/x86_64-pc-linux-gnu-library/3.6')
library(tidyverse)
library(ggplot2)
library(stringr)
library(stringi)
library(data.table)

mutant30 <- list.dirs("/hosts/linuxhome/mutant30/tmp/bramve/", full.names = TRUE, recursive = FALSE)

Time_vector <- seq(19000, 99000, by = 10000)

popfile_vector <- c()

for(files in mutant30){
  popfile <- paste0(files, "/popsize2.txt")
  popfile_vector <- c(popfile_vector, popfile)
}

ratio_frame <- data.frame(rearrangement_rate1 = numeric(0), rearrangement_rate2 = numeric(0), cost = character(0), production = numeric(0), outcome = numeric(0))

for(i in popfile_vector){
  temp <- fread(i)
  colnames(temp) <- c("Time", "popwt1", "popmu1", "popse", "popwt2", "popmu2")
  temp$Time <- as.numeric(temp$Time)
  ratios <- c()
  for(Timepoint in Time_vector){
    pop_ratio <- as.numeric(temp$popwt1[temp$Time == Timepoint]) / as.numeric(temp$popwt2[temp$Time == Timepoint])
    ratios <- c(ratios, pop_ratio)
  }
  outcome <- sum(ratios) / 9
  r1 <- str_match(i, 'r1(.*?)_r')
  stri_sub(r1[2],2,1) <- '.0'
  r1 <- as.numeric(r1[2])
  r2 <- str_match(i, '_r2(.*?)/p')
  r2 <- r2[2]
  if(str_detect(r2, ".0", negate = TRUE)){
    stri_sub(r2,2,1) <- '.0'
  }
  r2 <- as.numeric(r2)
  c <- str_match(i,'t_c(.*?)_pr')
  stri_sub(c, 2, 1) <- '.'
  cost <- as.character(c[2])
  pr <- str_match(i,'_pr(.*?)_r1')
  stri_sub(pr,2,1) <- '.'
  production <- as.numeric(pr[2])
  value_frame <- data.frame(r1, r2, cost, production, outcome)
  colnames(value_frame) <- c('rearrangement_rate1', 'rearrangement_rate2', 'cost', 'production', 'outcome')
  ratio_frame <- bind_rows(ratio_frame, value_frame)
}

print(ratio_frame)

ratio_frame %>% 
  filter(cost == "0.2") %>%
  ggplot(aes(x = rearrangement_rate2, y = rearrangement_rate1, fill = outcome)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'green', mid = 'white', midpoint = 1, limit = c(0.2, 1.8)) +
  theme_minimal() +
  geom_text(aes(x = rearrangement_rate2, y = rearrangement_rate1, label = round(outcome, digits = 3))) + 
  ggtitle('Ratio population sizes for strains with varying mutation rates (cost 0.2)') +
  
  
  ggsave('~/Documents/strepto/heatmap_competition_experiment/heatmap02.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)
  
ratio_frame %>% 
  filter(cost == "0.4") %>%
  ggplot(aes(x = rearrangement_rate2, y = rearrangement_rate1, fill = outcome)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'green', mid = 'white', midpoint = 1, limit = c(0.2, 1.8)) +
  theme_minimal() +
  geom_text(aes(x = rearrangement_rate2, y = rearrangement_rate1, label = round(outcome, digits =3))) + 
  ggtitle('Ratio population sizes for strains with varying mutation rates (cost 0.4)') +
    
    
  ggsave('~/Documents/strepto/heatmap_competition_experiment/heatmap04.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)  

ratio_frame %>% 
  filter(cost == "0.6") %>%
  ggplot(aes(x = rearrangement_rate2, y = rearrangement_rate1, fill = outcome)) +
  geom_tile(color = 'white') +
  scale_fill_gradient2(low = 'red', high = 'green', mid = 'white', midpoint = 1, limit = c(0.2, 1.8)) +
  theme_minimal() +
  geom_text(aes(x = rearrangement_rate2, y = rearrangement_rate1, label = round(outcome,digits = 3))) + 
  ggtitle('Ratio population sizes for strains with varying mutation rates (cost 0.6)') +
  
  
  ggsave('~/Documents/strepto/heatmap_competition_experiment/heatmap06.png',
         height = 210,
         width  = 297,
         units= 'mm',
         dpi= 150)  

