library(tidyverse)
library(microbenchmark)
library(profvis)
library(Rcpp)
library(splines)
library(xtable)
horse_data <- readr::read_csv("4_Horses.csv", col_types = cols(dead = col_logical()))
horse_missing <- horse_data %>%
  mutate(missing = ifelse(is.na(Temperature) == TRUE, TRUE,FALSE)) %>%
  group_by(dead,missing) %>%
  summarise(n = n())
