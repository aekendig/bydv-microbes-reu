##### info ####

# authors: Amy Kendig and Casey Easterday
# date last edited: 10/7/20
# goal: prepare data for statistical analyses


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
dat <- read_csv("data/REU Infection Data Final CAE 9-24-16.csv")


#### edit data ####

# inoculations
unique(dat$disease)

# make column for inoculation
# remove untested samples
# N addition columns
# long-term N treatment columns
dat2 <- dat %>%
  mutate(inoc_pav = ifelse(disease %in% c("PAV", "Co"), 1, 0),
         inoc_rpv = ifelse(disease %in% c("RPV", "Co"), 1, 0),
         nitrogen_added = ifelse(nutrient == "Ctrl", "low", "high"),
         N_added = ifelse(nitrogen_added == "low", 0, 1),
         soil = recode(soil, Sterile = "sterile", A = "low", D = "medium", H = "high"),
         soil_N = case_when(soil %in% c("sterile", "low") ~ 0,
                            soil == "medium" ~ 34,
                            soil == "high" ~ 272),
         microbes = ifelse(soil == "sterile", 0, 1)) %>%
  filter(!is.na(pav) & !is.na(rpv)) %>%
  select(-nutrient)
# 6 samples untested

# PCR values
unique(dat2$pav)
unique(dat2$rpv)
# half values when bands were visible, but not as strong as the controls

# round half values up
# remove contaminated samples
dat_up <- dat2 %>%
  mutate(pav = ceiling(pav),
         rpv = ceiling(rpv),
         co = ifelse(pav == 1 & rpv == 1, 1, 0)) %>%
  filter(!(inoc_rpv == 0 & rpv == 1) & !(inoc_pav == 0 & pav == 1))
# 31 samples contaminated

# round half values down
# remove contaminated samples
dat_dn <- dat2 %>%
  mutate(pav = floor(pav),
         rpv = floor(rpv),
         co = ifelse(pav == 1 & rpv == 1, 1, 0)) %>%
  filter(!(inoc_rpv == 0 & rpv == 1) & !(inoc_pav == 0 & pav == 1))
# 7 samples contaminated


#### output ####

write_csv(dat_up, "intermediate-data/bydv_microbes_data_rounded_up.csv")
write_csv(dat_dn, "intermediate-data/bydv_microbes_data_rounded_down.csv")
