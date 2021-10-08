##### info ####

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
         N_limit = ifelse(nitrogen_added == "low", 1, 0),
         soil = recode(soil, Sterile = "non-inoculated", A = "ambient N", D = "low N", H = "high N"),
         soil_N = case_when(soil %in% c("non-inoculated", "ambient N") ~ 0,
                            soil == "low N" ~ 34,
                            soil == "high N" ~ 272),
         microbes = ifelse(soil == "non-inoculated", 0, 1)) %>%
  filter(!is.na(pav) & !is.na(rpv)) %>%
  select(-nutrient)
# 6 samples untested

# PCR values
unique(dat2$pav)
unique(dat2$rpv)
# half values when bands were visible, but not as strong as the controls

# unintended inoculations
dat2 %>%
  filter(inoc_pav == 0 & pav %in% c(0.5, 1)) %>%
  group_by(disease, pav) %>%
  count()

dat2 %>%
  filter(inoc_rpv == 0 & rpv %in% c(0.5, 1)) %>%
  group_by(disease, rpv) %>%
  count()

# total inoculations
dat2 %>%
  group_by(disease) %>%
  count()

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

# sample sizes
dat_up %>%
  group_by(disease, soil, N_added) %>%
  count() %>%
  ungroup() %>%
  summarise(min = min(n),
            max = max(n))

dat_dn %>%
  group_by(disease, soil, N_added) %>%
  count() %>%
  ungroup() %>%
  summarise(min = min(n),
            max = max(n))


#### output ####

write_csv(dat_up, "intermediate-data/bydv_microbes_data_rounded_up.csv")
write_csv(dat_dn, "intermediate-data/bydv_microbes_data_rounded_down.csv")
