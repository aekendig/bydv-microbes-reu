#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
dat <- read_csv("data/Cedar_Creek_ELISA_111809.csv")


#### summarize ####

# edit data
dat2 <- dat %>%
  rename(PAV = `PAV call`,
         MAV = `MAV call`,
         RPV = `RPV call`) %>%
  mutate(PAV = case_when(PAV == "PAV" ~ 1,
                         PAV == "-" ~ 0),
         MAV = case_when(MAV == "MAV" ~ 1,
                         MAV == "-" ~ 0),
         RPV = case_when(RPV == "RPV" ~ 1,
                         RPV == "-" ~ 0))

# overall prevalence
sum(dat2$PAV)/nrow(dat2)
sum(dat2$MAV)/nrow(dat2)
sum(dat2$RPV)/nrow(dat2)

# sample size
nrow(dat2)

# plant species
unique(dat2$Spp)
length(unique(dat2$Spp))

# prevalence by species
