##### info ####

# authors: Amy Kendig and Casey Easterday
# date last edited: 10/7/20
# goal: analyze data when including all visible PCR bands as indicators of infection (rounding up)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
dat <- read_csv("intermediate-data/bydv_microbes_data_rounded_up.csv")


#### edit data ####

# reorganize factor levels
dat2 <- dat %>%
  mutate(soil = fct_relevel(soil, "sterile", "low", "medium", "high"),
         disease = fct_relevel(disease, "Healthy", "PAV", "RPV"),
         nitrogen_added = fct_relevel(nitrogen_added, "low", "high"))

# seperate by virus
pav_dat <- dat2 %>%
  filter(inoc_pav == 1)

rpv_dat <- dat2 %>%
  filter(inoc_rpv == 1)


#### visualize ####

# PAV infection prevalence
ggplot(pav_dat, aes(soil, pav, color = nitrogen_added)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.2)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  facet_wrap(~ disease)
# adding microbes reduces PAV infection unless N is high (because N reduces infection) or microbes have long-term exposure to high N
# adding microbes or N helps ameliorate negative effects of co-inoculation on PAV infection

# RPV infection prevalence
ggplot(rpv_dat, aes(soil, rpv, color = nitrogen_added)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.2)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  facet_wrap(~ disease)
# adding microbes reduces RPV infection, especially when N is high (because N increases infection)
# adding microbes or N helps ameliorate negative effects of co-inoculation on RPV infection


#### PAV model ####

# initial fit
pav_mod1 <- glm(pav ~ microbes * nitrogen_added * soil_N * inoc_rpv, data = pav_dat, family = binomial)
summary(pav_mod1)
# doesn't work because soil_N only applies to microbes == 1

# modify soil variable
pav_mod2 <- glm(pav ~ soil * nitrogen_added * inoc_rpv, data = pav_dat, family = binomial)
summary(pav_mod2)
# no significant effects


#### RPV model ####

# initial fit
rpv_mod1 <- glm(rpv ~ soil * nitrogen_added * inoc_pav, data = rpv_dat, family = binomial)
summary(rpv_mod1)
# no significant effects