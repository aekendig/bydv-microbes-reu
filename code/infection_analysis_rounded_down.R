##### info ####

# authors: Amy Kendig and Casey Easterday
# date last edited: 10/8/20
# goal: analyze data when only PCR bands at least as intense as controls are counted as indicators of infection (rounding down)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(MuMIn)

# import data
dat <- read_csv("intermediate-data/bydv_microbes_data_rounded_down.csv")


#### edit data ####

# reorganize factor levels
dat2 <- dat %>%
  mutate(soil = fct_relevel(soil, "sterile", "low", "medium", "high"),
         disease = fct_relevel(disease, "Healthy", "PAV", "RPV"),
         nitrogen_added = fct_relevel(nitrogen_added, "low", "high"))

# separate by virus
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
# adding microbes reduces PAV infection unless N is high (because N reduces infection)
# adding microbes helps ameliorate negative effects of co-inoculation on PAV infection

# RPV infection prevalence
ggplot(rpv_dat, aes(soil, rpv, color = nitrogen_added)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.2)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  facet_wrap(~ disease)
# adding microbes reduces RPV infection slightly when N is low or microbes come from a high N environment and N is high
# adding N helps ameliorate negative effects of co-inoculation on RPV infection


#### PAV model ####

# initial fit
pav_mod1 <- glm(pav ~ soil * N_added * inoc_rpv, data = pav_dat, family = binomial,
                na.action = na.fail)
summary(pav_mod1)
# all the standard errors are the same because the intercept is only 1's

# model average using AIC
pav_mod_avg <- model.avg(get.models(dredge(pav_mod1), subset = cumsum(weight) <= .95))
summary(pav_mod_avg)


#### RPV model ####

# initial fit
rpv_mod1 <- glm(rpv ~ soil * nitrogen_added * inoc_pav, data = rpv_dat, family = binomial,
                na.action = na.fail)
summary(rpv_mod1)
# no significant effects

# model average using AIC
rpv_mod_avg <- model.avg(get.models(dredge(rpv_mod1), subset = cumsum(weight) <= .95))
summary(rpv_mod_avg)
