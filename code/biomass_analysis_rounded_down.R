##### info ####

# authors: Amy Kendig and Casey Easterday
# date last edited: 11/25/20
# goal: analyze biomass data when only PCR bands at least as intense as controls are counted as indicators of infection (rounding down)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)

# import data
dat <- read_csv("intermediate-data/bydv_microbes_data_rounded_down.csv")


#### edit data ####

# reorganize factor levels
# select for successful infection
dat2 <- dat %>%
  mutate(soil = fct_relevel(soil, "sterile", "ambient N", "low N"),
         infection = case_when(disease == "PAV" & pav == 1 ~ "PAV infection",
                               disease == "RPV" & rpv == 1 ~ "RPV infection",
                               disease == "Co" & pav == 1 & rpv == 1 ~ "Co-infection",
                               disease == "Healthy" ~ "Mock inoculation",
                               TRUE ~ NA_character_) %>%
           fct_relevel("Mock inoculation", "PAV infection", "RPV infection"),
         nitrogen_added = fct_relevel(nitrogen_added, "low", "high"),
         log_biomass = log(biomass)) %>%
  filter(!is.na(infection))

# remove coinfection for stats
bio_dat <- dat2 %>%
  filter(infection != "Co-infection")

# soil N comparison
bio_soil_dat = filter(bio_dat, microbes == 1)


#### visualize ####

# small figures: 80 mm, large figures: 180 mm

# theme
theme_def <- theme_bw() +
  theme(axis.text = element_text(size = 8, color="black"),
        axis.title = element_text(size = 10, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, -10, -10, -10),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, color="black"))

# palettes
col_pal = c("white", "black")

# panel labels
pan_lab <- tibble(infection = unique(dat2$infection),
                  label = c("(a)", "(b)", "(c)", "(d)")) %>%
  mutate(soil = "sterile",
         biomass = 0.65,
         nitrogen_added = "low")

# biomass figure
pdf("output/biomass_figure_rounded_down.pdf", width = 6, height = 6)
ggplot(dat2, aes(soil, biomass, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2), shape = 21) +
  geom_text(data = pan_lab, aes(label = label), nudge_x = -0.4) +
  facet_wrap(~ infection) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  guides(fill = guide_legend(override.aes = list(shape = 21), direction = "horizontal", title.position = "top")) +
  xlab("Field soil N treatment") +
  ylab("Biomass (g)") +
  theme_def +
  theme(legend.position = c(0.88, 0.4))
dev.off()


#### all data model ####

# initial fit
bio_mod1 <- brm(log_biomass ~ soil * N_added * infection,
                     data = bio_dat, 
                     family = gaussian,
                     prior = c(prior(normal(0, 10), class = Intercept),
                               prior(normal(0, 10), class = b),
                               prior(cauchy(0, 1), class = sigma)),
                     iter = 6000, warmup = 1000, chains = 1)
summary(bio_mod1)

# increase chains
bio_mod2 <- update(bio_mod1, chains = 3)
summary(bio_mod2)
plot(bio_mod2)
pp_check(bio_mod2, nsamples = 50)


#### soil N model ####

# initial fit
bio_soil_mod1 <- brm(log_biomass ~ soil_N * N_added * infection,
                data = bio_soil_dat, 
                family = gaussian,
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 10), class = b),
                          prior(cauchy(0, 1), class = sigma)),
                iter = 6000, warmup = 1000, chains = 1)
summary(bio_soil_mod1)

# increase chains
bio_soil_mod2 <- update(bio_soil_mod1, chains = 3)
summary(bio_soil_mod2)
plot(bio_soil_mod2)
pp_check(bio_soil_mod2, nsamples = 50)


#### numbers for text ####

# all data model
(bio_samps <- posterior_samples(bio_mod2) %>%
  mutate(int = exp(b_Intercept),
         N_eff = exp(b_Intercept + b_N_added) - int,
         amb_eff = exp(b_Intercept + b_soilambientN) - exp(b_Intercept),
         PAV_eff = exp(b_Intercept + b_infectionPAVinfection) - exp(b_Intercept),
         RPV_eff = exp(b_Intercept + b_infectionRPVinfection) - exp(b_Intercept),
         N_eff_amb = exp(b_Intercept + b_soilambientN + b_N_added) - exp(b_Intercept + b_soilambientN)) %>%
  select(int:N_eff_amb) %>%
  pivot_longer(cols = int:N_eff_amb,
               names_to = "treatment",
               values_to = "biomass.g") %>%
  group_by(treatment) %>%
  mean_hdi(biomass.g))

# all data model
(bio_soil_samps <- posterior_samples(bio_soil_mod2) %>%
    mutate(int = exp(b_Intercept),
           N_eff = exp(b_Intercept + b_N_added) - int,
           soil_eff = exp(b_Intercept + b_soil_N) - exp(b_Intercept),
           PAV_eff = exp(b_Intercept + b_infectionPAVinfection) - exp(b_Intercept),
           RPV_eff = exp(b_Intercept + b_infectionRPVinfection) - exp(b_Intercept)) %>%
    select(int:RPV_eff) %>%
    pivot_longer(cols = int:RPV_eff,
                 names_to = "treatment",
                 values_to = "biomass.g") %>%
    group_by(treatment) %>%
    mean_hdi(biomass.g))



#### output ####

save(bio_mod2, file = "output/biomass_full_model_rounded_down.rda")
save(bio_soil_mod2, file = "output/biomass_soil_model_rounded_down.rda")
