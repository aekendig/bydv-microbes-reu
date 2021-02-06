##### info ####

# goal: analyze biomass data any visible PCR band is counted as indicators of infection (rounding up)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(broom)
library(brms)
library(tidybayes)

# import data
dat <- read_csv("intermediate-data/bydv_microbes_data_rounded_up.csv")


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
         log_biomass = log(biomass),
         microbes = ifelse(soil == "sterile", 0, 1),
         microbes_f = ifelse(microbes == 0, "sterile", "microbes") %>%
           fct_relevel("sterile")) %>%
  filter(!is.na(infection))

# replicates
dat2 %>%
  group_by(infection, N_added, soil) %>%
  count() %>%
  data.frame()
# co-infection only has 1 rep for some treatments

# dataset without coinfection
dat3 <- dat2 %>%
  filter(infection != "Co-infection")


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
pan_lab <- tibble(infection =levels(dat2$infection) %>%
                    fct_relevel("Mock inoculation", "PAV infection", "RPV infection"),
                  label = c("(a)", "(b)", "(c)", "(d)")) %>%
  mutate(microbes_f = "sterile",
         biomass = 0.64,
         nitrogen_added = "low")

# biomass figure
pdf("output/biomass_figure_microbes.pdf", width = 5, height = 5)
ggplot(dat2, aes(microbes_f, biomass, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2), shape = 21) +
  geom_text(data = pan_lab, aes(label = label), nudge_x = -0.45, fontface = "bold") +
  facet_wrap(~ infection) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  guides(fill = guide_legend(override.aes = list(shape = 21), direction = "horizontal", title.position = "top")) +
  xlab("Microbe inoculation") +
  ylab("Aboveground biomass [g/plant]") +
  theme_def +
  theme(legend.position = c(0.85, 0.38))
dev.off()


#### biomass model ####

# initial fit
bio_mod1 <- brm(log_biomass ~ soil * N_added * infection, 
                data = dat2,
                family = gaussian,
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 10), class = b)),
                iter = 6000, warmup = 1000, chains = 1)
summary(bio_mod1)

# increase chains
bio_mod2 <- update(bio_mod1, chains = 3)

# check model
summary(bio_mod2)
plot(bio_mod2)
pp_check(bio_mod2, nsamples = 50)

# microbes model
bio_mic_mod1 <- update(bio_mod2,
                       newdata = dat2,
                       formula = log_biomass ~ microbes * N_added * infection,
                       prior = c(prior(normal(0, 10), class = Intercept),
                                 prior(normal(0, 10), class = b)))

# check model
summary(bio_mic_mod1)
plot(bio_mic_mod1)
pp_check(bio_mic_mod1, nsamples = 50)

# compare with loo
bio_loo2 <- loo(bio_mod2, reloo = T)
bio_loo2 
# all k < 0.7 and elpd_loo <= 0.1
# good model fit
# reasonable to do model comparison
bio_mic_loo1 <- loo(bio_mic_mod1, reloo = T)
bio_mic_loo1 
# all k < 0.7 and elpd_loo <= 0.1
# good model fit
# reasonable to do model comparison
loo_compare(bio_loo2, bio_mic_loo1)
# microbes model is preferred


#### biomass model, no co-infection ####

# initial fit
bio_mod3 <- brm(log_biomass ~ soil * N_added * infection, 
                data = dat3,
                family = gaussian,
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 10), class = b)),
                iter = 6000, warmup = 1000, chains = 1)
summary(bio_mod3)

# increase chains
bio_mod4 <- update(bio_mod3, chains = 3)

# check model
summary(bio_mod4) # estimates are nearly identical to bio_mod2

# microbes model
bio_mic_mod2 <- update(bio_mod4,
                       newdata = dat3,
                       formula = log_biomass ~ microbes * N_added * infection,
                       prior = c(prior(normal(0, 10), class = Intercept),
                                 prior(normal(0, 10), class = b)))

# check model
summary(bio_mic_mod2) # estimates are nearly identical to bio_mic_mod1


#### values for text ####

bio_mic_post1 <- posterior_samples(bio_mic_mod1) %>%
  mutate(icp = exp(b_Intercept),
         high_N = exp(b_Intercept + b_N_added),
         N_effect = 100*(high_N - icp)/icp)

mean_hdi(bio_mic_post1$icp)
mean_hdi(bio_mic_post1$N_effect)


#### output ####
save(bio_mod2, file = "output/bio_bayesian_model_soil.rda")
save(bio_mic_mod1, file = "output/bio_bayesian_model_microbes.rda")

write_csv(tidy(summary(bio_mod2)$fixed), "output/bio_bayesian_model_soil.csv")
write_csv(tidy(summary(bio_mic_mod1)$fixed), "output/bio_bayesian_model_microbes.csv")
