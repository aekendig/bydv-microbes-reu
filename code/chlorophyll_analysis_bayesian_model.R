##### info ####

# goal: analyze chlorophyll data any visible PCR band is counted as indicators of infection (rounding up)


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
  rowwise() %>%
  mutate(chlorophyll = mean(c(chlorophyll_1, chlorophyll_2, chlorophyll_3), na.rm = T)) %>%
  ungroup() %>%
  mutate(soil = fct_relevel(soil, "sterile", "ambient N", "low N"),
         infection = case_when(disease == "PAV" & pav == 1 ~ "PAV infection",
                               disease == "RPV" & rpv == 1 ~ "RPV infection",
                               disease == "Co" & pav == 1 & rpv == 1 ~ "Co-infection",
                               disease == "Healthy" ~ "Mock inoculation",
                               TRUE ~ NA_character_) %>%
           fct_relevel("Mock inoculation", "PAV infection", "RPV infection"),
         infection_abb = recode(infection, 
                                "PAV infection" = "PAV",
                                "RPV infection" = "RPV",
                                "Co-infection" = "co",
                                "Mock inoculation" = "mock"),
         nitrogen_added = fct_relevel(nitrogen_added, "low", "high"),
         log_chlorophyll = log(chlorophyll),
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
        legend.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_blank())

# palettes
col_pal = c("white", "black")

# panel labels
pan_labs <- tibble(soil = levels(dat2$soil) %>% fct_relevel("sterile", "ambient N", "low N"),
                   label = c("(bold('a'))~sterile~soil",
                             "(bold('b'))~ambient~N~microbes",
                             "(bold('c'))~low~N~microbes",
                             "(bold('d'))~high~N~microbes")) %>%
  mutate(infection_abb = "mock",
         nitrogen_added = "low",
         chlorophyll = 32)

# sample sizes
chlor_samps <- dat2 %>%
  group_by(soil, nitrogen_added, infection_abb) %>%
  count() %>%
  mutate(chlorophyll = 14)

# figure
pdf("output/chlorophyll_figure.pdf", width = 4, height = 4.2)
ggplot(dat2, aes(infection_abb, chlorophyll, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.3), shape = 21) +
  geom_text(data = pan_labs, aes(label = label), hjust = 0, nudge_x = -0.55, parse = T, size = 3) +
  geom_text(data = chlor_samps, aes(label = n), size = 2.5, position = position_dodge(0.4)) +
  facet_wrap(~soil) +
  scale_fill_manual(values = col_pal, name = "Nitrogen supply") +
  xlab("Infection") +
  ylab("Leaf chlorophyll content (SPAD)") +
  coord_cartesian(ylim = c(14, 32.3)) +
  theme_def
dev.off()


#### chlorophyll model ####

# initial fit
chlor_mod1 <- brm(log_chlorophyll ~ soil * N_added * infection, 
                  data = dat3,
                  family = gaussian,
                  prior = c(prior(normal(0, 10), class = Intercept),
                            prior(normal(0, 10), class = b)),
                  iter = 6000, warmup = 1000, chains = 1)
summary(chlor_mod1)

# increase chains
chlor_mod2 <- update(chlor_mod1, chains = 3)

# check model
summary(chlor_mod2)
plot(chlor_mod2)
pp_check(chlor_mod2, nsamples = 50)

# microbes model
chlor_mic_mod1 <- update(chlor_mod2,
                         newdata = dat3,
                         formula = log_chlorophyll ~ microbes * N_added * infection,
                         prior = c(prior(normal(0, 10), class = Intercept),
                                   prior(normal(0, 10), class = b)))

# check model
summary(chlor_mic_mod1)
plot(chlor_mic_mod1)
pp_check(chlor_mic_mod1, nsamples = 50)

# compare with loo
chlor_loo2 <- loo(chlor_mod2, reloo = T)
chlor_loo2 
# all k < 0.7 and elpd_loo is 0.2, but smaller than other SE
# good model fit
# reasonable to do model comparison
chlor_mic_loo1 <- loo(chlor_mic_mod1, reloo = T)
chlor_mic_loo1 
# all k < 0.7 and elpd_loo is 0.2, but less than other SE
# good model fit
# reasonable to do model comparison
loo_compare(chlor_loo2, chlor_mic_loo1)
# microbes model is preferred


#### values for text ####

# soil model
chlor_post2 <- posterior_samples(chlor_mod2) %>%
  rename(soilA_pav_int = "b_soilambientN:infectionPAVinfection",
         soilA_N_int = "b_soilambientN:N_added",
         N_pav_int = "b_N_added:infectionPAVinfection",
         soilA_N_pav_int = "b_soilambientN:N_added:infectionPAVinfection") %>%
  mutate(icp = exp(b_Intercept),
         high_N = exp(b_Intercept + b_N_added),
         N_effect = 100*(high_N - icp)/icp,
         soilA = exp(b_Intercept + b_soilambientN),
         pav_soilA = exp(b_Intercept+ b_soilambientN + b_infectionPAVinfection + soilA_pav_int),
         pav_soilA_effect = 100*(pav_soilA - soilA)/soilA,
         soilA_N = exp(b_Intercept+ b_soilambientN + b_N_added + soilA_N_int),
         pav_soilA_N = exp(b_Intercept+ b_soilambientN + b_N_added + b_infectionPAVinfection + soilA_N_int + soilA_pav_int + N_pav_int + soilA_N_pav_int),
         pav_soilA_N_effect = 100*(pav_soilA_N - soilA_N)/soilA_N)

mean_hdi(chlor_post2$icp)
mean_hdi(chlor_post2$N_effect)
mean_hdi(chlor_post2$pav_soilA_effect)
mean_hdi(chlor_post2$pav_soilA_N_effect)

# microbes model
chlor_mic_post1 <- posterior_samples(chlor_mic_mod1) %>%
  mutate(icp = exp(b_Intercept),
         high_N = exp(b_Intercept + b_N_added),
         N_effect = 100*(high_N - icp)/icp)

mean_hdi(chlor_mic_post1$icp)
mean_hdi(chlor_mic_post1$N_effect)


#### output ####
save(chlor_mod2, file = "output/chlor_bayesian_model_soil.rda")
save(chlor_mic_mod1, file = "output/chlor_bayesian_model_microbes.rda")

write_csv(tidy(summary(chlor_mod2)$fixed), "output/chlor_bayesian_model_soil.csv")
write_csv(tidy(summary(chlor_mic_mod1)$fixed), "output/chlor_bayesian_model_microbes.csv")
