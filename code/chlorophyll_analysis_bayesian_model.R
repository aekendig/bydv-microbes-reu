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
library(scales)

# import data
dat <- read_csv("intermediate-data/bydv_microbes_data_rounded_up.csv")


#### edit data ####

# reorganize factor levels
# select for successful infection
dat2 <- dat %>%
  rowwise() %>%
  mutate(chlorophyll = mean(c(chlorophyll_1, chlorophyll_2, chlorophyll_3), na.rm = T)) %>%
  ungroup() %>%
  mutate(soil = fct_relevel(soil, "non-inoculated", "ambient N", "low N"),
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
         microbes = ifelse(soil == "non-inoculated", 0, 1),
         microbes_f = ifelse(microbes == 0, "non-inoculated", "microbes") %>%
           fct_relevel("non-inoculated")) %>%
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


#### initial visualization ####

ggplot(dat2, aes(infection_abb, chlorophyll, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.3), shape = 21)


#### chlorophyll model ####

# initial fit
chlor_mod1 <- brm(log_chlorophyll ~ soil * N_added * infection_abb, 
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

# # microbes model
# chlor_mic_mod1 <- update(chlor_mod2,
#                          newdata = dat3,
#                          formula = log_chlorophyll ~ microbes * N_added * infection,
#                          prior = c(prior(normal(0, 10), class = Intercept),
#                                    prior(normal(0, 10), class = b)))
# 
# # check model
# summary(chlor_mic_mod1)
# plot(chlor_mic_mod1)
# pp_check(chlor_mic_mod1, nsamples = 50)
# 
# # compare with loo
# chlor_loo2 <- loo(chlor_mod2, reloo = T)
# chlor_loo2 
# # all k < 0.7 and elpd_loo is 0.2, but smaller than other SE
# # good model fit
# # reasonable to do model comparison
# chlor_mic_loo1 <- loo(chlor_mic_mod1, reloo = T)
# chlor_mic_loo1 
# # all k < 0.7 and elpd_loo is 0.2, but less than other SE
# # good model fit
# # reasonable to do model comparison
# loo_compare(chlor_loo2, chlor_mic_loo1)
# # microbes model is preferred

save(chlor_mod2, file = "output/chlor_bayesian_model_soil.rda")
# save(chlor_mic_mod1, file = "output/chlor_bayesian_model_microbes.rda")

write_csv(tidy(summary(chlor_mod2)$fixed), "output/chlor_bayesian_model_soil.csv")
# write_csv(tidy(summary(chlor_mic_mod1)$fixed), "output/chlor_bayesian_model_microbes.csv")


#### visualize ####

# posterior dist
chlor_post <- as_draws_df(chlor_mod2) %>%
  rename_with(~ str_replace_all(., ":", "_")) %>%
  transmute(noninoculated_0_mock = exp(b_Intercept),
            noninoculated_1_mock = exp(b_Intercept + b_N_added),
            noninoculated_0_PAV = exp(b_Intercept + b_infection_abbPAV),
            noninoculated_1_PAV = exp(b_Intercept + b_N_added + b_infection_abbPAV + b_N_added_infection_abbPAV),
            noninoculated_0_RPV = exp(b_Intercept + b_infection_abbRPV),
            noninoculated_1_RPV = exp(b_Intercept + b_N_added + b_infection_abbRPV + b_N_added_infection_abbRPV),
            ambientN_0_mock = exp(b_Intercept + b_soilambientN),
            ambientN_1_mock = exp(b_Intercept + b_soilambientN + b_N_added + b_soilambientN_N_added),
            ambientN_0_PAV = exp(b_Intercept + b_soilambientN + b_infection_abbPAV + b_soilambientN_infection_abbPAV),
            ambientN_1_PAV = exp(b_Intercept + b_soilambientN + b_N_added + b_infection_abbPAV + b_N_added_infection_abbPAV + b_soilambientN_N_added + b_soilambientN_infection_abbPAV + b_soilambientN_N_added_infection_abbPAV),
            ambientN_0_RPV = exp(b_Intercept + b_soilambientN + b_infection_abbRPV + b_soilambientN_infection_abbRPV),
            ambientN_1_RPV = exp(b_Intercept + b_soilambientN + b_N_added + b_infection_abbRPV + b_N_added_infection_abbRPV + b_soilambientN_N_added + b_soilambientN_infection_abbRPV + b_soilambientN_N_added_infection_abbRPV),
            lowN_0_mock = exp(b_Intercept + b_soillowN),
            lowN_1_mock = exp(b_Intercept + b_soillowN + b_N_added + b_soillowN_N_added),
            lowN_0_PAV = exp(b_Intercept + b_soillowN + b_infection_abbPAV + b_soillowN_infection_abbPAV),
            lowN_1_PAV = exp(b_Intercept + b_soillowN + b_N_added + b_infection_abbPAV + b_N_added_infection_abbPAV + b_soillowN_N_added + b_soillowN_infection_abbPAV + b_soillowN_N_added_infection_abbPAV),
            lowN_0_RPV = exp(b_Intercept + b_soillowN + b_infection_abbRPV + b_soillowN_infection_abbRPV),
            lowN_1_RPV = exp(b_Intercept + b_soillowN + b_N_added + b_infection_abbRPV + b_N_added_infection_abbRPV + b_soillowN_N_added + b_soillowN_infection_abbRPV + b_soillowN_N_added_infection_abbRPV),
            highN_0_mock = exp(b_Intercept + b_soilhighN),
            highN_1_mock = exp(b_Intercept + b_soilhighN + b_N_added + b_soilhighN_N_added),
            highN_0_PAV = exp(b_Intercept + b_soilhighN + b_infection_abbPAV + b_soilhighN_infection_abbPAV),
            highN_1_PAV = exp(b_Intercept + b_soilhighN + b_N_added + b_infection_abbPAV + b_N_added_infection_abbPAV + b_soilhighN_N_added + b_soilhighN_infection_abbPAV + b_soilhighN_N_added_infection_abbPAV),
            highN_0_RPV = exp(b_Intercept + b_soilhighN + b_infection_abbRPV + b_soilhighN_infection_abbRPV),
            highN_1_RPV = exp(b_Intercept + b_soilhighN + b_N_added + b_infection_abbRPV + b_N_added_infection_abbRPV + b_soilhighN_N_added + b_soilhighN_infection_abbRPV + b_soilhighN_N_added_infection_abbRPV)) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "chlorophyll") %>%
  rowwise() %>%
  mutate(soil = str_split(treatment, "_")[[1]][1],
         N_added = str_split(treatment, "_")[[1]][2] %>% as.double(),
         infection_abb = str_split(treatment, "_")[[1]][3]) %>%
  ungroup() %>%
  mutate(soil = str_replace(soil, "N", " N"),
         soil = fct_relevel(soil, "noninoculated", "ambient N", "low N") %>%
           fct_recode("non-inoculated" = "noninoculated"),
         nitrogen_added = if_else(N_added == 1, "high", "low"),
         nitrogen_added = fct_relevel(nitrogen_added, "low"))

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
        strip.background = element_blank(),
        strip.text = element_blank())

col_pal <- viridis_pal(direction = -1)(4)

# figure
tiff("output/Figure_5.tiff", width = 180, height = 90, units = "mm", res = 300, compression = "lzw")
ggplot(chlor_post, aes(soil, chlorophyll, shape = nitrogen_added, color = infection_abb, fill = infection_abb, group = interaction(nitrogen_added, infection_abb))) +
  geom_point(data = dat3, size = 0.75, alpha = 0.5, position = position_jitterdodge(0.05, 0.05, 0.6)) +
  stat_pointinterval(.width = 0.95, position = position_dodge(0.6), alpha = 0.7, point_size = 2.5, interval_size = 0.75) +
  scale_color_manual(values = col_pal[c(1,2,4)], name = "Virus infection") +
  scale_fill_manual(values = col_pal[c(1,2,4)], name = "Virus infection") +
  scale_shape(name = "Nitrogen supply") +
  labs(x = "Field soil treatment", y = "Leaf chlorophyll content (SPAD)") +
  theme_def
dev.off()



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
mean_hdi(chlor_post2$high_N)
mean_hdi(chlor_post2$N_effect)

mean_hdi(chlor_post2$soilA)
mean_hdi(chlor_post2$pav_soilA)
mean_hdi(chlor_post2$pav_soilA_effect)

mean_hdi(chlor_post2$pav_soilA_N_effect)

# microbes model
chlor_mic_post1 <- posterior_samples(chlor_mic_mod1) %>%
  mutate(icp = exp(b_Intercept),
         high_N = exp(b_Intercept + b_N_added),
         N_effect = 100*(high_N - icp)/icp)

mean_hdi(chlor_mic_post1$icp)
mean_hdi(chlor_mic_post1$N_effect)
