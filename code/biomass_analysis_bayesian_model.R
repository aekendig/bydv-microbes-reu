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
         log_biomass = log(biomass),
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

ggplot(dat2, aes(infection_abb, biomass, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.3), shape = 21)


#### biomass model ####

# initial fit
bio_mod1 <- brm(log_biomass ~ soil * N_added * infection_abb, 
                data = dat3,
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

# # microbes model
# bio_mic_mod1 <- update(bio_mod2,
#                        newdata = dat3,
#                        formula = log_biomass ~ microbes * N_added * infection,
#                        prior = c(prior(normal(0, 10), class = Intercept),
#                                  prior(normal(0, 10), class = b)))
# 
# # check model
# summary(bio_mic_mod1)
# plot(bio_mic_mod1)
# pp_check(bio_mic_mod1, nsamples = 50)
# 
# # compare with loo
# bio_loo2 <- loo(bio_mod2)
# bio_loo2 
# # all k < 0.7 and elpd_loo <= 0.1
# # good model fit
# # reasonable to do model comparison
# bio_mic_loo1 <- loo(bio_mic_mod1)
# bio_mic_loo1 
# # all k < 0.7 and elpd_loo <= 0.1
# # good model fit
# # reasonable to do model comparison
# loo_compare(bio_loo2, bio_mic_loo1)
# # microbes model is preferred

# save models
save(bio_mod2, file = "output/bio_bayesian_model_soil.rda")
# save(bio_mic_mod1, file = "output/bio_bayesian_model_microbes.rda")

write_csv(tidy(summary(bio_mod2)$fixed), "output/bio_bayesian_model_soil.csv")
# write_csv(tidy(summary(bio_mic_mod1)$fixed), "output/bio_bayesian_model_microbes.csv")


#### visualize ####

# posterior dist
bio_post <- as_draws_df(bio_mod2) %>%
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
               values_to = "biomass") %>%
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
tiff("output/Figure_4.tiff", width = 180, height = 90, units = "mm", res = 300, compression = "lzw")
ggplot(bio_post, aes(soil, biomass, shape = nitrogen_added, color = infection_abb, fill = infection_abb, group = interaction(nitrogen_added, infection_abb))) +
  geom_point(data = dat3, size = 0.75, alpha = 0.5, position = position_jitterdodge(0.05, 0.05, 0.6)) +
  stat_pointinterval(.width = 0.95, position = position_dodge(0.6), alpha = 0.7, point_size = 2.5, interval_size = 0.75) +
  scale_color_manual(values = col_pal[c(1,2,4)], name = "Virus infection") +
  scale_fill_manual(values = col_pal[c(1,2,4)], name = "Virus infection") +
  scale_shape(name = "Nitrogen supply") +
  labs(x = "Field soil treatment", y = "Biomass (g)") +
  theme_def
dev.off()


#### values for text ####

# soil model
bio_post2 <- posterior_samples(bio_mod2) %>%
  mutate(icp = exp(b_Intercept),
         high_N = exp(b_Intercept + b_N_added),
         N_coef = 100*(exp(b_N_added) - 1),
         N_effect = 100*(high_N - icp)/icp,
         log_N_effect = 100*(b_Intercept + b_N_added - b_Intercept)/b_Intercept)

mean_hdi(bio_post2$icp)
mean_hdi(bio_post2$high_N)
mean_hdi(bio_post2$N_coef)
mean_hdi(bio_post2$N_effect)
mean_hdi(bio_post2$log_N_effect)

# microbes model
bio_mic_post1 <- posterior_samples(bio_mic_mod1) %>%
  mutate(icp = exp(b_Intercept),
         high_N = exp(b_Intercept + b_N_added),
         N_effect = 100*(high_N - icp)/icp)

mean_hdi(bio_mic_post1$icp)
mean_hdi(bio_mic_post1$N_effect)


#### power analysis ####

# source: https://solomonkurz.netlify.app/post/bayesian-power-analysis-part-i/

# use model output
summary(bio_mod2)

# parameters
mu_h <- -1.61
mu_i <- -1.61 - 0.23
sig <- 0.53
n_h <- nrow(dat3 %>%
              filter(soil == "sterile" & N_added == 0 & infection == "Mock inoculation"))
n_i <- nrow(dat3 %>%
              filter(soil == "sterile" & N_added == 0 & infection == "RPV infection"))
n_sim <- 1000
n_h2 <- n_h * 2
n_i2 <- n_i * 2
n_h3 <- n_h * 9
n_i3 <- n_i * 9
n_h4 <- n_h * 10
n_i4 <- n_i * 10

# first simulation
set.seed(1)

d <- tibble(group = c(rep("healthy", n_h), rep("infected", n_i))) %>%
  mutate(infection = ifelse(group == "healthy", 0, 1),
         y = ifelse(group == "healthy",
                    rnorm(n_h, mean = mu_h, sd = sig),
                    rnorm(n_i, mean = mu_i, sd = sig)))

fit <- brm(data = d,
           family = gaussian,
           y ~ infection,
           prior = c(prior(normal(-1.6, 2), class = Intercept),
                     prior(normal(0, 2), class = b)),
           seed = 1)

plot(fit)
summary(fit)
summary(fit)$fixed[2, c(1, 3, 4)]

# iterative function
sim_d_and_fit <- function(seed, n_c, n_t){
  
  set.seed(seed)
  
  d <- tibble(group = c(rep("healthy", n_c), rep("infected", n_t))) %>%
    mutate(infection = ifelse(group == "healthy", 0, 1),
           y = ifelse(group == "healthy",
                      rnorm(n_c, mean = mu_h, sd = sig),
                      rnorm(n_t, mean = mu_i, sd = sig)))
  
  fit_new <- update(fit,
                    newdata = d, 
                    seed = seed)
  
  return(summary(fit_new)$fixed[2, c(1, 3, 4)])
}

# run simulations
t1 <- Sys.time()

s <- tibble(seed = 1:n_sim) %>% 
  mutate(tidy = map(seed, sim_d_and_fit, n_c = n_h, n_t = n_i))

t2 <- Sys.time()

# format s
s2 <- s %>%
  unnest_wider(tidy) %>%
  rename(estimate = Estimate,
         lower = `l-95% CI`,
         upper = `u-95% CI`)

# visualize
ggplot(s2, aes(x = seed, y = estimate, ymin = lower, ymax = upper)) +
  geom_pointrange(fatten = 1/2) +
  geom_hline(yintercept = c(0, -0.23), color = "red") +
  labs(x = "seed (i.e., simulation index)",
       y = expression(beta[1]))

# summarize
s2 %>%
  mutate(sig = case_when(lower < 0 & upper < 0 ~ 1,
                         lower > 0 & upper > 0 ~ 1,
                         TRUE ~ 0)) %>%
  summarise(power = sum(sig) / n_sim)

# run simulations for 2x sample size
v <- tibble(seed = 1:n_sim) %>% 
  mutate(tidy = map(seed, sim_d_and_fit, n_c = n_h2, n_t = n_i2))

v2 <- v %>%
  unnest_wider(tidy) %>%
  rename(estimate = Estimate,
         lower = `l-95% CI`,
         upper = `u-95% CI`)

v2 %>%
  mutate(sig = case_when(lower < 0 & upper < 0 ~ 1,
                         lower > 0 & upper > 0 ~ 1,
                         TRUE ~ 0)) %>%
  summarise(power = sum(sig) / n_sim)

# run simulations for 9x sample size
z <- tibble(seed = 1:n_sim) %>% 
  mutate(tidy = map(seed, sim_d_and_fit, n_c = n_h3, n_t = n_i3))

z2 <- z %>%
  unnest_wider(tidy) %>%
  rename(estimate = Estimate,
         lower = `l-95% CI`,
         upper = `u-95% CI`)

z2 %>%
  mutate(sig = case_when(lower < 0 & upper < 0 ~ 1,
                         lower > 0 & upper > 0 ~ 1,
                         TRUE ~ 0)) %>%
  summarise(power = sum(sig) / n_sim)

# run simulations for 10x sample size
y <- tibble(seed = 1:n_sim) %>% 
  mutate(tidy = map(seed, sim_d_and_fit, n_c = n_h4, n_t = n_i4))

y2 <- y %>%
  unnest_wider(tidy) %>%
  rename(estimate = Estimate,
         lower = `l-95% CI`,
         upper = `u-95% CI`)

y2 %>%
  mutate(sig = case_when(lower < 0 & upper < 0 ~ 1,
                         lower > 0 & upper > 0 ~ 1,
                         TRUE ~ 0)) %>%
  summarise(power = sum(sig) / n_sim)

write_csv(s2, "output/simulated_datasets_1x_samp_size.csv")
write_csv(v2, "output/simulated_datasets_2x_samp_size.csv")
write_csv(z2, "output/simulated_datasets_9x_samp_size.csv")
write_csv(y2, "output/simulated_datasets_10x_samp_size.csv")
