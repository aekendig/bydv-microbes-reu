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

# function to convert logit to probability
logit2prob <- function(x){
  return(exp(x) / (1 + exp(x)))
}


#### edit data ####

# reorganize factor levels
dat2 <- dat %>%
  mutate(soil = fct_relevel(soil, "non-inoculated", "ambient N", "low N"),
         inoculation = case_when(disease %in% c("PAV", "RPV") ~ "single",
                                 disease == "Co" ~ "co-inoculation",
                                 disease == "Healthy" ~ "mock") %>%
           fct_relevel("mock", "single"),
         nitrogen_added = fct_relevel(nitrogen_added, "low", "high"),
         coinfection = case_when(disease == "Co" & pav == 1 & rpv == 1 ~ 1,
                                 TRUE ~ 0),
         microbes = ifelse(soil == "non-inoculated", 0, 1),
         microbes_f = ifelse(microbes == 0, "non-inoculated", "microbes") %>%
           fct_relevel("non-inoculated"))

# separate by virus
pav_dat <- dat2 %>%
  filter(inoc_pav == 1)

rpv_dat <- dat2 %>%
  filter(inoc_rpv == 1)

co_dat <- dat2 %>%
  filter(disease == "Co")

# sample sizes
dat2 %>%
  filter(disease != "Healthy") %>%
  group_by(disease, soil, N_added) %>%
  count() %>%
  data.frame()


#### initial visualization ####

ggplot(pav_dat, aes(inoculation, pav, color = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.3)) +
  facet_wrap(~soil)

ggplot(rpv_dat, aes(inoculation, rpv, color = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.3)) +
  facet_wrap(~soil)

ggplot(co_dat, aes(nitrogen_added, coinfection, color = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.3)) +
  stat_summary(geom = "point", fun = "mean", size = 2, position = position_dodge(0.3)) +
  facet_wrap(~soil)


#### PAV model ####

# initial fit
pav_mod1 <- brm(pav ~ soil * N_added * inoc_rpv, 
                data = pav_dat,
                family = bernoulli,
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 10), class = b)),
                iter = 6000, warmup = 1000, chains = 1)
summary(pav_mod1)

# increase chains
pav_mod2 <- update(pav_mod1, chains = 3)

# check model
summary(pav_mod2)
plot(pav_mod2)
pp_check(pav_mod2, ndraws = 50)

# # microbes model
# pav_mic_mod1 <- update(pav_mod2,
#                        newdata = pav_dat,
#                        formula = pav ~ microbes * N_added * inoc_rpv,
#                        prior = c(prior(normal(0, 10), class = Intercept),
#                                  prior(normal(0, 10), class = b)))
# 
# # check model
# summary(pav_mic_mod1)
# plot(pav_mic_mod1)
# pp_check(pav_mic_mod1, nsamples = 50)
# 
# # compare with loo
# pav_loo2 <- loo(pav_mod2, reloo = T)
# pav_loo2 
# # all k < 0.7 and elpd_loo <= 0.1
# # good model fit
# # reasonable to do model comparison
# pav_mic_loo1 <- loo(pav_mic_mod1)
# pav_mic_loo1 
# # all k < 0.7 and elpd_loo <= 0.1
# # good model fit
# # reasonable to do model comparison
# loo_compare(pav_loo2, pav_mic_loo1)
# # microbes model is slightly preferred


#### RPV model ####

# initial fit
rpv_mod1 <- brm(rpv ~ soil * N_added * inoc_pav, 
                data = rpv_dat,
                family = bernoulli,
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 10), class = b)),
                iter = 6000, warmup = 1000, chains = 1)
summary(rpv_mod1)

# increase chains
rpv_mod2 <- update(rpv_mod1, chains = 3)

# check model
summary(rpv_mod2)
plot(rpv_mod2)
pp_check(rpv_mod2, ndraws = 50)

# # microbes model
# rpv_mic_mod1 <- update(rpv_mod2,
#                        newdata = rpv_dat,
#                        formula = rpv ~ microbes * N_added * inoc_pav,
#                        prior = c(prior(normal(0, 10), class = Intercept),
#                                  prior(normal(0, 10), class = b)))
# 
# # check model
# summary(rpv_mic_mod1)
# plot(rpv_mic_mod1)
# pp_check(rpv_mic_mod1, nsamples = 50)
# 
# # compare with loo
# rpv_loo2 <- loo(rpv_mod2)
# rpv_loo2 
# # all k < 0.7 and elpd_loo <= 0.1
# # good model fit
# # reasonable to do model comparison
# rpv_mic_loo1 <- loo(rpv_mic_mod1, reloo = T)
# rpv_mic_loo1 
# # all k < 0.7
# # elpd_loo > 0.1, but still small compared to other SE
# # good model fit
# # reasonable to do model comparison
# loo_compare(rpv_loo2, rpv_mic_loo1)
# # microbes model is preferred


#### co-infection model ####

# initial fit
co_mod1 <- brm(coinfection ~ soil * N_added, 
               data = co_dat,
               family = bernoulli,
               prior = c(prior(normal(0, 10), class = Intercept),
                         prior(normal(0, 10), class = b)),
               iter = 6000, warmup = 1000, chains = 1)
summary(co_mod1)

# increase chains
co_mod2 <- update(co_mod1, chains = 3)

# check model
summary(co_mod2)
plot(co_mod2)
pp_check(co_mod2, ndraws = 50)

# # microbes model
# co_mic_mod1 <- update(co_mod2,
#                       newdata = co_dat,
#                       formula = coinfection ~ microbes * N_added,
#                       prior = c(prior(normal(0, 10), class = Intercept),
#                                 prior(normal(0, 10), class = b)))
# 
# # check model
# summary(co_mic_mod1)
# plot(co_mic_mod1)
# pp_check(co_mic_mod1, nsamples = 50)
# 
# # compare with loo
# co_loo2 <- loo(co_mod2, reloo = T)
# co_loo2 
# # all k < 0.7 and elpd_loo <= 0.1
# # elpd_loo > 0.1, but still small compared to other SE
# # good model fit
# # reasonable to do model comparison
# co_mic_loo1 <- loo(co_mic_mod1, reloo = T)
# co_mic_loo1 
# # all k < 0.7
# # elpd_loo > 0.1, but still small compared to other SE
# # good model fit
# # reasonable to do model comparison
# loo_compare(co_loo2, co_mic_loo1)
# # microbes model is preferred


#### model output ####
save(pav_mod2, file = "output/pav_bayesian_model_soil.rda")
# save(pav_mic_mod1, file = "output/pav_bayesian_model_microbes.rda")
save(rpv_mod2, file = "output/rpv_bayesian_model_soil.rda")
# save(rpv_mic_mod1, file = "output/rpv_bayesian_model_microbes.rda")
save(co_mod2, file = "output/co_bayesian_model_soil.rda")
# save(co_mic_mod1, file = "output/co_bayesian_model_microbes.rda")

write_csv(rownames_to_column(summary(pav_mod2)$fixed), "output/pav_bayesian_model_soil.csv")
# write_csv(tidy(summary(pav_mic_mod1)$fixed), "output/pav_bayesian_model_microbes.csv")
write_csv(rownames_to_column(summary(rpv_mod2)$fixed), "output/rpv_bayesian_model_soil.csv")
# write_csv(tidy(summary(rpv_mic_mod1)$fixed), "output/rpv_bayesian_model_microbes.csv")
write_csv(rownames_to_column(summary(co_mod2)$fixed), "output/co_bayesian_model_soil.csv")
# write_csv(tidy(summary(co_mic_mod1)$fixed), "output/co_bayesian_model_microbes.csv")


#### visualize ####

# posterior sample function
post_fun <- function(mod){
  
  out <- as_draws_df(mod) %>%
    rename_with(~ str_replace_all(., ":", "_")) %>%
    rename_with(~ str_replace_all(., "rpv", "co")) %>%
    rename_with(~ str_replace_all(., "pav", "co")) %>%
    transmute(noninoculated_0_0 = logit2prob(b_Intercept),
              noninoculated_1_0 = logit2prob(b_Intercept + b_N_added),
              noninoculated_0_1 = logit2prob(b_Intercept + b_inoc_co),
              noninoculated_1_1 = logit2prob(b_Intercept + b_N_added + b_inoc_co + b_N_added_inoc_co),
              ambientN_0_0 = logit2prob(b_Intercept + b_soilambientN),
              ambientN_1_0 = logit2prob(b_Intercept + b_soilambientN + b_N_added + b_soilambientN_N_added),
              ambientN_0_1 = logit2prob(b_Intercept + b_soilambientN + b_inoc_co + b_soilambientN_inoc_co),
              ambientN_1_1 = logit2prob(b_Intercept + b_soilambientN + b_N_added + b_inoc_co + b_N_added_inoc_co + b_soilambientN_N_added + b_soilambientN_inoc_co + b_soilambientN_N_added_inoc_co),
              lowN_0_0 = logit2prob(b_Intercept + b_soillowN),
              lowN_1_0 = logit2prob(b_Intercept + b_soillowN + b_N_added + b_soillowN_N_added),
              lowN_0_1 = logit2prob(b_Intercept + b_soillowN + b_inoc_co + b_soillowN_inoc_co),
              lowN_1_1 = logit2prob(b_Intercept + b_soillowN + b_N_added + b_inoc_co + b_N_added_inoc_co + b_soillowN_N_added + b_soillowN_inoc_co + b_soillowN_N_added_inoc_co),
              highN_0_0 = logit2prob(b_Intercept + b_soilhighN),
              highN_1_0 = logit2prob(b_Intercept + b_soilhighN + b_N_added + b_soilhighN_N_added),
              highN_0_1 = logit2prob(b_Intercept + b_soilhighN + b_inoc_co + b_soilhighN_inoc_co),
              highN_1_1 = logit2prob(b_Intercept + b_soilhighN + b_N_added + b_inoc_co + b_N_added_inoc_co + b_soilhighN_N_added + b_soilhighN_inoc_co + b_soilhighN_N_added_inoc_co)) %>%
    pivot_longer(cols = everything(),
                 names_to = "treatment",
                 values_to = "values") %>%
    rowwise() %>%
    mutate(soil = str_split(treatment, "_")[[1]][1],
           N_added = str_split(treatment, "_")[[1]][2] %>% as.double(),
           inoc = str_split(treatment, "_")[[1]][3] %>% as.double()) %>%
    ungroup() %>%
    mutate(soil = str_replace(soil, "N", " N"),
           soil = fct_relevel(soil, "noninoculated", "ambient N", "low N") %>%
             fct_recode("non-inoculated" = "noninoculated"),
           nitrogen_added = if_else(N_added == 1, "high", "low"),
           nitrogen_added = fct_relevel(nitrogen_added, "low"),
           inoculation = if_else(inoc == 1, "co-inoculation", "single"),
           inoculation = fct_relevel(inoculation, "single"))
  
  return(out)
}

# posterior samples
pav_post <- post_fun(pav_mod2) %>%
  rename(pav = values)

rpv_post <- post_fun(rpv_mod2) %>%
  rename(rpv = values)

co_post <- as_draws_df(co_mod2) %>%
  rename_with(~ str_replace_all(., ":", "_")) %>%
  transmute(noninoculated_0 = logit2prob(b_Intercept),
            noninoculated_1 = logit2prob(b_Intercept + b_N_added),
            ambientN_0 = logit2prob(b_Intercept + b_soilambientN),
            ambientN_1 = logit2prob(b_Intercept + b_soilambientN + b_N_added + b_soilambientN_N_added),
            lowN_0 = logit2prob(b_Intercept + b_soillowN),
            lowN_1 = logit2prob(b_Intercept + b_soillowN + b_N_added + b_soillowN_N_added),
            highN_0 = logit2prob(b_Intercept + b_soilhighN),
            highN_1 = logit2prob(b_Intercept + b_soilhighN + b_N_added + b_soilhighN_N_added)) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "coinfection") %>%
  rowwise() %>%
  mutate(soil = str_split(treatment, "_")[[1]][1],
         N_added = str_split(treatment, "_")[[1]][2] %>% as.double()) %>%
  ungroup() %>%
  mutate(soil = str_replace(soil, "N", " N"),
         soil = fct_relevel(soil, "noninoculated", "ambient N", "low N") %>%
           fct_recode("non-inoculated" = "noninoculated"),
         nitrogen_added = if_else(N_added == 1, "high", "low"),
         nitrogen_added = fct_relevel(nitrogen_added, "low"))

# theme
theme_def <- theme_bw() +
  theme(axis.text.y = element_text(size = 8, color="black"),
        axis.text.x = element_text(size = 8, color="black"),
        axis.title.y = element_text(size = 10, color="black"),
        axis.title.x = element_text(size = 10, color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin = margin(-10, -10, -10, -10),
        legend.background = element_blank(),
        # legend.position = "bottom",
        # legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_blank())

col_pal <- viridis_pal(direction = -1)(4)

# PAV infection prevalence
tiff("output/Figure_1.tiff", width = 160, height = 90, units = "mm", res = 300, compression = "lzw")
ggplot(pav_post, aes(soil, pav, shape = nitrogen_added, color = inoculation, fill = inoculation, group = interaction(nitrogen_added, inoculation))) +
  geom_point(data = pav_dat, size = 0.75, alpha = 0.5, position = position_jitterdodge(0.05, 0.05, 0.5)) +
  stat_pointinterval(.width = 0.95, point_interval = mean_hdi, position = position_dodge(0.5), alpha = 0.7, point_size = 2.5, interval_size = 0.75) +
  scale_color_manual(values = col_pal[2:3], name = "Virus inoculation",
                     labels = c("PAV", "co-inoculation")) +
  scale_fill_manual(values = col_pal[2:3], name = "Virus inoculation",
                    labels = c("PAV", "co-inoculation")) +
  scale_shape(name = "Nitrogen supply") +
  labs(x = "Field soil treatment", y = "BYDV-PAV incidence") +
  theme_def
dev.off()

# RPV infection prevalence
tiff("output/Figure_2.tiff", width = 160, height = 90, units = "mm", res = 300, compression = "lzw")
ggplot(rpv_post, aes(soil, rpv, shape = nitrogen_added, color = inoculation, fill = inoculation, group = interaction(nitrogen_added, inoculation))) +
  geom_point(data = rpv_dat, size = 0.75, alpha = 0.5, position = position_jitterdodge(0.05, 0.05, 0.5)) +
  stat_pointinterval(.width = 0.95, point_interval = mean_hdi, position = position_dodge(0.5), alpha = 0.7, point_size = 2.5, interval_size = 0.75) +
  scale_color_manual(values = col_pal[c(4,3)], name = "Virus inoculation",
                     labels = c("RPV", "co-inoculation")) +
  scale_fill_manual(values = col_pal[c(4,3)], name = "Virus inoculation",
                    labels = c("RPV", "co-inoculation")) +
  scale_shape(name = "Nitrogen supply") +
  labs(x = "Field soil treatment", y = "CYDV-RPV incidence") +
  theme_def +
  guides(shape = guide_legend(order = 2), col = guide_legend(order = 1), 
         fill = guide_legend(order = 1))
dev.off()

# Co-infection prevalence
tiff("output/Figure_3.tiff", width = 140, height = 90, units = "mm", res = 300, compression = "lzw")
ggplot(co_post, aes(x = soil, y = coinfection, shape = nitrogen_added)) +
  geom_point(data = co_dat, size = 0.75, alpha = 0.5, position = position_jitterdodge(0.05, 0.05, 0.25), color = col_pal[3]) +
  stat_pointinterval(.width = 0.95, point_interval = mean_hdi, position = position_dodge(0.25), alpha = 0.7, point_size = 2.5, interval_size = 0.75, color = col_pal[3], fill = col_pal[3]) +
  scale_shape(name = "Nitrogen supply") +
  labs(x = "Field soil treatment", y = "Co-infection incidence") +
  theme_def
dev.off()


#### values for text ####

# PAV incidence
pav_post2 <- posterior_samples(pav_mod2) %>%
  rename(b_N_rpv_int = "b_N_added:inoc_rpv",
         b_N_soilA_int = "b_soilambientN:N_added",
         b_N_soilL_int = "b_soillowN:N_added",
         b_N_soilH_int = "b_soilhighN:N_added",
         b_soilA_rpv_int = "b_soilambientN:inoc_rpv",
         b_soilL_rpv_int = "b_soillowN:inoc_rpv",
         b_soilH_rpv_int = "b_soilhighN:inoc_rpv",
         b_N_soilA_rpv_int = "b_soilambientN:N_added:inoc_rpv",
         b_N_soilL_rpv_int = "b_soillowN:N_added:inoc_rpv") %>%
  mutate(icp = exp(b_Intercept)/(1 + exp(b_Intercept)),
         high_N = exp(b_Intercept + b_N_added)/(1 + exp(b_Intercept + b_N_added)),
         N_effect = 100 * (high_N - icp)/icp,
         coinoc = exp(b_Intercept + b_inoc_rpv)/(1 + exp(b_Intercept + b_inoc_rpv)),
         co_effect = 100 * (coinoc - icp)/icp,
         coinoc_N = exp(b_Intercept + b_inoc_rpv + b_N_added + b_N_rpv_int)/(1 + exp(b_Intercept + b_inoc_rpv + b_N_added + b_N_rpv_int)),
         N_co_effect = 100 * (coinoc_N - coinoc)/coinoc,
         soilA = exp(b_Intercept + b_soilambientN)/(1 + exp(b_Intercept + b_soilambientN)),
         soilL = exp(b_Intercept + b_soillowN)/(1 + exp(b_Intercept + b_soillowN)),
         soilH = exp(b_Intercept + b_soilhighN)/(1 + exp(b_Intercept + b_soilhighN)),
         soilL_effect = 100 * (soilL - icp)/icp,
         soilA_N = exp(b_Intercept + b_soilambientN + b_N_added + b_N_soilA_int)/(1 + exp(b_Intercept + b_soilambientN + b_N_added + b_N_soilA_int)),
         soilL_N = exp(b_Intercept + b_soillowN + b_N_added + b_N_soilL_int)/(1 + exp(b_Intercept + b_soillowN + b_N_added + b_N_soilL_int)),
         soilH_N = exp(b_Intercept + b_soilhighN + b_N_added + b_N_soilH_int)/(1 + exp(b_Intercept + b_soilhighN + b_N_added + b_N_soilH_int)),
         N_soilA_effect = 100 * (soilA_N - soilA) / soilA,
         N_soilL_effect = 100 * (soilL_N - soilL) / soilL,
         N_soilH_effect = 100 * (soilH_N - soilH) / soilH,
         soilL_N_effect = 100 * (soilL_N - high_N) / high_N,
         soilA_co = exp(b_Intercept + b_soilambientN + b_inoc_rpv + b_soilA_rpv_int)/(1 + exp(b_Intercept + b_soilambientN + b_inoc_rpv + b_soilA_rpv_int)),
         soilL_co = exp(b_Intercept + b_soillowN + b_inoc_rpv + b_soilL_rpv_int)/(1 + exp(b_Intercept + b_soillowN + b_inoc_rpv + b_soilL_rpv_int)),
         soilH_co = exp(b_Intercept + b_soilhighN + b_inoc_rpv + b_soilH_rpv_int)/(1 + exp(b_Intercept + b_soilhighN + b_inoc_rpv + b_soilH_rpv_int)),
         co_soilA_effect = 100 * (soilA_co - soilA) / soilA,
         co_soilL_effect = 100 * (soilL_co - soilL) / soilL,
         co_soilH_effect = 100 * (soilH_co - soilH) / soilH,
         soilA_co_N = exp(b_Intercept + b_soilambientN + b_inoc_rpv + b_N_added + b_soilA_rpv_int + b_N_rpv_int + b_N_soilA_int + b_N_soilA_rpv_int)/(1 + exp(b_Intercept + b_soilambientN + b_inoc_rpv + b_N_added + b_soilA_rpv_int + b_N_rpv_int + b_N_soilA_int + b_N_soilA_rpv_int)),
         soilL_co_N = exp(b_Intercept + b_soillowN + b_inoc_rpv + b_N_added + b_soilL_rpv_int + b_N_rpv_int + b_N_soilL_int + b_N_soilL_rpv_int)/(1 + exp(b_Intercept + b_soillowN + b_inoc_rpv + b_N_added + b_soilL_rpv_int + b_N_rpv_int + b_N_soilL_int + b_N_soilL_rpv_int)),
         N_co_soilA_effect = 100 * (soilA_co_N - soilA_co) / soilA_co,
         N_co_soilL_effect = 100 * (soilL_co_N - soilL_co) / soilL_co)

# sterile soil paragraph
mean_hdi(pav_post2$icp)
mean_hdi(pav_post2$high_N)
mean_hdi(pav_post2$N_effect)

mean_hdi(pav_post2$coinoc)
mean_hdi(pav_post2$co_effect)

mean_hdi(pav_post2$N_co_effect)

# soil microbes paragraph
mean_hdi(pav_post2$soilL)
mean_hdi(pav_post2$soilL_effect)

mean_hdci(pav_post2$N_soilL_effect)
mean_hdi(pav_post2$N_soilA_effect)
mean_hdci(pav_post2$N_soilH_effect)

mean_hdi(pav_post2$co_soilL_effect)
mean_hdi(pav_post2$co_soilA_effect)

mean_hdci(pav_post2$soilH)
mean_hdi(pav_post2$soilH_co)
mean_hdi(pav_post2$co_soilH_effect)

# RPV incidence
rpv_post2 <- posterior_samples(rpv_mod2) %>%
  mutate(icp = exp(b_Intercept)/(1 + exp(b_Intercept)),
         coinoc = exp(b_Intercept + b_inoc_pav)/(1 + exp(b_Intercept + b_inoc_pav)),
         co_effect = 100 * (coinoc - icp)/icp)

# sterile soil paragraph
mean_hdi(rpv_post2$icp)
mean_hdi(rpv_post2$coinoc)
mean_hdi(rpv_post2$co_effect)

# co incidence
co_post2 <- posterior_samples(co_mod2) %>%
  mutate(icp = exp(b_Intercept)/(1 + exp(b_Intercept)))

mean_hdi(co_post2$icp)

# PAV incidence microbes model
pav_mic_post1 <- posterior_samples(pav_mic_mod1) %>%
  rename(b_N_rpv_int = "b_N_added:inoc_rpv",
         b_N_mic_int = "b_microbes:N_added",
         b_mic_rpv_int = "b_microbes:inoc_rpv",
         b_N_mic_rpv_int = "b_microbes:N_added:inoc_rpv") %>%
  mutate(icp = exp(b_Intercept)/(1 + exp(b_Intercept)),
         high_N = exp(b_Intercept + b_N_added)/(1 + exp(b_Intercept + b_N_added)),
         N_effect = high_N - icp,
         coinoc = exp(b_Intercept + b_inoc_rpv)/(1 + exp(b_Intercept + b_inoc_rpv)),
         co_effect = coinoc - icp,
         coinoc_N = exp(b_Intercept + b_inoc_rpv + b_N_added + b_N_rpv_int)/(1 + exp(b_Intercept + b_inoc_rpv + b_N_added + b_N_rpv_int)),
         N_co_effect = coinoc_N - coinoc,
         microbes = exp(b_Intercept + b_microbes)/(1 + exp(b_Intercept + b_microbes)),
         mic_effect = microbes - icp,
         microbes_N = exp(b_Intercept + b_microbes + b_N_added + b_N_mic_int)/(1 + exp(b_Intercept + b_microbes + b_N_added + b_N_mic_int)),
         N_mic_effect = microbes_N - microbes,
         microbes_co = exp(b_Intercept + b_microbes + b_inoc_rpv + b_mic_rpv_int)/(1 + exp(b_Intercept + b_microbes + b_inoc_rpv + b_mic_rpv_int)),
         co_mic_effect = microbes_co - microbes,
         microbes_co_N = exp(b_Intercept + b_microbes + b_inoc_rpv + b_N_added + b_mic_rpv_int + b_N_mic_int + b_N_rpv_int + b_N_mic_rpv_int)/(1 + exp(b_Intercept + b_microbes + b_inoc_rpv + b_N_added + b_mic_rpv_int + b_N_mic_int + b_N_rpv_int + b_N_mic_rpv_int)),
         N_co_mic_effect = microbes_co_N - microbes_co,
         mic_co_effect = microbes_co - coinoc)

ggplot(pav_mic_post1, aes(x = N_effect)) +
  geom_histogram(bins = 100)

mean_hdci(pav_mic_post1$icp)
mean_hdi(pav_mic_post1$N_effect)
mean_hdi(pav_mic_post1$co_effect)
mean_hdi(pav_mic_post1$N_co_effect)
mean_hdi(pav_mic_post1$mic_effect)
mean_hdi(pav_mic_post1$N_mic_effect)
mean_hdi(pav_mic_post1$co_mic_effect)
mean_hdi(pav_mic_post1$N_co_mic_effect)
mean_hdi(pav_mic_post1$mic_co_effect)

# RPV incidence microbes model
rpv_mic_post1 <- posterior_samples(rpv_mic_mod1) %>%
  mutate(icp = exp(b_Intercept)/(1 + exp(b_Intercept)),
         coinoc = exp(b_Intercept + b_inoc_pav)/(1 + exp(b_Intercept + b_inoc_pav)),
         co_effect = coinoc - icp)

mean_hdi(rpv_mic_post1$icp)
mean_hdi(rpv_mic_post1$co_effect)

# co incidence microbes model
co_mic_post1 <- posterior_samples(co_mic_mod1) %>%
  mutate(icp = exp(b_Intercept)/(1 + exp(b_Intercept)),
         microbes = exp(b_Intercept + b_microbes)/(1 + exp(b_Intercept + b_microbes)),
         mic_effect = microbes - icp)

mean_hdi(co_mic_post1$icp)
mean_hdi(co_mic_post1$mic_effect)


