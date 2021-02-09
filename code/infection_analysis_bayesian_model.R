##### info ####

# goal: analyze data when any visible PCR band is counted as indicators of infection (rounding up)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(cowplot)
library(broom)
library(brms)
library(tidybayes)

# import data
dat <- read_csv("intermediate-data/bydv_microbes_data_rounded_up.csv")


#### edit data ####

# reorganize factor levels
dat2 <- dat %>%
  mutate(soil = fct_relevel(soil, "sterile", "ambient N", "low N"),
         inoculation = case_when(disease %in% c("PAV", "RPV") ~ "Single inoculation",
                                 disease == "Co" ~ "Co-inoculation",
                                 disease == "Healthy" ~ "Mock inoculation") %>%
           fct_relevel("Mock inoculation", "Single inoculation"),
         nitrogen_added = fct_relevel(nitrogen_added, "low", "high"),
         coinfection = case_when(disease == "Co" & pav == 1 & rpv == 1 ~ 1,
                                 TRUE ~ 0),
         microbes = ifelse(soil == "sterile", 0, 1),
         microbes_f = ifelse(microbes == 0, "sterile", "microbes") %>%
           fct_relevel("sterile"))

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
pav_lab <- tibble(inoculation = c("Single inoculation", "Co-inoculation") %>% fct_relevel("Single inoculation"),
                  label = c("(a)", "(b)")) %>%
  mutate(microbes_f = "sterile",
         pav = 1,
         nitrogen_added = "low")

rpv_lab <- pav_lab %>%
  mutate(label = recode(label, "(a)" = "(c)", "(b)" = "(d)"),
         rpv = 1)

# PAV infection prevalence
pav_fig  <- ggplot(pav_dat, aes(microbes_f, pav, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2), shape = 21) +
  geom_text(data = pav_lab, aes(label = label), nudge_x = -0.45, fontface = "bold") +
  facet_wrap(~ inoculation) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  xlab("Microbe inoculation") +
  ylab("PAV incidence") +
  theme_def +
  theme(legend.position = "none",
        axis.title.x = element_blank())

# RPV infection prevalence
rpv_fig <- ggplot(rpv_dat, aes(microbes_f, rpv, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2), shape = 21) +
  geom_text(data = rpv_lab, aes(label = label), nudge_x = -0.45, fontface = "bold") +
  facet_wrap(~ inoculation) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  guides(fill = guide_legend(direction = "horizontal", title.position = "top")) +
  xlab("Microbe inoculation") +
  ylab("RPV incidence") +
  theme_def +
  theme(legend.position = c(0.13, 0.13),
        strip.text = element_blank(),
        legend.margin = margin(0, 0, 0, 0))

# combine figures
pdf("output/infection_figure_microbes.pdf", width = 5, height = 5)
plot_grid(pav_fig, rpv_fig,
          nrow = 2)
dev.off()


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
pp_check(pav_mod2, nsamples = 50)

# microbes model
pav_mic_mod1 <- update(pav_mod2,
                       newdata = pav_dat,
                       formula = pav ~ microbes * N_added * inoc_rpv,
                       prior = c(prior(normal(0, 10), class = Intercept),
                                 prior(normal(0, 10), class = b)))

# check model
summary(pav_mic_mod1)
plot(pav_mic_mod1)
pp_check(pav_mic_mod1, nsamples = 50)

# compare with loo
pav_loo2 <- loo(pav_mod2, reloo = T)
pav_loo2 
# all k < 0.7 and elpd_loo <= 0.1
# good model fit
# reasonable to do model comparison
pav_mic_loo1 <- loo(pav_mic_mod1)
pav_mic_loo1 
# all k < 0.7 and elpd_loo <= 0.1
# good model fit
# reasonable to do model comparison
loo_compare(pav_loo2, pav_mic_loo1)
# microbes model is slightly preferred


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
pp_check(rpv_mod2, nsamples = 50)

# microbes model
rpv_mic_mod1 <- update(rpv_mod2,
                       newdata = rpv_dat,
                       formula = rpv ~ microbes * N_added * inoc_pav,
                       prior = c(prior(normal(0, 10), class = Intercept),
                                 prior(normal(0, 10), class = b)))

# check model
summary(rpv_mic_mod1)
plot(rpv_mic_mod1)
pp_check(rpv_mic_mod1, nsamples = 50)

# compare with loo
rpv_loo2 <- loo(rpv_mod2)
rpv_loo2 
# all k < 0.7 and elpd_loo <= 0.1
# good model fit
# reasonable to do model comparison
rpv_mic_loo1 <- loo(rpv_mic_mod1, reloo = T)
rpv_mic_loo1 
# all k < 0.7
# elpd_loo > 0.1, but still small compared to other SE
# good model fit
# reasonable to do model comparison
loo_compare(rpv_loo2, rpv_mic_loo1)
# microbes model is preferred


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
pp_check(co_mod2, nsamples = 50)

# microbes model
co_mic_mod1 <- update(co_mod2,
                      newdata = co_dat,
                      formula = coinfection ~ microbes * N_added,
                      prior = c(prior(normal(0, 10), class = Intercept),
                                prior(normal(0, 10), class = b)))

# check model
summary(co_mic_mod1)
plot(co_mic_mod1)
pp_check(co_mic_mod1, nsamples = 50)

# compare with loo
co_loo2 <- loo(co_mod2, reloo = T)
co_loo2 
# all k < 0.7 and elpd_loo <= 0.1
# elpd_loo > 0.1, but still small compared to other SE
# good model fit
# reasonable to do model comparison
co_mic_loo1 <- loo(co_mic_mod1, reloo = T)
co_mic_loo1 
# all k < 0.7
# elpd_loo > 0.1, but still small compared to other SE
# good model fit
# reasonable to do model comparison
loo_compare(co_loo2, co_mic_loo1)
# microbes model is preferred


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
         N_effect = high_N - icp,
         coinoc = exp(b_Intercept + b_inoc_rpv)/(1 + exp(b_Intercept + b_inoc_rpv)),
         co_effect = coinoc - icp,
         coinoc_N = exp(b_Intercept + b_inoc_rpv + b_N_added + b_N_rpv_int)/(1 + exp(b_Intercept + b_inoc_rpv + b_N_added + b_N_rpv_int)),
         N_co_effect = coinoc_N - coinoc,
         soilA = exp(b_Intercept + b_soilambientN)/(1 + exp(b_Intercept + b_soilambientN)),
         soilL = exp(b_Intercept + b_soillowN)/(1 + exp(b_Intercept + b_soillowN)),
         soilH = exp(b_Intercept + b_soilhighN)/(1 + exp(b_Intercept + b_soilhighN)),
         soilL_effect = soilL - icp,
         soilA_N = exp(b_Intercept + b_soilambientN + b_N_added + b_N_soilA_int)/(1 + exp(b_Intercept + b_soilambientN + b_N_added + b_N_soilA_int)),
         soilL_N = exp(b_Intercept + b_soillowN + b_N_added + b_N_soilL_int)/(1 + exp(b_Intercept + b_soillowN + b_N_added + b_N_soilL_int)),
         soilH_N = exp(b_Intercept + b_soilhighN + b_N_added + b_N_soilH_int)/(1 + exp(b_Intercept + b_soilhighN + b_N_added + b_N_soilH_int)),
         N_soilA_effect = soilA_N - soilA,
         N_soilL_effect = soilL_N - soilL,
         N_soilH_effect = soilH_N - soilH,
         soilL_N_effect = soilL_N - high_N,
         soilA_co = exp(b_Intercept + b_soilambientN + b_inoc_rpv + b_soilA_rpv_int)/(1 + exp(b_Intercept + b_soilambientN + b_inoc_rpv + b_soilA_rpv_int)),
         soilL_co = exp(b_Intercept + b_soillowN + b_inoc_rpv + b_soilL_rpv_int)/(1 + exp(b_Intercept + b_soillowN + b_inoc_rpv + b_soilL_rpv_int)),
         soilH_co = exp(b_Intercept + b_soilhighN + b_inoc_rpv + b_soilH_rpv_int)/(1 + exp(b_Intercept + b_soilhighN + b_inoc_rpv + b_soilH_rpv_int)),
         co_soilA_effect = soilA_co - soilA,
         co_soilL_effect = soilL_co - soilL,
         co_soilH_effect = soilH_co - soilH,
         soilA_co_N = exp(b_Intercept + b_soilambientN + b_inoc_rpv + b_N_added + b_soilA_rpv_int + b_N_rpv_int + b_N_soilA_int + b_N_soilA_rpv_int)/(1 + exp(b_Intercept + b_soilambientN + b_inoc_rpv + b_N_added + b_soilA_rpv_int + b_N_rpv_int + b_N_soilA_int + b_N_soilA_rpv_int)),
         soilL_co_N = exp(b_Intercept + b_soillowN + b_inoc_rpv + b_N_added + b_soilL_rpv_int + b_N_rpv_int + b_N_soilL_int + b_N_soilL_rpv_int)/(1 + exp(b_Intercept + b_soillowN + b_inoc_rpv + b_N_added + b_soilL_rpv_int + b_N_rpv_int + b_N_soilL_int + b_N_soilL_rpv_int)),
         N_co_soilA_effect = soilA_co_N - soilA_co,
         N_co_soilL_effect = soilL_co_N - soilL_co)

mean_hdi(pav_post2$icp)
mean_hdi(pav_post2$N_effect)
mean_hdi(pav_post2$co_effect)
mean_hdi(pav_post2$N_co_effect)
mean_hdi(pav_post2$soilL_effect)
mean_hdi(pav_post2$soilA)
mean_hdi(pav_post2$soilA_N)
mean_hdi(pav_post2$N_soilA_effect)
mean_hdi(pav_post2$soilL)
mean_hdi(pav_post2$soilL_N)
mean_hdi(pav_post2$N_soilL_effect)
mean_hdci(pav_post2$soilH)
mean_hdci(pav_post2$soilH_N)
mean_hdci(pav_post2$N_soilH_effect)
mean_hdi(pav_post2$soilL_N_effect)
mean_hdi(pav_post2$co_soilA_effect)
mean_hdi(pav_post2$co_soilL_effect)
mean_hdi(pav_post2$co_soilH_effect)
mean_hdi(pav_post2$N_co_soilA_effect)
mean_hdi(pav_post2$N_co_soilL_effect)

# RPV incidence
rpv_post2 <- posterior_samples(rpv_mod2) %>%
  mutate(icp = exp(b_Intercept)/(1 + exp(b_Intercept)),
         coinoc = exp(b_Intercept + b_inoc_pav)/(1 + exp(b_Intercept + b_inoc_pav)),
         co_effect = coinoc - icp)

mean_hdi(rpv_post2$icp)
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


#### output ####
save(pav_mod2, file = "output/pav_bayesian_model_soil.rda")
save(pav_mic_mod1, file = "output/pav_bayesian_model_microbes.rda")
save(rpv_mod2, file = "output/rpv_bayesian_model_soil.rda")
save(rpv_mic_mod1, file = "output/rpv_bayesian_model_microbes.rda")
save(co_mod2, file = "output/co_bayesian_model_soil.rda")
save(co_mic_mod1, file = "output/co_bayesian_model_microbes.rda")

write_csv(tidy(summary(pav_mod2)$fixed), "output/pav_bayesian_model_soil.csv")
write_csv(tidy(summary(pav_mic_mod1)$fixed), "output/pav_bayesian_model_microbes.csv")
write_csv(tidy(summary(rpv_mod2)$fixed), "output/rpv_bayesian_model_soil.csv")
write_csv(tidy(summary(rpv_mic_mod1)$fixed), "output/rpv_bayesian_model_microbes.csv")
write_csv(tidy(summary(co_mod2)$fixed), "output/co_bayesian_model_soil.csv")
write_csv(tidy(summary(co_mic_mod1)$fixed), "output/co_bayesian_model_microbes.csv")