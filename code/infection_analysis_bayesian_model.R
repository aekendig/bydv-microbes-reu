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