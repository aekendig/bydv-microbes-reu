##### info ####

# authors: Amy Kendig and Casey Easterday
# date last edited: 10/14/20
# goal: analyze data when only PCR bands at least as intense as controls are counted as indicators of infection (rounding down)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(brms)
library(cowplot)
library(tidybayes)

# import data
dat <- read_csv("intermediate-data/bydv_microbes_data_rounded_down.csv")


#### edit data ####

# reorganize factor levels
dat2 <- dat %>%
  mutate(soil = fct_relevel(soil, "sterile", "low", "medium", "high"),
         inoculation = case_when(disease %in% c("PAV", "RPV") ~ "Single inoculation",
                                 disease == "Co" ~ "Co-inoculation",
                                 disease == "Healthy" ~ "Mock inoculation") %>%
           fct_relevel("Mock inoculation", "Single inoculation"),
         nitrogen_added = fct_relevel(nitrogen_added, "low", "high"))

# separate by virus
pav_dat <- dat2 %>%
  filter(inoc_pav == 1)

rpv_dat <- dat2 %>%
  filter(inoc_rpv == 1)

# microbe comparison
pav_microbe_dat = filter(pav_dat, soil_N == 0)
rpv_microbe_dat = filter(rpv_dat, soil_N == 0)


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
shape_pal = c(21, 22)

# panel labels
pav_lab <- tibble(inoculation = c("Single inoculation", "Co-inoculation") %>% fct_relevel("Single inoculation"),
                  label = c("(a)", "(b)")) %>%
  mutate(soil = "sterile",
         pav = 1,
         nitrogen_added = "low")

rpv_lab <- pav_lab %>%
  mutate(label = recode(label, "(a)" = "(c)", "(b)" = "(d)"),
         rpv = 1)

# PAV infection prevalence
pav_fig  <- ggplot(pav_dat, aes(soil, pav, fill = nitrogen_added, shape = inoculation)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2)) +
  geom_text(data = pav_lab, aes(label = label), nudge_x = -0.4) +
  facet_wrap(~ inoculation) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  scale_shape_manual(values = shape_pal, guide = F) +
  guides(fill = guide_legend(override.aes = list(shape = 21), direction = "horizontal", title.position = "top")) +
  xlab("Field soil N treatment") +
  ylab("PAV infection prevalence") +
  theme_def +
  theme(legend.position = c(0.13, 0.13),
        axis.title.x = element_blank())
# N reduces PAV unless microbes are added
# coinfection reduces PAV unless microbes are added

# RPV infection prevalence
rpv_fig <- ggplot(rpv_dat, aes(soil, rpv, fill = nitrogen_added, shape = inoculation)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2)) +
  geom_text(data = rpv_lab, aes(label = label), nudge_x = -0.4) +
  facet_wrap(~ inoculation) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  scale_shape_manual(values = shape_pal, guide = F) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab("Field soil N treatment") +
  ylab("RPV infection prevalence") +
  theme_def +
  theme(legend.position = "none",
        strip.text = element_blank())
# coinfection reduces RPV unless microbes are added


# combine figures
pdf("output/infection_figure_rounded_down.pdf", width = 6, height = 6)
plot_grid(pav_fig, rpv_fig,
          nrow = 2)
dev.off()


#### PAV microbe model ####

# intercept
pav_microbe_dat %>%
  filter(microbes == 0 & N_added == 0 & inoc_rpv == 0) %>%
  summarise(mean = mean(pav),
            zeros = sum(pav == 0))

# initial fit
pav_microbe_mod1 <- brm(pav ~ microbes * N_added * inoc_rpv,
                data = pav_microbe_dat, 
                family = bernoulli,
                prior = c(prior(normal(0, 100), class = Intercept),
                          prior(normal(0, 50), class = b)),
                iter = 6000, warmup = 1000, chains = 1)
summary(pav_microbe_mod1)

# add chains
pav_microbe_mod2 <- update(pav_microbe_mod1, chains = 3)

# check model
summary(pav_microbe_mod2)
plot(pav_microbe_mod2)
pp_check(pav_microbe_mod2, nsamples = 50)

# treatments
pav_microbe_trt <- pav_microbe_dat %>%
  select(microbes, N_added, inoc_rpv) %>%
  unique()

# posterior samples
pav_microbe_post <- as_tibble(posterior_samples(pav_microbe_mod2))

# prevalence
pav_microbe_prev <- pav_microbe_trt %>%
  expand_grid(pav_microbe_post) %>%
  rename_with(~ gsub(":", "_", .x, fixed = T)) %>%
  mutate(lin_exp = b_Intercept + 
           b_microbes * microbes + 
           b_N_added * N_added + 
           b_inoc_rpv * inoc_rpv +
           b_microbes_N_added * microbes * N_added +
           b_microbes_inoc_rpv * microbes * inoc_rpv +
           b_N_added_inoc_rpv * N_added * inoc_rpv +
           b_microbes_N_added_inoc_rpv * microbes * N_added * inoc_rpv,
         prev = exp(lin_exp) / (1 + exp(lin_exp)) * 100) %>%
  select(-c(b_Intercept:lin_exp))

# summary
pav_microbe_prev %>%
  group_by(microbes, N_added, inoc_rpv) %>%
  mean_hdci() %>%
  mutate(prev = round(prev),
         .lower = round(.lower),
         .upper = round(.upper))
# hdi splits the interval of the 0/0/0 treatment into three, but the discontinuities are on the scale of 0.001

# treatment effects
pav_microbe_prev %>%
  mutate(microbes = recode(microbes, "0" = "sterile", "1" = "microbes"),
         N_added = recode(N_added, "0" = "lowN", "1" = "hiN"),
         inoc_rpv = recode(inoc_rpv, "0" = "single", "1" = "co"),
         treatment = paste(microbes, N_added, inoc_rpv, sep = "_"),
         replicate = rep(seq(1:15000), 8)) %>%
  select(-c(microbes:inoc_rpv)) %>%
  pivot_wider(names_from = treatment, values_from = prev) %>%
  mutate(N_eff = sterile_hiN_single - sterile_lowN_single,
         mic_lowN_eff = microbes_lowN_single - sterile_lowN_single,
         mic_hiN_eff = microbes_hiN_single - sterile_hiN_single,
         co_lowN_eff = sterile_lowN_co - sterile_lowN_single,
         co_hiN_eff = sterile_hiN_co - sterile_hiN_single) %>%
  select(-c(replicate:sterile_hiN_single)) %>%
  pivot_longer(cols = N_eff:co_hiN_eff, names_to = "treatment", values_to = "effect") %>%
  group_by(treatment) %>%
  mean_hdi %>%
  mutate(effect = round(effect),
         .lower = round(.lower),
         .upper = round(.upper))
  

#### RPV microbe model ####

# initial fit
rpv_microbe_mod1 <- brm(rpv ~ microbes * N_added * inoc_pav,
                        data = rpv_microbe_dat, 
                        family = bernoulli,
                        prior = c(prior(normal(0, 100), class = Intercept),
                                  prior(normal(0, 50), class = b)),
                        iter = 6000, warmup = 1000, chains = 1)
summary(rpv_microbe_mod1)

# add chains
rpv_microbe_mod2 <- update(rpv_microbe_mod1, chains = 3)

# check model
summary(rpv_microbe_mod2)
plot(rpv_microbe_mod2)
pp_check(rpv_microbe_mod2, nsamples = 50)

# treatments
rpv_microbe_trt <- rpv_microbe_dat %>%
  select(microbes, N_added, inoc_pav) %>%
  unique()

# posterior samples
rpv_microbe_post <- as_tibble(posterior_samples(rpv_microbe_mod2))

# prevalence
rpv_microbe_prev <- rpv_microbe_trt %>%
  expand_grid(rpv_microbe_post) %>%
  rename_with(~ gsub(":", "_", .x, fixed = T)) %>%
  mutate(lin_exp = b_Intercept + 
           b_microbes * microbes + 
           b_N_added * N_added + 
           b_inoc_pav * inoc_pav +
           b_microbes_N_added * microbes * N_added +
           b_microbes_inoc_pav * microbes * inoc_pav +
           b_N_added_inoc_pav * N_added * inoc_pav +
           b_microbes_N_added_inoc_pav * microbes * N_added * inoc_pav,
         prev = exp(lin_exp) / (1 + exp(lin_exp)) * 100) %>%
  select(-c(b_Intercept:lin_exp))

# summary
rpv_microbe_prev %>%
  group_by(microbes, N_added, inoc_pav) %>%
  mean_hdci() %>%
  mutate(prev = round(prev),
         .lower = round(.lower),
         .upper = round(.upper))
# hdi splits the interval of the 0/0/0 treatment into three, but the discontinuities are on the scale of 0.001

# treatment effects
rpv_microbe_prev %>%
  mutate(microbes = recode(microbes, "0" = "sterile", "1" = "microbes"),
         N_added = recode(N_added, "0" = "lowN", "1" = "hiN"),
         inoc_pav = recode(inoc_pav, "0" = "single", "1" = "co"),
         treatment = paste(microbes, N_added, inoc_pav, sep = "_"),
         replicate = rep(seq(1:15000), 8)) %>%
  select(-c(microbes:inoc_pav)) %>%
  pivot_wider(names_from = treatment, values_from = prev) %>%
  mutate(mic_lowN_eff = microbes_lowN_single - sterile_lowN_single,
         mic_hiN_eff = microbes_hiN_single - sterile_hiN_single) %>%
  select(-c(replicate:sterile_hiN_single)) %>%
  pivot_longer(cols = mic_lowN_eff:mic_hiN_eff, names_to = "treatment", values_to = "effect") %>%
  group_by(treatment) %>%
  mean_hdi %>%
  mutate(effect = round(effect),
         .lower = round(.lower),
         .upper = round(.upper))


#### output ####

save(pav_microbe_mod2, file = "output/pav_microbe_model_rounded_down.rda")
save(rpv_microbe_mod2, file = "output/rpv_microbe_model_rounded_down.rda")
