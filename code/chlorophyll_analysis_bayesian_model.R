##### info ####

# goal: analyze chlorophyll data any visible PCR band is counted as indicators of infection (rounding up)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(broom)
library(brms)

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

dat2 %>%
  group_by(pav, rpv, N_added, soil) %>%
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
pan_lab <- tibble(infection =levels(dat2$infection) %>%
                    fct_relevel("Mock inoculation", "PAV infection", "RPV infection"),
                  label = c("(a)", "(b)", "(c)", "(d)")) %>%
  mutate(microbes_f = "sterile",
         chlorophyll = 32,
         nitrogen_added = "low")

# chlorophyll figure
pdf("output/chlorophyll_figure_microbes.pdf", width = 5, height = 5)
ggplot(dat2, aes(microbes_f, chlorophyll, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2), shape = 21) +
  geom_text(data = pan_lab, aes(label = label), nudge_x = -0.45, fontface = "bold") +
  facet_wrap(~ infection) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  guides(fill = guide_legend(override.aes = list(shape = 21), direction = "horizontal", title.position = "top")) +
  xlab("Microbe inoculation") +
  ylab("Leaf chlorophyll content [SPAD]") +
  theme_def +
  theme(legend.position = c(0.85, 0.92))
dev.off()


#### chlorophyll model ####

# initial fit
chlor_mod1 <- brm(log_chlorophyll ~ soil * N_added * infection, 
                  data = dat2,
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
                         newdata = dat2,
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
# all k < 0.7 and elpd_loo is 0.3, but smaller than other SE
# good model fit
# reasonable to do model comparison
chlor_mic_loo1 <- loo(chlor_mic_mod1, reloo = T)
chlor_mic_loo1 
# all k < 0.7 and elpd_loo is 0.2, but less than other SE
# good model fit
# reasonable to do model comparison
loo_compare(chlor_loo2, chlor_mic_loo1)
# microbes model is preferred


#### output ####
save(chlor_mod2, file = "output/chlor_bayesian_model_soil.rda")
save(chlor_mic_mod1, file = "output/chlor_bayesian_model_microbes.rda")

write_csv(tidy(summary(chlor_mod2)$fixed), "output/chlor_bayesian_model_soil.csv")
write_csv(tidy(summary(chlor_mic_mod1)$fixed), "output/chlor_bayesian_model_microbes.csv")
