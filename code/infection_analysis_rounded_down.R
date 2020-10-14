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
pav_dat_microbes = filter(pav_dat, soil_N == 0)
rpv_dat_microbes = filter(rpv_dat, soil_N == 0)


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


#### PAV model ####

# effect of microbes
pav_mod1 <- brm(pav ~ microbes * N_added * inoc_rpv,
                data = pav_dat_microbes, 
                family = bernoulli,
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 10), class = b)),
                iter = 6000, warmup = 1000, chains = 1)
summary(pav_mod1)

# add chains
pav_mod2 <- update(pav_mod1, chains = 3)
summary(pav_mod2)
plot(pav_mod2)
pp_check(pav_mod2, nsamples = 50)

#### left off here ####



#### RPV model ####

# initial fit
rpv_mod1 <- glm(rpv ~ soil * N_added * inoc_pav, data = rpv_dat, 
                family = binomial,
                na.action = na.fail)
summary(rpv_mod1)

# model average using AIC
rpv_mod_avg <- model.avg(dredge(rpv_mod1), subset = cumsum(weight) <= .95)
summary(rpv_mod_avg)
# copied and pasted output to Excel


#### output ####

save(pav_mod2, file = "output/pav_microbes_model_rounded_down.rda")
