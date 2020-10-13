##### info ####

# authors: Amy Kendig and Casey Easterday
# date last edited: 10/13/20
# goal: analyze data when including all visible PCR bands as indicators of infection (rounding up)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(MuMIn)
library(cowplot)

# import data
dat <- read_csv("intermediate-data/bydv_microbes_data_rounded_up.csv")


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
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  xlab("Field soil N treatment") +
  ylab("PAV infection prevalence") +
  theme_def +
  theme(legend.position = c(0.4, 0.2),
        axis.title.x = element_blank())
# adding microbes reduces PAV infection unless N is high (because N reduces infection) or microbes have long-term exposure to high N
# adding microbes or N helps ameliorate negative effects of co-inoculation on PAV infection

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
# adding microbes reduces RPV infection, especially when N is high (because N increases infection)
# adding microbes or N helps ameliorate negative effects of co-inoculation on RPV infection

# combine figures
pdf("output/infection_figure_rounded_up.pdf", width = 6, height = 6)
plot_grid(pav_fig, rpv_fig,
          nrow = 2)
dev.off()


#### PAV model ####

# initial fit
pav_mod1 <- glm(pav ~ soil * N_added * inoc_rpv, data = pav_dat, 
                family = binomial,
                na.action = na.fail)
summary(pav_mod1)
# no significant effects
# many of the standard errors are the same because the intercept is only 1's

# switch intercept
# pav_mod2 <- glm(pav ~ soil * N_limit * inoc_rpv, data = pav_dat, 
#                 family = binomial,
#                 na.action = na.fail)
# summary(pav_mod2)

# model average using AIC
pav_mod_avg <- model.avg(dredge(pav_mod1), cumsum(weight) <= 0.95)
summary(pav_mod_avg)
# component model results are the same whether model 1 or 2 is used
# summed weight is 1 -- because you have to include two 0.02 weighted models to get to 0.95?
# copied and pasted output to Excel


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