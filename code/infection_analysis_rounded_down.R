##### info ####

# goal: analyze data when only PCR bands at least as intense as controls are counted as indicators of infection (rounding down)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(cowplot)
library(broom)

# import data
dat <- read_csv("intermediate-data/bydv_microbes_data_rounded_down.csv")


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
                                 TRUE ~ 0))

# separate by virus
pav_dat <- dat2 %>%
  filter(inoc_pav == 1)

rpv_dat <- dat2 %>%
  filter(inoc_rpv == 1)

co_dat <- dat2 %>%
  filter(disease == "Co")


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
  mutate(soil = "sterile",
         pav = 1,
         nitrogen_added = "low")

rpv_lab <- pav_lab %>%
  mutate(label = recode(label, "(a)" = "(c)", "(b)" = "(d)"),
         rpv = 1)

# PAV infection prevalence
pav_fig  <- ggplot(pav_dat, aes(soil, pav, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2), shape = 21) +
  geom_text(data = pav_lab, aes(label = label), nudge_x = -0.4) +
  facet_wrap(~ inoculation) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  guides(fill = guide_legend(direction = "horizontal", title.position = "top")) +
  xlab("Field soil N treatment") +
  ylab("PAV incidence") +
  theme_def +
  theme(legend.position = c(0.13, 0.13),
        axis.title.x = element_blank())

# RPV infection prevalence
rpv_fig <- ggplot(rpv_dat, aes(soil, rpv, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2), shape = 21) +
  geom_text(data = rpv_lab, aes(label = label), nudge_x = -0.4) +
  facet_wrap(~ inoculation) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  xlab("Soil microbes") +
  ylab("RPV incidence") +
  theme_def +
  theme(legend.position = "none",
        strip.text = element_blank())

# combine figures
pdf("output/infection_figure_rounded_down.pdf", width = 6, height = 6)
plot_grid(pav_fig, rpv_fig,
          nrow = 2)
dev.off()


#### PAV model ####

# full model
pav_mod1 <- glm(pav ~ soil * N_added * inoc_rpv, 
                data = pav_dat,
                family = "binomial")
summary(pav_mod1)

# remove 3-way interaction?
pav_mod2 <- update(pav_mod1, ~. -soil:N_added:inoc_rpv)
summary(pav_mod2)
anova(pav_mod1, pav_mod2, test = "Chi") # yes

# remove 2-way interactions?
pav_mod3 <- update(pav_mod2, ~. -N_added:inoc_rpv)
summary(pav_mod3)
anova(pav_mod2, pav_mod3, test = "Chi") # yes

pav_mod4 <- update(pav_mod3, ~. -soil:N_added)
summary(pav_mod4)
anova(pav_mod3, pav_mod4, test = "Chi") # yes

pav_mod5 <- update(pav_mod4, ~. -soil:inoc_rpv)
summary(pav_mod5)
anova(pav_mod4, pav_mod5, test = "Chi") # no

# remove main effects?
pav_mod6 <- update(pav_mod4, ~. -N_added)
summary(pav_mod6)
anova(pav_mod4, pav_mod6, test = "Chi") # yes


#### PAV values ####

pav_dat %>%
  filter(soil == "sterile" & inoc_rpv == 0) %>%
  summarise(mean_pav = mean(pav),
            se_pav = sd(pav)/sqrt(n()))

pav_dat %>%
  group_by(soil, disease) %>%
  summarise(mean_pav = mean(pav)) %>%
  pivot_wider(names_from = disease,
              values_from = mean_pav) %>%
  mutate(change = PAV - Co)
  


#### RPV model ####

# full model
rpv_mod1 <- glm(rpv ~ soil * N_added * inoc_pav, 
                data = rpv_dat,
                family = "binomial")
summary(rpv_mod1)

# remove 3-way interaction?
rpv_mod2 <- update(rpv_mod1, ~. -soil:N_added:inoc_pav)
summary(rpv_mod2)
anova(rpv_mod1, rpv_mod2, test = "Chi") # yes

# remove 2-way interactions?
rpv_mod3 <- update(rpv_mod2, ~. -soil:N_added)
summary(rpv_mod3)
anova(rpv_mod2, rpv_mod3, test = "Chi") # yes

rpv_mod4 <- update(rpv_mod3, ~. -soil:inoc_pav)
summary(rpv_mod4)
anova(rpv_mod3, rpv_mod4, test = "Chi") # yes

rpv_mod5 <- update(rpv_mod4, ~. -N_added:inoc_pav)
summary(rpv_mod5)
anova(rpv_mod4, rpv_mod5, test = "Chi") # no

# remove main effects?
rpv_mod6 <- update(rpv_mod4, ~. -soil)
summary(rpv_mod6)
anova(rpv_mod4, rpv_mod6, test = "Chi") # yes


#### RPV values ####

rpv_dat %>%
  filter(N_added == 0 & inoc_pav == 0) %>%
  summarise(mean_rpv = mean(rpv),
            se_rpv = sd(rpv)/sqrt(n()))

rpv_dat %>%
  group_by(N_added, disease) %>%
  summarise(mean_rpv = mean(rpv)) %>%
  pivot_wider(names_from = disease,
              values_from = mean_rpv) %>%
  mutate(change = RPV - Co)


#### coinfection model ####

# full model
co_mod1 <- glm(coinfection ~ soil * N_added, 
                data = co_dat,
                family = "binomial")
summary(co_mod1)

# remove 2-way interaction?
co_mod2 <- update(co_mod1, ~. -soil:N_added)
summary(co_mod2)
anova(co_mod1, co_mod2, test = "Chi") # yes

# remove main effects?
co_mod3 <- update(co_mod2, ~. -soil)
summary(co_mod3)
anova(co_mod2, co_mod3, test = "Chi") # no

co_mod4 <- update(co_mod2, ~. -N_added)
summary(co_mod4)
anova(co_mod2, co_mod4, test = "Chi") # no

# high standard errors due to no cases at intercept
# alternative test

# switch factor levels
co_dat2 <- co_dat %>%
  mutate(soil = fct_rev(soil))

# full model
co_mod5 <- glm(coinfection ~ soil * N_limit, 
               data = co_dat2,
               family = "binomial")
summary(co_mod5) # still has some high SE values

# remove 2-way interaction?
co_mod6 <- update(co_mod5, ~. -soil:N_limit)
summary(co_mod6)
anova(co_mod5, co_mod6, test = "Chi") # yes

# remove main effects?
co_mod7 <- update(co_mod6, ~. -soil)
summary(co_mod7)
anova(co_mod6, co_mod7, test = "Chi") # no

co_mod8 <- update(co_mod6, ~. -N_limit)
summary(co_mod8)
anova(co_mod6, co_mod8, test = "Chi") # no


#### coinfection values ####

co_dat %>%
  filter(N_added == 0 & soil == "sterile") %>%
  summarise(mean_co = mean(co),
            se_co = sd(co)/sqrt(n()))

co_dat %>%
  group_by(soil) %>%
  summarise(mean_co = mean(co)) %>%
  filter(soil != "sterile") %>%
  summarise(mean_change = mean(mean_co))

co_dat %>%
  group_by(N_added) %>%
  summarise(mean_co = mean(co)) %>%
  pivot_wider(names_from = N_added,
              values_from = mean_co) %>%
  rename("low" = "0",
         "high" = "1") %>%
  mutate(change = high - low)



#### output ####
save(pav_mod6, file = "output/pav_model_rounded_down.rda")
save(rpv_mod6, file = "output/rpv_model_rounded_down.rda")
save(co_mod2, file = "output/coinfection_model_rounded_down.rda")

write_csv(tidy(pav_mod6), "output/pav_model_rounded_down.csv")
write_csv(tidy(rpv_mod6), "output/rpv_model_rounded_down.csv")
write_csv(tidy(co_mod2), "output/coinfection_model_rounded_down.csv")
