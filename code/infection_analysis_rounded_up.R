##### info ####

# goal: analyze data when any visible PCR band is counted as indicators of infection (rounding up)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(cowplot)
library(broom)

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
         soil_N_added = ifelse(soil %in% c("sterile", "ambient N"), 0, 1),
         soil_high_N = ifelse(soil == "high N", 1, 0))

# reset soil contrasts (Helmert coding)
# https://stats.idre.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/#User
# https://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html
#                          S   A   L   H
contrast_mat <- matrix(c(1/4, 1/4, 1/4, 1/4, # intercept
                         3/4, -1/4, -1/4, -1/4, # sterile vs others
                         0,  2/3, -1/3, -1/3, # ambient vs. N fert   
                         0,  0,  1/2,	-1/2), # low vs. high
                       ncol = 4)

# contrast_tmat <- solve(t(contrast_mat)) # inverse of transposed matrix (for user-defined)

# cor(contrast_tmat[, 2:4]) # test orthogonality
cor(contrast_mat[, 2:4])

contrasts(dat2$soil) <- contrast_mat[, 2:4]
# note that contr.helmert is reverse Helmert coding according to the ucla reference

# separate by virus
pav_dat <- dat2 %>%
  filter(inoc_pav == 1)

rpv_dat <- dat2 %>%
  filter(inoc_rpv == 1)

co_dat <- dat2 %>%
  filter(disease == "Co")

# sample sizes
pav_dat %>%
  group_by(soil, N_added, disease) %>%
  count()

rpv_dat %>%
  group_by(soil, N_added, disease) %>%
  count()


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
  geom_text(data = pav_lab, aes(label = label), nudge_x = -0.4, fontface = "bold") +
  facet_wrap(~ inoculation) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  guides(fill = guide_legend(direction = "horizontal", title.position = "top")) +
  xlab("Field soil N treatment") +
  ylab("PAV incidence") +
  theme_def +
  theme(legend.position = c(0.87, 0.13),
        axis.title.x = element_blank(),
        legend.margin = margin(0, 0, 0, 0))

# RPV infection prevalence
rpv_fig <- ggplot(rpv_dat, aes(soil, rpv, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2), shape = 21) +
  geom_text(data = rpv_lab, aes(label = label), nudge_x = -0.4, fontface = "bold") +
  facet_wrap(~ inoculation) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  xlab("Soil microbes") +
  ylab("RPV incidence") +
  theme_def +
  theme(legend.position = "none",
        strip.text = element_blank())

# combine figures
pdf("output/infection_figure_rounded_up.pdf", width = 6, height = 6)
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
anova(pav_mod1, pav_mod2, test = "Chi") # yes
summary(pav_mod2)

# remove 2-way interactions?
pav_mod3 <- update(pav_mod2, ~. -soil:inoc_rpv)
anova(pav_mod2, pav_mod3, test = "Chi") # no

pav_mod4 <- update(pav_mod2, ~. -soil:N_added)
anova(pav_mod2, pav_mod4, test = "Chi") # yes
summary(pav_mod4)

pav_mod5 <- update(pav_mod4, ~. -soil:inoc_rpv)
anova(pav_mod4, pav_mod5, test = "Chi") # no

pav_mod6 <- update(pav_mod4, ~. -N_added:inoc_rpv)
anova(pav_mod4, pav_mod6, test = "Chi") # yes
summary(pav_mod6)

pav_mod7 <- update(pav_mod6, ~. -soil:inoc_rpv)
anova(pav_mod6, pav_mod7, test = "Chi") # no

# remove main effects?
pav_mod8 <- update(pav_mod6, ~. -N_added)
anova(pav_mod6, pav_mod8, test = "Chi") # yes
summary(pav_mod8)

# check interaction again
pav_mod9 <- update(pav_mod8, ~. -soil:inoc_rpv)
anova(pav_mod8, pav_mod9, test = "Chi") # no

# final model
summary(pav_mod8)


#### PAV values ####

pav_dat %>%
  group_by(soil, disease) %>%
  summarise(mean_pav = mean(pav)) %>%
  pivot_wider(names_from = disease,
              values_from = mean_pav) %>%
  mutate(change = PAV - Co)

pav_dat %>%
  group_by(soil, disease) %>%
  summarise(mean_pav = mean(pav)) %>%
  ungroup() %>%
  mutate(logit = log(mean_pav/(1-mean_pav)))


#### PAV microbe model ####

# full model
pav_mic_mod1 <- glm(pav ~ microbes * N_added * inoc_rpv, 
                    data = pav_dat,
                    family = "binomial")

summary(pav_mic_mod1)

# remove 3-way interaction?
pav_mic_mod2 <- update(pav_mic_mod1, ~. -microbes:N_added:inoc_rpv)
anova(pav_mic_mod1, pav_mic_mod2, test = "Chi") # yes
summary(pav_mic_mod2)

# remove 2-way interactions?
pav_mic_mod3 <- update(pav_mic_mod2, ~. -N_added:inoc_rpv)
anova(pav_mic_mod2, pav_mic_mod3, test = "Chi") # yes
summary(pav_mic_mod3)

pav_mic_mod4 <- update(pav_mic_mod3, ~. -microbes:N_added)
anova(pav_mic_mod3, pav_mic_mod4, test = "Chi") # yes
summary(pav_mic_mod4)

pav_mic_mod5 <- update(pav_mic_mod4, ~. -microbes:inoc_rpv)
anova(pav_mic_mod4, pav_mic_mod5, test = "Chi") # no

# remove main effects?
pav_mic_mod6 <- update(pav_mic_mod4, ~. -N_added)
anova(pav_mic_mod4, pav_mic_mod6, test = "Chi") # yes
summary(pav_mic_mod6)

# re-test two-way interaction
pav_mic_mod7 <- update(pav_mic_mod6, ~. -microbes:inoc_rpv)
anova(pav_mic_mod6, pav_mic_mod7, test = "Chi") # no

# final model
summary(pav_mic_mod6)


#### RPV model ####

# full model
rpv_mod1 <- glm(rpv ~ soil * N_added * inoc_pav, 
                data = rpv_dat,
                family = "binomial")
summary(rpv_mod1)

# remove 3-way interaction?
rpv_mod2 <- update(rpv_mod1, ~. -soil:N_added:inoc_pav)
anova(rpv_mod1, rpv_mod2, test = "Chi") # yes
summary(rpv_mod2)

# remove 2-way interactions?
rpv_mod3 <- update(rpv_mod2, ~. -N_added:inoc_pav)
anova(rpv_mod2, rpv_mod3, test = "Chi") # yes
summary(rpv_mod3)

rpv_mod4 <- update(rpv_mod3, ~. -soil:N_added)
anova(rpv_mod3, rpv_mod4, test = "Chi") # yes
summary(rpv_mod4)

rpv_mod5 <- update(rpv_mod4, ~. -soil:inoc_pav)
anova(rpv_mod4, rpv_mod5, test = "Chi") # yes
summary(rpv_mod5)

# remove main effects?
rpv_mod6 <- update(rpv_mod5, ~. -soil)
anova(rpv_mod5, rpv_mod6, test = "Chi") # yes
summary(rpv_mod6)

rpv_mod7 <- update(rpv_mod6, ~. -N_added)
anova(rpv_mod6, rpv_mod7, test = "Chi") # yes
summary(rpv_mod7)

rpv_mod8 <- update(rpv_mod7, ~. -inoc_pav)
anova(rpv_mod7, rpv_mod8, test = "Chi") # no
summary(rpv_mod8)


#### RPV values ####

rpv_dat %>%
  group_by(disease) %>%
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
anova(co_mod2, co_mod3, test = "Chi") # yes

co_mod4 <- update(co_mod3, ~. -N_added)
summary(co_mod4)
anova(co_mod3, co_mod4, test = "Chi") # no


#### coinfection values ####

co_dat %>%
  group_by(soil) %>%
  summarise(mean_co = mean(co)) %>%
  pivot_wider(names_from = soil,
              values_from = mean_co) %>%
  rename("ambient" = "ambient N",
         "low" = "low N",
         "high" = "high N") %>%
  transmute(a_change = ambient - sterile,
            l_change = low - sterile,
            h_change = high - sterile) %>%
  pivot_longer(everything()) %>%
  summarise(mean_change = mean(value))

co_dat %>%
  group_by(N_added) %>%
  summarise(mean_co = mean(co)) %>%
  pivot_wider(names_from = N_added,
              values_from = mean_co) %>%
  rename("low" = "0",
         "high" = "1") %>%
  mutate(change = high - low)



#### output ####
save(pav_mod4, file = "output/pav_model_rounded_up.rda")
save(rpv_mod7, file = "output/rpv_model_rounded_up.rda")
save(co_mod3, file = "output/coinfection_model_rounded_up.rda")

write_csv(tidy(pav_mod4), "output/pav_model_rounded_up.csv")
write_csv(tidy(rpv_mod7), "output/rpv_model_rounded_up.csv")
write_csv(tidy(co_mod3), "output/coinfection_model_rounded_up.csv")
