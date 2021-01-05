##### info ####

# goal: analyze chlorophyll data when only PCR bands at least as intense as controls are counted as indicators of infection (rounding down)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(broom)

# import data
dat <- read_csv("intermediate-data/bydv_microbes_data_rounded_down.csv")


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
         log_chlorophyll = log(chlorophyll)) %>%
  filter(!is.na(infection))

# replicates
dat2 %>%
  group_by(infection, N_added, soil) %>%
  count() %>%
  data.frame()
# co-infection missing some treatments

# remove coinfection for stats
clr_dat <- dat2 %>%
  filter(infection != "Co-infection")


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
pan_lab <- tibble(infection = unique(dat2$infection),
                  label = c("(a)", "(b)", "(c)", "(d)")) %>%
  mutate(soil = "sterile",
         chlorophyll = 33,
         nitrogen_added = "low")

# chlorophyll figure
pdf("output/chlorophyll_figure_rounded_down.pdf", width = 6, height = 6)
ggplot(dat2, aes(soil, chlorophyll, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2), shape = 21) +
  geom_text(data = pan_lab, aes(label = label), nudge_x = -0.4, fontface = "bold") +
  facet_wrap(~ infection) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  guides(fill = guide_legend(override.aes = list(shape = 21), direction = "horizontal", title.position = "top")) +
  xlab("Field soil N treatment") +
  ylab("Leaf chlorophyll content [SPAD]") +
  theme_def +
  theme(legend.position = c(0.15, 0.6))
dev.off()


#### chlorophyll model ####

# full model
clr_mod1 <- lm(log_chlorophyll ~ soil * N_added * infection,
               data = clr_dat)
summary(clr_mod1)

# remove 3-way interaction?
clr_mod2 <- update(clr_mod1, ~. -soil:N_added:infection)
summary(clr_mod2)
anova(clr_mod1, clr_mod2, test = "F") # yes

# remove 2-way interactions?
clr_mod3 <- update(clr_mod2, ~. -N_added:infection)
summary(clr_mod3)
anova(clr_mod2, clr_mod3, test = "F") # yes

clr_mod4 <- update(clr_mod3, ~. -soil:infection)
summary(clr_mod4)
anova(clr_mod3, clr_mod4, test = "F") # yes

clr_mod5 <- update(clr_mod4, ~. -soil:N_added)
summary(clr_mod5)
anova(clr_mod4, clr_mod5, test = "F") # no

# remove main effects?
clr_mod6 <- update(clr_mod4, ~. -infection)
summary(clr_mod6)
anova(clr_mod4, clr_mod6, test = "F") # no


#### chlorophyll values ####
clr_dat %>%
  filter(soil == "sterile" & 
           nitrogen_added == "low" & 
           infection == "Mock inoculation") %>%
  summarise(mean_chlr = mean(chlorophyll),
            se_chlr = sd(chlorophyll)/n())

clr_dat %>%
  group_by(soil, nitrogen_added) %>%
  summarise(mean_chlr = mean(chlorophyll)) %>%
  pivot_wider(names_from = "nitrogen_added",
              values_from = "mean_chlr") %>%
  mutate(change = (high - low)/low)

clr_dat %>%
  group_by(soil, nitrogen_added) %>%
  summarise(mean_chlr = mean(chlorophyll)) %>%
  pivot_wider(names_from = "soil",
              values_from = "mean_chlr") %>%
  rename("ambient" = "ambient N",
         "low" = "low N",
         "high" = "high N") %>%
  mutate(ambient_change = (ambient - sterile)/sterile,
         low_change = (low - sterile)/sterile,
         high_change = (high - sterile)/sterile)

clr_dat %>%
  group_by(infection) %>%
  summarise(mean_chlr = mean(chlorophyll)) %>%
  pivot_wider(names_from = "infection",
              values_from = "mean_chlr") %>%
  rename("mock" = "Mock inoculation",
         "PAV" = "PAV infection",
         "RPV" = "RPV infection") %>%
  mutate(PAV_change = (PAV - mock)/mock,
         RPV_change = (RPV - mock)/mock)


#### output ####
save(clr_mod4, file = "output/chlorophyll_model_rounded_down.rda")
write_csv(tidy(clr_mod4), "output/chlorophyll_model_rounded_down.csv")
