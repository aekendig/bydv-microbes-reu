##### info ####

# goal: analyze biomass data any visible PCR band is counted as indicators of infection (rounding up)


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(broom)

# import data
dat <- read_csv("intermediate-data/bydv_microbes_data_rounded_up.csv")


#### edit data ####

# reorganize factor levels
# select for successful infection
dat2 <- dat %>%
  mutate(soil = fct_relevel(soil, "sterile", "ambient N", "low N"),
         infection = case_when(disease == "PAV" & pav == 1 ~ "PAV infection",
                               disease == "RPV" & rpv == 1 ~ "RPV infection",
                               disease == "Co" & pav == 1 & rpv == 1 ~ "Co-infection",
                               disease == "Healthy" ~ "Mock inoculation",
                               TRUE ~ NA_character_) %>%
           fct_relevel("Mock inoculation", "PAV infection", "RPV infection"),
         nitrogen_added = fct_relevel(nitrogen_added, "low", "high"),
         log_biomass = log(biomass)) %>%
  filter(!is.na(infection))

# replicates
dat2 %>%
  group_by(infection, N_added, soil) %>%
  count() %>%
  data.frame()
# co-infection only has 1 rep for some treatments


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
  mutate(soil = "sterile",
         biomass = 0.65,
         nitrogen_added = "low")

# biomass figure
pdf("output/biomass_figure_rounded_up.pdf", width = 6, height = 6)
ggplot(dat2, aes(soil, biomass, fill = nitrogen_added)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", size = 3, position = position_dodge(0.2), shape = 21) +
  geom_text(data = pan_lab, aes(label = label), nudge_x = -0.4, fontface = "bold") +
  facet_wrap(~ infection) +
  scale_fill_manual(values = col_pal, name = "N supply") +
  guides(fill = guide_legend(override.aes = list(shape = 21), direction = "horizontal", title.position = "top")) +
  xlab("Field soil N treatment") +
  ylab("Aboveground biomass [g/plant]") +
  theme_def +
  theme(legend.position = c(0.88, 0.4))
dev.off()


#### biomass model ####

# full model
bio_mod1 <- lm(log_biomass ~ soil * N_added * infection,
               data = dat2)
summary(bio_mod1)

# remove 3-way interaction?
bio_mod2 <- update(bio_mod1, ~. -soil:N_added:infection)
summary(bio_mod2)
anova(bio_mod1, bio_mod2, test = "F") # yes

# remove 2-way interactions?
bio_mod3 <- update(bio_mod2, ~. -soil:N_added)
summary(bio_mod3)
anova(bio_mod2, bio_mod3, test = "F") # yes

bio_mod4 <- update(bio_mod3, ~. -N_added:infection)
summary(bio_mod4)
anova(bio_mod3, bio_mod4, test = "F") # yes

bio_mod5 <- update(bio_mod4, ~. -soil:infection)
summary(bio_mod5)
anova(bio_mod4, bio_mod5, test = "F") # yes

# remove main effects?
bio_mod6 <- update(bio_mod5, ~. -soil)
summary(bio_mod6)
anova(bio_mod5, bio_mod6, test = "F") # yes

bio_mod7 <- update(bio_mod6, ~. -infection)
summary(bio_mod7)
anova(bio_mod6, bio_mod7, test = "F") # no

bio_mod8 <- update(bio_mod6, ~. -N_added)
summary(bio_mod8)
anova(bio_mod6, bio_mod8, test = "F") # no


#### biomass values ####
dat2 %>%
  group_by(nitrogen_added) %>%
  summarise(mean_bio = mean(biomass)) %>%
  pivot_wider(names_from = "nitrogen_added",
              values_from = "mean_bio") %>%
  mutate(change = (high - low)/low)

dat2 %>%
  group_by(infection) %>%
  summarise(mean_bio = mean(biomass)) %>%
  pivot_wider(names_from = "infection",
              values_from = "mean_bio") %>%
  rename("mock" = "Mock inoculation",
         "PAV" = "PAV infection",
         "RPV" = "RPV infection",
         "Co" = "Co-infection") %>%
  mutate(PAV_change = (PAV - mock)/mock,
         RPV_change = (RPV - mock)/mock,
         co_change = (Co - mock)/mock)


#### output ####
save(bio_mod6, file = "output/biomass_model_rounded_up.rda")
write_csv(tidy(bio_mod6), "output/biomass_model_rounded_up.csv")
