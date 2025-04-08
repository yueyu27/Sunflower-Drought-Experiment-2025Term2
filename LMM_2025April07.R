# LMM on Final A/C data
# Yue Yu
# 2025 April 07

# ============================
# Install and Load Required Libraries
# ============================
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(tidyverse)

# ============================
# Load Data
# ============================
setwd("/Users/yueyu/Desktop")
biomass <- read.table("Farbod_raw_data.csv", sep =",", header = TRUE, na.string = c("","NA"))
lapply(biomass,class)
head(biomass)

# ============================
# Log-transform Biomass
# ============================
biomass_log <- biomass %>%
  filter(!is.na(Total_shoot_bio_mass_g)) %>%
  mutate(log_biomass = log(Total_shoot_bio_mass_g)) %>%
  select(-Total_shoot_bio_mass_g)


biomass_log$Treatment <- as.factor(biomass_log$Treatment)
biomass_log$Mode <- as.factor(biomass_log$Mode)
biomass_log$Replicate <- as.factor(biomass_log$Replicate)
biomass_log$Sam_Line <- as.factor(biomass_log$Sam_Line)


head(biomass_log)
lapply(biomass_log,class)

# ============================
# LMM
# ============================

LMM <- lmer(log_biomass ~ Treatment + Mode + (1 | Sam_Line) + (1 | Replicate), data = biomass_log)
# warning: boundary (singular) fit: see help('isSingular')
# The warning means that lmer() found the variance of the random effect (Replicate) to be essentially zero, i.e., there is no measurable variation in log_biomass between replicates, once fixed effects are accounted for.
# Possible issue: fewer than 5 levels for random effect, variance estimation for random effects becomes unreliable. If you only have 2 or 3 replicates, that's likely causing the singularity.

summary(LMM)

# results

Estimate Std. Error t value
(Intercept)  0.64595    0.06906   9.353
TreatmentC   0.32062    0.07997   4.009
ModeT       -0.10351    0.07997  -1.294

# Fixed effects:
# Intercept: Mean log_biomass for Treatment A & Mode S (0.64595)
# TreatmentC: Increase compared to A (significant by T value)
# Mode T: Decrease small compared to Mode S (NON - significant by T value)

# END
