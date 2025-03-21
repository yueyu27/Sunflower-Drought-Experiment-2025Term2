# Drought experiment test code - LMM
# Yue Yu
# 2025 March 21st


library(lme4)
library(tidyverse)

# ------- Linear Model (LM): 
# This tests the effect of treatment while adjusting for genotype and block effects

# DO NOT USE FOR OUR EXPERIMENT SET UP
# because block should be incorporated as randome effect

model <- lm(biomass ~ treatment + genotype + (1 | block), data = mydata)
summary(model)



# ------- Linear Mixed Model (LMM) - Full model
# If you want to see if treatment affects biomass differently by genotype, you can include an interaction term in the model

LM_model_full <- lmer(biomass ~ treatment * genotype + (1 | replication) + (1 | replication:subblock), data = mydata)
summary(LM_model_full)


# Fixed effects:
#   - treatment * genotype: Tests whether treatment affects biomass differently by genotype (interaction).
#   - treatment: Tests for the main effect of treatment (overall difference between treatment and control).
#   - genotype: Tests for differences between genotypes across conditions.

# Random effects:
#   - (1 | replication): Accounts for variation between replication blocks (R1, R2, R3).
#   - (1 | replication:subblock): Accounts for variation within each sub-block (treatment and control within a block).





# ------- Linear Mixed Model (LMM) - Average among rep model

mydata_avg <- mydata %>%
  group_by(genotype, treatment, replication) %>%
  summarize(biomass_mean = mean(biomass), .groups = "drop")


LM_model_average <- lmer(biomass_mean ~ genotype * treatment + (1 | replication), data = mydata_avg)
summary(LM_model_average)


# Fixed effects:
#   - treatment * genotype: Tests whether treatment affects biomass differently by genotype (interaction).
#   - treatment: Tests for the main effect of treatment (overall difference between treatment and control).
#   - genotype: Tests for differences between genotypes across conditions.

# Random effects:
#   - (1 | replication): Accounts for variation between replication blocks (R1, R2, R3).







# END
