# Drought experiment test code
# Yue Yu
# 2025 March 11th

# -- getwd()
getwd()
setwd("/Users/yueyu/Desktop/Drought_R")

# --load pkgs
library(ggplot2)
library(tidyverse)
library(ggpubr)

# -- Load Treatment A/C biomass data
biomass <- read.delim("TreamentAC_Biomass.txt",header = T, na.string = c("","NA"))
lapply(biomass,class)
biomass$biomass <- as.numeric(biomass$biomass)

head(biomass)

# -- Calc Mean among all Replication for each SAM line per Treatment 
new_df <- biomass %>% 
          group_by(Treatment, SAM) %>% 
          summarise(mean_biomass = mean(biomass, na.rm = TRUE)) %>%
          ungroup()

dim(new_df)
head(new_df)

# -- Remove any line with NaN (no data recorded)
# -- In this case SAM281 Control does not have any data
new_df <- new_df[new_df$mean_biomass != "NaN",]


# -- Add Resistant and Susceptible to new column "Type"
Resistant_SAM <- c("SAM285","SAM083","SAM034","SAM254","SAM281","SAM259","SAM287","SAM264","SAM269","SAM075")
Susceptible_SAM <- c("SAM124","SAM051","SAM145","SAM035","SAM063","SAM172","SAM041","SAM139","SAM088","SAM163")

new_df$Type <- ifelse(new_df$SAM %in% Resistant_SAM, "Resistant",
                       ifelse(new_df$SAM %in% Susceptible_SAM, "Susceptible", NA))







# -- Subset Treatment A
a <- new_df[new_df$Treatment == "A",]
SAM_in_A <- unique(a$SAM) 
SAM_in_A

# -- Subset Treatment C 
c <- new_df[new_df$Treatment == "C",]
SAM_in_C <- unique(c$SAM) 
SAM_in_C

# -- Overlap # lines between A and C from R1-3:  17 lines 
overlap <- intersect(SAM_in_A,SAM_in_C)
overlap
length(overlap)




# -- Check mean biomass summary
summary(a$mean_biomass)
summary(c$mean_biomass)




# ---------------------------------------------------
# --  Perform 1st t-test: Overall Treatment A and C
# ---------------------------------------------------
t_test_result <- t.test(a$mean_biomass, c$mean_biomass, alternative = "two.sided", var.equal = FALSE, na.rm = TRUE)
print(t_test_result)
# YEAHHH!!! reject null hypotheses

# -- Extract p-value
p_value <- round(t_test_result$p.value, 5)
p_value


# -- Combine data into one dataframe (LONG format)
biomass_data <- data.frame(
  mean_biomass = c(a$mean_biomass, c$mean_biomass),
  Treatment = rep(c("Treatment_A", "Treatment_C"), times = c(length(a$mean_biomass), length(c$mean_biomass)))
)

# -- Plot boxplot
ggplot(biomass_data, aes(x = Treatment, y = mean_biomass, fill = Treatment)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Overall Mean Biomass compare between Treatment A and C",
       x = "Treatment",
       y = "Mean Biomass") +
  theme_minimal() +
  scale_fill_manual(values = c("Treatment_A" = "lightsalmon2", "Treatment_C" = "royalblue3")) +
  annotate("text", x = 1.5, y = max(biomass_data$mean_biomass, na.rm = TRUE) + 0.5, 
           label = paste("p =", p_value), size = 5, fontface = "bold")



# -------------------------------------------------------------------------
# --  Perform 2nd t-test: Resistant VS Susceptible within Treatment A 
# -------------------------------------------------------------------------

res_a <- a[a$Type == "Resistant",];res_a
sus_a <- a[a$Type == "Susceptible",];sus_a

t_test_result_2 <- t.test(res_a$mean_biomass, sus_a$mean_biomass, alternative = "two.sided", var.equal = FALSE, na.rm = TRUE)
print(t_test_result_2)

# -- Extract p-value
p_value_2 <- round(t_test_result_2$p.value, 3)
p_value_2


# -- Combine data into one dataframe (LONG format)
A_RES_SUS_biomass_data <- data.frame(
  mean_biomass = c(res_a$mean_biomass, sus_a$mean_biomass),
  Type = rep(c("Resistant", "Susceptible"), times = c(length(res_a$mean_biomass), length(sus_a$mean_biomass)))
)


# Ensure "Susceptible" is on the left and "Resistant" is on the right
A_RES_SUS_biomass_data$Type <- factor(A_RES_SUS_biomass_data$Type, levels = c("Susceptible", "Resistant"))

head(A_RES_SUS_biomass_data)


# -- Plot boxplot
ggplot(A_RES_SUS_biomass_data, aes(x = Type, y = mean_biomass, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Mean Biomass of predicted resistant and susceptible within Treatment A",
       x = "Predicted Type",
       y = "Mean Biomass") +
  theme_minimal() +
  scale_fill_manual(values = c("Resistant" = "lightsalmon4", "Susceptible" = "lightsalmon2")) +
  annotate("text", x = 1.5, y = max(A_RES_SUS_biomass_data$mean_biomass, na.rm = TRUE) + 0.5, 
           label = paste("p =", p_value_2), size = 5, fontface = "bold")



# -- Plot point plot
ggplot(A_RES_SUS_biomass_data, aes(x = Type, y = mean_biomass, fill = Type)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Mean Biomass of predicted resistant and susceptible within Treatment A",
       x = "Predicted Type",
       y = "Mean Biomass") +
  theme_minimal() +
  scale_fill_manual(values = c("Resistant" = "lightsalmon4", "Susceptible" = "lightsalmon2")) +
  annotate("text", x = 1.5, y = max(A_RES_SUS_biomass_data$mean_biomass, na.rm = TRUE) + 0.5, 
           label = paste("p =", p_value_2), size = 5, fontface = "bold")







# -------------------------------------------------------------------------
# --  Perform 3rd t-test: Resistant VS Susceptible within Treatment C 
# -------------------------------------------------------------------------

res_c <- c[c$Type == "Resistant",];res_c
sus_c <- c[c$Type == "Susceptible",];sus_c

t_test_result_3 <- t.test(res_c$mean_biomass, sus_c$mean_biomass, alternative = "two.sided", var.equal = FALSE, na.rm = TRUE)
print(t_test_result_3)

# -- Extract p-value
p_value_3 <- round(t_test_result_3$p.value, 3)
p_value_3


# -- Combine data into one dataframe (LONG format)
C_RES_SUS_biomass_data <- data.frame(
  mean_biomass = c(res_c$mean_biomass, sus_c$mean_biomass),
  Type = rep(c("Resistant", "Susceptible"), times = c(length(res_c$mean_biomass), length(sus_c$mean_biomass)))
)


# Ensure "Susceptible" is on the left and "Resistant" is on the right
C_RES_SUS_biomass_data$Type <- factor(C_RES_SUS_biomass_data$Type, levels = c("Susceptible", "Resistant"))

head(C_RES_SUS_biomass_data)


# -- Plot boxplot
ggplot(C_RES_SUS_biomass_data, aes(x = Type, y = mean_biomass, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Mean Biomass of predicted resistant and susceptible within Treatment C",
       x = "Predicted Type",
       y = "Mean Biomass") +
  theme_minimal() +
  scale_fill_manual(values = c("Resistant" = "royalblue4", "Susceptible" = "royalblue2")) +
  annotate("text", x = 1.5, y = max(C_RES_SUS_biomass_data$mean_biomass, na.rm = TRUE) + 0.5, 
           label = paste("p =", p_value_3), size = 5, fontface = "bold")



# -- Plot point plot
ggplot(C_RES_SUS_biomass_data, aes(x = Type, y = mean_biomass, fill = Type)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Mean Biomass of predicted resistant and susceptible within Treatment C",
       x = "Predicted Type",
       y = "Mean Biomass") +
  theme_minimal() +
  scale_fill_manual(values = c("Resistant" = "royalblue4", "Susceptible" = "royalblue2")) +
  annotate("text", x = 1.5, y = max(C_RES_SUS_biomass_data$mean_biomass, na.rm = TRUE) + 0.5, 
           label = paste("p =", p_value_3), size = 5, fontface = "bold")




# ------------------------------------------------------------------------------------
# --  Perform 4th t-test: Resistant VS Susceptible within Treatment A compared to C 
# ------------------------------------------------------------------------------------

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# Step 1: Combine dataframes a and c
combined <- bind_rows(a, c)

# Step 2: Create a new group variable
combined <- combined %>%
  mutate(Group = paste(Treatment, Type, sep = "_"))

# Step 3: Define the order of groups (optional but useful for consistent plots)
combined$Group <- factor(combined$Group, levels = c("A_Susceptible", "A_Resistant", "C_Susceptible", "C_Resistant"))


# Step 3.1: Define custom colors for each group
group_colors <- c(
  "A_Susceptible" = "#8093d7", # BLUE light
  "A_Resistant" = "#e7b8b8",   # RED light
  "C_Susceptible" = "#acacac", # Grey light
  "C_Resistant" = "#7a7575"    # Grey slightly darker
)


# Step 4: Create the boxplot with statistical comparisons
ggplot(combined, aes(x = Group, y = mean_biomass, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  stat_compare_means(
    comparisons = list(
      c("A_Susceptible", "A_Resistant"),
      c("C_Susceptible", "C_Resistant"),
      c("A_Susceptible", "C_Susceptible"),
      c("A_Resistant", "C_Resistant")
    ),
    method = "t.test",
    label = "p.signif", # show significance levels like *, **, *** 
    size = 5,             # font size of the asterisks
    color = "black",        # color of the label
    fontface = "bold" 
    ) +
  scale_fill_manual(values = group_colors) +
  ylim(c(1,6)) +
  labs(x = " ",
       y = " ") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),  # X axis labels
    axis.text.y = element_text(size = 14))


# END
