# Drought experiment test code - Days to wilt
# Yue Yu
# 2025 April09

# -- getwd()
getwd()
setwd("/Users/yueyu/Desktop")

# --load pkgs
library(tidyverse)
library(ggplot2)

# -- Load Treatment B days to full wilt data
wilt <- read.delim("TreatmentB_2025April09.txt",header = T, na.string = c("","NA"))
lapply(wilt,class)


# -- Adjust “days to wilt” using “V6_area"

# -------- First wilt ----
# Use linear model to adjust
model <- lm(Days_to_first_wilt ~ V6_area, data = wilt)

# Averge leaf area at V6
mean_v6 <- mean(wilt$V6_area, na.rm = TRUE); mean_v6

# Baseline days after accounting for V6_area
baseline_days <- coef(model)[1] + coef(model)[2] * mean_v6; baseline_days

# Adjusted first wilt (make biological sense)
wilt$adjusted_first_wilt <- residuals(model) + baseline_days



# -------- Full wilt ----
model2 <- lm(Days_to_full_wilt ~ V6_area, data = wilt)
mean_v6 <- mean(wilt$V6_area, na.rm = TRUE); mean_v6
baseline_days <- coef(model2)[1] + coef(model2)[2] * mean_v6; baseline_days
wilt$adjusted_full_wilt <- residuals(model2) + baseline_days


# -- Summarize adjusted data

summary(wilt$adjusted_first_wilt)
  Min.   1st Qu.  Median  Mean  3rd Qu.    Max. 
 0.3482  7.4425  9.2067  9.6545 10.8148 18.9710 

summary(wilt$adjusted_full_wilt)
  Min.   1st Qu.  Median  Mean   3rd Qu.    Max. 
  8.116  11.451  14.184  14.255  16.290  24.214 



# -- Calc Mean among all Replication per SAM line
new_df <- wilt %>% 
  group_by(SAM) %>% 
  summarise(mean_first_wilt = mean(adjusted_first_wilt, na.rm = TRUE),
  	        mean_full_wilt = mean(adjusted_full_wilt, na.rm = TRUE)) %>%
  ungroup()

dim(new_df)
head(new_df)


# -- Add Resistant and Susceptible to new column "Type"
Resistant_SAM <- c("SAM285","SAM083","SAM034","SAM254","SAM281","SAM259","SAM287","SAM264","SAM269","SAM075")
Susceptible_SAM <- c("SAM124","SAM051","SAM145","SAM035","SAM063","SAM172","SAM041","SAM139","SAM088","SAM163")

new_df$Type <- ifelse(new_df$SAM %in% Resistant_SAM, "Resistant",
                       ifelse(new_df$SAM %in% Susceptible_SAM, "Susceptible", NA))



# -------------------------------------------------------------------------
# --  Perform wilt t-test: Resistant VS Susceptible within Treatment B for WILT
# -------------------------------------------------------------------------

res_b <- new_df[new_df$Type == "Resistant",];res_b
sus_b <- new_df[new_df$Type == "Susceptible",];sus_b


t_first <- t.test(res_b$mean_first_wilt, sus_b$mean_first_wilt, alternative = "two.sided", var.equal = FALSE, na.rm = TRUE)
print(t_first)


t_full <- t.test(res_b$mean_full_wilt, sus_b$mean_full_wilt, alternative = "two.sided", var.equal = FALSE, na.rm = TRUE)
print(t_full)



# -- Extract p-value
p_value_b <- round(t_full$p.value, 3)
p_value_b

# 0.5 first wilt and  0.7 full wilt
# very sad indeed


# -- Combine data into one dataframe (LONG format)
FULL <- data.frame(
  Mean_Wilt = c(sus_b$mean_full_wilt,res_b$mean_full_wilt),
  Type = rep(c("Susceptible","Resistant"), times = c(length(sus_b$mean_full_wilt),length(res_b$mean_full_wilt)))
)

# Ensure "Susceptible" is on the left and "Resistant" is on the right
FULL$Type <- factor(FULL$Type, levels = c("Susceptible", "Resistant"))
head(FULL)


# -- Plot boxplot + points
ggplot(FULL, aes(x = Type, y = Mean_Wilt, fill = Type)) +
  geom_boxplot() +
  geom_point() +
  theme_minimal() +
  labs(title = "Mean full wilt of R and S within Treatment B (adjusted by V6 leaf area)",
       x = "Predicted Type",
       y = "Days to full wilt (mean)") +
  theme_minimal() +
  scale_fill_manual(values = c("Susceptible" = "lightseagreen","Resistant" = "aquamarine4")) +
  annotate("text", x = 1.5, y = max(FULL$Mean_Wilt, na.rm = TRUE) + 0.5, 
           label = paste("p =", p_value_b), size = 5, fontface = "bold")



# END
