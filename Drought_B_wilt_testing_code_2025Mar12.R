# Drought experiment test code - Days to wilt
# Yue Yu
# 2025 March 12th

# -- getwd()
getwd()
setwd("/Users/yueyu/Desktop")

# --load pkgs
library(ggplot2)

# -- Load Treatment B days to full wilt data
wilt <- read.delim("TreatmentB_2025March12.txt",header = T, na.string = c("","NA"))
lapply(wilt,class)
wilt$Days_to_full_wilt <- as.numeric(wilt$Days_to_full_wilt)


# -- Calc Mean among all Replication for each SAM line per Treatment 
new_df <- wilt %>% 
  group_by(Treatment, SAM) %>% 
  summarise(mean_wilt = mean(Days_to_full_wilt, na.rm = TRUE)) %>%
  ungroup()

dim(new_df)
head(new_df)
# 18 lines used for downstream analyses



# -- Add Resistant and Susceptible to new column "Type"
Resistant_SAM <- c("SAM285","SAM083","SAM034","SAM254","SAM281","SAM259","SAM287","SAM264","SAM269","SAM075")
Susceptible_SAM <- c("SAM124","SAM051","SAM145","SAM035","SAM063","SAM172","SAM041","SAM139","SAM088","SAM163")

new_df$Type <- ifelse(new_df$SAM %in% Resistant_SAM, "Resistant",
                       ifelse(new_df$SAM %in% Susceptible_SAM, "Susceptible", NA))


# -- Remove any line with NaN (no data recorded)
# -- In this case SAM281 Control does not have any data
new_df <- new_df[new_df$mean_wilt != "NaN",]
dim(new_df)
head(new_df)


# -------------------------------------------------------------------------
# --  Perform wilt t-test: Resistant VS Susceptible within Treatment B for WILT
# -------------------------------------------------------------------------

res_b <- new_df[new_df$Type == "Resistant",];res_b
sus_b <- new_df[new_df$Type == "Susceptible",];sus_b


t_b <- t.test(res_b$mean_wilt, sus_b$mean_wilt, alternative = "two.sided", var.equal = FALSE, na.rm = TRUE)
print(t_b)

# -- Extract p-value
p_value_b <- round(t_b$p.value, 3)
p_value_b


# -- Combine data into one dataframe (LONG format)
B_RES_SUS <- data.frame(
  Mean_Wilt = c(sus_b$mean_wilt,res_b$mean_wilt),
  Type = rep(c("Susceptible","Resistant"), times = c(length(sus_b$mean_wilt),length(res_b$mean_wilt)))
)

# Ensure "Susceptible" is on the left and "Resistant" is on the right
B_RES_SUS$Type <- factor(B_RES_SUS$Type, levels = c("Susceptible", "Resistant"))

head(B_RES_SUS)


# -- Plot boxplot
ggplot(B_RES_SUS, aes(x = Type, y = Mean_Wilt, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Mean full wilt days of predicted resistant and susceptible within Treatment B",
       x = "Predicted Type",
       y = "Days to full wilt (mean)") +
  theme_minimal() +
  scale_fill_manual(values = c("Susceptible" = "lightseagreen","Resistant" = "aquamarine4")) +
  annotate("text", x = 1.5, y = max(B_RES_SUS$Mean_Wilt, na.rm = TRUE) + 0.5, 
           label = paste("p =", p_value_b), size = 5, fontface = "bold")



# -- Plot Points
ggplot(B_RES_SUS, aes(x = Type, y = Mean_Wilt, fill = Type)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Mean full wilt days of predicted resistant and susceptible within Treatment B",
       x = "Predicted Type",
       y = "Days to full wilt (mean)") +
  theme_minimal() +
  scale_fill_manual(values = c("Susceptible" = "lightseagreen","Resistant" = "aquamarine4")) +
  annotate("text", x = 1.5, y = max(B_RES_SUS$Mean_Wilt, na.rm = TRUE) + 0.5, 
           label = paste("p =", p_value_b), size = 5, fontface = "bold")

