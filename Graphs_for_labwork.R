library(ggplot2)
library(openxlsx)

#load the data
data_m<-read.xlsx("C:/Users/Олеся/Desktop/graduation/lab/final_exp_results/meristem_len.xlsx")

data_h<-read.csv("C:/Users/Олеся/Desktop/graduation/lab/final_exp_results/stereo_results.csv")

# graph for meristem len
ggplot(data_m, aes(x = treatment, y = len, fill = treatment)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("CYCB1;1" = "#225ea8", "DMSO" = "#41b6c4")) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Meristem length distribution",
       x = "Treatment", y = "Length (μm)") +
  theme_minimal()

ggsave( filename = "C:/Users/Олеся/Desktop/graduation/lab/final_exp_results/Meristem length distributio.png",
        width = 4.5,
        height = 3.5,
        dpi = 300)

# graph for distance to root hairs
ggplot(data_h, aes(x = Sample_type, y = len, fill = Sample_type)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("GFP&RFP" = "#2c7fb8", "DMSO" = "#41b6c4", "GFP" = "#a1dab4")) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Distance to first root hairs",
       x = "Sample Type", y = "Length") +
  theme_minimal()
ggsave( filename = "C:/Users/Олеся/Desktop/graduation/lab/final_exp_results/distance to root hairs.png",
        width = 4.5,
        height = 3.5,
        dpi = 300)

#statistical analysis
# for meristem len 
dex_m <- data_m[data_m$treatment %in% c("DEX"),]
qqnorm(dex_m$len, main =  "Q-Q Plot for meristem lenght for CYCB1;1 samples")
qqline(dex_m$len, col = "red")

DMSO_m <- data_m[data_m$treatment %in% c("DMSO"),]
qqnorm(DMSO_m$len, main = "Q-Q Plot for meristem lenght for DMSO samples")
qqline(DMSO_m$len, col = "red")

anova_m <- aov(len ~ treatment, data = data_m)
summary(anova_m)

#for root hairs distanse
dex_h <- data_h[data_h$Sample_type %in% c("DEX"),]
qqnorm(dex_h$len, , main = "Q-Q Plot for distance to root hairs for RFP&GFP samples")
qqline(dex_h$len, col = "red")

DMSO_h <- data_h[data_h$Sample_type  %in% c("DMSO"),]
qqnorm(DMSO_h$len, main = "Q-Q Plot for distance to root hairs for DMSO samples")
qqline(DMSO_h$len, col = "red")

gfp_h <- data_h[data_h$Sample_type  %in% c("GFP"),]
qqnorm(gfp_h$len, main = "Q-Q Plot for distance to root hairs for GFP samples")
qqline(gfp_h$len, col = "red")


anova_h <- aov(len ~ Sample_type, data = data_h)
summary(anova_h)
TukeyHSD(anova_h)
plot(TukeyHSD(anova_h))


