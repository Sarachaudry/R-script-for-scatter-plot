# R scatter plot
This R script describe the steps to generate a scatterplot between Insulation Index ("IS") and Replication Timing signal ("RT) at three time points: Early, Mid, Late.

```R
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Set working directory
setwd("/scratch/99999/ha20be/Akram_Basslab_revised/EigenFiles_Compare_JuicevsCScoreTool")

# Read bedGraph files
IS_1kb <- read.table("fanc_RTallBRs_IS_1kb_2.5mbW.bedGraph", header = FALSE)
early <- read.table("E_ratio_2.smooth.1kbp.sorted.bedgraph", header = FALSE)
middle <- read.table("M_ratio_2.smooth.1kbp.sorted.bedgraph", header = FALSE)
late <- read.table("L_ratio_2.smooth.1kbp.sorted.bedgraph", header = FALSE)

# Assign column names
colnames(IS_1kb) <- c("chr", "start", "end", "value")
colnames(early)  <- c("chr", "start", "end", "value")
colnames(middle) <- c("chr", "start", "end", "value")
colnames(late)   <- c("chr", "start", "end", "value")

# Merge data
merged <- IS_1kb %>% rename(IS_1kb = value) %>% inner_join(early, by = c("chr", "start", "end")) %>%
rename(Early = value) %>% inner_join(middle, by = c("chr", "start", "end")) %>% rename(Middle = value) %>%
inner_join(late, by = c("chr", "start", "end")) %>% rename(Late = value)

# Pivot to long format
long_df <- merged %>% pivot_longer(cols = c(Early, Middle, Late), names_to = "RT_class", values_to = "RT_signal")

# Ensure correct factor order
long_df$RT_class <- factor(long_df$RT_class, levels = c("Early", "Middle", "Late"))

# Calculate Spearman correlations, ignoring missing values
cor_table <- long_df %>% group_by(RT_class) %>%
summarize(cor = cor(IS_1kb, RT_signal, method = "spearman", use = "complete.obs")) %>% arrange(RT_class)

# Set annotation positions and text (bottom right corner)
cor_labels <- data.frame(
RT_class = c("Early", "Middle", "Late"),
x = c(3.9, 3.9, 3.9),
y = c(0.95, 0.85, 0.75),
label = paste0("R = ", round(cor_table$cor, 2))
)

# Create plot and save as PDF
pdf("RT_vs_IS_1kb_CombinedPlot_EML_1kb_FINAL.pdf", width = 7, height = 6)

ggplot(long_df, aes(x = RT_signal, y = IS_1kb, color = RT_class)) +
geom_point(alpha = 0.3, size = 1) +
geom_smooth(method = "lm", se = FALSE, aes(group = RT_class), color = "black", size = 0.9) +
theme_minimal(base_size = 14) +
scale_color_manual(values = c("Early" = "blue", "Middle" = "green", "Late" = "red")) +
labs(
x = "RT_1kb",
y = "Hi-C IS_1kb",
title = "RT Class Correlation with Hi-C IS_1kb",
color = "RT Class"
) +
coord_cartesian(xlim = c(0, 5), ylim = c(-1, 1)) +
geom_text(data = cor_labels,
aes(x = x, y = y, label = label, color = RT_class),
inherit.aes = FALSE,
hjust = 0,
size = 4.2)

dev.off()
```
