# Fixed NMDS Plot Script
# This script fixes the "Too few points to calculate an ellipse" error
# by filtering groups with insufficient data before drawing ellipses

library(dplyr)
library(ggplot2)

# --- 1. Ensure Factors are Correct ---
# Make sure Age is a factor so we can get levels for the annotation
ord_data$Age <- factor(ord_data$Age, levels = c("Fetus", "Newborn", "7+days", "15+days", "20+days"))
ord_data$SampleType <- as.factor(ord_data$SampleType)

# --- 2. Create a Subset for Ellipses ---
# Filter to keep only groups that have at least 4 samples (required for stat_ellipse)
ord_data_ellipse <- ord_data %>%
  group_by(Age, SampleType) %>%
  filter(n() >= 4) %>%
  ungroup()

# --- 3. Define Annotation Text ---
# (Using the values you calculated)
permanova_text <- paste(
  "PERMANOVA (Bray-Curtis):",
  "\nNMDS Stress =", round(ord$stress, 3),
  "\nAge: R² = 0.211, P = 0.001***",
  "\nSample Type: R² = 0.087, P = 0.003**",
  sep = ""
)

# --- 4. Plot ---
p_nmds <- ggplot(ord_data, aes(x = NMDS1, y = NMDS2, color = SampleType)) +
    
    # POINTS: Plot ALL points (moved shape here to avoid warning in stat_ellipse)
    geom_point(aes(shape = Feed), size = 3, alpha = 0.8) +
    
    # ELLIPSES: Use the FILTERED data (only groups with n >= 4)
    stat_ellipse(data = ord_data_ellipse, 
                 aes(group = SampleType), 
                 linetype = 2, linewidth = 0.6) +
    
    # Facet by Age
    facet_wrap(~ Age) +
    
    # Labels and Theme
    labs(title = "NMDS Ordination of Microbial Communities",
         x = "NMDS1", y = "NMDS2",
         color = "Sample Type", shape = "Feed Status") +
    
    # Annotation (placed in the first available facet)
    geom_text(x = -Inf, y = Inf,
              label = permanova_text, 
              hjust = -0.05, vjust = 1.05, 
              size = 3.5, 
              data = data.frame(Age = levels(ord_data$Age)[1]), 
              inherit.aes = FALSE) +
    
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "gray95"))

print(p_nmds)

# Optional: Save the plot
ggsave("NMDS_plot_fixed.png", p_nmds, width = 12, height = 8, dpi = 300)
