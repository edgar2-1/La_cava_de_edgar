# =============================================================================
# COMPLETE WORKFLOW: From Data Loading to Publication Figures
# Kohler Microbiome Study
# Author: Edgardo L. Rosado Ramos
# Date: 2025-11-20
# =============================================================================

# This script runs the COMPLETE analysis from scratch
# It assumes you have the raw data files in your working directory

# =============================================================================
# STEP 0: INSTALL REQUIRED PACKAGES (RUN ONCE)
# =============================================================================

# Uncomment and run these lines if you need to install packages:
# install.packages(c("phyloseq", "vegan", "tidyverse", "ggplot2", "patchwork", "ggpubr"))
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "ANCOMBC"))

# =============================================================================
# STEP 1: LOAD LIBRARIES
# =============================================================================

library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(DESeq2)

# Optional (for advanced compositional analysis):
# library(ANCOMBC)

cat("✓ Libraries loaded\n")

# =============================================================================
# STEP 2: LOAD AND PREPARE DATA
# =============================================================================

# Read OTU table
otu_raw <- read.table("kohler.opti_mcc.0.03.subsample.shared", header = TRUE, sep = "\t")
otu_mat <- otu_raw[, -(1:3)]
rownames(otu_mat) <- otu_raw$Group 
otu_mat <- t(otu_mat)
otu_tab <- otu_table(otu_mat, taxa_are_rows = TRUE)

cat("✓ OTU table loaded:", nrow(otu_tab), "taxa,", ncol(otu_tab), "samples\n")

# Read taxonomy
tax_raw <- read.table("kohler.opti_mcc.0.03.cons.taxonomy", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
taxonomy_strings <- tax_raw$Taxonomy
names(taxonomy_strings) <- tax_raw$OTU

taxonomy_clean <- lapply(strsplit(taxonomy_strings, ";"), function(x) {
  x <- x[nchar(x) > 0]
  x <- gsub("\\(.*?\\)|\\\"", "", x)
  trimws(x)
})

max_ranks <- 7
taxonomy_mat <- do.call(rbind, lapply(taxonomy_clean, `length<-`, max_ranks))
colnames(taxonomy_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(taxonomy_mat) <- names(taxonomy_strings)
tax_tab <- tax_table(as.matrix(taxonomy_mat))

cat("✓ Taxonomy loaded\n")

# Read metadata
Kohler.design <- read.table("Kohlerdesign.txt",
                            header = TRUE,
                            sep = "\t",
                            row.names = 1,
                            stringsAsFactors = TRUE,
                            check.names = FALSE)

cat("✓ Metadata loaded:", nrow(Kohler.design), "samples\n")

# Create phyloseq object
common_samples <- intersect(sample_names(otu_tab), rownames(Kohler.design))
otu_tab_pruned <- prune_samples(common_samples, otu_tab)
sample_data_obj <- sample_data(Kohler.design[common_samples, ])
physeq <- phyloseq(otu_tab_pruned, tax_tab, sample_data_obj)

cat("✓ Phyloseq object created:", ntaxa(physeq), "taxa,", nsamples(physeq), "samples\n")

# Filter and rarefy
physeq_filt <- filter_taxa(physeq, function(x) sum(x > 3) > (0.2 * length(x)), TRUE)
min_depth <- min(sample_sums(physeq_filt))
set.seed(42)
physeq_rarefied <- rarefy_even_depth(physeq_filt, sample.size = min_depth, replace = FALSE)

cat("✓ Data rarefied to", min_depth, "reads per sample\n")

# Calculate alpha diversity
alpha_diversity_rarefied <- estimate_richness(physeq_rarefied, measures = c("Observed", "Shannon"))
sample_metadata_rarefied <- as.data.frame(sample_data(physeq_rarefied))
alpha_data_rarefied <- cbind(alpha_diversity_rarefied, sample_metadata_rarefied)

# Order factors
age_order <- c("Fetus", "Newborn", "7+days", "15+days", "20+days")
alpha_data_rarefied$Age <- factor(alpha_data_rarefied$Age, levels = age_order)

cat("✓ Alpha diversity calculated\n\n")

# =============================================================================
# STEP 3: SAMPLE SIZE ANALYSIS
# =============================================================================

cat("=== SAMPLE SIZE DISTRIBUTION ===\n")
sample_summary <- alpha_data_rarefied %>%
  group_by(Age, SampleType) %>%
  summarise(n = n(), .groups = "drop")

print(sample_summary)

# Create Figure 1: Sample Distribution
p_sample_dist <- ggplot(sample_summary, aes(x = Age, y = SampleType, fill = n)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = n), size = 6, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#d73027", high = "#1a9850", 
                      name = "Sample\nSize",
                      limits = c(0, 10),
                      breaks = c(0, 2, 4, 6, 8, 10)) +
  labs(title = "Figure 1: Sample Distribution Across Experimental Groups",
       x = "Developmental Stage", y = "Treatment") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Fig1_SampleDistribution.pdf", p_sample_dist, width = 8, height = 5)
cat("\n✓ Figure 1 saved: Fig1_SampleDistribution.pdf\n\n")

# =============================================================================
# STEP 4: BETA DIVERSITY ANALYSIS
# =============================================================================

cat("=== BETA DIVERSITY ANALYSIS ===\n")

# Calculate Bray-Curtis distance
bray_dist <- phyloseq::distance(physeq_rarefied, method = "bray")
metadata <- as.data.frame(sample_data(physeq_rarefied))

# Standard PERMANOVA
set.seed(42)
permanova_standard <- adonis2(
  bray_dist ~ SampleType + Age + Mother,
  data = metadata,
  permutations = 9999,
  by = "margin"
)

cat("\nPERMANOVA Results:\n")
print(permanova_standard)

# PERMDISP (check assumptions)
dispersion_sampletype <- betadisper(bray_dist, metadata$SampleType)
dispersion_age <- betadisper(bray_dist, metadata$Age)
dispersion_mother <- betadisper(bray_dist, metadata$Mother)

cat("\nPERMDISP Results (checking homogeneity of dispersion):\n")
cat("SampleType: p =", anova(dispersion_sampletype)$`Pr(>F)`[1], "\n")
cat("Age: p =", anova(dispersion_age)$`Pr(>F)`[1], "\n")
cat("Mother: p =", anova(dispersion_mother)$`Pr(>F)`[1], "\n")

# Create effect size table
effect_size_table <- data.frame(
  Factor = c("Mother", "SampleType", "Age"),
  R2 = c(
    permanova_standard$R2[which(rownames(permanova_standard) == "Mother")],
    permanova_standard$R2[which(rownames(permanova_standard) == "SampleType")],
    permanova_standard$R2[which(rownames(permanova_standard) == "Age")]
  ),
  F_stat = c(
    permanova_standard$F[which(rownames(permanova_standard) == "Mother")],
    permanova_standard$F[which(rownames(permanova_standard) == "SampleType")],
    permanova_standard$F[which(rownames(permanova_standard) == "Age")]
  ),
  p_value = c(
    permanova_standard$`Pr(>F)`[which(rownames(permanova_standard) == "Mother")],
    permanova_standard$`Pr(>F)`[which(rownames(permanova_standard) == "SampleType")],
    permanova_standard$`Pr(>F)`[which(rownames(permanova_standard) == "Age")]
  )
) %>%
  mutate(
    R2_percent = round(R2 * 100, 1),
    Effect_Size = case_when(
      R2 < 0.01 ~ "Negligible",
      R2 < 0.06 ~ "Small",
      R2 < 0.14 ~ "Medium",
      TRUE ~ "Large"
    ),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

cat("\nEffect Size Summary:\n")
print(effect_size_table)

write.csv(effect_size_table, "Table1_EffectSizes.csv", row.names = FALSE)
cat("\n✓ Table 1 saved: Table1_EffectSizes.csv\n\n")

# =============================================================================
# STEP 5: FIGURE 2 - MATERNAL EFFECTS (3-PANEL)
# =============================================================================

cat("=== CREATING FIGURE 2: MATERNAL EFFECTS ===\n")

# Panel A: PCoA colored by Mother
pcoa_res <- ordinate(physeq_rarefied, method = "PCoA", distance = bray_dist)

pcoa_mother <- plot_ordination(physeq_rarefied, pcoa_res, 
                               color = "Mother", shape = "SampleType") +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = Mother), linetype = 2, linewidth = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "A) Community Composition by Maternal Clutch",
       subtitle = paste0("PERMANOVA: R² = ", 
                        round(effect_size_table$R2[1], 3), 
                        ", p = ", 
                        format.pval(effect_size_table$p_value[1], digits = 3)),
       x = paste0("PCoA1 (", round(pcoa_res$values$Relative_eig[1] * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(pcoa_res$values$Relative_eig[2] * 100, 1), "%)")) +
  theme_bw(base_size = 12)

# Panel B: Dispersion within clutches
# NOTE: dispersion_mother was created in STEP 4 above (line with betadisper)
dispersion_df <- data.frame(
  Mother = metadata$Mother,
  Distance = dispersion_mother$distances
)

p_dispersion <- ggplot(dispersion_df, aes(x = Mother, y = Distance, fill = Mother)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "B) Within-Clutch Microbiome Variation",
       y = "Distance to Clutch Centroid (Bray-Curtis)",
       x = "Maternal ID") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Panel C: Shannon diversity by Mother
p_alpha_mother <- ggplot(alpha_data_rarefied, 
                         aes(x = Mother, y = Shannon, fill = Mother)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_brewer(palette = "Set1") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               fill = "white", color = "black") +
  labs(title = "C) Alpha Diversity by Maternal Clutch",
       y = "Shannon Diversity Index",
       x = "Maternal ID") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Combine panels
fig2 <- pcoa_mother / (p_dispersion | p_alpha_mother) +
  plot_annotation(
    title = "Figure 2: Maternal Identity Shapes Offspring Gut Microbiomes",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave("Fig2_MaternalEffects.pdf", fig2, width = 12, height = 10)
cat("✓ Figure 2 saved: Fig2_MaternalEffects.pdf\n\n")

# =============================================================================
# STEP 6: FIGURE 3 - ENVIRONMENTAL EFFECTS
# =============================================================================

cat("=== CREATING FIGURE 3: ENVIRONMENTAL EFFECTS ===\n")

# Focus on well-powered age groups (7+days and 15+days)
physeq_subset <- subset_samples(physeq_rarefied, 
                                Age %in% c("7+days", "15+days"))

bray_subset <- phyloseq::distance(physeq_subset, method = "bray")
pcoa_subset <- ordinate(physeq_subset, method = "PCoA", distance = bray_subset)

# Panel A: PCoA by SampleType
pcoa_env <- plot_ordination(physeq_subset, pcoa_subset, 
                            color = "SampleType", shape = "Age") +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = SampleType), linetype = 2, linewidth = 0.8) +
  scale_color_manual(values = c("Muck" = "#8c510a", "Sterile" = "#01665e")) +
  labs(title = "A) Community Composition by Environmental Exposure",
       subtitle = paste0("PERMANOVA: R² = ", 
                        round(effect_size_table$R2[2], 3), 
                        ", p = ", 
                        format.pval(effect_size_table$p_value[2], digits = 3),
                        " (7+days and 15+days only)"),
       x = paste0("PCoA1 (", round(pcoa_subset$values$Relative_eig[1] * 100, 1), "%)"),
       y = paste0("PCoA2 (", round(pcoa_subset$values$Relative_eig[2] * 100, 1), "%)")) +
  theme_bw(base_size = 12)

# Panel B: Alpha diversity by SampleType
alpha_subset <- alpha_data_rarefied %>%
  filter(Age %in% c("7+days", "15+days"))

p_alpha_env <- ggplot(alpha_subset, aes(x = SampleType, y = Shannon, 
                                        fill = SampleType)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("Muck" = "#8c510a", "Sterile" = "#01665e")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(title = "B) Shannon Diversity by Treatment",
       subtitle = "7+days and 15+days only",
       y = "Shannon Diversity Index",
       x = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

fig3 <- pcoa_env | p_alpha_env
ggsave("Fig3_EnvironmentalEffects.pdf", fig3, width = 14, height = 6)
cat("✓ Figure 3 saved: Fig3_EnvironmentalEffects.pdf\n\n")

# =============================================================================
# STEP 7: SUMMARY STATISTICS TABLE
# =============================================================================

cat("=== CREATING SUMMARY STATISTICS TABLE ===\n")

summary_stats <- alpha_data_rarefied %>%
  group_by(SampleType, Age) %>%
  summarise(
    n = n(),
    Shannon_mean = mean(Shannon, na.rm = TRUE),
    Shannon_sd = sd(Shannon, na.rm = TRUE),
    Shannon_se = sd(Shannon, na.rm = TRUE) / sqrt(n()),
    Observed_mean = mean(Observed, na.rm = TRUE),
    Observed_sd = sd(Observed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Shannon_summary = paste0(round(Shannon_mean, 2), " ± ", round(Shannon_se, 2)),
    Observed_summary = paste0(round(Observed_mean, 1), " ± ", round(Observed_sd, 1))
  )

print(summary_stats)
write.csv(summary_stats, "Table2_SummaryStats.csv", row.names = FALSE)
cat("\n✓ Table 2 saved: Table2_SummaryStats.csv\n\n")

# =============================================================================
# STEP 8: SAVE SESSION INFO
# =============================================================================

sink("SessionInfo.txt")
print(sessionInfo())
sink()

cat("✓ Session info saved: SessionInfo.txt\n\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║           ANALYSIS COMPLETE - FILES GENERATED              ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

cat("📊 FIGURES:\n")
cat("  ✓ Fig1_SampleDistribution.pdf\n")
cat("  ✓ Fig2_MaternalEffects.pdf (3-panel)\n")
cat("  ✓ Fig3_EnvironmentalEffects.pdf (2-panel)\n\n")

cat("📋 TABLES:\n")
cat("  ✓ Table1_EffectSizes.csv\n")
cat("  ✓ Table2_SummaryStats.csv\n\n")

cat("📄 OTHER:\n")
cat("  ✓ SessionInfo.txt\n\n")

cat("🎯 KEY FINDINGS:\n")
cat("  • Maternal effects: R² =", round(effect_size_table$R2[1], 3), 
    ", p =", format.pval(effect_size_table$p_value[1], digits = 3), "\n")
cat("  • SampleType effects: R² =", round(effect_size_table$R2[2], 3), 
    ", p =", format.pval(effect_size_table$p_value[2], digits = 3), "\n")
cat("  • Age effects: R² =", round(effect_size_table$R2[3], 3), 
    ", p =", format.pval(effect_size_table$p_value[3], digits = 3), "\n\n")

cat("📝 NEXT STEPS:\n")
cat("  1. Review the generated figures\n")
cat("  2. Check the effect size table\n")
cat("  3. Start writing your manuscript\n")
cat("  4. Target journals: Microorganisms, Frontiers, or PeerJ\n\n")

cat("═══════════════════════════════════════════════════════════════\n")
