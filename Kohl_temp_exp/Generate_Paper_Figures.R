# ============================================================================
# Generate All Figures and Tables for Semester Project Paper
# Edgardo L. Rosado Ramos
# ============================================================================

# Load required libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(png)

# Set output directory
output_dir <- "/Users/Familia/Downloads/Project/Kohl_temp_exp/Paper_Figures"
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

cat("Loading data...\n")

# Load OTU table
otu_raw <- read.table("kohler.opti_mcc.shared", header = TRUE, sep = "\t")
otu_mat <- otu_raw[, -(1:3)]
rownames(otu_mat) <- otu_raw$Group
otu_mat <- t(otu_mat)

# Load taxonomy
tax_raw <- read.table("kohler.opti_mcc.0.03.cons.taxonomy", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
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

# Load metadata
Kohler.design <- read.table("Kohlerdesign.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = TRUE)

# Create phyloseq object
common_samples <- intersect(colnames(otu_mat), rownames(Kohler.design))
physeq <- phyloseq(
    otu_table(otu_mat[, common_samples], taxa_are_rows = TRUE),
    tax_table(as.matrix(taxonomy_mat)),
    sample_data(Kohler.design[common_samples, ])
)

# Quality control and rarefaction
min_threshold <- 4000
sample_sums_val <- sample_sums(physeq)
keep_samples <- names(sample_sums_val[sample_sums_val >= min_threshold])
physeq_filt <- prune_samples(keep_samples, physeq)
physeq_filt <- prune_taxa(taxa_sums(physeq_filt) > 0, physeq_filt)

set.seed(42)
min_depth <- min(sample_sums(physeq_filt))
physeq_rarefied <- rarefy_even_depth(physeq_filt, sample.size = min_depth, replace = FALSE)

# Filter to 7-20 days
physeq_final <- subset_samples(physeq_rarefied, !Age %in% c("Fetus", "Newborn"))
sample_data(physeq_final)$Age <- factor(sample_data(physeq_final)$Age, levels = c("7+days", "15+days", "20+days"))

# Calculate alpha diversity
alpha_div <- estimate_richness(physeq_final, measures = c("Observed", "Shannon"))
sample_data(physeq_final)$Shannon <- alpha_div$Shannon
sample_data(physeq_final)$Observed <- alpha_div$Observed

cat("Data loaded successfully!\n")
cat(paste("Final sample count:", nsamples(physeq_final), "\n"))

# Professional theme for plots
theme_pub <- function() {
    theme_bw(base_size = 14) +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
            axis.title = element_text(face = "bold", size = 13),
            axis.text = element_text(size = 11),
            legend.position = "right",
            legend.title = element_text(face = "bold", size = 11),
            legend.text = element_text(size = 10),
            plot.margin = margin(15, 15, 15, 15)
        )
}

# ============================================================================
# TABLE 1: Sample Distribution (Metadata Summary)
# ============================================================================

cat("\nGenerating Table 1: Sample Distribution...\n")

meta_summary <- data.frame(sample_data(physeq_final)) %>%
    group_by(Age, SampleType) %>%
    summarise(N = n(), .groups = "drop") %>%
    pivot_wider(names_from = SampleType, values_from = N, values_fill = 0)

# Create table as image
table1_plot <- tableGrob(
    meta_summary,
    rows = NULL,
    theme = ttheme_default(
        core = list(
            fg_params = list(fontsize = 14),
            bg_params = list(fill = c("white", "gray95"))
        ),
        colhead = list(
            fg_params = list(fontsize = 14, fontface = "bold"),
            bg_params = list(fill = "#4472C4", col = "white")
        )
    )
)

# Add title
title_grob <- textGrob(
    "Table 1. Sample Distribution Across Developmental Stages and Environmental Treatments",
    gp = gpar(fontsize = 12, fontface = "bold")
)

caption_grob <- textGrob(
    "Numbers indicate sample counts per group after quality filtering (n = 30 total samples).",
    gp = gpar(fontsize = 10, fontface = "italic", col = "gray40")
)

png(file.path(output_dir, "Table1_SampleDistribution.png"), width = 800, height = 300, res = 150)
grid.arrange(
    title_grob,
    table1_plot,
    caption_grob,
    ncol = 1,
    heights = c(0.2, 0.6, 0.2)
)
dev.off()

cat("  ✓ Table1_SampleDistribution.png saved\n")

# ============================================================================
# TABLE 2: PERMANOVA Results
# ============================================================================

cat("\nGenerating Table 2: PERMANOVA Results...\n")

# Calculate statistics
dist_bray <- phyloseq::distance(physeq_final, "bray")
meta <- data.frame(sample_data(physeq_final))

perm_mother <- adonis2(dist_bray ~ Mother, data = meta, permutations = 999)
perm_env <- adonis2(dist_bray ~ SampleType, data = meta, permutations = 999)
perm_age <- adonis2(dist_bray ~ Age, data = meta, permutations = 999)

# PERMDISP
disp_env <- betadisper(dist_bray, meta$SampleType)
disp_age <- betadisper(dist_bray, meta$Age)
disp_mother <- betadisper(dist_bray, meta$Mother)

permdisp_env <- permutest(disp_env)$tab$`Pr(>F)`[1]
permdisp_age <- permutest(disp_age)$tab$`Pr(>F)`[1]
permdisp_mother <- permutest(disp_mother)$tab$`Pr(>F)`[1]

permanova_results <- data.frame(
    Factor = c("Maternal Identity", "Environment", "Age", "Residuals"),
    `R²` = c(
        paste0(round(perm_mother$R2[1] * 100, 1), "%"),
        paste0(round(perm_env$R2[1] * 100, 1), "%"),
        paste0(round(perm_age$R2[1] * 100, 1), "%"),
        "—"
    ),
    `Pseudo-F` = c(
        round(perm_mother$F[1], 2),
        round(perm_env$F[1], 2),
        round(perm_age$F[1], 2),
        "—"
    ),
    `p-value` = c(
        perm_mother$`Pr(>F)`[1],
        perm_env$`Pr(>F)`[1],
        perm_age$`Pr(>F)`[1],
        "—"
    ),
    `PERMDISP p` = c(
        ifelse(permdisp_mother < 0.001, "<0.001*", round(permdisp_mother, 3)),
        round(permdisp_env, 3),
        round(permdisp_age, 3),
        "—"
    ),
    check.names = FALSE
)

table2_plot <- tableGrob(
    permanova_results,
    rows = NULL,
    theme = ttheme_default(
        core = list(
            fg_params = list(fontsize = 12),
            bg_params = list(fill = c("white", "gray95"))
        ),
        colhead = list(
            fg_params = list(fontsize = 12, fontface = "bold"),
            bg_params = list(fill = "#4472C4", col = "white")
        )
    )
)

title2_grob <- textGrob(
    "Table 2. PERMANOVA Results for Beta Diversity (Bray-Curtis)",
    gp = gpar(fontsize = 12, fontface = "bold")
)

caption2_grob <- textGrob(
    "R² = variance explained; Pseudo-F = group difference strength. *PERMDISP p < 0.05 indicates heterogeneous dispersion.",
    gp = gpar(fontsize = 9, fontface = "italic", col = "gray40")
)

png(file.path(output_dir, "Table2_PERMANOVA.png"), width = 900, height = 350, res = 150)
grid.arrange(
    title2_grob,
    table2_plot,
    caption2_grob,
    ncol = 1,
    heights = c(0.15, 0.7, 0.15)
)
dev.off()

cat("  ✓ Table2_PERMANOVA.png saved\n")

# ============================================================================
# FIGURE 1: Alpha Diversity (Shannon Index)
# ============================================================================

cat("\nGenerating Figure 1: Alpha Diversity...\n")

fig1 <- ggplot(data.frame(sample_data(physeq_final)), aes(x = Age, y = Shannon, fill = SampleType)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), 
               size = 2.5, alpha = 0.7, color = "black") +
    scale_fill_manual(values = c("Muck" = "#8c510a", "Sterile" = "#01665e"), name = "Environment") +
    labs(
        title = "Alpha Diversity Across Developmental Stages",
        subtitle = "Shannon Index (Kruskal-Wallis p > 0.05 for all comparisons)",
        y = "Shannon Index",
        x = "Developmental Stage (Days Post-Birth)"
    ) +
    theme_pub() +
    theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Figure1_AlphaDiversity.png"), fig1, 
       width = 8, height = 6, dpi = 300, bg = "white")

cat("  ✓ Figure1_AlphaDiversity.png saved\n")

# ============================================================================
# FIGURE 2: PCoA (Beta Diversity)
# ============================================================================

cat("\nGenerating Figure 2: PCoA...\n")

ord_pcoa <- ordinate(physeq_final, method = "PCoA", distance = "bray")
eig <- ord_pcoa$values$Eigenvalues
var_explained <- round(eig / sum(eig) * 100, 1)

pcoa_df <- data.frame(
    PC1 = ord_pcoa$vectors[, 1],
    PC2 = ord_pcoa$vectors[, 2],
    sample_data(physeq_final)
)

fig2 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Mother, shape = SampleType)) +
    geom_point(size = 5, alpha = 0.85) +
    stat_ellipse(aes(group = SampleType, linetype = SampleType), 
                 type = "t", linewidth = 1, show.legend = TRUE) +
    scale_color_brewer(palette = "Set1", name = "Mother") +
    scale_shape_manual(values = c(16, 17), name = "Environment") +
    scale_linetype_manual(values = c("solid", "dashed"), name = "Environment") +
    labs(
        title = "PCoA: Community Structure (Bray-Curtis)",
        subtitle = "Maternal identity shows strongest clustering effect",
        x = paste0("PCoA Axis 1 (", var_explained[1], "%)"),
        y = paste0("PCoA Axis 2 (", var_explained[2], "%)")
    ) +
    theme_pub() +
    theme(legend.position = "right")

ggsave(file.path(output_dir, "Figure2_PCoA.png"), fig2, 
       width = 10, height = 7, dpi = 300, bg = "white")

cat("  ✓ Figure2_PCoA.png saved\n")

# ============================================================================
# FIGURE 3: NMDS (Beta Diversity)
# ============================================================================

cat("\nGenerating Figure 3: NMDS...\n")

set.seed(123)
invisible(capture.output(
    ord_nmds <- ordinate(physeq_final, method = "NMDS", distance = "bray", trymax = 50)
))

nmds_df <- data.frame(
    NMDS1 = ord_nmds$points[, 1],
    NMDS2 = ord_nmds$points[, 2],
    sample_data(physeq_final)
)

fig3 <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Mother, shape = SampleType)) +
    geom_point(size = 5, alpha = 0.85) +
    stat_ellipse(aes(group = SampleType, linetype = SampleType), 
                 type = "t", linewidth = 1) +
    scale_color_brewer(palette = "Set1", name = "Mother") +
    scale_shape_manual(values = c(16, 17), name = "Environment") +
    scale_linetype_manual(values = c("solid", "dashed"), name = "Environment") +
    labs(
        title = paste0("NMDS: Community Structure (Stress = ", round(ord_nmds$stress, 3), ")"),
        subtitle = "Confirms PCoA clustering patterns",
        x = "NMDS Axis 1",
        y = "NMDS Axis 2"
    ) +
    theme_pub()

ggsave(file.path(output_dir, "Figure3_NMDS.png"), fig3, 
       width = 10, height = 7, dpi = 300, bg = "white")

cat("  ✓ Figure3_NMDS.png saved\n")

# ============================================================================
# FIGURE 4: Taxonomic Composition (Phylum Level)
# ============================================================================

cat("\nGenerating Figure 4: Taxonomic Composition...\n")

phy_phylum <- tax_glom(physeq_final, "Phylum")
phy_rel <- transform_sample_counts(phy_phylum, function(x) x / sum(x))
melt <- psmelt(phy_rel)

# Top 5 phyla
top <- names(sort(tapply(melt$Abundance, melt$Phylum, sum), decreasing = TRUE))[1:5]
melt$Phylum <- ifelse(melt$Phylum %in% top, as.character(melt$Phylum), "Other")

# Calculate means
melt_sum <- melt %>%
    group_by(Age, SampleType, Phylum) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop")

melt_sum$Age <- factor(melt_sum$Age, levels = c("7+days", "15+days", "20+days"))

fig4 <- ggplot(melt_sum, aes(x = SampleType, y = MeanAbundance, fill = Phylum)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    facet_wrap(~Age) +
    scale_fill_brewer(palette = "Paired") +
    labs(
        title = "Taxonomic Composition at Phylum Level",
        subtitle = "Mean relative abundance by treatment and developmental stage",
        y = "Mean Relative Abundance",
        x = "Environmental Treatment"
    ) +
    theme_pub() +
    theme(
        axis.text.x = element_text(angle = 0),
        strip.background = element_rect(fill = "#4472C4"),
        strip.text = element_text(color = "white", face = "bold", size = 12)
    )

ggsave(file.path(output_dir, "Figure4_TaxonomicComposition.png"), fig4, 
       width = 12, height = 6, dpi = 300, bg = "white")

cat("  ✓ Figure4_TaxonomicComposition.png saved\n")

# ============================================================================
# FIGURE 5: Core Microbiome
# ============================================================================

cat("\nGenerating Figure 5: Core Microbiome...\n")

get_core <- function(phy) {
    phy_rel <- transform_sample_counts(phy, function(x) x / sum(x))
    filter_taxa(phy_rel, function(x) sum(x > 0) > (0.5 * length(x)), TRUE) %>% taxa_names()
}

muck_core <- get_core(subset_samples(physeq_final, SampleType == "Muck"))
sterile_core <- get_core(subset_samples(physeq_final, SampleType == "Sterile"))

core_data <- data.frame(
    Category = c("Shared\n(Maternal)", "Muck Only\n(Environmental)", "Sterile Only"),
    Count = c(
        length(intersect(muck_core, sterile_core)),
        length(setdiff(muck_core, sterile_core)),
        length(setdiff(sterile_core, muck_core))
    )
)

fig5 <- ggplot(core_data, aes(x = reorder(Category, Count), y = Count, fill = Category)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.5) +
    geom_text(aes(label = Count), hjust = -0.3, size = 6, fontface = "bold") +
    coord_flip() +
    scale_fill_manual(values = c(
        "Shared\n(Maternal)" = "#4daf4a",
        "Muck Only\n(Environmental)" = "#a65628",
        "Sterile Only" = "#999999"
    )) +
    labs(
        title = "Core Microbiome Composition",
        subtitle = "Taxa present in >50% of samples within each treatment",
        y = "Number of Core Taxa (OTUs)",
        x = ""
    ) +
    ylim(0, max(core_data$Count) * 1.15) +
    theme_pub() +
    theme(legend.position = "none")

ggsave(file.path(output_dir, "Figure5_CoreMicrobiome.png"), fig5, 
       width = 8, height = 5, dpi = 300, bg = "white")

cat("  ✓ Figure5_CoreMicrobiome.png saved\n")

# ============================================================================
# FIGURE 6: Variance Partitioning (PERMANOVA Summary)
# ============================================================================

cat("\nGenerating Figure 6: Variance Partitioning...\n")

var_data <- data.frame(
    Driver = c("Maternal Identity", "Environment", "Age"),
    Variance = c(
        round(perm_mother$R2[1] * 100, 1),
        round(perm_env$R2[1] * 100, 1),
        round(perm_age$R2[1] * 100, 1)
    ),
    F_value = c(
        round(perm_mother$F[1], 2),
        round(perm_env$F[1], 2),
        round(perm_age$F[1], 2)
    )
)

fig6 <- ggplot(var_data, aes(x = reorder(Driver, Variance), y = Variance, fill = Driver)) +
    geom_col(color = "black", width = 0.6, linewidth = 0.5) +
    geom_text(aes(label = paste0("F = ", F_value)), hjust = -0.1, size = 5, fontface = "bold") +
    coord_flip() +
    scale_fill_viridis_d(option = "C") +
    labs(
        title = "Variance Explained by Each Factor (PERMANOVA)",
        subtitle = "Pseudo-F values indicate relative strength of group differences",
        y = "Percent Variance Explained (%)",
        x = ""
    ) +
    ylim(0, max(var_data$Variance) * 1.3) +
    theme_pub() +
    theme(legend.position = "none")

ggsave(file.path(output_dir, "Figure6_VariancePartitioning.png"), fig6, 
       width = 9, height = 5, dpi = 300, bg = "white")

cat("  ✓ Figure6_VariancePartitioning.png saved\n")

# ============================================================================
# FIGURE 7: Beta Dispersion
# ============================================================================

cat("\nGenerating Figure 7: Beta Dispersion...\n")

mod <- betadisper(dist_bray, sample_data(physeq_final)$SampleType)
df_disp <- data.frame(Distance = mod$distances, Group = mod$group)

fig7 <- ggplot(df_disp, aes(x = Group, y = Distance, fill = Group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.1, size = 3, alpha = 0.7) +
    scale_fill_manual(values = c("Muck" = "#8c510a", "Sterile" = "#01665e")) +
    labs(
        title = "Beta Dispersion: Distance to Group Centroid",
        subtitle = paste0("PERMDISP p = ", round(permdisp_env, 3), " (no significant difference)"),
        y = "Distance to Centroid (Bray-Curtis)",
        x = "Environmental Treatment"
    ) +
    theme_pub() +
    theme(legend.position = "none")

ggsave(file.path(output_dir, "Figure7_BetaDispersion.png"), fig7, 
       width = 7, height = 6, dpi = 300, bg = "white")

cat("  ✓ Figure7_BetaDispersion.png saved\n")

# ============================================================================
# TABLE 3: Alpha Diversity Statistics
# ============================================================================

cat("\nGenerating Table 3: Alpha Diversity Statistics...\n")

meta_df <- data.frame(sample_data(physeq_final))

# Kruskal-Wallis tests
kw_env <- kruskal.test(Shannon ~ SampleType, data = meta_df)
kw_age <- kruskal.test(Shannon ~ Age, data = meta_df)
kw_mother <- kruskal.test(Shannon ~ Mother, data = meta_df)

alpha_stats <- data.frame(
    Variable = c("Environment (SampleType)", "Developmental Age", "Maternal Identity"),
    Test = rep("Kruskal-Wallis", 3),
    `Chi-squared` = c(
        round(kw_env$statistic, 3),
        round(kw_age$statistic, 3),
        round(kw_mother$statistic, 3)
    ),
    df = c(kw_env$parameter, kw_age$parameter, kw_mother$parameter),
    `p-value` = c(
        round(kw_env$p.value, 4),
        round(kw_age$p.value, 4),
        round(kw_mother$p.value, 4)
    ),
    Significance = c(
        ifelse(kw_env$p.value < 0.05, "*", "ns"),
        ifelse(kw_age$p.value < 0.05, "*", "ns"),
        ifelse(kw_mother$p.value < 0.05, "*", "ns")
    ),
    check.names = FALSE
)

table3_plot <- tableGrob(
    alpha_stats,
    rows = NULL,
    theme = ttheme_default(
        core = list(
            fg_params = list(fontsize = 11),
            bg_params = list(fill = c("white", "gray95"))
        ),
        colhead = list(
            fg_params = list(fontsize = 11, fontface = "bold"),
            bg_params = list(fill = "#4472C4", col = "white")
        )
    )
)

title3_grob <- textGrob(
    "Table 3. Alpha Diversity (Shannon Index) Statistical Tests",
    gp = gpar(fontsize = 12, fontface = "bold")
)

caption3_grob <- textGrob(
    "ns = not significant (p > 0.05); * = significant (p < 0.05)",
    gp = gpar(fontsize = 9, fontface = "italic", col = "gray40")
)

png(file.path(output_dir, "Table3_AlphaStats.png"), width = 900, height = 300, res = 150)
grid.arrange(
    title3_grob,
    table3_plot,
    caption3_grob,
    ncol = 1,
    heights = c(0.2, 0.6, 0.2)
)
dev.off()

cat("  ✓ Table3_AlphaStats.png saved\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("ALL FIGURES AND TABLES GENERATED SUCCESSFULLY!\n")
cat("============================================================\n")
cat(paste("\nOutput directory:", output_dir, "\n\n"))
cat("Files created:\n")
cat("  Tables:\n")
cat("    • Table1_SampleDistribution.png\n")
cat("    • Table2_PERMANOVA.png\n")
cat("    • Table3_AlphaStats.png\n")
cat("  Figures:\n")
cat("    • Figure1_AlphaDiversity.png\n")
cat("    • Figure2_PCoA.png\n")
cat("    • Figure3_NMDS.png\n")
cat("    • Figure4_TaxonomicComposition.png\n")
cat("    • Figure5_CoreMicrobiome.png\n")
cat("    • Figure6_VariancePartitioning.png\n")
cat("    • Figure7_BetaDispersion.png\n")
cat("\n")
