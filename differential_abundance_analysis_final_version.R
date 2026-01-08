library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)
library(lmerTest)

###############################################################################
folder_path <- "~/differential_abundance_analysis/"
setwd(folder_path)
###############################################################################

# 1. Read the species relative abundance data
species_data <- read.table(file = "metaphlan_results_new.tsv", sep = "\t", header = T, row.names = 1)
species <- species_data[grepl("s__", rownames(species_data)) & !grepl("t__", rownames(species_data)), ]

# ## Get taxonomy information:
# taxonomy <- data.frame("taxonomy" = rownames(species))
# 
# # Splitting into multiple columns
# taxonomy <- taxonomy %>% separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = "\\|", fill = "right")
# taxonomy_done <- as.data.frame(apply(taxonomy, 2, function(x) gsub("^[a-z]__", "", x)))
# write.csv(taxonomy_done, "/Users/M306307/Desktop/LBD/analysis/taxonomy_info.csv", row.names = F)

taxonomy <- read.csv(file = "taxonomy_info.csv")

# 2. Modify taxonomy names
new_rownames <- gsub(".*s__", "", rownames(species))
rownames(species) <- new_rownames

# 3. Modify sample names
col_names <- names(species)
new_colnames <- gsub("metaphlan_|_S.*", "", col_names)
colnames(species) <- new_colnames

# 4. Read the metadata
metadata <- read.csv(file = "metadata.csv", sep = ",", header = T, row.names = 1)
species_clean <- species[, c(rownames(metadata))]

# 5. Change the taxonomic data from percentage into proportion
colSums(species_clean)
species_prop <- data.frame(apply(species_clean, 2, function(x) x/sum(x)))
colSums(species_prop)

# 6. Find the appropriate presence cutoff for bacterial species data
relab_number <- unlist(species_prop, use.name = FALSE)
relab_number_ordered <- relab_number[order(relab_number, decreasing = T)]

pdf("rank_plot_with_threshold_microbial_species_done.pdf", height = 6, width = 8)
plot(x=c(1:length(relab_number_ordered)), y=log10(relab_number_ordered), pch = 20, xlab = "Rank", ylab = "Relative abundance of microbial species", main = "Rank plot of ordered relative abundance of microbial species", xlim = c(0, 14000))
abline(h=-4.7, col = "orange", lwd = 3)
dev.off()

## It seems 10^-4.5 is the good cutoff

## Number of low abundance species in each sample
cut_off <- 10^-4.7
species_cutoff <- species_prop
species_cutoff[species_cutoff <= cut_off] <- 0

# 7. Add metadata to the proportion data
species_t <- as.data.frame(t(species_cutoff))
species_all <- merge(species_t, metadata, by = "row.names", all = F)
rownames(species_all) <- species_all$Row.names
species_all <- species_all[, -1]
write.csv(species_all, "species_all.csv")

###############################################################################

## The following function performs a transformation on a dataset and merges it with metadata and group information to produce a final structured data frame.

trs <- function(dataset) {
  
  transformed <- asin(sqrt(dataset[, 1:(ncol(dataset)-33)]))

  return(transformed)
  
} 

# Find the bacterial species that are passing the prevalence cut-off (10%)

prevalence_cutoff <- function(transformed_data, metadata) {
  
  filtered_data <- transformed_data[, colSums(transformed_data > 0) >= 0.1*nrow(transformed_data)]
  result_df <- merge(filtered_data, metadata, by = "row.names", all = F)
  rownames(result_df) <- result_df$Row.names
  result_df <- result_df[, -1]
  
  return(result_df)
  
}

# The following function performs a mixed-effects model analysis to compare the abundance of species between two time points. It calculates statistical significance (p-value) and trend direction (difference) for each species across the given time points.

mixed_effect <- function(filtered_data, file_name) {
  
  p_diff_relab <- c()
  
  for (i in 1:(ncol(filtered_data)-33)) {
    p_diff_relab[i] <- tryCatch({summary(lmer(filtered_data[, i] ~ condition + age + BMI + (1|household_id), data = filtered_data, REML = F))[["coefficients"]][2,5]}, error = function(e) NA) # If error occurs, assign NA and continue
  }
  
  results <- data.frame("species" = names(filtered_data)[1:(ncol(filtered_data)-33)], "p_value" = p_diff_relab)
  results$q_value <- p.adjust(results$p_value, method = "BH")
  results_df <- merge(results, taxonomy, by = "species", all = F)
  
  write.csv(results_df, file = file_name, row.names = F)
  return(results_df)

}

log2_fc_results <- function(filtered_data, group1, group2, p_val_results, file_name) {
  
  # Re-name the column names so it contains the group labels
  fc_data <- filtered_data[, c(1:(ncol(filtered_data)-33), ncol(filtered_data)-30)]
  fc_data$new_id <- paste(fc_data$condition, rownames(fc_data), sep = "_")
  rownames(fc_data) <- fc_data$new_id
  fc_data <- fc_data[, 1:(ncol(filtered_data)-33)]
  fc_data_t <- as.data.frame(t(fc_data))
  
  # Compute medians for each group
  fc_data_t$median_g1 <- apply(fc_data_t[, grepl(group1, colnames(fc_data_t))], 1, function(x) median(x))
  fc_data_t$median_g2 <- apply(fc_data_t[, grepl(group2, colnames(fc_data_t))], 1, function(x) median(x))
  
  # Compute mean for each group
  fc_data_t$mean_g1 <- apply(fc_data_t[, grepl(group1, colnames(fc_data_t))], 1, function(x) mean(x))
  fc_data_t$mean_g2 <- apply(fc_data_t[, grepl(group2, colnames(fc_data_t))], 1, function(x) mean(x))
  
  # Compute log2 median fold change
  pseudocount <- cut_off
  fc_data_t$log2fc <- log2((fc_data_t$median_g1 + pseudocount)/(fc_data_t$median_g2 + pseudocount))
  
  # Compute log2 mean fold change
  pseudocount <- cut_off
  fc_data_t$log2fc_mean <- log2((fc_data_t$mean_g1 + pseudocount)/(fc_data_t$mean_g2 + pseudocount))
  
  # Compute mean difference
  fc_data_t$mean_difference <- fc_data_t$mean_g1 - fc_data_t$mean_g2
  
  # Ensure species names from p_val and log2fc results are the same
  if (!all(rownames(fc_data_t) %in% p_val_results$species)) {
    stop("Error: Some species in `fc_data_t` are missing in `p_results`.")
  }
  
  # Ensure there is a column named "p_value" in the p_val_results
  if (!"p_value" %in% colnames(p_val_results)) {
    stop("Error: `p_results` must contain a `p_value` column.")
  }
  
  # Merge p-values with log2fc results
  fc_data_t_p <- merge(fc_data_t, p_val_results, by.x = "row.names", by.y = "species", all = F)
  fc_data_t_p_clean <- fc_data_t_p[, c(1, (ncol(fc_data_t_p)-14):ncol(fc_data_t_p))]
  colnames(fc_data_t_p_clean)[1] <- "species"
  
  write.csv(fc_data_t_p_clean, file = file_name, row.names = F)
  return(list(fc_data_t_p_clean, fc_data_t_p))
  
}

## The following function processes the filtered dataset containing species abundance values for two conditions, computes the log2 fold-change, integrates p-values from mixed-effects linear model, and categorizes significant differences, and generates a volcano plot. 

volcano_plot <- function(fc_data_t_p, label1, label2, xlim1, xlim2) {
  
  # Define categories
  g1 <- paste0(label1, ">", label2)
  g2 <- paste0(label1, "<", label2)
  
  fc_data_t_p$category <- "All"
  fc_data_t_p$category[fc_data_t_p$log2fc > 0 & fc_data_t_p$p_value < 0.05] <- g1
  fc_data_t_p$category[fc_data_t_p$log2fc < 0 & fc_data_t_p$p_value < 0.05] <- g2
  fc_data_t_p$category <- factor(fc_data_t_p$category, levels = c("All", g1, g2))
  
  # Extract significant data and make them in an increasing order based on p_value
  positive_fold <- fc_data_t_p %>% filter(category == g1) %>% arrange(p_value)
  negative_fold <- fc_data_t_p %>% filter(category == g2) %>% arrange(p_value)
  
  write.csv(positive_fold, file = paste0("species_higher_in_", label1, "_than_", label2, ".csv"), row.names = F)
  write.csv(negative_fold, file = paste0("species_higher_in_", label2, "_than_", label1, ".csv"), row.names = F)
  
  plot <- ggplot(data = fc_data_t_p, aes(x = log2fc, y = -log10(p_value), color = category)) + 
    geom_point(size = 1, color = "lightgrey") + 
    xlim(as.numeric(xlim1), as.numeric(xlim2)) + 
    geom_hline(yintercept= -log10(0.05), color = "black", linetype = "dashed") +
    geom_vline(xintercept = log2(1), color = "black",linetype = "dashed") +
    
    ## Add points that have significant q-values in positive fold change
    geom_point(data = positive_fold, aes(x = log2fc, y = -log10(p_value)), size = 1.5) +
    geom_text_repel(data = positive_fold[1:(min(nrow(positive_fold), 5)), ], 
                    aes(x = log2fc, y = -log10(p_value), label= Row.names), 
                    max.overlaps = 15) +
    
    ## Add points that have significant p-values in negative fold change
    geom_point(data = negative_fold, aes(x = log2fc, y = -log10(p_value)), size = 1.5) +
    geom_text_repel(data = negative_fold[1:(min(nrow(negative_fold), 5)), ], 
                    aes(x = log2fc, y = -log10(p_value), label = Row.names), 
                    max.overlaps = 15) +
    
    labs(x = paste0("log2fold-change of relative abundance (", label1, "/", label2, ")"), 
         y = "-log10(p-value)") + ggtitle("") + 
    
    # Define color mapping
    scale_color_manual(values = setNames(c("brown2", "#45B1E9"), c(g1, g2)), 
                       labels = c(g1, g2)) +
    theme_bw()
  
  pdf(file = paste0("volcano_plot_", label1, "_vs_", label2, ".pdf"), height = 7, width = 8)
  print(plot)
  dev.off()
  
  return(list(positive_fold, negative_fold, fc_data_t_p))
  
}

###############################################################################

# LBD vs. Control#

# 1. Extract data for LBD and their controls
lbd_vs_control_relab <- species_all[species_all$condition %in% c("lbd", "lbd_control"), ]

# 2. Use arcsine square root transformation for the data
lbd_vs_control_relab_transformed <- trs(lbd_vs_control_relab)

# 3. Find the bacterial species that are passing the prevalence cut-off
lbd_vs_control_relab_transformed_filtered <- prevalence_cutoff(lbd_vs_control_relab_transformed, metadata)

# 4. Use the mixed-effect linear model to get p values condition variable. household_id will be random effect, and age and BMI are confounding variables. 
lbd_vs_control_relab_results <- mixed_effect(lbd_vs_control_relab_transformed_filtered, "lbd_vs_control_differential_abundance_p_values.csv")
lbd_vs_control_relab_results_adj <- read.csv(file = "lbd_vs_control_differential_abundance_p_values.csv")

sum(lbd_vs_control_relab_results$p_value < 0.01)
sum(lbd_vs_control_relab_results$p_value < 0.05)

# 5. Use median of each group to calculate the log2fc
lbd_vs_control_diff_abun_median <- log2_fc_results(lbd_vs_control_relab_transformed_filtered, "lbd_BIOME", "lbd_control", lbd_vs_control_relab_results, "lbd_vs_control_diff_abun_p_fc_done.csv")

# 6. Only focus on the differentially abundant species that have median relative abundance not 0, and p < 0.05
sig_lbd_vs_control_diff_abun_species <- lbd_vs_control_diff_abun_median[[1]][lbd_vs_control_diff_abun_median[[1]]$p_value < 0.05 & lbd_vs_control_diff_abun_median[[1]]$log2fc != 0, ]
nrow(sig_lbd_vs_control_diff_abun_species)
sig_lbd_vs_control_diff_abun_species$species
sig_lbd_vs_control_diff_abun_species$group <- ifelse(sig_lbd_vs_control_diff_abun_species$log2fc>0, "LBD>Control", "LBD<Control")
# write.csv(sig_lbd_vs_control_diff_abun_species, file = "sig_lbd_vs_control_diff_abun_species.csv", row.names = F)

# 7. Create the volcano plot
lbd_vs_control_fc_results <- volcano_plot(lbd_vs_control_diff_abun_median[[2]], "LBD", "Control", -12, 12) 

# 8. Create the overall boxplot
positive_fold <- lbd_vs_control_fc_results[[1]][1:4, ]
negative_fold <- lbd_vs_control_fc_results[[2]][1:5, ]
order_of_species_pos <- positive_fold[order(positive_fold$median_g1, decreasing = T), ]$Row.names
order_of_species_neg <- negative_fold[order(negative_fold$median_g2, decreasing = T), ]$Row.names
order_of_species <- c(order_of_species_pos, order_of_species_neg)
boxplot_data <- rbind(lbd_vs_control_fc_results[[1]], lbd_vs_control_fc_results[[2]])
boxplot_data_done <- boxplot_data[boxplot_data$Row.names %in% order_of_species, ]

boxplot_long_df <- melt(boxplot_data_done[, 1:51], id.vars = "Row.names")
boxplot_long_df$condition <- ifelse(grepl("control", boxplot_long_df$variable), "Control", "LBD")
colnames(boxplot_long_df)[1] <- "species"
boxplot_long_df$species <- factor(boxplot_long_df$species, levels = order_of_species)
boxplot_long_df$condition <- factor(boxplot_long_df$condition, levels = c("LBD", "Control"))

pdf(file = paste0("lbd_vs_control_boxplot_species_v2.pdf"), height = 8, width = 10)
plot <- ggplot(boxplot_long_df, 
                     aes(x = species, y = value, fill = condition)) + 
  geom_boxplot(outlier.color = "black", 
               outlier.fill = "white", 
               outlier.size = 2, 
               outlier.shape = 21) +
  scale_fill_manual(values = c("#ff9274", "#55b7e6")) +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Relative abundance (arcsine square root transformed)")
print(plot)
dev.off()

# 8.1 create boxplots for three species
order_of_species <- c("Roseburia_hominis", "Eubacterium_ramulus", "Intestinimonas_butyriciproducens")

boxplot_data <- rbind(lbd_vs_control_fc_results[[1]], lbd_vs_control_fc_results[[2]])
boxplot_data_done <- boxplot_data[boxplot_data$Row.names %in% order_of_species, ]

boxplot_long_df <- melt(boxplot_data_done[, 1:51], id.vars = "Row.names")
boxplot_long_df$condition <- ifelse(grepl("control", boxplot_long_df$variable), "Control", "LBD")
colnames(boxplot_long_df)[1] <- "species"
boxplot_long_df$species <- factor(boxplot_long_df$species, levels = order_of_species)
boxplot_long_df$condition <- factor(boxplot_long_df$condition, levels = c("LBD", "Control"))

pdf(file = paste0("lbd_vs_control_boxplot_3species.pdf"), height = 5, width = 3.5)
plot <- ggplot(boxplot_long_df, 
               aes(x = species, y = value, fill = condition)) + 
  geom_boxplot(outlier.color = "black", 
               outlier.fill = "white", 
               outlier.size = 2, 
               outlier.shape = 21) +
  scale_fill_manual(values = c("#ff9274", "#55b7e6")) +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Relative abundance (arcsine square root transformed)")
print(plot)
dev.off()

# 9. Calculate the difference between mean based on the relative abundance data before arcsine square-root transformed)
sig_lbd_vs_control_species_relab_og <- lbd_vs_control_relab[, colnames(lbd_vs_control_relab) %in% c(boxplot_data_done$Row.names, "condition")]

# Re-name the column names so it contains the group labels
sig_lbd_vs_control_species_relab_og$new_id <- paste(sig_lbd_vs_control_species_relab_og$condition, rownames(sig_lbd_vs_control_species_relab_og), sep = "_")
rownames(sig_lbd_vs_control_species_relab_og) <- sig_lbd_vs_control_species_relab_og$new_id
sig_lbd_vs_control_species_relab_og <- sig_lbd_vs_control_species_relab_og[, 1:(length(boxplot_data_done$Row.names))]
sig_lbd_vs_control_species_relab_og_t <- as.data.frame(t(sig_lbd_vs_control_species_relab_og))

# Compute mean for each group
sig_lbd_vs_control_species_relab_og_t$mean_lbd <- apply(sig_lbd_vs_control_species_relab_og_t[, grepl("lbd_BIOME", colnames(sig_lbd_vs_control_species_relab_og_t))], 1, function(x) mean(x))
sig_lbd_vs_control_species_relab_og_t$mean_control <- apply(sig_lbd_vs_control_species_relab_og_t[, grepl("lbd_control", colnames(sig_lbd_vs_control_species_relab_og_t))], 1, function(x) mean(x))

sig_lbd_vs_control_species_relab_og_t$mean_diff <- sig_lbd_vs_control_species_relab_og_t$mean_lbd - sig_lbd_vs_control_species_relab_og_t$mean_control
sig_lbd_vs_control_species_relab_og_t$species <- rownames(sig_lbd_vs_control_species_relab_og_t)

sig_lbd_vs_control_species_relab_og_t$species <- factor(sig_lbd_vs_control_species_relab_og_t$species, levels = order_of_species)

pdf(file = "mean_diff_bar_lbd_control_v2.pdf", height = 3.5, width = 8)
ggplot(sig_lbd_vs_control_species_relab_og_t, aes(x = as.factor(species), y = mean_diff)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "mean difference")
dev.off()

# 10. Create boxplot for each significant species
for (i in 1:length(order_of_species)) {
  plot <- ggplot(boxplot_long_df[boxplot_long_df$species == order_of_species[i], ], 
                 aes(x = condition, y = value, fill = condition)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(shape = 21, size = 2, position = position_jitter(width = 0.1)) +
    scale_fill_manual(values = c("#ff9274", "#55b7e6")) +
    theme_bw() + 
    theme(
      panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()) +
    labs(x = "", 
         y = "Microbial relative abundance (Arcsine square-root transformed)",
         title = order_of_species[i])
  
  pdf(file = paste0("individual_plots/lbd_vs_control_boxplot_", order_of_species[i], ".pdf"), height = 6, width = 4)
  print(plot)
  dev.off()
}

###############################################################################

# iRBD vs. Control#

# 1. Extract data for iRBD and their controls
irbd_vs_control_relab <- species_all[species_all$condition %in% c("irbd", "irbd_control"), ]

# 2. Use arcsine square root transformation for the data
irbd_vs_control_relab_transformed <- trs(irbd_vs_control_relab, metadata)

# 3. Find the bacterial species that are passing the prevalence cut-off
irbd_vs_control_relab_transformed_filtered <- prevalence_cutoff(irbd_vs_control_relab_transformed, metadata)

# 4. Use the mixed-effect linear model to get p values condition variable. household_id will be random effect, and age and BMI are confounding variables. 
irbd_vs_control_relab_results <- mixed_effect(irbd_vs_control_relab_transformed_filtered, "irbd_vs_control_differential_abundance_p_values.csv")
irbd_vs_control_relab_results <- read.csv(file = "irbd_vs_control/irbd_vs_control_differential_abundance_p_values.csv")

sum(irbd_vs_control_relab_results$p_value < 0.01)
sum(irbd_vs_control_relab_results$p_value < 0.05)

# 5. Use median of each group to calculate the log2fc
irbd_vs_control_diff_abun_median <- log2_fc_results(irbd_vs_control_relab_transformed_filtered, "irbd_BIOME", "irbd_control", irbd_vs_control_relab_results, "irbd_vs_control_diff_abun_p_fc_done.csv")

# 6. Only focus on the differentially abundant species that have median relative abundance not 0, and p < 0.05
sig_irbd_vs_control_diff_abun_species <- irbd_vs_control_diff_abun_median[[1]][irbd_vs_control_diff_abun_median[[1]]$p_value < 0.05 & irbd_vs_control_diff_abun_median[[1]]$log2fc != 0, ]
nrow(sig_irbd_vs_control_diff_abun_species)
sig_irbd_vs_control_diff_abun_species$group <- ifelse(sig_irbd_vs_control_diff_abun_species$log2fc>0, "irbd>Control", "irbd<Control")
# write.csv(sig_irbd_vs_control_diff_abun_species, file = "irbd_vs_control/sig_irbd_vs_control_diff_abun_species.csv", row.names = F)

# 7. Create the volcano plot
irbd_vs_control_fc_results <- volcano_plot(irbd_vs_control_diff_abun_median[[2]], "iRBD", "Control", -12, 12) 

# 8. Create the overall boxplot
positive_fold <- irbd_vs_control_fc_results[[1]][1:5, ]
negative_fold <- irbd_vs_control_fc_results[[2]][1:4, ]
order_of_species_pos <- positive_fold[order(positive_fold$median_g1, decreasing = T), ]$Row.names
order_of_species_neg <- negative_fold[order(negative_fold$median_g2, decreasing = T), ]$Row.names
order_of_species <- c(order_of_species_pos, order_of_species_neg)
boxplot_data <- irbd_vs_control_diff_abun_median[[2]][, 1:21]

boxplot_long_df <- melt(boxplot_data, id.vars = "Row.names")
colnames(boxplot_long_df)[1] <- "species"
boxplot_long_df$condition <- ifelse(grepl("control", boxplot_long_df$variable), "Control", "iRBD")
boxplot_long_df_done <- boxplot_long_df[boxplot_long_df$species %in% order_of_species, ]
boxplot_long_df_done$species <- factor(boxplot_long_df_done$species, levels = order_of_species)
boxplot_long_df_done$condition <- factor(boxplot_long_df_done$condition, levels = c("iRBD", "Control"))

pdf(file = paste0("irbd_vs_control_boxplot.pdf"), height = 6, width = 10)
plot <- ggplot(boxplot_long_df_done, 
               aes(x = species, y = value, fill = condition)) + 
  geom_boxplot(outlier.color = "black", 
               outlier.fill = "white", 
               outlier.size = 2, 
               outlier.shape = 21) +
  # Add jitter if desired:
  # geom_point(aes(fill = condition), shape = 21, size = 1.5, 
  #            position = position_jitter(width = 0.2), color = "black") +
  scale_fill_manual(values = c("#2DA248", "#fdc848")) +
  # ylim(as.numeric(0), as.numeric(0.5)) + 
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Microbial relative abundance (Arcsine square root transformed)")
print(plot)
dev.off()

# 9. Calculate the difference between mean based on the relative abundance data before arcsine square-root transformed)
sig_irbd_vs_control_species_relab_og <- irbd_vs_control_relab[, colnames(irbd_vs_control_relab) %in% c(order_of_species, "condition")]

# Re-name the column names so it contains the group labels
sig_irbd_vs_control_species_relab_og$new_id <- paste(sig_irbd_vs_control_species_relab_og$condition, rownames(sig_irbd_vs_control_species_relab_og), sep = "_")
rownames(sig_irbd_vs_control_species_relab_og) <- sig_irbd_vs_control_species_relab_og$new_id
sig_irbd_vs_control_species_relab_og <- sig_irbd_vs_control_species_relab_og[, 1:(length(boxplot_data_done$Row.names))]
sig_irbd_vs_control_species_relab_og_t <- as.data.frame(t(sig_irbd_vs_control_species_relab_og))

# Compute mean for each group
sig_irbd_vs_control_species_relab_og_t$mean_irbd <- apply(sig_irbd_vs_control_species_relab_og_t[, grepl("irbd_BIOME", colnames(sig_irbd_vs_control_species_relab_og_t))], 1, function(x) mean(x))
sig_irbd_vs_control_species_relab_og_t$mean_control <- apply(sig_irbd_vs_control_species_relab_og_t[, grepl("irbd_control", colnames(sig_irbd_vs_control_species_relab_og_t))], 1, function(x) mean(x))

sig_irbd_vs_control_species_relab_og_t$mean_diff <- sig_irbd_vs_control_species_relab_og_t$mean_irbd - sig_irbd_vs_control_species_relab_og_t$mean_control
sig_irbd_vs_control_species_relab_og_t$species <- rownames(sig_irbd_vs_control_species_relab_og_t)

sig_irbd_vs_control_species_relab_og_t$species <- factor(sig_irbd_vs_control_species_relab_og_t$species, levels = order_of_species)

pdf(file = "irbd_vs_control/mean_diff_bar_irbd_control_v2.pdf", height = 3.5, width = 8)
ggplot(sig_irbd_vs_control_species_relab_og_t, aes(x = as.factor(species), y = mean_diff)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "mean difference")
dev.off()

# 9. Create boxplot for each significant species
for (i in 1:nrow(sig_irbd_vs_control_diff_abun_species)) {
  plot <- ggplot(boxplot_long_df_done[boxplot_long_df_done$species == sig_irbd_vs_control_diff_abun_species$species[i], ], 
                 aes(x = condition, y = value, fill = condition)) + 
    geom_boxplot(outlier.color = "black", 
                 outlier.fill = "white", 
                 outlier.size = 2, 
                 outlier.shape = 21) +
    # geom_point(shape = 21, size = 2, position = position_jitter(width = 0.1)) +
    scale_fill_manual(values = c("#2DA248", "#fdc848")) +
    theme_bw() + 
    theme(
      panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "", 
         y = "Relative abundance (arcsine square-root transformed)",
         title = sig_irbd_vs_control_diff_abun_species$species[i])
  
  pdf(file = paste0("irbd_vs_control/individual_plots/irbd_vs_control_boxplot_", sig_irbd_vs_control_diff_abun_species$species[i], "_v1.pdf"), height = 6, width = 4)
  print(plot)
  dev.off()
}

###############################################################################

# LBD vs. iRBD #

# 1. Extract data for LBD and their controls
lbd_vs_irbd_relab <- species_all[species_all$condition %in% c("lbd", "irbd"), ]

# 2. Use arcsine square root transformation for the data
lbd_vs_irbd_relab_transformed <- trs(lbd_vs_irbd_relab)

# 3. Find the bacterial species that are passing the prevalence cut-off
lbd_vs_irbd_relab_transformed_filtered <- prevalence_cutoff(lbd_vs_irbd_relab_transformed, metadata)

# 4. Use the regular linear model to get p values condition variable. age and BMI are confounding variables. 
p <- c()

for (i in 1:(ncol(lbd_vs_irbd_relab_transformed_filtered)-33)) {
  p[i] <- summary(lm(lbd_vs_irbd_relab_transformed_filtered[, i] ~ condition + age + BMI, data = lbd_vs_irbd_relab_transformed_filtered))[["coefficients"]][2,4]
}

p_results <- data.frame("species" = names(lbd_vs_irbd_relab_transformed_filtered)[1:(ncol(lbd_vs_irbd_relab_transformed_filtered)-33)], "p_value" = p)

p_results$q_value <- p.adjust(p_results$p_value, method = "BH")

lbd_vs_irbd_relab_results <- merge(p_results, taxonomy, by = "species", all = F)

# write.csv(lbd_vs_irbd_relab_results, file = "lbd_vs_irbd/lbd_vs_irbd_differential_abundance_p_values.csv", row.names = F)

sum(lbd_vs_irbd_relab_results$p_value < 0.01)
sum(lbd_vs_irbd_relab_results$p_value < 0.05)

# 5. Use median of each group to calculate the log2fc
lbd_vs_irbd_diff_abun_median <- log2_fc_results(lbd_vs_irbd_relab_transformed_filtered, "lbd_BIOME", "irbd_BIOME", lbd_vs_irbd_relab_results, "lbd_vs_irbd_diff_abun_p_fc_done.csv")

# 6. Only focus on the differentially abundant species that have median relative abundance not 0, and p < 0.05
sig_lbd_vs_irbd_diff_abun_species <- lbd_vs_irbd_diff_abun_median[[1]][lbd_vs_irbd_diff_abun_median[[1]]$p_value < 0.05 & lbd_vs_irbd_diff_abun_median[[1]]$log2fc != 0, ]
nrow(sig_lbd_vs_irbd_diff_abun_species)
sig_lbd_vs_irbd_diff_abun_species$group <- ifelse(sig_lbd_vs_irbd_diff_abun_species$log2fc>0, "LBD>iRBD", "LBD<iRBD")
# write.csv(sig_lbd_vs_irbd_diff_abun_species, file = "lbd_vs_irbd/sig_lbd_vs_irbd_diff_abun_species.csv", row.names = F)

# 7. Create the volcano plot
lbd_vs_irbd_fc_results <- volcano_plot(lbd_vs_irbd_diff_abun_median[[2]], "LBD", "iRBD", -13, 13) 

# 8. Create the overall boxplot
positive_fold <- lbd_vs_irbd_fc_results[[1]][1:2, ]
negative_fold <- lbd_vs_irbd_fc_results[[2]][1:5, ]
order_of_species_pos <- positive_fold[order(positive_fold$median_g1, decreasing = T), ]$Row.names[1:(min(nrow(lbd_vs_irbd_fc_results[[1]]), 5))]
order_of_species_neg <- negative_fold[order(negative_fold$median_g2, decreasing = T), ]$Row.names[1:(min(nrow(lbd_vs_irbd_fc_results[[2]]), 5))]
order_of_species <- c(order_of_species_pos, order_of_species_neg)
boxplot_data <- lbd_vs_irbd_diff_abun_median[[2]][, 1:36]

boxplot_long_df <- melt(boxplot_data, id.vars = "Row.names")
colnames(boxplot_long_df)[1] <- "species"
boxplot_long_df$condition <- ifelse(grepl("lbd", boxplot_long_df$variable), "LBD", "iRBD")
boxplot_long_df_done <- boxplot_long_df[boxplot_long_df$species %in% order_of_species, ]
boxplot_long_df_done$species <- factor(boxplot_long_df_done$species, levels = order_of_species)
boxplot_long_df_done$condition <- factor(boxplot_long_df_done$condition, levels = c("LBD", "iRBD"))

pdf(file = paste0("lbd_vs_irbd_boxplot.pdf"), height = 6, width = 10)
plot <- ggplot(boxplot_long_df_done, 
               aes(x = species, y = value, fill = condition)) + 
  geom_boxplot(outlier.color = "black", 
               outlier.fill = "white", 
               outlier.size = 2, 
               outlier.shape = 21) +
  # Add jitter if desired:
  # geom_point(aes(fill = condition), shape = 21, size = 1.5, 
  #            position = position_jitter(width = 0.2), color = "black") +
  scale_fill_manual(values = c("#ff9274", "#fec849")) +
  # ylim(as.numeric(0), as.numeric(0.5)) + 
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Relative abundance (arcsine square-root transformed)")
print(plot)
dev.off()

# 9. Calculate the difference between mean based on the relative abundance data before arcsine square-root transformed)
sig_lbd_vs_irbd_species_relab_og <- lbd_vs_irbd_relab[, colnames(lbd_vs_irbd_relab) %in% c(order_of_species, "condition")]

# Re-name the column names so it contains the group labels
sig_lbd_vs_irbd_species_relab_og$new_id <- paste(sig_lbd_vs_irbd_species_relab_og$condition, rownames(sig_lbd_vs_irbd_species_relab_og), sep = "_")
rownames(sig_lbd_vs_irbd_species_relab_og) <- sig_lbd_vs_irbd_species_relab_og$new_id
sig_lbd_vs_irbd_species_relab_og <- sig_lbd_vs_irbd_species_relab_og[, 1:(length(order_of_species))]
sig_lbd_vs_irbd_species_relab_og_t <- as.data.frame(t(sig_lbd_vs_irbd_species_relab_og))

# Compute mean for each group
sig_lbd_vs_irbd_species_relab_og_t$mean_lbd <- apply(sig_lbd_vs_irbd_species_relab_og_t[, grepl("lbd_BIOME", colnames(sig_lbd_vs_irbd_species_relab_og_t))], 1, function(x) mean(x))
sig_lbd_vs_irbd_species_relab_og_t$mean_irbd <- apply(sig_lbd_vs_irbd_species_relab_og_t[, grepl("irbd_BIOME", colnames(sig_lbd_vs_irbd_species_relab_og_t))], 1, function(x) mean(x))

sig_lbd_vs_irbd_species_relab_og_t$mean_diff <- sig_lbd_vs_irbd_species_relab_og_t$mean_lbd - sig_lbd_vs_irbd_species_relab_og_t$mean_irbd
sig_lbd_vs_irbd_species_relab_og_t$species <- rownames(sig_lbd_vs_irbd_species_relab_og_t)

sig_lbd_vs_irbd_species_relab_og_t$species <- factor(sig_lbd_vs_irbd_species_relab_og_t$species, levels = order_of_species)

pdf(file = "lbd_vs_irbd/mean_diff_bar_lbd_irbd_v2.pdf", height = 3.5, width = 8)
ggplot(sig_lbd_vs_irbd_species_relab_og_t, aes(x = as.factor(species), y = mean_diff)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "mean difference")
dev.off()

# 9. Create boxplot for each significant species
for (i in 1:nrow(sig_lbd_vs_irbd_diff_abun_species)) {
  plot <- ggplot(boxplot_long_df[boxplot_long_df$species == sig_lbd_vs_irbd_diff_abun_species$species[i], ], 
                 aes(x = condition, y = value, fill = condition)) + 
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = c("#ff9274", "#2DA248")) +
    theme_bw() + 
    theme(
      panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()) +
    labs(x = "", 
         y = "Microbial relative abundance (Arcsine square root transformed)",
         title = order_of_species[i])
  
  pdf(file = paste0("lbd_vs_irbd_boxplot_", sig_lbd_vs_irbd_diff_abun_species$species[i], ".pdf"), height = 6, width = 4)
  print(plot)
  dev.off()
}

## Species that are higher in LBD compare to iRBD and control
intersect(lbd_vs_control_fc_results[[1]]$Row.names, lbd_vs_irbd_fc_results[[1]]$Row.names)
## Species that are higher in iRBD compare to LBD and control
intersect(irbd_vs_control_fc_results[[1]]$Row.names, lbd_vs_irbd_fc_results[[2]]$Row.names)
