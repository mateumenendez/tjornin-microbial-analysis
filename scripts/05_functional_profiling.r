# ======================================================================
# Script: 05_functional_profiling.R
# Supplementary Script — Nitrogen and Sulfur Metabolic Modules Analysis
# Author: Mateu Menendez Serra
# Affiliation: University of Copenhagen
# Date: October 2025
#
# Description:
#   This script quantifies the temporal dynamics of taxa associated to 
#   the nitrogen and sulfur metabolism by after the functional annotation
#   of the reference genomes based on KEGG.
#
#   The workflow includes:
#     1. Loading KEGG and metadata inputs
#     2. Defining functions for nitrogen and sulfur metabolism analysis
#     3. Constructing phyloseq objects with filtered taxa
#     4. Calculating module-level abundance changes (log fold change)
#     5. Visualizing temporal trends (Supplementary Figures SX8–SX14)
#
# Dependencies:
#   tidyverse, phyloseq,
#   ggh4x, egg, speedyseq, microbiome
# ======================================================================


# ----------------------------------------------------------------------
# (1) Load libraries and set working environment
# ----------------------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(ggdensity)
library(ggh4x)


setwd(dir = "")

# Load auxiliary functions
source("libs/lib.R")

# ----------------------------------------------------------------------
# (2) Load KEGG module completeness data
# ----------------------------------------------------------------------

kegg_data <- read_tsv("data/kegg-modules-summary-iceland.tsv.gz", col_types = cols(gene_caller_ids_in_module = "c"))


# ----------------------------------------------------------------------
# (3) Load metadata and define helper functions
# ----------------------------------------------------------------------

# Sample metadata
cdata <- read.table(file = "data/metadata.txt", sep = "\t", header = TRUE) %>%
  rename(year = y_bp)


# ----------------------------------------------------------------------
# (4) Define nitrogen metabolism analysis function
# ----------------------------------------------------------------------
nitrogen_metabolism <- function(taxa_data_phyloseq_df, kegg_data) {

  library(dplyr)
  library(forcats)

  process_module <- function(data, module_id, completeness_threshold, enzyme_pattern = NULL) {
    filtered_data <- data %>%
      filter(domain %in% c("d__Bacteria", "d__Archaea")) %>%
      mutate(reference = gsub("S__", "", OTU)) %>%
      rename(abundance = Abundance) %>%
      filter(abundance > 0) %>%
      inner_join(kegg_data %>% filter(module == module_id), by = "reference") %>%
      filter(pathwise_module_completeness >= completeness_threshold)

    if (!is.null(enzyme_pattern)) {
      filtered_data <- filtered_data %>%
        filter(grepl(enzyme_pattern, enzyme_hits_in_module))
    }

    filtered_data %>%
      group_by(label, year, depth_in_core_cm, module, module_name) %>%
      summarise(abundance = sum(abundance), .groups = "drop")
  }

  # Individual nitrogen pathways (KEGG modules)
  nfix <- process_module(taxa_data_phyloseq_df, "M00175", 0.6, "K02586|K02591")
  nitri <- process_module(taxa_data_phyloseq_df, "M00528", 0.5, "K10944")
  dnitri <- process_module(taxa_data_phyloseq_df, "M00529", 0.75, "K15864|K00368|K00376")
  d_nitr_red <- process_module(taxa_data_phyloseq_df, "M00530", 0.9)
  a_nitr_red <- process_module(taxa_data_phyloseq_df, "M00531", 0.9)
  commamox <- process_module(taxa_data_phyloseq_df, "M00804", 0.75)

  df1 <- bind_rows(nitri, nfix, dnitri, d_nitr_red, a_nitr_red, commamox) %>%
    group_by(module, module_name) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      log_fold_change = log2(abundance / lag(abundance)),
      log_fold_change_median = log2(abundance / median(abundance)),
      log_fold_change_oldest = log2(abundance / first(abundance))
    ) %>%
    ungroup()

  df1$log_fold_change[is.na(df1$log_fold_change)] <- 0

  df1 <- df1 %>%
    mutate(module_name = fct_reorder(module_name, -log_fold_change_oldest))

  # Custom color palette for nitrogen metabolism
  modules <- c(
    "Nitrogen fixation, nitrogen => ammonia",
    "Assimilatory nitrate reduction, nitrate => ammonia",
    "Dissimilatory nitrate reduction, nitrate => ammonia",
    "Denitrification, nitrate => nitrogen",
    "Nitrification, ammonia => nitrite",
    "Complete nitrification, comammox, ammonia => nitrite => nitrate"
  )
  colors_n <- c("#5c7d85", "#adb9a3", "#fff5c0", "#e6b265", "#b15555", "#763c3e")
  names(colors_n) <- modules

  return(list(df = df1, colors = colors_n))
}


# ----------------------------------------------------------------------
# (5) Define sulfur metabolism analysis function
# ----------------------------------------------------------------------
sulfur_metabolism <- function(taxa_data_phyloseq_df, kegg_data) {

  library(dplyr)
  library(forcats)

  process_module <- function(data, module_id, completeness_threshold, enzyme_pattern = NULL) {
    filtered_data <- data %>%
      filter(domain %in% c("d__Bacteria", "d__Archaea")) %>%
      mutate(reference = gsub("S__", "", OTU)) %>%
      rename(abundance = Abundance) %>%
      filter(abundance > 0) %>%
      inner_join(kegg_data %>% filter(module == module_id), by = "reference") %>%
      filter(pathwise_module_completeness >= completeness_threshold)

    if (!is.null(enzyme_pattern)) {
      filtered_data <- filtered_data %>%
        filter(grepl(enzyme_pattern, enzyme_hits_in_module))
    }

    filtered_data %>%
      group_by(label, year, depth_in_core_cm, module, module_name) %>%
      summarise(abundance = sum(abundance), .groups = "drop")
  }

  # Individual sulfur pathways (KEGG modules)
  a_sul_red <- process_module(taxa_data_phyloseq_df, "M00176", 0.9, "K00390")
  d_sul_red <- process_module(taxa_data_phyloseq_df, "M00596", 0.75)
  thiosul_ox <- process_module(taxa_data_phyloseq_df, "M00595", 1)
  sulf_sulf_ass <- process_module(taxa_data_phyloseq_df, "M00616", 1)

  df1 <- bind_rows(a_sul_red, d_sul_red, thiosul_ox, sulf_sulf_ass) %>%
    group_by(module, module_name) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      log_fold_change = log2(abundance / lag(abundance)),
      log_fold_change_median = log2(abundance / median(abundance)),
      log_fold_change_oldest = log2(abundance / first(abundance))
    ) %>%
    ungroup()

  df1$log_fold_change[is.na(df1$log_fold_change)] <- 0

  df1 <- df1 %>%
    mutate(module_name = fct_reorder(module_name, -log_fold_change_oldest))

  # Custom color palette for sulfur metabolism
  modules <- c(
    "Assimilatory sulfate reduction, sulfate => H2S",
    "Thiosulfate oxidation by SOX complex, thiosulfate => sulfate",
    "Dissimilatory sulfate reduction, sulfate => H2S",
    "Sulfate-sulfur assimilation"
  )
  colors_s <- c("#0a9396", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226")
  names(colors_s) <- modules

  return(list(df = df1, colors = colors_s))
}


# ----------------------------------------------------------------------
# (6) Load community taxonomic profiles
# ----------------------------------------------------------------------

tax_data_stats <- read_tsv("results/dmg-summary-ssp.tsv.gz") %>%
  rename(year = y_bp)

phylofilt <- readRDS("results/filtered_ps_nmds.rds")

tax_data_stats <- tax_data_stats %>%
  filter(species %in% taxa_names(phylofilt))


# ----------------------------------------------------------------------
# (7) Construct phyloseq object
# ----------------------------------------------------------------------

tax_data_df <- tax_data_stats %>%
  dplyr::select(label, name, abundance) %>%
  distinct() %>%
  pivot_wider(names_from = "label", values_from = abundance, values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("name")

taxonomy_df <- tax_data_stats %>%
  dplyr::select(domain, phylum, class, order, family, genus, species, subspecies) %>%
  distinct() %>%
  filter(subspecies %in% rownames(tax_data_df)) %>%
  filter(domain != "d__Eukaryota") %>%
  mutate(name = subspecies) %>%
  as.data.frame() %>%
  column_to_rownames("name")

cdata_df <- cdata %>%
  distinct() %>%
  filter(label %in% colnames(tax_data_df)) %>%
  mutate(name = label) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames("name")

tax_data_phyloseq <- phyloseq::phyloseq(
  otu_table(tax_data_df, taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy_df)),
  sample_data(cdata_df)
)

# Relaxed filtering to remove low-abundance taxa
tax_data_phyloseq_filt <- filt_physeq(
  tax_data_phyloseq,
  min_counts = 1000, ncounts = 3, nsites = 0, vcoeff = 0, mean_prop_thresh = 1e-5
) %>%
  subset_taxa(domain != "d__Viruses")

# Transform to relative abundances
tax_data_phyloseq_df <- tax_data_phyloseq_filt %>%
  microbiome::transform("compositional") %>%
  speedyseq::psmelt() %>%
  as_tibble()

# ----------------------------------------------------------------------
# (8) Run nitrogen metabolism analysis
# ----------------------------------------------------------------------

df1_n <- nitrogen_metabolism(tax_data_phyloseq_df, kegg_data)

# Supplementary Figures SX8–SX10
# SX8: Combined area plot, SX9: Faceted by module, SX10: Combined layout
# ----------------------------------------------------------------------

plotAn <- ggplot() +
  geom_area(data = df1_n$df, aes(year, log_fold_change_oldest, fill = module_name),
            color = "black", alpha = 1) +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  scale_x_continuous(breaks = seq(200, 2000, 100)) +
  theme(text = element_text(size = 14),
        legend.position = "bottom",
        strip.text.x = element_text(hjust = 0),
        aspect.ratio = 1.5 / 1) +
  guides(fill = guide_legend(title = NULL, ncol = 1)) +
  scale_fill_manual(values = df1_n$colors) +
  ylab("Fold change") +
  xlab("") +
  coord_flip()

png("plots/fig_SX8.png", width = 7, height = 10, units = "in", res = 600)
plotAn
dev.off()


plotBn <- ggplot() +
  geom_area(data = df1_n$df, aes(year, log_fold_change_oldest, fill = module_name),
            color = "black", alpha = 1, position = "stack") +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  facet_nested(~module_name) +
  scale_x_continuous(breaks = seq(200, 2000, 100)) +
  guides(fill = guide_legend(title = NULL, nrow = 2)) +
  scale_fill_manual(values = df1_n$colors) +
  ylab("Fold change") +
  coord_flip() +
  xlab("") +
  theme(text = element_text(size = 14),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_blank(),
        aspect.ratio = 3 / 1)

png("plots/fig_SX9.png", width = 15, height = 10, units = "in", res = 600)
plotBn
dev.off()

png("plots/fig_SX10.png", width = 19, height = 10, units = "in", res = 600)
egg::ggarrange(plotAn, plotBn + theme(legend.position = "none"), nrow = 1, widths = c(1, 3))
dev.off()


# ----------------------------------------------------------------------
# (9) Run sulfur metabolism analysis
# ----------------------------------------------------------------------

df1_s <- sulfur_metabolism(tax_data_phyloseq_df, kegg_data)

# Supplementary Figures SX11–SX14
# SX11: Combined area plot, SX12: Faceted by module,
# SX13: Combined layout, SX14: Comparative nitrogen–sulfur
# ----------------------------------------------------------------------

plotAs <- ggplot() +
  geom_area(data = df1_s$df, aes(year, log_fold_change_oldest, fill = module_name),
            color = "black", alpha = 1, position = "stack") +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  scale_x_continuous(breaks = seq(200, 2000, 100)) +
  theme(text = element_text(size = 14),
        legend.position = "bottom",
        strip.text.x = element_text(hjust = 0),
        aspect.ratio = 1.5 / 1) +
  guides(fill = guide_legend(title = NULL, ncol = 1)) +
  scale_fill_manual(values = df1_s$colors) +
  ylab("Fold change") +
  xlab("") +
  coord_flip()

png("plots/fig_SX11.png", width = 7, height = 10, units = "in", res = 600)
plotAs
dev.off()

plotBs <- ggplot() +
  geom_area(data = df1_s$df, aes(year, log_fold_change_oldest, fill = module_name),
            color = "black", alpha = 1, position = "stack") +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  facet_nested(~module_name) +
  scale_x_continuous(breaks = seq(200, 2000, 100)) +
  guides(fill = guide_legend(title = NULL, nrow = 2)) +
  scale_fill_manual(values = df1_s$colors) +
  ylab("Fold change") +
  coord_flip() +
  xlab("") +
  theme(text = element_text(size = 14),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_blank(),
        aspect.ratio = 3 / 1)

png("plots/fig_SX12.png", width = 15, height = 10, units = "in", res = 600)
plotBs
dev.off()

png("plots/fig_SX13.png", width = 19, height = 10, units = "in", res = 600)
egg::ggarrange(plotAs, plotBs + theme(legend.position = "none"), nrow = 1, widths = c(1, 3))
dev.off()

png("plots/fig_SX14.png", width = 12, height = 10, units = "in", res = 600)
egg::ggarrange(
  plotAn + labs(title = "Nitrogen metabolism"),
  plotAs + labs(title = "Sulfur metabolism"),
  nrow = 1
)
dev.off()


# ----------------------------------------------------------------------
# End of script
# ----------------------------------------------------------------------

