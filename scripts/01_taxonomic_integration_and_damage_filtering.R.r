# ======================================================================
# Script: 01_taxonomic_integration_and_damage_filtering.R
# Author: Mateu Menendez Serra
# Affiliation: University of Copenhagen
# Contact: mateu.menendez@sund.ku.dk
# Date: October 2025
#
# Description:
#   This script performs the integration of the microbial taxonomic profiling 
# 	lca and damage estimation at subspecies and species levels using 
# 	metaDMG and filterBAM outputs.
#   It computes Lin’s concordance correlation coefficient (CCC) to
#   evaluate the fit of the observed DNA damage patterns to the expected ones,
#	classifies taxa as “Damaged” or “Non-damaged”, and generates abundance tables that
#   integrate taxonomic, abundance, and damage information.
#   Finally, it builds phyloseq objects at subspecies and species levels
#   for downstream diversity and ordination analyses.
#
# Input:
#   - data/metadata.txt
#   - data/metadmg/tp-lca.summary.tsv.gz
#   - data/metadmg/tp-damage.tsv.gz
#   - data/filterbam/tp-mapping-filtered.summary.tsv.gz
#
# Output:
#   - results/dmg-subspecies-lite.tsv.gz
#   - results/aggregated/dmg-species-lite.tsv.gz
#   - results/dmg-summary-ssp.tsv.gz
#   - results/dmg-summary-sp.tsv.gz
#   - results/ps-unfilt-ssp.rds
#   - results/ps-unfilt-sp.rds
#
# Dependencies:
#   tidyverse, phyloseq
#
# Note:
#   - dmg_fwd_CCC(): function for computing Lin’s CCC, defined in libs/dmg.R
#   - get_st(): optional normalization/filtering function, defined in libs/lib.R
#
# ======================================================================


# ----------------------------------------------------------------------
# (1) Setup environment and load libraries
# ----------------------------------------------------------------------
# Load required packages and source custom functions.

setwd(dir = "")
source("libs/dmg.R")
source("libs/lib.R")

library(tidyverse)
library(phyloseq)

# ----------------------------------------------------------------------
# (2) Load metadata and LCA taxonomic assignments
# ----------------------------------------------------------------------
# Load core sample metadata and metaDMG LCA summary, keeping only
# subspecies and species ranks, and excluding Eukaryotic hits.

cdata <- read.table(file = "data/metadata.txt", sep = "\t", header = TRUE)

tax_data <- read_tsv("data/metadmg/tp-lca.summary.tsv.gz") %>%
    filter(rank %in% c("subspecies", "species")) %>%
    separate(
        tax_path,
        into = c("root", "domain", "lineage", "kingdom", "phylum", "class",
                 "order", "family", "genus", "species", "subspecies"),
        sep = ";"
    ) %>%
    filter(domain != "d__Eukaryota") %>%
    inner_join(cdata, by = "label")



# ----------------------------------------------------------------------
# (3) Evaluate DNA damage model fit at subspecies level
# ----------------------------------------------------------------------
# Assess Lin’s concordance correlation coefficient (CCC) for taxa
# with >100 reads at the subspecies level, and classify fits as
# “good” or “bad”.

dmg <- read_tsv("data/metadmg/tp-damage.tsv.gz") %>%
    filter(label %in% cdata$label) %>%
    inner_join(tax_data %>% select(label, taxid, name, rank), by = c("label", "taxid"))

dmg_ssp <- dmg %>%
    filter(rank == "subspecies", nreads >= 100) %>%
    rename(tax_name = taxid, n_reads = nreads)

samples <- unique(cdata$label)

dat <- dmg_fwd_CCC(dmg_ssp, samples, ci = "asymptotic", nperm = 100, nproc = 24)

dat1 <- dat %>%
    ungroup() %>%
    inner_join(dmg_ssp, by = c("label", "tax_name")) %>%
    mutate(
        fit = ifelse(rho_c >= 0.85 & C_b > 0.9 &
                     round(rho_c_perm_pval, 3) < 0.1 &
                     !is.na(rho_c), "good", "bad"),
        fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit)
    )

dat1 %>%
    select(label, tax_name, A_b, q_b, c_b, phi_b, fit) %>%
    write_tsv("results/dmg-subspecies-lite.tsv.gz")



# ----------------------------------------------------------------------
# (4) Evaluate DNA damage model fit at species level
# ----------------------------------------------------------------------
# Repeat CCC computation and fit classification for species-level taxa.

dmg_sp <- dmg %>%
    filter(rank == "species", nreads >= 100) %>%
    rename(tax_name = taxid, n_reads = nreads)

samples <- unique(cdata$label)

dat <- dmg_fwd_CCC(dmg_sp, samples, ci = "asymptotic", nperm = 100, nproc = 24)

dat1 <- dat %>%
    ungroup() %>%
    inner_join(dmg_sp, by = c("label", "tax_name")) %>%
    mutate(
        fit = ifelse(rho_c >= 0.85 & C_b > 0.9 &
                     round(rho_c_perm_pval, 3) < 0.1 &
                     !is.na(rho_c), "good", "bad"),
        fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit)
    )

dat1 %>%
    select(label, tax_name, A_b, q_b, c_b, phi_b, fit) %>%
    write_tsv("results/aggregated/dmg-species-lite.tsv.gz")



# ----------------------------------------------------------------------
# (5) Combine taxonomic assignments and damage results (subspecies)
# ----------------------------------------------------------------------
# Integrate mapping stats and damage classification.
# Taxa with bad fits or A_b < 0.1 are marked as “Non-damaged”.

lca_ssp <- read_tsv("data/metadmg/tp-lca.summary.tsv.gz") %>%
    filter(rank == "subspecies") %>%
    separate(
        tax_path,
        into = c("root", "domain", "lineage", "kingdom", "phylum", "class",
                 "order", "family", "genus", "species", "subspecies"),
        sep = ";"
    ) %>%
    filter(domain != "d__Eukaryota") %>%
    rename(nreads = n_reads) %>%
    left_join(cdata)

mapping_stats <- read_tsv("data/filterbam/tp-mapping-filtered.summary.tsv.gz") %>%
    filter(norm_entropy > 0.6, breadth > 0.01, breadth_exp_ratio >= 0.75) %>%
    mutate(
        name = paste0("S__", reference),
        abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad)
    )

tax_data_ssp <- lca_ssp %>%
    select(-nreads, -abundance) |>
    inner_join(mapping_stats)

dmg_ssp <- read_tsv("results/dmg-subspecies-lite.tsv.gz") %>%
    rename(taxid = tax_name) %>%
    mutate(damage = ifelse(fit == "bad", 0, A_b))

tax_data_dmg <- tax_data_ssp %>%
    inner_join(dmg_ssp, by = c("label", "taxid")) %>%
    left_join(cdata, by = "label") %>%
    mutate(
        is_dmg = ifelse(damage > 0.1, "Damaged", "Non-damaged"),
        is_dmg = ifelse(fit == "bad", "Non-damaged", is_dmg)
    )

tax_data_dmg %>%
    write_tsv("results/dmg-summary-ssp.tsv.gz")



# ----------------------------------------------------------------------
# (6) Create phyloseq object at subspecies level
# ----------------------------------------------------------------------
# Build count, taxonomy, and metadata tables, and assemble a phyloseq object.

tax_data_dmg <- read_tsv("results/dmg-summary-ssp.tsv.gz")

tax_data_df <- tax_data_dmg %>%
    inner_join(cdata, by = "label") %>%
    select(label, subspecies, abundance) %>%
    pivot_wider(names_from = "label", values_from = abundance, values_fill = 0) %>%
    as.data.frame() %>%
    column_to_rownames("subspecies")

taxonomy_df <- tax_data_dmg %>%
    select(domain, phylum, class, order, family, genus, species, subspecies) %>%
    mutate(name = subspecies) %>%
    distinct() %>%
    as.data.frame() %>%
    column_to_rownames("name")

cdata_df <- cdata %>%
    filter(label %in% colnames(tax_data_df)) %>%
    mutate(name = label) %>%
    distinct() %>%
    as.data.frame()
rownames(cdata_df) <- cdata_df$name

ps_subspecies <- phyloseq(
    otu_table(tax_data_df, taxa_are_rows = TRUE),
    tax_table(as.matrix(taxonomy_df)),
    sample_data(cdata_df)
)

saveRDS(ps_subspecies, "results/ps-unfilt-ssp.rds")



# ----------------------------------------------------------------------
# (7) Combine taxonomic assignments and damage results (species)
# ----------------------------------------------------------------------
# Repeat integration and filtering for species-level taxa.

lca_sp <- read_tsv("data/metadmg/tp-lca.summary.tsv.gz") %>%
    filter(rank == "species") %>%
    separate(
        tax_path,
        into = c("root", "domain", "lineage", "kingdom", "phylum", "class",
                 "order", "family", "genus", "species", "subspecies"),
        sep = ";"
    ) %>%
    filter(domain != "d__Eukaryota") %>%
    rename(nreads = n_reads) %>%
    left_join(cdata)

dmg_sp <- read_tsv("results/dmg-species-lite.tsv.gz") %>%
    rename(taxid = tax_name) %>%
    mutate(damage = ifelse(fit == "bad", 0, A_b))

tax_data_dmg <- lca_sp %>%
    inner_join(dmg_sp) %>%
    ungroup() %>%
    left_join(cdata) %>%
    mutate(
        is_dmg = ifelse(damage > 0.1, "Damaged", "Non-damaged"),
        is_dmg = ifelse(fit == "bad", "Non-damaged", is_dmg)
    )

tax_data_dmg %>%
    write_tsv("results/dmg-summary-sp.tsv.gz")



# ----------------------------------------------------------------------
# (8) Create phyloseq object at species level
# ----------------------------------------------------------------------
# Assemble count, taxonomy, and metadata tables for species-level dataset.

tax_data_dmg_sp <- read_tsv("results/dmg-summary-sp.tsv.gz")

tax_data_df <- tax_data_dmg_sp %>%
    inner_join(cdata, by = "label") %>%
    select(label, species, abundance) %>%
    pivot_wider(names_from = "label", values_from = abundance, values_fill = 0) %>%
    as.data.frame() %>%
    column_to_rownames("species")

taxonomy_df <- tax_data_dmg_sp %>%
    select(domain, phylum, class, order, family, genus, species) %>%
    mutate(name = species) %>%
    distinct() %>%
    as.data.frame() %>%
    column_to_rownames("name")

cdata_df <- cdata %>%
    filter(label %in% colnames(tax_data_df)) %>%
    mutate(name = label) %>%
    distinct() %>%
    as.data.frame()
rownames(cdata_df) <- cdata_df$name

ps_species <- phyloseq(
    otu_table(tax_data_df, taxa_are_rows = TRUE),
    tax_table(as.matrix(taxonomy_df)),
    sample_data(cdata_df)
)

saveRDS(ps_species, "results/ps-unfilt-sp.rds")



# ----------------------------------------------------------------------
# End of script
# ----------------------------------------------------------------------