library(tidyverse)
library(knitr)
library(readxl)
library(measurements)
library(ggmap)
library(ggrepel)
library(lvplot)
library(cowplot)
library(phyloseq)


setwd(dir = "/projects/caeg/people/ngm902/apps/repos/tjornin-microbial-analysis")
source("/projects/fernandezguerra/people/ngm902/scripts/r-miscellaneous.R")

# Load tax combining damage and abundances
tax_data_dmg <- read_tsv("results/aggregated/dmg-summary-ssp.tsv.gz")

# Load labels data and taxonomy
cdata <- read.table(file = "data/metadata.txt", sep = "\t", header = T)

################################################
#### Data filtering ############################
################################################

tax_data_filt <- tax_data_dmg |>
    inner_join(cdata) |>
    select(label, subspecies, abundance)

# let's create a phyloseq object
tax_data_df <- tax_data_filt |>
    pivot_wider(names_from = "label", values_from = abundance, values_fill = 0) |>
    as.data.frame() |>
    column_to_rownames("subspecies")

taxonomy_df <- tax_data_dmg |>
    ungroup() |>
    select(domain, phylum, class, order, family, genus, species, subspecies) |>
    mutate(name = subspecies) |>
    distinct() |>
    as.data.frame() |>
    column_to_rownames("name")

cdata_df <- cdata |>
    filter(label %in% colnames(tax_data_df)) |>
    mutate(name = label) |>
    distinct() |>
    as.data.frame() 

rownames(cdata_df) <- cdata_df$name


# phyloseq object
tax_data_phyloseq <- phyloseq::phyloseq(
    otu_table(tax_data_df, taxa_are_rows = TRUE),
    tax_table(as.matrix(taxonomy_df)),
    sample_data(cdata_df)
)

# Let's save the phyloseq object in case we need it later
saveRDS(tax_data_phyloseq, "results/aggregated/ps-unfilt-ssp.rds")


ps_filt <- get_st(tax_data_phyloseq, nsites = 0.05, vcoeff = 0, mean_prop_thresh = 0, norm = "median", glom = NULL)
saveRDS(ps_filt, "results/aggregated/ps-stand-filt-ssp.rds")





# ## Gloom at species level
# tax_data_phyloseq_sp <- tax_glom(tax_data_phyloseq, taxrank = "species")
# # taxa_names(tax_data_phyloseq_sp)

# # Extraer los nuevos nombres de taxa
# new_taxa_names <- tax_table(tax_data_phyloseq_sp)[, "species"]  # Ajusta según la taxonomía disponible

# # Convertir los valores a caracteres y asignar a `taxa_names`
# taxa_names(tax_data_phyloseq_sp) <- as.character(new_taxa_names)
# saveRDS(tax_data_phyloseq_sp, "./results/aggregated/ps-unfilt-sp_glom.rds")









# Load tax combining damage and abundances
tax_data_dmg_sp <- read_tsv("results/aggregated/dmg-summary-sp.tsv.gz")

################################################
#### Data filtering ############################
################################################

tax_data_filt <- tax_data_dmg_sp |>
    inner_join(cdata) |>
    select(label, species, abundance)

# let's create a phyloseq object
tax_data_df <- tax_data_filt |>
    pivot_wider(names_from = "label", values_from = abundance, values_fill = 0) |>
    as.data.frame() |>
    column_to_rownames("species")

taxonomy_df <- tax_data_dmg_sp |>
    ungroup() |>
    select(domain, phylum, class, order, family, genus, species) |>
    mutate(name = species) |>
    distinct() |>
    as.data.frame() |>
    column_to_rownames("name")

cdata_df <- cdata |>
    filter(label %in% colnames(tax_data_df)) |>
    mutate(name = label) |>
    distinct() |>
    as.data.frame()

rownames(cdata_df) <- cdata_df$name


# phyloseq object
tax_data_phyloseq <- phyloseq::phyloseq(
    otu_table(tax_data_df, taxa_are_rows = TRUE),
    tax_table(as.matrix(taxonomy_df)),
    sample_data(cdata_df)
)

# Let's save the phyloseq object in case we need it later
saveRDS(tax_data_phyloseq, "./results/aggregated/ps-unfilt-sp.rds")

ps_filt <- get_st(tax_data_phyloseq, nsites = 0.05, vcoeff = 0, mean_prop_thresh = 0, norm = "median", glom = NULL)
saveRDS(ps_filt, "./results/aggregated/ps-stand-filt-sp.rds")
