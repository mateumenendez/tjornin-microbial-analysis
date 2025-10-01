library(tidyverse)
library(janitor)
library(ggpubr)
# install.packages("phytools")
# library(phytools)
library(lvplot)
library(phyloseq)
library(ggdensity)
library(showtext)
library(ggh4x)
showtext_auto()

setwd(dir = "/projects/caeg/people/ngm902/apps/repos/tjornin-microbial-analysis")

source("libs/lib.R")
# source("libs/get-metadata.R")
# source("libs/get-gtdb-tree.R")
source("/projects/fernandezguerra/people/ngm902/scripts/r-miscellaneous.R")

# showtext::showtext_auto()
# showtext::showtext_opts(dpi = 300)
showtext::showtext_auto(FALSE)



# Let's read the kegg data
# kegg_data_gtdb <- read_tsv("/projects/fernandezguerra/people/ngm902/functionalDB/kegg-modules-summary-gtdb_207-rep.tsv.gz", col_types = cols(gene_caller_ids_in_module="c")) |>
kegg_data_gtdb <- read_tsv("/projects/caeg/people/ngm902/apps/repos/iceland-microbial-analysis/data/cdata/kegg-modules-summary-gtdb_207-rep.tsv.gz", col_types = cols(gene_caller_ids_in_module="c")) |>
    mutate(reference = genome_name) 
    # filter(pathwise_module_completeness >= 0.9) |>

# kegg_data_nongtdb <- read_tsv("/projects/fernandezguerra/people/ngm902/functionalDB/kegg-modules-summary-non-gtdb_207-rep.tsv.gz", col_types = cols(gene_caller_ids_in_module = "c")) |>
kegg_data_nongtdb <- read_tsv("/projects/caeg/people/ngm902/apps/repos/iceland-microbial-analysis/data/cdata/kegg-modules-summary-non-gtdb_207-rep-iceland.tsv.gz", col_types = cols(gene_caller_ids_in_module = "c")) |>
    mutate(reference = genome_name) 

kegg_data <- kegg_data_gtdb |> 
    bind_rows(kegg_data_nongtdb) |>
    dplyr::select(-label)


# Get all metadata that we might need
# Samples metadata
# Let's load the cdata
cdata <- read.table(file = "data/metadata.txt", sep = "\t", header = T) %>%
    rename(year = y_bp)



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
        filtered_data <- filtered_data %>% filter(grepl(enzyme_pattern, enzyme_hits_in_module))
    }

    filtered_data %>%
        group_by(label, year, depth_in_core_cm, module, module_name) %>%
        summarise(abundance = sum(abundance), .groups = "drop")
    }

    nfix <- process_module(taxa_data_phyloseq_df, "M00175", 0.6, "K02586|K02591")
    nitri <- process_module(taxa_data_phyloseq_df, "M00528", 0.5, "K10944")
    dnitri <- process_module(taxa_data_phyloseq_df, "M00529", 0.75, "K15864|K00368|K00376")
    d_nitr_red <- process_module(taxa_data_phyloseq_df, "M00530", 0.9)
    a_nitr_red <- process_module(taxa_data_phyloseq_df, "M00531", 0.9)
    commamox <- process_module(taxa_data_phyloseq_df, "M00804", 0.75)
    anammox <- process_module(taxa_data_phyloseq_df, "M00973", 0.75)

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

    # Six color palette
    modules <- c(
        "Nitrogen fixation, nitrogen => ammonia", # M00175
        "Assimilatory nitrate reduction, nitrate => ammonia", # M00531
        "Dissimilatory nitrate reduction, nitrate => ammonia", # M00530
        "Denitrification, nitrate => nitrogen", # M00529
        "Nitrification, ammonia => nitrite", # M00528
        "Complete nitrification, comammox, ammonia => nitrite => nitrate" # M00804
        # "Anammox, nitrite + ammonia => nitrogen" # M00973
        # "Methanogenesis"
    )

    colors_n <- c("#5c7d85", "#adb9a3", "#fff5c0", "#e6b265", "#b15555", "#763c3e")
    names(colors_n) <- modules

    return(list(df = df1, colors = colors_n))
    # return(df1, colors_n)

}


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
        filtered_data <- filtered_data %>% filter(grepl(enzyme_pattern, enzyme_hits_in_module))
    }

    filtered_data %>%
        group_by(label, year, depth_in_core_cm, module, module_name) %>%
        summarise(abundance = sum(abundance), .groups = "drop")
    }

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
  
    # Six color palette
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




### FILTERED COMMUNITY ###
# Let's get the taxonomic profiles, we will use the tad80 estimated abundances when available
tax_data <- read_tsv("results/aggregated/dmg-summary-ssp.tsv.gz") %>%
    rename(year = y_bp)

tax_data_stats <- tax_data |>
    filter(norm_entropy > 0.6, breadth > 0.01, breadth_exp_ratio >= 0.75) |>
    mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad))

tax_data_stats %>% pull(species) %>% unique() %>% length()

phylofilt <- readRDS("results/aggregated/filtered_ps_nmds.rds")

tax_data_stats <- tax_data_stats %>%
    filter(species %in% taxa_names(phylofilt))



# let's create a phyloseq object
tax_data_df <- tax_data_stats %>%
    # filter(domain != "d__Eukaryota", domain != "d__Viruses") %>%
    dplyr::select(label, name, abundance) %>%
    distinct() %>%
    pivot_wider(names_from = "label", values_from = abundance, values_fill = 0) %>%
    as.data.frame() %>%
    column_to_rownames("name")

taxonomy_df <- tax_data_stats %>%
    dplyr::select(domain, phylum, class, order, family, genus, species, subspecies) %>%
    distinct() %>%
    filter(subspecies %in% rownames(tax_data_df)) %>%
    filter(domain != "d__Eukaryota") |> 
    mutate(name = subspecies) %>%
    distinct() %>%
    as.data.frame() %>%
    column_to_rownames("name")

cdata_df <- cdata %>%
    distinct() %>%
    filter(label %in% colnames(tax_data_df)) %>%
    mutate(name = label) %>%
    distinct() %>%
    as.data.frame() %>%
    remove_rownames() %>%  # Ensure no row names exist
    column_to_rownames("name")

# phyloseq object
tax_data_phyloseq <- phyloseq::phyloseq(
    otu_table(tax_data_df, taxa_are_rows = TRUE),
    tax_table(as.matrix(taxonomy_df)),
    sample_data(cdata_df)
)

# Let's do a very relaxed filtering to remove low abundant taxa
tax_data_phyloseq_filt <- filt_physeq(tax_data_phyloseq, min_counts = 1000, ncounts = 3, nsites = 0, vcoeff = 0, mean_prop_thresh = 1e-5) %>%
    subset_taxa(., domain != "d__Viruses")

# Let's get proportions
tax_data_phyloseq_df <- tax_data_phyloseq_filt %>%
    # tax_glom(taxrank = "phylum") %>%
    microbiome::transform("compositional") %>%
    speedyseq::psmelt() %>%
    as_tibble()



#### RUN IT ####
df1_n <- nitrogen_metabolism(tax_data_phyloseq_df, kegg_data)



plotAn <- ggplot() +
    geom_area(data = df1_n$df, aes(year, log_fold_change_oldest, fill = module_name), color = "black", alpha = 1,  ) +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    scale_x_continuous(breaks = seq(200, 2000, 100))+
    theme(
        text = element_text(size = 14),
        legend.position = "bottom",
        strip.text.x = element_text(hjust = 0),
        aspect.ratio = 1.5/1
    ) +
    guides(fill = guide_legend(title = NULL, ncol = 1)) +
    scale_fill_manual(values = df1_n$colors) +
    ylab("Fold change") +
    xlab("")+
    coord_flip()

png(file = "plots/fig_SX8.png", width = 7, height = 10, units = "in", res = 600)
plotAn
dev.off()

plotBn <- ggplot() +
    geom_area(data = df1_n$df, aes(year, log_fold_change_oldest, fill = module_name), color = "black", alpha = 1, position = "stack") +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    facet_nested(~module_name)+
    scale_x_continuous(breaks = seq(200, 2000, 100))+
    guides(fill = guide_legend(title = NULL, nrow = 2)) +
    scale_fill_manual(values = df1_n$colors) +
    ylab("Fold change") +
    coord_flip()+
    xlab("") +
    theme(
        text = element_text(size = 14),
        legend.position = "bottom",
        # panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        # strip.text.x = element_text(hjust = 0),
        strip.text.x = element_blank(), # This removes facet titles
        aspect.ratio = 3/1,
    )

png(file = "plots/fig_SX9.png", width = 15, height = 10, units = "in", res = 600)
plotBn
dev.off()

png(file = "plots/fig_SX10.png", width = 19, height = 10, units = "in", res = 600)
egg::ggarrange(plotAn, plotBn + theme(legend.position = "none"), nrow = 1, widths = c(1,3))
dev.off()




df1_s <- sulfur_metabolism(tax_data_phyloseq_df, kegg_data)


plotAs <- ggplot() +
    geom_area(data = df1_s$df, aes(year, log_fold_change_oldest, fill = module_name), color = "black", alpha = 1, position = "stack") +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    scale_x_continuous(breaks = seq(200, 2000, 100))+
    theme(
        text = element_text(size = 14),
        legend.position = "bottom",
        # panel.grid = element_blank(),
        # strip.background = element_blank(),
        # strip.placement = "outside",
        strip.text.x = element_text(hjust = 0),
        aspect.ratio = 1.5/1
    ) +
    guides(fill = guide_legend(title = NULL, ncol = 1)) +
    scale_fill_manual(values = df1_s$colors) +
    ylab("Fold change") +
    xlab("")+
    coord_flip()

png(file = "plots/fig_SX11.png", width = 7, height = 10, units = "in", res = 600)
plotAs
dev.off()

plotBs <- ggplot() +
    geom_area(data = df1_s$df, aes(year, log_fold_change_oldest, fill = module_name), color = "black", alpha = 1, position = "stack") +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    facet_nested(~module_name)+
    scale_x_continuous(breaks = seq(200, 2000, 100))+
    guides(fill = guide_legend(title = NULL, nrow = 2)) +
    scale_fill_manual(values = df1_s$colors) +
    ylab("Fold change") +
    coord_flip()+
    xlab("") +
    theme(
        text = element_text(size = 14),
        legend.position = "bottom",
        # panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        # strip.text.x = element_text(hjust = 0),
        strip.text.x = element_blank(), # This removes facet titles
        aspect.ratio = 3/1,
    )

png(file = "plots/fig_SX12.png", width = 15, height = 10, units = "in", res = 600)
plotBs
dev.off()

png(file = "plots/fig_SX13.png", width = 19, height = 10, units = "in", res = 600)
egg::ggarrange(plotAs, plotBs + theme(legend.position = "none"), nrow = 1, widths = c(1,3))
dev.off()

png(file = "plots/fig_SX14.png", width = 12, height = 10, units = "in", res = 600)
egg::ggarrange(plotAn + labs(title = "Nitrogen metabolism"), plotAs + labs(title = "Sulfur metabolism"), nrow = 1)
dev.off()

