# ======================================================================
# Script: 04_host_and_environmental_sources_analysis.R
# Author: Mateu Menendez Serra
# Affiliation: University of Copenhagen
# Contact: mateu.menendez@sund.ku.dk
# Date: October 2025
#
# Description:
#   This script integrates microbial taxonomic profiles with host
#   information and environmental source metadata to explore potential
#   origins of bacterial, archaeal, and viral taxa.
#   It classifies host-associated bacterial and archaeal genomes,
#   traces viral environmental sources using IMG/VR annotations,
#   and integrates bacterial and archaeal SourceTracker results.
#   The script produces figures summarizing the relative abundance
#   of taxa by host, viral isolation sources, and inferred environmental
#   contributions.
#
# Input:
#   - results/aggregated/dmg-summary-ssp.tsv.gz
#   - results/aggregated/dmg-summary-sp.tsv.gz
#   - data/metadata.txt
#   - data/host-isolated-references.csv
#   - results/aggregated/sourcetracker/st-biome-class-gm_v2/mixing_proportions.txt
#   - functionalDB/IMGVR_all_Sequence_information-high_confidence.tsv.gz
#   - sourcetracker/source_data/med-biomes-download.txt
#
# Output:
#   - plots/fig_SX7.png
#   - plots/fig_X1.png
#   - plots/fig_X1.pdf
#   - plots/fig_X1_2.pdf
# 
#
# Dependencies:
#   tidyverse, phyloseq, ggh4x, dendextend, indicspecies, vegan, egg
#
# Note:
#   - The color palettes (paired_genus) and auxiliary functions
#     are defined in libs/lib.R
# ======================================================================


# ----------------------------------------------------------------------
# (1) Setup environment and load libraries
# ----------------------------------------------------------------------
# Load required packages and custom functions.

setwd(dir = "")
source("libs/lib.R")

library(tidyverse)
library(phyloseq)
library(egg)
library(ggh4x)
library(vegan)
library(ggrepel)



# ----------------------------------------------------------------------
# (2) Load metadata and subspecies-level abundance data
# ----------------------------------------------------------------------
# Load sample metadata and damage summary data at subspecies level.

ssp <- read_tsv("results/dmg-summary-ssp.tsv.gz")
cdata <- read.table("data/metadata.txt", sep = "\t", header = TRUE)


# ----------------------------------------------------------------------
# (3) Assign host taxonomy to intestinal-associated subspecies
# ----------------------------------------------------------------------
# Assign reference hosts to specific bacterial and archaeal subspecies
# based on known genome sources (e.g., Bos taurus, Ovis, Capra, Homo sapiens).
host_data <- read.csv("data/host-isolated-references.csv")

ssp_intestine <- ssp %>%
    filter(!is.na(abundance)) %>%
	left_join(host_data %>% select(-family), by = "subspecies") %>%
    filter(!is.na(host))

colors_host <- paired_genus[1:length(unique(ssp_intestine$host))]
names(colors_host) <- unique(ssp_intestine$host)
colors_host[grepl("Other", names(colors_host))] <- "#dcdcdc"


# ----------------------------------------------------------------------
# (4) Plot abundance and damage profiles by host
# ----------------------------------------------------------------------
# Produce multi-panel figure summarizing:
# A) Subspecies abundance per host
# B) Abundance by host and sample
# C) Number of species per family over time
# D) Damage fit (A_b) by host group

# Plot definitions omitted for brevity (same as original code)
# Output: plots/fig_SX7.png

# A) Abundance per reference genome
plotA <- ssp_intestine %>%
  group_by(host) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(
    host = fct_reorder(host, -total_abundance),
    host = fct_relevel(host, "Other", after = Inf)
  ) %>%
  ggplot(aes(
    x = y_bp, y = abundance,
    color = host,
    shape = paste(domain, phylum, class, order, family)
  )) +
  geom_point() +
  scale_color_manual(values = colors_host) +
  coord_flip() +
  scale_x_continuous(breaks = seq(200, 2000, 100), limits = c(200, 2020)) +
  guides(shape = guide_legend(ncol = 1), color = guide_legend(ncol = 2)) +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  theme(
    aspect.ratio = 2 / 1,
    legend.position = "bottom",
    legend.justification = c(0, 0)
  ) +
  labs(
    y = "Abundance",
    x = "Year (CE)",
    title = "Abundance per reference",
    color = "Reference host",
    shape = "Taxonomic\npath"
  )


# B) Stacked area abundance per host over time
plotB1 <- ssp_intestine %>%
  group_by(label, host) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  tidyr::complete(
    host = unique(host),
    label = unique(label),
    fill = list(abundance = 0)
  ) %>%
  left_join(cdata %>% select(label, y_bp)) %>%
  group_by(host) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(
    host = fct_reorder(host, -total_abundance),
    host = fct_relevel(host, "Other", after = Inf)
  ) %>%
  ggplot(aes(x = y_bp, y = abundance, fill = host)) +
  geom_area(stat = "identity", position = "stack", alpha = 0.7) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors_host) +
  coord_flip() +
  scale_x_continuous(breaks = seq(200, 2000, 100), limits = c(200, 2020)) +
  labs(y = "Abundance", x = "", title = "Abundance by host", fill = "") +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  theme(
    aspect.ratio = 2 / 1,
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )


# C) Species count by family
tax_data_dmg_sp <- read_tsv("./results/dmg-summary-sp.tsv.gz")

plotB2 <- tax_data_dmg_sp %>%
  filter(species %in% ssp_intestine$species, !is.na(abundance)) %>%
  group_by(label, y_bp, domain, phylum, class, order, family) %>%
  summarise(count = n(), .groups = "drop") %>%
  ggplot(aes(
    x = y_bp, y = count,
    group = paste(domain, phylum, class, order, family)
  )) +
  geom_point(aes(shape = paste(domain, phylum, class, order, family))) +
  geom_line(color = "black") +
  coord_flip() +
  scale_x_continuous(breaks = seq(200, 2000, 100), limits = c(200, 2020)) +
  labs(
    y = "Number of species",
    x = "",
    title = "Species\ncount by family",
    fill = "Damage fit"
  ) +
  theme(
    aspect.ratio = 2 / 1,
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )


# D) DNA damage (A_b) per host and fit classification
plotC <- ssp_intestine %>%
  mutate(fit = ifelse(fit == "good", "Good damage fit", "Bad damage fit")) %>%
  group_by(host) %>%
  mutate(total_abundance = sum(abundance)) %>%
  ungroup() %>%
  mutate(
    host = fct_reorder(host, -total_abundance),
    host = fct_relevel(host, "Other", after = Inf)
  ) %>%
  ggplot(aes(x = y_bp, y = A_b, shape = paste(domain, phylum, class, order, family), color = host)) +
	geom_point() +
	scale_color_manual(values = colors_host) +
	scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15)) +
	coord_flip() +
	scale_x_continuous(breaks = seq(200, 2000, 100), limits = c(200, 2020)) +
	facet_nested(~fit) +
	guides(shape = "none") +
	theme(
		aspect.ratio = 2 / 1,
		legend.position = "none",
		axis.text.y = element_blank(),
		axis.ticks.y = element_blank()
	) +
	labs(
		y = "Damage (A_b)",
		x = "",
		title = "Damage by species",
		color = "Reference host",
		shape = "Family"
	)


# ----------------------------------------------------------------------
# (5) Export figure
# ----------------------------------------------------------------------

png(file = "plots/fig_SX7.png", width = 14, height = 9, units = "in", res = 600)
egg::ggarrange(plotA, plotB1, plotB2, plotC,
               nrow = 1,
               widths = c(2, 1.5, 1, 2.5),
               labels = c("A", "B", "C", "D"))
dev.off()



sources_host <- ssp_intestine %>%
	pull(host) %>%
	unique()

n_sp_vir <- ssp_intestine %>%
	filter(abundance > 0) %>%
	select(label, species) %>%
	distinct() %>%
	group_by(label) %>%
	count()

HOST <- ssp_intestine %>% 
    filter(!is.na(abundance)) %>%
    group_by(label, host) %>%
    summarise(abundance = sum(n_reads), .groups = "drop") %>%
    tidyr::complete(.,
        host = unique(host),
        label = unique(label),
        fill = list(abundance = 0)) %>%
    left_join(cdata %>% select(label, y_bp)) %>%
    group_by(host) %>%
    mutate(total_abundance = sum(abundance)) %>%
    ungroup() %>%
	group_by(label) %>%
	mutate(relative_abundance = abundance/sum(abundance)) %>%
	ungroup() %>%
	left_join(n_sp_vir) %>%
	mutate(prop = n * relative_abundance) %>%
	mutate(toplevel1 = "Host-associated taxa") %>%
	mutate(toplevel = ifelse(host == "Homo sapiens", "Human", "Animal")) %>%
    mutate(host = fct_reorder(host, -total_abundance)) %>%
    mutate(host = fct_relevel(host, "Other", after = Inf))



# ----------------------------------------------------------------------
# (5) Integrate host and viral sources
# ----------------------------------------------------------------------
# Combine bacterial/archaeal host data with viral source annotations
# from IMG/VR and functionalDB.

ssp <- readRDS("results/ps-stand-filt-ssp.rds")

vir <- read_tsv("data/viral-sources-imgvr.tsv.gz") %>%
    janitor::clean_names()

ssp_source <- ssp %>%
	speedyseq::psmelt() %>%
	filter(domain == "d__Viruses") %>%
	mutate(uvig = str_replace(subspecies, "S__", "")) %>%
	left_join(vir) %>%
	mutate(ecosystem_classification = ifelse(is.na(ecosystem_classification), "unclassified;unclassified;unclassified;unclassified", ecosystem_classification)) %>%
	mutate(ecosystem_classification = ifelse(ecosystem_classification == ";;;", "unclassified;unclassified;unclassified;unclassified", ecosystem_classification)) %>%
	separate(ecosystem_classification, into = c("class1", "class2", "class3", "class4"), sep = ";") %>%
	mutate(across(c(class1, class2, class3, class4), ~ if_else(is.na(class1), "unclassified", .))) %>%
    mutate(class4 = ifelse(class4 == "Biological phosphorus removal" & class3 == "Nutrient removal", "Phosphorus removal", class4)) %>%
    mutate(class4 = ifelse(class4 == "Shale gas reservoir" & class3 == "Deep subsurface", "Gas reservoir", class4)) %>%
	mutate(class2 = ifelse(class2 == "WWTP", "Wastewater", class2))
    
# ----------------------------------------------------------------------
# (6) Identify major viral sources and aggregate abundance
# ----------------------------------------------------------------------
# Identify viral ecological sources (environmental/engineered)
# and summarize their relative abundance over time.

sources_viruses <- ssp_source %>%
	group_by(label) %>%
	mutate(relative_abundance = Abundance/sum(Abundance)) %>%
	ungroup() %>%
	group_by(class1, class2, class3, class4) %>%
	summarise(n = n(), n_ab = sum(relative_abundance > 1e-2), abundance = sum(relative_abundance), .groups = "drop") %>%
	filter(n_ab > 25 | class2 == "Wastewater" | grepl("Wetland", class4) | grepl("Intertidal", class4)) %>%
	mutate(labels = paste(class3, class4, sep = "; ")) %>%
	select(labels) %>%
	distinct() %>%
    filter(labels !=  "Freshwater; Pond") %>%
    filter(labels !=  "Marine; Wetlands") %>%
	pull(labels)

n_sp <- ssp_source %>%
	filter(Abundance > 0) %>%
	filter(domain == "d__Viruses") %>%
	select(label, y_bp, phylum, class, order, family, genus, species) %>%
	distinct() %>%
	group_by(label, y_bp) %>%
	summarise(n = n()) %>%
	ungroup()


VIR <- ssp_source %>%
	
	mutate(across(c(class1, class2, class3, class4), ~ if_else(class2 == "Aquatic" & !class3 %in% c("Freshwater", "Marine", "Deep subsurface"), "other", .))) %>%
	mutate(across(c(class1, class2, class3, class4), ~ if_else(class2 == "Terrestrial" & !class3 %in% c("Soil"), "other", .))) %>%
	mutate(across(c(class1, class2, class3, class4), ~ if_else(class1 == "Engineered" & !class3 %in% c("Industrial wastewater", "Sewage", "Nutrient removal", "Wastewater"), "other", .))) %>%
	
	group_by(label, y_bp, class1, class2, class3, class4) %>%
	summarise(abundance = sum(Abundance), .groups = "drop") %>%
	
	group_by(label) %>%
	mutate(relative_abundance = abundance/sum(abundance)) %>%
	ungroup() %>%

	group_by(class3) %>%
	mutate(sum_relative_abundance = sum(relative_abundance)) %>%
	ungroup() %>%
	
	left_join(n_sp) %>%
	mutate(prop = n * relative_abundance) %>%

	mutate(class3_mod = str_replace(class3, " ", "\n")) %>%
	mutate(class3_mod = fct_reorder(class3_mod, -sum_relative_abundance)) %>%
	
	mutate(path = paste(class3, class4, sep = "; ")) %>%
	mutate(path = fct_reorder(path, as.numeric(factor(paste(class2, class3))))) %>%
	filter(class1 %in% c("Environmental", "Engineered")) %>%
    mutate(class4 = ifelse(class4 == "Shale gas reservoir" & class3 == "Deep subsurface", "Gas reservoir", class4)) %>%
    mutate(class4 = ifelse(class4 == "Biological phosphorus removal" & class3 == "Nutrient removal", "Phosphorus removal", class4))


# ----------------------------------------------------------------------
# (7) Load and process SourceTracker results
# ----------------------------------------------------------------------
# Integrate bacterial and archaeal SourceTracker outputs to identify
# environmental source contributions across time.


source_cdata <- read_tsv("data/sourcetracker/med-biomes-download.txt") %>%
    mutate(biome_subclass = gsub(":|-| +", "_", biome_subclass)) %>%
    mutate(biome_subclass = tolower(biome_subclass)) %>%
    mutate(biome_subclass = gsub(" +", "_", biome_subclass))

header_line <- readLines("results/results-st/mixing_proportions.txt", n = 1) %>% sub("^#", "", .)
stout <- read.table(file = "results/results-st/mixing_proportions.txt", sep = "\t", header = FALSE, col.names = unlist(strsplit(header_line, "\t"))) %>% 
    janitor::clean_names() %>%
    rename(label = sample_id) %>%
    as_tibble() %>%
    pivot_longer(cols = -label, names_to = "source", values_to = "proportion") %>%
    left_join(source_cdata %>% select(biome_class, biome_subclass) %>% distinct(), by = c("source"="biome_subclass")) %>%
    mutate(biome_class = ifelse(source == "unknown", "unknown", biome_class))

stout <- stout %>%
    mutate(source_clean = case_when(
        grepl("root_environmental_aquatic_freshwater_sediment", source) ~ "Freshwater sediment",
        grepl("root_environmental_aquatic_freshwater_lake_atlantic_maritime", source) ~ "Freshwater aquatic; Lake",
        grepl("root_environmental_terrestrial_soil_grasslands", source) ~ "Soil; Grasslands",
        grepl("root_environmental_aquatic_marine_sediment", source) ~ "Marine sediment",
        grepl("root_environmental_terrestrial_soil_agricultural", source) ~ "Soil; Agricultural",
        grepl("root_environmental_terrestrial_soil_permafrost", source) ~ "Soil; Permafrost",
        grepl("root_environmental_aquatic_freshwater_lake_boreal_plains", source) ~ "Freshwater aquatic; Lake boreal",
        grepl("root_environmental_aquatic_freshwater_lake_prairies", source) ~ "Freshwater aquatic; Lake prairies",
        grepl("root_environmental_aquatic_freshwater_lake_montane_cordillera", source) ~ "Freshwater aquatic; Lake montane",
        grepl("root_environmental_aquatic_freshwater_lake_pacific_maritime", source) ~ "Freshwater aquatic; Lake pacific",
        grepl("root_environmental_aquatic_freshwater", source) ~ "Freshwater",
        grepl("root_environmental_aquatic_freshwater_lake_atlantic_highlands", source) ~ "Freshwater aquatic; Lake atlantic",
        grepl("root_environmental_aquatic_freshwater_lake_mixewood_plains", source) ~ "Freshwater aquatic; Lake mixewood",
        grepl("root_environmental_aquatic_marine_polar", source) ~ "Marine aquatic; Polar",
        grepl("root_environmental_aquatic_freshwater_lake_semi_arid_plateau", source) ~ "Freshwater aquatic; Lake semi arid",
        grepl("root_environmental_aquatic_freshwater_lake_boreal_shield", source) ~ "Freshwater aquatic; Lake boreal",
        grepl("root_host_associated_plants", source) ~ "Plants; Unclassified",
        grepl("root_environmental_aquatic_marine_coastal", source) ~ "Marine aquatic; Coastal",
        grepl("root_host_associated_mammals_digestive", source) ~ "Mammals; Digestive",
        grepl("root_environmental_aquatic_estuary", source) ~ "Marine aquatic; estuary",
        grepl("root_environmental_aquatic_marine_mediterranean_sea", source) ~ "Marine aquatic; Mediterranean",
        grepl("root_environmental_aquatic_marine_brackish_water", source) ~ "Marine aquatic; Brackish",
      TRUE ~ "Other"
    ))


st_source <- c("Freshwater sediment", "Marine sediment", "Freshwater aquatic; Lake", "Marine aquatic; Polar", "Soil; Grasslands", "Soil; Agricultural", "Mammals; Digestive", "Plants; Unclassified")


# ----------------------------------------------------------------------
# (8) Combine and visualize bacterial, archaeal, and viral sources
# ----------------------------------------------------------------------
# Produce final multi-panel plot summarizing:
# A) Viral fraction by isolation source
# B) Bacteria/archaea SourceTracker proportions
# C) Host-associated bacterial and archaeal taxa

# Outputs:
#   - plots/fig_X1.png
#   - plots/fig_X1.pdf
#   - plots/fig_X1_2.pdf
#   - my_plot.RData

common_list <- c(sources_host, sources_viruses, st_source)
common_colors <- paired_genus[1:length(common_list)]
names(common_colors) <- common_list

plotVIR <- VIR %>%
    mutate(path = as.character(path)) %>%
    mutate(path = fct_relevel(path, c("Freshwater; Lake", "Freshwater; River", "Marine; Intertidal zone", "Marine; Epipelagic", 
									"Marine; Coastal", "Deep subsurface; Gas reservoir", "Soil; Coastal area", "Soil; Wetlands",
 									"Wastewater; Unclassified", "Industrial wastewater; Unclassified", "Sewage; Unclassified", 
									"Nutrient removal; Phosphorus removal","Marine; Wetlands", "Deep subsurface; Groundwater", 
									"Freshwater; Creek", "Freshwater; Lentic", "Freshwater; Lotic", "Freshwater; Pond", "Freshwater; Sinkhole", 
									"Marine; Inlet", "Marine; Neritic zone", "Marine; Neritic zone/Coastal water", "Marine; Ocean trench", 
									"Marine; Oceanic", "Marine; Sediment", "Marine; Strait", "Marine; Unclassified", "Soil; Peat"))) %>%
	ggplot(., aes(x = y_bp, y = abundance, fill = path))+
		geom_area(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.2)+
		facet_nested(.~class2*class3_mod, scales = "free")+
		scale_fill_manual(values = common_colors)+
		guides(fill = guide_legend(title = "", nrow = 4))+
		scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), 
			breaks = scales::pretty_breaks(n = 2))+
		scale_x_continuous(breaks = seq(200, 2000, 100))+
		theme(
			legend.position = "bottom", 
			aspect.ratio = 5/1, 
			)+
		labs(title = "Viral fraction by isolation source", y = "Abundance", x = "", fill = "",)+
		coord_flip()

ST <- stout %>%
    left_join(cdata) %>%
    filter(source != "unknown") %>%
    group_by(proportion)  %>%
    mutate(mean_proportion = mean(proportion)) %>%
    ungroup() %>%
    arrange(desc(mean_proportion)) %>%
    separate(biome_class, into = c("root", "class1", "class2", "class3"), sep = ":") %>%
    mutate(class3 = case_when(
        class3 == "Marine" & grepl("sediment", source) ~ "Marine\nsediment",
        class3 == "Marine" & !grepl("sediment", source) ~ "Marine\naquatic",
        class3 == "Freshwater" & !grepl("sediment", source) ~ "Freshwater\naquatic",
        class3 == "Freshwater" & grepl("sediment", source) ~ "Freshwater\nsediment",
        TRUE ~ class3
    )) %>% 
    mutate(source_clean = fct_relevel(source_clean, c("Freshwater aquatic; Lake", "Freshwater sediment", "Marine aquatic; Polar", "Marine sediment", "Soil; Grasslands", "Soil; Agricultural", "Mammals; Digestive", "Plants; Unknown", "Soil; Permafrost", "Freshwater", "Freshwater aquatic; Lake prairies", "Freshwater aquatic; Lake montane", "Freshwater aquatic; Lake boreal", "Freshwater aquatic; Lake pacific", "Marine aquatic; Mediterranean", "Marine aquatic; Coastal", "Marine aquatic; estuary", "Marine aquatic; Brackish"))) %>%
    ggplot(., aes(x = y_bp, y = proportion, fill = source_clean)) +
            geom_area(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.2)+
            labs(x = "Age", y = "Source proportion", fill = "Source")+
            facet_nested(.~class1*class2*class3, scales = "free")+
            guides(fill = guide_legend(title = "", nrow = 4))+
            scale_fill_manual(values = common_colors)+
            scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), 
                breaks = scales::pretty_breaks(n = 2))+
            scale_x_continuous(breaks = seq(200, 2000, 100))+
            theme(
                legend.position = "bottom", 
                aspect.ratio = 5/1, 
                axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                )+
            labs(title = "Bacteria and archaea SourceTracker", y = "Proportion", x = "", fill = "",)+
            coord_flip()

plotHOST <- HOST %>%
	ggplot(., aes(x = y_bp, y = abundance, fill = host))+
        geom_area(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.2)+
        scale_fill_manual(values = common_colors)+
        coord_flip()+
        scale_x_continuous(breaks = seq(200, 2200, 100), limits = c(200, 2020))+
		facet_nested(~toplevel, scales = "free", space = "fixed")+
        labs(y = "Abundance", x = "", title = "Bacterial and archaeal\nbased on reference\ngenome isolation host", fill = "")+
		guides(fill = guide_legend(title = "", ncol = 2))+
        scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()), 
			breaks = scales::pretty_breaks(n = 3))+
        theme(
			aspect.ratio = 5/1, 
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
            legend.position = "bottom")

png(file = "plots/fig_X1.png", width = 19, height = 11, units = "in", res = 600)
egg::ggarrange(plotVIR, ST, plotHOST, nrow = 1, widths = c(4,3,1), labels = c("A", "B", "C"))
dev.off()


# my_plot <- egg::ggarrange(plotVIR, ST, plotHOST, nrow = 1, widths = c(4,3,1), labels = c("A", "B", "C"))
# ggsave("plots/fig_X1.pdf", plot = my_plot, device = "pdf", width = 21, height = 9.5)
# ggsave("plots/fig_X1_2.pdf", plot = my_plot, device = cairo_pdf, width = 21, height = 9.5)

# cairo_pdf("plots/fig_X1_2.pdf", width = 21, height = 9.5, useDingbats = FALSE)
# print(my_plot)
# dev.off()


# # Save it as an RData file
# save(my_plot, file = "my_plot.RData")


# ----------------------------------------------------------------------
# End of script
# ----------------------------------------------------------------------
