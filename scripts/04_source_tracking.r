library(tidyverse)
library(knitr)
library(readxl)
library(measurements)
library(ggmap)
library(ggrepel)
library(lvplot)
library(cowplot)
library(phyloseq)
library(gridExtra)
library(dendextend)
library(ggdendro)
library(ggh4x)
library(indicspecies)
library(vegan)

setwd(dir = "/projects/caeg/people/ngm902/apps/repos/tjornin-microbial-analysis")
source("/projects/caeg/people/ngm902/scripts/r-miscellaneous.R")

ssp <- read_tsv("./results/aggregated/dmg-summary-ssp.tsv.gz")
cdata <- read.table(file = "data/metadata.txt", sep = "\t", header = T)

ssp_intestine <- ssp %>%
    filter(!is.na(abundance)) %>%
    mutate(host = case_when(
            # Atopobiaceae
            subspecies == "S__GCA_902793195.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_017510875.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902758085.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902780305.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902764955.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902793465.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_002449555.1" ~ "Ovis (sheep)",
            subspecies == "S__GCA_017460265.1" ~ "Other", #"Roe deer",
            subspecies == "S__GCA_017941705.1" ~ "Other", #"Water deer",
            subspecies == "S__GCA_902788705.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__3300028048_46  " ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902786015.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_018376165.1" ~ "Homo sapiens",
            subspecies == "S__GCA_900540955.1" ~ "Homo sapiens",
            subspecies == "S__GCA_902772085.1" ~ "Bos taurus (Cattle)",

            # Eggerthellaceae
            subspecies == "S__GCA_017511465.1" ~ "Ovis (sheep)",
            subspecies == "S__GCA_902780685.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902781245.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_015062605.1" ~ "Capra (Goat)",
            subspecies == "S__3300021326_23  " ~ "Capra (Goat)",
            subspecies == "S__GCA_902764575.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902788545.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_905234785.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_015062635.1" ~ "Other", #"Cervus elaphus",
            subspecies == "S__GCA_017434045.1" ~ "Capra (Goat)",

            # Lachnospiraceae
            # subspecies == "S__GCF_900101015.1" ~ "Other/unknown", #"Unknown",
            subspecies == "S__GCA_902777355.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902801415.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902788035.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902801585.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902798065.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902779635.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_017510065.1" ~ "Ovis (sheep)",
            subspecies == "S__GCA_900321785.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902786735.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCF_900111235.1" ~ "Ovis (sheep)",
            subspecies == "S__GCF_003609635.1" ~ "Ovis (sheep)",
            subspecies == "S__GCA_902760435.1" ~ "Bos taurus (Cattle)",
            # subspecies == "S__GCF_002251225.2" ~ "Environmental (aquatic)",
            # subspecies == "S__GCF_003201285.1" ~ "Other/unknown", #"Unknown",

            # Methanobrevibacter
            subspecies == "S__GCA_016291275.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__3300010282_39" ~ "Other", # Capybara
            subspecies == "S__GCA_015062985.1" ~ "Capra (Goat)",
            subspecies == "S__GCA_902756745.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902800455.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_900318035.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902765105.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902764015.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902801725.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_015062905.1" ~ "Capra (Goat)",
            subspecies == "S__GCA_017432245.1" ~ "Capra (Goat)",
            subspecies == "S__GCA_902774685.1" ~ "Bos taurus (Cattle)",
            # subspecies == "S__GCF_900103415.1" ~ "Other/unknown",
            subspecies == "S__GCF_001548675.1" ~ "Bos (other)", # indicus,
            subspecies == "S__GCA_017553445.1" ~ "Other", # Water buffalo
            subspecies == "S__GCA_002509015.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902764425.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902781165.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_900319985.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_002498445.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902771435.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902777885.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902789505.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_900314695.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCF_003111625.1" ~ "Bos (other)", #unkown
            subspecies == "S__GCA_902799005.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_900320955.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902794475.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_015656755.1" ~ "Bos (other)", #Bos grunniens,
            subspecies == "S__GCA_902800515.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_902764455.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCF_001477655.1" ~ "Ovis (sheep)",
            subspecies == "S__GCA_902798335.1" ~ "Bos taurus (Cattle)",
            subspecies == "S__GCA_015062915.1" ~ "Capra (Goat)",
            subspecies == "S__GCF_000404165.1" ~ "Ovis (sheep)",
            subspecies == "S__GCF_002813085.1" ~ "Homo sapiens",
            subspecies == "S__GCF_000016525.1" ~ "Homo sapiens",

            TRUE ~ NA
        )) %>%
        filter(!is.na(host))

colors_host <- paired_genus[1:length(ssp_intestine %>% pull(host) %>% unique())]
names(colors_host) <- ssp_intestine %>% pull(host) %>% unique()
colors_host[grepl("Other", names(colors_host))] <- "#dcdcdc"


plotA <- ssp_intestine %>%
    group_by(host) %>%
    mutate(total_abundance = sum(abundance)) %>%
    ungroup() %>%
    mutate(host = fct_reorder(host, -total_abundance)) %>%
    mutate(host = fct_relevel(host, "Other", after = Inf)) %>%
    ggplot(., aes(x = y_bp, y = abundance, color = host, shape = paste(domain, phylum, class, order, family)))+
        geom_point()+
        scale_color_manual(values = colors_host)+
        coord_flip()+
        scale_x_continuous(breaks = seq(200, 2000, 100), limits = c(200, 2020))+
        guides(shape = guide_legend(ncol = 1), color = guide_legend(ncol = 2))+
        scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
        theme(aspect.ratio = 2/1, legend.position = "bottom", legend.justification = c(0, 0))+
        labs(y = "Abundance", x = "Year (CE)", title = "Abundance per reference", color = "Reference host", shape = "Taxonomic\npath")

plotB1 <- ssp_intestine %>%
    filter(!is.na(abundance)) %>%
    group_by(label, host) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    tidyr::complete(.,
        host = unique(host),
        label = unique(label),
        fill = list(abundance = 0)) %>%
    left_join(cdata %>% select(label, y_bp)) %>%
    group_by(host) %>%
    mutate(total_abundance = sum(abundance)) %>%
    ungroup() %>%
    mutate(host = fct_reorder(host, -total_abundance)) %>%
    mutate(host = fct_relevel(host, "Other", after = Inf)) %>% 
    ggplot(., aes(x = y_bp, y = abundance, fill = host))+
        geom_area(stat = "identity", position = "stack", alpha = 0.7)+
        geom_bar(stat = "identity", position = "stack")+
        scale_fill_manual(values = colors_host)+
        coord_flip()+
        scale_x_continuous(breaks = seq(200, 2000, 100), limits = c(200, 2020))+
        labs(y = "Abundance", x = "", title = "Abundance by host", fill = "")+
        scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
        theme(aspect.ratio = 2/1, 
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
            legend.position = "none")

tax_data_dmg_sp <- read_tsv("./results/aggregated/dmg-summary-sp.tsv.gz")

plotB2 <- tax_data_dmg_sp %>%
    filter(species %in% ssp_intestine$species) %>%
    filter(!is.na(abundance)) %>%
    group_by(label, y_bp, domain, phylum, class, order, family) %>%
    summarise(count = n(), .groups = "drop") %>%
    ggplot(., aes(x = y_bp, y = count, group = paste(domain, phylum, class, order, family)))+
        geom_point(aes(shape = paste(domain, phylum, class, order, family)))+
        geom_line(color = "black")+
        coord_flip()+
        scale_x_continuous(breaks = seq(200, 2000, 100), limits = c(200, 2020))+
        labs(y = "Number of species", x = "", title = "Species\ncount by family", fill = "Damage fit")+
        theme(aspect.ratio = 2/1, axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")

plotC <- ssp_intestine %>%
    mutate(fit = ifelse(fit == "good", "Good damage fit", "Bad damage fit")) %>%
    group_by(host) %>%
    mutate(total_abundance = sum(abundance)) %>%
    ungroup() %>%
    mutate(host = fct_reorder(host, -total_abundance)) %>%
    mutate(host = fct_relevel(host, "Other", after = Inf)) %>%
    ggplot(., aes(x = y_bp, y = A_b, shape = paste(domain, phylum, class, order, family), color = host))+
        geom_point()+
        scale_color_manual(values = colors_host)+
        scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15)) +
        coord_flip()+
        scale_x_continuous(breaks = seq(200, 2000, 100), limits = c(200, 2020))+
        facet_nested(~fit)+
        guides(shape = "none")+
        theme(aspect.ratio = 2/1, legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())+
        labs(y = "Damage (A_b)", x = "", title = "Damage by species", color = "Reference host", shape = "Family")

png(file = "plots/fig_SX7.png", width = 14, height = 9, units = "in", res = 600)
egg::ggarrange(plotA, plotB1, plotB2, plotC, nrow = 1, widths = c(2,1.5,1,2.5), labels = c("A", "B", "C", "D"))
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




########## Viruses sources ##########
ssp <- readRDS("./results/aggregated/ps-stand-filt-ssp.rds")

vir <- read_tsv(file = "/projects/caeg/people/ngm902/apps/functionalDB/IMGVR_all_Sequence_information-high_confidence.tsv.gz") %>%
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
    

# Get the main sources
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
	
	# filter(path %in% names(colors_genus)) %>%
	# filter(class1 != "unclassified") %>%
	# filter(class2 != "other") %>%
	# mutate(path = paste(class2, class3, class4, sep = "; ")) %>%
	mutate(path = paste(class3, class4, sep = "; ")) %>%
	mutate(path = fct_reorder(path, as.numeric(factor(paste(class2, class3))))) %>%
	filter(class1 %in% c("Environmental", "Engineered")) %>%
    # mutate(class1 = fct_relevel(class1, "Environmental", "Engineered")) %>%
    mutate(class4 = ifelse(class4 == "Shale gas reservoir" & class3 == "Deep subsurface", "Gas reservoir", class4)) %>%
    mutate(class4 = ifelse(class4 == "Biological phosphorus removal" & class3 == "Nutrient removal", "Phosphorus removal", class4))


 
    
######### Sourcetracker ##########

source_cdata <- read_tsv("/projects/caeg/people/ngm902/apps/sourcetracker/source_data/med-biomes-download.txt") %>%
    mutate(biome_subclass = gsub(":|-| +", "_", biome_subclass)) %>%
    mutate(biome_subclass = tolower(biome_subclass)) %>%
    mutate(biome_subclass = gsub(" +", "_", biome_subclass))
    
header_line <- readLines("/projects/caeg/people/ngm902/tjornen/results/aggregated/sourcetracker/st-biome-class-gm_v2/mixing_proportions.txt", n = 1) %>% sub("^#", "", .)
stout <- read.table(file = "/projects/caeg/people/ngm902/tjornen/results/aggregated/sourcetracker/st-biome-class-gm_v2/mixing_proportions.txt", sep = "\t", header = FALSE, col.names = unlist(strsplit(header_line, "\t"))) %>% 
# header_line <- readLines("/projects/caeg/people/ngm902/tjornen/results/aggregated/sourcetracker/st-sp_min/mixing_proportions.txt", n = 1) %>% sub("^#", "", .)
# stout <- read.table(file = "/projects/caeg/people/ngm902/tjornen/results/aggregated/sourcetracker/st-sp_min/mixing_proportions.txt", sep = "\t", header = FALSE, col.names = unlist(strsplit(header_line, "\t"))) %>% 
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
stout %>% 
    group_by(source_clean) %>%
    summarise(mean_proportion = mean(proportion),
        median_proportion = median(proportion),
        sum_proportion = sum(proportion),
        max_proportion = max(proportion), .groups = "drop") %>%
    mutate(is_selected = ifelse(source_clean %in% st_source, "Selected", "Other")) %>%
    View()

    
## Create a common color palete
common_list <- c(sources_host, sources_viruses, st_source)
common_colors <- paired_genus[1:length(common_list)]
names(common_colors) <- common_list


plotVIR <- VIR %>%
    mutate(path = as.character(path)) %>%
    mutate(path = fct_relevel(path, c("Freshwater; Lake", "Freshwater; River", "Marine; Intertidal zone", "Marine; Epipelagic", "Marine; Coastal", "Deep subsurface; Gas reservoir", "Soil; Coastal area", "Soil; Wetlands",
 "Wastewater; Unclassified", "Industrial wastewater; Unclassified", "Sewage; Unclassified", "Nutrient removal; Phosphorus removal","Marine; Wetlands",
"Deep subsurface; Groundwater", "Freshwater; Creek", "Freshwater; Lentic", "Freshwater; Lotic", "Freshwater; Pond", "Freshwater; Sinkhole", "Marine; Inlet", 
 "Marine; Neritic zone", "Marine; Neritic zone/Coastal water", "Marine; Ocean trench", "Marine; Oceanic", "Marine; Sediment", "Marine; Strait", 
"Marine; Unclassified", "Soil; Peat"))) %>%
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
        # geom_vline(xintercept = c(877, 1104, 1362, 1477, 1783, 1845, 2014), linetype = 2)+
        theme(
			aspect.ratio = 5/1, 
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
            legend.position = "bottom")

png(file = "plots/fig_X1.png", width = 19, height = 11, units = "in", res = 600)
egg::ggarrange(plotVIR, ST, plotHOST, nrow = 1, widths = c(4,3,1), labels = c("A", "B", "C"))
dev.off()


my_plot <- egg::ggarrange(plotVIR, ST, plotHOST, nrow = 1, widths = c(4,3,1), labels = c("A", "B", "C"))
ggsave("plots/fig_X1.pdf", plot = my_plot, device = "pdf", width = 21, height = 9.5)
ggsave("plots/fig_X1_2.pdf", plot = my_plot, device = cairo_pdf, width = 21, height = 9.5)

cairo_pdf("plots/fig_X1_2.pdf", width = 21, height = 9.5, useDingbats = FALSE)
print(my_plot)
dev.off()


# Save it as an RData file
save(my_plot, file = "my_plot.RData")