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
source("/projects/fernandezguerra/people/ngm902/scripts/r-miscellaneous.R")

# tax_data_dmg_sp <- read_tsv("results/aggregated/dmg-summary-ssp.tsv.gz")

# png(file = "presentation0.png", width = 8, height = 10, units = "in", res = 600)
# tax_data_dmg_sp %>%
#     group_by(y_bp) %>%
#     summarise(median_A_b = median(A_b), sd = sd(A_b, na.rm = TRUE), n = n(), .groups = "drop") %>%
#     ggplot(., aes(x = y_bp, y = median_A_b))+
#         geom_errorbar(aes(ymin = median_A_b - sd, ymax = median_A_b + sd), color = "darkgrey", width = 10)+
#         geom_point()+
#         # geom_boxplot(aes(group = y_bp), outlier.size = 0.1)+
#         coord_flip()+
#         labs(x = "Year (CE)", y = "Median A_b")
# dev.off()



# Load labels data and taxonomy
cdata <- read.table(file = "data/metadata.txt", sep = "\t", header = T)

tax_data_dmg_sp <- read_tsv("./results/aggregated/dmg-summary-sp.tsv.gz")

# tax_data_dmg_sum <- tax_data_dmg_sp %>%
#     group_by(species) %>%
#     summarise(A_b = max(A_b), .groups = "drop")


# # Initial (derep) read number by sample
stats_derep <- read.table(file = "data/aggregated/all.stats-derep-summary.tsv.gz", sep = "\t", header = T) %>%
    mutate(step = "derep") %>%
    select(label, num_seqs, step)

stats_initial <- read.table(file = "data/aggregated/all.stats-initial-summary.tsv.gz", sep = "\t", header = T) %>%
    mutate(step = "initial") %>%
    select(label, num_seqs, step)

plotA <- rbind(stats_derep, stats_initial) %>%
    left_join(cdata) %>%
    ggplot(., aes(x = y_bp, y = num_seqs, color = step))+
        geom_point()+
        geom_line()+
        coord_flip()+
        scale_x_continuous(breaks = seq(200, 2000, 100))+
        scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
        labs(x = "Year (CE)", y = "Number of reads", title = "Nº of reads", color = "Reads")+
        theme(legend.position = "bottom")

plotB <- readRDS("./results/aggregated/ps-unfilt-sp.rds") %>%
    speedyseq::psmelt() %>%
    filter(Abundance > 0) %>%
    group_by(label) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    left_join(cdata %>% select(label, y_bp)) %>%
    ggplot(., aes(x = y_bp, y = n))+
        geom_area(alpha = 0.7)+
        geom_bar(stat = "identity")+
        coord_flip()+
        scale_x_continuous(breaks = seq(200, 2000, 100))+
        labs(x = "", y = "Number of species", title = "Nº of species")+
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

plotC <- readRDS("./results/aggregated/ps-unfilt-sp.rds") %>%
    speedyseq::psmelt() %>%
    filter(Abundance > 0) %>%
    group_by(label, domain) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    left_join(cdata %>% select(label, y_bp)) %>%
    mutate(domain = fct_relevel(domain, "d__Bacteria", "d__Archaea", "d__Viruses")) %>%
    ggplot(., aes(x = y_bp, y = n, fill = domain))+
        geom_area(stat = "identity", alpha = 0.7)+
        geom_bar(stat = "identity")+
        scale_fill_manual(values= c("d__Viruses"="#f3c96d", "d__Bacteria"="#74afb9", "d__Archaea"="#9f443d"))+
        facet_nested(.~domain, scales = "free")+
        coord_flip()+
        scale_x_continuous(breaks = seq(200, 2000, 100))+
        scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
        labs(x = "", y = "Number of species", title = "Nº of species by domain")+
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")

plotD <- readRDS("./results/aggregated/ps-unfilt-sp.rds") %>%
    speedyseq::psmelt() %>%
    rename(Abundance = Abundance) %>%
    group_by(label, domain) %>%
    summarise(Abundance = sum(Abundance)) %>%
    ungroup() %>%
    left_join(cdata %>% select(label, y_bp)) %>%
    mutate(domain = fct_relevel(domain, "d__Bacteria", "d__Archaea", "d__Viruses")) %>%
    ggplot(., aes(x = y_bp, y = Abundance, fill = domain))+
        geom_area(stat = "identity", position = "stack", alpha = 0.8)+
        geom_bar(stat = "identity", position = "stack")+
        scale_fill_manual(values= c("d__Viruses"="#f3c96d", "d__Bacteria"="#74afb9", "d__Archaea"="#9f443d"))+
        labs(x = "", y = "Total abundance", title = "Domain abundance", fill = "Domain")+
        scale_x_continuous(breaks = seq(200, 2000, 100))+
        scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
        coord_flip()+
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")

n_sp <- readRDS("./results/aggregated/ps-unfilt-sp.rds") %>%
    speedyseq::psmelt() %>%
	filter(Abundance > 0) %>%
	select(label, y_bp, domain, phylum, class, order, family, genus, species) %>%
	distinct() %>%
	group_by(label) %>%
	summarise(n = n()) %>%
	ungroup()

plotE <- readRDS("./results/aggregated/ps-unfilt-sp.rds") %>%
    speedyseq::psmelt() %>%
    rename(Abundance = Abundance) %>%
    group_by(label, domain) %>%
    summarise(Abundance = sum(Abundance)) %>%
    ungroup() %>%
    group_by(label) %>%
    mutate(proportion = Abundance/sum(Abundance)) %>%
    ungroup() %>%
    left_join(n_sp) %>%
    mutate(relative_abundance = n*proportion) %>%
    left_join(cdata %>% select(label, y_bp)) %>%
    mutate(domain = fct_relevel(domain, "d__Bacteria", "d__Archaea", "d__Viruses")) %>%
    ggplot(., aes(x = y_bp, y = relative_abundance, fill = domain))+
        geom_area(stat = "identity", position = "stack", alpha = 0.8)+
        geom_bar(stat = "identity", position = "stack")+
        scale_fill_manual(values= c("d__Viruses"="#f3c96d", "d__Bacteria"="#74afb9", "d__Archaea"="#9f443d"))+
        labs(x = "", y = "Number of species", title = "Abundance weighted\nnumber of species", fill = "")+
        scale_x_continuous(breaks = seq(200, 2000, 100))+
        scale_y_continuous()+
        coord_flip()+
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")

png(file = "plots/fig_SX1.png", width = 15, height = 9, units = "in", res = 600)
egg::ggarrange(plotA, plotB, plotC, plotD, plotE, nrow = 1, widths = c(1,0.6,1.5,1.5,1.5), labels = c("A", "B", "C", "D", "E"))
dev.off()








to_remove <- tax_data_dmg_sp %>%
	filter(domain != "d__Viruses") %>%
	group_by(species) %>%
	mutate(
		n = n(),
		max_A_b = max(A_b)) %>%
	ungroup() %>%
	filter(n > 30 & max_A_b < 0.05) %>%
	pull(species) %>% unique()


# ps <- readRDS("./results/aggregated/taxonomy/ps-unfilt-sp_glom.rds") %>%
#     subset_taxa(., domain != "d__Viruses") %>%
#     get_st(., nsites = 0.05, vcoeff = 1.5, mean_prop_thresh = 1e-04, norm = "none", glom = NULL, norm_out = FALSE) %>%
#     prune_taxa(!(taxa_names(.) %in% to_remove), .) 


phylo_filt <- readRDS("results/aggregated/ps-unfilt-sp.rds") %>%
    subset_taxa(., domain != "d__Viruses") %>%
    get_st(., nsites = 0.05, vcoeff = 1.5, mean_prop_thresh = 1e-04, norm = "median", glom = NULL, norm_out = FALSE) %>%
    prune_taxa(!(taxa_names(.) %in% to_remove), .)  # Remove specified taxa

saveRDS(phylo_filt, "./results/aggregated/filtered_ps_nmds.rds")

alpha_div <- estimate_richness(phylo_filt, measures = c("Shannon", "Simpson", "Observed"))


plot <- alpha_div %>%
    rownames_to_column("label") %>%
    pivot_longer(cols = -label, names_to = "index", values_to = "value") %>%
    left_join(cdata %>% select(label, y_bp)) %>%
    ggplot(., aes(x = y_bp, y = value))+
        geom_point()+
        geom_line()+
        scale_x_continuous(breaks = seq(200, 2000, 100))+
        facet_nested(.~index, scales = "free_x")+
        coord_flip()

png(file = "plots/fig_SX2.png", width = 10, height = 10, units = "in", res = 600)
plot
dev.off()




n_sp <- phylo_filt %>%
    speedyseq::psmelt() %>%
	filter(Abundance > 0) %>%
	select(label, y_bp, domain, phylum, class, order, family, genus, species) %>%
	distinct() %>%
	group_by(label) %>%
	summarise(n = n()) %>%
	ungroup()

# Get distances
dist <- phylo_filt %>%
    otu_table() %>%
    t() %>%
    as.data.frame() %>%
    vegan::decostand(., method = "hellinger") %>%
    vegan::vegdist(method = "bray")

# Permanova amb el dereplicated-reads
nreads <- dist %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames() %>%
    as.data.frame() %>%
    rename(label = 1) %>%
    left_join(stats_derep) %>%
    mutate(num_seqs = ifelse(is.na(num_seqs), mean(num_seqs, na.rm = TRUE), num_seqs))

permanova <- adonis2(dist ~ nreads$num_seqs, permutations = 999, dist = "bray")
permanova


# Get nmds points
set.seed(123456789)
nmdspoints <- dist %>%
    vegan::metaMDS(try = 900, k = 2) %>%
    {stress_value <<- .$stress; .} %>%
    .$points %>%
    data.frame() %>% 
    rownames_to_column("label") %>%
    left_join(cdata)
cat(paste0("stress value: ", round(stress_value, 3)))

# Get clustering information
k_values <- 2:30
avg_sil_values <- sapply(k_values, function(k) avg_silhouette(hclust(sqrt(dist), method = "ward.D2"), k, dist))
optimal_k <- k_values[which.max(avg_sil_values)]
optimal_k

plot(k_values, avg_sil_values, xlab="k",ylab="av. silhouette",type="b", pch=19)

clusters_opt <- cutree(hclust(sqrt(dist), method = "ward.D2"), k = optimal_k) %>%
    as.data.frame() %>% 
    rownames_to_column("label") %>% 
    rename_at(vars(2), ~ c("cluster_opt")) 

clusters_opt <- clusters_opt %>%
    left_join(cdata %>% select(label, y_bp)) %>%

    group_by(cluster_opt) %>%
    mutate(early_year_cluster_opt = min(y_bp)) %>%
    ungroup() %>%
    
    
    arrange(early_year_cluster_opt) %>% 
    mutate(cluster_order = dense_rank(early_year_cluster_opt)) %>%
    
    select(label, cluster_order) %>%
    rename(cluster_opt = cluster_order)

nmdspoints_clus <- nmdspoints %>%
    left_join(clusters_opt)



# Set cluster colors to be consistent across plots
clusts <- nmdspoints_clus %>% pull(cluster_opt) %>% unique()
clust_cols <- paired_genus[1:length(clusts)]
names(clust_cols) <- clusts




hc <- hclust(sqrt(dist), method = "ward.D2")
dend <- as.dendrogram(hc)
labels(dend) <- labels(dend) %>% data.frame() %>% rename(label = 1) %>% left_join(nmdspoints_clus %>% mutate(y_bp = round(y_bp, 0))) %>% pull(y_bp)

labels_original <- labels(dend)  # Extrae los labels originales
years <- as.numeric(labels_original)  # Convertir a valores numéricos

# Ordenar los labels de izquierda a derecha por el valor medio de los años en cada rama
order_new <- order(years)  # Ordenar en orden creciente

# Aplicar la rotación en el dendrograma
dend <- dendextend::rotate(dend, order = order_new)

Dendplot <- dend %>%
   set("branches_k_color", k=10, value = clust_cols[clust_cols %>% names() %>% as.numeric()]) %>%
   set("branches_lwd", 1.7) %>%
   set("leaves_pch", 21) %>% 
   set("leaves_cex", 1) %>%
   set("leaves_col", "white") %>%
   set("leaves_bg", "white") 



plotA <- nmdspoints_clus %>%
    left_join(stats_derep %>% select(label, num_seqs)) %>%
    mutate(num_seqs = ifelse(label == "CGG3017468", 1e8, num_seqs)) %>% # parche temporal
    ggplot(., aes(MDS1, MDS2))+
        geom_point(aes(size = num_seqs, fill = as.factor(cluster_opt)), shape = 21)+
        geom_text_repel(aes(label = round(y_bp, 0)), size = 3)+
        scale_fill_manual(values = clust_cols)+
        guides(fill = "none")+
        scale_size_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
        annotate("text", x = Inf, y = Inf, label = paste("stress:", round(stress_value, 3)), hjust = 1.2, vjust = 1.8, size = 4, color = "black")+
        labs(title = "Non-metric MultiDimensional Scaling (NMDS) samples representation",
            size = "Derreplicated\nreads",
            fill = "Date (CE)",
            shape = "Extraction and flowcell")+
        theme(aspect.ratio = 2/2.5, legend.position = "bottom")

plotB <- nmdspoints_clus %>%
    ggplot(., aes(y = y_bp))+
        geom_point(aes(x = MDS1, color = as.factor(cluster_opt)))+
        scale_color_manual(values = clust_cols)+
        geom_segment(data = . %>% group_by(cluster_opt) %>% summarise(y_min = min(y_bp), y_max = max(y_bp), .groups = "drop"), 
                     aes(x = 2.4, xend = 2.4, y = y_min, yend = y_max, color = as.factor(cluster_opt)),
                     linewidth = 1, arrow = arrow(length = unit(0.05, "inches"), ends = "both", type = "open"))+
        geom_text(data = . %>% group_by(cluster_opt) %>% summarise(mval = median(y_bp, .groups = "drop")), aes(x = 2.9, y = mval, label = cluster_opt), size = 3)+
        geom_segment(data = . %>% group_by(episode) %>% summarise(y_min = min(y_bp), y_max = max(y_bp), .groups = "drop"), 
                     aes(x = 3.5, xend = 3.5, y = y_min, yend = y_max, color = as.factor(episode)),
                     linewidth = 1, arrow = arrow(length = unit(0.05, "inches"), ends = "both", type = "open"))+
        geom_text(data = . %>% group_by(episode) %>% summarise(mval = mean(y_bp, .groups = "drop")), aes(x = 4, y = mval, label = episode))+
        labs(color = "Cluster\ngroup", y = "Age (CE)", title = "MDS1 vs age with groups and episodes")+
        xlim(c(-2, 4.3))+
        theme(aspect.ratio = 2/1.2, legend.position = "bottom")

png(file = "plots/fig_SX3.png", width = 12, height = 7, units = "in", res = 600)
egg::ggarrange(plotA, plotB, nrow = 1, widths = c(2, 1), labels = c("A", "B"))
dev.off()

png(file = "plots/fig_SX4.png", width = 12, height = 4.5, units = "in", res = 600)
plot(Dendplot, xlab = "Samples (year)", ylab = "Height")
title("C    Hierarchical clustering dendrogram of samples based on sqrt-transformed Bray-Curtis dissimilarity", adj = 0)
dev.off()







#  Indicator species analysis
data_ind <- phylo_filt %>%
    # microbiome::transform("compositional") %>%
    otu_table() %>%
    as.data.frame() %>%
    t()  %>%
    as.data.frame() %>%
    rownames_to_column(var = "label") %>% 
    full_join(clusters_opt) %>%
    column_to_rownames("label")

data_ind_obs <- data_ind %>% select(-cluster_opt)

indval <- multipatt(data_ind_obs, data_ind$cluster_opt, func = "IndVal.g", duleg = TRUE)

p_value_threshold <- 0.01
indva_stat <- 0.5

indv_group <- indval$sign %>%
    as.data.frame() %>% 
    filter(p.value <= p_value_threshold & stat > indva_stat) %>%
    rownames_to_column(var = "taxa") %>%
    select(taxa, index, stat) %>%
    mutate(index = paste0("Cluster ", index))%>%
    rename(species = taxa) 

plotA <- phylo_filt %>%
    speedyseq::psmelt() %>%
    as_tibble() %>%
    rename(taxa = OTU, n_reads = Abundance) %>%
    group_by(label) %>%
    mutate(relative_abundance = n_reads/sum(n_reads)) %>%
    ungroup() %>% 
    left_join(indv_group) %>% 
    mutate(index = ifelse(is.na(index), "No indval", index)) %>% 
    group_by(label, index) %>%
    summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
    left_join(n_sp) %>%
    mutate(relative_abundance = n*relative_abundance) %>%
    left_join(cdata %>% select(label, y_bp)) %>%
    left_join(clusters_opt) %>%
    mutate(y_indval = ifelse(index == paste0("Cluster ", cluster_opt), y_bp, NA)) %>%
    mutate(index = fct_relevel(index, "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10", "No indval")) %>%
    ggplot(., aes(x = y_bp, y = relative_abundance, fill = index))+
        geom_bar(stat = "identity", position = "stack", alpha = 1, width = 1)+
        geom_area(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.2)+
        scale_fill_manual(values = c(unname(clust_cols[clust_cols %>% names() %>% as.numeric()]), "grey"))+
        guides(fill = guide_legend(ncol = 3))+
        labs(x = "Year (CE)", title = "Abundance weighted\nnumber of species", y = "Number of species", fill = "")+
        scale_x_continuous(breaks = seq(200, 2000, 100))+
        coord_flip()+
        theme(legend.position = "bottom")

plotB <- phylo_filt %>%
    speedyseq::psmelt() %>%
    as_tibble() %>%
    rename(taxa = OTU, n_reads = Abundance) %>%
    group_by(label) %>%
    mutate(relative_abundance = n_reads/sum(n_reads)) %>%
    ungroup() %>%
    left_join(indv_group) %>% 
    mutate(index = ifelse(is.na(index), "No indval", index)) %>% 
    group_by(across(-c(n_reads, relative_abundance, class, order, family, genus, species, taxa, stat))) %>%
    summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>% 
    group_by(phylum) %>%
    mutate(mean_relative_abundance = mean(relative_abundance)) %>%
    ungroup() %>%
    mutate(across(c(phylum),
            ~ if_else(mean_relative_abundance < quantile(mean_relative_abundance, 0.5), "other", .))) %>%
    ungroup() %>% 
    select(-mean_relative_abundance) %>%   
    group_by(across(-c(relative_abundance))) %>%
    summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
    left_join(clusters_opt) %>%
    mutate(y_indval = ifelse(index == paste0("Cluster ", cluster_opt), y_bp, NA)) %>%
    mutate(index = fct_relevel(index, "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10", "No indval")) %>%
    mutate(episode = case_when(
        index == "Cluster 1" ~ "Episode\n~I",
        index == "Cluster 2" ~ "Episode\n~II",
        index == "Cluster 3" ~ "Episode\n~III",
        index %in% c("Cluster 4", "Cluster 5", "Cluster 6") ~ "Episode\n~IV",
        index == "Cluster 7" ~ "Episode\n~V",
        index %in% c("Cluster 8", "Cluster 9", "Cluster 10") ~ "Episode\n~VI",
        TRUE ~ "No indval"
    )) %>% 
    mutate(episode = fct_relevel(episode, "Episode\n~I", "Episode\n~II", "Episode\n~III", "Episode\n~IV", "Episode\n~V", "Episode\n~VI", "No indval")) %>%
    ggplot(., aes(x = y_bp, y = relative_abundance, fill = paste(domain, phylum)))+
        geom_rect(data = . %>% select(label, y_indval, index, episode) %>% distinct(), aes(xmin = y_indval-5, xmax = y_indval+5, ymin = 0, ymax = Inf), 
            fill = "darkgrey", alpha = 0.4, inherit.aes = FALSE)+
        geom_area(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.1)+
        geom_bar(stat = "identity", position = "stack", alpha = 0.9)+
        scale_fill_manual(values = paired_genus[11:50])+
        scale_y_continuous(labels = scales::percent_format(), breaks = scales::pretty_breaks(n = 2))+
        scale_x_continuous(breaks = seq(200, 2000, 100))+
        labs(x = "", y = "Relative abundance", fill = "", title = "Main phyla relative abundance for each cluster group")+
        facet_nested(.~episode*index, scales = "fixed", space = "fixed")+
        theme(legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank())+
        coord_flip()


png(file = "plots/fig_SX5.png", width = 17, height = 10, units = "in", res = 600)
egg::ggarrange(plotA, plotB, nrow = 1, widths = c(1,5), labels = c("A", "B"))
dev.off()




plotNon <- phylo_filt %>%
    speedyseq::psmelt() %>%
    as_tibble() %>%
    rename(taxa = OTU, n_reads = Abundance) %>%
    group_by(label) %>%
    mutate(relative_abundance = n_reads/sum(n_reads)) %>%
    ungroup() %>%
    left_join(indv_group) %>% 
    mutate(index = ifelse(is.na(index), "No indval", index)) %>% 
    group_by(across(-c(n_reads, relative_abundance, class, order, family, genus, species, taxa, stat))) %>%
    summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>% 
    group_by(phylum) %>%
    mutate(mean_relative_abundance = mean(relative_abundance)) %>%
    ungroup() %>%
    mutate(across(c(phylum),
            ~ if_else(mean_relative_abundance < quantile(mean_relative_abundance, 0.5), "other", .))) %>%
    ungroup() %>% 
    select(-mean_relative_abundance) %>%   
    group_by(across(-c(relative_abundance))) %>%
    summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
    left_join(clusters_opt) %>%
    mutate(y_indval = ifelse(index == paste0("Cluster ", cluster_opt), y_bp, NA)) %>%
    mutate(index = fct_relevel(index, "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10", "No indval")) %>%
    mutate(episode = case_when(
        index == "Cluster 1" ~ "Episode\n~I",
        index == "Cluster 2" ~ "Episode\n~II",
        index == "Cluster 3" ~ "Episode\n~III",
        index %in% c("Cluster 4", "Cluster 5", "Cluster 6") ~ "Episode\n~IV",
        index == "Cluster 7" ~ "Episode\n~V",
        index %in% c("Cluster 8", "Cluster 9", "Cluster 10") ~ "Episode\n~VI",
        TRUE ~ "No indval"
    )) %>% 
    mutate(episode = fct_relevel(episode, "Episode\n~I", "Episode\n~II", "Episode\n~III", "Episode\n~IV", "Episode\n~V", "Episode\n~VI", "No indval")) %>%
    filter(episode == "No indval") %>%

    ggplot(., aes(x = y_bp, y = relative_abundance, fill = paste(domain, phylum)))+
        geom_rect(data = . %>% select(label, y_indval, index, episode) %>% distinct(), aes(xmin = y_indval-5, xmax = y_indval+5, ymin = 0, ymax = Inf), 
            fill = "darkgrey", alpha = 0.4, inherit.aes = FALSE)+
        geom_area(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.1)+
        geom_bar(stat = "identity", position = "stack", alpha = 0.9)+
        scale_fill_manual(values = paired_genus[30:50])+
        scale_y_continuous(labels = scales::percent_format(), breaks = scales::pretty_breaks(n = 2))+
        scale_x_continuous(breaks = seq(200, 2000, 100))+
        labs(x = "", y = "Relative abundance", fill = "", title = "Taxa not associated to any cluster of samples")+
        facet_nested(.~domain, scales = "fixed", space = "fixed")+
        theme(
            legend.position = "right", 
            aspect.ratio = 3/1)+
        coord_flip()
    
png(file = "plots/fig_SX6.png", width = 10, height = 10, units = "in", res = 600)
plotNon
dev.off()







tax_paths_unique <- c()

cluster <- 2
for(cluster in 1:10){

    plot_reads <- phylo_filt %>%
        speedyseq::psmelt() %>%
        as_tibble() %>%
        rename(taxa = OTU, n_reads = Abundance) %>%
        mutate(relative_abundance = n_reads) %>%
        group_by(label) %>%
        mutate(relative_abundance = n_reads/sum(n_reads)) %>%
        ungroup() %>%
        left_join(indv_group) %>%
        filter(index == paste0("Cluster ", cluster)) %>%        
        filter(!is.na(index)) %>%
        group_by(across(-c(n_reads, relative_abundance, species, taxa, stat))) %>%
        summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
        group_by(genus) %>%
        mutate(mean_relative_abundance = mean(relative_abundance)) %>%
        ungroup()
    
    phylo_filt %>%
        speedyseq::psmelt() %>%
        as_tibble() %>%
        rename(taxa = OTU, n_reads = Abundance) %>%
        group_by(label) %>%
        mutate(relative_abundance = n_reads/sum(n_reads)) %>%
        ungroup() %>%
        left_join(indv_group) %>%
        filter(index == paste0("Cluster ", cluster)) %>%
        group_by(domain, phylum, order, family, genus, species) %>%
        summarise(mean_relative_abundance = mean(relative_abundance), .groups = "drop") %>%
        arrange(desc(mean_relative_abundance)) %>%
        # slice(1:30) %>%
        {species_number <<- nrow(.); .} %>%
        {species_list <<- .; .} %>%
        knitr::kable()
    
    phylo_filt %>%
        speedyseq::psmelt() %>%
        as_tibble() %>%
        rename(taxa = OTU, n_reads = Abundance) %>%
        group_by(label) %>%
        mutate(relative_abundance = n_reads/sum(n_reads)) %>%
        ungroup() %>%
        left_join(indv_group) %>%
        filter(index == paste0("Cluster ", cluster)) %>%
        group_by(domain, phylum, order, family, genus, species) %>%
        summarise(mean_relative_abundance = mean(relative_abundance), .groups = "drop") %>%
        arrange(desc(mean_relative_abundance)) %>%
        slice(1:20) %>%
        mutate(species = fct_reorder(species, mean_relative_abundance)) %>%
        ggplot(., aes(x = species, y = mean_relative_abundance))+
            geom_point()+
            coord_flip()
        
    cat(paste0(species_number, " species identified as Indicator species"))

    kept_family <- plot_reads %>% select(genus, mean_relative_abundance) %>% distinct() %>% arrange(desc(mean_relative_abundance)) %>% slice(1:15) %>% pull(genus)

    plot_reads_filt <- plot_reads %>%
        mutate(
            across(c(phylum, class, order, family, genus),
                ~ if_else(genus %in% kept_family, ., "other"))
        ) %>%
        mutate(
            across(c(class, order, family, genus),
                ~ if_else(genus == "other", "",.))
        ) %>%
        ungroup() %>% 
        select(-mean_relative_abundance) %>%   
        group_by(across(-c(relative_abundance))) %>%
        summarise(relative_abundance = sum(relative_abundance), .groups = "drop") %>%
        left_join(clusters_opt) %>%
        mutate(y_indval = ifelse(index == paste0("Group_", cluster_opt), y_bp, NA))

    # Set color palette
    new_tax <- plot_reads_filt %>% mutate(taxa = paste(domain, phylum, class, order, family, genus)) %>% pull(taxa) %>% unique()

    for(tx in new_tax){
        if(!tx %in% tax_paths_unique){
            tax_paths_unique <- c(tax_paths_unique, tx)
        } else{NULL}
    }
    
    colors_genus <- paired_genus[!grepl("gray", paired_genus)]
    colors_genus <- colors_genus[1:length(tax_paths_unique)]
    names(colors_genus) <- tax_paths_unique
    colors_genus[grepl("other", names(colors_genus))] <- c("#dcdcdc", "#c9c9c9")

    plotA <- ggplot(plot_reads_filt, aes(x = y_bp, y = relative_abundance, fill = paste(domain, phylum, class, order, family, genus)))+
                geom_rect(data = plot_reads_filt %>% select(label, y_indval, index) %>% distinct(), aes(xmin = y_indval-5, xmax = y_indval+5, ymin = 0, ymax = Inf), 
                fill = "#1F78B4", alpha = 0.2, inherit.aes = FALSE)+
                geom_bar(stat = "identity", position = "stack", alpha = 1, width = 1)+
                geom_area(stat = "identity", position = "stack", alpha = 0.7, color = "black", linewidth = 0.2)+
                scale_fill_manual(values = colors_genus)+
                guides(fill = guide_legend(ncol = 1))+
                scale_y_continuous(labels = scales::percent_format())+
                scale_x_continuous(breaks = seq(200, 2000, 100))+
                labs(x = "Age", y = "Relative abundance", fill = "", title = paste0("Cluster group ", cluster))+
                theme(legend.position = "none", aspect.ratio = 2/1)+
                coord_flip()

    plotB <- plot_reads_filt %>% 
        filter(index == paste0("Group_", cluster_opt)) %>%
        group_by(across(c(domain, phylum, class, order, family, genus, index))) %>%
        summarise(relative_abundance = mean(relative_abundance), .groups = "drop") %>%
        ggplot(., aes(x = index, y = relative_abundance, fill = paste(domain, phylum, class, order, family, genus)))+
            geom_bar(stat = "identity", position = "stack", alpha = 0.7, color = "black", width = 1, linewidth = 0.2)+
            scale_fill_manual(values = colors_genus)+
            scale_y_continuous(labels = scales::percent_format())+
            guides(fill = guide_legend(ncol = 1))+
            labs(x = "", y = "Relative abundance", fill = "Taxonomic path")+
            theme(legend.position = "bottom")+
            facet_nested(.~domain, scales = "free_y")+
            theme(legend.position = "right", strip.text.y = element_text(size = 8))
    
    plotC <- phylo_filt %>%
            speedyseq::psmelt() %>%
            as_tibble() %>%
            rename(taxa = OTU, n_reads = Abundance) %>%
            mutate(relative_abundance = n_reads) %>%
            group_by(label) %>%
            mutate(relative_abundance = n_reads/sum(n_reads)) %>%
            ungroup() %>%
            left_join(indv_group) %>%
            filter(index == paste0("Group_", cluster)) %>%        
            filter(!is.na(index)) %>%
            left_join(clusters_opt) %>%
            filter(index == paste0("Group_", cluster_opt)) %>%
            group_by(domain, phylum, class, order, family, genus, species, index) %>%
            summarise(relative_abundance = mean(relative_abundance), .groups = "drop") %>%
            left_join(indv_group) %>%
            left_join(tax_data_dmg_sum, by = "species") %>%
            mutate(color_name = ifelse(paste(domain, phylum, class, order, family, genus) %in% tax_paths_unique, paste(domain, phylum, class, order, family, genus), NA)) %>%
            mutate(color_name = ifelse(is.na(color_name) , paste0(domain, " other    "), color_name)) %>%
            ggplot(., aes(x = stat, y = relative_abundance, size = A_b))+
                geom_point(aes(color = color_name), alpha = 1)+
                guides(color = "none")+
                geom_text_repel(data = . %>% filter(!grepl("other", color_name)), aes(label = species), size = 3)+
                scale_color_manual(values = colors_genus)+
                scale_size_continuous(limits = c(0, 0.3))+
                scale_y_sqrt(labels = scales::percent_format())+
                labs(x = "Inval statistic (specificity to group)", y = "Mean relative abundance in group samples", color = "Taxonomic path", size = "Maximum A_b", title = "Relative abundance and Indval statistic")+
                theme(legend.position = "right", aspect.ratio = 1/1.9)


    png(file = paste0("plots/aggregated/figures/indval_clusters/cluster_", cluster, "_genus.png"), width = 17, height = 12, units = "in", res = 600)
    ggpubr::ggarrange(egg::ggarrange(plotA, plotB, nrow = 1, widths = c(1,1)), ggpubr::ggarrange(plotC, NULL, widths = c(1.4,1)), ncol = 1, heights = c(1.6, 1))
    dev.off()


}

