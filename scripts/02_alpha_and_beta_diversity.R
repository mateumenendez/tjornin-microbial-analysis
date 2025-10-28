# ======================================================================
# Script: 02_alpha_and_beta_diversity.R
# Author: Mateu Menendez Serra
# Affiliation: University of Copenhagen
# Contact: mateu.menendez@sund.ku.dk
# Date: October 2025
#
# Description:
#   This script performs alpha-diversity analysis, NMDS ordination, 
#   hierarchical clustering based on Bray-Curtis dissimilarity, and 
#   indicator species analysis to detect taxa associated with specific 
#   sample clusters.
#
# Input:
#   - data/metadata.txt
#   - data/all.stats-derep-summary.tsv.gz
#   - data/all.stats-initial-summary.tsv.gz
#   - results/ps-unfilt-sp.rds
#   - results/dmg-summary-sp.tsv.gz
#
# Output:
#   - plots/fig_SX1.png
#   - plots/fig_SX2.png
#   - plots/fig_SX3.png
#   - plots/fig_SX4.png
#   - plots/fig_SX5.png
#   - plots/fig_SX6.png
#   - results/filtered_ps_nmds.rds
#
# Dependencies:
#   tidyverse, phyloseq, vegan, indicspecies, dendextend, ggh4x, ggrepel, egg
#
#
# Note:
#   - get_st(): custom normalization function defined in libs/lib.R
#   - avg_silhouette(): custom function for average silhouette width defined in libs/lib.R
#
# ======================================================================


# ----------------------------------------------------------------------
# (1) Setup environment and load libraries
# ----------------------------------------------------------------------
# Set working directory
# Load required R packages and source custom functions.

setwd(dir = "")
source("libs/lib.R")

# Uncomment to install required packages
# install.packages(c("tidyverse", "phyloseq", "vegan", "indicspecies", "dendextend", "ggh4x", "ggrepel", "egg"))
library(tidyverse)
library(ggrepel)
library(phyloseq)
library(dendextend)
library(ggh4x)
library(egg)
library(indicspecies)
library(vegan)



# ----------------------------------------------------------------------
# (2) Load metadata and sequencing statistics
# ----------------------------------------------------------------------
# Load core metadata (sample labels, depths, ages, etc.)
# and sequencing summary data (initial and dereplicated read counts).

cdata <- read.table(file = "data/metadata.txt", sep = "\t", header = T)

# ----------------------------------------------------------------------
# (3) Exploratory plots of read counts and species richness
# ----------------------------------------------------------------------
# Generate summary plots showing read numbers per sample,
# species richness, and relative abundance by domain.
# Outputs: plots/fig_SX1.png

plotA <- cdata %>%
    pivot_longer(cols = c("initial", "derep"), names_to = "step", values_to = "num_seqs") %>%
    ggplot(., aes(x = y_bp, y = num_seqs, color = step))+
        geom_point()+
        geom_line()+
        coord_flip()+
        scale_x_continuous(breaks = seq(200, 2000, 100))+
        scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()))+
        labs(x = "Year (CE)", y = "Number of reads", title = "Nº of reads", color = "Reads")+
        theme(legend.position = "bottom")

plotB <- readRDS("./results/ps-unfilt-sp.rds") %>%
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

plotC <- readRDS("./results/ps-unfilt-sp.rds") %>%
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

plotD <- readRDS("./results/ps-unfilt-sp.rds") %>%
    speedyseq::psmelt() %>%
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

n_sp_all <- readRDS("./results/ps-unfilt-sp.rds") %>%
    speedyseq::psmelt() %>%
	filter(Abundance > 0) %>%
	select(label, y_bp, domain, phylum, class, order, family, genus, species) %>%
	distinct() %>%
	group_by(label) %>%
	summarise(n = n()) %>%
	ungroup()

plotE <- readRDS("./results/ps-unfilt-sp.rds") %>%
    speedyseq::psmelt() %>%
    group_by(label, domain) %>%
    summarise(Abundance = sum(Abundance)) %>%
    ungroup() %>%
    group_by(label) %>%
    mutate(proportion = Abundance/sum(Abundance)) %>%
    ungroup() %>%
    left_join(n_sp_all) %>%
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



# ----------------------------------------------------------------------
# (4) Filtering potential modern taxa
# ----------------------------------------------------------------------
# Filter out taxa showing low DNA damage (A_b < 0.05)
# and widespread occurrence (n > 30 samples).
# Save filtered phyloseq object for downstream analyses.
# Outputs: results/filtered_ps_nmds.rds

to_remove <- read_tsv("results/dmg-summary-sp.tsv.gz") %>%
	filter(domain != "d__Viruses") %>%
	group_by(species) %>%
	mutate(
		n = n(),
		max_A_b = max(A_b)) %>%
	ungroup() %>%
	filter(n > 30 & max_A_b < 0.05) %>%
	pull(species) %>% unique()


phylo_filt <- readRDS("results/ps-unfilt-sp.rds") %>%
    subset_taxa(., domain != "d__Viruses") %>%
    get_st(., nsites = 0.05, vcoeff = 1.5, mean_prop_thresh = 1e-04, norm = "median", glom = NULL, norm_out = FALSE) %>%
    prune_taxa(!(taxa_names(.) %in% to_remove), .)

saveRDS(phylo_filt, "./results/filtered_ps_nmds.rds")



# ----------------------------------------------------------------------
# (5) Alpha diversity estimation
# ----------------------------------------------------------------------
# Compute Observed, Shannon, and Simpson diversity indices
# for each sample using the filtered dataset.
# Outputs: plots/fig_SX2.png

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


# ----------------------------------------------------------------------
# (6) Compute Bray-Curtis distances and test sequencing depth effect
# ----------------------------------------------------------------------
# Transform OTU table with Hellinger transformation,
# compute Bray-Curtis dissimilarities, and perform PERMANOVA
# to test for sequencing depth bias.

distance_matrix <- phylo_filt %>%
    otu_table() %>%
    t() %>%
    as.data.frame() %>%
    vegan::decostand(., method = "hellinger") %>%
    vegan::vegdist(method = "bray")

nreads <- distance_matrix %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames() %>%
    as.data.frame() %>%
    rename(label = 1) %>%
    left_join(cdata %>% select(label, derep)) 

permanova <- adonis2(distance_matrix ~ nreads$derep, permutations = 999, dist = "bray")
permanova
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# adonis2(formula = distance_matrix ~ nreads$derep, permutations = 999, dist = "bray")
#                 Df SumOfSqs      R2      F Pr(>F)    
# nreads$derep    1   1.6789 0.09227 6.0986  0.001 ***
# Residual        60  16.5170 0.90773                  
# Total           61  18.1958 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




# ----------------------------------------------------------------------
# (7) Non-metric MultiDimensional Scaling (NMDS)
# ----------------------------------------------------------------------
# Perform NMDS ordination based on Bray-Curtis distances,
# extract sample coordinates, and calculate stress value.
# Outputs: intermediate NMDS results and plots.

set.seed(123456789)
nmdspoints <- distance_matrix %>%
    vegan::metaMDS(try = 900, k = 2) %>%
    {stress_value <<- .$stress; .} %>%
    .$points %>%
    data.frame() %>% 
    rownames_to_column("label") %>%
    left_join(cdata)
cat(paste0("stress value: ", round(stress_value, 3)))


# ----------------------------------------------------------------------
# (8) Optimal cluster number estimation
# ----------------------------------------------------------------------
# Determine the optimal number of clusters using the average silhouette method.
# This section computes silhouette width for k = 2:30.
# Merge clustering information with NMDS output for plotting.
# Outputs: plots/fig_SX3.png and fig_SX4.png

k_values <- 2:30
avg_sil_values <- sapply(k_values, function(k) avg_silhouette(hclust(sqrt(distance_matrix), method = "ward.D2"), k, distance_matrix))
optimal_k <- k_values[which.max(avg_sil_values)]
cat(paste0("Optimal number of clusters: ", optimal_k, "\n"))

plot(k_values, avg_sil_values, xlab="k",ylab="av. silhouette",type="b", pch=19)

# Add clustering information to NMDS output and extract clusters temporal ranges.
clusters_opt <- cutree(hclust(sqrt(distance_matrix), method = "ward.D2"), k = optimal_k) %>%
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

# Plot NMDS ordination and relationship of MDS1 axis with sample age
plotA <- nmdspoints_clus %>%
    left_join(cdata %>% select(label, derep)) %>% 
    ggplot(., aes(MDS1, MDS2))+
        geom_point(aes(size = derep, fill = as.factor(cluster_opt)), shape = 21)+
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



# ----------------------------------------------------------------------
# (9) Hierarchical clustering and dendrogram visualization
# ----------------------------------------------------------------------
# Cluster samples using Ward’s method (on sqrt-transformed distances),
# reorder dendrogram by sample age, and visualize cluster structure.
# Outputs: plots/fig_SX4.png

hc <- hclust(sqrt(distance_matrix), method = "ward.D2")
dend <- as.dendrogram(hc)
labels(dend) <- labels(dend) %>% data.frame() %>% rename(label = 1) %>% left_join(nmdspoints_clus %>% mutate(y_bp = round(y_bp, 0))) %>% pull(y_bp)

labels_original <- labels(dend)
years <- as.numeric(labels_original)
order_new <- order(years)
dend <- dendextend::rotate(dend, order = order_new)

Dendplot <- dend %>%
   set("branches_k_color", k=10, value = clust_cols[clust_cols %>% names() %>% as.numeric()]) %>%
   set("branches_lwd", 1.7) %>%
   set("leaves_pch", 21) %>% 
   set("leaves_cex", 1) %>%
   set("leaves_col", "white") %>%
   set("leaves_bg", "white") 

png(file = "plots/fig_SX4.png", width = 12, height = 4.5, units = "in", res = 600)
plot(Dendplot, xlab = "Samples (year)", ylab = "Height")
title("C    Hierarchical clustering dendrogram of samples based on sqrt-transformed Bray-Curtis dissimilarity", adj = 0)
dev.off()



# ----------------------------------------------------------------------
# (10) Indicator Species Analysis (IndVal)
# ----------------------------------------------------------------------
# Identify taxa significantly associated with each cluster
# using the IndVal.g method (indicspecies::multipatt).
# Apply thresholds: p < 0.01, stat > 0.5.

# Prepare data for Indval analysis
data_ind <- phylo_filt %>%
    otu_table() %>%
    as.data.frame() %>%
    t()  %>%
    as.data.frame() %>%
    rownames_to_column(var = "label") %>% 
    full_join(clusters_opt) %>%
    column_to_rownames("label")

data_ind_obs <- data_ind %>% select(-cluster_opt)

# Perform Indval analysis using multipatt function from indicspecies package
indval <- multipatt(data_ind_obs, data_ind$cluster_opt, func = "IndVal.g", duleg = TRUE)

# Set thresholds for significance and indicator value
p_value_threshold <- 0.01
indval_stat <- 0.5

# Extract significant indicator species for each cluster
indv_group <- indval$sign %>%
    as.data.frame() %>% 
    filter(p.value <= p_value_threshold & stat > indval_stat) %>%
    rownames_to_column(var = "taxa") %>%
    select(taxa, index, stat) %>%
    mutate(index = paste0("Cluster ", index))%>%
    rename(species = taxa) 


# ----------------------------------------------------------------------
# (11) Visualize abundance-weighted species richness per cluster
# ----------------------------------------------------------------------
# Plot the relative contribution of indicator species 
# to total species richness across time and indicator species relative abundance.
# Outputs: plots/fig_SX5.png

# Extract number of species per sample in the filtered phyloseq object
# to later represent the abundance-weighted number of species per sample.
n_sp <- phylo_filt %>%
    speedyseq::psmelt() %>%
	filter(Abundance > 0) %>%
	select(label, y_bp, domain, phylum, class, order, family, genus, species) %>%
	distinct() %>%
	group_by(label) %>%
	summarise(n = n()) %>%
	ungroup()


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



# ----------------------------------------------------------------------
# (12) Relative abundance of non-indicator taxa
# ----------------------------------------------------------------------
# Show the temporal distribution of taxa not assigned to any cluster group.
# Outputs: plots/fig_SX6.png

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


# ----------------------------------------------------------------------
# End of script
# ----------------------------------------------------------------------