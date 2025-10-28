# ======================================================================
# Script: 03_sourcetracker_data_preparation.R
# Author: Mateu Menendez Serra
# Affiliation: University of Copenhagen
# Date: October 2025
#
# Description:
#   This script prepares the "sink" (sample) and "source" (environmental)
#   datasets used for SourceTracker analysis. It integrates phyloseq
#   objects, filters taxa and samples, harmonizes metadata, identifies
#   ubiquitous genera shared across biomes, and exports standardized
#   BIOM and MAP files for downstream SourceTracker runs.
#
# Dependencies:
#   tidyverse, biomformat, phyloseq, speedyseq, microbiome, showtext
# ======================================================================


# ----------------------------------------------------------------------
# (1) Load libraries and set working directory
# ----------------------------------------------------------------------
library(tidyverse)
library(biomformat)
library(phyloseq)
library(showtext)

setwd(dir = "")


# ----------------------------------------------------------------------
# (2) Load sink data (target samples)
# ----------------------------------------------------------------------

# Load pre-filtered phyloseq object (excluding viruses)
phylo_filt <- readRDS("results/filtered_ps_nmds.rds")

# Extract metadata and taxonomy
cdata <- sample_data(phylo_filt)
tax_info <- tax_table(phylo_filt) %>%
  as_tibble() %>%
  rename(reference = 1) %>%
  mutate(reference = sub("S__", "", reference))

sink_data <- phylo_filt

# Extract taxonomy table (bacteria + archaea only)
sink_taxonomy <- speedyseq::psmelt(sink_data) %>%
  as_tibble() %>%
  filter(domain %in% c("d__Bacteria", "d__Archaea")) %>%
  select(reference = OTU, domain, phylum, class, order, family, genus, species) %>%
  distinct()

# Extract sink sample metadata
sink_cdata <- speedyseq::psmelt(sink_data) %>%
  as_tibble() %>%
  select(label) %>%
  distinct() %>%
  mutate(SourceSink = "sink")

# Convert OTU table to data.frame
st_df <- phylo_filt %>%
  otu_table() %>%
  t() %>%
  as.data.frame()

# Build taxonomy and metadata data.frames
st_taxonomy <- tax_info %>%
  filter(species %in% colnames(st_df)) %>%
  select(domain, phylum, class, order, family, genus, species) %>%
  distinct() %>%
  mutate(reference = species) %>%
  as.data.frame() %>%
  column_to_rownames("reference")

st_cdata <- sink_cdata %>%
  as.data.frame() %>%
  mutate(clabel = label) %>%
  column_to_rownames("clabel")

# Create phyloseq object for sink data
st_sink <- phyloseq(
  otu_table(st_df, taxa_are_rows = FALSE),
  tax_table(as.matrix(st_taxonomy)),
  sample_data(st_cdata)
)

sink_names <- sample_names(st_sink)


# ----------------------------------------------------------------------
# (3) Load source data (environmental samples)
# ----------------------------------------------------------------------

st_source <- readRDS("data/sourcetracker/st_source.rds")

# ----------------------------------------------------------------------
# (4) Source data cleaning and normalization
# ----------------------------------------------------------------------

# Remove samples dominated by a single species (>90% abundance)
to_rm <- microbiome::transform(st_source, "compositional") %>%
  speedyseq::psmelt() %>%
  as_tibble() %>%
  filter(Abundance > 0.9) %>%
  select(Sample) %>%
  distinct()

st_source <- prune_samples(!(sample_names(st_source) %in% to_rm$Sample), st_source)

# Equalize sequencing depth between source and sink
cut_s <- st_sink %>%
  sample_sums() %>%
  min()

st_source <- prune_samples(sample_sums(st_source) >= cut_s, st_source)
st_source <- prune_taxa(taxa_sums(st_source) > 0, st_source)


# ----------------------------------------------------------------------
# (5) Identify ubiquitous taxa shared across biomes
# ----------------------------------------------------------------------

# Merge metadata of source and sink
s_data <- st_source %>%
  sample_data() %>%
  as_tibble() %>%
  rbind(
    st_sink %>%
      sample_data() %>%
      as_tibble() %>%
      mutate(biome = NA, biome_class = NA, biome_subclass = NA)
  ) %>%
  select(-1)

med_st_source <- st_source
med_st_sink <- st_sink

biomes <- unique(s_data %>% filter(SourceSink == "source") %>% pull(biome_subclass))

ubiquitous_taxa <- list()

for (b in biomes) {
  print(b)

  source_samples_biome <- s_data %>%
    filter(biome_subclass == b, SourceSink == "source") %>%
    pull(label)
  sink_samples_biome <- s_data %>%
    filter(SourceSink == "sink") %>%
    pull(label)

  st_biome_source <- prune_samples(sample_names(med_st_source) %in% source_samples_biome, med_st_source)
  st_biome_sink <- prune_samples(sample_names(med_st_sink) %in% sink_samples_biome, med_st_sink)

  st_biome_source <- prune_taxa(taxa_sums(st_biome_source) > 0, st_biome_source)
  st_biome_sink <- prune_taxa(taxa_sums(st_biome_sink) > 0, st_biome_sink)

  st_biome_source_genus <- tax_table(st_biome_source) %>% as_tibble() %>% pull(genus) %>% unique()
  st_biome_sink_genus <- tax_table(st_biome_sink) %>% as_tibble() %>% pull(genus) %>% unique()

  ubiquitous_taxa[[b]] <- intersect(st_biome_source_genus, st_biome_sink_genus)
}

df <- enframe(ubiquitous_taxa, name = "habitat", value = "taxa") %>%
  unnest(cols = c(taxa))

ubiquitous_taxa <- unique(unlist(ubiquitous_taxa))


# ----------------------------------------------------------------------
# (6) Merge source and sink, filter and normalize
# ----------------------------------------------------------------------

med_st_ps <- merge_phyloseq(med_st_source, med_st_sink)
med_st_ps <- speedyseq::tax_glom(med_st_ps, "genus")

taxa_names(med_st_ps) <- tax_table(med_st_ps) %>%
  as_tibble() %>%
  pull(genus)

med_st_ps <- subset_taxa(med_st_ps, taxa_names(med_st_ps) %in% ubiquitous_taxa)

# Standard filtering thresholds
nseq_filt <- sample_sums(med_st_ps) %>%
  enframe(name = "Sample", value = "Abundance") %>%
  filter(Sample %in% sink_names) %>%
  pull(Abundance) %>%
  min()

med_st_ps <- prune_samples(sample_sums(med_st_ps) >= nseq_filt, med_st_ps)
med_st_ps <- prune_taxa(taxa_sums(med_st_ps) > 0, med_st_ps)

# Remove rare taxa based on mean relative abundance
mean_prop_thresh <- 1e-4
tax_filt <- speedyseq::psmelt(med_st_ps) %>%
  select(Sample, OTU, Abundance) %>%
  group_by(Sample) %>%
  mutate(prop = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  group_by(OTU) %>%
  summarise(mean_prop = sum(prop) / nsamples(med_st_ps)) %>%
  filter(mean_prop >= mean_prop_thresh) %>%
  pull(OTU)

med_st_ps <- prune_taxa(tax_filt, med_st_ps)
med_st_ps <- prune_taxa(taxa_sums(med_st_ps) > 0, med_st_ps)


# ----------------------------------------------------------------------
# (7) Final filtering by sample completeness
# ----------------------------------------------------------------------

taxa_counts <- apply(otu_table(med_st_ps), 1, function(x) sum(x > 0))
samples_to_keep <- names(taxa_counts)[taxa_counts >= 10]
med_st_ps <- prune_samples(samples_to_keep, med_st_ps)
med_st_ps <- prune_taxa(taxa_sums(med_st_ps) > 0, med_st_ps)


# ----------------------------------------------------------------------
# (8) Standardize counts and export SourceTracker input files
# ----------------------------------------------------------------------

total <- min(sample_sums(med_st_ps))
standf <- function(x, t = total) round(t * (x / sum(x)))
med_st_ps <- transform_sample_counts(med_st_ps, standf)

# Export BIOM file
med_st_biom <- biomformat::make_biom(
  data = t(as((otu_table(med_st_ps, taxa_are_rows = FALSE)), "matrix")),
  matrix_element_type = "int"
)
biomformat::write_biom(
  med_st_biom,
  biom_file = "./results/st-biome-class.biom"
)

# Export MAP file
med_st_ps %>%
  speedyseq::psmelt() %>%
  as_tibble() %>%
  select(Sample, SourceSink, Env = biome_subclass) %>%
  distinct() %>%
  rename("#SampleID" = Sample) %>%
  write_tsv("./results/st-biome-class.map")
