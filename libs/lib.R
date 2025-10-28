custom_theme <- theme(strip.background = element_rect(fill = "white", color = "black") , strip.text.y.right = element_text(angle = 0))
theme_set(theme_bw()+custom_theme)

paired_genus <- c("#a6cee3", "#1f78b4", "#33a02c", "#b2df8a", "#fb9a99", "#e31a1c", "#fdbf6f", "#cab2d6", "#6a3d9a", #"#ffff99", 
  "tan3", "darkslategray2", "darkslategray4", "#8dd593", "#009E73", "violet", "#d33f6a", "darkgoldenrod1", "#b5bbe3", #"khaki1",
  "purple1", "deepskyblue2", "deepskyblue4", "darkolivegreen3", "olivedrab", "coral2", "brown", "sienna1", "orangered", "thistle2", 
  "mediumpurple3", "tan1", "chocolate3", "lightsteelblue1", "lightsteelblue3", "springgreen4", "lightsalmon", "indianred1", 
  "gold", "orange", "darkorchid4", "darkorchid1", "#c6dec7", "antiquewhite", "gray", "green4", "yellowgreen", "#fdbf6f", 
  "#ef9708", "#6a3d9a", "#cab2d6", "#ffff99", "tan3", "darkslategray2", "darkslategray4", "#CC6666", "#ff7f00", "#c6dec7", 
  "seashell4", "sienna", "sienna1", "sienna2", "sienna3", "sienna4", "skyblue", "skyblue1", "skyblue2", "skyblue3", "skyblue4", 
  "slateblue", "slateblue1", "slateblue2", "slateblue3", "slateblue4", "slategray", "slategray1", "slategray2", "slategray3", 
  "slategray4", "slategrey", "snow", "snow1", "snow2", "snow3", "snow4", "springgreen", "springgreen1", "springgreen2", 
  "springgreen3", "springgreen4", "steelblue", "steelblue1", "steelblue2", "steelblue3", "steelblue4", "tan", "tan1", "tan2", 
  "tan3", "tan4", "thistle", "thistle1", "thistle2", "thistle3", "thistle4", "tomato", "tomato1", "tomato2", "tomato3", 
  "tomato4", "turquoise", "turquoise1", "turquoise2", "turquoise3", "turquoise4", "violet", "violetred", "violetred1", 
  "violetred2", "violetred3", "violetred4", "wheat", "wheat1", "wheat2", "wheat3", "wheat4", "whitesmoke", "yellow", 
  "yellow1", "yellow2", "yellow3", "yellow4", "yellowgreen", "#1f78b4", "darkslategray2", "darkolivegreen3", "#33a02c", 
  "coral2", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "green4", "yellowgreen", "lightgoldenrod1", "tan3", 
  "violet", "#d33f6a", "#8dd593", "#009E73", "sienna1", "chocolate3", "#ffff99", "orange", "#b2df8a", "limegreen", "tan1", 
  "orangered", "deepskyblue2", "#fb9a99", "brown", "gold", "darkgoldenrod1", "mediumpurple1", "purple1", "deepskyblue3", 
  "#a6cee3", "lightsalmon", "indianred1", "olivedrab", "yellowgreen", "saddlebrown", "peru", "gray48", "linen", "springgreen")



filt_physeq <- function(physeq, glom = NULL, min_counts = 10000, ncounts = 3, nsites = 0.2, vcoeff = 3, mean_prop_thresh = 1e-5) {
  suppressMessages({
    if (!is.null(glom)) {
      cat(paste("Tax glomming to:", glom, "\n"))
      physeq <- speedyseq::tax_glom(physeq, taxrank = glom)
    }

    cat(paste("Filtering samples with at least", min_counts, "counts...\n"))
    physeq_filt <- prune_samples(sample_sums(physeq) >= min_counts, physeq)

    cat(paste("Filtering taxa not seen at least", ncounts, "in", nsites, "of the samples...\n"))
    physeq_filt <- prune_taxa(taxa_sums(physeq_filt) > 0, physeq_filt)
    physeq_filt <- filter_taxa(physeq_filt, function(x) sum(x > ncounts) > (nsites * length(x)), TRUE)
    physeq_filt <- prune_samples(sample_sums(physeq_filt) >= min_counts, physeq_filt)

    # Scale to the median of the sample depth
    cat("Scaling to median sample depth", "\n")
    total <- median(sample_sums(physeq_filt))
    standf <- function(x, t = total) round(t * (x / sum(x)))
    physeq_filt_m <- transform_sample_counts(physeq_filt, standf)
    # Filter ASVs with a CV <= 3.0
    # osd2014_dada2_phyloseq_alpha_norare = filter_taxa(osd2014_dada2_phyloseq_alpha_filt, function(x) sd(x)/mean(x) > 3.0, TRUE)

    if (is.null(vcoeff)) {
      vcoeff <- speedyseq::psmelt(physeq_filt_m) |>
        as_tibble() |>
        group_by(OTU) |>
        summarise(cv = sd(Abundance) / mean(Abundance)) |>
        pull(cv) |>
        median()
      cat(paste("WARNING: Using median CV of all taxa\n"))
      cat(paste("Filtering taxa with a lower CV than", vcoeff, "...\n"))
    } else {
      cat(paste("Filtering taxa with a lower CV than", vcoeff, "...\n"))
    }

    physeq_filt_v <- filter_taxa(physeq_filt_m, function(x) sd(x) / mean(x) > vcoeff, TRUE)
    physeq_filt <- prune_taxa(taxa_names(physeq_filt) %in% taxa_names(physeq_filt_v), physeq_filt)

    # remove those with very low abundance across samples
    cat(paste("Filtering taxa with a mean proportion <=", mean_prop_thresh, "\n"))
    tax_filt <- speedyseq::psmelt(physeq_filt) |>
      dplyr::select(Sample, OTU, Abundance) |>
      group_by(Sample) |>
      mutate(prop = Abundance / sum(Abundance)) |>
      ungroup() |>
      group_by(OTU) |>
      summarise(mean_prop = sum(prop) / nsamples(physeq_filt)) |>
      ungroup() |>
      filter(mean_prop >= mean_prop_thresh) |>
      pull(OTU)
    physeq_filt <- prune_taxa(tax_filt, physeq_filt)
    physeq_filt <- prune_taxa(taxa_sums(physeq_filt) > 0, physeq_filt)
    physeq_filt <- prune_samples(sample_sums(physeq_filt) >= 0, physeq_filt)
    cat(paste0("Kept ", ntaxa(physeq_filt), " taxa (Original:", ntaxa(physeq), ") summing a total of ", sum(taxa_sums(physeq_filt)), " counts (Original:", sum(taxa_sums(physeq)), ")\n"))
    return(physeq_filt)
  })
}

