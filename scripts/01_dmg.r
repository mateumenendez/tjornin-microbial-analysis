library(tidyverse)
library(knitr)
library(readxl)
library(measurements)
library(ggmap)
library(ggrepel)
library(lvplot)


setwd(dir = "/projects/caeg/people/ngm902/apps/repos/tjornin-microbial-analysis")
source("libs/dmg.R")


source("/projects/fernandezguerra/people/ngm902/scripts/r-miscellaneous.R")
# source("/projects/fernandezguerra/people/ngm902/scripts/get-metadata.R")




cdata <- read.table(file = "data/metadata.txt", sep = "\t", header = T)

tax_data <- read_tsv("data/aggregated/tp-lca.summary.tsv.gz") %>%
    filter(rank %in% c("subspecies", "species", "genus")) %>% 
    separate(tax_path, into = c("root", "domain", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "specie", "subspecie"), sep = ";") %>%
    filter(domain != "d__Eukaryota") %>%
    inner_join(cdata)

# Let's get the damage at each node
dmg <- read_tsv("data/aggregated/tp-damage.tsv.gz") %>%
    filter(label %in% cdata$label) |>
    inner_join(tax_data |> select(label, taxid, name, rank))

# For subspecies
# Let's identify all those taxa with a bad fit
dmg_ssp <- dmg |>
    filter(rank == "subspecies") |>
    filter(nreads >= 100) |>
    rename(tax_name = taxid, n_reads = nreads)

samples <- cdata$label |> unique()

dat <- dmg_fwd_CCC(dmg_ssp, samples, ci = "asymptotic", nperm = 100, nproc = 24)

dat1 <- dat |>
    ungroup() |>
    inner_join(dmg_ssp) |>
    mutate(fit = ifelse(rho_c >= 0.85 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) |>
    mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))

dat1 %>%
    select(label, tax_name, A_b, q_b, c_b, phi_b, fit) |>
    write_tsv("./results/aggregated/dmg-subspecies-lite.tsv.gz")





# For species
# Let's identify all those taxa with a bad fit
dmg_sp <- dmg |>
    filter(rank == "species") |>
    filter(nreads >= 100) |>
    rename(tax_name = taxid, n_reads = nreads)

samples <- cdata$label |> unique()

dat <- dmg_fwd_CCC(dmg_sp, samples, ci = "asymptotic", nperm = 100, nproc = 24)

dat1 <- dat |>
    ungroup() |>
    inner_join(dmg_sp) |>
    mutate(fit = ifelse(rho_c >= 0.85 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) |>
    mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))

dat1 %>%
    select(label, tax_name, A_b, q_b, c_b, phi_b, fit) |>
    write_tsv("./results/aggregated/dmg-species-lite.tsv.gz")








# plotA <- dat1 |>
#     group_by(label, fit) |>
#     count() |>
#     ungroup() |>
#     inner_join(cdata) |>
#     mutate(label2 = paste(label, " (", y_bp, ")", sep = "")) |>
#     arrange(desc(depth_in_core_cm)) %>%
#     mutate(label2 = factor(label2, levels = unique(label2))) |>
#     ggplot(aes(n, label2, fill = fit)) +
#         geom_col() +
#         theme(
#             legend.position = "bottom",
#             legend.title = element_blank(),
#             text = element_text(size = 10),
#         ) +
#         labs(x = "Number of references", y = "Samples sorted by age", title = "Number of references exhibiting good and bad fit to the expected damage pattern") 
        
# plotB <- dat1 |>
#     group_by(label, fit) |>
#     count() |>
#     ungroup() |> 
#     group_by(label) %>%
#     summarize(
#         total = sum(n),
#         good_count = sum(n[fit == "good"]),
#         percentage_good = good_count / total
#     ) %>% 
#     ungroup() %>% 
#     inner_join(cdata) |> 
#     mutate(label2 = paste(label, " (", y_bp, ")", sep = "")) |>
#     arrange(desc(depth_in_core_cm)) %>%
#     mutate(label2 = factor(label2, levels = unique(label2))) |>
#     ggplot(aes(percentage_good, label2)) +
#         geom_point() +
#         geom_line()+
#         theme(
#             axis.text.y =element_blank(), 
#             axis.ticks.y = element_blank(), 
#             legend.position = "top",
#             legend.title = element_blank(),
#             text = element_text(size = 10),
#         ) +
#         labs(x = "Proportion of good fits", y = "", title = "Proportion of references exhibiting good and bad fit") +
#         scale_x_continuous(labels = scales::percent_format())+
#         guides(fill = guide_legend(nrow = 1))


# png(file = "plots/aggregated/preliminary_assessment/damage_fit_subspecies_ss-dmg.png", width = 13, height = 10, units = "in", res = 600)
# egg::ggarrange(plotA, plotB, nrow = 1, widths = c(2, 1))
# dev.off()



# ## Smiley plots
# ## Example of good fit (9 references)
# dat_filt <- dat1 %>%
#     left_join(cdata) %>%
#     filter(fit == "good")

# samp <- dat_filt %>% pull(label) %>% unique()
# select <- dat_filt %>% filter(label == samp[1]) %>% pull(name) %>% unique() %>% sample(., 9)

# plots100 <- purrr::map(.x = samp[1], dat = dat_filt %>% filter(name %in% select), .f = function(x, dat, ncol = 5, nrow = 2, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
#     data <- dat |>
#             filter(label == x) %>%
#             arrange(name)
#     grid_size <- calculate_plot_grid(length(data$tax_name))
#     l <- lapply(data$name, function(X) {
#     df1 <- data |>
#                 filter(name == X)
#     p <- get_dmg_decay_fit(df1, orient = orient, pos = pos, p_breaks = p_breaks)
#     p <- p + ggtitle(X)
#     return(p)
#         })
#     plot <- ggpubr::ggarrange(plotlist = l, ncol = grid_size$cols, nrow = grid_size$rows, align = "hv")
#     }, , .progress = TRUE)

# png(file = "plots/aggregated/preliminary_assessment/smiley_plot_example_good_fit_example.png", width = 20, height = 15, units = "in", res = 600)
# plots100
# dev.off()


# ## Example of bad fit (9 references)
# dat_filt <- dat1 %>%
#     left_join(cdata) %>%
#     filter(fit == "bad")

# samp <- dat_filt %>% pull(label) %>% unique()
# select <- dat_filt %>% filter(label == samp[1]) %>% pull(name) %>% unique() %>% sample(., 9)

# plots100 <- purrr::map(.x = samp[1], dat = dat_filt %>% filter(name %in% select), .f = function(x, dat, ncol = 5, nrow = 2, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
#     data <- dat |>
#             filter(label == x) %>%
#             arrange(name)
#     grid_size <- calculate_plot_grid(length(data$tax_name))
#     l <- lapply(data$name, function(X) {
#     df1 <- data |>
#                 filter(name == X)
#     p <- get_dmg_decay_fit(df1, orient = orient, pos = pos, p_breaks = p_breaks)
#     p <- p + ggtitle(X)
#     return(p)
#         })
#     plot <- ggpubr::ggarrange(plotlist = l, ncol = grid_size$cols, nrow = grid_size$rows, align = "hv")
#     }, , .progress = TRUE)


# png(file = "plots/aggregated/preliminary_assessment/smiley_plot_example_bad_fit_example.png", width = 20, height = 15, units = "in", res = 600)
# plots100
# dev.off()










## Subspecies level
lca_ssp <- read_tsv("data/aggregated/tp-lca.summary.tsv.gz") %>%
    filter(rank == "subspecies") |>
    separate(tax_path, into = c("root", "domain", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies"), sep = ";") %>%
    filter(domain != "d__Eukaryota") %>%
    rename(nreads = n_reads) %>%
    left_join(cdata)

# Mapping stats from filterBAM
mapping_stats <- read_tsv("data/aggregated/tp-mapping-filtered.summary.tsv.gz") |>
    filter(norm_entropy > 0.6, breadth > 0.01, breadth_exp_ratio >= 0.75) |> 
    mutate(name = paste0("S__", reference)) |>
    mutate(abundance = ifelse(tax_abund_tad == 0, tax_abund_read, tax_abund_tad))

tax_data_ssp <- lca_ssp |>
    select(-nreads, -abundance) |>
    inner_join(mapping_stats)

# Add dmage. If the fit is bad, the damage will be 0
dmg_ssp <- read_tsv("results/aggregated/dmg-subspecies-lite.tsv.gz") |>
    rename(taxid = tax_name) |>
    mutate(damage = ifelse(fit == "bad", 0, A_b))

# Combine damage and taxonomic assignments
tax_data_dmg <- tax_data_ssp %>% 
    inner_join(dmg_ssp)

tax_data_dmg <- tax_data_dmg |>
    ungroup() |>
    left_join(cdata) %>%
    mutate(is_dmg = ifelse(damage > 0.1, "Damaged", "Non-damaged")) |>
    mutate(is_dmg = ifelse(fit == "bad", "Non-damaged", is_dmg))

tax_data_dmg |>
    write_tsv("results/aggregated/dmg-summary-ssp.tsv.gz")







## Species level
lca_sp <- read_tsv("data/aggregated/tp-lca.summary.tsv.gz") %>%
    filter(rank == "species") |>
    separate(tax_path, into = c("root", "domain", "lineage", "kingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies"), sep = ";") %>%
    filter(domain != "d__Eukaryota") %>%
    rename(nreads = n_reads) %>%
    left_join(cdata)

# Add dmage. If the fit is bad, the damage will be 0
dmg_sp <- read_tsv("./results/aggregated/dmg-species-lite.tsv.gz") |>
    rename(taxid = tax_name) |>
    mutate(damage = ifelse(fit == "bad", 0, A_b))

# Combine damage and taxonomic assignments
tax_data_dmg <- lca_sp %>%
    inner_join(dmg_sp)

tax_data_dmg <- tax_data_dmg |>
    ungroup() |>
    left_join(cdata) %>%
    mutate(is_dmg = ifelse(damage > 0.1, "Damaged", "Non-damaged")) |>
    mutate(is_dmg = ifelse(fit == "bad", "Non-damaged", is_dmg))

tax_data_dmg |>
    write_tsv("results/aggregated/dmg-summary-sp.tsv.gz")

