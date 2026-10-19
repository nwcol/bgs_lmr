
#################
# 
# This script loads the tables generated with moments++ and bgshr
# for different models (constrained elements, mutation maps)
# and prepares them for downstream SEM analyses and visualization
# NOTE: this script should be run from within each current_models dir
#
#################

source("~/Devel/bgs_lmr/tools/R_scripts/HelperFunctions.R")

constrained_models <- c("merged_cds_phastcons", "split_cds_phastcons", "split_cds_regulatory", "merged_cds_regulatory")
mutation_maps <- c("carlson", "gnomad", "roulette")

human_maps_10kb <- vector("list", length=22 * 4 * length(constrained_models)) 
pos <- 1
for(u in mutation_maps) {
  for(m in constrained_models) {
    print(paste("Retrieving:", m, "|", u, "..."))
    pb <- txtProgressBar(min=0, max=22, style=3)
    for(c in 1:22) {
      setTxtProgressBar(pb, c)
      tmp <- fread(paste0(m, "/", u, "/B_tbl_YRI_chr", c, "_10kb.csv.gz"))
      
      tmp <- dplyr::select(tmp, -exp_del_pi) %>% setDT()
      
      tmp[, mut_map := u]
      tmp[, constrained_model := m]
      tmp[, demography := "equilibrium"]
      
      human_maps_10kb[[pos]] <- tmp
      pos <- pos + 1
    }
    close(pb)
  }
}

human_maps_10kb <- data.table::rbindlist(human_maps_10kb)
human_maps_10kb <- dplyr::arrange(human_maps_10kb, mut_map, constrained_model, demography, chrom, chromStart) %>% as.data.table() %>% setDT()
human_maps_10kb[, scale := 1e4]

# introducing mutation noise in expected pi
S <- 25 # diploid sample size
window_size <- 1e+4

simulate_pi <- function(exp_pi, S, window_size) {
  n_windows <- length(exp_pi)
  
  pi_sim <- matrix(
    rbinom(n_windows * S, size = window_size, prob = exp_pi) / window_size,
    nrow = n_windows,
    ncol = S)
  
  rowMeans(pi_sim)
}

human_maps_10kb[, pi_sim_tot := simulate_pi(exp_pi, S = S, window_size = 1e4),
                  by = .(mut_map, constrained_model)]

human_maps_10kb[, coverage := num_sites / 1e+4]
human_maps_10kb[, prop_del_sites := del_sites / num_sites]

human_maps_10kb <- human_maps_10kb %>% group_by(constrained_model, mut_map, chrom) %>%
  mutate(bin_100kb=chromEnd %/% 1e+5, bin_1Mb=chromEnd %/% 1e+6) %>% setDT()

# 100 kb binning
setnames(human_maps_10kb, old="num_sites", new="weight")
human_maps_100kb <- human_maps_10kb %>%
  group_by(constrained_model, mut_map, chrom, bin_100kb) %>% 
  summarise(
    num_sites = sum(weight),
    del_sites = sum(del_sites),
    avg_mut = weighted.mean(avg_mut, weight, na.rm = TRUE),
    avg_rec = weighted.mean(avg_rec, weight, na.rm = TRUE),
    B = weighted.mean(B, weight, na.rm = TRUE),
    exp_pi = weighted.mean(exp_pi, weight, na.rm = TRUE),
    avg_pi = weighted.mean(avg_pi, weight, na.rm = TRUE),
    pi_sim_tot = weighted.mean(pi_sim_tot, weight, na.rm = TRUE),
    prop_del_sites = weighted.mean(prop_del_sites, weight, na.rm = TRUE),
    .groups = "drop") %>% setDT()

human_maps_100kb[, scale := 1e5]
human_maps_100kb[, chromStart := as.integer(bin_100kb * 1e5)]
human_maps_100kb[, chromEnd := as.integer((bin_100kb + 1) * 1e+5)]
human_maps_100kb[, coverage:= num_sites / 1e+5]
human_maps_100kb[, prop_del_sites := del_sites / num_sites]
human_maps_100kb[, bin_100kb := NULL]

setcolorder(human_maps_100kb, c("chrom", "chromStart", "chromEnd", "constrained_model", "mut_map", "scale", "coverage", "num_sites", "del_sites",
                                "prop_del_sites", "avg_mut", "avg_rec", "B", "exp_pi", "avg_pi", "pi_sim_tot"))
fwrite(human_maps_100kb, "maps/maps_100kb.csv.gz")

# 1 Mb binning
human_maps_1Mb <- human_maps_10kb %>%
  group_by(constrained_model, mut_map, chrom, bin_1Mb) %>% 
  summarise(
    num_sites = sum(weight),
    del_sites = sum(del_sites),
    avg_mut = weighted.mean(avg_mut, weight, na.rm = TRUE),
    avg_rec = weighted.mean(avg_rec, weight, na.rm = TRUE),
    B = weighted.mean(B, weight, na.rm = TRUE),
    exp_pi = weighted.mean(exp_pi, weight, na.rm = TRUE),
    avg_pi = weighted.mean(avg_pi, weight, na.rm = TRUE),
    pi_sim_tot = weighted.mean(pi_sim_tot, weight, na.rm = TRUE),
    prop_del_sites = weighted.mean(prop_del_sites, weight, na.rm = TRUE),
    .groups = "drop") %>% setDT()

human_maps_1Mb[, scale := 1e6]
human_maps_1Mb[, chromStart := as.integer(bin_1Mb * 1e6)]
human_maps_1Mb[, chromEnd := as.integer((bin_1Mb + 1) * 1e+6)]
human_maps_1Mb[, coverage:= num_sites / 1e+6]
human_maps_1Mb[, prop_del_sites := del_sites / num_sites]
human_maps_1Mb[, bin_1Mb := NULL]

setcolorder(human_maps_1Mb, c("chrom", "chromStart", "chromEnd", "constrained_model", "mut_map", "scale", "coverage", "num_sites", "del_sites",
                              "prop_del_sites", "avg_mut", "avg_rec", "B", "exp_pi", "avg_pi", "pi_sim_tot"))
fwrite(human_maps_1Mb, "maps/maps_1Mb.csv.gz")

human_maps_10kb[, bin_100kb := NULL]
human_maps_10kb[, bin_1Mb := NULL]
setnames(human_maps_10kb, old="weight", new="num_sites")

setcolorder(human_maps_10kb, c("chrom", "chromStart", "chromEnd", "constrained_model", "mut_map", "scale", "coverage", "num_sites", "del_sites",
                               "prop_del_sites", "avg_mut", "avg_rec", "B", "exp_pi", "avg_pi", "pi_sim_tot"))
fwrite(human_maps_10kb, "maps/maps_10kb.csv.gz")