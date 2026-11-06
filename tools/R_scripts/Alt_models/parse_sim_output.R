library(data.table)
library(tidyverse)

setwd("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/10Mb/")
dirs <- c("demog_predictions", "fit_Ne_predictions_good_DFE", "fit_Ne_predictions_bad_DFE")

for(d in dirs) {
  sim_B_tbls <- list.files(d, pattern="*csv.gz")
  sim_pi <- list.files("observed_diversity")
  
  sim_maps_1kb <- vector("list", length=length(sim_B_tbls)) 
  pb <- txtProgressBar(min=0, max=length(sim_B_tbls), style=3)
  for(i in 1:length(sim_B_tbls)) {
  
    setTxtProgressBar(pb, i)
    B <- fread(paste0(d, "/", sim_B_tbls[i]))
    pi <- fread(paste0("observed_diversity/", sim_pi[i]))
    
    B[, avg_pi := pi$ave_pi]
    sim_maps_1kb[[i]] <- B
  }
  close(pb)
  
  sim_maps_1kb <- data.table::rbindlist(sim_maps_1kb)
  sim_maps_1kb <- dplyr::arrange(sim_maps_1kb, chrom, chromStart) %>% as.data.table() %>% setDT()
  
  sim_maps_1kb <- sim_maps_1kb %>% 
    mutate(bin_10kb=(chromEnd - 1e3) %/% 1e+4, bin_100kb=(chromEnd - 1e3) %/% 1e+5, bin_1Mb=(chromEnd - 1e3) %/% 1e+6) %>% setDT()
  
  # 10 kb binning
  sim_maps_10kb <- sim_maps_1kb %>%
    group_by(chrom, bin_10kb) %>% 
    summarise(
      chrom=mean(chrom),
      avg_mut=mean(avg_mut),
      avg_rec=mean(avg_rec),
      B=mean(B),
      exp_pi=mean(exp_pi),
      avg_pi=mean(avg_pi),
      prop_del_sites=sum(del_sites) / sum(num_sites),
      .groups="drop") %>% setDT()
  
  sim_maps_10kb[, scale := 1e4]
  sim_maps_10kb[, chromStart := as.integer(bin_10kb * 1e4)]
  sim_maps_10kb[, chromEnd := as.integer((bin_10kb + 1) * 1e+4)]
  sim_maps_10kb[, bin_10kb := NULL]
  
  setcolorder(sim_maps_10kb, c("chrom", "chromStart", "chromEnd", "scale", "prop_del_sites",
                               "avg_mut", "avg_rec", "B", "exp_pi", "avg_pi"))
  fwrite(sim_maps_10kb, paste0(d, "/maps/maps_10kb.csv.gz"))
  
  # 100 kb binning
  sim_maps_100kb <- sim_maps_1kb %>%
    group_by(chrom, bin_100kb) %>% 
    summarise(
      chrom=mean(chrom),
      avg_mut=mean(avg_mut),
      avg_rec=mean(avg_rec),
      B=mean(B),
      exp_pi=mean(exp_pi),
      avg_pi=mean(avg_pi),
      prop_del_sites=sum(del_sites) / sum(num_sites),
      .groups="drop") %>% setDT()
  
  sim_maps_100kb[, scale := 1e5]
  sim_maps_100kb[, chromStart := as.integer(bin_100kb * 1e5)]
  sim_maps_100kb[, chromEnd := as.integer((bin_100kb + 1) * 1e+5)]
  sim_maps_100kb[, bin_100kb := NULL]
  
  setcolorder(sim_maps_100kb, c("chrom", "chromStart", "chromEnd", "scale", "prop_del_sites",
                               "avg_mut", "avg_rec", "B", "exp_pi", "avg_pi"))
  fwrite(sim_maps_100kb, paste0(d, "/maps/maps_100kb.csv.gz"))
  
  # 1 Mb binning
  sim_maps_1Mb <- sim_maps_1kb %>%
    group_by(chrom, bin_1Mb) %>% 
    summarise(
      chrom=mean(chrom),
      avg_mut=mean(avg_mut),
      avg_rec=mean(avg_rec),
      B=mean(B),
      exp_pi=mean(exp_pi),
      avg_pi=mean(avg_pi),
      prop_del_sites=sum(del_sites) / sum(num_sites),
      .groups="drop") %>% setDT()
  
  sim_maps_1Mb[, scale := 1e6]
  sim_maps_1Mb[, chromStart := as.integer(bin_1Mb * 1e6)]
  sim_maps_1Mb[, chromEnd := as.integer((bin_1Mb + 1) * 1e+6)]
  sim_maps_1Mb[, bin_1Mb := NULL]
  
  setcolorder(sim_maps_1Mb, c("chrom", "chromStart", "chromEnd", "scale", "prop_del_sites",
                                "avg_mut", "avg_rec", "B", "exp_pi", "avg_pi"))
  fwrite(sim_maps_1Mb, paste0(d, "/maps/maps_1Mb.csv.gz"))
}
