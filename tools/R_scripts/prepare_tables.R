
#################
# 
# This script loads the tables generated with moments++ and bgshr
# for different models (constrained elements, mutation maps)
# and prepares them for downstream SEM analyses and visualization
# NOTE: this script should be run from within each current_models dir
#
#################

source("~/Data/bgs_lmr/human_data/tools/HelperFunctions.R")

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

######### Bonus 2D plots ######### 
# Finding the parameter combination which results in simulated maps most similar to the real ones 
plot_2d <- F
if(plot_2d) { 
  
  # summarizing in terms of AC and CV (not for uniform mut maps as they are flat)
  human_summaries_eq <- human_maps_eq %>% filter(., constrained_model=="phastcons") %>% # arbitrary, they are identical
    group_by(chrom, mut_map, scale) %>%
    mutate(cv_mut=sd(avg_mut)/mean(avg_mut)) %>%
    mutate(cv_rec=sd(avg_rec)/mean(avg_rec)) %>%
    mutate(ac_mut=as.numeric(unlist(acf(avg_mut, lag=1, pl=F))[2])) %>%
    mutate(ac_rec=as.numeric(unlist(acf(avg_rec, lag=1, pl=F))[2])) %>%
    distinct(chrom, .keep_all=T) %>% ungroup() %>% setDT()
  
  human_summaries_eq[, chr := paste0("chr", chrom)]
  
  human_summaries_eq <- dplyr::select(human_summaries_eq,
                                      c(chr, cv_mut, cv_rec, ac_mut, ac_rec, scale, mut_map))
  
  # summarizing in terms of AC and CV (not for uniform mut maps as they are flat)
  human_summaries_4e <- human_maps_4e %>% filter(., constrained_model=="phastcons") %>% # arbitrary, they are identical
    group_by(chrom, mut_map, scale) %>%
    mutate(cv_mut=sd(avg_mut)/mean(avg_mut)) %>%
    mutate(cv_rec=sd(avg_rec)/mean(avg_rec)) %>%
    mutate(ac_mut=as.numeric(unlist(acf(avg_mut, lag=1, pl=F))[2])) %>%
    mutate(ac_rec=as.numeric(unlist(acf(avg_rec, lag=1, pl=F))[2])) %>%
    distinct(chrom, .keep_all=T) %>% ungroup() %>% setDT()
  
  human_summaries_4e[, chr := paste0("chr", chrom)]
  
  human_summaries_4e <- dplyr::select(human_summaries_4e,
                                      c(chr, cv_mut, cv_rec, ac_mut, ac_rec, scale, mut_map))
  
  # collecting simulated maps
  muts_10kb <- as.data.frame(matrix(nrow=num_models*num_reps, ncol=1e+4))
  recs_10kb <- as.data.frame(matrix(nrow=num_models*num_reps, ncol=1e+4))
  
  muts_100kb <- as.data.frame(matrix(nrow=num_models*num_reps, ncol=1e+3))
  recs_100kb <- as.data.frame(matrix(nrow=num_models*num_reps, ncol=1e+3))
  
  muts_1Mb <- as.data.frame(matrix(nrow=num_models*num_reps, ncol=1e+2))
  recs_1Mb <- as.data.frame(matrix(nrow=num_models*num_reps, ncol=1e+2))
  
  pb <- txtProgressBar(min=0, max=num_models, style=3)
  for(m in 1:num_models) {
    setTxtProgressBar(pb, m)
    for(rep in 1:num_reps) {
      maps_1kb <- fread(paste("model_", m,
                              "/rep_", rep, 
                              "/maps_1kb.csv.gz", sep=""))
      
      maps_1kb[, bin_1kb   := .I]
      maps_1kb[, bin_10kb  := (bin_1kb - 1) %/% 10]
      maps_1kb[, bin_100kb := (bin_1kb - 1) %/% 100]
      maps_1kb[, bin_1Mb   := (bin_1kb - 1) %/% 1000]
      
      maps_10kb <- maps_1kb %>% group_by(bin_10kb) %>% 
        summarise_at(c("avg_mut", "avg_rec"), mean)
      maps_100kb <- maps_1kb %>% group_by(bin_100kb) %>% 
        summarise_at(c("avg_mut", "avg_rec"), mean)
      maps_1Mb <- maps_1kb %>% group_by(bin_1Mb) %>% 
        summarise_at(c("avg_mut", "avg_rec"), mean)
      
      idx <- rep + (m - 1) * num_reps
      
      muts_10kb[idx,] <- maps_10kb$avg_mut
      recs_10kb[idx,] <- maps_10kb$avg_rec
      
      muts_100kb[idx,] <- maps_100kb$avg_mut
      recs_100kb[idx,] <- maps_100kb$avg_rec
      
      muts_1Mb[idx,] <- maps_1Mb$avg_mut
      recs_1Mb[idx,] <- maps_1Mb$avg_rec
    }
  }
  close(pb)
  
  # mutation maps
  names(muts_10kb) <- paste("bin", 1:ncol(muts_10kb), sep="_")
  muts_10kb$model <- sort(rep(1:num_models, num_reps))
  muts_10kb$replicate <- rep(1:10, num_models)
  
  names(muts_100kb) <- paste("bin", 1:ncol(muts_100kb), sep="_")
  muts_100kb$model <- sort(rep(1:num_models, num_reps))
  muts_100kb$replicate <- rep(1:10, num_models)
  
  names(muts_1Mb) <- paste("bin", 1:ncol(muts_1Mb), sep="_")
  muts_1Mb$model <- sort(rep(1:num_models, num_reps))
  muts_1Mb$replicate <- rep(1:10, num_models)
  
  # recombination maps
  names(recs_10kb) <- paste("bin", 1:ncol(recs_10kb), sep="_")
  recs_10kb$model <- sort(rep(1:num_models, num_reps))
  recs_10kb$replicate <- rep(1:10, num_models)
  
  names(recs_100kb) <- paste("bin", 1:ncol(recs_100kb), sep="_")
  recs_100kb$model <- sort(rep(1:num_models, num_reps))
  recs_100kb$replicate <- rep(1:10, num_models)
  
  names(recs_1Mb) <- paste("bin", 1:ncol(recs_1Mb), sep="_")
  recs_1Mb$model <- sort(rep(1:num_models, num_reps))
  recs_1Mb$replicate <- rep(1:10, num_models)
  
  fwrite(muts_10kb, "mut_maps_10kb_models.csv.gz")
  fwrite(muts_100kb, "mut_maps_100kb_models.csv.gz")
  fwrite(muts_1Mb, "mut_maps_1Mb_models.csv.gz")
  
  fwrite(recs_10kb, "rec_maps_10kb_models.csv.gz")
  fwrite(recs_100kb, "rec_maps_100kb_models.csv.gz")
  fwrite(recs_1Mb, "rec_maps_1Mb_models.csv.gz")
  
  # Plotting 2D maps (AC x CV), simulated and real
  
  muts_10kb <- fread("mut_maps_10kb_models.csv.gz")
  muts_100kb <- fread("mut_maps_100kb_models.csv.gz")
  muts_1Mb <- fread("mut_maps_1Mb_models.csv.gz")
  
  recs_10kb <- fread("rec_maps_10kb_models.csv.gz")
  recs_100kb <- fread("rec_maps_100kb_models.csv.gz")
  recs_1Mb <- fread("rec_maps_1Mb_models.csv.gz")
  
  # mutation maps
  ## 10 kb
  muts_10kb$replicate <- rep(1:10, num_models)
  muts_10kb_models <- left_join(muts_10kb, models, by="model") %>%
    pivot_longer(cols=starts_with("bin"), 
                 values_to="mu",
                 names_to="bin")
  
  muts_10kb_models <- muts_10kb_models %>% group_by(model, replicate) %>%
    mutate(cv_mut=sd(mu)/mean(mu)) %>%
    distinct(model, 
             replicate, 
             .keep_all=T) %>%
    dplyr::select(., -bin)
  
  muts_10kb_models$ac_mut <- apply(muts_10kb[,1:100], 1, # 1st order autocorrelation
                                   function(x) as.numeric(unlist(acf(x, lag=1, pl=F))[2]))
  muts_10kb_models$scale <- 1e+4
  
  ## 100 kb
  muts_100kb$replicate <- rep(1:10, num_models)
  muts_100kb_models <- left_join(muts_100kb, models, by="model") %>%
    pivot_longer(cols=starts_with("bin"), 
                 values_to="mu",
                 names_to="bin")
  
  muts_100kb_models <- muts_100kb_models %>% group_by(model, replicate) %>%
    mutate(cv_mut=sd(mu)/mean(mu)) %>%
    distinct(model, replicate, .keep_all=T) %>%
    dplyr::select(., -bin)
  
  muts_100kb_models$ac_mut <- apply(muts_100kb[,1:100], 1, # 1st order autocorrelation
                                    function(x) as.numeric(unlist(acf(x, lag=1, pl=F))[2]))
  muts_100kb_models$scale <- 1e+5
  
  ## 1 Mb
  muts_1Mb$replicate <- rep(1:10, num_models)
  muts_1Mb_models <- left_join(muts_1Mb, models, by="model") %>%
    pivot_longer(cols=starts_with("bin"), 
                 values_to="mu",
                 names_to="bin")
  
  muts_1Mb_models <- muts_1Mb_models %>% group_by(model, replicate) %>%
    mutate(cv_mut=sd(mu)/mean(mu)) %>%
    distinct(model, replicate, .keep_all=T) %>%
    dplyr::select(., -bin)
  
  muts_1Mb_models$ac_mut <- apply(muts_1Mb[,1:100], 1, # 1st order autocorrelation
                                  function(x) as.numeric(unlist(acf(x, lag=1, pl=F))[2]))
  muts_1Mb_models$scale <- 1e+6
  
  muts <- rbind.data.frame(muts_10kb_models, muts_100kb_models, muts_1Mb_models)
  
  p <- ggplot(data=muts, 
              aes(y=cv_mut, 
                  x=ac_mut, 
                  color=as.factor(shapes_rate_u),
                  shape=as.factor(avg_mut_spans))) +
    geom_point(size=3) + theme_bw() +
    geom_point(data=human_summaries, 
               aes(x=ac_mut, y=cv_mut, color=mut_map),
               shape=4) +
    #geom_text(data=human_summaries, aes(label=chr), vjust=0.8, size=10) +
    facet_grid(scale~mult_r_u) +
    scale_shape_manual(values=c(0, 1, 2), name="Avg. Mut. Spans") + 
    # https://stats.stackexchange.com/questions/185425/how-to-determine-the-critical-values-of-acf
    geom_vline(xintercept=-1.96 / sqrt(100-1), color="grey", linetype="dashed") + 
    geom_vline(xintercept=1.96 / sqrt(100-1), color="grey", linetype="dashed") + 
    scale_color_manual(values=c("brown1", "green4", "purple3", "cyan3", 
                                "grey1", "blue", "orange"),
                       name="Shape=Rate Gamma") + 
    scale_y_continuous(breaks=pretty_breaks()) +
    scale_x_continuous(breaks=pretty_breaks()) +
    labs(title=NULL, x="AC (mutation rates)", y="CV (mutation rates)") +
    theme(strip.text.y=element_text(size=16),        
          strip.text.x=element_text(size=16),        
          axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          axis.text.x=element_text(size=12),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16),
          legend.position="bottom")
  save_plot("2D_plot_mut.svg", p, base_height=12, base_width=18)
  
  # recombination maps
  ## 10 kb
  recs_10kb$replicate <- rep(1:10, num_models)
  recs_10kb_models <- left_join(recs_10kb, models, by="model") %>%
    pivot_longer(cols=starts_with("bin"), 
                 values_to="rec",
                 names_to="bin")
  
  recs_10kb_models <- recs_10kb_models %>% group_by(model, replicate) %>%
    mutate(cv_rec=sd(rec)/mean(rec)) %>%
    distinct(model, replicate, .keep_all=T) %>%
    dplyr::select(., -bin)
  
  recs_10kb_models$ac_rec <- apply(recs_10kb[,1:100], 1, # 1st order autocorrelation
                                   function(x) as.numeric(unlist(acf(x, lag=1, pl=F))[2]))
  recs_10kb_models$scale <- 1e+4
  
  ## 100 kb
  recs_100kb$replicate <- rep(1:10, num_models)
  recs_100kb_models <- left_join(recs_100kb, models, by="model") %>%
    pivot_longer(cols=starts_with("bin"), 
                 values_to="rec",
                 names_to="bin")
  
  recs_100kb_models <- recs_100kb_models %>% group_by(model, replicate) %>%
    mutate(cv_rec=sd(rec)/mean(rec)) %>%
    distinct(model, replicate, .keep_all=T) %>%
    dplyr::select(., -bin)
  
  recs_100kb_models$ac_rec <- apply(recs_100kb[,1:100], 1, # 1st order autocorrelation
                                    function(x) as.numeric(unlist(acf(x, lag=1, pl=F))[2]))
  recs_100kb_models$scale <- 1e+5
  
  ## 1 Mb
  recs_1Mb$replicate <- rep(1:10, num_models)
  recs_1Mb_models <- left_join(recs_1Mb, models, by="model") %>%
    pivot_longer(cols=starts_with("bin"), 
                 values_to="rec",
                 names_to="bin")
  
  recs_1Mb_models <- recs_1Mb_models %>% group_by(model, replicate) %>%
    mutate(cv_rec=sd(rec)/mean(rec)) %>%
    distinct(model, replicate, .keep_all=T) %>%
    dplyr::select(., -bin)
  
  recs_1Mb_models$ac_rec <- apply(recs_1Mb[,1:100], 1, # 1st order autocorrelation
                                  function(x) as.numeric(unlist(acf(x, lag=1, pl=F))[2]))
  recs_1Mb_models$scale <- 1e+6
  
  recs <- rbind.data.frame(recs_10kb_models, recs_100kb_models, recs_1Mb_models)
  
  q <- ggplot(data=recs, 
              aes(y=cv_rec, 
                  x=ac_rec, 
                  color=as.factor(shapes_rate_r),
                  shape=as.factor(avg_rec_spans))) +
    geom_point(size=3) + theme_bw() +
    geom_point(data=human_summaries, 
               aes(x=ac_rec, y=cv_rec),
               shape=4, color="black") +
    #geom_text(data=human_summaries, aes(label=chr), vjust=0.8, size=10) +
    facet_grid(scale~mult_r_u) +
    scale_shape_manual(values=c(0, 1, 2), name="Avg. Mut. Spans") + 
    # https://stats.stackexchange.com/questions/185425/how-to-determine-the-critical-values-of-acf
    geom_vline(xintercept=-1.96 / sqrt(100-1), color="grey", linetype="dashed") + 
    geom_vline(xintercept=1.96 / sqrt(100-1), color="grey", linetype="dashed") + 
    scale_color_manual(values=c("brown1", "green4", "purple3", "cyan3"),
                       name="Shape=Rate Gamma") + 
    scale_y_continuous(breaks=pretty_breaks()) +
    scale_x_continuous(breaks=pretty_breaks()) +
    labs(title=NULL, x="AC (recombination rates)", y="CV (recombination rates)") +
    theme(strip.text.y=element_text(size=16),       
          strip.text.x=element_text(size=16),        
          axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          axis.text.x=element_text(size=12),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16),
          legend.position="bottom")
  save_plot("2D_plot_rec.svg", q, base_height=12, base_width=18)
}