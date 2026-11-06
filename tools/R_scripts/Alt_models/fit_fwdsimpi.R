
#################
# 
# This script uses prepared tables to build SEMs 
# NOTE: this script should be run from within each current_models dir
#
#################

source("~/Data/bgs_lmr/human_data/tools/HelperFunctions.R")

setwd("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/10Mb/fit_Ne_predictions_good_DFE/") # adjust
num_chrom <- 100 # adjust

#################
#
# Loading data
#
#################

sim_maps_10kb <- fread("maps/maps_10kb.csv.gz")
sim_maps_100kb <- fread("maps/maps_100kb.csv.gz") 
sim_maps_1Mb <- fread("maps/maps_1Mb.csv.gz") 

#################
#
# Specifying causal models (DAGs)
#
#################

# these are the four candidate model (fit by chr) in lavaan format
model_1 <- "
    avg_pi ~ start(0.25)*c1*avg_mut
    avg_pi ~ start(0.5)*c2*B
    B ~ start(-0.25)*c3*avg_mut
    "

model_2 <- "
    avg_pi ~ start(0.25)*c1*avg_mut
    avg_pi ~ start(0.5)*c2*B
    avg_pi ~ start(-0.25)*c3*prop_del_sites
    B ~ start(-0.25)*c4*prop_del_sites
    B ~ start(-0.25)*c5*avg_mut
    mu_B_pi:=c1*c5
    "

model_3 <- "
    avg_pi ~ start(0.25)*c1*avg_mut
    avg_pi ~ start(0.5)*c2*B
    avg_pi ~ start(-0.25)*c3*prop_del_sites
    B ~ start(-0.25)*c4*prop_del_sites
    B ~ start(-0.25)*c5*avg_mut
    B ~ start(0.25)*c6*avg_rec
    avg_mut ~ start(0.25)*c7*avg_rec
    avg_mut ~~ 0*prop_del_sites
    mu_B_pi:=c1*c5
    rec_B_pi:=c2*c6
    rec_mu_pi:=c1*c7
    rec_mu_B_pi:=c2*c5*c7
    "

model_4 <- "
    avg_pi ~ start(0.25)*c1*avg_mut
    avg_pi ~ start(0.5)*c2*B
    avg_pi ~ start(-0.25)*c3*prop_del_sites
    B ~ start(-0.25)*c4*prop_del_sites
    B ~ start(-0.25)*c5*avg_mut
    B ~ start(0.25)*c6*avg_rec
    avg_mut ~ start(0.25)*c7*avg_rec
    avg_mut ~~ prop_del_sites
    mu_B_pi:=c1*c5
    rec_B_pi:=c2*c6
    rec_mu_pi:=c1*c7
    rec_mu_B_pi:=c2*c5*c7
    "

#################
#
# Fitting models 
#
#################

################## 10 kb #################

tbl_10kb <- sim_maps_10kb
table_10kb_g <- tbl_10kb[,1:3]
table_10kb_g[, avg_pi := log(tbl_10kb$avg_pi)]
table_10kb_g[, avg_rec := log(tbl_10kb$avg_rec)] 
table_10kb_g[, avg_mut := log(tbl_10kb$avg_mut)] 
table_10kb_g[, B := log(tbl_10kb$B)] 
table_10kb_g[, prop_del_sites := log1p(tbl_10kb$prop_del_sites)]

table_10kb_g <- add_genome_offsets_and_midpoint(table_10kb_g, chrom_lengths=NULL)
dt <- copy(table_10kb_g) %>% arrange(chrom, midpoint) %>% as.data.table() %>% setDT()

# simple correlations
data <- dplyr::select(table_10kb_g, avg_pi, prop_del_sites, B, avg_rec, avg_mut)
n <- ncol(data)

cor_mat <- matrix(NA, n, n)
p_mat <- matrix(NA, n, n)
colnames(cor_mat) <- rownames(cor_mat) <- colnames(data)
colnames(p_mat) <- rownames(p_mat) <- colnames(data)

for(i in 1:n) {
  for(j in 1:n) {
    test <- cor.test(data[[i]], data[[j]], method = "pearson", use = "complete.obs")
    cor_mat[i, j] <- test$estimate
    p_mat[i, j] <- test$p.value
  }
}

combined_mat <- matrix(NA_real_, n, n)
for(i in 1:n) {
  for(j in 1:n) {
    if(i > j) {
      combined_mat[i, j] <- round(cor_mat[i, j], 4)
    } else if (i < j) {
      combined_mat[i, j] <- signif(p_mat[i, j], 3)
    }
  }
}

combined_mat <- setDT(as.data.frame(combined_mat))
rownames(combined_mat) <- colnames(combined_mat) <- colnames(data)
fwrite(combined_mat, paste("sem/pearson_gw_10kb.csv", sep=""), row.names=T, quote=F)

# model 3
edges <- list(
  avg_pi = list(
    list(p = "avg_mut", label = "c1", start = 0.25),
    list(p = "B", label = "c2", start = 0.5),
    list(p = "prop_del_sites",label = "c3", start = -0.25)
  ),
  B = list(
    list(p = "prop_del_sites",label = "c4", start = -0.25),
    list(p = "avg_mut", label = "c5", start = -0.25),
    list(p = "avg_rec", label = "c6", start = 0.25)
  ),
  avg_mut = list(
    list(p = "avg_rec", label = "c7", start = 0.25)
  ),
  avg_rec = list(),
  prop_del_sites = list()
)

covs <- list()

defined <- c(
  mu_B_pi = "c1*c5",
  rec_B_pi = "c2*c6",
  rec_mu_pi = "c1*c7",
  rec_mu_B_pi  = "c2*c5*c7"
)

model3_constrained <- make_lavaan_model(edges,
                                        covs = covs,
                                        defined = defined,
                                        groups = num_chrom,
                                        free_by = "groups",
                                        label_prefix = "x")

model3_constrained <- paste0(model3_constrained, "\n", "avg_mut ~~ 0*prop_del_sites")

out <- preprocess_sem_from_lavaan(model3_constrained,
                                  df = dt,
                                  taper_radius = 5e5,
                                  pos = "midpoint",
                                  chrom = "chrom")
full_df <- out$full_df
sem_df_std <- out$sem_df_std
sem_df_std$chrom <- full_df$chrom
sem_df_std$chrom <- as.factor(sem_df_std$chrom)

m3_gw <- sem(model3_constrained,
             data = sem_df_std,
             estimator = "ML",
             meanstructure = FALSE,
             std.ov = TRUE,
             fixed.x = TRUE,
             optim.method = "default",
             group="chrom",
             group.equal=c("regressions", "residuals"),
             control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

if(lavInspect(m3_gw, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m3_gw, "r2"))),
                          as.numeric(unlist(inspect(m3_gw, "r2")))),
         "sem/lavaan_out/fwd_simpi_DAG3_r2_10kb_genomewide.csv",
         col.names=F)
  
  fwrite(parameterestimates(m3_gw, ci=T),
         "sem/lavaan_out/fwd_simpi_DAG3_params_10kb_genomewide.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m3_gw)),
                          as.numeric(fitMeasures(m3_gw))),
         "sem/lavaan_out/fwd_simpi_DAG3_fitmeasures_10kb_genomewide.csv",
         col.names=F)
} else {
  print("Model m3_gw did not converge (10 kb)")
}

# build a per_group_map that marks every lhs~rhs as group-specific
per_group_map <- character()
for (lhs in names(edges)) {
  for (r in edges[[lhs]]) {
    rhs <- if (is.list(r)) r$p else r
    per_group_map[paste0(lhs, "~", rhs)] <- "groups"
  }
}

model_free_with_group_defs <- make_lavaan_model_free_with_grouped_defined(
  edges = edges, covs = covs, defined = defined,
  groups = num_chrom, per_group_map = per_group_map, label_prefix = "x"
)

model_free_with_group_defs <- paste0(model_free_with_group_defs, "\n", "avg_mut ~~ 0*prop_del_sites")

out <- preprocess_sem_from_lavaan(model_free_with_group_defs,
                                  df = dt,
                                  taper_radius = 5e5,
                                  pos = "midpoint",
                                  chrom = "chrom")
full_df <- out$full_df
sem_df_std <- out$sem_df_std
sem_df_std$chrom <- full_df$chrom
sem_df_std$chrom <- as.factor(sem_df_std$chrom)

m3_chr <- sem(model_free_with_group_defs,
              data = sem_df_std,
              estimator = "ML",
              meanstructure = FALSE,
              std.ov = TRUE,
              fixed.x = TRUE,
              optim.method = "default",
              group="chrom",
              group.equal="none",
              control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

lrt_res <- tryCatch(
  lavTestLRT(m3_gw, m3_chr, method = "satorra.bentler.2001"),
  error = function(e) {
    message("lavTestLRT failed: ", e$message)
    NULL
  }
)

out_path <- "sem/lavaan_out/fwd_simpi_LRT_DAG3_gw-equal-freechrom_10kb.csv"

if(!is.null(lrt_res)) {
  df <- as.data.table(lrt_res, keep.rownames = "comparison")
  num_cols <- setdiff(names(df), "comparison")
  df[, (num_cols) := lapply(.SD, function(x) {
    xn <- as.numeric(x)
    if (all(is.na(xn)) && is.character(x)) return(x) else return(xn)
  }), .SDcols = num_cols]
  fwrite(df, out_path)
  message("Wrote lavTestLRT results to: ", out_path)
} else {
  # Fallback: AIC/BIC 
  aic_free <- AIC(m3_chr); aic_constr <- AIC(m3_gw)
  bic_free <- BIC(m3_chr); bic_constr <- BIC(m3_gw)
  df_free <- tryCatch(fitMeasures(m3_chr, "df"), error = function(e) NA)
  df_constr<- tryCatch(fitMeasures(m3_gw, "df"), error = function(e) NA)
  nobs_free <- tryCatch(fitMeasures(m3_chr, "nobs"), error = function(e) NA)
  nobs_constr <- tryCatch(fitMeasures(m3_gw, "nobs"), error = function(e) NA)
  # approximate p-value with AIC
  delta_AIC <- aic_free - aic_constr  
  LR <- 2 * df_constr - delta_AIC
  p_value <- pchisq(LR, df = df_constr, lower.tail = FALSE)
  
  fallback_df <- data.table(
    comparison = c("free", "constrained"),
    AIC = c(aic_free, aic_constr),
    BIC = c(bic_free, bic_constr),
    df = c(df_free, df_constr),
    nobs = c(nobs_free, nobs_constr),
    `Pr(>Chisq)` = p_value
  )
  fwrite(fallback_df, out_path)
  message("lavTestLRT failed; wrote AIC/BIC fallback to: ", out_path)
}

if(lavInspect(m3_chr, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m3_chr, "r2"))),
                          as.numeric(unlist(inspect(m3_chr, "r2")))),
         "sem/lavaan_out/fwd_simpi_DAG3_r2_10kb_genomewide_freechrom.csv",
         col.names=F)
  
  fwrite(parameterestimates(m3_chr, ci=T),
         "sem/lavaan_out/fwd_simpi_DAG3_params_10kb_genomewide_freechrom.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m3_chr)),
                          as.numeric(fitMeasures(m3_chr))),
         "sem/lavaan_out/fwd_simpi_DAG3_fitmeasures_10kb_genomewide_freechrom.csv",
         col.names=F)

} else {
  print("Model m3_chr did not converge (10 kb)")
}

# model 4
edges <- list(
  avg_pi = list(
    list(p = "avg_mut", label = "c1", start = 0.25),
    list(p = "B", label = "c2", start = 0.5),
    list(p = "prop_del_sites", label = "c3", start = -0.25)
  ),
  B = list(
    list(p = "prop_del_sites", label = "c4", start = -0.25),
    list(p = "avg_mut", label = "c5", start = -0.25),
    list(p = "avg_rec", label = "c6", start = 0.25)
  ),
  avg_mut = list(
    list(p = "avg_rec", label = "c7", start = 0.25)
  ),
  avg_rec = list(),
  prop_del_sites = list()
)

covs <- list(
  c("avg_mut", "prop_del_sites")
)

defined <- c(
  mu_B_pi = "c1*c5",
  rec_B_pi = "c2*c6",
  rec_mu_pi = "c1*c7",
  rec_mu_B_pi = "c2*c5*c7"
)

model4_constrained <- make_lavaan_model(edges,
                                        covs = covs,
                                        defined = defined,
                                        groups = num_chrom,
                                        free_by = "groups",
                                        label_prefix = "x")

out <- preprocess_sem_from_lavaan(model4_constrained,
                                  df = dt,
                                  taper_radius = 5e5,
                                  pos = "midpoint",
                                  chrom = "chrom")
full_df <- out$full_df
sem_df_std <- out$sem_df_std
sem_df_std$chrom <- full_df$chrom
sem_df_std$chrom <- as.factor(sem_df_std$chrom)

m4_gw <- sem(model4_constrained,
             data = sem_df_std,
             estimator = "ML",
             meanstructure = FALSE,
             std.ov = TRUE,
             fixed.x = TRUE,
             optim.method = "default",
             group="chrom",
             group.equal=c("regressions", "residuals"),
             control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

if(lavInspect(m4_gw, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m4_gw, "r2"))),
                          as.numeric(unlist(inspect(m4_gw, "r2")))),
         "sem/lavaan_out/fwd_simpi_DAG4_r2_10kb_genomewide.csv",
         col.names=F)
  
  fwrite(parameterestimates(m4_gw, ci=T),
         "sem/lavaan_out/fwd_simpi_DAG4_params_10kb_genomewide.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m4_gw)),
                          as.numeric(fitMeasures(m4_gw))),
         "sem/lavaan_out/fwd_simpi_DAG4_fitmeasures_10kb_genomewide.csv",
         col.names=F)
  
} else {
  print("Model m4_gw did not converge (10 kb)")
}

# build a per_group_map that marks every lhs~rhs as group-specific
per_group_map <- character()
for (lhs in names(edges)) {
  for (r in edges[[lhs]]) {
    rhs <- if (is.list(r)) r$p else r
    per_group_map[paste0(lhs, "~", rhs)] <- "groups"
  }
}

model_free_with_group_defs <- make_lavaan_model_free_with_grouped_defined(
  edges = edges, covs = covs, defined = defined,
  groups = num_chrom, per_group_map = per_group_map, label_prefix = "x"
)

out <- preprocess_sem_from_lavaan(model_free_with_group_defs,
                                  df = dt,
                                  taper_radius = 5e5,
                                  pos = "midpoint",
                                  chrom = "chrom")
full_df <- out$full_df
sem_df_std <- out$sem_df_std
sem_df_std$chrom <- full_df$chrom
sem_df_std$chrom <- as.factor(sem_df_std$chrom)

m4_chr <- sem(model_free_with_group_defs,
              data = sem_df_std,
              estimator = "ML",
              meanstructure = FALSE,
              std.ov = TRUE,
              fixed.x = TRUE,
              optim.method = "default",
              group="chrom",
              group.equal="none",
              control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

lrt_res <- tryCatch(
  lavTestLRT(m4_gw, m4_chr, method = "satorra.bentler.2001"),
  error = function(e) {
    message("lavTestLRT failed: ", e$message)
    NULL
  }
)

out_path <- "sem/lavaan_out/fwd_simpi_LRT_DAG4_gw-equal-freechrom_10kb.csv"

if (!is.null(lrt_res)) {
  df <- as.data.table(lrt_res, keep.rownames = "comparison")
  num_cols <- setdiff(names(df), "comparison")
  df[, (num_cols) := lapply(.SD, function(x) {
    xn <- as.numeric(x)
    if (all(is.na(xn)) && is.character(x)) return(x) else return(xn)
  }), .SDcols = num_cols]
  fwrite(df, out_path)
  message("Wrote lavTestLRT results to: ", out_path)
} else {
  # Fallback: AIC/BIC and sample-size-adjusted comparisons
  aic_free <- AIC(m4_chr); aic_constr <- AIC(m4_gw)
  bic_free <- BIC(m4_chr); bic_constr <- BIC(m4_gw)
  df_free <- tryCatch(fitMeasures(m4_chr, "df"), error = function(e) NA)
  df_constr<- tryCatch(fitMeasures(m4_gw, "df"), error = function(e) NA)
  nobs_free <- tryCatch(fitMeasures(m4_chr, "nobs"), error = function(e) NA)
  nobs_constr <- tryCatch(fitMeasures(m4_gw, "nobs"), error = function(e) NA)
  # approximate p-value with AIC
  delta_AIC <- aic_free - aic_constr  
  LR <- 2 * df_constr - delta_AIC
  p_value <- pchisq(LR, df = df_constr, lower.tail = FALSE)
  
  fallback_df <- data.table(
    comparison = c("free", "constrained"),
    AIC = c(aic_free, aic_constr),
    BIC = c(bic_free, bic_constr),
    df = c(df_free, df_constr),
    nobs = c(nobs_free, nobs_constr),
    `Pr(>Chisq)` = p_value
  )
  fwrite(fallback_df, out_path)
  message("lavTestLRT failed; wrote AIC/BIC fallback to: ", out_path)
}

if(lavInspect(m4_chr, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m4_chr, "r2"))),
                          as.numeric(unlist(inspect(m4_chr, "r2")))),
         "sem/lavaan_out/fwd_simpi_DAG4_r2_10kb_genomewide_freechrom.csv",
         col.names=F)
  
  fwrite(parameterestimates(m4_chr, ci=T),
         "sem/lavaan_out/fwd_simpi_DAG4_params_10kb_genomewide_freechrom.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m4_chr)),
                          as.numeric(fitMeasures(m4_chr))),
        "sem/lavaan_out/fwd_simpi_DAG4_fitmeasures_10kb_genomewide_freechrom.csv",
         col.names=F)
} else {
  print("Model m4_chr did not converge (10 kb)")
}

# one model per chrom
for(c in 0:(num_chrom-1)) {
  
  m3_fit <- sem(model_3,
                data=filter(dt, chrom==c),
                estimator = "ML",
                meanstructure = FALSE,
                std.ov = TRUE,
                fixed.x = TRUE,
                optim.method = "default",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m3_fit, "converged")) {
    fwrite(rbind.data.frame(names(inspect(m3_fit, "r2")),
                            as.numeric(inspect(m3_fit, "r2"))),
           paste("sem/lavaan_out/fwd_simpi_DAG3_r2_10kb_chr",
                 c, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m3_fit, ci=T),
           paste("sem/lavaan_out/fwd_simpi_DAG3_params_10kb_chr",
                 c, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m3_fit)),
                            as.numeric(fitMeasures(m3_fit))),
           paste("sem/lavaan_out/fwd_simpi_DAG3_fitmeasures_10kb_chr",
                 c, ".csv", sep=""),
           col.names=F)
  } else {
    print(paste("Model 3 did not converge for chr", c, "(10 kb)"))
  }
  
  m4_fit <- sem(model_4,
                data=filter(dt, chrom==c),
                estimator = "ML",
                meanstructure = FALSE,
                std.ov = TRUE,
                fixed.x = TRUE,
                optim.method = "default",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m4_fit, "converged")) {
    fwrite(rbind.data.frame(names(inspect(m4_fit, "r2")),
                            as.numeric(inspect(m4_fit, "r2"))),
           paste("sem/lavaan_out/fwd_simpi_DAG4_r2_10kb_chr",
                 c, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m4_fit, ci=T),
           paste("sem/lavaan_out/fwd_simpi_DAG4_params_10kb_chr",
                 c, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m4_fit)),
                            as.numeric(fitMeasures(m4_fit))),
           paste("sem/lavaan_out/fwd_simpi_DAG4_fitmeasures_10kb_chr",
                 c, ".csv", sep=""),
           col.names=F)
    
    fwrite(lavTestLRT(m3_fit, m4_fit),
           paste("sem/lavaan_out/fwd_simpi_LRT_DAG3_DAG4_10kb_chr",
                 c, ".csv", sep=""))
  } else {
    print(paste("Model 4 did not converge for chr", c, "(10 kb)"))
  }
}

################## bootstrapping 10 kb #################

nboot <- 5e3
chr_list <- vector("list", length=num_chrom)
for(c in 0:(num_chrom-1)) {
  
  cat(paste("chrom", c + 1, "\n"))
  
  cXX <- dplyr::filter(dt, chrom==c)
  cXX_r2_m3 <- numeric(length=nboot)
  
  pb <- txtProgressBar(min=0, max=nboot, style=3)
  for(b in 1:nboot) {
    
    setTxtProgressBar(pb, b)
    cXX_boot <- cXX %>% dplyr::slice_sample(prop=1, replace=TRUE)
    
    m3_fit <- sem(model_3,
                  data=cXX_boot,
                  estimator = "ML",
                  meanstructure = FALSE,
                  std.ov = TRUE,
                  fixed.x = TRUE,
                  optim.method = "default",
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    cXX_r2_m3[b] <- inspect(m3_fit, "r2")[1]
  }
  close(pb)
  
  cXX_boot_r2 <- cbind.data.frame(cXX_r2_m3)
  names(cXX_boot_r2) <- "model_3"
  
  cXX_boot_r2$boot <- 1:nrow(cXX_boot_r2)
  cXX_boot_r2$chrom <- c+1
  #cXX_boot_r2_m <- pivot_longer(cXX_boot_r2, cols=starts_with("model"), names_to="model", values_to="r-sqrd")
  chr_list[[c+1]] <- cXX_boot_r2
}

boot_r2 <- data.table::rbindlist(chr_list)
names(boot_r2)[1] <- "r-sqrd"
boot_r2$chrom <- factor(boot_r2$chrom)

b_10kb <- ggplot(boot_r2, aes(x = chrom, y = `r-sqrd`)) +
  geom_boxplot(alpha = 0.8) + 
  labs(x = "Chromosome", y = expression("Bootstrapped " * R[pi]^2)) + theme_classic() +
  scale_y_continuous(breaks=pretty_breaks(), limits=c(0, 1)) +
  theme(strip.text=element_text(size=16),
        axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_blank(),
        legend.position="none")
save_plot("sem/boot_r2_chroms_10kb.pdf", b_10kb, base_height=8, base_width=17)

################# 100 kb #################

tbl_100kb <- sim_maps_100kb
table_100kb_g <- tbl_100kb[,1:3]
table_100kb_g[, avg_pi := log(tbl_100kb$avg_pi)]
table_100kb_g[, avg_rec := log(tbl_100kb$avg_rec)] 
table_100kb_g[, avg_mut := log(tbl_100kb$avg_mut)] 
table_100kb_g[, B := log(tbl_100kb$B)] 
table_100kb_g[, prop_del_sites := log1p(tbl_100kb$prop_del_sites)]

table_100kb_g <- add_genome_offsets_and_midpoint(table_100kb_g, chrom_lengths=NULL)
dt <- copy(table_100kb_g) %>% arrange(chrom, midpoint) %>% as.data.table() %>% setDT()

data <- dplyr::select(table_100kb_g, avg_pi, prop_del_sites, B, avg_rec, avg_mut)
n <- ncol(data)

cor_mat <- matrix(NA, n, n)
p_mat <- matrix(NA, n, n)
colnames(cor_mat) <- rownames(cor_mat) <- colnames(data)
colnames(p_mat) <- rownames(p_mat) <- colnames(data)

for(i in 1:n) {
  for(j in 1:n) {
    test <- cor.test(data[[i]], data[[j]], method = "pearson", use = "complete.obs")
    cor_mat[i, j] <- test$estimate
    p_mat[i, j] <- test$p.value
  }
}

combined_mat <- matrix(NA_real_, n, n)
for(i in 1:n) {
  for(j in 1:n) {
    if(i > j) {
      combined_mat[i, j] <- round(cor_mat[i, j], 4)
    } else if (i < j) {
      combined_mat[i, j] <- signif(p_mat[i, j], 3)
    }
  }
}
combined_mat <- setDT(as.data.frame(combined_mat))
rownames(combined_mat) <- colnames(combined_mat) <- colnames(data)
fwrite(combined_mat, paste("sem/pearson_gw_100kb.csv", sep=""), row.names=T, quote=F)

# model 3
edges <- list(
  avg_pi = list(
    list(p = "avg_mut", label = "c1", start = 0.25),
    list(p = "B", label = "c2", start = 0.5),
    list(p = "prop_del_sites", label = "c3", start = -0.25)
  ),
  B = list(
    list(p = "prop_del_sites",label = "c4", start = -0.25),
    list(p = "avg_mut", label = "c5", start = -0.25),
    list(p = "avg_rec", label = "c6", start = 0.25)
  ),
  avg_mut = list(
    list(p = "avg_rec", label = "c7", start = 0.25)
  ),
  avg_rec = list(),
  prop_del_sites = list()
)

covs <- list()

defined <- c(
  mu_B_pi = "c1*c5",
  rec_B_pi = "c2*c6",
  rec_mu_pi = "c1*c7",
  rec_mu_B_pi = "c2*c5*c7"
)

model3_constrained <- make_lavaan_model(edges,
                                        covs = covs,
                                        defined = defined,
                                        groups = num_chrom,
                                        free_by = "groups",
                                        label_prefix = "x")

model3_constrained <- paste0(model3_constrained, "\n", "avg_mut ~~ 0*prop_del_sites")

m3_gw <- sem(model3_constrained,
             data = dt,         
             estimator = "ML",              
             meanstructure = FALSE,            
             std.ov = TRUE,
             fixed.x = TRUE,
             optim.method = "default",         
             group="chrom",
             group.equal = c("regressions", "residuals"),
             control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

if(lavInspect(m3_gw, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m3_gw, "r2"))),
                          as.numeric(unlist(inspect(m3_gw, "r2")))), 
         "sem/lavaan_out/fwd_simpi_DAG3_r2_100kb_genomewide.csv",
         col.names=F)
  
  fwrite(parameterestimates(m3_gw, ci=T), 
         "sem/lavaan_out/fwd_simpi_DAG3_params_100kb_genomewide.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m3_gw)),
                          as.numeric(fitMeasures(m3_gw))), 
         "sem/lavaan_out/fwd_simpi_DAG3_fitmeasures_100kb_genomewide.csv",
         col.names=F)
} else {
  print("Model m3_gw did not converge (100 kb)")
}

# build a per_group_map that marks every lhs~rhs as group-specific
per_group_map <- character()
for (lhs in names(edges)) {
  for (r in edges[[lhs]]) {
    rhs <- if (is.list(r)) r$p else r
    per_group_map[paste0(lhs, "~", rhs)] <- "groups"
  }
}

model_free_with_group_defs <- make_lavaan_model_free_with_grouped_defined(
  edges = edges, covs = covs, defined = defined,
  groups = num_chrom, per_group_map = per_group_map, label_prefix = "x"
)

model_free_with_group_defs <- paste0(model_free_with_group_defs, "\n", "avg_mut ~~ 0*prop_del_sites")

m3_chr <- sem(model_free_with_group_defs,
              data = dt,         
              estimator = "ML",              
              meanstructure = FALSE,            
              std.ov = TRUE,    
              fixed.x = TRUE,
              optim.method = "default",         
              group="chrom",
              group.equal="none",
              control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

lrt_res <- tryCatch(
  lavTestLRT(m3_gw, m3_chr, method = "satorra.bentler.2001"),
  error = function(e) {
    message("lavTestLRT failed: ", e$message)
    NULL
  }
)

out_path <- "sem/lavaan_out/fwd_simpi_LRT_DAG3_gw-equal-freechrom_100kb.csv"

if(!is.null(lrt_res)) {
  df <- as.data.table(lrt_res, keep.rownames = "comparison")
  num_cols <- setdiff(names(df), "comparison")
  df[, (num_cols) := lapply(.SD, function(x) {
    xn <- as.numeric(x)
    if (all(is.na(xn)) && is.character(x)) return(x) else return(xn)
  }), .SDcols = num_cols]
  fwrite(df, out_path)
  message("Wrote lavTestLRT results to: ", out_path)
} else {
  # Fallback: AIC/BIC and sample-size-adjusted comparisons
  aic_free <- AIC(m3_chr); aic_constr <- AIC(m3_gw)
  bic_free <- BIC(m3_chr); bic_constr <- BIC(m3_gw)
  df_free  <- tryCatch(fitMeasures(m3_chr, "df"), error = function(e) NA)
  df_constr<- tryCatch(fitMeasures(m3_gw, "df"), error = function(e) NA)
  nobs_free<- tryCatch(fitMeasures(m3_chr, "nobs"), error = function(e) NA)
  nobs_constr<- tryCatch(fitMeasures(m3_gw, "nobs"), error = function(e) NA)
  # approximate p-value with AIC
  delta_AIC <- aic_free - aic_constr  
  LR <- 2 * df_constr - delta_AIC
  p_value <- pchisq(LR, df = df_constr, lower.tail = FALSE)
  
  fallback_df <- data.table(
    comparison = c("free", "constrained"),
    AIC = c(aic_free, aic_constr),
    BIC = c(bic_free, bic_constr),
    df = c(df_free, df_constr),
    nobs = c(nobs_free, nobs_constr),
    `Pr(>Chisq)` = p_value
  )
  fwrite(fallback_df, out_path)
  message("lavTestLRT failed; wrote AIC/BIC fallback to: ", out_path)
}

if(lavInspect(m3_chr, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m3_chr, "r2"))),
                          as.numeric(unlist(inspect(m3_chr, "r2")))), 
         "sem/lavaan_out/fwd_simpi_DAG3_r2_100kb_genomewide_freechrom.csv",
         col.names=F)
  
  fwrite(parameterestimates(m3_chr, ci=T), 
         "sem/lavaan_out/fwd_simpi_DAG3_params_100kb_genomewide_freechrom.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m3_chr)),
                          as.numeric(fitMeasures(m3_chr))), 
         "sem/lavaan_out/fwd_simpi_DAG3_fitmeasures_100kb_genomewide_freechrom.csv",
         col.names=F)
} else {
  print("Model m3_chr did not converge (100 kb)")
}

# model 4
edges <- list(
  avg_pi = list(
    list(p = "avg_mut", label = "c1", start = 0.25),
    list(p = "B", label = "c2", start = 0.5),
    list(p = "prop_del_sites", label = "c3", start = -0.25)
  ),
  B = list(
    list(p = "prop_del_sites", label = "c4", start = -0.25),
    list(p = "avg_mut", label = "c5", start = -0.25),
    list(p = "avg_rec", label = "c6", start = 0.25)
  ),
  avg_mut = list(
    list(p = "avg_rec", label = "c7", start = 0.25)
  ),
  avg_rec = list(),
  prop_del_sites = list()
)

covs <- list(
  c("avg_mut", "prop_del_sites")
)

defined <- c(
  mu_B_pi = "c1*c5",
  rec_B_pi = "c2*c6",
  rec_mu_pi = "c1*c7",
  rec_mu_B_pi = "c2*c5*c7"
)

model4_constrained <- make_lavaan_model(edges,
                                        covs = covs,
                                        defined = defined,
                                        groups = num_chrom,
                                        free_by = "groups",
                                        label_prefix = "x")

m4_gw <- sem(model4_constrained,
             data = dt,         
             estimator = "ML",              
             meanstructure = FALSE,            
             std.ov = TRUE,
             fixed.x = TRUE,
             optim.method = "default",         
             group="chrom",
             group.equal=c("regressions", "residuals"),
             control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

if(lavInspect(m4_gw, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m4_gw, "r2"))),
                          as.numeric(unlist(inspect(m4_gw, "r2")))), 
        "sem/lavaan_out/fwd_simpi_DAG4_r2_100kb_genomewide.csv",
         col.names=F)
  
  fwrite(parameterestimates(m4_gw, ci=T), 
         "sem/lavaan_out/fwd_simpi_DAG4_params_100kb_genomewide.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m4_gw)),
                          as.numeric(fitMeasures(m4_gw))), 
         "sem/lavaan_out/fwd_simpi_DAG4_fitmeasures_100kb_genomewide.csv",
         col.names=F)
  
} else {
  print("Model m4_gw did not converge (100 kb)")
}

# build a per_group_map that marks every lhs~rhs as group-specific
per_group_map <- character()
for(lhs in names(edges)) {
  for(r in edges[[lhs]]) {
    rhs <- if (is.list(r)) r$p else r
    per_group_map[paste0(lhs, "~", rhs)] <- "groups"
  }
}

model_free_with_group_defs <- make_lavaan_model_free_with_grouped_defined(
  edges = edges, covs = covs, defined = defined,
  groups = num_chrom, per_group_map = per_group_map, label_prefix = "x"
)

m4_chr <- sem(model_free_with_group_defs,
              data = dt,         
              estimator = "ML",              
              meanstructure = FALSE,            
              std.ov = TRUE,    
              fixed.x = TRUE,
              optim.method = "default",         
              group="chrom",
              group.equal="none",
              control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

lrt_res <- tryCatch(
  lavTestLRT(m4_gw, m4_chr, method = "satorra.bentler.2001"),
  error = function(e) {
    message("lavTestLRT failed: ", e$message)
    NULL
  }
)

out_path <- "sem/lavaan_out/fwd_simpi_LRT_DAG4_gw-equal-freechrom_100kb.csv"

if(!is.null(lrt_res)) {
  df <- as.data.table(lrt_res, keep.rownames = "comparison")
  num_cols <- setdiff(names(df), "comparison")
  df[, (num_cols) := lapply(.SD, function(x) {
    xn <- as.numeric(x)
    if (all(is.na(xn)) && is.character(x)) return(x) else return(xn)
  }), .SDcols = num_cols]
  fwrite(df, out_path)
  message("Wrote lavTestLRT results to: ", out_path)
} else {
  # Fallback: AIC/BIC and sample-size-adjusted comparisons
  aic_free <- AIC(m4_chr); aic_constr <- AIC(m4_gw)
  bic_free <- BIC(m4_chr); bic_constr <- BIC(m4_gw)
  df_free  <- tryCatch(fitMeasures(m4_chr, "df"), error = function(e) NA)
  df_constr<- tryCatch(fitMeasures(m4_gw, "df"), error = function(e) NA)
  nobs_free<- tryCatch(fitMeasures(m4_chr, "nobs"), error = function(e) NA)
  nobs_constr<- tryCatch(fitMeasures(m4_gw, "nobs"), error = function(e) NA)
  # approximate p-value with AIC
  delta_AIC <- aic_free - aic_constr  
  LR <- 2 * df_constr - delta_AIC
  p_value <- pchisq(LR, df = df_constr, lower.tail = FALSE)
  
  fallback_df <- data.table(
    comparison = c("free", "constrained"),
    AIC = c(aic_free, aic_constr),
    BIC = c(bic_free, bic_constr),
    df = c(df_free, df_constr),
    nobs = c(nobs_free, nobs_constr),
    `Pr(>Chisq)` = p_value
  )
  fwrite(fallback_df, out_path)
  message("lavTestLRT failed; wrote AIC/BIC fallback to: ", out_path)
}

if(lavInspect(m4_chr, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m4_chr, "r2"))),
                          as.numeric(unlist(inspect(m4_chr, "r2")))), 
         "sem/lavaan_out/fwd_simpi_DAG4_r2_100kb_genomewide_freechrom.csv",
         col.names=F)
  
  fwrite(parameterestimates(m4_chr, ci=T), 
         "sem/lavaan_out/fwd_simpi_DAG4_params_100kb_genomewide_freechrom.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m4_chr)),
                          as.numeric(fitMeasures(m4_chr))), 
         "sem/lavaan_out/fwd_simpi_DAG4_fitmeasures_100kb_genomewide_freechrom.csv",
         col.names=F)
} else {
  print("Model m4_chr did not converge (100 kb)")
}

# one model per chrom
for(c in 0:(num_chrom-1)) {
  
  m3_fit <- sem(model_3,
                data=filter(dt, chrom==c),
                estimator = "ML",
                meanstructure = FALSE,
                std.ov = TRUE,
                fixed.x = TRUE,
                optim.method = "default",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m3_fit, "converged")) {
    fwrite(rbind.data.frame(names(inspect(m3_fit, "r2")),
                            as.numeric(inspect(m3_fit, "r2"))),
           paste("sem/lavaan_out/fwd_simpi_DAG3_r2_100kb_chr",
                 c, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m3_fit, ci=T),
           paste("sem/lavaan_out/fwd_simpi_DAG3_params_100kb_chr",
                 c, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m3_fit)),
                            as.numeric(fitMeasures(m3_fit))),
           paste("sem/lavaan_out/fwd_simpi_DAG3_fitmeasures_100kb_chr",
                 c, ".csv", sep=""),
           col.names=F)
  } else {
    print(paste("Model 3 did not converge for chr", c, "(100 kb)"))
  }
  
  m4_fit <- sem(model_4,
                data=filter(dt, chrom==c),
                estimator = "ML",
                meanstructure = FALSE,
                std.ov = TRUE,
                fixed.x = TRUE,
                optim.method = "default",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m4_fit, "converged")) {
    fwrite(rbind.data.frame(names(inspect(m4_fit, "r2")),
                            as.numeric(inspect(m4_fit, "r2"))),
           paste("sem/lavaan_out/fwd_simpi_DAG4_r2_100kb_chr",
                 c, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m4_fit, ci=T),
           paste("sem/lavaan_out/fwd_simpi_DAG4_params_100kb_chr",
                 c, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m4_fit)),
                            as.numeric(fitMeasures(m4_fit))),
           paste("sem/lavaan_out/fwd_simpi_DAG4_fitmeasures_100kb_chr",
                 c, ".csv", sep=""),
           col.names=F)
    
    fwrite(lavTestLRT(m3_fit, m4_fit),
           paste("sem/lavaan_out/fwd_simpi_LRT_DAG3_DAG4_100kb_chr",
                 c, ".csv", sep=""))
  } else {
    print(paste("Model 4 did not converge for chr", c, "(100 kb)"))
  }
}

################## bootstrapping 100 kb ################## 

nboot <- 5e3
chr_list <- vector("list", length=num_chrom)
for(c in 0:(num_chrom-1)) {
  
  cat(paste("chrom", c + 1, "\n"))
  
  cXX <- dplyr::filter(dt, chrom==c)
  cXX_r2_m3 <- numeric(length=nboot)
  
  pb <- txtProgressBar(min=0, max=nboot, style=3)
  for(b in 1:nboot) {
    
    setTxtProgressBar(pb, b)
    cXX_boot <- cXX %>% dplyr::slice_sample(prop=1, replace=TRUE)
    
    m3_fit <- sem(model_3,
                  data=cXX_boot,
                  estimator = "ML",
                  meanstructure = FALSE,
                  std.ov = TRUE,
                  fixed.x = TRUE,
                  optim.method = "default",
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    cXX_r2_m3[b] <- inspect(m3_fit, "r2")[1]
  }
  close(pb)
  
  cXX_boot_r2 <- cbind.data.frame(cXX_r2_m3)
  names(cXX_boot_r2) <- "model_3"
  
  cXX_boot_r2$boot <- 1:nrow(cXX_boot_r2)
  cXX_boot_r2$chrom <- c+1
  #cXX_boot_r2_m <- pivot_longer(cXX_boot_r2, cols=starts_with("model"), names_to="model", values_to="r-sqrd")
  chr_list[[c+1]] <- cXX_boot_r2
}

boot_r2 <- data.table::rbindlist(chr_list)
names(boot_r2)[1] <- "r-sqrd"
boot_r2$chrom <- factor(boot_r2$chrom)

b_100kb <- ggplot(boot_r2, aes(x = chrom, y = `r-sqrd`)) +
  geom_boxplot(alpha = 0.8) + 
  labs(x=NULL, y = expression("Bootstrapped " * R[pi]^2)) + theme_classic() +
  scale_y_continuous(breaks=pretty_breaks(), limits=c(0, 1)) +
  theme(strip.text=element_text(size=16),
        axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_blank(),
        legend.position="none")
save_plot("sem/boot_r2_chroms_100kb.pdf", b_100kb, base_height=8, base_width=17)

################# 1 Mb #################

tbl_1Mb <- sim_maps_1Mb
table_1Mb_g <- tbl_1Mb[,1:3]
table_1Mb_g[, avg_pi := log(tbl_1Mb$avg_pi)]
table_1Mb_g[, avg_rec := log(tbl_1Mb$avg_rec)] 
table_1Mb_g[, avg_mut := log(tbl_1Mb$avg_mut)] 
table_1Mb_g[, B := log(tbl_1Mb$B)] 
table_1Mb_g[, prop_del_sites := log1p(tbl_1Mb$prop_del_sites)]

table_1Mb_g <- add_genome_offsets_and_midpoint(table_1Mb_g, chrom_lengths=NULL)
dt <- copy(table_1Mb_g) %>% arrange(chrom, midpoint) %>% as.data.table() %>% setDT()

data <- dplyr::select(table_1Mb_g, avg_pi, prop_del_sites, B, avg_rec, avg_mut)
n <- ncol(data)

cor_mat <- matrix(NA, n, n)
p_mat <- matrix(NA, n, n)
colnames(cor_mat) <- rownames(cor_mat) <- colnames(data)
colnames(p_mat) <- rownames(p_mat) <- colnames(data)

for(i in 1:n) {
  for(j in 1:n) {
    test <- cor.test(data[[i]], data[[j]], method = "pearson", use = "complete.obs")
    cor_mat[i, j] <- test$estimate
    p_mat[i, j] <- test$p.value
  }
}

combined_mat <- matrix(NA_real_, n, n)
for(i in 1:n) {
  for(j in 1:n) {
    if(i > j) {
      combined_mat[i, j] <- round(cor_mat[i, j], 4)
    } else if (i < j) {
      combined_mat[i, j] <- signif(p_mat[i, j], 3)
    }
  }
}
combined_mat <- setDT(as.data.frame(combined_mat))
rownames(combined_mat) <- colnames(combined_mat) <- colnames(data)
fwrite(combined_mat, paste("sem/pearson_gw_1Mb.csv", sep=""), row.names=T, quote=F)

# model 3
edges <- list(
  avg_pi = list(
    list(p = "avg_mut", label = "c1", start = 0.25),
    list(p = "B", label = "c2", start = 0.5),
    list(p = "prop_del_sites", label = "c3", start = -0.25)
  ),
  B = list(
    list(p = "prop_del_sites",label = "c4", start = -0.25),
    list(p = "avg_mut", label = "c5", start = -0.25),
    list(p = "avg_rec", label = "c6", start = 0.25)
  ),
  avg_mut = list(
    list(p = "avg_rec", label = "c7", start = 0.25)
  ),
  avg_rec = list(),
  prop_del_sites = list()
)

covs <- list()

defined <- c(
  mu_B_pi = "c1*c5",
  rec_B_pi = "c2*c6",
  rec_mu_pi = "c1*c7",
  rec_mu_B_pi = "c2*c5*c7"
)

model3_constrained <- make_lavaan_model(edges,
                                        covs = covs,
                                        defined = defined,
                                        groups = num_chrom,
                                        free_by = "groups",
                                        label_prefix = "x")

model3_constrained <- paste0(model3_constrained, "\n", "avg_mut ~~ 0*prop_del_sites")

m3_gw <- sem(model3_constrained,
             data = dt,         
             estimator = "ML",              
             meanstructure = FALSE,            
             std.ov = TRUE,
             fixed.x = TRUE,
             optim.method = "default",         
             group="chrom",
             group.equal = c("regressions", "residuals"),
             control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

if(lavInspect(m3_gw, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m3_gw, "r2"))),
                          as.numeric(unlist(inspect(m3_gw, "r2")))), 
         "sem/lavaan_out/fwd_simpi_DAG3_r2_1Mb_genomewide.csv",
         col.names=F)
  
  fwrite(parameterestimates(m3_gw, ci=T), 
         "sem/lavaan_out/fwd_simpi_DAG3_params_1Mb_genomewide.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m3_gw)),
                          as.numeric(fitMeasures(m3_gw))), 
         "sem/lavaan_out/fwd_simpi_DAG3_fitmeasures_1Mb_genomewide.csv",
         col.names=F)
} else {
  print("Model m3_gw did not converge (1 Mb)")
}

# build a per_group_map that marks every lhs~rhs as group-specific
per_group_map <- character()
for (lhs in names(edges)) {
  for (r in edges[[lhs]]) {
    rhs <- if (is.list(r)) r$p else r
    per_group_map[paste0(lhs, "~", rhs)] <- "groups"
  }
}

model_free_with_group_defs <- make_lavaan_model_free_with_grouped_defined(
  edges = edges, covs = covs, defined = defined,
  groups = num_chrom, per_group_map = per_group_map, label_prefix = "x"
)

model_free_with_group_defs <- paste0(model_free_with_group_defs, "\n", "avg_mut ~~ 0*prop_del_sites")

m3_chr <- sem(model_free_with_group_defs,
              data = dt,         
              estimator = "ML",              
              meanstructure = FALSE,            
              std.ov = TRUE,    
              fixed.x = TRUE,
              optim.method = "default",         
              group="chrom",
              group.equal="none",
              control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

lrt_res <- tryCatch(
  lavTestLRT(m3_gw, m3_chr, method = "satorra.bentler.2001"),
  error = function(e) {
    message("lavTestLRT failed: ", e$message)
    NULL
  }
)

out_path <- "sem/lavaan_out/fwd_simpi_LRT_DAG3_gw-equal-freechrom_1Mb.csv"

if(!is.null(lrt_res)) {
  df <- as.data.table(lrt_res, keep.rownames = "comparison")
  num_cols <- setdiff(names(df), "comparison")
  df[, (num_cols) := lapply(.SD, function(x) {
    xn <- as.numeric(x)
    if (all(is.na(xn)) && is.character(x)) return(x) else return(xn)
  }), .SDcols = num_cols]
  fwrite(df, out_path)
  message("Wrote lavTestLRT results to: ", out_path)
} else {
  # Fallback: AIC/BIC and sample-size-adjusted comparisons
  aic_free <- AIC(m3_chr); aic_constr <- AIC(m3_gw)
  bic_free <- BIC(m3_chr); bic_constr <- BIC(m3_gw)
  df_free  <- tryCatch(fitMeasures(m3_chr, "df"), error = function(e) NA)
  df_constr<- tryCatch(fitMeasures(m3_gw, "df"), error = function(e) NA)
  nobs_free<- tryCatch(fitMeasures(m3_chr, "nobs"), error = function(e) NA)
  nobs_constr<- tryCatch(fitMeasures(m3_gw, "nobs"), error = function(e) NA)
  # approximate p-value with AIC
  delta_AIC <- aic_free - aic_constr  
  LR <- 2 * df_constr - delta_AIC
  p_value <- pchisq(LR, df = df_constr, lower.tail = FALSE)
  
  fallback_df <- data.table(
    comparison = c("free", "constrained"),
    AIC = c(aic_free, aic_constr),
    BIC = c(bic_free, bic_constr),
    df = c(df_free, df_constr),
    nobs = c(nobs_free, nobs_constr),
    `Pr(>Chisq)` = p_value
  )
  fwrite(fallback_df, out_path)
  message("lavTestLRT failed; wrote AIC/BIC fallback to: ", out_path)
}

if(lavInspect(m3_chr, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m3_chr, "r2"))),
                          as.numeric(unlist(inspect(m3_chr, "r2")))), 
         "sem/lavaan_out/fwd_simpi_DAG3_r2_1Mb_genomewide_freechrom.csv",
         col.names=F)
  
  fwrite(parameterestimates(m3_chr, ci=T), 
         "sem/lavaan_out/fwd_simpi_DAG3_params_1Mb_genomewide_freechrom.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m3_chr)),
                          as.numeric(fitMeasures(m3_chr))), 
         "sem/lavaan_out/fwd_simpi_DAG3_fitmeasures_1Mb_genomewide_freechrom.csv",
         col.names=F)
} else {
  print("Model m3_chr did not converge (1 Mb)")
}

# model 4
edges <- list(
  avg_pi = list(
    list(p = "avg_mut", label = "c1", start = 0.25),
    list(p = "B", label = "c2", start = 0.5),
    list(p = "prop_del_sites", label = "c3", start = -0.25)
  ),
  B = list(
    list(p = "prop_del_sites", label = "c4", start = -0.25),
    list(p = "avg_mut", label = "c5", start = -0.25),
    list(p = "avg_rec", label = "c6", start = 0.25)
  ),
  avg_mut = list(
    list(p = "avg_rec", label = "c7", start = 0.25)
  ),
  avg_rec = list(),
  prop_del_sites = list()
)

covs <- list(
  c("avg_mut", "prop_del_sites")
)

defined <- c(
  mu_B_pi = "c1*c5",
  rec_B_pi = "c2*c6",
  rec_mu_pi = "c1*c7",
  rec_mu_B_pi = "c2*c5*c7"
)

model4_constrained <- make_lavaan_model(edges,
                                        covs = covs,
                                        defined = defined,
                                        groups = num_chrom,
                                        free_by = "groups",
                                        label_prefix = "x")

m4_gw <- sem(model4_constrained,
             data = dt,         
             estimator = "ML",              
             meanstructure = FALSE,            
             std.ov = TRUE,
             fixed.x = TRUE,
             optim.method = "default",         
             group="chrom",
             group.equal=c("regressions", "residuals"),
             control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

if(lavInspect(m4_gw, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m4_gw, "r2"))),
                          as.numeric(unlist(inspect(m4_gw, "r2")))), 
         "sem/lavaan_out/fwd_simpi_DAG4_r2_1Mb_genomewide.csv",
         col.names=F)
  
  fwrite(parameterestimates(m4_gw, ci=T), 
         "sem/lavaan_out/fwd_simpi_DAG4_params_1Mb_genomewide.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m4_gw)),
                          as.numeric(fitMeasures(m4_gw))), 
         "sem/lavaan_out/fwd_simpi_DAG4_fitmeasures_1Mb_genomewide.csv",
         col.names=F)
  
} else {
  print("Model m4_gw did not converge (1 Mb)")
}

# build a per_group_map that marks every lhs~rhs as group-specific
per_group_map <- character()
for(lhs in names(edges)) {
  for(r in edges[[lhs]]) {
    rhs <- if (is.list(r)) r$p else r
    per_group_map[paste0(lhs, "~", rhs)] <- "groups"
  }
}

model_free_with_group_defs <- make_lavaan_model_free_with_grouped_defined(
  edges = edges, covs = covs, defined = defined,
  groups = num_chrom, per_group_map = per_group_map, label_prefix = "x"
)

m4_chr <- sem(model_free_with_group_defs,
              data = dt,         
              estimator = "ML",              
              meanstructure = FALSE,            
              std.ov = TRUE,    
              fixed.x = TRUE,
              optim.method = "default",         
              group="chrom",
              group.equal="none",
              control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

lrt_res <- tryCatch(
  lavTestLRT(m4_gw, m4_chr, method = "satorra.bentler.2001"),
  error = function(e) {
    message("lavTestLRT failed: ", e$message)
    NULL
  }
)

out_path <- "sem/lavaan_out/fwd_simpi_LRT_DAG4_gw-equal-freechrom_1Mb.csv"

if(!is.null(lrt_res)) {
  df <- as.data.table(lrt_res, keep.rownames = "comparison")
  num_cols <- setdiff(names(df), "comparison")
  df[, (num_cols) := lapply(.SD, function(x) {
    xn <- as.numeric(x)
    if (all(is.na(xn)) && is.character(x)) return(x) else return(xn)
  }), .SDcols = num_cols]
  fwrite(df, out_path)
  message("Wrote lavTestLRT results to: ", out_path)
} else {
  # Fallback: AIC/BIC and sample-size-adjusted comparisons
  aic_free <- AIC(m4_chr); aic_constr <- AIC(m4_gw)
  bic_free <- BIC(m4_chr); bic_constr <- BIC(m4_gw)
  df_free  <- tryCatch(fitMeasures(m4_chr, "df"), error = function(e) NA)
  df_constr<- tryCatch(fitMeasures(m4_gw, "df"), error = function(e) NA)
  nobs_free<- tryCatch(fitMeasures(m4_chr, "nobs"), error = function(e) NA)
  nobs_constr<- tryCatch(fitMeasures(m4_gw, "nobs"), error = function(e) NA)
  # approximate p-value with AIC
  delta_AIC <- aic_free - aic_constr  
  LR <- 2 * df_constr - delta_AIC
  p_value <- pchisq(LR, df = df_constr, lower.tail = FALSE)
  
  fallback_df <- data.table(
    comparison = c("free", "constrained"),
    AIC = c(aic_free, aic_constr),
    BIC = c(bic_free, bic_constr),
    df = c(df_free, df_constr),
    nobs = c(nobs_free, nobs_constr),
    `Pr(>Chisq)` = p_value
  )
  fwrite(fallback_df, out_path)
  message("lavTestLRT failed; wrote AIC/BIC fallback to: ", out_path)
}

if(lavInspect(m4_chr, "converged")) {
  fwrite(rbind.data.frame(names(unlist(inspect(m4_chr, "r2"))),
                          as.numeric(unlist(inspect(m4_chr, "r2")))), 
         "sem/lavaan_out/fwd_simpi_DAG4_r2_1Mb_genomewide_freechrom.csv",
         col.names=F)
  
  fwrite(parameterestimates(m4_chr, ci=T), 
         "sem/lavaan_out/fwd_simpi_DAG4_params_1Mb_genomewide_freechrom.csv")
  
  fwrite(rbind.data.frame(names(fitMeasures(m4_chr)),
                          as.numeric(fitMeasures(m4_chr))), 
         "sem/lavaan_out/fwd_simpi_DAG4_fitmeasures_1Mb_genomewide_freechrom.csv",
         col.names=F)
} else {
  print("Model m4_chr did not converge (1 Mb)")
}

# one model per chrom
for(c in 0:(num_chrom-1)) {
  
  m3_fit <- sem(model_3,
                data=filter(dt, chrom==c),
                estimator = "ML",
                meanstructure = FALSE,
                std.ov = TRUE,
                fixed.x = TRUE,
                optim.method = "default",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m3_fit, "converged")) {
    fwrite(rbind.data.frame(names(inspect(m3_fit, "r2")),
                            as.numeric(inspect(m3_fit, "r2"))),
           paste("sem/lavaan_out/fwd_simpi_DAG3_r2_1Mb_chr",
                 c, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m3_fit, ci=T),
           paste("sem/lavaan_out/fwd_simpi_DAG3_params_1Mb_chr",
                 c, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m3_fit)),
                            as.numeric(fitMeasures(m3_fit))),
           paste("sem/lavaan_out/fwd_simpi_DAG3_fitmeasures_1Mb_chr",
                 c, ".csv", sep=""),
           col.names=F)
  } else {
    print(paste("Model 3 did not converge for chr", c, "(1 Mb)"))
  }
  
  m4_fit <- sem(model_4,
                data=filter(dt, chrom==c),
                estimator = "ML",
                meanstructure = FALSE,
                std.ov = TRUE,
                fixed.x = TRUE,
                optim.method = "default",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m4_fit, "converged")) {
    fwrite(rbind.data.frame(names(inspect(m4_fit, "r2")),
                            as.numeric(inspect(m4_fit, "r2"))),
           paste("sem/lavaan_out/fwd_simpi_DAG4_r2_1Mb_chr",
                 c, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m4_fit, ci=T),
           paste("sem/lavaan_out/fwd_simpi_DAG4_params_1Mb_chr",
                 c, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m4_fit)),
                            as.numeric(fitMeasures(m4_fit))),
           paste("sem/lavaan_out/fwd_simpi_DAG4_fitmeasures_1Mb_chr",
                 c, ".csv", sep=""),
           col.names=F)
    
    fwrite(lavTestLRT(m3_fit, m4_fit),
           paste("sem/lavaan_out/fwd_simpi_LRT_DAG3_DAG4_1Mb_chr",
                 c, ".csv", sep=""))
  } else {
    print(paste("Model 4 did not converge for chr", c, "(1 Mb)"))
  }
}

################## bootstrapping 1 Mb ################## 

nboot <- 5e3
chr_list <- vector("list", length=num_chrom)
for(c in 0:(num_chrom-1)) {
  
  cat(paste("chrom", c + 1, "\n"))
  
  cXX <- dplyr::filter(dt, chrom==c)
  cXX_r2_m3 <- numeric(length=nboot)
  
  pb <- txtProgressBar(min=0, max=nboot, style=3)
  for(b in 1:nboot) {
    
    setTxtProgressBar(pb, b)
    cXX_boot <- cXX %>% dplyr::slice_sample(prop=1, replace=TRUE)
    
    m3_fit <- sem(model_3,
                  data=cXX_boot,
                  estimator = "ML",
                  meanstructure = FALSE,
                  std.ov = TRUE,
                  fixed.x = TRUE,
                  optim.method = "default",
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    cXX_r2_m3[b] <- inspect(m3_fit, "r2")[1]
  }
  close(pb)
  
  cXX_boot_r2 <- cbind.data.frame(cXX_r2_m3)
  names(cXX_boot_r2) <- "model_3"
  
  cXX_boot_r2$boot <- 1:nrow(cXX_boot_r2)
  cXX_boot_r2$chrom <- c+1
  #cXX_boot_r2_m <- pivot_longer(cXX_boot_r2, cols=starts_with("model"), names_to="model", values_to="r-sqrd")
  chr_list[[c+1]] <- cXX_boot_r2
}

boot_r2 <- data.table::rbindlist(chr_list)
names(boot_r2)[1] <- "r-sqrd"
boot_r2$chrom <- factor(boot_r2$chrom)

b_1Mb <- ggplot(boot_r2, aes(x = chrom, y = `r-sqrd`)) +
  geom_boxplot(alpha = 0.8) + 
  labs(x=NULL, y = expression("Bootstrapped " * R[pi]^2)) + theme_classic() +
  scale_y_continuous(breaks=pretty_breaks(), limits=c(0, 1)) +
  theme(strip.text=element_text(size=16),
        axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_blank(),
        legend.position="none")
save_plot("sem/boot_r2_chroms_1Mb.pdf", b_1Mb, base_height=8, base_width=17)

boot_fwdsimpi <- plot_grid(b_1Mb, b_100kb, b_10kb, ncol=1, labels="AUTO")
save_plot("sem/boot_r2_fwd_simpi.pdf", boot_fwdsimpi, base_height=8, base_width=16)

