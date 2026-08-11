
#################
# 
# This script uses prepared tables to build SEMs 
# NOTE: this script should be run from within each models dir
#
#################

source("~/Devel/bgs_lmr/tools/Rscripts/HelperFunctions.R")

constrained_models <- c("merged_cds_phastcons", "split_cds_phastcons", "split_cds_regulatory", "merged_cds_regulatory")
mutation_maps <- c("carlson", "gnomad", "roulette")

#################
#
# Loading data
#
#################

# arbitrarily subsetting to the roulette mut map
# (irrelevant since u and B will be replaced by latent constructs)
human_maps_10kb <- fread("maps/maps_10kb.csv.gz") %>% filter(., mut_map=="roulette")
human_maps_100kb <- fread("maps/maps_100kb.csv.gz") %>% filter(., mut_map=="roulette") 
human_maps_1Mb <- fread("maps/maps_1Mb.csv.gz") %>% filter(., mut_map=="roulette") 

# recombination latent (omitted for now, using pyrho_YRI instead)
# rec_latent_10kb <- fread("maps/rec_latent_10kb.csv.gz")
# rec_latent_100kb <- fread("maps/rec_latent_100kb.csv.gz")
# rec_latent_1Mb <- fread("maps/rec_latent_1Mb.csv.gz")

# mutation latent
mut_latent_10kb <- fread("maps/mut_latent_10kb.csv.gz")
mut_latent_100kb <- fread("maps/mut_latent_100kb.csv.gz")
mut_latent_1Mb <- fread("maps/mut_latent_1Mb.csv.gz")

# B-map latent
B_latent_10kb <- data.frame()
B_latent_100kb <- data.frame()
B_latent_1Mb <- data.frame()

for(m in constrained_models) {
  
  tmp_10kb<- fread(paste0("maps/B_latent_", m, "_10kb.csv.gz"))
  tmp_100kb <- fread(paste0("maps/B_latent_", m, "_100kb.csv.gz"))
  tmp_1Mb <- fread(paste0("maps/B_latent_", m, "_1Mb.csv.gz"))
  
  B_latent_10kb <- rbind.data.frame(B_latent_10kb, tmp_10kb)
  B_latent_100kb <- rbind.data.frame(B_latent_100kb, tmp_100kb)
  B_latent_1Mb <- rbind.data.frame(B_latent_1Mb, tmp_1Mb)
}

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

# whether to assess robustness of r^2 by shuffling windows (computationally expensive)
permut <- F
bootstrap <- F

################## 10 kb #################
for(m in constrained_models) {

  print(paste("10kb", m, Sys.Date(), Sys.time()), sep="\t")

  tmp_B <- filter(B_latent_10kb, constrained_model==m) %>% dplyr::select(., -constrained_model)
  
  tbl_10kb <- filter(human_maps_10kb, constrained_model==m, coverage >= 0.75)
  table_10kb_g <- tbl_10kb[,1:3]
  table_10kb_g[, avg_pi := log(tbl_10kb$avg_pi)]
  table_10kb_g[, avg_rec := log(tbl_10kb$avg_rec)] # pyrho YRI pop.
  table_10kb_g[, prop_del_sites := log1p(tbl_10kb$prop_del_sites)]

  #table_10kb_g <- left_join(table_10kb_g, rec_latent_10kb, by=c("chrom", "chromStart", "chromEnd")) 
  table_10kb_g <- left_join(table_10kb_g, tmp_B, by=c("chrom", "chromStart", "chromEnd"))
  table_10kb_g <- left_join(table_10kb_g, mut_latent_10kb, by=c("chrom", "chromStart", "chromEnd"))
  table_10kb_g <- add_genome_offsets_and_midpoint(table_10kb_g, chrom_lengths=NULL)
  
  setnames(table_10kb_g, old=c("B_latent_log", "mut_latent_log"), new=c("B", "avg_mut"))
  
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
  fwrite(combined_mat, paste("sem/pearson_gw_10kb_", m, ".csv", sep=""), row.names=T, quote=F)

  # model 1
  edges <- list(
    avg_pi = list(
      list(p = "avg_mut", label = "c1", start = 0.25),
      list(p = "B", label = "c2", start = 0.5)
    ),
    B = list(
      list(p = "avg_mut", label = "c3", start = -0.25)
    )
  )

  covs <- list()

  defined <- c(
    direct_mut_pi = "c1",
    indirect_mut_pi = "c3*c2",
    total_mut_pi = "c1 + c3*c2",
    direct_B_pi = "c2"
  )

  model1_constrained <- make_lavaan_model(edges,
                                          covs = covs,
                                          defined = defined,
                                          groups = 22,
                                          free_by = "groups",
                                          label_prefix = "x")

  out <- preprocess_sem_from_lavaan(model1_constrained,
                                    df = dt,
                                    taper_radius = 5e5,
                                    pos = "midpoint",
                                    chrom = "chrom")
  full_df <- out$full_df
  sem_df_std <- out$sem_df_std
  sem_df_std$chrom <- full_df$chrom
  sem_df_std$chrom <- as.factor(sem_df_std$chrom)

  m1_gw <- sem(model1_constrained,
               data = sem_df_std,
               estimator = "MLR",
               missing = "fiml",
               meanstructure = TRUE ,
               std.ov = TRUE,
               fixed.x = FALSE,
               optim.method = "nlminb",
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

  summary(m1_gw, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
  inspect(m1_gw, "coef")

  if(lavInspect(m1_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m1_gw, "r2"))),
                            as.numeric(unlist(inspect(m1_gw, "r2")))),
           paste("sem/lavaan_out/obspi_YRI_DAG1_r2_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""),
           col.names=F)

    fwrite(parameterestimates(m1_gw, ci=T),
           paste("sem/lavaan_out/obspi_YRI_DAG1_params_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""))

    fwrite(rbind.data.frame(names(fitMeasures(m1_gw)),
                            as.numeric(fitMeasures(m1_gw))),
           paste("sem/lavaan_out/obspi_YRI_DAG1_fitmeasures_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""),
           col.names=F)
  } else {
    print("Model m1_gw did not converge (10 kb)")
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
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

  m1_chr <- sem(model_free_with_group_defs,
                data = sem_df_std,
                estimator = "MLR",
                missing = "fiml",
                meanstructure = TRUE ,
                std.ov = TRUE,
                fixed.x = FALSE,
                optim.method = "nlminb",
                group="chrom",
                group.equal="none",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

  lrt_res <- tryCatch(
    lavTestLRT(m1_gw, m1_chr, method = "satorra.bentler.2001"),
    error = function(e) {
      message("lavTestLRT failed: ", e$message)
      NULL
    }
  )

  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG1_gw-equal-freechrom_10kb_latents_", m, ".csv")

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
    aic_free <- AIC(m1_chr); aic_constr <- AIC(m1_gw)
    bic_free <- BIC(m1_chr); bic_constr <- BIC(m1_gw)
    df_free <- tryCatch(fitMeasures(m1_chr, "df"), error = function(e) NA)
    df_constr <- tryCatch(fitMeasures(m1_gw, "df"), error = function(e) NA)
    nobs_free <- tryCatch(fitMeasures(m1_chr, "nobs"), error = function(e) NA)
    nobs_constr <- tryCatch(fitMeasures(m1_gw, "nobs"), error = function(e) NA)
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

  if(lavInspect(m1_chr, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m1_chr, "r2"))),
                            as.numeric(unlist(inspect(m1_chr, "r2")))),
           paste("sem/lavaan_out/obspi_YRI_DAG1_r2_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)

    fwrite(parameterestimates(m1_chr, ci=T),
           paste("sem/lavaan_out/obspi_YRI_DAG1_params_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""))

    fwrite(rbind.data.frame(names(fitMeasures(m1_chr)), as.numeric(fitMeasures(m1_chr))),
           paste("sem/lavaan_out/obspi_YRI_DAG1_fitmeasures_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
  } else {
    print("Model m1_chr did not converge (10 kb)")
  }

  # model 2
  edges <- list(
    avg_pi = list(
      list(p = "avg_mut", label = "c1", start = 0.25),
      list(p = "B", label = "c2", start = 0.5),
      list(p = "prop_del_sites", label = "c3", start = -0.25)
    ),
    B = list(
      list(p = "prop_del_sites", label = "c4", start = -0.25),
      list(p = "avg_mut", label = "c5", start = -0.25)
    ),
    avg_mut = list(),
    prop_del_sites = list()
  )

  covs <- list()

  defined <- c(
    mu_B_pi = "c1*c5"
  )

  model2_constrained <- make_lavaan_model(edges,
                                          covs = covs,
                                          defined = defined,
                                          groups = 22,
                                          free_by = "groups",
                                          label_prefix = "x")

  out <- preprocess_sem_from_lavaan(model2_constrained,
                                    df = dt,
                                    taper_radius = 5e5,
                                    pos = "midpoint",
                                    chrom = "chrom")
  full_df <- out$full_df
  sem_df_std <- out$sem_df_std
  sem_df_std$chrom <- full_df$chrom
  sem_df_std$chrom <- as.factor(sem_df_std$chrom)

  m2_gw <- sem(model2_constrained,
               data = sem_df_std,
               estimator = "MLR",
               missing = "fiml",
               fixed.x = FALSE,
               meanstructure = TRUE ,
               std.ov = TRUE,
               optim.method = "nlminb",
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

  if(lavInspect(m2_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m2_gw, "r2"))),
                            as.numeric(unlist(inspect(m2_gw, "r2")))),
           paste("sem/lavaan_out/obspi_YRI_DAG2_r2_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""),
           col.names=F)

    fwrite(parameterestimates(m2_gw, ci=T),
           paste("sem/lavaan_out/obspi_YRI_DAG2_params_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""))

    fwrite(rbind.data.frame(names(fitMeasures(m2_gw)),
                            as.numeric(fitMeasures(m2_gw))),
           paste("sem/lavaan_out/obspi_YRI_DAG2_fitmeasures_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""),
           col.names=F)
  } else {
    print("Model m2_gw did not converge (10 kb)")
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
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

  m2_chr <- sem(model_free_with_group_defs,
                data = sem_df_std,
                estimator = "MLR",
                missing = "fiml",
                meanstructure = TRUE ,
                std.ov = TRUE,
                fixed.x = FALSE,
                optim.method = "nlminb",
                group="chrom",
                group.equal="none",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

  lrt_res <- tryCatch(
    lavTestLRT(m2_gw, m2_chr, method = "satorra.bentler.2001"),
    error = function(e) {
      message("lavTestLRT failed: ", e$message)
      NULL})

  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG2_gw-equal-freechrom_10kb_latents_", m, ".csv")

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
    aic_free <- AIC(m2_chr); aic_constr <- AIC(m2_gw)
    bic_free <- BIC(m2_chr); bic_constr <- BIC(m2_gw)
    df_free <- tryCatch(fitMeasures(m2_chr, "df"), error = function(e) NA)
    df_constr <- tryCatch(fitMeasures(m2_gw, "df"), error = function(e) NA)
    nobs_free <- tryCatch(fitMeasures(m2_chr, "nobs"), error = function(e) NA)
    nobs_constr <- tryCatch(fitMeasures(m2_gw, "nobs"), error = function(e) NA)
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

  if(lavInspect(m2_chr, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m2_chr, "r2"))), as.numeric(unlist(inspect(m2_chr, "r2")))),
           paste("sem/lavaan_out/obspi_YRI_DAG2_r2_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)

    fwrite(parameterestimates(m2_chr, ci=T),
           paste("sem/lavaan_out/obspi_YRI_DAG2_params_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""))

    fwrite(rbind.data.frame(names(fitMeasures(m2_chr)), as.numeric(fitMeasures(m2_chr))),
           paste("sem/lavaan_out/obspi_YRI_DAG2_fitmeasures_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
  } else {
    print("Model m2_chr did not converge (10 kb)")
  }

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
                                          groups = 22,
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
               estimator = "MLR",
               missing = "fiml",
               meanstructure = TRUE ,
               std.ov = TRUE,
               fixed.x = FALSE,
               optim.method = "nlminb",
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

  if(lavInspect(m3_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m3_gw, "r2"))),
                            as.numeric(unlist(inspect(m3_gw, "r2")))),
           paste("sem/lavaan_out/obspi_YRI_DAG3_r2_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""),
           col.names=F)

    fwrite(parameterestimates(m3_gw, ci=T),
           paste("sem/lavaan_out/obspi_YRI_DAG3_params_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""))

    fwrite(rbind.data.frame(names(fitMeasures(m3_gw)),
                            as.numeric(fitMeasures(m3_gw))),
           paste("sem/lavaan_out/obspi_YRI_DAG3_fitmeasures_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""),
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
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
                estimator = "MLR",
                missing = "fiml",
                meanstructure = TRUE ,
                std.ov = TRUE,
                fixed.x = FALSE,
                optim.method = "nlminb",
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

  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG3_gw-equal-freechrom_10kb_latents_", m, ".csv")

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
           paste("sem/lavaan_out/obspi_YRI_DAG3_r2_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)

    fwrite(parameterestimates(m3_chr, ci=T),
           paste("sem/lavaan_out/obspi_YRI_DAG3_params_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""))

    fwrite(rbind.data.frame(names(fitMeasures(m3_chr)),
                            as.numeric(fitMeasures(m3_chr))),
           paste("sem/lavaan_out/obspi_YRI_DAG3_fitmeasures_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""),
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
                                          groups = 22,
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
               estimator = "MLR",
               missing = "fiml",
               meanstructure = TRUE ,
               std.ov = TRUE,
               fixed.x = FALSE,
               optim.method = "nlminb",
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

  if(lavInspect(m4_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m4_gw, "r2"))),
                            as.numeric(unlist(inspect(m4_gw, "r2")))),
           paste("sem/lavaan_out/obspi_YRI_DAG4_r2_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""),
           col.names=F)

    fwrite(parameterestimates(m4_gw, ci=T),
           paste("sem/lavaan_out/obspi_YRI_DAG4_params_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""))

    fwrite(rbind.data.frame(names(fitMeasures(m4_gw)),
                            as.numeric(fitMeasures(m4_gw))),
           paste("sem/lavaan_out/obspi_YRI_DAG4_fitmeasures_10kb_genomewide_",
                 "latents_", m, ".csv", sep=""),
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
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
                estimator = "MLR",
                missing = "fiml",
                meanstructure = TRUE ,
                std.ov = TRUE,
                fixed.x = FALSE,
                optim.method = "nlminb",
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

  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG4_gw-equal-freechrom_10kb_latents_", m, ".csv")

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
           paste("sem/lavaan_out/obspi_YRI_DAG4_r2_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)

    fwrite(parameterestimates(m4_chr, ci=T),
           paste("sem/lavaan_out/obspi_YRI_DAG4_params_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""))

    fwrite(rbind.data.frame(names(fitMeasures(m4_chr)),
                            as.numeric(fitMeasures(m4_chr))),
           paste("sem/lavaan_out/obspi_YRI_DAG4_fitmeasures_10kb_genomewide_",
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
  } else {
    print("Model m4_chr did not converge (10 kb)")
  }

  print("chromosome-specific models...")

  # one model per chrom
  for(c in 1:22) {

    data <- filter(table_10kb_g, chrom==c) %>%
      dplyr::select(., avg_pi, prop_del_sites, B, avg_rec, avg_mut)
    n <- ncol(data)

    cor_mat <- matrix(NA, n, n)
    p_mat <- matrix(NA, n, n)
    colnames(cor_mat) <- rownames(cor_mat) <- colnames(data)
    colnames(p_mat) <- rownames(p_mat) <- colnames(data)

    for (i in 1:n) {
      for (j in 1:n) {
        test <- cor.test(data[[i]], data[[j]], method = "pearson", use = "complete.obs")
        cor_mat[i, j] <- test$estimate
        p_mat[i, j] <- test$p.value
      }
    }

    combined_mat <- matrix(NA_real_, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i > j) {
          combined_mat[i, j] <- round(cor_mat[i, j], 4)
        } else if (i < j) {
          combined_mat[i, j] <- signif(p_mat[i, j], 3)
        }
      }
    }

    combined_mat <- setDT(as.data.frame(combined_mat))
    rownames(combined_mat) <- colnames(combined_mat) <- colnames(data)
    fwrite(combined_mat, paste("sem/pearson_chr", c, "_", m, "_10kb.csv", sep=""),
           row.names=T, quote=F)

    m1_fit <- sem(model_1,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",
                  missing = "fiml",
                  meanstructure = TRUE ,
                  std.ov = TRUE,
                  fixed.x = TRUE,
                  optim.method = "nlminb",
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

    if(lavInspect(m1_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m1_fit, "r2")),
                              as.numeric(inspect(m1_fit, "r2"))),
             paste("sem/lavaan_out/obspi_YRI_DAG1_r2_10kb_chr",
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)

      fwrite(parameterestimates(m1_fit, ci=T),
             paste("sem/lavaan_out/obspi_YRI_DAG1_params_10kb_chr",
                   c, "_latents_", m, ".csv", sep=""))
    } else {
      print(paste("Model 1 did not converge for chr", c, "(10 kb)"))
    }

    m2_fit <- sem(model_2,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",
                  missing = "fiml",
                  meanstructure = TRUE ,
                  std.ov = TRUE,
                  fixed.x = TRUE,
                  optim.method = "nlminb",
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

    if(lavInspect(m2_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m2_fit, "r2")),
                              as.numeric(inspect(m2_fit, "r2"))),
             paste("sem/lavaan_out/obspi_YRI_DAG2_r2_10kb_chr",
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)

      fwrite(parameterestimates(m2_fit, ci=T),
             paste("sem/lavaan_out/obspi_YRI_DAG2_params_10kb_chr",
                   c, "_latents_", m, ".csv", sep=""))
    } else {
      print(paste("Model 2 did not converge for chr", c, "(10 kb)"))
    }

    m3_fit <- sem(model_3,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",
                  missing = "fiml",
                  meanstructure = TRUE ,
                  std.ov = TRUE,
                  fixed.x = TRUE,
                  optim.method = "nlminb",
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

    if(lavInspect(m3_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m3_fit, "r2")),
                              as.numeric(inspect(m3_fit, "r2"))),
             paste("sem/lavaan_out/obspi_YRI_DAG3_r2_10kb_chr",
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)

      fwrite(parameterestimates(m3_fit, ci=T),
             paste("sem/lavaan_out/obspi_YRI_DAG3_params_10kb_chr",
                   c, "_latents_", m, ".csv", sep=""))

      fwrite(rbind.data.frame(names(fitMeasures(m3_fit)),
                              as.numeric(fitMeasures(m3_fit))),
             paste("sem/lavaan_out/obspi_YRI_DAG3_fitmeasures_10kb_chr",
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
    } else {
      print(paste("Model 3 did not converge for chr", c, "(10 kb)"))
    }

    m4_fit <- sem(model_4,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",
                  missing = "fiml",
                  meanstructure = TRUE ,
                  std.ov = TRUE,
                  fixed.x = TRUE,
                  optim.method = "nlminb",
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

    if(lavInspect(m4_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m4_fit, "r2")),
                              as.numeric(inspect(m4_fit, "r2"))),
             paste("sem/lavaan_out/obspi_YRI_DAG4_r2_10kb_chr",
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)

      fwrite(parameterestimates(m4_fit, ci=T),
             paste("sem/lavaan_out/obspi_YRI_DAG4_params_10kb_chr",
                   c, "_latents_", m, ".csv", sep=""))

      fwrite(rbind.data.frame(names(fitMeasures(m4_fit)),
                              as.numeric(fitMeasures(m4_fit))),
             paste("sem/lavaan_out/obspi_YRI_DAG4_fitmeasures_10kb_chr",
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)

      fwrite(safeLRT(m3_fit, m4_fit),
             paste("sem/lavaan_out/obspi_LRT_DAG3_DAG4_10kb_chr",
                   c, "_latents_", m, ".csv", sep=""))
    } else {
      print(paste("Model 4 did not converge for chr", c, "(10 kb)"))
    }
  }
}

if(bootstrap) { # bootstrapping windows
   nboot <- 5e3
   chr_list <- vector("list", length=num_chrom)
   for(c in 1:22) {
  
     cat(paste("chrom", c, "\n"))
  
     cXX <- filter(dt, chrom==c)
     cXX_r2_m1 <- numeric(length=nboot)
     cXX_r2_m2 <- numeric(length=nboot)
     cXX_r2_m3 <- numeric(length=nboot)
     cXX_r2_m4 <- numeric(length=nboot)
  
     pb <- txtProgressBar(min=0, max=nboot, style=3)
     for(b in 1:nboot) {
  
       setTxtProgressBar(pb, b)
       cXX_boot <- cXX %>% dplyr::slice_sample(prop=1, replace=TRUE)
  
       m1_fit <- sem(model_1,
                     data=cXX_boot,
                     estimator = "MLR",
                     missing = "fiml",
                     meanstructure = TRUE ,
                     std.ov = TRUE,
                     fixed.x = FALSE,
                     optim.method = "nlminb",
                     control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
       cXX_r2_m1[b] <- inspect(m1_fit, "r2")[1]
  
       m2_fit <- sem(model_2,
                     data=cXX_boot,
                     estimator = "MLR",
                     missing = "fiml",
                     meanstructure = TRUE,
                     std.ov = TRUE,
                     fixed.x = FALSE,
                     optim.method = "nlminb",
                     control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
       cXX_r2_m2[b] <- inspect(m2_fit, "r2")[1]
  
       m3_fit <- sem(model_3,
                     data=cXX_boot,
                     estimator = "MLR",
                     missing = "fiml",
                     meanstructure = TRUE ,
                     std.ov = TRUE,
                     fixed.x = FALSE,
                     optim.method = "nlminb",
                     control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
       cXX_r2_m3[b] <- inspect(m3_fit, "r2")[1]
  
       m4_fit <- sem(model_4,
                     data=cXX_boot,
                     estimator = "MLR",
                     missing = "fiml",
                     meanstructure = TRUE ,
                     std.ov = TRUE,
                     fixed.x = FALSE,
                     optim.method = "nlminb",
                     control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
       cXX_r2_m4[b] <- inspect(m4_fit, "r2")[1]
     }
     close(pb)
  
     cXX_boot_r2 <- cbind.data.frame(cXX_r2_m1, cXX_r2_m2, cXX_r2_m3, cXX_r2_m4)
     names(cXX_boot_r2) <- paste("model", 1:4, sep="_")
  
     cXX_boot_r2$boot <- 1:nrow(cXX_boot_r2)
     cXX_boot_r2$chrom <- c
     #cXX_boot_r2_m <- pivot_longer(cXX_boot_r2, cols=starts_with("model"), names_to="model", values_to="r-sqrd")
     chr_list[[c]] <- cXX_boot_r2
   }
  
   samp_size <- as.data.frame(table(dt$chrom))
   names(samp_size) <- c("chrom", "windows")
   samp_size$chrom <- as.integer(as.character(samp_size$chrom))
  
   boot_r2 <- data.table::rbindlist(chr_list)
   names(boot_r2)[1:4] <- paste("DAG", 1:4)
   boot_extended <- left_join(boot_r2, samp_size)
  
   cXX_boot_r2_m <- pivot_longer(boot_extended, cols=starts_with("DAG"), names_to="DAG", values_to="r-sqrd")
   cXX_boot_r2_m$chrom <- factor(cXX_boot_r2_m$chrom)
  
   b <- ggplot(cXX_boot_r2_m, aes(x = chrom, y = `r-sqrd`, fill=windows)) +
     geom_boxplot(alpha = 0.8) + facet_wrap(~DAG, ncol=2) +
     labs(x = "Chromosome", y = expression("Bootstrapped " * R[pi]^2)) + theme_classic() +
     scale_y_continuous(breaks=pretty_breaks(), limits=c(0, 1)) +
     theme(strip.text=element_text(size=16),
           axis.title=element_text(size=16),
           axis.text=element_text(size=12),
           legend.position="bottom")
   save_plot("sem/boot_r2_chroms_10kb.pdf", b, base_height=8, base_width=17)
}

if(permut) { # permuting chr labels
   npermut <- 1e3
   cXX_r2_permuts <- vector("list", length=npermut)
   pb <- txtProgressBar(min=0, max=npermut, style=3)
   for(p in 1:npermut) {
  
     setTxtProgressBar(pb, p)
  
     cXX <- dt
     cXX$chrom <- sample(cXX$chrom)
     chr_list <- vector("list", length=22)
  
     for(c in 1:22) {
       cXX_permut <- filter(cXX, chrom==c)
  
       m1_fit <- sem(model_1,
                     data=cXX_permut,
                     estimator = "MLR",
                     missing = "fiml",
                     meanstructure = TRUE ,
                     std.ov = TRUE,
                     fixed.x = FALSE,
                     optim.method = "nlminb",
                     control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
       m2_fit <- sem(model_2,
                     data=cXX_permut,
                     estimator = "MLR",
                     missing = "fiml",
                     meanstructure = TRUE ,
                     std.ov = TRUE,
                     fixed.x = FALSE,
                     optim.method = "nlminb",
                     control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
       m3_fit <- sem(model_3,
                     data=cXX_permut,
                     estimator = "MLR",
                     missing = "fiml",
                     meanstructure = TRUE ,
                     std.ov = TRUE,
                     fixed.x = FALSE,
                     optim.method = "nlminb",
                     control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
       m4_fit <- sem(model_4,
                     data=cXX_permut,
                     estimator = "MLR",
                     missing = "fiml",
                     meanstructure = TRUE ,
                     std.ov = TRUE,
                     fixed.x = FALSE,
                     optim.method = "nlminb",
                     control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
       chr_list[[c]] <- data.table(t(c(inspect(m1_fit, "r2")[1], inspect(m2_fit, "r2")[1], inspect(m3_fit, "r2")[1], inspect(m4_fit, "r2")[1], c)))
     }
  
     block <- data.table::rbindlist(chr_list)
     names(block) <- c(paste("DAG", 1:4, sep="_"), "chrom")
     block$permut <- p
     cXX_r2_permuts[[p]] <- block
   }
   close(pb)
  
   cXX_permut_r2 <- data.table::rbindlist(cXX_r2_permuts)
   cXX_permut_r2_m <- pivot_longer(cXX_permut_r2, cols=starts_with("DAG"), names_to="DAG", values_to="r-sqrd")
   cXX_permut_r2_m$chrom <- factor(cXX_permut_r2_m$chrom)
  
   p <- ggplot(cXX_permut_r2_m, aes(x = chrom, y = `r-sqrd`)) +
     geom_boxplot(alpha = 0.8) + facet_wrap(~DAG, ncol=2) +
     labs(x = "Chromosome", y = expression("Permuted " * R[pi]^2)) + theme_classic() +
     scale_y_continuous(breaks=pretty_breaks(), limits=c(0, 1)) +
     theme(strip.text=element_text(size=16),
           axis.title=element_text(size=16),
           axis.text=element_text(size=12),
           legend.position="bottom")
   save_plot("sem/permut_r2_chroms_10kb.pdf", p, base_height=8, base_width=17)
}

################# 100 kb #################
for(m in constrained_models) { 
  
  print(paste("100kb", m, Sys.Date(), Sys.time()), sep="\t")
  
  tmp_B <- filter(B_latent_100kb, constrained_model==m) %>% dplyr::select(., -constrained_model)
  
  tbl_100kb <- filter(human_maps_100kb, constrained_model==m, coverage >= 0.5)  
  table_100kb_g <- tbl_100kb[,1:3]
  table_100kb_g[, avg_pi := log(tbl_100kb$avg_pi)]
  table_100kb_g[, avg_rec := log(tbl_100kb$avg_rec)] # pyrho YRI pop.
  table_100kb_g[, prop_del_sites := log1p(tbl_100kb$prop_del_sites)]
  
  #table_100kb_g <- left_join(table_100kb_g, rec_latent_100kb, by=c("chrom", "chromStart", "chromEnd"))
  table_100kb_g <- left_join(table_100kb_g, tmp_B, by=c("chrom", "chromStart", "chromEnd"))
  table_100kb_g <- left_join(table_100kb_g, mut_latent_100kb, by=c("chrom", "chromStart", "chromEnd"))
  table_100kb_g <- add_genome_offsets_and_midpoint(table_100kb_g, chrom_lengths=NULL)
  
  setnames(table_100kb_g, old=c("B_latent_log", "mut_latent_log"), new=c("B", "avg_mut"))
  
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
  fwrite(combined_mat, paste("sem/pearson_gw_100kb_", m, ".csv", sep=""), row.names=T, quote=F)
  
  # model 1
  edges <- list(
    avg_pi = list(
      list(p = "avg_mut", label = "c1", start = 0.25),
      list(p = "B", label = "c2", start = 0.5)
    ),
    B = list(
      list(p = "avg_mut", label = "c3", start = -0.25)
    )
  )
  
  covs <- list()
  
  defined <- c(
    direct_mut_pi   = "c1",            
    indirect_mut_pi = "c3*c2",       
    total_mut_pi = "c1 + c3*c2", 
    direct_B_pi = "c2" 
  )
  
  model1_constrained <- make_lavaan_model(edges,
                                          covs = covs,
                                          defined = defined,
                                          groups = 22,
                                          free_by = "groups",
                                          label_prefix = "x")
  
  m1_gw <- sem(model1_constrained,
               data = dt,         
               estimator = "MLR",              
               missing = "fiml",             
               meanstructure = TRUE ,            
               std.ov = TRUE,
               fixed.x = FALSE,
               optim.method = "nlminb",
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m1_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m1_gw, "r2"))),
                            as.numeric(unlist(inspect(m1_gw, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_r2_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m1_gw, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_params_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m1_gw)),
                            as.numeric(fitMeasures(m1_gw))), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_fitmeasures_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
  } else {
    print("Model m1_gw did not converge (100 kb)")
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
  )
  
  m1_chr <- sem(model_free_with_group_defs,
                data = dt,         
                estimator = "MLR",              
                missing = "fiml",             
                meanstructure = TRUE ,            
                std.ov = TRUE,  
                fixed.x = FALSE,
                optim.method = "nlminb",         
                group="chrom",
                group.equal="none",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  lrt_res <- tryCatch(
    lavTestLRT(m1_gw, m1_chr, method = "satorra.bentler.2001"),
    error = function(e) {
      message("lavTestLRT failed: ", e$message)
      NULL
    }
  )
  
  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG1_gw-equal-freechrom_100kb_latents_", m, ".csv")
  
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
    aic_free <- AIC(m1_chr); aic_constr <- AIC(m1_gw)
    bic_free <- BIC(m1_chr); bic_constr <- BIC(m1_gw)
    df_free  <- tryCatch(fitMeasures(m1_chr, "df"), error = function(e) NA)
    df_constr<- tryCatch(fitMeasures(m1_gw, "df"), error = function(e) NA)
    nobs_free<- tryCatch(fitMeasures(m1_chr, "nobs"), error = function(e) NA)
    nobs_constr<- tryCatch(fitMeasures(m1_gw, "nobs"), error = function(e) NA)
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
  
  if(lavInspect(m1_chr, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m1_chr, "r2"))),
                            as.numeric(unlist(inspect(m1_chr, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_r2_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m1_chr, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_params_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m1_chr)),
                            as.numeric(fitMeasures(m1_chr))), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_fitmeasures_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
  } else {
    print("Model m1_chr did not converge (100 kb)")
  }
  
  # model 2
  edges <- list(
    avg_pi = list(
      list(p = "avg_mut", label = "c1", start = 0.25),
      list(p = "B", label = "c2", start = 0.5),
      list(p = "prop_del_sites", label = "c3", start = -0.25)
    ),
    B = list(
      list(p = "prop_del_sites", label = "c4", start = -0.25),
      list(p = "avg_mut", label = "c5", start = -0.25)
    ),
    avg_mut = list(), 
    prop_del_sites = list()
  )
  
  covs <- list()
  
  defined <- c(
    mu_B_pi = "c1*c5"
  )
  
  model2_constrained <- make_lavaan_model(edges,
                                          covs = covs,
                                          defined = defined,
                                          groups = 22,
                                          free_by = "groups",
                                          label_prefix = "x")
  
  m2_gw <- sem(model2_constrained,
               data = dt,         
               estimator = "MLR",              
               missing = "fiml",
               fixed.x = FALSE,
               meanstructure = TRUE ,            
               std.ov = TRUE,               
               optim.method = "nlminb",         
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m2_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m2_gw, "r2"))),
                            as.numeric(unlist(inspect(m2_gw, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_r2_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m2_gw, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_params_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m2_gw)),
                            as.numeric(fitMeasures(m2_gw))), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_fitmeasures_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
  } else {
    print("Model m2_gw did not converge (100 kb)")
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
  )
  
  m2_chr <- sem(model_free_with_group_defs,
                data = dt,         
                estimator = "MLR",              
                missing = "fiml",             
                meanstructure = TRUE ,            
                std.ov = TRUE,    
                fixed.x = FALSE,
                optim.method = "nlminb",         
                group="chrom",
                group.equal="none",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  lrt_res <- tryCatch(
    lavTestLRT(m2_gw, m2_chr, method = "satorra.bentler.2001"),
    error = function(e) {
      message("lavTestLRT failed: ", e$message)
      NULL
    }
  )
  
  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG2_gw-equal-freechrom_100kb_latents_", m, ".csv")
  
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
    aic_free <- AIC(m2_chr); aic_constr <- AIC(m2_gw)
    bic_free <- BIC(m2_chr); bic_constr <- BIC(m2_gw)
    df_free  <- tryCatch(fitMeasures(m2_chr, "df"), error = function(e) NA)
    df_constr<- tryCatch(fitMeasures(m2_gw, "df"), error = function(e) NA)
    nobs_free<- tryCatch(fitMeasures(m2_chr, "nobs"), error = function(e) NA)
    nobs_constr<- tryCatch(fitMeasures(m2_gw, "nobs"), error = function(e) NA)
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
  
  if(lavInspect(m2_chr, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m2_chr, "r2"))),
                            as.numeric(unlist(inspect(m2_chr, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_r2_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m2_chr, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_params_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m2_chr)),
                            as.numeric(fitMeasures(m2_chr))), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_fitmeasures_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
  } else {
    print("Model m2_chr did not converge (100 kb)")
  }
  
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
                                          groups = 22,
                                          free_by = "groups",
                                          label_prefix = "x")
  
  model3_constrained <- paste0(model3_constrained, "\n", "avg_mut ~~ 0*prop_del_sites")
  
  m3_gw <- sem(model3_constrained,
               data = dt,         
               estimator = "MLR",              
               missing = "fiml",             
               meanstructure = TRUE,            
               std.ov = TRUE,
               fixed.x = FALSE,
               optim.method = "nlminb",         
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m3_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m3_gw, "r2"))),
                            as.numeric(unlist(inspect(m3_gw, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG3_r2_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m3_gw, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG3_params_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m3_gw)),
                            as.numeric(fitMeasures(m3_gw))), 
           paste("sem/lavaan_out/obspi_YRI_DAG3_fitmeasures_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
  )
  
  model_free_with_group_defs <- paste0(model_free_with_group_defs, "\n", "avg_mut ~~ 0*prop_del_sites")
  
  m3_chr <- sem(model_free_with_group_defs,
                data = dt,         
                estimator = "MLR",              
                missing = "fiml",             
                meanstructure = TRUE ,            
                std.ov = TRUE,    
                fixed.x = FALSE,
                optim.method = "nlminb",         
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
  
  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG3_gw-equal-freechrom_100kb_latents_", m, ".csv")
  
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
           paste("sem/lavaan_out/obspi_YRI_DAG3_r2_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m3_chr, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG3_params_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m3_chr)),
                            as.numeric(fitMeasures(m3_chr))), 
           paste("sem/lavaan_out/obspi_YRI_DAG3_fitmeasures_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
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
                                          groups = 22,
                                          free_by = "groups",
                                          label_prefix = "x")
  
  m4_gw <- sem(model4_constrained,
               data = dt,         
               estimator = "MLR",              
               missing = "fiml",             
               meanstructure = TRUE ,            
               std.ov = TRUE,
               fixed.x = FALSE,
               optim.method = "nlminb",         
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m4_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m4_gw, "r2"))),
                            as.numeric(unlist(inspect(m4_gw, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG4_r2_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m4_gw, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG4_params_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m4_gw)),
                            as.numeric(fitMeasures(m4_gw))), 
           paste("sem/lavaan_out/obspi_YRI_DAG4_fitmeasures_100kb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
  )
  
  m4_chr <- sem(model_free_with_group_defs,
                data = dt,         
                estimator = "MLR",              
                missing = "fiml",             
                meanstructure = TRUE ,            
                std.ov = TRUE,    
                fixed.x = FALSE,
                optim.method = "nlminb",         
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
  
  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG4_gw-equal-freechrom_100kb_latents_", m, ".csv")
  
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
           paste("sem/lavaan_out/obspi_YRI_DAG4_r2_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m4_chr, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG4_params_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m4_chr)),
                            as.numeric(fitMeasures(m4_chr))), 
           paste("sem/lavaan_out/obspi_YRI_DAG4_fitmeasures_100kb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
  } else {
    print("Model m4_chr did not converge (100 kb)")
  }
  
  # one model per chrom
  for(c in 1:22) {
    
    data <- filter(table_100kb_g, chrom==c) %>%
      dplyr::select(., avg_pi, prop_del_sites, B, avg_rec, avg_mut)
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
    fwrite(combined_mat, paste("sem/pearson_chr", c, "_", m, "_100kb.csv", sep=""),
           row.names=T, quote=F)

    m1_fit <- sem(model_1,
                  data=filter(dt, chrom==c),
                  estimator  = "MLR",              
                  missing = "fiml",             
                  meanstructure = TRUE ,            
                  std.ov = TRUE,     
                  fixed.x = TRUE,
                  optim.method = "nlminb",    
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m1_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m1_fit, "r2")),
                              as.numeric(inspect(m1_fit, "r2"))), 
             paste("sem/lavaan_out/obspi_YRI_DAG1_r2_100kb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(parameterestimates(m1_fit, ci=T), 
             paste("sem/lavaan_out/obspi_YRI_DAG1_params_100kb_chr", 
                   c, "_latents_", m, ".csv", sep=""))
    } else {
      print(paste("Model 1 did not converge for chr", c, "(100 kb)"))
    }
    
    m2_fit <- sem(model_2,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",              
                  missing = "fiml",             
                  meanstructure = TRUE ,            
                  std.ov = TRUE,     
                  fixed.x = TRUE,
                  optim.method = "nlminb",    
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m2_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m2_fit, "r2")),
                              as.numeric(inspect(m2_fit, "r2"))), 
             paste("sem/lavaan_out/obspi_YRI_DAG2_r2_100kb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(parameterestimates(m2_fit, ci=T), 
             paste("sem/lavaan_out/obspi_YRI_DAG2_params_100kb_chr", 
                   c, "_latents_", m, ".csv", sep=""))
    } else {
      print(paste("Model 2 did not converge for chr", c, "(100 kb)"))
    }
    
    m3_fit <- sem(model_3,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",              
                  missing = "fiml",             
                  meanstructure = TRUE ,            
                  std.ov = TRUE,     
                  fixed.x = TRUE,
                  optim.method = "nlminb",    
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m3_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m3_fit, "r2")),
                              as.numeric(inspect(m3_fit, "r2"))), 
             paste("sem/lavaan_out/obspi_YRI_DAG3_r2_100kb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(parameterestimates(m3_fit, ci=T), 
             paste("sem/lavaan_out/obspi_YRI_DAG3_params_100kb_chr", 
                   c, "_latents_", m, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m3_fit)),
                              as.numeric(fitMeasures(m3_fit))), 
             paste("sem/lavaan_out/obspi_YRI_DAG3_fitmeasures_100kb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
    } else {
      print(paste("Model 3 did not converge for chr", c, "(100 kb)"))
    }
    
    m4_fit <- sem(model_4,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",              
                  missing = "fiml",             
                  meanstructure = TRUE ,            
                  std.ov = TRUE,     
                  fixed.x = TRUE,
                  optim.method = "nlminb",    
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m4_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m4_fit, "r2")),
                              as.numeric(inspect(m4_fit, "r2"))), 
             paste("sem/lavaan_out/obspi_YRI_DAG4_r2_100kb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(parameterestimates(m4_fit, ci=T), 
             paste("sem/lavaan_out/obspi_YRI_DAG4_params_100kb_chr", 
                   c, "_latents_", m, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m4_fit)),
                              as.numeric(fitMeasures(m4_fit))), 
             paste("sem/lavaan_out/obspi_YRI_DAG4_fitmeasures_100kb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(safeLRT(m3_fit, m4_fit),
             paste("sem/lavaan_out/obspi_LRT_DAG3_DAG4_100kb_chr",
                   c, "_latents_", m, ".csv", sep=""))
    } else {
      print(paste("Model 4 did not converge for chr", c, "(100 kb)"))
    }
  }
}

if(bootstrap) { # bootstrapping windows
  nboot <- 5e3
  chr_list <- vector("list", length=num_chrom)
  for(c in 1:22) {
    
    cat(paste("chrom", c, "\n"))
    
    cXX <- filter(dt, chrom==c)
    cXX_r2_m1 <- numeric(length=nboot)
    cXX_r2_m2 <- numeric(length=nboot)
    cXX_r2_m3 <- numeric(length=nboot)
    cXX_r2_m4 <- numeric(length=nboot)
    
    pb <- txtProgressBar(min=0, max=nboot, style=3)
    for(b in 1:nboot) {
      
      setTxtProgressBar(pb, b)
      cXX_boot <- cXX %>% dplyr::slice_sample(prop=1, replace=TRUE)
      
      m1_fit <- sem(model_1,
                    data=cXX_boot,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      cXX_r2_m1[b] <- inspect(m1_fit, "r2")[1]
      
      m2_fit <- sem(model_2,
                    data=cXX_boot,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      cXX_r2_m2[b] <- inspect(m2_fit, "r2")[1]
      
      m3_fit <- sem(model_3,
                    data=cXX_boot,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      cXX_r2_m3[b] <- inspect(m3_fit, "r2")[1]
      
      m4_fit <- sem(model_4,
                    data=cXX_boot,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      cXX_r2_m4[b] <- inspect(m4_fit, "r2")[1]
    }
    close(pb)
    
    cXX_boot_r2 <- cbind.data.frame(cXX_r2_m1, cXX_r2_m2, cXX_r2_m3, cXX_r2_m4)
    names(cXX_boot_r2) <- paste("model", 1:4, sep="_")
    
    cXX_boot_r2$boot <- 1:nrow(cXX_boot_r2)
    cXX_boot_r2$chrom <- c
    #cXX_boot_r2_m <- pivot_longer(cXX_boot_r2, cols=starts_with("model"), names_to="model", values_to="r-sqrd")
    chr_list[[c]] <- cXX_boot_r2
  }
  
  samp_size <- as.data.frame(table(dt$chrom))
  names(samp_size) <- c("chrom", "windows")
  samp_size$chrom <- as.integer(as.character(samp_size$chrom))
  
  boot_r2 <- data.table::rbindlist(chr_list)
  names(boot_r2)[1:4] <- paste("DAG", 1:4)
  boot_extended <- left_join(boot_r2, samp_size)
  
  cXX_boot_r2_m <- pivot_longer(boot_extended, cols=starts_with("DAG"), names_to="DAG", values_to="r-sqrd")
  cXX_boot_r2_m$chrom <- factor(cXX_boot_r2_m$chrom)
  
  b <- ggplot(cXX_boot_r2_m, aes(x = chrom, y = `r-sqrd`, fill=windows)) +
    geom_boxplot(alpha = 0.8) + facet_wrap(~DAG, ncol=2) +
    labs(x = "Chromosome", y = expression("Bootstrapped " * R[pi]^2)) + theme_classic() +
    scale_y_continuous(breaks=pretty_breaks(), limits=c(0, 1)) +
    theme(strip.text=element_text(size=16),
          axis.title=element_text(size=16),
          axis.text=element_text(size=12),
          legend.position="bottom")
  
  save_plot("sem/boot_r2_chroms_100kb.pdf", b, base_height=8, base_width=17)
}

if(permut) { # permuting chr labels
  npermut <- 1e3
  cXX_r2_permuts <- vector("list", length=npermut)
  pb <- txtProgressBar(min=0, max=npermut, style=3)
  for(p in 1:npermut) {
    
    setTxtProgressBar(pb, p)
    
    cXX <- dt
    cXX$chrom <- sample(cXX$chrom)
    
    chr_list <- vector("list", length=22)
    
    for(c in 1:22) {
      cXX_permut <- filter(cXX, chrom==c)
      
      m1_fit <- sem(model_1,
                    data=cXX_permut,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      m2_fit <- sem(model_2,
                    data=cXX_permut,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      m3_fit <- sem(model_3,
                    data=cXX_permut,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      m4_fit <- sem(model_4,
                    data=cXX_permut,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      chr_list[[c]] <- data.table(t(c(inspect(m1_fit, "r2")[1], inspect(m2_fit, "r2")[1], inspect(m3_fit, "r2")[1], inspect(m4_fit, "r2")[1], c)))
    }
    
    block <- data.table::rbindlist(chr_list)
    names(block) <- c(paste("DAG", 1:4, sep="_"), "chrom")
    block$permut <- p
    cXX_r2_permuts[[p]] <- block
  }
  close(pb)
  
  cXX_permut_r2 <- data.table::rbindlist(cXX_r2_permuts)
  cXX_permut_r2_m <- pivot_longer(cXX_permut_r2, cols=starts_with("DAG"), names_to="DAG", values_to="r-sqrd")
  cXX_permut_r2_m$chrom <- factor(cXX_permut_r2_m$chrom)
  
  p <- ggplot(cXX_permut_r2_m, aes(x = chrom, y = `r-sqrd`)) +
    geom_boxplot(alpha = 0.8) + facet_wrap(~DAG, ncol=2) +
    labs(x = "Chromosome", y = expression("Permuted " * R[pi]^2)) + theme_classic() +
    scale_y_continuous(breaks=pretty_breaks(), limits=c(0, 1)) +
    theme(strip.text=element_text(size=16),
          axis.title=element_text(size=16),
          axis.text=element_text(size=12),
          legend.position="bottom")
  save_plot("sem/permut_r2_chroms_100kb.pdf", p, base_height=8, base_width=17)
}

################# 1 Mb #################
for(m in constrained_models) { 
  
  print(paste("1Mb", m, Sys.Date(), Sys.time()), sep="\t")
  
  tmp_B <- filter(B_latent_1Mb, constrained_model==m) %>% dplyr::select(., -constrained_model)
  
  tbl_1Mb <- filter(human_maps_1Mb, constrained_model==m, coverage >= 0.2)  
  table_1Mb_g <- tbl_1Mb[,1:3]
  table_1Mb_g[, avg_pi := log(tbl_1Mb$avg_pi)]
  table_1Mb_g[, avg_rec := log(tbl_1Mb$avg_rec)] # pyrho YRI pop.
  table_1Mb_g[, prop_del_sites := log1p(tbl_1Mb$prop_del_sites)]
  
  #table_1Mb_g <- left_join(table_1Mb_g, rec_latent_1Mb, by=c("chrom", "chromStart", "chromEnd"))
  table_1Mb_g <- left_join(table_1Mb_g, tmp_B, by=c("chrom", "chromStart", "chromEnd"))
  table_1Mb_g <- left_join(table_1Mb_g, mut_latent_1Mb, by=c("chrom", "chromStart", "chromEnd"))
  table_1Mb_g <- add_genome_offsets_and_midpoint(table_1Mb_g, chrom_lengths=NULL)
  
  setnames(table_1Mb_g, old=c("B_latent_log", "mut_latent_log"), new=c("B", "avg_mut"))
  
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
  fwrite(combined_mat, paste("sem/pearson_gw_1Mb_", m, ".csv", sep=""), row.names=T, quote=F)
  
  # model 1
  edges <- list(
    avg_pi = list(
      list(p = "avg_mut", label = "c1", start = 0.25),
      list(p = "B", label = "c2", start = 0.5)
    ),
    B = list(
      list(p = "avg_mut", label = "c3", start = -0.25)
    )
  )
  
  covs <- list()
  
  defined <- c(
    direct_mut_pi   = "c1",            
    indirect_mut_pi = "c3*c2",       
    total_mut_pi = "c1 + c3*c2", 
    direct_B_pi = "c2" 
  )
  
  model1_constrained <- make_lavaan_model(edges,
                                          covs = covs,
                                          defined = defined,
                                          groups = 22,
                                          free_by = "groups",
                                          label_prefix = "x")
  
  m1_gw <- sem(model1_constrained,
               data = dt,         
               estimator = "MLR",              
               missing = "fiml",             
               meanstructure = TRUE ,            
               std.ov = TRUE,
               fixed.x = FALSE,
               optim.method = "nlminb",
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m1_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m1_gw, "r2"))),
                            as.numeric(unlist(inspect(m1_gw, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_r2_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m1_gw, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_params_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m1_gw)),
                            as.numeric(fitMeasures(m1_gw))), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_fitmeasures_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
  } else {
    print("Model m1_gw did not converge (100 Mb)")
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
  )
  
  m1_chr <- sem(model_free_with_group_defs,
                data = dt,         
                estimator = "MLR",              
                missing = "fiml",             
                meanstructure = TRUE ,            
                std.ov = TRUE,  
                fixed.x = FALSE,
                optim.method = "nlminb",         
                group="chrom",
                group.equal="none",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  lrt_res <- tryCatch(
    lavTestLRT(m1_gw, m1_chr, method = "satorra.bentler.2001"),
    error = function(e) {
      message("lavTestLRT failed: ", e$message)
      NULL
    }
  )
  
  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG1_gw-equal-freechrom_1Mb_latents_", m, ".csv")
  
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
    aic_free <- AIC(m1_chr); aic_constr <- AIC(m1_gw)
    bic_free <- BIC(m1_chr); bic_constr <- BIC(m1_gw)
    df_free <- tryCatch(fitMeasures(m1_chr, "df"), error = function(e) NA)
    df_constr<- tryCatch(fitMeasures(m1_gw, "df"), error = function(e) NA)
    nobs_free <- tryCatch(fitMeasures(m1_chr, "nobs"), error = function(e) NA)
    nobs_constr <- tryCatch(fitMeasures(m1_gw, "nobs"), error = function(e) NA)
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
  
  if(lavInspect(m1_chr, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m1_chr, "r2"))),
                            as.numeric(unlist(inspect(m1_chr, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_r2_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m1_chr, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_params_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m1_chr)),
                            as.numeric(fitMeasures(m1_chr))), 
           paste("sem/lavaan_out/obspi_YRI_DAG1_fitmeasures_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
  } else {
    print("Model m1_chr did not converge (1 Mb)")
  }
  
  # model 2
  edges <- list(
    avg_pi = list(
      list(p = "avg_mut", label = "c1", start = 0.25),
      list(p = "B", label = "c2", start = 0.5),
      list(p = "prop_del_sites", label = "c3", start = -0.25)
    ),
    B = list(
      list(p = "prop_del_sites", label = "c4", start = -0.25),
      list(p = "avg_mut", label = "c5", start = -0.25)
    ),
    avg_mut = list(), 
    prop_del_sites = list()
  )
  
  covs <- list()
  
  defined <- c(
    mu_B_pi = "c1*c5"
  )
  
  model2_constrained <- make_lavaan_model(edges,
                                          covs = covs,
                                          defined = defined,
                                          groups = 22,
                                          free_by = "groups",
                                          label_prefix = "x")
  
  m2_gw <- sem(model2_constrained,
               data = dt,         
               estimator = "MLR",              
               missing = "fiml",
               fixed.x = FALSE,
               meanstructure = TRUE ,            
               std.ov = TRUE,               
               optim.method = "nlminb",         
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m2_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m2_gw, "r2"))),
                            as.numeric(unlist(inspect(m2_gw, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_r2_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m2_gw, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_params_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m2_gw)),
                            as.numeric(fitMeasures(m2_gw))), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_fitmeasures_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
  } else {
    print("Model m2_gw did not converge (1 Mb)")
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
  )
  
  m2_chr <- sem(model_free_with_group_defs,
                data = dt,         
                estimator = "MLR",              
                missing = "fiml",             
                meanstructure = TRUE ,            
                std.ov = TRUE,    
                fixed.x = FALSE,
                optim.method = "nlminb",         
                group="chrom",
                group.equal="none",
                control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  lrt_res <- tryCatch(
    lavTestLRT(m2_gw, m2_chr, method = "satorra.bentler.2001"),
    error = function(e) {
      message("lavTestLRT failed: ", e$message)
      NULL
    }
  )
  
  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG2_gw-equal-freechrom_1Mb_latents_", m, ".csv")
  
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
    aic_free <- AIC(m2_chr); aic_constr <- AIC(m2_gw)
    bic_free <- BIC(m2_chr); bic_constr <- BIC(m2_gw)
    df_free  <- tryCatch(fitMeasures(m2_chr, "df"), error = function(e) NA)
    df_constr<- tryCatch(fitMeasures(m2_gw, "df"), error = function(e) NA)
    nobs_free<- tryCatch(fitMeasures(m2_chr, "nobs"), error = function(e) NA)
    nobs_constr<- tryCatch(fitMeasures(m2_gw, "nobs"), error = function(e) NA)
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
  
  if(lavInspect(m2_chr, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m2_chr, "r2"))),
                            as.numeric(unlist(inspect(m2_chr, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_r2_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m2_chr, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_params_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m2_chr)),
                            as.numeric(fitMeasures(m2_chr))), 
           paste("sem/lavaan_out/obspi_YRI_DAG2_fitmeasures_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
  } else {
    print("Model m2_chr did not converge (1 Mb)")
  }
  
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
                                          groups = 22,
                                          free_by = "groups",
                                          label_prefix = "x")
  
  model3_constrained <- paste0(model3_constrained, "\n", "avg_mut ~~ 0*prop_del_sites")
  
  m3_gw <- sem(model3_constrained,
               data = dt,         
               estimator = "MLR",              
               missing = "fiml",             
               meanstructure = TRUE ,            
               std.ov = TRUE,
               fixed.x = FALSE,
               optim.method = "nlminb",         
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m3_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m3_gw, "r2"))),
                            as.numeric(unlist(inspect(m3_gw, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG3_r2_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m3_gw, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG3_params_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m3_gw)),
                            as.numeric(fitMeasures(m3_gw))), 
           paste("sem/lavaan_out/obspi_YRI_DAG3_fitmeasures_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
  )
  
  model_free_with_group_defs <- paste0(model_free_with_group_defs, "\n", "avg_mut ~~ 0*prop_del_sites")
  
  m3_chr <- sem(model_free_with_group_defs,
                data = dt,         
                estimator = "MLR",              
                missing = "fiml",             
                meanstructure = TRUE ,            
                std.ov = TRUE,    
                fixed.x = FALSE,
                optim.method = "nlminb",         
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
  
  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG3_gw-equal-freechrom_1Mb_latents_", m, ".csv")
  
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
           paste("sem/lavaan_out/obspi_YRI_DAG3_r2_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m3_chr, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG3_params_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m3_chr)),
                            as.numeric(fitMeasures(m3_chr))), 
           paste("sem/lavaan_out/obspi_YRI_DAG3_fitmeasures_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
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
    rec_mu_B_pi  = "c2*c5*c7"
  )
  
  model4_constrained <- make_lavaan_model(edges,
                                          covs = covs,
                                          defined = defined,
                                          groups = 22,
                                          free_by = "groups",
                                          label_prefix = "x")
  
  m4_gw <- sem(model4_constrained,
               data = dt,         
               estimator = "MLR",              
               missing = "fiml",             
               meanstructure = TRUE ,            
               std.ov = TRUE,
               fixed.x = FALSE,
               optim.method = "nlminb",         
               group="chrom",
               group.equal="regressions",
               control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
  if(lavInspect(m4_gw, "converged")) {
    fwrite(rbind.data.frame(names(unlist(inspect(m4_gw, "r2"))),
                            as.numeric(unlist(inspect(m4_gw, "r2")))), 
           paste("sem/lavaan_out/obspi_YRI_DAG4_r2_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m4_gw, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG4_params_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m4_gw)),
                            as.numeric(fitMeasures(m4_gw))), 
           paste("sem/lavaan_out/obspi_YRI_DAG4_fitmeasures_1Mb_genomewide_", 
                 "latents_", m, ".csv", sep=""),
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
    groups = 22, per_group_map = per_group_map, label_prefix = "x"
  )
  
  m4_chr <- sem(model_free_with_group_defs,
                data = dt,         
                estimator = "MLR",              
                missing = "fiml",             
                meanstructure = TRUE ,            
                std.ov = TRUE,    
                fixed.x = FALSE,
                optim.method = "nlminb",         
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
  
  out_path <- paste0("sem/lavaan_out/obspi_LRT_DAG4_gw-equal-freechrom_1Mb_latents_", m, ".csv")
  
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
           paste("sem/lavaan_out/obspi_YRI_DAG4_r2_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
    
    fwrite(parameterestimates(m4_chr, ci=T), 
           paste("sem/lavaan_out/obspi_YRI_DAG4_params_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""))
    
    fwrite(rbind.data.frame(names(fitMeasures(m4_chr)),
                            as.numeric(fitMeasures(m4_chr))), 
           paste("sem/lavaan_out/obspi_YRI_DAG4_fitmeasures_1Mb_genomewide_", 
                 "latents_", m, "_freechrom.csv", sep=""),
           col.names=F)
  } else {
    print("Model m4_chr did not converge (1 Mb)")
  }
  
  # one model per chrom
  for(c in 1:22) {
    
    data <- filter(table_1Mb_g, chrom==c) %>%
      dplyr::select(., avg_pi, prop_del_sites, B, avg_rec, avg_mut)
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
    fwrite(combined_mat, paste("sem/pearson_chr", c, "_", m, "_1Mb.csv", sep=""),
           row.names=T, quote=F)
    
    m1_fit <- sem(model_1,
                  data=filter(dt, chrom==c),
                  estimator  = "MLR",              
                  missing = "fiml",             
                  meanstructure = TRUE ,            
                  std.ov = TRUE,     
                  fixed.x = TRUE,
                  optim.method = "nlminb",    
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m1_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m1_fit, "r2")),
                              as.numeric(inspect(m1_fit, "r2"))), 
             paste("sem/lavaan_out/obspi_YRI_DAG1_r2_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(parameterestimates(m1_fit, ci=T), 
             paste("sem/lavaan_out/obspi_YRI_DAG1_params_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""))
    } else {
      print(paste("Model 1 did not converge for chr", c, "(1 Mb)"))
    }
    
    m2_fit <- sem(model_2,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",              
                  missing = "fiml",             
                  meanstructure = TRUE ,            
                  std.ov = TRUE,     
                  fixed.x = TRUE,
                  optim.method = "nlminb",    
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m2_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m2_fit, "r2")),
                              as.numeric(inspect(m2_fit, "r2"))), 
             paste("sem/lavaan_out/obspi_YRI_DAG2_r2_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(parameterestimates(m2_fit, ci=T), 
             paste("sem/lavaan_out/obspi_YRI_DAG2_params_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""))
    } else {
      print(paste("Model 2 did not converge for chr", c, "(1 Mb)"))
    }
    
    m3_fit <- sem(model_3,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",              
                  missing = "fiml",             
                  meanstructure = TRUE ,            
                  std.ov = TRUE,     
                  fixed.x = TRUE,
                  optim.method = "nlminb",    
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m3_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m3_fit, "r2")),
                              as.numeric(inspect(m3_fit, "r2"))), 
             paste("sem/lavaan_out/obspi_YRI_DAG3_r2_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(parameterestimates(m3_fit, ci=T), 
             paste("sem/lavaan_out/obspi_YRI_DAG3_params_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m3_fit)),
                              as.numeric(fitMeasures(m3_fit))), 
             paste("sem/lavaan_out/obspi_YRI_DAG3_fitmeasures_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
    } else {
      print(paste("Model 3 did not converge for chr", c, "(1 Mb)"))
    }
    
    m4_fit <- sem(model_4,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",              
                  missing = "fiml",             
                  meanstructure = TRUE ,            
                  std.ov = TRUE,     
                  fixed.x = TRUE,
                  optim.method = "nlminb",    
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m4_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m4_fit, "r2")),
                              as.numeric(inspect(m4_fit, "r2"))), 
             paste("sem/lavaan_out/obspi_YRI_DAG4_r2_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(parameterestimates(m4_fit, ci=T), 
             paste("sem/lavaan_out/obspi_YRI_DAG4_params_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m4_fit)),
                              as.numeric(fitMeasures(m4_fit))), 
             paste("sem/lavaan_out/obspi_YRI_DAG4_fitmeasures_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(safeLRT(m3_fit, m4_fit),
             paste("sem/lavaan_out/obspi_LRT_DAG3_DAG4_1Mb_chr",
                   c, "_latents_", m, ".csv", sep=""))
    } else {
      print(paste("Model 4 did not converge for chr", c, "(1 Mb)"))
    }
  }
}

if(bootstrap) { # bootstrapping windows
  nboot <- 5e3
  chr_list <- vector("list", length=num_chrom)
  for(c in 1:22) {
    
    cat(paste("chrom", c, "\n"))
    
    cXX <- filter(dt, chrom==c)
    cXX_r2_m1 <- numeric(length=nboot)
    cXX_r2_m2 <- numeric(length=nboot)
    cXX_r2_m3 <- numeric(length=nboot)
    cXX_r2_m4 <- numeric(length=nboot)
    
    pb <- txtProgressBar(min=0, max=nboot, style=3)
    for(b in 1:nboot) {
      
      setTxtProgressBar(pb, b)
      cXX_boot <- cXX %>% dplyr::slice_sample(prop=1, replace=TRUE)
      
      m1_fit <- sem(model_1,
                    data=cXX_boot,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      cXX_r2_m1[b] <- inspect(m1_fit, "r2")[1]
      
      m2_fit <- sem(model_2,
                    data=cXX_boot,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      cXX_r2_m2[b] <- inspect(m2_fit, "r2")[1]
      
      m3_fit <- sem(model_3,
                    data=cXX_boot,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      cXX_r2_m3[b] <- inspect(m3_fit, "r2")[1]
      
      m4_fit <- sem(model_4,
                    data=cXX_boot,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      cXX_r2_m4[b] <- inspect(m4_fit, "r2")[1]
    }
    close(pb)
    
    cXX_boot_r2 <- cbind.data.frame(cXX_r2_m1, cXX_r2_m2, cXX_r2_m3, cXX_r2_m4)
    names(cXX_boot_r2) <- paste("model", 1:4, sep="_")
    
    cXX_boot_r2$boot <- 1:nrow(cXX_boot_r2)
    cXX_boot_r2$chrom <- c
    #cXX_boot_r2_m <- pivot_longer(cXX_boot_r2, cols=starts_with("model"), names_to="model", values_to="r-sqrd")
    chr_list[[c]] <- cXX_boot_r2
  }
  
  samp_size <- as.data.frame(table(dt$chrom))
  names(samp_size) <- c("chrom", "windows")
  samp_size$chrom <- as.integer(as.character(samp_size$chrom))
  
  boot_r2 <- data.table::rbindlist(chr_list)
  names(boot_r2)[1:4] <- paste("DAG", 1:4)
  boot_extended <- left_join(boot_r2, samp_size)
  
  cXX_boot_r2_m <- pivot_longer(boot_extended, cols=starts_with("DAG"), names_to="DAG", values_to="r-sqrd")
  cXX_boot_r2_m$chrom <- factor(cXX_boot_r2_m$chrom)
  
  b <- ggplot(cXX_boot_r2_m, aes(x = chrom, y = `r-sqrd`, fill=windows)) +
    geom_boxplot(alpha = 0.8) + facet_wrap(~DAG, ncol=2) +
    labs(x = "Chromosome", y = expression("Bootstrapped " * R[pi]^2)) + theme_classic() +
    scale_y_continuous(breaks=pretty_breaks(), limits=c(0, 1)) +
    theme(strip.text=element_text(size=16),
          axis.title=element_text(size=16),
          axis.text=element_text(size=12),
          legend.position="bottom")
  
  save_plot("sem/boot_r2_chroms_1Mb.pdf", b, base_height=8, base_width=17)
}

if(permut) { # permuting chr labels
  npermut <- 1e3
  cXX_r2_permuts <- vector("list", length=npermut)
  pb <- txtProgressBar(min=0, max=npermut, style=3)
  for(p in 1:npermut) {
    
    setTxtProgressBar(pb, p)
    
    cXX <- dt
    cXX$chrom <- sample(cXX$chrom)
    
    chr_list <- vector("list", length=22)
    
    for(c in 1:22) {
      cXX_permut <- filter(cXX, chrom==c)
      
      m1_fit <- sem(model_1,
                    data=cXX_permut,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      m2_fit <- sem(model_2,
                    data=cXX_permut,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      m3_fit <- sem(model_3,
                    data=cXX_permut,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      m4_fit <- sem(model_4,
                    data=cXX_permut,
                    estimator = "MLR",
                    missing = "fiml",
                    meanstructure = TRUE ,
                    std.ov = TRUE,
                    fixed.x = FALSE,
                    optim.method = "nlminb",
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      chr_list[[c]] <- data.table(t(c(inspect(m1_fit, "r2")[1], inspect(m2_fit, "r2")[1], inspect(m3_fit, "r2")[1], inspect(m4_fit, "r2")[1], c)))
    }
    
    block <- data.table::rbindlist(chr_list)
    names(block) <- c(paste("DAG", 1:4, sep="_"), "chrom")
    block$permut <- p
    cXX_r2_permuts[[p]] <- block
  }
  close(pb)
  
  cXX_permut_r2 <- data.table::rbindlist(cXX_r2_permuts)
  cXX_permut_r2_m <- pivot_longer(cXX_permut_r2, cols=starts_with("DAG"), names_to="DAG", values_to="r-sqrd")
  cXX_permut_r2_m$chrom <- factor(cXX_permut_r2_m$chrom)
  
  p <- ggplot(cXX_permut_r2_m, aes(x = chrom, y = `r-sqrd`)) +
    geom_boxplot(alpha = 0.8) + facet_wrap(~DAG, ncol=2) +
    labs(x = "Chromosome", y = expression("Permuted " * R[pi]^2)) + theme_classic() +
    scale_y_continuous(breaks=pretty_breaks(), limits=c(0, 1)) +
    theme(strip.text=element_text(size=16),
          axis.title=element_text(size=16),
          axis.text=element_text(size=12),
          legend.position="bottom")
  save_plot("sem/permut_r2_chroms_1Mb.pdf", p, base_height=8, base_width=17)
}
