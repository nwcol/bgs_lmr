
#################
# 
# This script uses prepared tables to build SEMs 
# NOTE: this script should be run from within each models dir
#
#################

source("~/Devel/bgs_lmr/tools/Rscripts/HelperFunctions.R")

constrained_models <- c("split_cds_phastcons", "split_cds_regulatory")
mutation_maps <- c("roulette")

#################
#
# Fitting models with E[pi] + noise
#
#################

##################  10 kb ################## 
human_maps_10kb <- fread("maps/maps_10kb.csv.gz")

for(u in mutation_maps) {
  for(m in constrained_models) { 
    
    print(paste(u, m, Sys.Date(), Sys.time()), sep="\t")

    tbl_10kb <- filter(human_maps_10kb,
                       coverage >= 0.75,
                       mut_map==u,
                       constrained_model==m)
    
    tbl_10kb_log <- tbl_10kb[,1:3]
    tbl_10kb_log[, B := log(tbl_10kb$B)]
    tbl_10kb_log[, avg_mut := log(tbl_10kb$avg_mut)]
    tbl_10kb_log[, avg_rec := log(tbl_10kb$avg_rec)]
    tbl_10kb_log[, pi_sim_tot := log(tbl_10kb$pi_sim_tot)]
    tbl_10kb_log[, prop_del_sites := log1p(tbl_10kb$prop_del_sites)]
    
    table_10kb_g <- add_genome_offsets_and_midpoint(tbl_10kb_log, chrom_lengths=NULL)
    dt <- copy(table_10kb_g) %>% arrange(chrom, midpoint) %>% as.data.table() %>% setDT()
    
    # model 1
    edges <- list(
      pi_sim_tot = list(
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
                 std.ov = TRUE,
                 fixed.x = FALSE,
                 optim.method = "nlminb",         
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m1_gw, "converged")) {
      fwrite(parameterestimates(m1_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_params_10kb_genomewide_", 
                   m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m1_gw)),
                              as.numeric(fitMeasures(m1_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_fitmeasures_10kb_genomewide_", 
                   m, "_", u, ".csv", sep=""),
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
                  data  = sem_df_std,         
                  estimator = "MLR",              
                  missing = "fiml",             
                  std.ov = TRUE,  
                  fixed.x = FALSE,
                  optim.method = "nlminb",         
                  group="chrom",
                  group.equal="none",
                  control=list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    lrt_res <- tryCatch(
      lavTestLRT(m1_gw, m1_chr, method = "satorra.bentler.2001"),
      error = function(e) {
        message("lavTestLRT failed: ", e$message)
        NULL
      }
    )
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG1_gw-equal-freechrom_10kb_", m, "_", u, ".csv")
    
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
      fwrite(parameterestimates(m1_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_params_10kb_genomewide_", 
                  m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m1_chr)),
                              as.numeric(fitMeasures(m1_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_fitmeasures_10kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    } else {
      print("Model m1_chr did not converge (10 kb)")
    }
    
    # model 2
    edges <- list(
      pi_sim_tot = list(
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
                 std.ov = TRUE,               
                 optim.method = "nlminb",         
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m2_gw, "converged")) {
      fwrite(parameterestimates(m2_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_params_10kb_genomewide_", 
                   m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m2_gw)),
                              as.numeric(fitMeasures(m2_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_fitmeasures_10kb_genomewide_", 
                   m, "_", u, ".csv", sep=""),
             col.names=F)
    } else {
      print("Model m2_gw did not converge (10 kb)")
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
    
    m2_chr <- sem(model_free_with_group_defs,
                  data = sem_df_std,         
                  estimator = "MLR",              
                  missing = "fiml",             
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
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG2_gw-equal-freechrom_10kb_", m, "_", u, ".csv")
    
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
      # Fallback: AIC/BIC 
      aic_free <- AIC(m2_chr); aic_constr <- AIC(m2_gw)
      bic_free <- BIC(m2_chr); bic_constr <- BIC(m2_gw)
      df_free <- tryCatch(fitMeasures(m2_chr, "df"), error = function(e) NA)
      df_constr<- tryCatch(fitMeasures(m2_gw, "df"), error = function(e) NA)
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
      fwrite(parameterestimates(m2_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_params_10kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m2_chr)),
                              as.numeric(fitMeasures(m2_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_fitmeasures_10kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    } else {
      print("Model m2_chr did not converge (10 kb)")
    }
    
    # model 3
    edges <- list(
      pi_sim_tot = list(
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
                 estimator  = "MLR",              
                 missing  = "fiml",             
                 std.ov = TRUE,
                 fixed.x = FALSE,
                 optim.method = "nlminb",         
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))

    if(lavInspect(m3_gw, "converged")) {
      fwrite(parameterestimates(m3_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_params_10kb_genomewide_", 
                   m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m3_gw)),
                              as.numeric(fitMeasures(m3_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_fitmeasures_10kb_genomewide_", 
                   m, "_", u, ".csv", sep=""),
             col.names=F)
    } else {
      print("Model m3_gw did not converge (10 kb)")
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
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG3_gw-equal-freechrom_10kb_", m, "_", u, ".csv")
    
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
      # Fallback: AIC/BIC 
      aic_free <- AIC(m3_chr); aic_constr <- AIC(m3_gw)
      bic_free <- BIC(m3_chr); bic_constr <- BIC(m3_gw)
      df_free <- tryCatch(fitMeasures(m3_chr, "df"), error = function(e) NA)
      df_constr <- tryCatch(fitMeasures(m3_gw, "df"), error = function(e) NA)
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
      fwrite(parameterestimates(m3_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_params_10kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m3_chr)),
                              as.numeric(fitMeasures(m3_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_fitmeasures_10kb_genomewide_", 
                  m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    } else {
      print("Model m3_chr did not converge (10 kb)")
    }
    
    # model 4
    edges <- list(
      pi_sim_tot = list(
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
                 estimator  = "MLR",              
                 missing = "fiml",             
                 std.ov = TRUE,
                 fixed.x = FALSE,
                 optim.method = "nlminb",         
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m4_gw, "converged")) {
      fwrite(parameterestimates(m4_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_params_10kb_genomewide_", 
                   m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m4_gw)),
                              as.numeric(fitMeasures(m4_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_fitmeasures_10kb_genomewide_", 
                  m, "_", u, ".csv", sep=""),
             col.names=F)
    } else {
      print("Model m4_gw did not converge (10 kb)")
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
    
    m4_chr <- sem(model_free_with_group_defs,
                  data = sem_df_std,         
                  estimator = "MLR",              
                  missing = "fiml",             
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
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG4_gw-equal-freechrom_10kb_", m, "_", u, ".csv")
    
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
      # Fallback: AIC/BIC 
      aic_free <- AIC(m4_chr); aic_constr <- AIC(m4_gw)
      bic_free <- BIC(m4_chr); bic_constr <- BIC(m4_gw)
      df_free <- tryCatch(fitMeasures(m4_chr, "df"), error = function(e) NA)
      df_constr <- tryCatch(fitMeasures(m4_gw, "df"), error = function(e) NA)
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
      fwrite(parameterestimates(m4_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_params_10kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m4_chr)),
                              as.numeric(fitMeasures(m4_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_fitmeasures_10kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    } else {
      print("Model m4_chr did not converge (10 kb)")
    }
    
    # one model per chrom
    for(c in 1:22) {

      model_1 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      B ~ start(-0.25)*c3*avg_mut
      "
      
      m1_fit <- sem(model_1,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      if(lavInspect(m1_fit, "converged")) {
        fwrite(parameterestimates(m1_fit, ci=T), 
               paste("sem/lavaan_out/simpi_YRI_DAG1_params_10kb_chr", 
                     c, "_", m, "_", u, ".csv", sep=""))
      } else {
        print(paste("Model 1 did not converge for chr", c, "(10 kb)"))
      }
      
      model_2 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      pi_sim_tot ~ start(-0.25)*c3*prop_del_sites
      B ~ start(-0.25)*c4*prop_del_sites
      B ~ start(-0.25)*c5*avg_mut
      mu_B_pi:=c1*c5
      "
      
      m2_fit <- sem(model_2,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      if(lavInspect(m2_fit, "converged")) {
        fwrite(parameterestimates(m2_fit, ci=T), 
               paste("sem/lavaan_out/simpi_YRI_DAG2_params_10kb_chr", 
                     c, "_", m, "_", u, ".csv", sep=""))
      } else {
        print(paste("Model 2 did not converge for chr", c, "(10 kb)"))
      }
      
      model_3 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      pi_sim_tot ~ start(-0.25)*c3*prop_del_sites
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
      
      m3_fit <- sem(model_3,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      if(lavInspect(m3_fit, "converged")) {
        fwrite(parameterestimates(m3_fit, ci=T), 
               paste("sem/lavaan_out/simpi_YRI_DAG3_params_10kb_chr", 
                     c, "_", m, "_", u, ".csv", sep=""))
        
        fwrite(rbind.data.frame(names(fitMeasures(m3_fit)),
                                as.numeric(fitMeasures(m3_fit))), 
               paste("sem/lavaan_out/simpi_YRI_DAG3_fitmeasures_10kb_chr", 
                     c, "_", m, "_", u, ".csv", sep=""),
               col.names=F)
      } else {
        print(paste("Model 3 did not converge for chr", c, "(10 kb)"))
      }
      
      model_4 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      pi_sim_tot ~ start(-0.25)*c3*prop_del_sites
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
      
      m4_fit <- sem(model_4,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      if(lavInspect(m4_fit, "converged")) {
        fwrite(parameterestimates(m4_fit, ci=T), 
               paste("sem/lavaan_out/simpi_YRI_DAG4_params_10kb_chr", 
                     c, "_", m, "_", u, ".csv", sep=""))
        
        fwrite(rbind.data.frame(names(fitMeasures(m4_fit)),
                                as.numeric(fitMeasures(m4_fit))), 
               paste("sem/lavaan_out/simpi_YRI_DAG4_fitmeasures_10kb_chr", 
                     c, "_", m, "_", u, ".csv", sep=""),
               col.names=F)
        
        fwrite(safeLRT(m3_fit, m4_fit),
               paste("sem/lavaan_out/simpi_LRT_DAG3_DAG4_10kb_chr",
                     c, "_", m, "_", u, ".csv", sep=""))
      } else {
        print(paste("Model 4 did not converge for chr", c, "(10 kb)"))
      }
    }
  }
}

##################  100 kb ################## 
human_maps_100kb <- fread("maps/maps_100kb.csv.gz") 

for(u in mutation_maps) {
  for(m in constrained_models) { 
    
    print(paste(u, m, Sys.Date(), Sys.time()), sep="\t")
    
    tbl_100kb <- filter(human_maps_100kb,
                       coverage >= 0.5,
                       mut_map==u,
                       constrained_model==m)
    
    tbl_100kb_log <- tbl_100kb[,1:3]
    tbl_100kb_log[, B := log(tbl_100kb$B)]
    tbl_100kb_log[, avg_mut := log(tbl_100kb$avg_mut)]
    tbl_100kb_log[, avg_rec := log(tbl_100kb$avg_rec)]
    tbl_100kb_log[, pi_sim_tot := log(tbl_100kb$pi_sim_tot)]
    tbl_100kb_log[, prop_del_sites := log1p(tbl_100kb$prop_del_sites)]
    
    table_100kb_g <- add_genome_offsets_and_midpoint(tbl_100kb_log, chrom_lengths=NULL)
    dt <- copy(table_100kb_g) %>% arrange(chrom, midpoint) %>% as.data.table() %>% setDT()
    
    # model 1
    edges <- list(
      pi_sim_tot = list(
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
    
    m1_gw <- sem(model1_constrained,
                 data  = dt,         
                 estimator = "MLR",              
                 missing = "fiml",             
                 std.ov = TRUE,
                 fixed.x = FALSE,
                 optim.method = "nlminb", 
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m1_gw, "converged")) {
      fwrite(parameterestimates(m1_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_params_100kb_genomewide_", 
                   m, "_", u,".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m1_gw)),
                              as.numeric(fitMeasures(m1_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_fitmeasures_100kb_genomewide_", 
                   m, "_", u, ".csv", sep=""),
             col.names=F)
    } else {
      print("Model m1_gw did not converge (100 kb)")
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
    
    m1_chr <- sem(model_free_with_group_defs,
                  data = dt,         
                  estimator = "MLR",              
                  missing = "fiml",             
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
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG1_gw-equal-freechrom_100kb_", m, "_", u, ".csv")
    
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
      fwrite(parameterestimates(m1_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_params_100kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m1_chr)),
                              as.numeric(fitMeasures(m1_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_fitmeasures_100kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    } else {
      print("Model m1_chr did not converge (100 kb)")
    }
    
    # model 2
    edges <- list(
      pi_sim_tot = list(
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
                 data  = dt,         
                 estimator = "MLR",              
                 missing = "fiml",
                 fixed.x = FALSE,
                 std.ov = TRUE,               
                 optim.method = "nlminb",         
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m2_gw, "converged")) {
      fwrite(parameterestimates(m2_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_params_100kb_genomewide_", 
                   m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m2_gw)),
                              as.numeric(fitMeasures(m2_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_fitmeasures_100kb_genomewide_", 
                   m, "_", u, ".csv", sep=""),
             col.names=F)
    } else {
      print("Model m2_gw did not converge (100 kb)")
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
    
    m2_chr <- sem(model_free_with_group_defs,
                  data = dt,         
                  estimator = "MLR",              
                  missing = "fiml",             
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
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG2_gw-equal-freechrom_100kb_", m, "_", u, ".csv")
    
    if(!is.null(lrt_res)) {
      df <- as.data.table(lrt_res, keep.rownames = "comparison")
      num_cols <- setdiff(names(df), "comparison")
      df[, (num_cols) := lapply(.SD, function(x) {
        xn <- as.numeric(x)
        if(all(is.na(xn)) && is.character(x)) return(x) else return(xn)
      }), .SDcols = num_cols]
      fwrite(df, out_path)
      message("Wrote lavTestLRT results to: ", out_path)
    } else {
      # Fallback: AIC/BIC 
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
      fwrite(parameterestimates(m2_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_params_100kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m2_chr)),
                              as.numeric(fitMeasures(m2_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_fitmeasures_100kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    }
    
    # model 3
    edges <- list(
      pi_sim_tot = list(
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
                 std.ov = TRUE,
                 fixed.x = FALSE,
                 optim.method = "nlminb",         
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
  
    if(lavInspect(m3_gw, "converged")) {
      fwrite(parameterestimates(m3_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_params_100kb_genomewide_", 
                   m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m3_gw)),
                              as.numeric(fitMeasures(m3_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_fitmeasures_100kb_genomewide_", 
                   m, "_", u, ".csv", sep=""),
             col.names=F)
    } else {
      print("Model m3_gw did not converge (100 kb)")
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
    
    model_free_with_group_defs <- paste0(model_free_with_group_defs, "\n", "avg_mut ~~ 0*prop_del_sites")
    
    m3_chr <- sem(model_free_with_group_defs,
                  data = dt,         
                  estimator = "MLR",              
                  missing = "fiml",             
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
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG3_gw-equal-freechrom_100kb_", m, "_", u, ".csv")
    
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
      df_constr <- tryCatch(fitMeasures(m3_gw, "df"), error = function(e) NA)
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
      fwrite(parameterestimates(m3_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_params_100kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m3_chr)),
                              as.numeric(fitMeasures(m3_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_fitmeasures_100kb_genomewide_", 
                    m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    } else {
      print("Model m3_chr did not converge (100 kb)")
    }
    
    # model 4
    edges <- list(
      pi_sim_tot = list(
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
                 std.ov = TRUE,
                 fixed.x = FALSE,
                 optim.method = "nlminb",         
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m4_gw, "converged")) {
      fwrite(parameterestimates(m4_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_params_100kb_genomewide_", 
                    m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m4_gw)),
                              as.numeric(fitMeasures(m4_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_fitmeasures_100kb_genomewide_", 
                   m, "_", u, ".csv", sep=""),
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
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG4_gw-equal-freechrom_100kb_", m, "_", u, ".csv")
    
    if(!is.null(lrt_res)) {
      df <- as.data.table(lrt_res, keep.rownames = "comparison")
      num_cols <- setdiff(names(df), "comparison")
      df[, (num_cols) := lapply(.SD, function(x) {
        xn <- as.numeric(x)
        if(all(is.na(xn)) && is.character(x)) return(x) else return(xn)
      }), .SDcols = num_cols]
      fwrite(df, out_path)
      message("Wrote lavTestLRT results to: ", out_path)
    } else {
      # Fallback: AIC/BIC
      aic_free <- AIC(m4_chr); aic_constr <- AIC(m4_gw)
      bic_free <- BIC(m4_chr); bic_constr <- BIC(m4_gw)
      df_free <- tryCatch(fitMeasures(m4_chr, "df"), error = function(e) NA)
      df_constr <- tryCatch(fitMeasures(m4_gw, "df"), error = function(e) NA)
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
      fwrite(parameterestimates(m4_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_params_100kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m4_chr)),
                              as.numeric(fitMeasures(m4_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_fitmeasures_100kb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    } else {
      print("Model m4_chr did not converge (100 kb)")
    }
    
    # one model per chrom
    for(c in 1:22) {
      
      model_1 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      B ~ start(-0.25)*c3*avg_mut
      "
      
      m1_fit <- sem(model_1,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      fwrite(parameterestimates(m1_fit, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_params_100kb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""))
      
      model_2 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      pi_sim_tot ~ start(-0.25)*c3*prop_del_sites
      B ~ start(-0.25)*c4*prop_del_sites
      B ~ start(-0.25)*c5*avg_mut
      mu_B_pi:=c1*c5
      "
      
      m2_fit <- sem(model_2,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      fwrite(parameterestimates(m2_fit, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_params_100kb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""))
      
      model_3 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      pi_sim_tot ~ start(-0.25)*c3*prop_del_sites
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
      
      m3_fit <- sem(model_3,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      fwrite(parameterestimates(m3_fit, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_params_100kb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m3_fit)),
                              as.numeric(fitMeasures(m3_fit))), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_fitmeasures_100kb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""),
             col.names=F)
      
      model_4 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      pi_sim_tot ~ start(-0.25)*c3*prop_del_sites
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
      
      m4_fit <- sem(model_4,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      fwrite(parameterestimates(m4_fit, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_params_100kb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m4_fit)),
                              as.numeric(fitMeasures(m4_fit))), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_fitmeasures_100kb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""),
             col.names=F)
      
      fwrite(safeLRT(m3_fit, m4_fit),
             paste("sem/lavaan_out/simpi_LRT_DAG3_DAG4_100kb_chr",
                   c, "_", m, "_", u, ".csv", sep=""))
    }
  }
}

##################  1 Mb ################## 
human_maps_1Mb <- fread("maps/maps_1Mb.csv.gz")

for(u in mutation_maps) {
  for(m in constrained_models) { 
    
    print(paste(u, m, Sys.Date(), Sys.time()), sep="\t")
    
    tbl_1Mb <- filter(human_maps_1Mb,
                        coverage >= 0.25,
                        mut_map==u,
                        constrained_model==m)
    
    tbl_1Mb_log <- tbl_1Mb[,1:3]
    tbl_1Mb_log[, B := log(tbl_1Mb$B)]
    tbl_1Mb_log[, avg_mut := log(tbl_1Mb$avg_mut)]
    tbl_1Mb_log[, avg_rec := log(tbl_1Mb$avg_rec)]
    tbl_1Mb_log[, pi_sim_tot := log(tbl_1Mb$pi_sim_tot)]
    tbl_1Mb_log[, prop_del_sites := log1p(tbl_1Mb$prop_del_sites)]
    
    table_1Mb_g <- add_genome_offsets_and_midpoint(tbl_1Mb_log, chrom_lengths=NULL)
    dt <- copy(table_1Mb_g) %>% arrange(chrom, midpoint) %>% as.data.table() %>% setDT()
    
    # model 1
    edges <- list(
      pi_sim_tot = list(
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
    
    m1_gw <- sem(model1_constrained,
                 data  = dt,         
                 estimator = "MLR",              
                 missing = "fiml",             
                 std.ov = TRUE,
                 fixed.x = FALSE,
                 optim.method = "nlminb", 
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m1_gw, "converged")) {
      fwrite(parameterestimates(m1_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_params_1Mb_genomewide_", 
                   m, "_", u,".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m1_gw)),
                              as.numeric(fitMeasures(m1_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_fitmeasures_1Mb_genomewide_", 
                   m, "_", u, ".csv", sep=""),
             col.names=F)
    } else {
      print("Model m1_gw did not converge (1 Mb)")
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
    
    m1_chr <- sem(model_free_with_group_defs,
                  data = dt,         
                  estimator = "MLR",              
                  missing = "fiml",             
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
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG1_gw-equal-freechrom_1Mb_", m, "_", u, ".csv")
    
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
      fwrite(parameterestimates(m1_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_params_1Mb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m1_chr)),
                              as.numeric(fitMeasures(m1_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_fitmeasures_1Mb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    } else {
      print("Model m1_chr did not converge (100 kb)")
    }
    
    # model 2
    edges <- list(
      pi_sim_tot = list(
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
                 data  = dt,         
                 estimator = "MLR",              
                 missing = "fiml",
                 fixed.x = FALSE,
                 std.ov = TRUE,               
                 optim.method = "nlminb",         
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m2_gw, "converged")) {
      fwrite(parameterestimates(m2_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_params_1Mb_genomewide_", 
                   m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m2_gw)),
                              as.numeric(fitMeasures(m2_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_fitmeasures_1Mb_genomewide_", 
                   m, "_", u, ".csv", sep=""),
             col.names=F)
    } else {
      print("Model m2_gw did not converge (100 kb)")
    }
    
    # build a per_group_map that marks every lhs~rhs as group-specific
    per_group_map <- character()
    for(lhs in names(edges)) {
      for(r in edges[[lhs]]) {
        rhs <- if(is.list(r)) r$p else r
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
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG2_gw-equal-freechrom_1Mb_", m, "_", u, ".csv")
    
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
      # Fallback: AIC/BIC 
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
      fwrite(parameterestimates(m2_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_params_1Mb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m2_chr)),
                              as.numeric(fitMeasures(m2_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_fitmeasures_1Mb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    }
    
    # model 3
    edges <- list(
      pi_sim_tot = list(
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
                 std.ov = TRUE,
                 fixed.x = FALSE,
                 optim.method = "nlminb",         
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m3_gw, "converged")) {
      fwrite(parameterestimates(m3_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_params_1Mb_genomewide_", 
                   m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m3_gw)),
                              as.numeric(fitMeasures(m3_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_fitmeasures_1Mb_genomewide_", 
                   m, "_", u, ".csv", sep=""),
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
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG3_gw-equal-freechrom_1Mb_", m, "_", u, ".csv")
    
    if(!is.null(lrt_res)) {
      df <- as.data.table(lrt_res, keep.rownames = "comparison")
      num_cols <- setdiff(names(df), "comparison")
      df[, (num_cols) := lapply(.SD, function(x) {
        xn <- as.numeric(x)
        if(all(is.na(xn)) && is.character(x)) return(x) else return(xn)
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
      fwrite(parameterestimates(m3_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_params_1Mb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m3_chr)),
                              as.numeric(fitMeasures(m3_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_fitmeasures_1Mb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    } else {
      print("Model m3_chr did not converge (100 kb)")
    }
    
    # model 4
    edges <- list(
      pi_sim_tot = list(
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
                 std.ov = TRUE,
                 fixed.x = FALSE,
                 optim.method = "nlminb",         
                 group="chrom",
                 group.equal="regressions",
                 control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m4_gw, "converged")) {
      fwrite(parameterestimates(m4_gw, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_params_1Mb_genomewide_", 
                   m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m4_gw)),
                              as.numeric(fitMeasures(m4_gw))), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_fitmeasures_1Mb_genomewide_", 
                   m, "_", u, ".csv", sep=""),
             col.names=F)
    } else {
      print("Model m4_gw did not converge (100 kb)")
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
    
    m4_chr <- sem(model_free_with_group_defs,
                  data = dt,         
                  estimator = "MLR",              
                  missing = "fiml",             
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
    
    out_path <- paste0("sem/lavaan_out/simpi_LRT_DAG4_gw-equal-freechrom_1Mb_", m, "_", u, ".csv")
    
    if(!is.null(lrt_res)) {
      df <- as.data.table(lrt_res, keep.rownames = "comparison")
      num_cols <- setdiff(names(df), "comparison")
      df[, (num_cols) := lapply(.SD, function(x) {
        xn <- as.numeric(x)
        if(all(is.na(xn)) && is.character(x)) return(x) else return(xn)
      }), .SDcols = num_cols]
      fwrite(df, out_path)
      message("Wrote lavTestLRT results to: ", out_path)
    } else {
      # Fallback: AIC/BIC 
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
      fwrite(parameterestimates(m4_chr, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_params_1Mb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m4_chr)),
                              as.numeric(fitMeasures(m4_chr))), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_fitmeasures_1Mb_genomewide_", 
                   m, "_", u, "_freechrom.csv", sep=""),
             col.names=F)
    } else {
      print("Model m4_chr did not converge (100 kb)")
    }
    
    # one model per chrom
    for(c in 1:22) {
      
      model_1 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      B ~ start(-0.25)*c3*avg_mut
      "
      
      m1_fit <- sem(model_1,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      fwrite(parameterestimates(m1_fit, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG1_params_1Mb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""))
      
      model_2 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      pi_sim_tot ~ start(-0.25)*c3*prop_del_sites
      B ~ start(-0.25)*c4*prop_del_sites
      B ~ start(-0.25)*c5*avg_mut
      mu_B_pi:=c1*c5
      "
      
      m2_fit <- sem(model_2,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      fwrite(parameterestimates(m2_fit, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG2_params_1Mb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""))
      
      model_3 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      pi_sim_tot ~ start(-0.25)*c3*prop_del_sites
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
      
      m3_fit <- sem(model_3,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      fwrite(parameterestimates(m3_fit, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_params_1Mb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m3_fit)),
                              as.numeric(fitMeasures(m3_fit))), 
             paste("sem/lavaan_out/simpi_YRI_DAG3_fitmeasures_1Mb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""),
             col.names=F)
      
      model_4 <- "
      pi_sim_tot ~ start(0.25)*c1*avg_mut
      pi_sim_tot ~ start(0.5)*c2*B
      pi_sim_tot ~ start(-0.25)*c3*prop_del_sites
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
      
      m4_fit <- sem(model_4,
                    data=filter(dt, chrom==c),
                    estimator = "MLR",              
                    missing = "fiml",             
                    std.ov = TRUE,     
                    fixed.x = TRUE,
                    optim.method = "nlminb",    
                    control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
      
      fwrite(parameterestimates(m4_fit, ci=T), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_params_1Mb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m4_fit)),
                              as.numeric(fitMeasures(m4_fit))), 
             paste("sem/lavaan_out/simpi_YRI_DAG4_fitmeasures_1Mb_chr", 
                   c, "_", m, "_", u, ".csv", sep=""),
             col.names=F)
      
      fwrite(safeLRT(m3_fit, m4_fit),
             paste("sem/lavaan_out/simpi_LRT_DAG3_DAG4_1Mb_chr",
                   c, "_", m, "_", u, ".csv", sep=""))
    }
  }
}

