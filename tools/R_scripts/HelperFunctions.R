#################
#
# Loading packages
#
#################

library(data.table)
library(tidyverse)
library(dplyr)
library(purrr)
library(tibble)
library(cowplot)
library(scales)
library(GenomicRanges)
library(R.utils)
library(RColorBrewer)
library(rtracklayer) # genonomic lift-over
library(patchwork)

library(stats)
library(MASS)
library(lmtest)
library(nlme)
library(car)
library(interactions)
library(ppcor)
library(Hmisc)
library(MVN)
library(Matrix)
options(rgl.useNULL=TRUE) # "solves" issue with rgl package in MacOS BigSur
library(heplots)
library(parallel) # mcapply
library(pbmcapply)

library(survcomp) # for fisher's C test, combine.test()
library(piecewiseSEM)
library(semTools)
library(graph)
library(ggm)
library(lavaan)
library(dagitty)
library(ggdag)
library(magrittr)
library(ggpubr)
library(GGally)

#################
#
# Helper for plotting
#
#################

scale.4d <- function(x) sprintf("%.4f", x) # for adjusting precision in plots

# custom function to mutate a subset of rows (a few chunks use it)
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

elem_labels <- c("split_cds_regulatory" = "CDS + Regulatory", 
                 "split_cds_phastcons" = "CDS + phastCons",
                 "merged_cds_regulatory" = "CDS + Regulatory", 
                 "merged_cds_phastcons" = "CDS + phastCons")

mut_labels <- c("roulette" = "Roulette", 
                "carlson" = "Carlson",
                "gnomad" = "Gnomad")

model_labels <- c("split_cds_regulatory" = "CDS + Regulatory", 
                  "split_cds_phastcons" = "CDS + phastCons",
                  "merged_cds_regulatory" = "CDS + Regulatory", 
                  "merged_cds_phastcons" = "CDS + phastCons")

num_chrom <- 22

#################
#
# Helper Functions
#
#################


# ---------- lavaan model selection  ----------
safeLRT <- function(fitA, fitB, robust_method = "satorra.bentler.2010", verbose = FALSE) {
  dfA <- fitMeasures(fitA, "df")
  dfB <- fitMeasures(fitB, "df")
  chisqA <- fitMeasures(fitA, "chisq")
  chisqB <- fitMeasures(fitB, "chisq")
  chisqA_scaled <- fitMeasures(fitA, "chisq.scaled")
  chisqB_scaled <- fitMeasures(fitB, "chisq.scaled")
  
  if(any(is.na(c(dfA, dfB, chisqA, chisqB)))) {
    stop("One or both models lack standard df or chisq. Check fitMeasures(fit, c('df','chisq')).")
  }
  
  if(dfA > dfB) {
    fit_small <- fitA; fit_large <- fitB
    name_small <- deparse(substitute(fitA)); name_large <- deparse(substitute(fitB))
  } else if (dfB > dfA) {
    fit_small <- fitB; fit_large <- fitA
    name_small <- deparse(substitute(fitB)); name_large <- deparse(substitute(fitA))
  }
  
  # robust LRT if both scaled chisq available
  scaled_small <- fitMeasures(fit_small, "chisq.scaled")
  scaled_large <- fitMeasures(fit_large, "chisq.scaled")
  if(!is.na(scaled_small) && !is.na(scaled_large)) {
    if (verbose) message("Using robust LRT (", robust_method, ").")
    return(lavTestLRT(fit_small, fit_large, method = robust_method))
  }
  
  # fallback: standard LRT using unscaled chisq
  chisq_small <- fitMeasures(fit_small, "chisq")
  chisq_large <- fitMeasures(fit_large, "chisq")
  df_small <- fitMeasures(fit_small, "df")
  df_large <- fitMeasures(fit_large, "df")
  
  chi_diff <- chisq_small - chisq_large
  df_diff  <- df_small - df_large
  
  if(is.na(chi_diff) || is.na(df_diff)) {
    stop("Cannot compute fallback LRT: chi-square or df is NA for one of the models.")
  }
  if(df_diff <= 0) {
    stop("Invalid df difference in fallback LRT (df_diff <= 0). Models may not be nested or ordering failed.")
  }
  
  p_value <- pchisq(chi_diff, df_diff, lower.tail = FALSE)
  out <- data.table(chisq_diff = chi_diff, df_diff = df_diff, `Pr(>Chisq)` = p_value)
  rownames(out) <- "chisq"
  if(verbose) message("Using standard (unscaled) LRT fallback.")
  return(out)
}

# ---------- Genomic SEM  ----------
add_genome_offsets_and_midpoint <- function(tbl,
                                            chrom_lengths = NULL,
                                            chrom_col = "chrom",
                                            start_col = "chromStart",
                                            end_col = "chromEnd",
                                            out_midpoint_name = "midpoint",
                                            order_chrom_by = c("natural", "provided", "lexicographic"),
                                            keep_integer_midpoint = FALSE) {
  
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")
  order_chrom_by <- match.arg(order_chrom_by)
  
  dt <- data.table::as.data.table(tbl)
  
  # validate input columns
  if (!all(c(chrom_col, start_col, end_col) %in% names(dt))) {
    stop("tbl must contain columns: ", paste(c(chrom_col, start_col, end_col), collapse = ", "))
  }
  
  # normalize types
  dt[[chrom_col]] <- as.character(dt[[chrom_col]])
  dt[[start_col]] <- as.numeric(dt[[start_col]])
  dt[[end_col]] <- as.numeric(dt[[end_col]])
  
  # build chrom_lengths if not supplied
  if (is.null(chrom_lengths)) {
    chrom_lengths <- dt[, .(length = max(get(end_col), na.rm = TRUE)), by = chrom_col]
    data.table::setnames(chrom_lengths, old = chrom_col, new = "chrom", skip_absent = TRUE)
  } else {
    chrom_lengths <- data.table::as.data.table(chrom_lengths)
    # normalize expected column names
    if (!("chrom" %in% names(chrom_lengths))) {
      if (chrom_col %in% names(chrom_lengths)) data.table::setnames(chrom_lengths, chrom_col, "chrom")
    }
    if (!("length" %in% names(chrom_lengths))) {
      if ("chromEnd" %in% names(chrom_lengths)) data.table::setnames(chrom_lengths, "chromEnd", "length")
      else if ("end" %in% names(chrom_lengths)) data.table::setnames(chrom_lengths, "end", "length")
      else stop("chrom_lengths must provide columns: chrom and length (or chromEnd/end)")
    }
    chrom_lengths[, chrom := as.character(chrom)]
    chrom_lengths[, length := as.numeric(length)]
  }
  
  # determine chromosome order
  chroms <- unique(chrom_lengths$chrom)
  if (order_chrom_by == "provided") {
    chrom_lengths <- chrom_lengths[match(chroms, chrom_lengths$chrom)]
  } else if (order_chrom_by == "natural") {
    extract_num <- function(x) {
      x2 <- sub("^chr", "", x, ignore.case = TRUE)
      suppressWarnings(as.integer(sub("^.*?_?([0-9]+)$", "\\1", x2)))
    }
    nums <- extract_num(chroms)
    has_num <- !is.na(nums)
    # numeric chromosomes first sorted, then the remainder lexicographically
    ord <- c(order(nums[has_num], na.last = NA), order(chroms[!has_num]))
    chrom_ordered <- chroms[ord]
    chrom_lengths <- chrom_lengths[match(chrom_ordered, chrom_lengths$chrom)]
  } else {
    chrom_lengths <- chrom_lengths[order(chrom)]
  }
  
  # compute cumulative offsets such that the first chromosome has offset 0
  chrom_lengths[, chrom_offset := as.numeric(data.table::shift(cumsum(length), fill = 0))]
  
  # join offsets into dt
  setkeyv(chrom_lengths, "chrom")
  dt_temp <- copy(dt)
  data.table::setnames(dt_temp, old = chrom_col, new = "chrom", skip_absent = TRUE)
  dt_temp <- chrom_lengths[dt_temp, on = .(chrom)]
  
  # compute genome-wide coordinates; first chromosome receives offset 0 by design
  dt_temp[, gstart := chrom_offset + get(start_col)]
  dt_temp[, gend   := chrom_offset + get(end_col)]
  dt_temp[, (out_midpoint_name) := (gstart + gend) / 2]
  
  if (isTRUE(keep_integer_midpoint)) {
    dt_temp[, (out_midpoint_name) := floor(get(out_midpoint_name))]
  }
  
  # restore original chrom column name if different
  if (chrom_col != "chrom") data.table::setnames(dt_temp, "chrom", chrom_col)
  
  # assemble output columns (preserve original columns order then append new columns)
  new_cols <- c("chrom_offset", "gstart", "gend", out_midpoint_name)
  # avoid duplicating columns if input already had same names
  new_cols <- new_cols[!new_cols %in% names(tbl)]
  out_cols <- c(names(tbl), new_cols)
  out <- dt_temp[, ..out_cols]
  return(out[])
}

make_lavaan_model <- function(edges,
                              covs = NULL,
                              defined = NULL,
                              groups = NULL,
                              free_by = c("shared", "groups", "named"),
                              per_group_map = NULL,
                              label_prefix = "p",
                              sep = "  ") {
  free_by <- match.arg(free_by)
  stopifnot(is.list(edges) && length(edges) > 0)
  
  # helper to normalize an RHS spec into a list(predictor, label, start)
  normalize_rhs <- function(r) {
    if (is.character(r) && length(r) == 1) return(list(p = r, label = NULL, start = NULL))
    if (is.list(r) && !is.null(r$p)) {
      return(list(p = as.character(r$p),
                  label = if (!is.null(r$label)) as.character(r$label) else NULL,
                  start = if (!is.null(r$start)) as.numeric(r$start) else NULL))
    }
    stop("RHS must be a string or list with element 'p'")
  }
  
  # collect labeled edges in deterministic order (for mapping if needed)
  collect_labels <- function(edges_list) {
    labs <- character(0)
    for (lhs in names(edges_list)) {
      rhss <- edges_list[[lhs]]
      for (r in rhss) {
        nr <- normalize_rhs(r)
        if (!is.null(nr$label)) labs <- c(labs, nr$label)
      }
    }
    unique(labs)
  }
  original_labels <- collect_labels(edges)
  label_to_index <- setNames(seq_along(original_labels), original_labels)
  
  lines <- character(0)
  edge_i <- 0L
  for (lhs in names(edges)) {
    rhss <- edges[[lhs]]
    for (r in rhss) {
      edge_i <- edge_i + 1L
      nr <- normalize_rhs(r)
      p <- nr$p; lab <- nr$label; st <- nr$start
      
      start_part <- if (!is.null(st)) paste0("start(", st, ")*") else ""
      
      if (free_by == "shared") {
        # label or unlabeled (shared across groups)
        if (!is.null(lab)) lab_part <- paste0(lab, "*") else lab_part <- ""
        line <- paste0(lhs, " ~ ", start_part, lab_part, p)
        lines <- c(lines, line)
        
      } else if (free_by == "groups") {
        if (is.null(groups) || groups < 1) {
          # we still allow groups to be NULL; labels remain per-edge
          lab_name <- paste0(label_prefix, edge_i)
        } else {
          lab_name <- paste0(label_prefix, edge_i)
        }
        lab_part <- paste0(lab_name, "*")
        line <- paste0(lhs, " ~ ", start_part, lab_part, p)
        lines <- c(lines, line)
        
      } else { # free_by == "named"
        # per_group_map allows specifying which particular edges should be group-specific
        key <- paste0(lhs, "~", p)
        if (!is.null(per_group_map) && !is.na(per_group_map[key]) && per_group_map[key] == "groups") {
          if (is.null(groups) || groups < 1) stop("groups must be provided for per-group labels")
          lab_names <- paste0(label_prefix, edge_i, "_g", seq_len(groups))
          lab_vec <- paste0("c(", paste(lab_names, collapse = ", "), ")*")
          line <- paste0(lhs, " ~ ", start_part, lab_vec, p)
          lines <- c(lines, line)
        } else {
          # otherwise give a unique per-edge label (shared across groups)
          lab_name <- if (!is.null(lab)) lab else paste0(label_prefix, edge_i)
          lab_part <- paste0(lab_name, "*")
          line <- paste0(lhs, " ~ ", start_part, lab_part, p)
          lines <- c(lines, line)
        }
      }
    }
  }
  
  # covariances
  if (!is.null(covs) && length(covs) > 0) {
    for (cv in covs) {
      if (is.character(cv) && length(cv) == 2) {
        lines <- c(lines, paste0(cv[1], " ~~ ", cv[2]))
      } else if (is.list(cv) && length(cv) == 2) {
        lines <- c(lines, paste0(cv[[1]], " ~~ ", cv[[2]]))
      } else {
        stop("covs must be pairs of variable names")
      }
    }
  }
  
  # defined parameters
  if (!is.null(defined) && length(defined) > 0) {
    defined_lines <- character(0)
    for (nm in names(defined)) {
      expr <- defined[[nm]]
      # if original labels (like c1,c2,...) exist and we've generated per-edge labels,
      # remap occurrences of original labels to the new per-edge label names where possible.
      for (orig_lab in original_labels) {
        idx <- label_to_index[[orig_lab]]
        if (!is.null(idx) && !is.na(idx)) {
          new_lab <- paste0(label_prefix, idx)
          pattern <- paste0("\\b", orig_lab, "\\b")
          expr <- gsub(pattern, new_lab, expr, perl = TRUE)
        }
      }
      defined_lines <- c(defined_lines, paste0(nm, " := ", expr))
    }
    lines <- c(lines, defined_lines)
  }
  
  paste(lines, collapse = paste0("\n", sep))
}


make_lavaan_model_free_with_grouped_defined <- function(edges,
                                                        covs = NULL,
                                                        defined = NULL,
                                                        groups,
                                                        per_group_map = NULL,
                                                        label_prefix = "x",
                                                        defined_suffix_sep = "_g",
                                                        sep = "  ") {
  stopifnot(is.list(edges) && length(edges) > 0)
  stopifnot(!missing(groups) && is.numeric(groups) && groups >= 1)
  
  normalize_rhs <- function(r) {
    if (is.character(r) && length(r) == 1) return(list(p = r, label = NULL, start = NULL))
    if (is.list(r) && !is.null(r$p)) {
      return(list(p = as.character(r$p),
                  label = if (!is.null(r$label)) as.character(r$label) else NULL,
                  start = if (!is.null(r$start)) as.numeric(r$start) else NULL))
    }
    stop("RHS must be a string or list with element 'p'")
  }
  
  # collect original labeled tokens in deterministic order
  orig_labels <- character(0)
  for (lhs in names(edges)) for (r in edges[[lhs]]) {
    nr <- normalize_rhs(r)
    if (!is.null(nr$label)) orig_labels <- c(orig_labels, nr$label)
  }
  orig_labels <- unique(orig_labels)
  label_index <- setNames(seq_along(orig_labels), orig_labels)
  
  lines <- character(0)
  edge_i <- 0L
  
  # build regression and covariance lines
  for (lhs in names(edges)) {
    rhss <- edges[[lhs]]
    for (r in rhss) {
      edge_i <- edge_i + 1L
      nr <- normalize_rhs(r)
      p <- nr$p; lab <- nr$label; st <- nr$start
      start_part <- if (!is.null(st)) paste0("start(", st, ")*") else ""
      key <- paste0(lhs, "~", p)
      is_group_edge <- !is.null(per_group_map) && !is.na(per_group_map[key]) && per_group_map[key] == "groups"
      
      if (is_group_edge) {
        lab_names <- paste0(label_prefix, edge_i, defined_suffix_sep, seq_len(groups))
        lab_vec <- paste0("c(", paste(lab_names, collapse = ", "), ")*")
        line <- paste0(lhs, " ~ ", start_part, lab_vec, p)
      } else {
        lab_name <- if (!is.null(lab)) lab else paste0(label_prefix, edge_i)
        line <- paste0(lhs, " ~ ", start_part, lab_name, "*", p)
      }
      lines <- c(lines, line)
    }
  }
  
  # covariances
  if (!is.null(covs) && length(covs) > 0) {
    for (cv in covs) {
      if (is.character(cv) && length(cv) == 2) {
        lines <- c(lines, paste0(cv[1], " ~~ ", cv[2]))
      } else if (is.list(cv) && length(cv) == 2) {
        lines <- c(lines, paste0(cv[[1]], " ~~ ", cv[[2]]))
      } else stop("covs must be pairs")
    }
  }
  
  # expand defined expressions into per-group scalar definitions
  if (!is.null(defined) && length(defined) > 0) {
    # helper: build per-group vector tokens for an original label if it was emitted as group-specific
    build_vec_tokens <- function(orig_lab, edge_index) {
      paste0(label_prefix, edge_index, defined_suffix_sep, seq_len(groups))
    }
    
    # we need to know which edge index corresponds to each original label
    # walk edges in same order to map label -> edge_index (first occurrence)
    label_to_edge_index <- list()
    edge_i2 <- 0L
    for (lhs in names(edges)) {
      for (r in edges[[lhs]]) {
        edge_i2 <- edge_i2 + 1L
        nr <- normalize_rhs(r)
        if (!is.null(nr$label) && !(nr$label %in% names(label_to_edge_index))) {
          label_to_edge_index[[nr$label]] <- edge_i2
        }
      }
    }
    
    for (nm in names(defined)) {
      expr <- defined[[nm]]
      # find occurrences of original labels and replace with per-group token names
      # but we produce explicit scalar defined lines per group
      for (g in seq_len(groups)) {
        expr_g <- expr
        for (orig in orig_labels) {
          if (orig %in% names(label_to_edge_index)) {
            ei <- label_to_edge_index[[orig]]
            token_g <- paste0(label_prefix, ei, defined_suffix_sep, g)
            expr_g <- gsub(paste0("\\b", orig, "\\b"), token_g, expr_g, perl = TRUE)
          } else {
            # if original label not found, leave token unchanged (rare)
            expr_g <- gsub(paste0("\\b", orig, "\\b"), orig, expr_g, perl = TRUE)
          }
        }
        def_name <- paste0(nm, defined_suffix_sep, g)
        lines <- c(lines, paste0(def_name, " := ", expr_g))
      }
    }
  }
  
  paste(lines, collapse = paste0("\n", sep))
}
# ---------- Generic ----------
`%||%` <- function(a, b) if (!is.null(a)) a else b

std_by_sd <- function(x, eps = 1e-8) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s < eps) s <- eps
  (x - m) / s
}

safe_std <- function(x) {
  m <- median(x, na.rm=TRUE)
  s <- stats::mad(x, constant=1, na.rm=TRUE)
  if (is.na(s) || s < 1e-8) s <- sd(x, na.rm=TRUE) + 1e-8
  (x - m) / s
}

quant_norm <- function(x) {
  r <- rank(x, na.last = "keep")
  q <- qnorm((r - 0.5) / sum(!is.na(x)))
  (q - median(q, na.rm = TRUE)) / (stats::mad(q, constant = 1, na.rm = TRUE) + 1e-8)
}

# ---------- Blocking ----------
define_blocks <- function(pos, chrom = NULL, max_gap = 1e5, max_span = 2e6, min_size = 1L) {
  if (is.null(pos) || length(pos) == 0) return(integer(0))
  ok <- !is.na(pos)
  out <- rep(NA_integer_, length(pos))
  if (!is.null(chrom)) {
    chroms <- unique(chrom[ok])
    offset <- 0L
    for (c in chroms) {
      idx <- which(chrom == c & ok)
      if (length(idx) == 0) next
      subpos <- pos[idx]
      subout <- define_blocks(subpos, chrom = NULL, max_gap = max_gap, max_span = max_span, min_size = min_size)
      subout2 <- as.integer(as.factor(subout)) + offset
      out[idx] <- subout2
      offset <- max(offset, max(subout2, na.rm = TRUE))
    }
    return(out)
  }
  
  ord <- order(pos)
  pos_o <- pos[ord]
  blk <- integer(length(pos_o))
  cur <- 1L; blk[1] <- cur; start <- pos_o[1]
  if (length(pos_o) > 1) {
    for (i in 2:length(pos_o)) {
      if (is.na(pos_o[i]) || is.na(pos_o[i - 1])) {
        cur <- cur + 1L
        start <- pos_o[i]
        blk[i] <- cur
        next
      }
      if ((pos_o[i] - pos_o[i - 1] > max_gap) || ((pos_o[i] - start) > max_span)) {
        cur <- cur + 1L
        start <- pos_o[i]
      }
      blk[i] <- cur
    }
  }
  out[ord] <- blk
  
  if (!is.null(min_size) && min_size > 1L) {
    tb <- table(out[ok])
    small_blks <- as.integer(names(tb)[tb < min_size])
    if (length(small_blks) > 0) {
      for (sb in small_blks) {
        idsb <- which(out == sb)
        if (length(idsb) == 0) next
        pos_sb <- pos[idsb]
        neigh_idx <- which(out != sb & !is.na(out))
        if (length(neigh_idx) == 0) next
        dists <- abs(pos[neigh_idx] - median(pos_sb, na.rm = TRUE))
        pick <- neigh_idx[which.min(dists)]
        out[idsb] <- out[pick]
      }
      out[ok] <- as.integer(as.factor(out[ok]))
    }
  }
  out
}

# ---------- Decay estimation ----------
estimate_decay <- function(resid, pos, n_bins = 50, model = "exponential") {
  if (length(resid) == 0 || length(pos) == 0) return(NULL)
  ok <- !is.na(resid) & !is.na(pos)
  if (sum(ok) < 4) return(NULL)
  tryCatch({
    samp_idx <- if (sum(ok) > 2000) sample(which(ok), 2000) else which(ok)
    pos_s <- pos[samp_idx]
    resid_s <- resid[samp_idx]
    dmat <- abs(outer(pos_s, pos_s, "-"))
    sill <- var(resid[ok], na.rm = TRUE)
    if (is.na(sill) || sill <= 0) sill <- stats::mad(resid[ok], constant = 1, na.rm = TRUE)^2
    if (is.na(sill) || sill <= 0) sill <- 1e-8
    med_d <- median(dmat[lower.tri(dmat)], na.rm = TRUE)
    phi_est <- max(med_d, 1)
    nugget <- 0
    cov_fun <- function(h) sill * exp(- (h) / phi_est)
    list(method = "heuristic", phi = phi_est, sill = sill, nugget = nugget, cov_fun = cov_fun, converged = FALSE)
  }, error = function(e) NULL)
}

# ---------- Covariance matrix builder ----------
build_cov_matrix <- function(pos, cov_fun, nugget = 0, ridge = 1e-8, taper_radius = NULL, as_sparse = TRUE) {
  if (length(pos) == 0) return(NULL)
  ok <- !is.na(pos)
  if (sum(ok) == 0) return(NULL)
  tryCatch({
    idx <- which(ok)
    pos_ok <- pos[idx]
    n <- length(pos_ok)
    D <- as.matrix(stats::dist(pos_ok))
    C <- cov_fun(D)
    if (!is.null(nugget) && nugget > 0) diag(C) <- diag(C) + nugget
    if (!is.null(taper_radius) && is.numeric(taper_radius) && taper_radius > 0) {
      Tmat <- pmax(0, 1 - D / taper_radius)
      C <- C * Tmat
    }
    C <- (C + t(C)) / 2
    diag(C) <- diag(C) + ridge
    if (as_sparse) {
      thr <- max(diag(C)) * 1e-12
      C[abs(C) < thr] <- 0
      C_sparse <- Matrix::Matrix(C, sparse = TRUE)
      list(C = C_sparse, idx = idx, n = n, pos_range = range(pos_ok, na.rm = TRUE))
    } else {
      list(C = C, idx = idx, n = n, pos_range = range(pos_ok, na.rm = TRUE))
    }
  }, error = function(e) NULL)
}

# ---------- GLS solver ----------
gls_solve <- function(y, X, C_sp) {
  if (length(y) == 0) return(NULL)
  ok <- !is.na(y)
  if (sum(ok) == 0) return(NULL)
  tryCatch({
    y2 <- y[ok]
    X2 <- if (is.matrix(X)) X[ok, , drop = FALSE] else as.matrix(X[ok, , drop = FALSE])
    if (!is.null(C_sp$idx)) {
      map_idx <- match(which(ok), C_sp$idx)
      if (any(is.na(map_idx))) {
        covered_rows <- which(!is.na(map_idx))
        if (length(covered_rows) < 4) return(NULL)
        y2 <- y2[covered_rows]; X2 <- X2[covered_rows, , drop = FALSE]
        map_idx2 <- map_idx[covered_rows]
        C_sub <- C_sp$C[map_idx2, map_idx2, drop = FALSE]
      } else {
        C_sub <- C_sp$C
      }
    } else {
      C_sub <- C_sp$C
    }
    
    if (inherits(C_sub, "sparseMatrix")) {
      chol_try <- try(Matrix::Cholesky(C_sub, LDL = FALSE), silent = TRUE)
      if (inherits(chol_try, "try-error")) return(NULL)
      CinvX <- try(Matrix::solve(chol_try, X2), silent = TRUE)
      Cinvy <- try(Matrix::solve(chol_try, y2), silent = TRUE)
      if (inherits(CinvX, "try-error") || inherits(Cinvy, "try-error")) return(NULL)
      XtCinvX <- crossprod(as.matrix(CinvX))
      if (det(XtCinvX) == 0 || any(is.na(XtCinvX))) {
        beta <- try(MASS::ginv(XtCinvX) %*% crossprod(as.matrix(CinvX), as.numeric(Cinvy)), silent = TRUE)
        if (inherits(beta, "try-error")) return(NULL)
      } else {
        beta <- try(solve(XtCinvX, crossprod(as.matrix(CinvX), as.numeric(Cinvy))), silent = TRUE)
        if (inherits(beta, "try-error")) return(NULL)
      }
      fitted <- as.numeric(X2 %*% beta)
      resid <- as.numeric(y2 - fitted)
      fitted_full <- rep(NA_real_, length(y)); resid_full <- rep(NA_real_, length(y))
      if (exists("covered_rows")) {
        ok_idx <- which(ok)[covered_rows]
        fitted_full[ok_idx] <- fitted; resid_full[ok_idx] <- resid
      } else {
        fitted_full[ok] <- fitted; resid_full[ok] <- resid
      }
      diagnostics <- list(p = ncol(X2), used_n = sum(!is.na(fitted_full)), rankX = qr(X2)$rank)
      list(fitted = fitted_full, resid = resid_full, beta = as.numeric(beta), chol = chol_try, diagnostics = diagnostics)
    } else {
      C_dense <- as.matrix(C_sub)
      Cinv <- try(solve(C_dense), silent = TRUE)
      if (inherits(Cinv, "try-error")) return(NULL)
      XtCinvX <- crossprod(X2, Cinv %*% X2)
      beta <- try(solve(XtCinvX, crossprod(X2, Cinv %*% y2)), silent = TRUE)
      if (inherits(beta, "try-error")) return(NULL)
      fitted <- as.numeric(X2 %*% beta)
      resid <- as.numeric(y2 - fitted)
      fitted_full <- rep(NA_real_, length(y)); resid_full <- rep(NA_real_, length(y))
      fitted_full[ok] <- fitted; resid_full[ok] <- resid
      diagnostics <- list(p = ncol(X2), used_n = sum(!is.na(fitted_full)), rankX = qr(X2)$rank)
      list(fitted = fitted_full, resid = resid_full, beta = as.numeric(beta), diagnostics = diagnostics)
    }
  }, error = function(e) NULL)
}

# ---------- Process a single block ----------
process_block <- function(df_block, response, pos_col = "pos", covariates = NULL,
                          model = "exponential", bin_n = 50, taper_radius = NULL,
                          ridge = 1e-8, nugget_override = NULL, min_block_size = 6,
                          persist_dir = NULL, verbose = FALSE, block_id = NULL) {
  res <- list(n = NA_integer_, method = NA_character_, success = FALSE, warnings = character(0))
  if (!is.data.frame(df_block)) { res$warnings <- c(res$warnings, "df_block not a data.frame"); return(res) }
  if (!(response %in% names(df_block))) { res$warnings <- c(res$warnings, "response missing"); return(res) }
  if (!(pos_col %in% names(df_block))) { res$warnings <- c(res$warnings, "pos_col missing"); return(res) }
  pos <- df_block[[pos_col]]; y <- df_block[[response]]
  n <- length(y); res$n <- n
  if (all(is.na(y))) {
    res$warnings <- c(res$warnings, "response all NA in block")
    res$method <- "no_data_block"
    res$fitted <- rep(NA_real_, n); res$resid <- rep(NA_real_, n)
    return(res)
  }
  
  X0 <- NULL
  if (!is.null(covariates)) {
    present <- covariates[covariates %in% names(df_block)]
    if (length(present) > 0) X0 <- df_block[, present, drop = FALSE]
    if (!is.null(X0) && all(sapply(X0, function(col) all(is.na(col))))) X0 <- NULL
  }
  
  ols <- tryCatch({
    if (is.null(X0)) stats::lm(y ~ 1) else {
      mm <- try(stats::model.matrix(~ ., data = X0), silent = TRUE)
      if (inherits(mm, "try-error") || ncol(mm) == 0) stats::lm(y ~ 1) else stats::lm.fit(mm, y)
    }
  }, error = function(e) NULL)
  if (is.null(ols)) {
    res$warnings <- c(res$warnings, "OLS fit failed; aborting block")
    res$method <- "OLS_failed"
    res$fitted <- rep(NA_real_, n); res$resid <- rep(NA_real_, n)
    return(res)
  }
  fitted0 <- tryCatch(if (inherits(ols, "lm")) stats::fitted(ols) else as.numeric(ols$fitted.values), error = function(e) rep(NA_real_, n))
  resid0  <- tryCatch(if (inherits(ols, "lm")) stats::resid(ols) else as.numeric(ols$residuals), error = function(e) rep(NA_real_, n))
  
  if (sum(!is.na(resid0)) < min_block_size) {
    res$warnings <- c(res$warnings, "small block; using OLS residuals")
    res$method <- "OLS_small_block"; res$fitted <- fitted0; res$resid <- resid0
    return(res)
  }
  
  dec <- estimate_decay(resid0, pos, n_bins = bin_n, model = model)
  if (is.null(dec)) {
    res$warnings <- c(res$warnings, "decay estimation failed; using pooled fallback")
    dec <- list(method = "fallback", nugget = 0, sill = var(resid0, na.rm = TRUE),
                cov_fun = function(h) var(resid0, na.rm = TRUE) * exp(-h / max(diff(range(pos, na.rm = TRUE)), 1)),
                converged = FALSE)
  }
  res$decay <- dec
  nugget <- if (!is.null(nugget_override)) nugget_override else dec$nugget %||% 0
  
  C_sp <- build_cov_matrix(pos, cov_fun = dec$cov_fun, nugget = nugget, ridge = ridge, taper_radius = taper_radius, as_sparse = TRUE)
  if (is.null(C_sp)) {
    res$warnings <- c(res$warnings, "covariance construction failed; returning OLS residuals")
    res$method <- "OLS_cov_failed"; res$fitted <- fitted0; res$resid <- resid0
    return(res)
  }
  
  Xmat <- if (is.null(X0)) matrix(1, n, 1) else cbind(1, as.matrix(X0))
  gls_out <- gls_solve(y = y, X = Xmat, C_sp = C_sp)
  if (is.null(gls_out)) {
    res$warnings <- c(res$warnings, "GLS solve failed; falling back to OLS residuals")
    res$method <- "GLS_failed"; res$fitted <- fitted0; res$resid <- resid0
    return(res)
  }
  
  resid_gls <- gls_out$resid; fitted_gls <- gls_out$fitted
  if (all(is.na(resid_gls))) {
    res$warnings <- c(res$warnings, "GLS produced all-NA residuals")
    res$method <- "GLS_allNA"; res$fitted <- fitted_gls; res$resid <- resid_gls
    return(res)
  }
  
  res$method <- "GLS"
  res$fitted <- fitted_gls
  res$resid <- resid_gls
  res$gls_diagnostics <- gls_out$diagnostics %||% list()
  res$success <- TRUE
  
  res
}

# ---------- Per-variable pipeline ----------
gls_pipeline_full <- function(df, response, pos = "midpoint", chrom = "chrom",
                              covariates = NULL, block_max_gap = 1e5, block_max_span = 2e6,
                              model = "exponential", bin_n = 50, taper_radius = NULL,
                              ridge = 1e-8, nugget_override = NULL, min_block_size = 6,
                              persist_dir = NULL, parallel = FALSE, verbose = TRUE) {
  if (!is.data.frame(df)) stop("df must be a data.frame")
  if (!(response %in% names(df))) stop("response not found")
  if (!(pos %in% names(df))) stop("pos not found")
  
  chrom_vec <- if (!is.null(chrom) && chrom %in% names(df)) df[[chrom]] else NULL
  blocks <- define_blocks(df[[pos]], chrom = chrom_vec, max_gap = block_max_gap, max_span = block_max_span)
  if(length(blocks) != nrow(df) || all(is.na(blocks))) {
    if(verbose) message("define_blocks returned invalid length or all NA; using single block")
    blocks <- rep(1L, nrow(df))
  }
  
  fitted_gls <- rep(NA_real_, nrow(df))
  resid_gls <- rep(NA_real_, nrow(df))
  block_fits <- list()
  
  unique_blocks <- sort(unique(blocks[!is.na(blocks)]))
  if(verbose) message("Processing ", length(unique_blocks), " blocks for response '", response, "'")
  
  for(b in unique_blocks) {
    idx <- which(blocks == b)
    df_block <- df[idx, , drop = FALSE]
    pb <- tryCatch(
      process_block(
        df_block, response = response, pos_col = pos, covariates = covariates,
        model = model, bin_n = bin_n, taper_radius = taper_radius,
        ridge = ridge, nugget_override = nugget_override, min_block_size = min_block_size,
        persist_dir = NULL, verbose = verbose, block_id = b
      ),
      error = function(e) {
        if(verbose) message("Error in block ", b, ": ", conditionMessage(e))
        list(n = length(idx), method = "error", warnings = c(conditionMessage(e)),
             fitted = rep(NA_real_, length(idx)), resid = rep(NA_real_, length(idx)))
      }
    )
    fitted_gls[idx] <- if(!is.null(pb$fitted)) pb$fitted else rep(NA_real_, length(idx))
    resid_gls[idx]  <- if(!is.null(pb$resid)) pb$resid else rep(NA_real_, length(idx))
    block_fits[[as.character(b)]] <- pb
  }
  
  out <- as.data.frame(df)
  out$fitted_gls <- fitted_gls
  out$resid_gls <- resid_gls
  attr(out, "block_fits") <- block_fits
  out
}

# ---------- Wrapper ----------
gls_preprocess_variables <- function(df, vars, pos = "midpoint", chrom = "chrom",
                                     covariates = NULL,
                                     block_max_gap = 1e5, block_max_span = 2e6,
                                     model = "exponential", bin_n = 50,
                                     taper_radius = NULL, ridge = 1e-8, nugget_override = NULL,
                                     min_block_size = 6, persist_dir = NULL,
                                     parallel = FALSE, future_workers = 4, verbose = TRUE) {
  if (!is.data.frame(df)) stop("df must be a data.frame or data.table")
  if (!is.character(vars) || length(vars) == 0) stop("vars must be a non-empty character vector")
  missing_vars <- setdiff(vars, names(df))
  if (length(missing_vars) > 0) stop("Some vars not found in df: ", paste(missing_vars, collapse = ", "))
  if (!pos %in% names(df)) stop("pos column not found: ", pos)
  if (!is.null(chrom) && !chrom %in% names(df)) stop("chrom column not found: ", chrom)
  
  nrow_orig <- nrow(df)
  full_df <- df
  diagnostics <- list(global = list(call = match.call(),
                                    pos = pos, chrom = chrom, vars = vars,
                                    block_max_gap = block_max_gap, block_max_span = block_max_span,
                                    model = model, taper_radius = taper_radius, ridge = ridge,
                                    min_block_size = min_block_size,
                                    parallel = parallel, future_workers = future_workers))
  diagnostics$variables <- vector("list", length(vars))
  names(diagnostics$variables) <- vars
  
  proc_one_var <- function(varname) {
    if(verbose) message("Processing variable: ", varname)
    valid_count <- sum(!is.na(df[[varname]]))
    if(valid_count == 0) {
      warning("Variable '", varname, "' has zero non-NA values; skipping and producing NA outputs")
      n <- nrow_orig
      full_cols <- list()
      full_cols[[paste0(varname, "_fitted_gls")]] <- rep(NA_real_, n)
      full_cols[[paste0(varname, "_resid_gls")]]   <- rep(NA_real_, n)
      diag <- list(var = varname, skipped = TRUE, reason = "all NA", method_table = NA, fallback_count = NA, acf_summary = list(), block_fits = list())
      return(list(success = TRUE, full_cols = full_cols, diagnostics = diag))
    }
    
    out_try <- try(gls_pipeline_full(df = df, response = varname, pos = pos, chrom = chrom,
                                     covariates = covariates, block_max_gap = block_max_gap, block_max_span = block_max_span,
                                     model = model, bin_n = bin_n, taper_radius = taper_radius, ridge = ridge,
                                     nugget_override = nugget_override, min_block_size = min_block_size,
                                     persist_dir = NULL, parallel = FALSE, verbose = verbose), silent = TRUE)
    
    if (inherits(out_try, "try-error")) {
      warning("gls_pipeline_full failed for variable '", varname, "': ", conditionMessage(attr(out_try, "condition")))
      n <- nrow_orig
      full_cols <- list()
      full_cols[[paste0(varname, "_fitted_gls")]] <- rep(NA_real_, n)
      full_cols[[paste0(varname, "_resid_gls")]]   <- rep(NA_real_, n)
      diag <- list(var = varname, error = conditionMessage(attr(out_try, "condition")), skipped = TRUE, block_fits = list())
      return(list(success = FALSE, full_cols = full_cols, diagnostics = diag))
    }
    
    out <- out_try
    n <- nrow(out)
    get_or_na <- function(x, name) if (name %in% names(x)) x[[name]] else rep(NA_real_, n)
    fitted_col <- get_or_na(out, "fitted_gls")
    resid_col  <- get_or_na(out, "resid_gls")
    
    full_cols <- list()
    full_cols[[paste0(varname, "_fitted_gls")]] <- fitted_col
    full_cols[[paste0(varname, "_resid_gls")]]   <- resid_col
    
    block_fits <- attr(out, "block_fits")
    if (is.null(block_fits)) block_fits <- list()
    methods_vec <- tryCatch({
      v <- sapply(block_fits, function(x) if (!is.null(x$method)) as.character(x$method) else NA_character_)
      if (length(v) == 0) character(0) else v
    }, error = function(e) character(0))
    method_table <- if (length(methods_vec) == 0) NA else table(methods_vec)
    fallback_count <- NA_integer_
    if (!is.na(method_table)[1]) {
      mt_names <- names(method_table)
      sel <- grepl("OLS|fallback", mt_names, ignore.case = TRUE)
      if (any(sel)) fallback_count <- sum(method_table[sel]) else fallback_count <- 0L
    }
    
    acf_summary <- lapply(seq_along(block_fits), function(i) {
      bf <- block_fits[[i]]
      if (!is.null(bf$resid)) {
        acf_vals <- try(stats::acf(bf$resid, plot = FALSE, na.action = na.pass)$acf, silent = TRUE)
        if (inherits(acf_vals, "try-error")) return(list(block = i, ok = FALSE))
        return(list(block = i, ok = TRUE, acf1 = as.numeric(acf_vals[2]), n = if (!is.null(bf$n)) bf$n else NA))
      } else {
        return(list(block = i, ok = FALSE, reason = "no resid"))
      }
    })
    
    diag <- list(var = varname, method_table = method_table, fallback_count = fallback_count,
                 acf_summary = acf_summary, block_fits = block_fits)
    list(success = TRUE, full_cols = full_cols, diagnostics = diag)
  }
  
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) stop("future.apply required for parallel = TRUE; install it or set parallel = FALSE")
    res_list <- future.apply::future_lapply(vars, proc_one_var, future.seed = TRUE)
  } else {
    res_list <- lapply(vars, proc_one_var)
  }
  names(res_list) <- vars
  
  for (v in vars) {
    r <- res_list[[v]]
    cols <- r$full_cols
    for (nm in names(cols)) full_df[[nm]] <- cols[[nm]]
    diagnostics$variables[[v]] <- r$diagnostics
  }
  
  sem_df <- data.frame(matrix(NA_real_, nrow = nrow_orig, ncol = length(vars)))
  colnames(sem_df) <- vars
  for (i in seq_along(vars)) {
    v <- vars[i]
    resid_col_name <- paste0(v, "_resid_gls")
    if (resid_col_name %in% names(full_df)) {
      sem_df[[v]] <- full_df[[resid_col_name]]
    } else {
      sem_df[[v]] <- rep(NA_real_, nrow_orig)
    }
  }
  
  attr(full_df, "gls_preproc_diagnostics") <- diagnostics
  list(sem_df = sem_df, full_df = full_df, diagnostics = diagnostics)
}

preprocess_sem_from_lavaan <- function(model,
                                       df,
                                       pos = "midpoint",
                                       chrom = "chrom",
                                       covariates = NULL,
                                       taper_radius = 5e5,
                                       block_max_gap = 2e5,
                                       block_max_span = 2e6,
                                       ridge = 1e-8,
                                       parallel = FALSE,
                                       verbose = TRUE,
                                       apply_std = TRUE,
                                       gls_extra_args = list()) {
  # model: lavaan model string
  # df: data.frame / data.table with SEM variables and pos/chrom columns
  # apply_std: if TRUE, create mean/SD standardized sem_df_std (returned)
  # gls_extra_args: extra named args forwarded to gls_preprocess_variables
  if(!is.character(model) || length(model) == 0) stop("model must be a lavaan model string")
  if(!is.data.frame(df)) stop("df must be a data.frame or data.table")
  
  # helper: simple mean/SD standardizer with tiny-EPS guard
  std_by_sd <- function(x, eps = 1e-8) {
    if(all(is.na(x))) return(rep(NA_real_, length(x)))
    m <- mean(x, na.rm = TRUE)
    s <- stats::sd(x, na.rm = TRUE)
    if(is.na(s) || s < eps) s <- eps
    (x - m) / s
  }
  
  # extract variable names from the lavaan model
  vars_in_model <- NULL
  if(requireNamespace("lavaan", quietly = TRUE)) {
    parsed <- tryCatch(lavaan::lavaanify(model), error = function(e) NULL)
    if (!is.null(parsed) && is.data.frame(parsed)) {
      sel_ops <- parsed$op %in% c("~", "=~", "~~", ":=")
      lhs_vars <- unique(as.character(parsed$lhs[sel_ops]))
      rhs_vars <- unique(as.character(parsed$rhs[sel_ops]))
      vars_in_model <- unique(c(lhs_vars, rhs_vars))
    }
  }
  
  # fallback regex-based extraction if lavaan parse failed
  if(is.null(vars_in_model) || length(vars_in_model) == 0) {
    mlines <- unlist(strsplit(model, "\n"))
    mlines <- trimws(gsub("#.*", "", mlines))
    mlines <- mlines[mlines != ""]
    toks <- character(0)
    for(ln in mlines) {
      parts <- strsplit(ln, split = "(~|=~|~~|:=)", perl = TRUE)[[1]]
      for(p in parts) {
        p2 <- gsub("[A-Za-z0-9_\\.]+\\s*\\([^\\)]*\\)", " ", p)
        p2 <- gsub("[\\+\\-\\*:\\(\\)\\,\\=\\>\\<\\?\\!\\;\\|/\\\\]", " ", p2)
        p2 <- gsub("\\s+", " ", trimws(p2))
        if (nchar(p2) > 0) toks <- c(toks, unlist(strsplit(p2, "\\s+")))
      }
    }
    toks <- toks[nzchar(toks)]
    toks <- toks[!grepl("^[0-9\\.-]+$", toks)]
    toks <- unique(toks)
    vars_in_model <- toks
  }
  
  # intersect with available columns in df
  vars_found <- intersect(vars_in_model, names(df))
  if(length(vars_found) == 0) stop("No model variables found in df columns. Check model and df names.")
  if(verbose) message("Identified SEM variables (present in df): ", paste(sort(vars_found), collapse = ", "))
  
  # call gls_preprocess_variables on the found variables
  call_args <- c(list(df = df,
                      vars = vars_found,
                      pos = pos,
                      chrom = chrom,
                      covariates = covariates,
                      taper_radius = taper_radius,
                      block_max_gap = block_max_gap,
                      block_max_span = block_max_span,
                      ridge = ridge,
                      parallel = parallel,
                      verbose = verbose),
                 gls_extra_args)
  
  res <- do.call(gls_preprocess_variables, call_args)
  
  if (!is.list(res) || !all(c("sem_df", "full_df") %in% names(res))) {
    stop("gls_preprocess_variables did not return expected structure (sem_df, full_df).")
  }
  
  sem_df_raw <- res$sem_df
  full_df <- res$full_df
  diagnostics <- res$diagnostics %||% attr(full_df, "gls_preproc_diagnostics")
  
  # Optionally standardize each sem variable by mean/SD
  sem_df_std <- NULL
  if (isTRUE(apply_std)) {
    sem_df_std <- as.data.frame(lapply(sem_df_raw, function(x) std_by_sd(x)))
    names(sem_df_std) <- names(sem_df_raw)
    attr(sem_df_std, "std_method") <- "mean_sd"
  }
  
  # Attach metadata and return
  attr(sem_df_raw, "model_vars") <- vars_found
  attr(full_df, "model_vars") <- vars_found
  result <- list(sem_df_raw = sem_df_raw,
                 full_df = full_df,
                 diagnostics = diagnostics)
  if (!is.null(sem_df_std)) result$sem_df_std <- sem_df_std
  
  return(result)
}
