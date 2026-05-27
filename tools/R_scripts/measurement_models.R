#################
# 
# This script tables with mutation, recombination and B-values
# then builds measurement models and writes latent predictors to files
# NOTE: this script should be run from within each current_models dir
#
#################

source("~/Data/bgs_lmr/human_data/tools/HelperFunctions.R")

constrained_models <- c("merged_cds_phastcons", "split_cds_phastcons", "split_cds_regulatory", "merged_cds_regulatory")

num_chrom <- 22

# loading data
human_maps_10kb <- fread("maps/maps_10kb.csv.gz")
human_maps_100kb <- fread("maps/maps_100kb.csv.gz")
human_maps_1Mb <- fread("maps/maps_1Mb.csv.gz")

#################
#
# Mutation 
#
#################

# we build the measurement models from the (filtered) 10 kb tables instead of
# from the raw mutation maps to match the sites used to compute diversity
mut_maps_10kb <- filter(human_maps_10kb, constrained_model==constrained_models[1]) %>% # arbitrary constrained model
  dplyr::select(., c(chrom, chromStart, chromEnd, avg_mut, mut_map)) %>%
  pivot_wider(., values_from=avg_mut, names_from=mut_map) 

mut_maps_100kb <- filter(human_maps_100kb, constrained_model==constrained_models[1]) %>% # arbitrary constrained model
  dplyr::select(., c(chrom, chromStart, chromEnd, avg_mut, mut_map)) %>%
  pivot_wider(., values_from=avg_mut, names_from=mut_map) 

mut_maps_1Mb <- filter(human_maps_1Mb, constrained_model==constrained_models[1]) %>% # arbitrary constrained model
  dplyr::select(., c(chrom, chromStart, chromEnd, avg_mut, mut_map)) %>%
  pivot_wider(., values_from=avg_mut, names_from=mut_map) 

################# testing measurement model on chr1, 1 Mb scale ################# 
test_chr1 <- FALSE
if(test_chr1) {
  names(mut_maps_1Mb)[1] <- "chrom"
  human_maps_1Mb <- fread("maps/maps_1Mb.csv.gz")
  
  obs_pi_1Mb <- filter(human_maps_1Mb, constrained_model=="phastcons", mut_map=="roulette") %>% # arbitrary
    dplyr::select(., c(chrom, chromStart, chromEnd, coverage, avg_pi))
  
  df_1Mb <- left_join(mut_maps_1Mb, obs_pi_1Mb, by=c("chrom", "chromStart", "chromEnd")) %>%
    distinct(., chrom, chromStart, .keep_all=T)
  
  # ismc chr 1, from the book chapter Barroso & Dutheil 2026
  ismc_chr1 <- fread("~/Data/bgs_lmr/human_data/data/Mbuti_chr1.theta.1Mb.bedgraph")
  ismc_chr1$chrom <- 1
  ismc_chr1$chrom <- as.integer(ismc_chr1$chrom)
  names(ismc_chr1)[4:9] <- c(paste("ismc", 1:4, sep="_"), "ismc_mean", "ismc_joint")
  
  ismc_measurement_1 <- "
  ismc=~ismc_1 + ismc_2 + ismc_3 + ismc_4
  "
  
  m1 <- sem(model=ismc_measurement_1,
            data=ismc_chr1,
            estimator="MLR",
            optim.method="nlminb",
            std.ov=T,
            std.lv=T)
  
  summary(m1)
  inspect(m1, "r2")
  parameterEstimates(m1)
  fitMeasures(m1)
  # p-value = 0.206
  # CFI = 0.998
  # GFI = 0.994
  
  ismc_measurement_2 <- "
  ismc=~ismc_1 + ismc_2 + ismc_3 + ismc_4 + ismc_joint
  "
  
  m2 <- sem(model=ismc_measurement_2,
            data=ismc_chr1,
            estimator="MLR",
            optim.method="nlminb",
            std.ov=T,
            std.lv=T)
  
  lavTestLRT(m1, m2)
  # p = 0.5689
  
  ismc_chr1$ismc_latent <- as.numeric(predict(m1))
  
  q <- ggpairs(ismc_chr1[,4:ncol(ismc_chr1)])
  save_plot("sem/ismc_mut_Mbuti_chr1_1Mb.pdf", q, base_height=7, base_width=9)
  
  # preparing to fit measurement model with other mutation maps
  df_chr1_1Mb <- left_join(filter(df_1Mb, chrom==1),
                           ismc_chr1,
                           by=c("chrom", "chromStart", "chromEnd"))
  
  mut_chr1_1Mb <-  filter(df_chr1_1Mb, coverage >= 0.25) 
  
  mut_chr1_1Mb_log <- cbind.data.frame(mut_chr1_1Mb[,1:3],
                                       apply(mut_chr1_1Mb[,4:(ncol(mut_chr1_1Mb)-1)],
                                             2, log))
  mut_chr1_1Mb_log$ismc_latent <- mut_chr1_1Mb$ismc_latent
  
  mut_measurement_1 <- "
  mu=~mean_mut_masked_roulette + mean_mut_masked_carlson + mean_mut_masked_gnomad + ismc_mean 
  "
  
  m3 <- sem(model=mut_measurement_1,
            data=mut_chr1_1Mb_log,
            estimator="MLR",
            optim.method="nlminb",
            std.ov=T,
            std.lv=T,
            missing="fiml")
  
  summary(m3)
  fitMeasures(m3)
  modificationIndices(m3)
  
  mut_measurement_2 <- "
  mu=~mean_mut_masked_roulette + mean_mut_masked_carlson + mean_mut_masked_gnomad + ismc_mean 
  mean_mut_masked_carlson ~~ ismc_mean
  " 
  
  m4 <- sem(model=mut_measurement_2,
            data=mut_chr1_1Mb_log,
            estimator="MLR",
            optim.method="nlminb",
            std.ov=T,
            std.lv=T,
            missing="fiml")
  
  summary(m4)
  fitmeasures(m4)
  # p=0.162
  inspect(m4, "r2")
  parameterEstimates(m4)
  lavInspect(m4, what="free")
  vcov(m4)
  
  mut_chr1_1Mb_log$mut_latent <- as.numeric(predict(m4))
  
  u <- mut_chr1_1Mb_log %>% dplyr::select(., mean_mut_masked_roulette, 
                                             mean_mut_masked_carlson, 
                                             mean_mut_masked_gnomad,
                                             ismc_mean,
                                          mut_latent)
  
  names(u) <- c("roulette", "carlson", "gnomad", "ismc", "latent")
  p <- ggpairs(u)
  
  save_plot("sem/mut_measurement_chr1_1Mb.pdf", p, base_height=7, base_width=9)
}

################# all chr models #################
# (3 indicators => cannot test model fit, we want to extract the latent)

setDT(mut_maps_10kb)
mut_10kb_log <- mut_maps_10kb[,1:3]
mut_10kb_log[, roulette := log(mut_maps_10kb$roulette)]
mut_10kb_log[, carlson := log(mut_maps_10kb$carlson)]
mut_10kb_log[, gnomad := log(mut_maps_10kb$gnomad)]

setDT(mut_maps_100kb)
mut_100kb_log <- mut_maps_100kb[,1:3]
mut_100kb_log[, roulette := log(mut_maps_100kb$roulette)]
mut_100kb_log[, carlson := log(mut_maps_100kb$carlson)]
mut_100kb_log[, gnomad := log(mut_maps_100kb$gnomad)]

setDT(mut_maps_1Mb)
mut_1Mb_log <- mut_maps_1Mb[,1:3]
mut_1Mb_log[, roulette := log(mut_maps_1Mb$roulette)]
mut_1Mb_log[, carlson := log(mut_maps_1Mb$carlson)]
mut_1Mb_log[, gnomad := log(mut_maps_1Mb$gnomad)]

mut_measurement <- "
  mu =~ a*roulette + a*carlson + a*gnomad
  mu ~~ 1*mu # fix latent variance for identification
"

vars <- lavNames(mut_measurement, type = "ov")
print("Measurement model mutation...")

mut_10kb_log <- mut_10kb_log[rowSums(!is.na(mut_10kb_log[, ..vars])) > 0, ]
m5 <- sem(model=mut_measurement,
          data=mut_10kb_log,
          estimator="MLR",
          optim.method="nlminb",
          std.ov=T,
          missing="fiml")
  
mut_10kb_log$mut_latent_log <- lavPredict(m5)[,1]   
fwrite(dplyr::select(mut_10kb_log, c(chrom, chromStart, chromEnd, mut_latent_log)), "maps/mut_latent_10kb.csv.gz")
  
mut_100kb_log <- mut_100kb_log[rowSums(!is.na(mut_100kb_log[, ..vars])) > 0, ]
m5 <- sem(model=mut_measurement,
          data=mut_100kb_log,
          estimator="MLR",
          optim.method="nlminb",
          std.ov=T,
          missing="fiml")

mut_100kb_log$mut_latent_log <- lavPredict(m5)[,1]   
fwrite(dplyr::select(mut_100kb_log, c(chrom, chromStart, chromEnd, mut_latent_log)), "maps/mut_latent_100kb.csv.gz")
  
mut_1Mb_log <- mut_1Mb_log[rowSums(!is.na(mut_1Mb_log[, ..vars])) > 0, ]
m5 <- sem(model=mut_measurement,
          data=mut_1Mb_log,
          estimator="MLR",
          optim.method="nlminb",
          std.ov=T,
          missing="fiml")

mut_1Mb_log$mut_latent_log <- lavPredict(m5)[,1]  
fwrite(dplyr::select(mut_1Mb_log, c(chrom, chromStart, chromEnd, mut_latent_log)), "maps/mut_latent_1Mb.csv.gz")

# u <- mut_1Mb_log %>% dplyr::select(., mean_mut_masked_roulette, 
#                                       mean_mut_masked_carlson, 
#                                       mean_mut_masked_gnomad,
#                                       mut_latent_log)
# 
# names(u) <- c("roulette", "carlson", "gnomad", "latent")
# p <- ggpairs(u)

#################
#
# Recombination 
#
#################

################# loading data ################# 
# pyrho YRI
pyrho <- vector("list", length=num_chrom)
for(c in 1:num_chrom) {
  pyrho[[c]] <- fread(paste0("~/Data/bgs_lmr/human_data/data/",
                            "recombination_tables/pyrho/",
                            "pyrho_YRI_tbl_chr", c, ".csv.gz"))
}
pyrho_YRI_1kb <- data.table::rbindlist(pyrho)

# pyrho CEU
pyrho <- vector("list", length=num_chrom)
for(c in 1:num_chrom) {
  pyrho[[c]] <- fread(paste0("~/Data/bgs_lmr/human_data/data/",
                            "recombination_tables/pyrho/",
                            "pyrho_CEU_tbl_chr", c, ".csv.gz"))
}
pyrho_CEU_1kb <- data.table::rbindlist(pyrho)

# decode
decode <- vector("list", length=num_chrom)
for(c in 1:num_chrom) {
  decode[[c]] <- fread(paste0("~/Data/bgs_lmr/human_data/data/",
                             "recombination_tables/deCODE/",
                             "deCODE_tbl_chr", c, ".csv.gz"))
}
decode_1kb <- data.table::rbindlist(decode)

pyrho_CEU_1kb$rate[pyrho_CEU_1kb$num_sites_map==0] <- NA
pyrho_YRI_1kb$rate[pyrho_YRI_1kb$num_sites_map==0] <- NA
decode_1kb$rate[decode_1kb$num_sites_map==0] <- NA

names(pyrho_CEU_1kb)[5] <- "pyrho_CEU"
names(pyrho_YRI_1kb)[5] <- "pyrho_YRI"
names(decode_1kb)[5] <- "deCODE"

pyrho_CEU_gr <- makeGRangesFromDataFrame(dplyr::select(pyrho_CEU_1kb, -num_sites_map), keep.extra.columns=T)
pyrho_YRI_gr <- makeGRangesFromDataFrame(dplyr::select(pyrho_YRI_1kb, -num_sites_map), keep.extra.columns=T)
decode_gr <- makeGRangesFromDataFrame(dplyr::select(decode_1kb, -num_sites_map), keep.extra.columns=T)

# transform column names before merging
dt_CEU <- as.data.table(pyrho_CEU_gr)[, pyrho_CEU := mcols(pyrho_CEU_gr)$pyrho_CEU][, .(seqnames, start, end, pyrho_CEU)]
dt_YRI <- as.data.table(pyrho_YRI_gr)[, pyrho_YRI := mcols(pyrho_YRI_gr)$pyrho_YRI][, .(seqnames, start, end, pyrho_YRI)]
dt_de <- as.data.table(decode_gr)[, deCODE := mcols(decode_gr)$deCODE][, .(seqnames, start, end, deCODE)]

rec_maps_1kb <- Reduce(function(a, b) merge(a, b, by = c("seqnames","start","end"), all = TRUE),
                       list(dt_CEU, dt_YRI, dt_de))

names(rec_maps_1kb)[1:3] <- c("chrom", "chromStart", "chromEnd")

rec_maps_1kb[, chrom := as.integer(sub("^chr", "", chrom))]

rec_maps_1kb <- rec_maps_1kb %>% group_by(chrom) %>% 
  mutate(bin_10kb=(1:n() - 1) %/% 10,
         bin_100kb=(1:n() - 1) %/% 100,
         bin_1Mb=(1:n() - 1) %/% 1000)

rec_maps_10kb <- rec_maps_1kb %>%
  group_by(chrom, bin_10kb) %>% 
  mutate(mean_pyrho_CEU=mean(pyrho_CEU),
         mean_pyrho_YRI=mean(pyrho_YRI),
         mean_deCODE=mean(deCODE),
         chromStart=min(chromStart),
         chromEnd=max(chromEnd)) %>%
  distinct(bin_10kb, .keep_all=T) %>%
  ungroup() %>%
  dplyr::select(., c(chrom, chromStart, chromEnd, 
                     mean_pyrho_CEU,
                     mean_pyrho_YRI,
                     mean_deCODE))

rec_maps_100kb <- rec_maps_1kb %>%
  group_by(chrom, bin_100kb) %>% 
  mutate(mean_pyrho_CEU=mean(pyrho_CEU),
         mean_pyrho_YRI=mean(pyrho_YRI),
         mean_deCODE=mean(deCODE),
         chromStart=min(chromStart),
         chromEnd=max(chromEnd)) %>%
  distinct(bin_100kb, .keep_all=T) %>%
  ungroup() %>%
  dplyr::select(., c(chrom, chromStart, chromEnd, 
                     mean_pyrho_CEU,
                     mean_pyrho_YRI,
                     mean_deCODE))

rec_maps_1Mb <- rec_maps_1kb %>%
  group_by(chrom, bin_1Mb) %>% 
  mutate(mean_pyrho_CEU=mean(pyrho_CEU),
         mean_pyrho_YRI=mean(pyrho_YRI),
         mean_deCODE=mean(deCODE),
         chromStart=min(chromStart),
         chromEnd=max(chromEnd)) %>%
  distinct(bin_1Mb, .keep_all=T) %>%
  ungroup() %>%
  dplyr::select(., c(chrom, chromStart, chromEnd, 
                     mean_pyrho_CEU,
                     mean_pyrho_YRI,
                     mean_deCODE))

################# testing measurement model on chr1, 1 Mb scale ################# 
if(interactive()) {
  # ismc chr 1, from the book chapter Barroso & Dutheil 2026
  ismc_chr1 <- fread("/media/gvbarroso/extradrive1/iSMC/iSMC_Mbuti/Mbuti_chr1.rho.1Mb.bedgraph")
  names(ismc_chr1)[4:9] <- c(paste("ismc", 1:4, sep="_"), "ismc_mean", "ismc_joint")
  
  # w/ rec we must filter by hand since there's no "coverage" info in the tables
  ismc_missing <-  fread("/media/gvbarroso/extradrive1//iSMC/iSMC_Mbuti/Mbuti_chr1.missing.prop.1Mb.bedgraph")
  
  ismc_chr1$missing_prop <- ismc_missing$`LP6005441-DNA_A08`
  ismc_chr1$ismc_1[ismc_chr1$missing_prop > 0.75] <- NA
  ismc_chr1$ismc_2[ismc_chr1$missing_prop > 0.75] <- NA
  ismc_chr1$ismc_3[ismc_chr1$missing_prop > 0.75] <- NA
  ismc_chr1$ismc_4[ismc_chr1$missing_prop > 0.75] <- NA
  ismc_chr1$ismc_mean[ismc_chr1$missing_prop > 0.75] <- NA
  ismc_chr1$ismc_joint[ismc_chr1$missing_prop > 0.75] <- NA
  
  ismc_measurement_1 <- "
  ismc=~ismc_1 + ismc_2 + ismc_3 + ismc_4
  "
  
  m1 <- sem(model=ismc_measurement_1,
            data=ismc_chr1,
            estimator="MLR",
            optim.method="nlminb",
            std.ov=T,
            std.lv=T,
            missing="fiml")
  
  summary(m1, fit.measures = TRUE)
  inspect(m1, "r2")
  parameterEstimates(m1)
  fitMeasures(m1)
  
  ismc_measurement_2 <- "
  ismc=~ismc_1 + ismc_2 + ismc_3 + ismc_4 + ismc_joint
  "
  
  m2 <- sem(model=ismc_measurement_2,
            data=ismc_chr1,
            estimator="MLR",
            optim.method="nlminb",
            std.ov=T,
            std.lv=T,
            missing="fiml")
  
  lavTestLRT(m1, m2)
  
  ismc_chr1$ismc_latent <- as.numeric(predict(m1))
  
  q <- ggpairs(ismc_chr1[,4:ncol(ismc_chr1)])
  save_plot("real_data_maps/ismc_rec_pairs_chr1_1Mb.svg", q, base_height=10, base_width=12)
  
  df_chr1_1Mb <- left_join(filter(rec_maps_1Mb, chrom=="chr1"),
                           ismc_chr1,
                           by=c("chrom", "chromStart", "chromEnd"))
  
  p1 <- ggpairs(dplyr::select(df_chr1_1Mb, c(mean_pyrho_CEU,
                                             mean_pyrho_YRI,
                                             mean_deCODE,
                                             ismc_1,
                                             ismc_2,
                                             ismc_3,
                                             ismc_4,
                                             ismc_mean,
                                             ismc_latent)))
  
  save_plot("real_data_maps/rec_pairs_chr1_1Mb.svg", p1, base_width=12, base_height=10)
  
  rec_chr1_log <- cbind.data.frame(df_chr1_1Mb[,1:3],
                                   apply(df_chr1_1Mb[,4:ncol(df_chr1_1Mb)],
                                         2, log))
  rec_chr1_log$ismc_latent <- df_chr1_1Mb$ismc_latent
  
  rec_measurement_1 <- "
  rec=~mean_pyrho_CEU + mean_pyrho_YRI + mean_deCODE + ismc_latent 
  "
  
  m3 <- sem(model=rec_measurement_1,
            data=rec_chr1_log,
            estimator="MLR",
            optim.method="nlminb",
            std.ov=T,
            std.lv=T,
            missing="fiml")
  
  summary(m3)
  inspect(m3, "r2")
  parameterEstimates(m4)
  lavInspect(m4, what="free")
  vcov(m4)
  
  m4 <- sem(model=rec_measurement_1,
            data=na.omit(rec_chr1_log),
            estimator="MLR",
            optim.method="nlminb",
            std.ov=T,
            std.lv=T,
            estimator="GLS")
  
  summary(m4)
  inspect(m4, "r2")
  parameterEstimates(m4)
  lavInspect(m4, what="free")
  vcov(m4)
  
  rec_chr1_log$rec_latent <- as.numeric(predict(m3))
}

################# all chr models ################# 
# cannot test model fit, we want to extract the latent

setDT(rec_maps_10kb)
rec_10kb_log <- rec_maps_10kb[,1:3]
rec_10kb_log[, mean_pyrho_CEU := log1p(rec_maps_10kb$mean_pyrho_CEU)]
rec_10kb_log[, mean_pyrho_YRI := log1p(rec_maps_10kb$mean_pyrho_YRI)]
rec_10kb_log[, mean_deCODE := log1p(rec_maps_10kb$mean_deCODE)]

setDT(rec_maps_100kb)
rec_100kb_log <- rec_maps_100kb[,1:3]
rec_100kb_log[, mean_pyrho_CEU := log1p(rec_maps_100kb$mean_pyrho_CEU)]
rec_100kb_log[, mean_pyrho_YRI := log1p(rec_maps_100kb$mean_pyrho_YRI)]
rec_100kb_log[, mean_deCODE := log1p(rec_maps_100kb$mean_deCODE)]

setDT(rec_maps_1Mb)
rec_1Mb_log <- rec_maps_1Mb[,1:3]
rec_1Mb_log[, mean_pyrho_CEU := log1p(rec_maps_1Mb$mean_pyrho_CEU)]
rec_1Mb_log[, mean_pyrho_YRI := log1p(rec_maps_1Mb$mean_pyrho_YRI)]
rec_1Mb_log[, mean_deCODE := log1p(rec_maps_1Mb$mean_deCODE)]

rec_measurement <- "
  rec =~ a*mean_pyrho_CEU +
         a*mean_pyrho_YRI +
         a*mean_deCODE

  rec ~~ 1*rec # fix latent variance for identification
"

vars <- lavNames(rec_measurement, type = "ov")

rec_10kb_log[, (vars) := lapply(.SD, scale), .SDcols = vars]
rec_100kb_log[, (vars) := lapply(.SD, scale), .SDcols = vars]
rec_1Mb_log[, (vars) := lapply(.SD, scale), .SDcols = vars]

print("Measurement model recombination...")

rec_10kb_log <- rec_10kb_log[rowSums(!is.na(rec_10kb_log[, ..vars])) > 0]
m5 <- sem(model=rec_measurement,
          data=rec_10kb_log,
          estimator="MLR",
          optim.method="nlminb",
          std.ov=T,
          missing="fiml")

rec_10kb_log$rec_latent_log <- lavPredict(m5)[,1]   
df <- dplyr::select(rec_10kb_log, c(chrom, chromStart, chromEnd, rec_latent_log)) %>% setDT()
fwrite(df, paste0("maps/rec_latent_10kb.csv.gz"))

rec_100kb_log <- rec_100kb_log[rowSums(!is.na(rec_100kb_log[, ..vars])) > 0]
m5 <- sem(model=rec_measurement,
          data=rec_100kb_log,
          estimator="MLR",
          optim.method="nlminb",
          std.ov=T,
          missing="fiml")

rec_100kb_log$rec_latent_log <- lavPredict(m5)[,1]   
df <- dplyr::select(rec_100kb_log, c(chrom, chromStart, chromEnd, rec_latent_log)) %>% setDT()
fwrite(df, paste0("maps/rec_latent_100kb.csv.gz"))

rec_1Mb_log <- rec_1Mb_log[rowSums(!is.na(rec_1Mb_log[, ..vars])) > 0]
m5 <- sem(model=rec_measurement,
          data=rec_1Mb_log,
          estimator="MLR",
          optim.method="nlminb",
          std.ov=T,
          missing="fiml")

rec_1Mb_log$rec_latent_log <- lavPredict(m5)[,1]   
df <- dplyr::select(rec_1Mb_log, c(chrom, chromStart, chromEnd, rec_latent_log)) %>% setDT()
fwrite(df, paste0("maps/rec_latent_1Mb.csv.gz"))

#################
#
# B-values (all highly correlated)
#
#################

################# all chr models ################# 
# B maps from each mutation map 
# B-maps from mut maps are highly correlated -> we use PCA instead of SEM for numerical stability

for(m in constrained_models) {
  
  print(paste("Measurement model B-values:", m))
  
  tmp_10kb <- filter(human_maps_10kb, constrained_model==m) %>% 
    dplyr::select(., c(chrom, chromStart, chromEnd, constrained_model, B, mut_map)) %>%
    pivot_wider(., values_from=B, names_from=mut_map) 
  
  #b <- na.omit(tmp_10kb) %>% dplyr::select(., carlson, gnomad, roulette) %>% cor() # very high
  
  setDT(tmp_10kb)
  tmp_10kb_log <- tmp_10kb[,1:4]
  tmp_10kb_log[, carlson := log(tmp_10kb$carlson)]
  tmp_10kb_log[, gnomad := log(tmp_10kb$gnomad)]
  tmp_10kb_log[, roulette := log(tmp_10kb$roulette)]
  
  tmp_100kb <- filter(human_maps_100kb, constrained_model==m) %>%
    dplyr::select(., c(chrom, chromStart, chromEnd, constrained_model, B, mut_map)) %>%
    pivot_wider(., values_from=B, names_from=mut_map) 
  
  setDT(tmp_100kb)
  tmp_100kb_log <- tmp_100kb[,1:4]
  tmp_100kb_log[, carlson := log(tmp_100kb$carlson)]
  tmp_100kb_log[, gnomad := log(tmp_100kb$gnomad)]
  tmp_100kb_log[, roulette := log(tmp_100kb$roulette)]
  
  tmp_1Mb <- filter(human_maps_1Mb, constrained_model==m) %>%
    dplyr::select(., c(chrom, chromStart, chromEnd, constrained_model, B, mut_map)) %>%
    pivot_wider(., values_from=B, names_from=mut_map) 
  
  setDT(tmp_1Mb)
  tmp_1Mb_log <- tmp_1Mb[,1:4]
  tmp_1Mb_log[, carlson := log(tmp_1Mb$carlson)]
  tmp_1Mb_log[, gnomad := log(tmp_1Mb$gnomad)]
  tmp_1Mb_log[, roulette := log(tmp_1Mb$roulette)]
  
  vars <- c("carlson", "gnomad", "roulette")
  
  # PCA on the covariance matrix not to miss windows due to lower Carlson coverage
  S <- cov(tmp_10kb_log[, ..vars], use = "pairwise.complete.obs")
  pc <- eigen(S)
  loadings <- pc$vectors[, 1]
  Z <- scale(tmp_10kb_log[, ..vars])
  
  pc1_raw <- as.vector(Z %*% loadings)
  ref <- tmp_10kb_log$roulette
  sgn <- sign(cor(pc1_raw, ref, use = "complete.obs"))
  tmp_10kb_log$B_latent_log <- as.numeric(pc1_raw) * sgn
  
  print(cor(na.omit(tmp_10kb_log[, 5:8])))
  
  fwrite(dplyr::select(tmp_10kb_log, c(chrom, chromStart, chromEnd, constrained_model, B_latent_log)),
         paste0("maps/B_latent_", m, "_10kb.csv.gz"))
  
  # PCA on the covariance matrix not to miss windows due to lower Carlson coverage
  S <- cov(tmp_100kb_log[, ..vars], use = "pairwise.complete.obs")
  pc <- eigen(S)
  loadings <- pc$vectors[, 1]
  Z <- scale(tmp_100kb_log[, ..vars])
  
  pc1_raw <- as.vector(Z %*% loadings)
  ref <- tmp_100kb_log$roulette
  sgn <- sign(cor(pc1_raw, ref, use = "complete.obs"))
  tmp_100kb_log$B_latent_log <- as.numeric(pc1_raw) * sgn
  
  print(cor(na.omit(tmp_100kb_log[, 5:8])))
  
  fwrite(dplyr::select(tmp_100kb_log, c(chrom, chromStart, chromEnd, constrained_model, B_latent_log)),
         paste0("maps/B_latent_", m, "_100kb.csv.gz"))
  
  # PCA on the covariance matrix not to miss windows due to lower Carlson coverage
  S <- cov(tmp_1Mb_log[, ..vars], use = "pairwise.complete.obs")
  pc <- eigen(S)
  loadings <- pc$vectors[, 1]
  Z <- scale(tmp_1Mb_log[, ..vars])
  
  pc1_raw <- as.vector(Z %*% loadings)
  ref <- tmp_1Mb_log$roulette
  sgn <- sign(cor(pc1_raw, ref, use = "complete.obs"))
  tmp_1Mb_log$B_latent_log <- as.numeric(pc1_raw) * sgn
  
  print(cor(na.omit(tmp_1Mb_log[, 5:8])))
  
  fwrite(dplyr::select(tmp_1Mb_log, c(chrom, chromStart, chromEnd, constrained_model, B_latent_log)),
         paste0("maps/B_latent_", m, "_1Mb.csv.gz"))
}
