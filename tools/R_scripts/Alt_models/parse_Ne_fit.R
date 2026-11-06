library(data.table)
library(tidyverse)
library(cowplot)

setwd("~/Data/bgs_lmr/PNAS_reviews/")

files <- list.files(pattern="*log.txt")

Ne_opt <- vector("list", length=5 * 22)
x <- 1
for(t in c(5, 10, 20, 30, 40)) {
  for(m in c("carlson", "gnomad", "roulette")) {
    for(c in 1:22) {
      line <- readLines(paste0("cutoff_phastcons/fit_Ne/top", t, "/split_cds_phastcons_", m, "_chr", c, "_log.txt"))[1]
      Ne <- as.numeric(unlist(strsplit(unlist(strsplit(line, ","))[1], ": "))[2])
      
      Ne_opt[[x]] <- cbind.data.frame(Ne, m, t, c)
      x <- x + 1
    }
  }
}

Ne_opt <- data.table::rbindlist(Ne_opt)
names(Ne_opt) <- c("Ne", "map", "cut-off", "chrom")

fwrite(Ne_opt[`cut-off`==5,], "cutoff_phastcons/fit_Ne/top5/Ne_estimates.csv")
fwrite(Ne_opt[`cut-off`==10,], "cutoff_phastcons/fit_Ne/top10/Ne_estimates.csv")
fwrite(Ne_opt[`cut-off`==20,], "cutoff_phastcons/fit_Ne/top20/Ne_estimates.csv")
fwrite(Ne_opt[`cut-off`==30,], "cutoff_phastcons/fit_Ne/top30/Ne_estimates.csv")
fwrite(Ne_opt[`cut-off`==40,], "cutoff_phastcons/fit_Ne/top40/Ne_estimates.csv")

p <- ggplot(data=Ne_opt, aes(x=chrom, y=Ne, shape=as.factor(`cut-off`))) +
  geom_point(size=3) + theme_classic() +
  facet_grid(map~`cut-off`) +
  labs(title=NULL, y="Estimated Ne", x="Chromosome") +
  scale_shape_manual(values=0:4) +
  theme(strip.text=element_text(size=16),
        axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        legend.position="none")
save_plot("cutoff_phastcons/fit_Ne/CDS+phast_Ne_cut-offs.pdf", p, base_height=6, base_width=12)

# using Ne estimates to edit the variables files for submitting predic_B jobs
vars_top5 <- fread("cutoff_phastcons/top5/split_cds_phastcons_variables.txt")

uBGS_carlson <- filter(vars_top5, V3=="carlson")$V2
uBGS_gnomad <- filter(vars_top5, V3=="gnomad")$V2
uBGS_roulette <- uBGS_gnomad # proxy

df_carlson <- cbind.data.frame(chr, uBGS_carlson, "carlson", "32GB")
df_gnomad <- cbind.data.frame(chr, uBGS_gnomad, "gnomad", "32GB")
df_roulette <- cbind.data.frame(chr, uBGS_roulette, "roulette", "32GB")

names(df_carlson) <- paste0("V", 1:ncol(df_carlson))
names(df_gnomad) <- paste0("V", 1:ncol(df_carlson))
names(df_roulette) <- paste0("V", 1:ncol(df_carlson))

df_vars <- rbind.data.frame(df_carlson, df_gnomad, df_roulette) %>% setDT()
df_vars[, fnames := paste0("split_cds_phastcons_", V3, "_chr", V1, ".csv")]

Ne5 <- Ne_opt[`cut-off`==5,]
setnames(Ne5, old=c("chrom", "map"), new=c("V1", "V3"))
df_co5 <- left_join(df_vars, Ne5, by=c("V1", "V3"))
df_co5[, `cut-off` := NULL]
setcolorder(df_co5, c("V1", "V2", "Ne", "V3", "V4", "fnames"))
setorder(df_co5, "V1", "V3")
fwrite(df_co5, "cutoff_phastcons/split_cds_phastcons_top5_variables.txt", col.names=F)

Ne10 <- Ne_opt[`cut-off`==10,]
setnames(Ne10, old=c("chrom", "map"), new=c("V1", "V3"))
df_co10 <- left_join(df_vars, Ne10, by=c("V1", "V3"))
df_co10[, `cut-off` := NULL]
setcolorder(df_co10, c("V1", "V2", "Ne", "V3", "V4", "fnames"))
setorder(df_co10, "V1", "V3")
fwrite(df_co10, "cutoff_phastcons/split_cds_phastcons_top10_variables.txt", col.names=F)

Ne20 <- Ne_opt[`cut-off`==20,]
setnames(Ne20, old=c("chrom", "map"), new=c("V1", "V3"))
df_co20 <- left_join(df_vars, Ne20, by=c("V1", "V3"))
df_co20[, `cut-off` := NULL]
setcolorder(df_co20, c("V1", "V2", "Ne", "V3", "V4", "fnames"))
setorder(df_co20, "V1", "V3")
fwrite(df_co20, "cutoff_phastcons/split_cds_phastcons_top20_variables.txt", col.names=F)

Ne30 <- Ne_opt[`cut-off`==30,]
setnames(Ne30, old=c("chrom", "map"), new=c("V1", "V3"))
df_co30 <- left_join(df_vars, Ne30, by=c("V1", "V3"))
df_co30[, `cut-off` := NULL]
setcolorder(df_co30, c("V1", "V2", "Ne", "V3", "V4", "fnames"))
setorder(df_co30, "V1", "V3")
fwrite(df_co30, "cutoff_phastcons/split_cds_phastcons_top30_variables.txt", col.names=F)

Ne40 <- Ne_opt[`cut-off`==40,]
setnames(Ne40, old=c("chrom", "map"), new=c("V1", "V3"))
df_co40 <- left_join(df_vars, Ne40, by=c("V1", "V3"))
df_co40[, `cut-off` := NULL]
setcolorder(df_co40, c("V1", "V2", "Ne", "V3", "V4", "fnames"))
setorder(df_co40, "V1", "V3")
fwrite(df_co40, "cutoff_phastcons/split_cds_phastcons_top40_variables.txt", col.names=F)
