library(data.table)
library(tidyverse)
library(cowplot)

setwd("~/Data/bgs_lmr/PNAS_reviews/deCODE/fit_Ne/")

Ne_opt <- vector("list", length=2 * 3 * 22)
for(j in 0:65) {
  opt <- fread(paste0("split_cds_phastcons/reports/Ne_j", j, ".txt"))
  Ne_opt[[j+1]] <- cbind.data.frame(opt, j + 1)
}

Ne_opt <- data.table::rbindlist(Ne_opt)
names(Ne_opt) <- c("Ne", "job")
models_phast <- fread("split_cds_phastcons/split_cds_variables.txt") 
names(models_phast) <- c("chrom", "unlinkedB", "map", "memory", "name")
models_phast$job <- 1:nrow(models_phast)

Ne_phast <- merge(models_phast, Ne_opt, by="job") %>% setDT()
Ne_phast[, job := NULL]
setcolorder(Ne_phast, c("chrom", "unlinkedB", "Ne", "map", "memory", "name"))

fwrite(Ne_phast, "decode_phastcons_variables.txt", col.names=F)

# regulatory
Ne_opt <- vector("list", length=2 * 3 * 22)
for(j in 0:65) {
  opt <- fread(paste0("split_cds_regulatory/reports/Ne_j", j, ".txt"))
  Ne_opt[[j+1]] <- cbind.data.frame(opt, j + 1)
}

Ne_opt <- data.table::rbindlist(Ne_opt)
names(Ne_opt) <- c("Ne", "job")
models_reg <- fread("split_cds_regulatory/split_cds_variables.txt")
names(models_reg) <- c("chrom", "unlinkedB", "map", "memory", "name")
models_reg$job <- 1:nrow(models_reg)

Ne_reg <- merge(models_reg, Ne_opt, by="job") %>% setDT()
Ne_reg[, job := NULL]
setcolorder(Ne_reg, c("chrom", "unlinkedB", "Ne", "map", "memory", "name"))

fwrite(Ne_reg, "decode_regulatory_variables.txt", col.names=F)

Ne_reg_tbl <- dplyr::select(Ne_reg, c("chrom", "Ne", "map"))
Ne_reg_tbl$constrained_model <- "split CDS+regulatory"

Ne_phast_tbl <- dplyr::select(Ne_phast, c("chrom", "Ne", "map"))
Ne_phast_tbl$constrained_model <- "split CDS+phastCons"

Ne_tbl <- rbind.data.frame(Ne_reg_tbl, Ne_phast_tbl) %>% setDT()
chr_lengths <- fread("~/Devel/bgs_lmr/data/h38.genome")
Ne_tbl <- merge(Ne_tbl, chr_lengths, by="chrom")

mean_lines <- Ne_tbl |> dplyr::group_by(map, constrained_model) |> dplyr::summarise(mean_Ne = mean(Ne))

p <- ggplot(data=Ne_tbl, aes(x=length/1e6, y=Ne, color=map)) +
  geom_text(aes(label=chrom), vjust=2, size=3, show.legend=F) +
  geom_hline(data=mean_lines, aes(yintercept=mean_Ne, color=map), linewidth=0.8, inherit.aes=F) +
  facet_wrap(~constrained_model, ncol=1) + theme_classic() + 
  scale_y_continuous(breaks=seq(18e3, 38e3, 2e3), limits=c(18e3, 38e3)) +
  scale_x_continuous(breaks=pretty_breaks()) +
  scale_color_manual(values=c("steelblue2", "orange2", "seagreen"), name="Mutation map") +
  labs(title=NULL, x="Length (Mb)", y="Estimated Ne") +
  theme(strip.text.y=element_text(size=16),  
        strip.text.x=element_text(size=16),
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
p

save_plot("~/Data/bgs_lmr/PNAS_reviews/deCODE/Ne_est.pdf", p, base_height=8, base_width=12)
