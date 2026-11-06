
library(data.table)
library(tidyverse)
library(cowplot)

#################
# 
# comparison between phastCons and regulatory models
#
#################

file_list <- list.files("data/pi_tables/roulette/")
pis_1kb <- data.table::rbindlist(lapply(paste0("data/pi_tables/roulette/", file_list), fread))
pis_1kb[, mut_map :="roulette"]
pis_1kb[, pos := chromStart + 500]

file_list <- list.files("models/B_maps_1kb/equilibrium_granular_Ne/split_cds_phastcons/roulette/", full.names=T)
bs_1kb_phast <- data.table::rbindlist(lapply(file_list, fread))
bs_1kb_phast[, mut_map :="roulette"]
bs_1kb_phast[, element := "CDS + phastCons"]

file_list <- list.files("models/B_maps_1kb/equilibrium_granular_Ne/split_cds_regulatory/roulette/", , full.names=T)
bs_1kb_reg <- data.table::rbindlist(lapply(file_list, fread))
bs_1kb_reg[, mut_map :="roulette"]
bs_1kb_reg[, element := "CDS + regulatory"]

# pause to find intersection
inter <- merge(bs_1kb_reg, bs_1kb_phast, by=c("chrom", "pos"))
inter <- merge(inter, pis_1kb, by=c("chrom", "pos"))
inter[B.x > 0.9 & B.y > 0.9,]
setnames(inter, old=c("B.x", "B.y"), new=c("B_reg", "B_phast"))
inter <- dplyr::select(inter, c(chrom, pos, B_reg, B_phast, num_sites))
fwrite(inter, "B_intersection.csv.gz")

tbl_1kb_reg <- merge(pis_1kb, bs_1kb_reg) %>% filter(., num_sites >= 500) %>% setDT()
tbl_1kb_phast <- merge(pis_1kb, bs_1kb_phast) %>% filter(., num_sites >= 500) %>% setDT()

# assigning to 1-percentile bins of B
tbl_1kb_reg[, B_bin := cut(B, breaks = quantile(B, probs = seq(0, 1, 0.01), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]
tbl_1kb_phast[, B_bin := cut(B, breaks = quantile(B, probs = seq(0, 1, 0.01), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)]

tbl_1kb_reg[, B_percentile := frank(B, ties.method = "average") / .N]
tbl_1kb_phast[, B_percentile := frank(B, ties.method = "average") / .N]

cor.test(tbl_1kb_reg$B_percentile, tbl_1kb_phast$B_percentile)

tbl_1kb <- rbind.data.frame(tbl_1kb_reg, tbl_1kb_phast) %>% as.data.frame()
tbl_1kb_fixed <- data.table::as.data.table(tbl_1kb)

hist <- ggplot(data=tbl_1kb_fixed, aes(x=B, fill=element)) +
  geom_histogram(alpha=0.75, binwidth=0.005, color="black", position="identity") +
  scale_fill_manual(values=c("firebrick", "seagreen4")) +
  labs(x="B", y="Density", fill=NULL) + theme_classic() +
  scale_x_continuous(breaks = seq(0, 1, length.out=5), expand = expansion(mult = 0.01)) +
  theme(panel.grid.minor=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        axis.text.y=element_text(hjust=1),
        legend.position = c(0.2, 0.75),
        legend.text=element_text(size = 16))
save_plot("plots/B_1kb_hist.pdf", hist, base_height=5, base_width=9)

dt_wide <- dcast(tbl_1kb_fixed[element %in% c("CDS + phastCons", "CDS + regulatory")],
                 chrom + pos + mut_map ~ element, value.var = "B")

sp <- ggplot(dt_wide[sample(.N, .N * 0.05)], aes(x = `CDS + phastCons`, y = `CDS + regulatory`)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_smooth(method="lm", se=FALSE, color="red") +
  stat_cor(method = "pearson") +
  labs(x="B (CDS + phastCons)", y="B (CDS + regulatory)") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))
save_plot("plots/B_scatter.pdf", sp, base_height=5, base_width=9)

cp <- plot_grid(hist, sp, nrow=1, align="h", labels=c("C", "D"))
save_plot("plots/B_reg_phast.pdf", cp, base_height=5, base_width=15)

tbl_1kb_reg[, B_1pB := mean(B), by=.(B_bin)] # mean B within each 1-percentile of B
tbl_1kb_reg[, pi_1pB := mean(avg_pi), by=.(B_bin)] # mean pi within each 1-percentile of B
tbl_1kb_phast[, B_1pB := mean(B), by=.(B_bin)] # mean B within each 1-percentile of B
tbl_1kb_phast[, pi_1pB := mean(avg_pi), by=.(B_bin)] # mean pi within each 1-percentile of B

clean_reg <- tbl_1kb_reg[!duplicated(B_bin)]
clean_phast <- tbl_1kb_phast[!duplicated(B_bin)]

plot(clean_reg$B_1pB, clean_reg$pi_1pB)
plot(clean_phast$B_1pB, clean_phast$pi_1pB)

dat <- inter[B_reg > 0.9 & B_phast > 0.9 & num_sites > 500,]
  
dat_m <- pivot_longer(dat, cols=starts_with("B"), names_to="model", values_to="B") %>% setDT()
dat_m$model <- factor(dat_m$model, levels=c("B_reg", "B_phast"), labels=c("Regulatory", "PhastCons"))

dat_m[, gap := pos - data.table::shift(pos), by=chrom]
dat_m[!duplicated(chrom), gap := 0]

dat_m[, contiguous := gap == 1e+3 | gap == 0]
dat_m[, run_id := rleid(contiguous)]
dat_m[, plot_group := paste0("r", run_id, "_cm_", model)]

bs_1kb_gw <- ggplot(data=filter(dat_m, chrom==21), aes(x = pos / 1e6)) + theme_classic() +
  geom_path(data=filter(dat_m, chrom==21), aes(y = B, color = model, group = plot_group), linewidth = 0.75) +
  geom_point(data=filter(dat_m, chrom==21), shape=16, aes(y = B, color = model, group = plot_group), size = 1.5) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_shape_manual(values=c(0, 1)) +
  labs(title = NULL, x = "Position (Mb)", y = "B") +
  theme(axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.position="bottom")
