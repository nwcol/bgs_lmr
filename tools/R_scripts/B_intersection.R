

dat <- fread("B_intersection.csv.gz")

dat <- dat[B_reg > 0.9 & B_phast > 0.9 & num_sites > 500,]
  
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
