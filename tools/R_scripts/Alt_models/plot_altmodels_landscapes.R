library(ggplotify)
library(GGally)
library(tidyverse)
library(data.table)
library(cowplot)
library(gridExtra)

# Load 1 Mb
cousins_1Mb <- fread("~/Data/bgs_lmr/PNAS_reviews/cousins_demo/maps/maps_1Mb.csv.gz")
singleNe_1Mb <- fread("~/Data/bgs_lmr/human_data/current_models/equilibrium/maps/maps_1Mb.csv.gz")
decode_1Mb <- fread("~/Data/bgs_lmr/PNAS_reviews/deCODE/maps/maps_1Mb.csv.gz")

cousins_1Mb$model <- "Cousins"
singleNe_1Mb$model <- "Single-Ne"
decode_1Mb$model <- "deCODE"

maps_1Mb <- rbind.data.frame(cousins_1Mb, singleNe_1Mb, decode_1Mb) %>% setDT()

# Load 100 kb
cousins_100kb <- fread("~/Data/bgs_lmr/PNAS_reviews/cousins_demo/maps/maps_100kb.csv.gz")
singleNe_100kb <- fread("~/Data/bgs_lmr/human_data/current_models/equilibrium/maps/maps_100kb.csv.gz")
decode_100kb <- fread("~/Data/bgs_lmr/PNAS_reviews/deCODE/maps/maps_100kb.csv.gz")

cousins_100kb$model <- "Cousins"
singleNe_100kb$model <- "Single-Ne"
decode_100kb$model <- "deCODE"

maps_100kb <- rbind.data.frame(cousins_100kb, singleNe_100kb, decode_100kb) %>% setDT()

c <- 8
dat_1Mb <- maps_1Mb %>%
  filter(mut_map=="roulette", 
         constrained_model %in% constrained_models,
         chrom == c, 
         coverage >= 0.25) %>%
  arrange(chromStart) %>%
  mutate(gap = chromStart - lag(chromStart, default = chromStart[1]))

dat_1Mb[1, gap := 0]

dat_1Mb[, contiguous := gap == 1e+6 | gap == 0]
dat_1Mb[, run_id := rleid(contiguous)]
dat_1Mb <- dat_1Mb[contiguous == TRUE]

dat_1Mb[, constrained_model := as.character(constrained_model)]
dat_1Mb[, mut_map := as.character(mut_map)]

dat_1Mb$constrained_model <- factor(dat_1Mb$constrained_model, 
                                    levels = names(elem_labels), 
                                    labels = unname(elem_labels)) 

dat_1Mb[, plot_group := paste0("r", run_id, "_cm_", model)]

p_1Mb <- ggplot(dat_1Mb, aes(x = chromStart / 1e6)) + theme_classic() +
  geom_path(data = dat_1Mb, 
            aes(y = avg_pi,
                color = "obs_pi", 
                group = paste0("r", run_id)), linewidth = 0.75) +
  geom_point(data = dat_1Mb,
             aes(y = avg_pi,
                 color = "obs_pi", 
                 group = paste0("r", run_id)), size = 1.5) +
  geom_path(data = dat_1Mb,
            aes(y = exp_pi,
                color = model,
                group = plot_group), linewidth = 0.75) +
  geom_point(data = dat_1Mb, shape=16,
             aes(y = exp_pi, 
                 color = model, 
                 group = plot_group), size = 1.5) +
  facet_wrap(~constrained_model, nrow=2) +
  scale_x_continuous(breaks = seq(from=0, to=140, by=10), expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values = c("obs_pi" = "black",
                                "deCODE" = "cyan3",
                                "Single-Ne" = "brown1",
                                "Cousins" = "purple3"),
                     breaks = c("obs_pi", "deCODE", "Single-Ne", "Cousins"),
                     labels = c("Observed",
                                "deCODE",
                                "Single-Ne", 
                                "Cousins"),
                     name=NULL) +
  #annotate("text", label = "Chromosome 8", size = 5,
  #         x = mean(range(dat_1Mb$chromStart / 1e6)),
  #         y = 0.95 * max(range(c(dat_1Mb$avg_pi, dat_1Mb$exp_pi)))) +
  labs(title = NULL, x = "Position chromosome 8 (1-Mb scale)", y = expression(pi)) +
  theme(axis.title.x = element_blank(),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 22),
        #legend.position = c(0.85, 0.85),
        legend.position = "bottom",
        legend.text=element_text(size = 16))
save_plot("~/Data/bgs_lmr/PNAS_reviews/alt_models_landscapes_chr8.pdf", p_1Mb, base_height=8, base_width=15)

c <- 22
dat_100kb <- maps_100kb %>%
  filter(mut_map=="roulette",
         constrained_model %in% constrained_models,
         chrom == c, 
         coverage >= 0.5) %>%
  arrange(chromStart) %>%
  mutate(gap = chromStart - lag(chromStart, default = chromStart[1]))

dat_100kb$constrained_model <- factor(dat_100kb$constrained_model, 
                                      levels = names(elem_labels), 
                                      labels = unname(elem_labels)) 
dat_100kb[1, gap := 0]

dat_100kb[, contiguous := gap == 1e+5 | gap == 0]
dat_100kb[, run_id := rleid(contiguous)]
dat_100kb <- dat_100kb[contiguous == TRUE]

dat_100kb[, constrained_model := as.character(constrained_model)]
dat_100kb[, mut_map := as.character(mut_map)]

dat_100kb[, plot_group := paste0("r", run_id, "_cm_", model)]

p_100kb <- ggplot(dat_100kb, aes(x = chromStart / 1e6)) + theme_classic() +
  geom_path(data = dat_100kb, 
            aes(y = avg_pi,
                color = "obs_pi", 
                group = paste0("r", run_id)), linewidth = 0.75) +
  geom_point(data = dat_100kb,
             aes(y = avg_pi,
                 color = "obs_pi", 
                 group = paste0("r", run_id)), size = 1.5) +
  geom_path(data = dat_100kb,
            aes(y = exp_pi, 
                color = model,
                group = plot_group), linewidth = 0.75) +
  geom_point(data = dat_100kb,
             shape=16,
             aes(y = exp_pi,
                 color = model, 
                 group = plot_group), size = 1.5) +
  #annotate("text", label = "Chromosome 22", size = 5,
  #         x = mean(range(dat_100kb$chromStart / 1e6)),
  #         y = 0.95 * max(range(c(dat_100kb$avg_pi, dat_100kb$exp_pi)))) +
  facet_wrap(~constrained_model, nrow=2) +
  scale_x_continuous(breaks = seq(from=20, to=50, by=5), expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values = c("obs_pi" = "black",
                                "deCODE" = "cyan3",
                                "Single-Ne" = "brown1",
                                "Cousins" = "purple3"),
                     breaks = c("obs_pi", "deCODE", "Single-Ne", "Cousins"),
                     labels = c("Observed",
                                "deCODE",
                                "Single-Ne", 
                                "Cousins"),
                     name=NULL) +
  labs(title = NULL, x = "Position chromosome 22 (100-kb scale)", y = expression(pi)) +
  theme(strip.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.position="bottom")
save_plot("~/Data/bgs_lmr/PNAS_reviews/alt_models_landscapes_chr22.pdf", p_100kb, base_height=8, base_width=15)

# looking at B-map correlations between alt models
maps_100kb[, group:=paste(chrom, chromStart, chromEnd, sep="-")]
dt_100kb <- maps_100kb %>%
  dplyr::filter(mut_map=="roulette",
         constrained_model==constrained_models[1],
         coverage >= 0.5) %>% dplyr::select(., c("B", "avg_rec", "model", "group")) %>% setDT()

wide <- dcast(dt_100kb, group ~ model, value.var = c("B", "avg_rec"))

all_windows <- wide %>% dplyr::select(., starts_with("B_")) %>% ggpairs()
deCODE_gtX <- wide %>% dplyr::filter(., avg_rec_deCODE > 1e-9) %>%
  dplyr::select(., starts_with("B_")) %>% ggpairs()

save_plot("~/Data/bgs_lmr/PNAS_reviews/B_alt_models_all.pdf", all_windows, base_height=7, base_width=10)
save_plot("~/Data/bgs_lmr/PNAS_reviews/B_alt_models_deCODEgt1e-9.pdf", deCODE_gtX, base_height=7, base_width=10)


a <- ggplot(wide, aes(x = avg_rec_Cousins, y = avg_rec_deCODE)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_smooth(method="lm", se=FALSE, color="red") +
  stat_cor(method = "pearson") +
  labs(x="pyrho", y="deCODE") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))

b <- ggplot(dplyr::filter(wide, avg_rec_deCODE > 1e-9), aes(x = avg_rec_Cousins, y = avg_rec_deCODE)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_smooth(method="lm", se=FALSE, color="red") +
  stat_cor(method = "pearson") +
  labs(x="pyrho", y=NULL) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.text.y=element_blank(),
        axis.title=element_text(size=16))

q <- plot_grid(a, b, labels="AUTO", nrow=1)
save_plot("~/Data/bgs_lmr/PNAS_reviews/rec_maps.pdf", q, base_height=5, base_width=10)
