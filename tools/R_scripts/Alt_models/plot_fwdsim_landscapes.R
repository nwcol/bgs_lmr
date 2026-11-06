library(ggplotify)
library(GGally)
library(tidyverse)
library(data.table)
library(cowplot)
library(gridExtra)


# Load 1 Mb
true_1Mb <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/100Mb/demog_predictions/maps/maps_1Mb.csv.gz")
fit_good_1Mb <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/100Mb/fit_Ne_predictions_good_DFE/maps/maps_1Mb.csv.gz")
fit_bad_1Mb <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/100Mb/fit_Ne_predictions_bad_DFE/maps/maps_1Mb.csv.gz")

true_1Mb$model <- "True model"
fit_good_1Mb$model <- "Eq. good-DFE"
fit_bad_1Mb$model <- "Eq. bad-DFE"

maps_1Mb <- rbind.data.frame(true_1Mb, fit_good_1Mb, fit_bad_1Mb) %>% setDT()

# Load 100 kb
true_100kb <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/10Mb/demog_predictions/maps/maps_100kb.csv.gz")
fit_good_100kb <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/10Mb/fit_Ne_predictions_good_DFE/maps/maps_100kb.csv.gz")
fit_bad_100kb <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/10Mb/fit_Ne_predictions_bad_DFE/maps/maps_100kb.csv.gz")

true_100kb$model <- "True model"
fit_good_100kb$model <- "Eq. good-DFE"
fit_bad_100kb$model <- "Eq. bad-DFE"

maps_100kb <- rbind.data.frame(true_100kb, fit_good_100kb, fit_bad_100kb) %>% setDT()

dat_1Mb <- filter(maps_1Mb, chrom==1)
p_1Mb <- ggplot(dat_1Mb, aes(x = chromStart / 1e6)) + theme_classic() +
  geom_path(data = dat_1Mb, aes(y = avg_pi, color = "obs_pi", group=model), linewidth = 0.75) +
  geom_point(data = dat_1Mb, aes(y = avg_pi, color = "obs_pi", group=model), size = 1.5) +
  geom_path(data = dat_1Mb, aes(y = exp_pi, color = model, group=model), linewidth = 0.75) +
  geom_point(data = dat_1Mb, shape=16, aes(y = exp_pi, color = model, group=model), size = 1.5) +
  scale_x_continuous(breaks = pretty_breaks(), expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values = c("obs_pi" = "black",
                                "Eq. good-DFE" = "seagreen4",
                                "Eq. bad-DFE" = "purple2",
                                "True model" = "brown1"),
                     breaks = c("obs_pi", "Eq. good-DFE", "Eq. bad-DFE", "True model"),
                     labels=c(expression(pi[obs]), "Eq. good-DFE", "Eq. bad-DFE", "True model"),
                     name=NULL) +
  labs(title = NULL, x = "Position (Mb)", y = expression(pi)) +
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        legend.position = "bottom", #c(0.85, 0.85),
        legend.text=element_text(size = 18))


dat_100kb <- filter(maps_100kb, chrom==1)
p_100kb <- ggplot(dat_100kb, aes(x = chromStart / 1e6)) + theme_classic() +
  geom_path(data = dat_100kb, aes(y = avg_pi, color = "obs_pi", group=model), linewidth = 0.75) +
  geom_point(data = dat_100kb, aes(y = avg_pi, color = "obs_pi", group=model), size = 1.5) +
  geom_path(data = dat_100kb, aes(y = exp_pi, color = model, group=model), linewidth = 0.75) +
  geom_point(data = dat_100kb, shape=16, aes(y = exp_pi, color = model, group=model), size = 1.5) +
  scale_x_continuous(breaks = pretty_breaks(), expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values = c("obs_pi" = "black",
                                "Eq. good-DFE" = "seagreen4",
                                "Eq. bad-DFE" = "purple2",
                                "True model" = "brown1"),
                     breaks = c("obs_pi", "Eq. good-DFE", "Eq. bad-DFE", "True model"),
                     labels=c(expression(pi[obs]), "Eq. good-DFE", "Eq. bad-DFE", "True model"),
                     name=NULL) +
  labs(title = NULL, x = NULL, y = expression(pi)) +
  theme(axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.position="none")

p <- plot_grid(p_100kb, p_1Mb, nrow=2, rel_heights=c(1, 1.15), labels="AUTO", align="v")
save_plot("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/fwdsim_pi.pdf", p, base_height=8, base_width=15)

# 100 Mb chrom
vars_100Mb <-  dplyr::select(maps_1Mb, "prop_del_sites", "avg_mut", "avg_rec", "B", "exp_pi", "avg_pi") %>% sample_frac(., 0.2)
names(vars_100Mb) <- c("E", "u", "r", "B", "pi_exp", "pi_obs")

vars_10Mb <- dplyr::select(maps_100kb, "prop_del_sites", "avg_mut", "avg_rec", "B", "exp_pi", "avg_pi") %>% sample_frac(., 0.1)
names(vars_10Mb) <- c("E", "u", "r", "B", "pi_exp", "pi_obs")

cors_100Mb_lands <- ggpairs(vars_100Mb)
cors_10Mb_lands <- ggpairs(vars_10Mb)

save_plot("100Mb/cors_100Mb_lands.pdf", cors_100Mb_lands, base_heigh=7, base_width=15)
save_plot("10Mb/cors_10Mb_lands.pdf", cors_10Mb_lands, base_heigh=7, base_width=15)





