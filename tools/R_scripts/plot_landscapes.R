
#################
# 
# This script uses prepared tables to make some landscape plots
# NOTE: this script should be run from within each current_models dir
#
#################

source("~/Data/bgs_lmr/human_data/tools/HelperFunctions.R")

constrained_models <- c("split_cds_phastcons", "split_cds_regulatory")
mutation_maps <- c("carlson", "gnomad", "roulette")

#################
# 
# Genomewide correlations, part I
#
#################

maps_1Mb <- fread("maps/maps_1Mb.csv.gz")
maps_100kb <- fread("maps/maps_100kb.csv.gz")
maps_10kb <- fread("maps/maps_10kb.csv.gz") 

cor_list <- vector("list", length=3*3*2)
x <- 1
for(u in mutation_maps) {
  for(m in constrained_models) {

    gw_maps_1Mb <- maps_1Mb %>% filter(mut_map==u, constrained_model==m, coverage >= 0.25)
    piexp <- cor.test(gw_maps_1Mb$avg_pi, gw_maps_1Mb$exp_pi)$estimate
    B <- cor.test(gw_maps_1Mb$avg_pi, gw_maps_1Mb$B)$estimate
    mut <- cor.test(gw_maps_1Mb$avg_pi, gw_maps_1Mb$avg_mut)$estimate
    rec <- cor.test(gw_maps_1Mb$avg_pi, gw_maps_1Mb$avg_rec)$estimate
    Elem <- cor.test(gw_maps_1Mb$avg_pi, gw_maps_1Mb$del_sites)$estimate
    
    cor_list[[x]] <- cbind.data.frame(piexp, B, mut, rec, Elem, m, u)
    x <- x + 1
  }
}

cor_dt_1Mb <- data.table::rbindlist(cor_list)
names(cor_dt_1Mb) <- c("pi_exp", "B", "mut_rate", "rec_rate", "E", "elements", "mutation_map")
cor_dt_1Mb$scale <- 1e+6

cor_list <- vector("list", length=3*3*2)
x <- 1
for(u in mutation_maps) {
  for(m in constrained_models) {

    gw_maps_100kb <- maps_100kb %>% filter(mut_map==u, constrained_model==m, coverage >= 0.5)
    piexp <- cor.test(gw_maps_100kb$avg_pi, gw_maps_100kb$exp_pi)$estimate
    B <- cor.test(gw_maps_100kb$avg_pi, gw_maps_100kb$B)$estimate
    mut <- cor.test(gw_maps_100kb$avg_pi, gw_maps_100kb$avg_mut)$estimate
    rec <- cor.test(gw_maps_100kb$avg_pi, gw_maps_100kb$avg_rec)$estimate
    Elem <- cor.test(gw_maps_100kb$avg_pi, gw_maps_100kb$del_sites)$estimate
    
    cor_list[[x]] <- cbind.data.frame(piexp, B, mut, rec, Elem, m, u)
    x <- x + 1
  }
}

cor_dt_100kb <- data.table::rbindlist(cor_list)
names(cor_dt_100kb) <- c("pi_exp", "B", "mut_rate", "rec_rate", "E", "elements", "mutation_map")
cor_dt_100kb$scale <- 1e+5

cor_list <- vector("list", length=3*3*2)
x <- 1
for(u in mutation_maps) {
  for(m in constrained_models) {

    gw_maps_10kb <- maps_10kb %>% filter(mut_map==u, constrained_model==m, coverage >= 0.75)
    piexp <- cor.test(gw_maps_10kb$avg_pi, gw_maps_10kb$exp_pi)$estimate
    B <- cor.test(gw_maps_10kb$avg_pi, gw_maps_10kb$B)$estimate
    mut <- cor.test(gw_maps_10kb$avg_pi, gw_maps_10kb$avg_mut)$estimate
    rec <- cor.test(gw_maps_10kb$avg_pi, gw_maps_10kb$avg_rec)$estimate
    Elem <- cor.test(gw_maps_10kb$avg_pi, gw_maps_10kb$del_sites)$estimate
    
    cor_list[[x]] <- cbind.data.frame(piexp, B, mut, rec, Elem, m, u)
    x <- x + 1
  }
}

cor_dt_10kb <- data.table::rbindlist(cor_list)
names(cor_dt_10kb) <- c("pi_exp", "B", "mut_rate", "rec_rate", "E", "elements", "mutation_map")
cor_dt_10kb$scale <- 1e+4

cor_dt <- data.table(rbind.data.frame(cor_dt_10kb, cor_dt_100kb, cor_dt_1Mb))
m_cor_dt <- pivot_longer(cor_dt, cols=c("pi_exp", "B", "mut_rate", "rec_rate", "E"),
                         names_to="predictor", values_to="correlation")

m_cor_dt$elements <- factor(m_cor_dt$elements, 
                            levels = names(elem_labels), 
                            labels = unname(elem_labels)) 

m_cor_dt$mutation_map <- factor(m_cor_dt$mutation_map, 
                                levels = names(mut_labels), 
                                labels = unname(mut_labels))

p <- ggplot(m_cor_dt, aes(x=scale/1e3, y=correlation^2, color=predictor)) +
  facet_grid(mutation_map~elements) + theme_classic() +
  geom_line(data = m_cor_dt, linewidth = 0.75) +
  geom_point(size = 3) +
  scale_color_manual(name = NULL, 
                     values = c( B = "brown1", 
                                 E = "green4",
                                 mut_rate = "gold3",
                                 pi_exp = "purple3",
                                 rec_rate = "grey2"), 
                     labels = c( B = "B",
                                 E = "E", 
                                 mut_rate = expression(mu), 
                                 pi_exp = expression(pi[exp]),
                                 rec_rate = "r")) +
  scale_x_log10(breaks = c(10, 100, 1000)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  ylab(expression(R^2)) + xlab("Scale (kb)") +
  theme(strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom")
save_plot("sem/R2_pi_gw.pdf", p, base_height=8, base_width=15)

#################
# 
# Landscapes (main figures)
#
#################

# diversity (target) under equilibrium (fitted Ne)
c <- 22
dat_100kb <- maps_100kb %>%
  filter(mut_map=="roulette",
         constrained_model=="split_cds_phastcons",
         chrom == c, 
         coverage >= 0.5) %>%
  arrange(chromStart) %>%
  mutate(gap = chromStart - lag(chromStart, default = chromStart[1]))

dat_100kb[1, gap := 0]

dat_100kb[, contiguous := gap == 1e+5 | gap == 0]
dat_100kb[, run_id := rleid(contiguous)]
dat_100kb <- dat_100kb[contiguous == TRUE]

dat_100kb[, constrained_model := as.character(constrained_model)]
dat_100kb[, mut_map := as.character(mut_map)]

dat_100kb[, plot_group := paste0("r", run_id, "_cm_", constrained_model)]

pi_100kb <- ggplot(dat_100kb, aes(x = chromStart / 1e6)) +
  geom_path(data = dat_100kb, 
            aes(y = avg_pi, group = paste0("r", run_id)), linewidth = 0.75) +
  geom_point(data = dat_100kb,
             aes(y = avg_pi, group = paste0("r", run_id)), size = 1.5) +
  theme_classic() +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(title = NULL, x = "Position (Mb)", y = expression(pi)) +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 22))
save_plot("plots/c22pi100k.pdf", pi_100kb, base_width=15, base_height=5)

# plots landscapes (main figure)
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

dat_1Mb[, plot_group := paste0("r", run_id, "_cm_", constrained_model)]

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
                color = constrained_model,
                group = plot_group), linewidth = 0.75) +
  geom_point(data = dat_1Mb, shape=16,
             aes(y = exp_pi, 
                 color = constrained_model, 
                 group = plot_group), size = 1.5) +
  scale_x_continuous(breaks = seq(from=0, to=140, by=10), expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values = c("obs_pi" = "black",
                                "CDS + phastCons" = "cyan4",
                                "CDS + Regulatory" = "brown1"),
                     breaks = c("obs_pi", "CDS + phastCons", "CDS + Regulatory"),
                     labels = c("Observed",
                                "Expected (CDS + phastCons)",
                                "Expected (CDS + Regulatory)"),
                     name=NULL) +
  annotate("text", label = "Chromosome 8", size = 5,
           x = mean(range(dat_1Mb$chromStart / 1e6)),
           y = 0.95 * max(range(c(dat_1Mb$avg_pi, dat_1Mb$exp_pi)))) +
  labs(title = NULL, x = NULL, y = expression(pi)) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 22),
        legend.position = c(0.85, 0.85),
        legend.text=element_text(size = 16))

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

dat_100kb[, plot_group := paste0("r", run_id, "_cm_", constrained_model)]

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
                color = constrained_model,
                group = plot_group), linewidth = 0.75) +
  geom_point(data = dat_100kb,
             shape=16,
             aes(y = exp_pi,
                 color = constrained_model, 
                 group = plot_group), size = 1.5) +
  annotate("text", label = "Chromosome 22", size = 5,
           x = mean(range(dat_100kb$chromStart / 1e6)),
           y = 0.95 * max(range(c(dat_100kb$avg_pi, dat_100kb$exp_pi)))) +
  scale_x_continuous(breaks = seq(from=20, to=50, by=5), expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values = c("obs_pi" = "black",
                                "CDS + Regulatory" = "brown1",
                                "CDS + phastCons" = "cyan4"), name=NULL,
                     labels=c("Observed", "Expected (CDS + Regulatory)", "Expected (CDS + phastCons)")) +
  labs(title = NULL, x = "Position (Mb)", y = expression(pi)) +
  theme(strip.text.y = element_blank(),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.position="bottom")
save_plot("plots/chr22_100kb.pdf", p_100kb, base_height=6, base_width=15)
p_100kb <- p_100kb + theme(legend.position="none")
p <- plot_grid(p_1Mb, p_100kb, nrow=2, rel_heights=c(1, 1.15), labels="AUTO", align="v")
save_plot("plots/human_diversity_chr8_chr22.pdf", p, base_height=9, base_width=15)

#################
# 
# Landscapes (supplemental figures)
#
#################

maps_1Mb <- fread("maps/maps_1Mb.csv.gz")
for(c in 1:22) {
  dat_1Mb <- maps_1Mb %>%
    filter(constrained_model %in% constrained_models,
           mut_map=="roulette",
           chrom==c, coverage >= 0.25) %>%
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
  
  dat_1Mb[, plot_group := paste0("r", run_id, "_cm_", constrained_model)]
  
  q_1Mb <- ggplot(dat_1Mb, aes(x = chromStart / 1e6)) + theme_classic() +
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
                  color = constrained_model,
                  group = plot_group), linewidth = 0.75) +
    geom_point(data = dat_1Mb, shape=16,
               aes(y = exp_pi, 
                   color = constrained_model, 
                   group = plot_group), size = 1.5) +
    scale_x_continuous(breaks = pretty_breaks(), expand = expansion(mult = 0.01)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = c("obs_pi" = "black",
                                  "CDS + phastCons" = "cyan4",
                                  "CDS + Regulatory" = "brown1"),
                       breaks = c("obs_pi", "CDS + phastCons", "CDS + Regulatory"),
                       labels = c("Observed",
                                  "Expected (CDS + phastCons)",
                                  "Expected (CDS + Regulatory)"),
                       name=NULL) +
    labs(title = NULL, x = NULL, y = expression(pi)) +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 16),
          axis.title.y = element_text(size = 22),
          legend.position ="bottom",
          legend.text=element_text(size = 16))
  save_plot(paste("plots/pi_c", c, "_1Mb.pdf", sep=""), q_1Mb, base_height=6, base_width=12)
}

#################
# 
# Map comparison at 10 kb resolution
#
#################

c <- 22
dat_10kb <- maps_10kb %>%
  filter(chrom == c, coverage >= 0.75, chromStart >= 39.5e6, chromEnd <= 40.5e6) %>%
  arrange(chromStart) %>%
  mutate(gap = chromStart - lag(chromStart, default = chromStart[1]))

dat_10kb$constrained_model <- factor(dat_10kb$constrained_model, 
                                      levels = names(elem_labels), 
                                      labels = unname(elem_labels)) 
dat_10kb[1, gap := 0]

dat_10kb[, contiguous := gap == 1e+4 | gap == 0]
dat_10kb[, run_id := rleid(contiguous)]
dat_10kb <- dat_10kb[contiguous == TRUE]

dat_10kb[, constrained_model := as.character(constrained_model)]
dat_10kb[, mut_map := as.character(mut_map)]

focm <- c("CDS + Regulatory", "CDS + phastCons")
dat_10kb[, plot_group := paste0("r", run_id, "_cm_", mut_map)]
q_10kb <- ggplot(data=filter(dat_10kb, constrained_model %in% focm),
                 aes(x = chromStart / 1e6)) + theme_classic() +
  geom_path(data=filter(dat_10kb, constrained_model %in% focm),
            aes(y = exp_pi, 
                color = mut_map,
                group = plot_group), linewidth = 0.75) +
  geom_point(data=filter(dat_10kb, constrained_model %in% focm),
             shape=16,
             aes(y = exp_pi,
                 color = mut_map, 
                 group = plot_group), size = 1.5) +
  scale_x_continuous(breaks = pretty_breaks(), expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values = c("gnomad" = "seagreen4",
                                "roulette" = "brown1",
                                "carlson" = "cyan3"), name=NULL,
                     labels=c("Carlson",
                              "Gnomad",
                              "Roulette")) +
  labs(title = NULL, x = "Position (Mb)", y = expression(E(pi))) +
  facet_wrap(~constrained_model, ncol=1, strip.position = "right") + 
  theme(strip.text=element_text(size=16),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.position="bottom")
save_plot("plots/pi_c22_10kb.pdf", q_10kb, base_height=12, base_width=10)

dat_10kb[, plot_group := paste0("r", run_id, "_cm_", constrained_model)]
r_10kb <- ggplot(dat_10kb, aes(x = chromStart / 1e6)) + theme_classic() +
  geom_path(data = dat_10kb,
            aes(y = B, 
                color = constrained_model,
                group = plot_group), linewidth = 0.75) +
  geom_point(data = dat_10kb,
             shape=16,
             aes(y = B,
                 color = constrained_model, 
                 group = plot_group), size = 1.5) +
  scale_x_continuous(breaks = pretty_breaks(), expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values = c("CDS + Regulatory" = "purple3",
                                "CDS + phastCons" = "seagreen4",
                                "CDS + Regulatory" = "brown1",
                                "CDS + phastCons" = "cyan3"), name=NULL) +
  labs(title = NULL, x = "Position (Mb)", y = "B") +
  facet_wrap(~mut_map, ncol=1, strip.position = "right") + 
  theme(strip.text=element_text(size= 16),
        axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.position="bottom")
save_plot("plots/B_c22_10kb.pdf", r_10kb, base_height=12, base_width=10)

#################
# 
# simple scatterplots for slides
#
#################

r <- ggplot(filter(maps_100kb,
                   mut_map=="roulette",
                   constrained_model=="split_cds_phastcons",
                   chrom==22,
                   coverage >= 0.55), aes(y = avg_pi, x = avg_rec)) +
  theme_classic() + geom_point(size = 3) +
  scale_x_continuous(breaks = pretty_breaks()) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  ylab(expression(pi)) + xlab(expression(r)) +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 12))
save_plot("plots/rec_pi.pdf", r, base_height=6, base_width=10)

u <- ggplot(filter(maps_100kb,
                   mut_map=="roulette",
                   constrained_model=="split_cds_phastcons",
                   chrom==22,
                   coverage >= 0.5), 
            aes(y = avg_pi, x = avg_mut)) +
  theme_classic() + geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  ylab(expression(pi)) + xlab(expression(mu)) +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 12))
save_plot("plots/mut_pi.pdf", u, base_height=6, base_width=10)

b <- ggplot(filter(maps_100kb,
                   mut_map=="roulette",
                   constrained_model=="split_cds_phastcons",
                   chrom==22,
                   coverage >= 0.5), 
            aes(y = avg_pi, x = B)) +
  theme_classic() + geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  ylab(expression(pi)) + xlab(expression(B)) +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 12))
save_plot("plots/b_pi.pdf", b, base_height=6, base_width=10)

e <- ggplot(filter(maps_100kb,
                   mut_map=="roulette",
                   constrained_model=="split_cds_phastcons",
                   chrom==22,
                   coverage > 0.5), 
            aes(y = avg_pi, x = prop_del_sites)) +
  theme_classic() + geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  ylab(expression(pi)) + xlab("Density of constrained sites") +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 12))
save_plot("plots/e_pi.pdf", e, base_height=6, base_width=10)

#################
# 
# comparative correlations plot 
#
#################

tbl_1Mb <- maps_1Mb
cor_list <- vector("list", length=22*2*3) # chr * elements * mut maps
x <- 1
for(u in unique(tbl_1Mb$mut_map)) {
  for(m in constrained_models) {
    
    dat_pp <- filter(tbl_1Mb, mut_map==u, constrained_model==m, coverage >= 0.25) 
    
    # genome-wide 
    piexp <- cor.test(dat_pp$avg_pi, dat_pp$exp_pi)$estimate^2
    B <- cor.test(dat_pp$avg_pi, dat_pp$B)$estimate^2
    mut <- cor.test(dat_pp$avg_pi, dat_pp$avg_mut)$estimate^2
    rec <- cor.test(dat_pp$avg_pi, dat_pp$avg_rec)$estimate^2
    Elem <- cor.test(dat_pp$avg_pi, dat_pp$del_sites)$estimate^2
    
    gw <- cbind.data.frame(piexp, B, mut, rec, Elem, m, u)
    gw$data <- "genome-wide"
    
    # per chromosome
    chr_list <- vector("list", length=22)
    for(c in 1:22) {
  
      dat_c <- filter(dat_pp, chrom==c) 
      
      piexp <- cor.test(dat_c$avg_pi, dat_c$exp_pi)$estimate^2
      B <- cor.test(dat_c$avg_pi, dat_c$B)$estimate^2
      mut <- cor.test(dat_c$avg_pi, dat_c$avg_mut)$estimate^2
      rec <- cor.test(dat_c$avg_pi, dat_c$avg_rec)$estimate^2
      Elem <- cor.test(dat_c$avg_pi, dat_c$del_sites)$estimate^2
      
      chr <- cbind.data.frame(piexp, B, mut, rec, Elem, m, u)
      chr$data <- paste("chr", c, sep="")
      
      chr_list[[c]] <- chr
    }
    
    chr <- data.table::rbindlist(chr_list)
    cor_list[[x]] <- rbind.data.frame(gw, chr)
    x <- x + 1
  }
}

cor_dt_1Mb <- data.table::rbindlist(cor_list)
names(cor_dt_1Mb) <- c("pi_exp", "B", "mut_rate", "rec_rate", "E", "elements", "mutation_map", "data")
cor_dt_1Mb$scale <- 1e+6

cor_dt_1Mb$elements <- factor(cor_dt_1Mb$elements, 
                              levels = names(elem_labels), 
                              labels = unname(elem_labels)) 

cor_dt_1Mb$mutation_map <- factor(cor_dt_1Mb$mutation_map, 
                                  levels = names(mut_labels), 
                                  labels = unname(mut_labels))

m_r2_dt <- pivot_longer(cor_dt_1Mb, 
                        cols=c("pi_exp", "B", "mut_rate", "rec_rate", "E"),
                        names_to="variable",
                        values_to="rsqrd")

new_B <- c(0.6, 0.67) 
new_colors <- c("red", "purple") 
new_labels <- c("Murphy2022", "Buffalo2024") 

levels_var <- unique(m_r2_dt$variable) 
other_studies <- data.frame(variable = factor(rep("B", length(new_B)), levels = levels_var),
                            rsqrd = new_B, label = new_labels, color = new_colors)

dat_gw <- filter(m_cor_dt, scale==1e6, predictor!="pi_exp")
dat_gw$rsqrd <- dat_gw$correlation^2
dat_gw$variable <- dat_gw$predictor
dat_gw <- dplyr::select(dat_gw, c(variable, rsqrd, elements, mutation_map)) 
dat_gw$label <- "This study"
dat_gw$color <- "black"
dat_gw$rsqrd <- as.numeric(dat_gw$rsqrd)

piexp_gw <- filter(m_cor_dt, predictor=="pi_exp", scale==1e6) %>%
  dplyr::select(., c("elements", "mutation_map", "correlation"))
piexp_gw$variable <- "pi_exp"
piexp_gw$rsqrd <- piexp_gw$correlation^2
piexp_gw <- dplyr::select(piexp_gw, -correlation)

dat_gw <- full_join(dat_gw, piexp_gw)

dat_gw$label <- "This study"
dat_gw$color <- "black"
dat_gw$rsqrd <- as.numeric(dat_gw$rsqrd)

# fake-applying elements and mutation maps to other studies
other_studies_extended <- crossing(other_studies, 
                                   "CDS + phastCons",
                                   c("Roulette", "Gnomad", "Carlson"))
names(other_studies_extended)[5] <- "elements"
names(other_studies_extended)[6] <- "mutation_map"

new_points <- full_join(other_studies_extended, dat_gw) %>% setDT()

cp <- ggplot(data=m_r2_dt, aes(y = rsqrd, x = variable)) +
  facet_grid(mutation_map~elements) + theme_classic() + 
  geom_point(size = 3, shape=2, color="cyan3") +
  scale_x_discrete(labels = c( "mut_rate" = expression(mu),
                               "rec_rate" = expression(r),
                               "pi_exp" = expression(pi[exp] == B %*% mu),
                               "B" = expression(B),
                               "prop_del_sites" = expression(E))) +
  scale_y_continuous(breaks = pretty_breaks(), limits=c(0, 1)) +
  ylab(expression(R[pi]^2)) + xlab("Variable") +
  theme(strip.text = element_text(size = 16),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 14),
        legend.position="bottom") +
  geom_point(data = new_points, aes(x = variable, y = rsqrd, color = label), 
             size = 5, shape = 8, inherit.aes = FALSE) +
  scale_color_manual(name="Genome-wide", 
                     values = setNames(unique(new_points$color),
                                       unique(new_points$label)))
save_plot("sem/detailed_r2.pdf", cp, base_height=7, base_width=10)

#################
# 
# comparison between phastCons and regulatory models
#
#################

file_list <- list.files("../../data/pi_tables/roulette/")
pis_1kb <- data.table::rbindlist(lapply(paste0("../../data/pi_tables/roulette/", file_list), fread))
pis_1kb[, mut_map :="roulette"]
pis_1kb[, pos := chromStart + 500]

file_list <- list.files("split_cds_phastcons/roulette/maps_1kb/")
bs_1kb_phast <- data.table::rbindlist(lapply(paste0("split_cds_phastcons/roulette/maps_1kb/", file_list), fread))
bs_1kb_phast[, mut_map :="roulette"]
bs_1kb_phast[, element := "CDS + phastCons"]

file_list <- list.files("split_cds_regulatory/roulette/maps_1kb/")
bs_1kb_reg <- data.table::rbindlist(lapply(paste0("split_cds_regulatory/roulette/maps_1kb/", file_list), fread))
bs_1kb_reg[, mut_map :="roulette"]
bs_1kb_reg[, element := "CDS + regulatory"]

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

# 1 Mb scatterplots
maps_1Mb_m <- pivot_longer(maps_1Mb, cols=c("B", "exp_pi"), names_to="Var", values_to="Val")

x1 <- ggplot(filter(maps_1Mb,
                    mut_map=="roulette",
                    constrained_model %in% constrained_models,
                    coverage > 0.25), 
             aes(y = avg_pi, x = B, color=constrained_model)) +
  theme_classic() + geom_point(size = 3) +
  geom_smooth(data = ~ subset(.x, constrained_model == "split_cds_regulatory"),
              method = "lm", se = FALSE, color = "red") +
  geom_smooth(data = ~ subset(.x, constrained_model == "split_cds_phastcons"),
              method = "lm", se = FALSE, color = "green") +
  scale_color_manual(values = c("split_cds_regulatory" = "firebrick",
                                "split_cds_phastcons" = "seagreen4"), 
                     labels=c("CDS + phastCons",
                              "CDS + regulatory"),
                     name=NULL) +
  stat_cor(method = "pearson", 
           label.x.npc = "left", label.y.npc = "top",
           size = 6, show.legend = FALSE,
           output.type = "expression",
           aes(label = ..rr.label..)) +
  scale_x_continuous(breaks = pretty_breaks(), expand = expansion(mult = 0.05)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  ylab(expression(pi[obs])) + xlab(expression(B)) +
  theme(plot.margin = margin(r = 20),
        axis.title = element_text(size = 26),
        axis.text = element_text(size = 16),
        legend.position = c(0.2, 0.75),
        legend.text=element_text(size = 16),
        legend.box="vertical",
        legend.key.size=unit(12, "pt"),
        legend.title=element_text(size=16, face="bold"))

x2 <- ggplot(filter(maps_1Mb,
                    mut_map=="roulette",
                    constrained_model %in% constrained_models,
                    coverage > 0.25), 
             aes(y = avg_pi, x = exp_pi, color=constrained_model)) +
  theme_classic() + geom_point(size = 3) +
  geom_smooth(data = ~ subset(.x, constrained_model == "split_cds_regulatory"),
              method = "lm", se = FALSE, color = "red") +
  geom_smooth(data = ~ subset(.x, constrained_model == "split_cds_phastcons"),
              method = "lm", se = FALSE, color = "green") +
  scale_color_manual(values = c("split_cds_regulatory" = "firebrick",
                                "split_cds_phastcons" = "seagreen4"), 
                     labels=c("CDS + phastCons",
                              "CDS + regulatory"),
                     name=NULL) +
  stat_cor(method = "pearson", 
           label.x.npc = "left", label.y.npc = "top",
           size = 6, show.legend = FALSE,
           output.type = "expression",
           aes(label = ..rr.label..)) +
  scale_x_continuous(breaks = pretty_breaks(), expand = expansion(mult = 0.05)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  ylab(expression(pi[obs])) + xlab(expression(pi[exp])) +
  theme(plot.margin = margin(r = 20, l=20),
        axis.title = element_text(size = 26),
        axis.text = element_text(size = 16),
        legend.position = c(0.2, 0.75),
        legend.text=element_text(size = 16),
        legend.box="vertical",
        legend.key.size=unit(12, "pt"),
        legend.title=element_text(size=16, face="bold"))

mut_latent_1Mb <- fread("maps/mut_latent_1Mb.csv.gz") # mut latent
tab_1Mb <- merge(filter(maps_1Mb, 
                        mut_map=="roulette",
                        constrained_model == "split_cds_regulatory", # arbitrary, won't be used
                        coverage > 0.25),
                 mut_latent_1Mb) %>% dplyr::select(., c(chrom, avg_pi, exp_pi, avg_mut, mut_latent_log))
tab_1Mb[, c("Roulette", "Latent") := lapply(.SD, scale), .SDcols = c("avg_mut", "mut_latent_log")]
tab_1Mb_m <- pivot_longer(tab_1Mb, cols=c("Roulette", "Latent"), names_to="Map", values_to="Rate")

x3 <- ggplot(tab_1Mb_m, aes(y = avg_pi, x = Rate, color=Map)) +
  theme_classic() + geom_point(size = 3) +
  geom_smooth(data = ~ subset(.x, Map == "Roulette"),
              method = "lm", se = FALSE, color = "cyan3") +
  geom_smooth(data = ~ subset(.x, Map == "Latent"),
              method = "lm", se = FALSE, color = "gold2") +
  scale_color_manual(values = c("Roulette" = "steelblue4", "Latent" = "orange2"), name=NULL) +
  stat_cor(method = "pearson", 
           label.x.npc = "left", label.y.npc = "top",
           size = 6, show.legend = FALSE,
           output.type = "expression",
           aes(label = ..rr.label..)) +
  scale_x_continuous(breaks = pretty_breaks(), expand = expansion(mult = 0.05)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  ylab(expression(pi[obs])) + xlab(expression("Standardized " * mu)) +
  theme(plot.margin = margin(l = 20),
        axis.title = element_text(size = 26),
        axis.text = element_text(size = 16),
        legend.position = c(0.125, 0.75),
        legend.text=element_text(size = 16),
        legend.box="vertical",
        legend.key.size=unit(12, "pt"),
        legend.title=element_text(size=16, face="bold"))
xp <- plot_grid(x1, x2, x3, nrow=1, align="h", labels=c("A", "B", "C"))
save_plot("plots/genome_wide_cors_2.pdf", xp, base_height=7, base_width=24)

# output by fit_obspi.R
r2_sem <- fread("sem/obspi_chr_r2pi.csv") %>% filter(., DAG==4) %>% 
  dplyr::select(., -DAG) %>% setDT()
r2_sem$method <- "sem"

# Pearson correlations per chr at different scales
elem_list <- vector("list", length=2) # elements 
x <- 1
for(m in constrained_models) {

  chr_list <- vector("list", length=22)
  for(c in 1:22) {
    
    dat_1Mb <- filter(maps_1Mb, mut_map=="roulette", constrained_model==m, coverage >= 0.25, chrom==c) 
    dat_100kb <- filter(maps_100kb, mut_map=="roulette", constrained_model==m, coverage >= 0.5, chrom==c) 
    dat_10kb <- filter(maps_10kb, mut_map=="roulette", constrained_model==m, coverage >= 0.75, chrom==c) 
    
    piexp_1Mb <- cor.test(dat_1Mb$avg_pi, dat_1Mb$exp_pi)$estimate^2
    piexp_100kb <- cor.test(dat_100kb$avg_pi, dat_100kb$exp_pi)$estimate^2
    piexp_10kb <- cor.test(dat_10kb$avg_pi, dat_10kb$exp_pi)$estimate^2

    tmp <- rbind.data.frame(piexp_1Mb, piexp_100kb, piexp_10kb)
    names(tmp) <- "avg_pi"
    tmp$chrom <- c
    tmp$constrained_model <- m
    tmp$scale <- c(1e6, 1e5, 1e4)
    
    chr_list[[c]] <- tmp
  }
  
  elem_list[[x]] <- data.table::rbindlist(chr_list)
  x <- x + 1
}

r2_pearson <- data.table::rbindlist(elem_list)
r2_pearson$method <- "pearson"

r2_pearson$constrained_model <- factor(r2_pearson$constrained_model,
                                       levels = names(elem_labels), 
                                       labels = unname(elem_labels)) 

r2_methods <- rbind.data.frame(r2_sem, r2_pearson)
r2_methods <- r2_methods %>% 
  mutate(scale_offset = case_when(method == "sem" ~ scale * 1.2, 
                                  method == "pearson" ~ scale * 0.8,
                                  TRUE ~ scale))

p1 <- ggplot(data=filter(r2_methods, constrained_model=="CDS + phastCons"),
            aes(x=scale_offset/1e3, y=avg_pi, shape=method)) +
  geom_point(size=3) + theme_classic() + 
  scale_shape_manual(values=c(0, 1),
                     labels=c(expression(R[pi]^2 * "  " * "Pearson"), expression(R[pi]^2 * "  " * "SEM")), 
                     name=NULL) +
  geom_line(aes(group = interaction(chrom, scale, constrained_model)),
            linewidth = 0.6, alpha = 0.4) +
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_x_log10(breaks=c(10, 100, 1000), expand = expansion(mult = 0.05)) +
  labs(title=NULL, x="Scale (kb)", y=expression(R[pi]^2 * "  " * "(CDS + phastCons)")) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=14), 
        axis.text.x = element_text(margin=margin(t=10)),
        legend.position = c(0.15, 0.8),
        legend.text=element_text(size = 16),
        legend.box="vertical",
        legend.key.size=unit(12, "pt"),
        legend.title=element_text(size=16, face="bold"))

p2 <- ggplot(data=filter(r2_methods, constrained_model=="CDS + Regulatory"),
             aes(x=scale_offset/1e3, y=avg_pi, shape=method)) +
  geom_point(size=3) + theme_classic() + 
  scale_shape_manual(values=c(0, 1), 
                     labels=c(expression(R[pi]^2 * "  " * "Pearson"), expression(R[pi]^2 * "  " * "SEM")), 
                     name=NULL) +
  geom_line(aes(group = interaction(chrom, scale, constrained_model)),
            linewidth = 0.6, alpha = 0.4) +
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_x_log10(breaks=c(10, 100, 1000), expand = expansion(mult = 0.05)) +
  labs(title=NULL, x="Scale (kb)", y=expression(R[pi]^2 * "  " * "(CDS + Regulatory)")) +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=14), 
        axis.text.x = element_text(margin=margin(t=10)),
        legend.position = c(0.15, 0.8),
        legend.text=element_text(size = 16),
        legend.box="vertical",
        legend.key.size=unit(12, "pt"),
        legend.title=element_text(size=16, face="bold"))

y <- plot_grid(hist, p1, p2, x1, x2, x3, nrow=2, label_size = 18, align="hv", labels="AUTO")
save_plot("plots/collage.pdf", y, base_width=24, base_height=12)

p3 <- ggplot(data=r2_methods, aes(x=scale_offset/1e3, y=avg_pi, shape=method)) +
  geom_point(size=3) + theme_classic() + facet_wrap(~constrained_model) +
  scale_shape_manual(values=c(0, 1), 
                     labels=c(expression(R^2 * "(Pearson)"), expression(R[pi]^2 * "(SEM)")), 
                     name=NULL) +
  geom_line(aes(group = interaction(chrom, scale, constrained_model)),
            linewidth = 0.6, alpha = 0.4) +
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_x_log10(breaks=c(10, 100, 1000), expand = expansion(mult = 0.05)) +
  labs(title=NULL, x="Scale (kb)", y="Variance Explained") +
  theme(strip.text=element_text(size=16),
        axis.title=element_text(size=16), 
        axis.text=element_text(size=14), 
        axis.text.x = element_text(margin=margin(t=10)),
        legend.position = "bottom",
        legend.text=element_text(size = 16),
        legend.box="vertical",
        legend.key.size=unit(12, "pt"),
        legend.title=element_text(size=16, face="bold"))
save_plot("plots/r2_comparison.pdf", p3, base_height=5, base_width=10)

#################
# 
# AOV
#
#################

# 1Mb
maps_1Mb_reg <- fread("maps/maps_1Mb.csv.gz") %>%
  filter(., mut_map=="roulette", constrained_model=="split_cds_regulatory") %>%
  dplyr::select(., c("chrom", "chromStart", "chromEnd",  "avg_pi",
                     "avg_mut", "avg_rec", "B", "exp_pi", "prop_del_sites", "coverage"))

maps_1Mb_phast <- fread("maps/maps_1Mb.csv.gz") %>%
  filter(., mut_map=="roulette", constrained_model=="split_cds_phastcons") %>%
  dplyr::select(., c("chrom", "chromStart", "chromEnd", "avg_pi",
                     "avg_mut", "avg_rec", "B", "exp_pi", "prop_del_sites", "coverage"))

# lifting over AOV regions from hg19 to 38
chain <- import.chain("../../data/hg19ToHg38.over.chain")
aov <- fread("../../data/Gilbert2020_AOV.csv") %>% 
  dplyr::select(., c("Chromosome", "START", "END"))
colnames(aov) <- c("seqnames", "start", "end")
aov$seqnames <- paste0("chr", aov$seqnames)

aov_gr <- makeGRangesFromDataFrame(aov)

lifted <- liftOver(aov_gr, chain)
lifted_aov_gr <- unlist(lifted)

maps_1Mb_reg$chromStart <- maps_1Mb_reg$chromStart + 1
maps_1Mb_phast$chromStart <- maps_1Mb_phast$chromStart + 1

maps_1Mb_reg$chrom <- paste0("chr", maps_1Mb_reg$chrom)
maps_1Mb_phast$chrom<- paste0("chr", maps_1Mb_phast$chrom)

maps_1Mb_reg_gr <- makeGRangesFromDataFrame(maps_1Mb_reg, keep.extra.columns=T)
maps_1Mb_phast_gr <- makeGRangesFromDataFrame(maps_1Mb_phast, keep.extra.columns=T)

# independent of constrained element type
hits <- findOverlaps(query=lifted_aov_gr, subject=maps_1Mb_reg_gr) 
hit_list <- as.data.frame(hits)
names(hit_list) <- c("aov_bin", "bin_1Mb")

aov$aov_bin <- 1:nrow(aov)
maps_1Mb_reg$bin_1Mb <- 1:nrow(maps_1Mb_reg)
maps_1Mb_phast$bin_1Mb <- 1:nrow(maps_1Mb_phast)

aov_genes <- merge(maps_1Mb_reg, hit_list) %>% 
  left_join(., dplyr::select(aov, -c("seqnames", "start", "end")), by="aov_bin")

aov_phast <- merge(maps_1Mb_phast, hit_list) %>% 
  left_join(., dplyr::select(aov, -c("seqnames", "start", "end")), by="aov_bin")

tbl_aov_genes <- left_join(maps_1Mb_reg, aov_genes)
tbl_aov_phast <- left_join(maps_1Mb_phast, aov_phast)

setkey(tbl_aov_genes, c(chrom, chromStart, chromEnd))
tbl_aov_genes <- foverlaps(
  tbl_aov_genes[, .(chrom, chromStart, chromEnd = chromStart + 1e6 - 1)],
  tbl_aov_genes[, .(chrom, chromStart, chromEnd = chromStart + 1e6 - 1)],
  type = "any",
  mult = "first"
)

tbl_aov_genes <- tbl_aov_genes %>% filter(., coverage >= 0.25) %>%
  mutate(gap = chromStart - lag(chromStart, default = chromStart[1]))

tbl_aov_genes[1, gap := 0]

tbl_aov_genes[, contiguous := gap == 1e+6 | gap == 0]
tbl_aov_genes[, run_id := rleid(contiguous)]
tbl_aov_genes <- tbl_aov_genes[contiguous == TRUE]

tbl_aov_genes$constrained_model <- "Regulatory"

tbl_aov_phast <- tbl_aov_phast %>% filter(., coverage >= 0.25) %>%
  mutate(gap = chromStart - lag(chromStart, default = chromStart[1]))

tbl_aov_phast[1, gap := 0]

tbl_aov_phast[, contiguous := gap == 1e+6 | gap == 0]
tbl_aov_phast[, run_id := rleid(contiguous)]
tbl_aov_phast <- tbl_aov_phast[contiguous == TRUE]

tbl_aov_phast$constrained_model <- "phastCons"

tbl_aov <- setDT(rbind.data.frame(tbl_aov_genes, tbl_aov_phast))
tbl_aov[, plot_group := paste0("r", run_id, "_cm_", constrained_model)]

for(c in unique(tbl_aov$chrom)) {
  
  tbl_aov_c <- filter(tbl_aov, chrom==c) |> arrange(chromStart)
  
  p_rec <- p_rec <- ggplot(tbl_aov_c, aes(x = chromStart/1e6, y = avg_rec * 1e8)) +
    geom_step(color = "burlywood", aes(group = paste0("r", run_id))) +
    theme_classic() +
    labs(x = NULL, y = "r") +
    geom_rect(data = tbl_aov_c |> filter(!is.na(aov_bin)),
              aes(xmin = chromStart / 1e6, xmax = chromEnd / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          legend.position="none")
  
  p_mut <- ggplot(tbl_aov_c, aes(x = chromStart / 1e6, y = avg_mut * 1e8)) +
    geom_step(color = "orchid", aes(group = paste0("r", run_id))) +
    theme_classic() +
    geom_rect(data = tbl_aov_c |> filter(!is.na(aov_bin)),
              aes(xmin = chromStart / 1e6, xmax = chromEnd / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
    labs(x = NULL, y = expression(mu)) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          legend.position="none")
  
  p_sel <- ggplot(tbl_aov_c, aes(x = chromStart / 1e6, y = prop_del_sites * 100)) +
    theme_classic() +
    geom_step(aes(color = constrained_model, group = plot_group)) +
    scale_color_manual(values=c("Genes" = "brown1", "phastCons" = "green4")) +
    geom_rect(data = tbl_aov_c |> filter(!is.na(aov_bin)),
              aes(xmin = chromStart / 1e6, xmax = chromEnd / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
    labs(x = NULL, y = " % del. sites") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          legend.position="none")
  
  p_B <- ggplot(tbl_aov_c, aes(x = chromStart / 1e6, y = B)) +
    theme_classic() +
    geom_step(aes(color = constrained_model, group = plot_group)) +
    scale_color_manual(values=c("Genes" = "brown1", "phastCons" = "green4")) +
    geom_rect(data = tbl_aov_c |> filter(!is.na(aov_bin)),
              aes(xmin = chromStart / 1e6, xmax = chromEnd / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
    labs(x = NULL, y = "B") +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          legend.position="none")
  
  p_pi <- ggplot(tbl_aov_c, aes(x = chromStart / 1e6)) + theme_classic() +
    geom_rect(data = tbl_aov_c |> filter(!is.na(aov_bin)),
              aes( xmin = chromStart / 1e6, xmax = chromEnd / 1e6, ymin = -Inf, ymax = Inf ),
              fill = "grey80", alpha = 0.4, inherit.aes = FALSE ) +
    geom_line(data = tbl_aov_c, 
              aes(y = avg_pi * 1e3, 
                  color = "obs_pi", 
                  group = paste0("r", run_id)), linewidth = 0.75) +
    geom_point(data = tbl_aov_c,
               aes(y = avg_pi * 1e3,
                   color = "obs_pi", 
                   group = paste0("r", run_id)), size = 1.5) +
    geom_line(data = tbl_aov_c,
              aes(y = exp_pi * 1e3, 
                  color = constrained_model,
                  group = plot_group), linewidth = 0.75) +
    geom_point(data = tbl_aov_c,
               shape=2,
               aes(y = exp_pi * 1e3,
                   color = constrained_model, 
                   group = plot_group), size = 2) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = c("obs_pi" = "black",
                                  "Genes" = "brown1",
                                  "phastCons" = "green4"),
                       name=NULL) +
    labs(title = NULL, x = "Position (Mb)", y = expression(pi)) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 22),
          legend.position = "bottom")
  
  final_plot <- p_rec / p_mut / p_sel / p_B / p_pi + 
    plot_layout(heights = c(1, 1, 1, 1, 6)) &
    theme(plot.margin = unit(c(0,0,0,0), "pt"), panel.spacing = unit(0, "pt"))
  
  save_plot(paste("plots/aov_1Mb_", c, ".pdf", sep=""), final_plot, base_height=13, base_width=16)
}
