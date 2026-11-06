
#################
# 
# This script uses prepared tables to make some landscape plots
# NOTE: this script should be run from within each current_models dir
#
#################

source("~/Devel/bgs_lmr/tools/R_scripts/HelperFunctions.R")

constrained_models <- c("split_cds_phastcons", "split_cds_regulatory")
mutation_maps <- c("carlson", "gnomad", "roulette")

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
  annotate("text", label = "Chromosome 22", size = 5,
           x = mean(range(dat_100kb$chromStart / 1e6)),
           y = 0.95 * max(range(c(dat_100kb$avg_pi, dat_100kb$exp_pi)))) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(title = NULL, x = "Position (Mb)", y = expression(pi)) +
  theme(axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 22))
save_plot("plots/c22pi100k.pdf", pi_100kb, base_width=15, base_height=6)

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
                     labels=c("Expected (CDS + phastCons)",  "Expected (CDS + Regulatory)", "Observed")) +
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
