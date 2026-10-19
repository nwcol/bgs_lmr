
#################
# 
# This script analyses differences in genomic features across chromosomes
#
#
#################

source("~/Devel/bgs_lmr/tools/R_scripts/HelperFunctions.R")

maps_10kb <- fread("models/tables_10kb/equilibrium_granular_Ne/maps/maps_10kb.csv.gz")
maps_100kb <- fread("models/tables_10kb/equilibrium_granular_Ne/maps/maps_100kb.csv.gz")
maps_1Mb <- fread("models/tables_10kb/equilibrium_granular_Ne/maps/maps_1Mb.csv.gz")

################### mutation
# summarizing in terms of AC and CV (not for uniform mut maps as they are flat)
## 10 kb
human_summaries_10kb <- maps_10kb %>% filter(., constrained_model=="split_cds_phastcons") %>% # arbitrary, features are identical wrt elements
  group_by(chrom, mut_map, scale) %>%
  mutate(cv_mut=sd(avg_mut)/mean(avg_mut)) %>%
  mutate(cv_rec=sd(avg_rec)/mean(avg_rec)) %>%
  mutate(ac_mut=as.numeric(unlist(acf(avg_mut, lag=1, pl=F))[2])) %>%
  mutate(ac_rec=as.numeric(unlist(acf(avg_rec, lag=1, pl=F))[2])) %>%
  distinct(chrom, .keep_all=T) %>% ungroup() %>% setDT()
human_summaries_10kb <- dplyr::select(human_summaries_10kb, c(chrom, cv_mut, cv_rec, ac_mut, ac_rec, scale, mut_map))

## 100 kb
human_summaries_100kb <- maps_100kb %>% filter(., constrained_model=="split_cds_phastcons") %>% # arbitrary, features are identical wrt elements
  group_by(chrom, mut_map, scale) %>%
  mutate(cv_mut=sd(avg_mut)/mean(avg_mut)) %>%
  mutate(cv_rec=sd(avg_rec)/mean(avg_rec)) %>%
  mutate(ac_mut=as.numeric(unlist(acf(avg_mut, lag=1, pl=F))[2])) %>%
  mutate(ac_rec=as.numeric(unlist(acf(avg_rec, lag=1, pl=F))[2])) %>%
  distinct(chrom, .keep_all=T) %>% ungroup() %>% setDT()
human_summaries_100kb <- dplyr::select(human_summaries_100kb, c(chrom, cv_mut, cv_rec, ac_mut, ac_rec, scale, mut_map))

## 1 Mb
human_summaries_1Mb <- maps_1Mb %>% filter(., constrained_model=="split_cds_phastcons") %>% # arbitrary, features are identical wrt elements
  group_by(chrom, mut_map, scale) %>%
  mutate(cv_mut=sd(avg_mut)/mean(avg_mut)) %>%
  mutate(cv_rec=sd(avg_rec)/mean(avg_rec)) %>%
  mutate(ac_mut=as.numeric(unlist(acf(avg_mut, lag=1, pl=F))[2])) %>%
  mutate(ac_rec=as.numeric(unlist(acf(avg_rec, lag=1, pl=F))[2])) %>%
  distinct(chrom, .keep_all=T) %>% ungroup() %>% setDT()
human_summaries_1Mb <- dplyr::select(human_summaries_1Mb, c(chrom, cv_mut, cv_rec, ac_mut, ac_rec, scale, mut_map))

human_summaries <- rbind.data.frame(human_summaries_10kb, human_summaries_100kb, human_summaries_1Mb)

# Plotting 2D mutation maps (AC x CV)
p <- ggplot(data=human_summaries, aes(y=cv_mut, x=ac_mut)) + theme_bw() +
  facet_grid(mut_map~scale/1e3, scale="free") +
  ggrepel::geom_text_repel(aes(label=chrom), vjust=0.3, size=3, max.overlaps=Inf) +
  # https://stats.stackexchange.com/questions/185425/how-to-determine-the-critical-values-of-acf
  #geom_vline(xintercept=-1.96 / sqrt(100-1), color="grey", linetype="dashed") + 
  #geom_vline(xintercept=1.96 / sqrt(100-1), color="grey", linetype="dashed") + 
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_x_continuous(breaks=pretty_breaks()) +
  labs(title=NULL, x="Autocorrelation", y="Coefficient of Variation") +
  theme(strip.text.y=element_text(size=16),        
        strip.text.x=element_text(size=16),        
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("2D_plot_mut.pdf", p, base_height=6, base_width=12)

################### recombination
# decode map
maps_decode <- list.files("~/Data/bgs_lmr/PNAS_reviews/deCODE/rec_maps_formatted/", full.names=T)
maps_decode <- data.table::rbindlist(lapply(maps_decode, fread))

decode_10kb <- maps_decode[, .(mean_rate = mean(`Rate(cM/Mb)`, na.rm = TRUE)),
                             by = .(Chromosome, start=((`Position(bp)` - 1) %/% 1e4) * 1e4, end=(((`Position(bp)` - 1) %/% 1e4) + 1) * 1e4)]
decode_10kb[, scale := 1e4]

decode_100kb <- maps_decode[, .(mean_rate = mean(`Rate(cM/Mb)`, na.rm = TRUE)),
                           by = .(Chromosome, start=((`Position(bp)` - 1) %/% 1e5) * 1e5, end=(((`Position(bp)` - 1) %/% 1e5) + 1) * 1e5)]
decode_100kb[, scale := 1e5]

decode_1Mb <- maps_decode[, .(mean_rate = mean(`Rate(cM/Mb)`, na.rm = TRUE)),
                            by = .(Chromosome, start=((`Position(bp)` - 1) %/% 1e6) * 1e6, end=(((`Position(bp)` - 1) %/% 1e6) + 1) * 1e6)]
decode_1Mb[, scale := 1e6]

decode <- rbind.data.frame(decode_10kb, decode_100kb, decode_1Mb)
decode[, Chromosome := sub("^chr", "", Chromosome)]
setnames(decode, old="Chromosome", new="chrom")

decode_summaries <- decode %>%
  group_by(chrom, scale) %>%
  mutate(cv_rec=sd(mean_rate)/mean(mean_rate)) %>%
  mutate(ac_rec=as.numeric(unlist(acf(mean_rate, lag=1, pl=F))[2])) %>%
  distinct(chrom, .keep_all=T) %>% ungroup() %>% setDT()
decode_summaries <- dplyr::select(decode_summaries, c(chrom, cv_rec, ac_rec, scale))
decode_summaries[, rec_map := "deCODE"]

pyrho_summaries <- dplyr::filter(human_summaries, mut_map=="roulette") %>% dplyr::select(c(chrom, cv_rec, ac_rec, scale)) %>% setDT()
pyrho_summaries[, rec_map := "pyrho"]

rec_summaries <- rbind.data.frame(decode_summaries, pyrho_summaries)

# Plotting 2D recombination maps (AC x CV)
p <- ggplot(data=rec_summaries, aes(y=cv_rec, x=ac_rec)) + theme_bw() +
  facet_grid(rec_map~scale/1e3, scale="free") +
  ggrepel::geom_text_repel(aes(label=chrom), vjust=0.3, size=3, max.overlaps=Inf) +
  # https://stats.stackexchange.com/questions/185425/how-to-determine-the-critical-values-of-acf
  #geom_vline(xintercept=-1.96 / sqrt(100-1), color="grey", linetype="dashed") + 
  #geom_vline(xintercept=1.96 / sqrt(100-1), color="grey", linetype="dashed") + 
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_x_continuous(breaks=pretty_breaks()) +
  labs(title=NULL, x="Autocorrelation", y="Coefficient of Variation") +
  theme(strip.text.y=element_text(size=16),        
        strip.text.x=element_text(size=16),        
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("2D_plot_rec.pdf", p, base_height=6, base_width=12)
