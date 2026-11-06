
setwd("~/Data/bgs_lmr/PNAS_reviews/DAGs5-6/")

source("~/Devel/bgs_lmr/tools/R_scripts/HelperFunctions.R")

constrained_models <- c("split_cds_phastcons", "split_cds_regulatory")

#################
#
# Drawing DAGs
#
#################

coords_df <- data.frame( x = c(pi = 1, B = 2, E = 0, u = 0.5, r = 1.5),
                         y = c(pi = 1, B = 2, E = 2, u = 3, r = 4), 
                         stringsAsFactors = FALSE )

DAG5 <- dagify(
  pi ~ E,
  pi ~ B,
  pi ~ u,
  pi ~ r,
  #B ~ u,
  #B ~ r,
  B ~ E
  #u ~ r
)

td <- tidy_dagitty(DAG5) 
nodes <- unique(td$data$name)

coords_df_clean <- coords_df[nodes, , drop = FALSE ]

coords_list <- list(x = setNames(as.numeric(coords_df_clean$x), rownames(coords_df_clean)), 
                    y = setNames(as.numeric(coords_df_clean$y), rownames(coords_df_clean)))

coordinates(DAG5) <- coords_list
coordinates(DAG5)

DD5p <- DAG5 %>% 
  tidy_dagitty() %>%
  arrange(name) %>% 
  ggplot(aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_dag_node() +
  geom_dag_edges() +
  theme_dag() + 
  geom_dag_text(parse=TRUE, 
                label = c("B", "E", "pi", "r", "mu"))

DAG6 <- dagify(
  pi ~ E,
  pi ~ B,
  pi ~ u,
  #B ~ u,
  r ~ B,
  B ~ E,
  r ~ u
)

td <- tidy_dagitty(DAG6) 
nodes <- unique(td$data$name)

coords_df_clean <- coords_df[nodes, , drop = FALSE ]

coords_list <- list(x = setNames(as.numeric(coords_df_clean$x), rownames(coords_df_clean)), 
                    y = setNames(as.numeric(coords_df_clean$y), rownames(coords_df_clean)))

coordinates(DAG6) <- coords_list
coordinates(DAG6)

DD6p <- DAG6 %>% 
  tidy_dagitty() %>%
  arrange(name) %>% 
  ggplot(aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_dag_node() +
  geom_dag_edges() +
  theme_dag() + 
  geom_dag_text(parse=TRUE, 
                label = c("B", "E", "pi", "r", "mu"))

x <- plot_grid(DD5p, DD6p, nrow=1, labels=5:6)
save_plot("DAGs5-6.pdf", x, base_width=10, base_height=5)

#################
#
# Specifying causal models (DAGs)
#
#################

model_5 <- "
    avg_pi ~ start(0.25)*c1*avg_rec
    avg_pi ~ start(0.25)*c2*avg_mut
    avg_pi ~ start(0.5)*c3*B
    avg_pi ~ start(-0.25)*c4*prop_del_sites
    B ~ start(-0.25)*c5*prop_del_sites
    "

model_6 <- "
    avg_pi ~ start(0.25)*c1*avg_mut
    avg_pi ~ start(0.5)*c2*B
    avg_pi ~ start(-0.25)*c3*prop_del_sites
    B ~ start(-0.25)*c4*prop_del_sites
    avg_rec ~ start(0.25)*c6*B
    avg_rec ~ start(0.25)*c7*avg_mut
    "
#B ~ start(-0.25)*c5*avg_mut
#################
#
# Loading real data
#
#################

path <- "~/Devel/bgs_lmr/models/tables_10kb/equilibrium_granular_Ne/"

human_maps_1Mb <- fread(paste0(path, "maps/maps_1Mb.csv.gz")) %>% filter(., mut_map=="roulette") 
mut_latent_1Mb <- fread(paste0(path, "maps/mut_latent_1Mb.csv.gz"))
B_latent_1Mb <- data.frame()

for(m in constrained_models) {
  tmp_1Mb <- fread(paste0(path, "maps/B_latent_", m, "_1Mb.csv.gz"))
  B_latent_1Mb <- rbind.data.frame(B_latent_1Mb, tmp_1Mb)
}

#################
#
# Fitting to real data (1 Mb)
#
#################

for(m in constrained_models) { 
  
  print(paste("1Mb", m, Sys.Date(), Sys.time()), sep="\t")
  
  tmp_B <- filter(B_latent_1Mb, constrained_model==m) %>% dplyr::select(., -constrained_model)
  
  tbl_1Mb <- filter(human_maps_1Mb, constrained_model==m, coverage >= 0.2)  
  table_1Mb_g <- tbl_1Mb[,1:3]
  table_1Mb_g[, avg_pi := log(tbl_1Mb$avg_pi)]
  table_1Mb_g[, avg_rec := log(tbl_1Mb$avg_rec)] # pyrho YRI pop.
  table_1Mb_g[, prop_del_sites := log1p(tbl_1Mb$prop_del_sites)]
  
  table_1Mb_g <- left_join(table_1Mb_g, tmp_B, by=c("chrom", "chromStart", "chromEnd"))
  table_1Mb_g <- left_join(table_1Mb_g, mut_latent_1Mb, by=c("chrom", "chromStart", "chromEnd"))
  table_1Mb_g$chrom <- as.integer(table_1Mb_g$chrom)
  table_1Mb_g <- add_genome_offsets_and_midpoint(table_1Mb_g, chrom_lengths=NULL)
  table_1Mb_g$chrom <- as.integer(table_1Mb_g$chrom)
  
  setnames(table_1Mb_g, old=c("B_latent_log", "mut_latent_log"), new=c("B", "avg_mut"))
  
  dt <- copy(table_1Mb_g) %>% arrange(chrom, midpoint) %>% as.data.table() %>% setDT()
  
  for(c in 1:22) { # one model per chrom
    
    m5_fit <- sem(model_5,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",              
                  missing = "fiml",             
                  meanstructure = FALSE,            
                  std.ov = TRUE,     
                  fixed.x = TRUE,
                  optim.method = "nlminb",    
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m5_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m5_fit, "r2")),
                              as.numeric(inspect(m5_fit, "r2"))), 
             paste("sem/lavaan_out/obspi_YRI_DAG5_r2_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(parameterestimates(m5_fit, ci=T), 
             paste("sem/lavaan_out/obspi_YRI_DAG5_params_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m5_fit)),
                              as.numeric(fitMeasures(m5_fit))), 
             paste("sem/lavaan_out/obspi_YRI_DAG5_fitmeasures_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
    } else {
      print(paste("Model 5 did not converge for chr", c, "(1 Mb)"))
    }
    
    m6_fit <- sem(model_6,
                  data=filter(dt, chrom==c),
                  estimator = "MLR",              
                  missing = "fiml",             
                  meanstructure = FALSE,            
                  std.ov = TRUE,     
                  fixed.x = TRUE,
                  optim.method = "nlminb",    
                  control = list(trace = 0, rel.tol = 1e-8, iter.max = 1000))
    
    if(lavInspect(m6_fit, "converged")) {
      fwrite(rbind.data.frame(names(inspect(m6_fit, "r2")),
                              as.numeric(inspect(m6_fit, "r2"))), 
             paste("sem/lavaan_out/obspi_YRI_DAG6_r2_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
      
      fwrite(parameterestimates(m6_fit, ci=T), 
             paste("sem/lavaan_out/obspi_YRI_DAG6_params_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""))
      
      fwrite(rbind.data.frame(names(fitMeasures(m6_fit)),
                              as.numeric(fitMeasures(m6_fit))), 
             paste("sem/lavaan_out/obspi_YRI_DAG6_fitmeasures_1Mb_chr", 
                   c, "_latents_", m, ".csv", sep=""),
             col.names=F)
    } else {
      print(paste("Model 6 did not converge for chr", c, "(1 Mb)"))
    }
  }
}

#################
#
# Plotting real data (1 Mb)
#
#################

chr_fitmeasures_1Mb <- data.frame()
for(dag in 5:6) {
  for(m in constrained_models) { 
    chr_list <- list(length=num_chrom)
    for(c in 1:num_chrom) {
      
      tmp <- fread(paste("sem/lavaan_out/obspi_YRI_DAG", dag, 
                         "_fitmeasures_1Mb_chr", 
                         c, "_latents_", m, ".csv", sep=""))
      
      tmp$chrom <- c
      tmp$DAG <- dag
      tmp$constrained_model <- m
      
      chr_list[[c]] <- tmp
    }
    
    df <- do.call(rbind.data.frame, chr_list)
    chr_fitmeasures_1Mb <- rbind.data.frame(chr_fitmeasures_1Mb, df)
  }
}

chr_fit_1Mb <- dplyr::select(chr_fitmeasures_1Mb, c("pvalue", "aic", "bic", "cfi", "gfi", "chrom", "DAG", "constrained_model"))

chr_fit_1Mb$constrained_model <- factor(chr_fit_1Mb$constrained_model, 
                                        levels = names(elem_labels), 
                                        labels = unname(elem_labels)) 

chr_fit_1Mb <- chr_fit_1Mb %>% 
  mutate(DAG_offset = case_when(constrained_model == "CDS + phastCons" ~ DAG * 0.975, 
                                constrained_model == "CDS + Regulatory" ~ DAG * 1.025,
                                TRUE ~ DAG)) %>% setDT()

p2 <- ggplot(data=chr_fit_1Mb, aes(x=DAG_offset, y=-log10(pvalue), color=constrained_model)) +
  geom_point(size=3, shape=2) + theme_classic() + 
  geom_hline(yintercept=-log10(0.05/num_chrom), linetype="dashed", color="darkgrey") +
  scale_y_continuous(breaks=seq(0, 15, by=5)) +
  scale_x_continuous(breaks=c(5, 6)) +
  scale_color_manual(values=c("orange2", "magenta1"), name="DAG") +
  labs(title=NULL, x=NULL, y=expression(-log[10](italic(p)))) +
  theme(strip.text.y=element_text(size=16),       
        strip.text.x=element_text(size=16),       
        axis.title=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.text.x=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="none")

# fake plot for legend
fake <- ggplot(chr_fit_1Mb, aes(x=1, y=1, color=constrained_model)) +
  geom_point(size = 4) +   
  scale_color_manual(values = c("orange2", "magenta1"), name = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 16),   
        legend.title = element_text(size = 16))
legend <- cowplot::get_legend(fake)

q2 <- ggplot(data=chr_fit_1Mb, aes(x=DAG_offset, y=cfi, color=constrained_model, shape=cfi==0.49)) +
  geom_point(size=3, shape=2) + theme_classic() + 
  geom_hline(yintercept=0.9, linetype="dashed", color="darkgrey") +
  scale_y_continuous(breaks=pretty_breaks(), limits=c(0.49, 1)) +
  scale_x_continuous(breaks=c(5, 6)) +
  scale_color_manual(values=c("orange2", "magenta1"), name="DAG") +
  labs(title=NULL, x="DAG", y="CFI") +
  theme(strip.text.x=element_blank(),       
        axis.title=element_text(size=18), 
        axis.text=element_text(size=16), 
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="none")
y <- plot_grid(plot_grid(p2, q2, ncol=1, align="v", labels="AUTO", rel_heights=c(1, 1.175)), legend, nrow=2, rel_heights=c(1, 0.05))
save_plot("sem/chr_fit_compressed_DAGs5-6.pdf", y, base_height=6, base_width=7)

z <- plot_grid(x, y, ncol=1, align="v", rel_heights=c(0.8, 1))
save_plot("fivesix.pdf", z, base_height=9, base_width=9)
