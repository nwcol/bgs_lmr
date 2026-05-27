
source("~/Data/bgs_lmr/human_data/tools/HelperFunctions.R")

constrained_models <- c("split_cds_phastcons", "split_cds_regulatory")
mutation_maps <- c("roulette")

#################
#
# Fit measures & model selection
#
#################

# (Nested) model selection, genome-wide models, free or equal-chr regressions (LRT)
## 10 kb
LRT_gw_freechrom_10kb <- data.frame()
for(u in mutation_maps) {
  for(m in constrained_models) { 
    for(d in 1:4) {
      
      path <- paste0("sem/lavaan_out/simpi_LRT_DAG", d, "_gw-equal-freechrom_10kb_", m, "_", u, ".csv")
      
      tmp <- if(file.exists(path)) {
        fread(path)
      } else {
        NULL
      } 
      
      df <- as.data.frame(tmp$`Pr(>Chisq)`[2],)  
      df$constrained_model <- m
      df$dag <- d
      
      LRT_gw_freechrom_10kb <- rbind.data.frame(LRT_gw_freechrom_10kb, df)
    }
  }
}
LRT_gw_freechrom_10kb$scale <- 1e+4

## 100 kb
LRT_gw_freechrom_100kb <- data.frame()
for(u in mutation_maps) {
  for(m in constrained_models) { 
    for(d in 1:4) { 
      
      path <- paste0("sem/lavaan_out/simpi_LRT_DAG", d, "_gw-equal-freechrom_100kb_", m, "_", u, ".csv")
      
      tmp <- if(file.exists(path)) {
        fread(path)
      } else {
        NULL
      } 
      
      df <- as.data.frame(tmp$`Pr(>Chisq)`[2],)  
      df$constrained_model <- m
      df$dag <- d
      
      LRT_gw_freechrom_100kb <- rbind.data.frame(LRT_gw_freechrom_100kb, df)
    }
  }
}
LRT_gw_freechrom_100kb$scale <- 1e+5

## 1 Mb
LRT_gw_freechrom_1Mb <- data.frame()
for(u in mutation_maps) {
  for(m in constrained_models) { 
    for(d in 2:4) { 
      
      path <- paste0("sem/lavaan_out/simpi_LRT_DAG", d, "_gw-equal-freechrom_1Mb_", m, "_", u, ".csv")
      
      tmp <- if(file.exists(path)) {
        fread(path)
      } else {
        NULL
      } 
      
      df <- as.data.frame(tmp$`Pr(>Chisq)`[2],)  
      df$constrained_model <- m
      df$dag <- d
      
      LRT_gw_freechrom_1Mb <- rbind.data.frame(LRT_gw_freechrom_1Mb, df)
    }
  }
}
LRT_gw_freechrom_1Mb$scale <- 1e+6

gw_LRT <- rbind.data.frame(LRT_gw_freechrom_10kb, LRT_gw_freechrom_100kb, LRT_gw_freechrom_1Mb)
names(gw_LRT)[1] <- "pvalue"

p <- ggplot(data=filter(gw_LRT, constrained_model!="functional_merged_exons"), 
            aes(x=scale/1e3,  y=-log10(pvalue), shape=as.factor(dag), color=constrained_model)) +
  geom_point(size=3) + theme_classic() + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="darkgrey") +
  scale_x_log10(breaks=c(10, 100, 1000)) +
  scale_shape_manual(values=c(0, 1, 2, 3), name="DAG") +
  scale_color_manual(values=c("purple3", "seagreen4", "brown1"),
                     labels = c("functional_split_exons" = "Genes", "phastcons" = "phastCons"),
                     name=NULL) +
  labs(title=NULL, x="Scale (kb)", y=expression(-log[10](italic(p)))) +
  theme(strip.text.y=element_text(size=16),       
        strip.text.x=element_text(size=16),      
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("sem/simpi_gw_LRT_freechrom.pdf", p, base_height=4, base_width=12)

# per chromosome SEM models
## 10 kb
chr_fitmeasures_10kb <- data.frame()
for(dag in 3:4) {
  for(u in mutation_maps) {
    for(m in constrained_models) { 
      chr_list <- list(length=22)
      for(c in 1:22) {
        tmp <- fread(paste("sem/lavaan_out/simpi_YRI_DAG", dag, 
                           "_fitmeasures_10kb_chr", c, "_", m, "_", u, ".csv", sep=""))
        
        tmp$chrom <- c
        tmp$DAG <- dag
        tmp$mut_map <- u
        tmp$constrained_model <- m
        
        chr_list[[c]] <- tmp
      }
      
      df <- do.call(rbind.data.frame, chr_list)
      chr_fitmeasures_10kb <- rbind.data.frame(chr_fitmeasures_10kb, df)
    }
  }
}
chr_fitmeasures_10kb$scale <- 1e+4

## 100 kb
chr_fitmeasures_100kb <- data.frame()
for(dag in 3:4) {
  for(u in mutation_maps) {
    for(m in constrained_models) { 
      chr_list <- list(length=22)
      for(c in 1:22) {
        tmp <- fread(paste("sem/lavaan_out/simpi_YRI_DAG", dag, 
                           "_fitmeasures_100kb_chr", c, "_", m, "_", u, ".csv", sep=""))
        
        tmp$chrom <- c
        tmp$DAG <- dag
        tmp$mut_map <- u
        tmp$constrained_model <- m
        
        chr_list[[c]] <- tmp
      }
      
      df <- do.call(rbind.data.frame, chr_list)
      chr_fitmeasures_100kb <- rbind.data.frame(chr_fitmeasures_100kb, df)
    }
  }
}
chr_fitmeasures_100kb$scale <- 1e+5

## 1 Mb
chr_fitmeasures_1Mb <- data.frame()
for(dag in 3:4) {
  for(u in mutation_maps) {
    for(m in constrained_models) { 
      chr_list <- list(length=22)
      for(c in 1:22) {
        
        fname <- file.path(paste0("sem/lavaan_out/simpi_YRI_DAG", dag,
                                  "_fitmeasures_1Mb_chr", c, "_", m, "_", u, ".csv"))
        
        if(file.exists(fname)) {
          
          tmp <- fread(paste("sem/lavaan_out/simpi_YRI_DAG", dag, 
                             "_fitmeasures_1Mb_chr", c, "_", m, "_", u, ".csv", sep=""))
          
          tmp$chrom <- c
          tmp$DAG <- dag
          tmp$mut_map <- u
          tmp$constrained_model <- m
          
          chr_list[[c]] <- tmp
        }
      }
      
      df <- do.call(rbind.data.frame, chr_list)
      chr_fitmeasures_1Mb <- rbind.data.frame(chr_fitmeasures_1Mb, df)
    }
  }
}
chr_fitmeasures_1Mb$scale <- 1e+6

chr_fitmeasures <- rbind.data.frame(chr_fitmeasures_10kb, chr_fitmeasures_100kb, chr_fitmeasures_1Mb)
chr_fitmeasures <- filter(chr_fitmeasures, constrained_model!="functional_merged_exons")

chr_fitmeasures$constrained_model <- factor(chr_fitmeasures$constrained_model, 
                                            levels = names(elem_labels), 
                                            labels = unname(elem_labels)) 

chr_fitmeasures$scale_kb <- chr_fitmeasures$scale / 1e3

chr_fitmeasures <- chr_fitmeasures %>% 
  mutate(scale_offset = case_when(DAG == 3 ~ scale * 0.9, 
                                  DAG == 4 ~ scale * 1.1,
                                  TRUE ~ scale))

fwrite(dplyr::select(chr_fitmeasures, c(pvalue, cfi, chrom, DAG, mut_map, constrained_model, scale)),
       "sem/simpi_chr_fitmeasures.csv")

for(mut in mutation_maps) {
  p1 <- ggplot(data=filter(chr_fitmeasures, mut_map==mut), 
               aes(x=scale_offset / 1e3,  y=-log10(pvalue), color=as.factor(DAG))) +
    geom_point(size=2, shape=2) + 
    facet_wrap(~constrained_model) + theme_classic() + 
    geom_hline(yintercept=-log10(0.05 / 22), linetype="dashed", color="darkgrey") +
    scale_y_continuous(breaks=pretty_breaks()) +
    scale_x_log10(breaks=c(10, 100, 1000)) +
    scale_color_manual(values=c("brown1", "cyan3"), name="DAG") +
    labs(title=NULL, x=NULL, y=expression(-log[10](italic(p)))) +
    theme(strip.text.y=element_text(size=16),     
          strip.text.x=element_text(size=16),       
          axis.title=element_text(size=16), 
          axis.text=element_text(size=12), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16),
          legend.position="none")
  
  p2 <- ggplot(data=filter(chr_fitmeasures, mut_map==mut), 
               aes(x=scale_offset / 1e3, y=cfi, color=as.factor(DAG))) +
    geom_point(size=2, shape=2) + 
    facet_wrap(~constrained_model) + theme_classic() + 
    geom_hline(yintercept=0.9, linetype="dashed", color="darkgrey") +
    scale_y_continuous(breaks=pretty_breaks()) +
    scale_x_log10(breaks=c(10, 100, 1000)) +
    scale_color_manual(values=c("brown1", "cyan3"), name="DAG") +
    labs(title=NULL, x=NULL, y="CFI") +
    theme(#strip.text.y=element_blank(),     
      strip.text.x=element_blank(),
      axis.title=element_text(size=16), 
      axis.text=element_text(size=12), 
      axis.text.x=element_text(size=12),
      legend.text=element_text(size=16),
      legend.title=element_text(size=16),
      legend.position="bottom")
  sim_gfit <- plot_grid(p1, p2, ncol=1, rel_heights=c(1, 1.25), align="v")
  save_plot(paste0("sem/simpi_fitmeasures_", mut, ".pdf"), sim_gfit, base_height=8, base_width=12)
}

# (Nested) model selection, DAGs 3 and 4 (LRT)
# (DAG 4 adds a covariance term between constrained sites and mutation rates)
## 10 kb
chr_LRT_10kb <- data.frame()
for(u in mutation_maps) {
  for(m in constrained_models) { 
    chr_list <- list(length=22)
    for(c in 1:22) {
      
      tmp <- fread(paste0("sem/lavaan_out/simpi_LRT_DAG3_DAG4_10kb_chr",
                          c, "_", m, "_", u, ".csv"))
      
      tmp <- as.data.frame(tmp$`Pr(>Chisq)`[2],)  
      tmp$chrom <- c
      tmp$mut_map <- u
      tmp$constrained_model <- m
      
      chr_list[[c]] <- tmp
    }
    
    df <- do.call(rbind.data.frame, chr_list)
    chr_LRT_10kb <- rbind.data.frame(chr_LRT_10kb, df)
  }
}
chr_LRT_10kb$scale <- 1e+4

## 100 kb
chr_LRT_100kb <- data.frame()
for(u in mutation_maps) {
  for(m in constrained_models) { 
    chr_list <- list(length=22)
    for(c in 1:22) {
      
      tmp <- fread(paste0("sem/lavaan_out/simpi_LRT_DAG3_DAG4_100kb_chr",
                          c, "_", m, "_", u, ".csv"))
      
      tmp <- as.data.frame(tmp$`Pr(>Chisq)`[2],)  
      tmp$chrom <- c
      tmp$mut_map <- u
      tmp$constrained_model <- m
      
      chr_list[[c]] <- tmp
    }
    
    df <- do.call(rbind.data.frame, chr_list)
    chr_LRT_100kb <- rbind.data.frame(chr_LRT_100kb, df)
  }
}
chr_LRT_100kb$scale <- 1e+5

## 1 Mb
chr_LRT_1Mb <- data.frame()
for(u in mutation_maps) {
  for(m in constrained_models) { 
    chr_list <- list(length=22)
    for(c in 1:22) {
      
      tmp <- fread(paste0("sem/lavaan_out/simpi_LRT_DAG3_DAG4_1Mb_chr",
                          c, "_", m, "_", u, ".csv"))
      
      tmp <- as.data.frame(tmp$`Pr(>Chisq)`[2],)  
      tmp$chrom <- c
      tmp$mut_map <- u
      tmp$constrained_model <- m
      
      chr_list[[c]] <- tmp
    }
    
    df <- do.call(rbind.data.frame, chr_list)
    chr_LRT_1Mb <- rbind.data.frame(chr_LRT_1Mb, df)
  }
}
chr_LRT_1Mb$scale <- 1e+6

chr_LRT <- rbind.data.frame(chr_LRT_10kb, chr_LRT_100kb, chr_LRT_1Mb)
names(chr_LRT)[1] <- "pvalue"

chr_LRT <- filter(chr_LRT, constrained_model!="functional_split_exons")

chr_LRT$constrained_model <- factor(chr_LRT$constrained_model, 
                                    levels = names(elem_labels), 
                                    labels = unname(elem_labels)) 

chr_LRT <- chr_LRT %>% 
  mutate(scale_offset = case_when(mut_map == "gnomad" ~ scale * 0.85, 
                                  mut_map == "carlson" ~ scale * 1, 
                                  mut_map == "roulette" ~ scale * 1.2, 
                                  TRUE ~ scale))

p <- ggplot(data=chr_LRT,
            aes(x=scale_offset/1e3, 
                y=-log10(pvalue), 
                color=mut_map)) +
  geom_point(size=3, shape=2) +
  geom_hline(yintercept=-log10(0.05 / 22), linetype="dashed", color="darkgrey") +
  facet_grid(~constrained_model, scale="free_y") + theme_classic() + 
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_x_log10(breaks=c(10, 100, 1000)) +
  scale_color_manual(values=c("purple3", "seagreen4", "brown1"), name="Mut. Map") +
  labs(title=NULL, x="Scale (kb)", y=expression(-log[10](italic(p)))) +
  theme(strip.text.y=element_text(size=16),      
        strip.text.x=element_text(size=16),       
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("sem/simpi_chr_LRT_DAGs3vs4.pdf", p, base_height=5, base_width=12)

#################
#
# Parameter estimates
#
#################

# per chromosome SEM models
## 10 kb
chr_params_10kb <- data.frame()
for(dag in 3:4) {
  for(u in mutation_maps) {
    for(m in constrained_models) { 
      chr_list <- list(length=22)
      for(c in 1:22) {
        
        tmp <- fread(paste0("sem/lavaan_out/simpi_YRI_DAG", dag, "_params_10kb_chr", 
                            c, "_", m, "_", u, ".csv"))
        
        tmp$chrom <- c
        tmp$DAG <- dag
        tmp$mut_map <- u
        tmp$constrained_model <- m
        
        chr_list[[c]] <- tmp
      }
      
      df <- do.call(rbind.data.frame, chr_list)
      chr_params_10kb <- rbind.data.frame(chr_params_10kb, df)
    }
  }
}
chr_params_10kb$scale <- 1e+4

## 100 kb
chr_params_100kb <- data.frame()
for(dag in 3:4) {
  for(u in mutation_maps) {
    for(m in constrained_models) { 
      chr_list <- list(length=22)
      for(c in 1:22) {
        
        tmp <- fread(paste0("sem/lavaan_out/simpi_YRI_DAG", dag, "_params_100kb_chr", 
                            c, "_", m, "_", u, ".csv"))
        
        tmp$chrom <- c
        tmp$DAG <- dag
        tmp$mut_map <- u
        tmp$constrained_model <- m
        
        chr_list[[c]] <- tmp
      }
      
      df <- do.call(rbind.data.frame, chr_list)
      chr_params_100kb <- rbind.data.frame(chr_params_100kb, df)
    }
  }
}
chr_params_100kb$scale <- 1e+5

## 1 Mb
chr_params_1Mb <- data.frame()
for(dag in 3:4) {
  for(u in mutation_maps) {
    for(m in constrained_models) { 
      chr_list <- list(length=22)
      for(c in 1:22) {
        
        tmp <- fread(paste0("sem/lavaan_out/simpi_YRI_DAG", dag, "_params_1Mb_chr", 
                            c, "_", m, "_", u, ".csv"))
        
        tmp$chrom <- c
        tmp$DAG <- dag
        tmp$mut_map <- u
        tmp$constrained_model <- m
        
        chr_list[[c]] <- tmp
      }
      
      df <- do.call(rbind.data.frame, chr_list)
      chr_params_1Mb <- rbind.data.frame(chr_params_1Mb, df)
    }
  }
}
chr_params_1Mb$scale <- 1e+6

chr_params <- rbind.data.frame(chr_params_10kb, chr_params_100kb, chr_params_1Mb)  %>%
  filter(., constrained_model!="functional_merged_exons") %>%
  filter(., op=="~" | op==":=" | (lhs=="avg_mut" & op=="~~" & rhs=="prop_del_sites")) 

chr_params$sig <- chr_params$pvalue < 0.05 / 22
chr_params$constrained_model <- factor(chr_params$constrained_model, 
                                       levels = names(elem_labels), 
                                       labels = unname(elem_labels)) 

chr_params$scale_kb <- chr_params$scale / 1e3
chr_params[op=="~~",]$rhs <- "covar"
chr_params$is_19 <- chr_params$chrom==19

fwrite(dplyr::select(chr_params, -scale_kb), "sem/simpi_chr_PC_DAGs.csv")

for(u in mutation_maps) {
  p <- ggplot(data=filter(chr_params,
                          DAG==4,
                          mut_map==u,
                          lhs=="pi_sim_tot" | lhs=="rec_B_pi" | lhs=="rec_mu_pi" | op=="~~"),
              aes(x=rhs, y=est, shape=sig, color=is_19, group=constrained_model)) +
    geom_point(size=3, position = position_dodge(width = 0.3)) + 
    facet_wrap(~scale_kb, ncol=1) + theme_classic() + 
    scale_color_manual(values=c("cyan3", "brown1"), name="Chromosome",
                       labels = c("FALSE" = "Others", "TRUE" = "19")) +
    scale_shape_manual(values=c(8, 2), name=NULL,
                       labels = c("FALSE" = expression(p >= 0.05/22), "TRUE" = expression(p < 0.05/22))) +
    scale_y_continuous(breaks=pretty_breaks()) +
    #scale_fill_manual(values = c("Genes" = "green", "phastCons" = "grey80"), name = "Model") +
    scale_x_discrete(labels = c( "avg_mut" = expression(mu %->% pi),
                                 "c1*c7" = expression(r %->% mu %->% pi),
                                 "c2*c6" = expression(r %->% B %->% pi),
                                 "B" = expression(B %->% pi),
                                 "prop_del_sites" = expression(E %->% pi),
                                 "covar" = expression(mu %~% E))) +
    geom_hline(yintercept=0, color="black", linetype="dashed") +
    labs(title=NULL, x="Path", y="Coefficient") +
    theme(strip.text.y=element_text(size=16),       
          strip.text.x=element_text(size=16),       
          axis.title=element_text(size=16), 
          axis.text=element_text(size=14), 
          legend.text=element_text(size=16),
          legend.title=element_text(size=16),
          legend.position="bottom")
  save_plot(paste0("sem/simpi_chr_params_pi_DAG4_", u, ".pdf"), p, base_height=10, base_width=10)
}

chr_params_ext <- left_join(chr_params, 
                            chr_fitmeasures,
                            by=c("chrom", "DAG", "mut_map",
                                 "constrained_model", "scale"))

names(chr_params_ext)[8] <- "pvalue_est"
names(chr_params_ext)[23] <- "pvalue_DAG"

p1 <- ggplot(data=filter(chr_params_ext,
                         lhs=="pi_sim_tot", 
                         DAG==3), 
             aes(x=chrom, 
                 y=est, 
                 color=mut_map,
                 shape=rhs)) +
  geom_point(size=3) + geom_hline(yintercept=0, linetype="dashed") +
  facet_grid(scale/1e3~constrained_model) + theme_classic() + 
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_x_continuous(breaks=1:22) + 
  scale_shape_manual(values=c(0, 1, 2), name="Variable") +
  scale_color_manual(values=c("purple3", "seagreen4", "brown1"), name="Mut. Map") +
  labs(title=NULL, x="Chromosome", y="Direct Std. PC to Diversity") +
  theme(strip.text.y=element_text(size=16),     
        strip.text.x=element_text(size=16),     
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("sem/simpi_chr_params_pi_DAG3.pdf", p1, base_height=8, base_width=18)


p2 <- ggplot(data=filter(chr_params_ext,
                         lhs=="pi_sim_tot", 
                         DAG==4), 
             aes(x=chrom, 
                 y=est, 
                 color=mut_map,
                 shape=rhs)) +
  geom_point(size=3) + geom_hline(yintercept=0, linetype="dashed") +
  facet_grid(scale/1e3~constrained_model) + theme_classic() + 
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_x_continuous(breaks=1:22) + 
  scale_shape_manual(values=c(0, 1, 2), name="Variable") +
  scale_color_manual(values=c("purple3", "seagreen4", "brown1"), name="Mut. Map") +
  labs(title=NULL, x="Chromosome", y="Direct Std. PC to Diversity") +
  theme(strip.text.y=element_text(size=16),     
        strip.text.x=element_text(size=16),     
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("sem/simpi_chr_params_pi_DAG4.pdf", p2, base_height=8, base_width=18)

