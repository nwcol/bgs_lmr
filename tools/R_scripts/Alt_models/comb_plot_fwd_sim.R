# set 1: 100x 10-Mb chr
demog_set1 <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/10Mb/demog_predictions/summary_SEM.csv")
fit_Ne_good_set1 <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/10Mb/fit_Ne_predictions_good_DFE/summary_SEM.csv")
fit_Ne_bad_set1 <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/10Mb/fit_Ne_predictions_bad_DFE/summary_SEM.csv")

demog_set1[, model := "True model"]
fit_Ne_good_set1[, model := "Eq. good-DFE"]
fit_Ne_bad_set1[, model := "Eq. bad-DFE"]

fwd_sim_set1 <- rbind.data.frame(demog_set1, fit_Ne_good_set1, fit_Ne_bad_set1) %>%
                dplyr::filter(., scale != 1e+6) %>% setDT()
fwd_sim_set1[, sim := "100x 10-Mb chr"]

# set 2: 10x 100-Mb chr
demog_set2 <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/100Mb/demog_predictions/summary_SEM.csv")
fit_Ne_good_set2 <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/100Mb/fit_Ne_predictions_good_DFE/summary_SEM.csv")
fit_Ne_bad_set2 <- fread("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/100Mb/fit_Ne_predictions_bad_DFE/summary_SEM.csv")

demog_set2[, model := "True model"]
fit_Ne_good_set2[, model := "Eq. good-DFE"]
fit_Ne_bad_set2[, model := "Eq. bad-DFE"]

fwd_sim_set2 <- rbind.data.frame(demog_set2, fit_Ne_good_set2, fit_Ne_bad_set2) %>% setDT()
fwd_sim_set2[, sim := "10x 100-Mb chr"]

fwd_sim <- rbind.data.frame(fwd_sim_set1, fwd_sim_set2) %>% setDT()

x1 <- ggplot(data=filter(fwd_sim, color_var!="metric", DAG==4), aes(x=model, y=value, color=color_var, group=sim, shape=sim)) +
  geom_point(size=3, position=position_dodge(width=0.6)) + theme_bw() +
  facet_grid(scale/1e3 ~ variable, labeller=label_parsed) +
  scale_shape_manual(name=NULL, values=c(1, 2)) +
  scale_color_manual(name=NULL,
                     values=c("sigTRUE"="brown1", "sigFALSE"="cyan3"),
                     labels=c("sigTRUE"="Significant", "sigFALSE"="Non-significant"),
                     breaks=c("sigTRUE","sigFALSE")) +
  geom_hline(yintercept=0, linetype="dashed") +
  guides(color=guide_legend(override.aes=list(shape=16))) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=NULL, x="B-value parameters", y="Estimate") +
  theme(strip.text.y=element_text(size=16),
        strip.text.x=element_text(size=16),
        axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/fwd_sim_est.pdf", x1, base_height=8, base_width=12)

x2 <- ggplot(data=filter(fwd_sim, color_var=="metric", DAG==4), aes(x=model, y=value, group=sim, shape=sim)) +
  geom_point(size=3, position=position_dodge(width=0.6)) + theme_bw() +
  facet_grid(scale/1e3~variable, labeller=label_parsed) + 
  scale_shape_manual(name=NULL, values=c(1, 2)) +
  scale_y_continuous(breaks=pretty_breaks()) + 
  labs(title=NULL, x="B-value parameters", y="Estimate") +
  theme(strip.text.y=element_text(size=16),    
        strip.text.x=element_text(size=16),       
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("~/Data/bgs_lmr/PNAS_reviews/fwd_sims/fwd_sims_to_share_pnas/fwd_sim_fit.pdf", x2, base_height=8, base_width=12)


# LRT, set 1
LRT_DAG3vs4_set1_demog <- fread("10Mb/demog_predictions/sem/chr_LRT_DAGs3-4.csv")
LRT_DAG3vs4_set1_fitbad <- fread("10Mb/fit_Ne_predictions_bad_DFE/sem/chr_LRT_DAGs3-4.csv")
LRT_DAG3vs4_set1_fitgood <- fread("10Mb/fit_Ne_predictions_good_DFE/sem/chr_LRT_DAGs3-4.csv")

LRT_DAG3vs4_set1_demog[, model := "True model"]
LRT_DAG3vs4_set1_fitgood[, model := "Eq. good-DFE"]
LRT_DAG3vs4_set1_fitbad[, model := "Eq. bad-DFE"]

fwd_sim_set1 <- rbind.data.frame(LRT_DAG3vs4_set1_demog, LRT_DAG3vs4_set1_fitgood, LRT_DAG3vs4_set1_fitbad) %>% setDT()
fwd_sim_set1[, sim := "100x 10-Mb chr"]

# LRT, set 2
LRT_DAG3vs4_set2_demog <- fread("100Mb/demog_predictions/sem/chr_LRT_DAGs3-4.csv")
LRT_DAG3vs4_set2_fitbad <- fread("100Mb/fit_Ne_predictions_bad_DFE/sem/chr_LRT_DAGs3-4.csv")
LRT_DAG3vs4_set2_fitgood <- fread("100Mb/fit_Ne_predictions_good_DFE/sem/chr_LRT_DAGs3-4.csv")

LRT_DAG3vs4_set2_demog[, model := "True model"]
LRT_DAG3vs4_set2_fitgood[, model := "Eq. good-DFE"]
LRT_DAG3vs4_set2_fitbad[, model := "Eq. bad-DFE"]

fwd_sim_set2 <- rbind.data.frame(LRT_DAG3vs4_set2_demog, LRT_DAG3vs4_set2_fitgood, LRT_DAG3vs4_set2_fitbad) %>% setDT()
fwd_sim_set2[, sim := "10x 100-Mb chr"]

fwd_sim <- rbind.data.frame(fwd_sim_set1, fwd_sim_set2) %>% setDT()

x3 <- ggplot(data=fwd_sim, aes(x=model, y=-log10(pvalue), group=sim, shape=sim)) +
  geom_point(size=3, position=position_dodge(width=0.6)) + theme_bw() +
  facet_wrap(~scale/1e3) + 
  geom_hline(yintercept=-log10(0.05/110), linetype="dashed") +
  scale_shape_manual(name=NULL, values=c(1, 2)) +
  scale_y_continuous(breaks=pretty_breaks()) + 
  labs(title=NULL, x="Model", y=expression(-log[10](italic(p)))) +
  theme(strip.text.y=element_text(size=16),    
        strip.text.x=element_text(size=16),       
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("fwdsim_LRT_DAG3vsDAG4.pdf", x3, base_height=5, base_width=12)

