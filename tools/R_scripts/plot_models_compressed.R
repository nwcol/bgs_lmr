
chr_params_ext <- fread("sem/obspi_chr_PC_DAGs.csv")
chr_r2_ext_DAGs34 <- fread("sem/chr_r2_ext_DAGs34.csv")

# compressing plot
PC_pi <- filter(chr_params_ext, lhs=="avg_pi") %>% 
  dplyr::select(., c("rhs", "est", "chrom", "DAG", "constrained_model", "scale", "sig")) %>% setDT()

setnames(PC_pi, old="rhs", new="Var")
PC_pi <- pivot_wider(PC_pi, names_from=Var, values_from=c(est, sig)) %>% setDT()

chr_r2_pi <- filter(chr_r2_ext_DAGs34) %>%
  dplyr::select(., c("avg_pi", "chrom", "DAG", "scale", "pvalue", "cfi")) %>% setDT()

setnames(chr_r2_pi, old="avg_pi", new="R2_pi")

pc_dt <- PC_pi[chr_r2_pi, on = .(chrom, DAG, scale)]
pc_dt_m <- pivot_longer(pc_dt, cols=c(starts_with("est"), "R2_pi", "pvalue", "cfi"), names_to="variable") %>% setDT()
# map significance to the correct estimates
pc_dt_m[, sig := fifelse(variable == "est_avg_mut",  sig_avg_mut,
                         fifelse(variable == "est_B", sig_B,
                                 fifelse(variable == "est_prop_del_sites", sig_prop_del_sites, NA)))]
# CFI and R2 don't have p-values
pc_dt_m[, variable := factor(variable, levels=c("est_avg_mut", "est_B", "est_prop_del_sites", "R2_pi", "pvalue", "cfi"),
                             labels=c("u", "B", "E", "R^2 * '(' * pi * ')'", "p-value", "CFI"))]
pc_dt_m[, color_var := fifelse(variable %in% c("u","B","E"), ifelse(sig, "sigTRUE", "sigFALSE"), "metric")]
pc_dt_m <- dplyr::select(pc_dt_m, c("chrom", "DAG", "constrained_model", "scale", "variable", "value", "color_var"))
fwrite(pc_dt_m, "summary_SEM.csv")

x1 <- ggplot(data=filter(pc_dt_m, color_var!="metric", DAG==4), aes(x=constrained_model, y=value, color=color_var)) +
  geom_point(size=3, shape=1) + facet_grid(scale/1e3~variable, labeller=label_parsed) + 
  scale_color_manual(name=NULL, values=c("sigTRUE"="brown1", "sigFALSE"="cyan3"),
                     labels = c("sigTRUE"="Significant", "sigFALSE"="Non-significant"),
                     breaks=c("sigTRUE", "sigFALSE")) +
  geom_hline(yintercept=0, linetype="dashed") +
  guides(color=guide_legend(override.aes=list(shape=16))) +
  scale_y_continuous(breaks=pretty_breaks()) + theme_bw() + 
  labs(title=NULL, x="Constraint model", y="Estimate") +
  theme(strip.text.y=element_text(size=16),    
        strip.text.x=element_text(size=16),       
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("sem/PC_compressed.pdf", x1, base_height=8, base_width=10)

x2 <- ggplot(data=filter(pc_dt_m, color_var=="metric", DAG==4), aes(x=constrained_model, y=value)) +
  geom_point(size=3, shape=1) + facet_grid(scale/1e3~variable, labeller=label_parsed) + 
  guides(color=guide_legend(override.aes=list(shape=16))) +
  scale_y_continuous(breaks=pretty_breaks()) + theme_bw() + 
  labs(title=NULL, x="Constraint model", y="Estimate") +
  theme(strip.text.y=element_text(size=16),    
        strip.text.x=element_text(size=16),       
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("sem/fit_compressed.pdf", x2, base_height=8, base_width=8)

