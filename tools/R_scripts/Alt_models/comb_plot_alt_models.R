library(tidyverse)
library(data.table)

cousins <- fread("~/Data/bgs_lmr/PNAS_reviews/cousins_demo/summary_SEM.csv")
singleNe <- fread("~/Data/bgs_lmr/human_data/current_models/equilibrium/summary_SEM.csv")
decode <- fread("~/Data/bgs_lmr/PNAS_reviews/deCODE/summary_SEM.csv")

cousins$model <- "Cousins"
singleNe$model <- "Single-Ne"
decode$model <- "deCODE"

alt_models <- rbind.data.frame(cousins, singleNe, decode)

alt_models$variable <- factor(alt_models$variable,
                              levels=c("u->pi", "B->pi", "E->pi", "r->B->pi", "r->u->pi", "R^2 * '(' * pi * ')'"))

my_labeller <- as_labeller(c(
  "u->pi" = "u %->% pi",
  "B->pi" = "B %->% pi",
  "E->pi" = "E %->% pi",
  "r->B->pi" = "r %->% B %->% pi",
  "r->u->pi" = "r %->% u %->% pi",
  "R^2 * '(' * pi * ')'" = "R^2 * '(' * pi * ')'"
), default = label_parsed)

x1 <- ggplot(data=filter(alt_models, color_var != "metric", DAG == 4),
  aes(x=constrained_model, y=value, color=color_var, group=model, shape=model)) +
  geom_point(size=3, position=position_dodge(width=0.6)) + theme_bw() +
  facet_grid(scale/1e3 ~ variable, labeller=labeller(variable=my_labeller)) + 
  scale_shape_manual(name=NULL, values=c(0, 1, 2)) +
  scale_color_manual(name=NULL,
    values=c("sigTRUE"="brown1", "sigFALSE"="cyan3"),
    labels=c("sigTRUE"="Significant", "sigFALSE"="Non-significant"),
    breaks=c("sigTRUE", "sigFALSE")) +
  geom_hline(yintercept=0, linetype="dashed") +
  guides(color=guide_legend(override.aes=list(shape=18))) +
  scale_y_continuous(breaks=pretty_breaks()) +
  labs(title=NULL, x="Constraint model", y="Estimate") +
  theme(strip.text.y=element_text(size=16),
        strip.text.x=element_text(size=16),
        axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("~/Data/bgs_lmr/PNAS_reviews/alt_models_est.pdf", x1, base_height=10, base_width=17)

x2 <- ggplot(data=filter(alt_models, color_var=="metric", DAG==4),
             aes(x=constrained_model, y=value, group=model, shape=model)) +
  geom_point(size=3, position=position_dodge(width=0.6)) + theme_bw() +
  facet_grid(scale/1e3~variable, labeller=label_parsed) + 
  scale_shape_manual(name=NULL, values=c(0, 1, 2)) +
  scale_y_continuous(breaks=pretty_breaks()) + 
  guides(color=guide_legend(override.aes=list(shape=18))) +
  labs(title=NULL, x="Constraint model", y="Estimate") +
  theme(strip.text.y=element_text(size=16),    
        strip.text.x=element_text(size=16),       
        axis.title=element_text(size=16), 
        axis.text=element_text(size=12), 
        axis.text.x=element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("~/Data/bgs_lmr/PNAS_reviews/alt_models_fit.pdf", x2, base_height=8, base_width=12)

