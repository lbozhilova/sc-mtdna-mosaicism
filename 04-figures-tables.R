###################
##### Figures #####
###################

# Last edited: 11/04/23 by LVB

# Description: Manuscript and SI figures and tables. 

#----- Packages
library("tidyverse")
library("ggbeeswarm")
library("kimura")
library("reshape2")

#----- Code
source("00-plot-setup.R")

#----- Data
load("data/parsed/01-parsed-data.RData")
load("data/parsed/01-parsed-replicates.RData")
load("data/parsed/01-parsed-copy-numbers.RData")
load("data/parsed/02-kimura-fits.RData")
load("data/parsed/03-models.RData")

#----- General functionality
# Tissue labels
get_tissue_labs <- function(x) {
  x <- x %>% unique %>% as.character
  x[order(match(x, tissue_names))]
}

#----- Main: Figure 1
# (A) Pipeline
fig_1a <- ggplot() + 
  theme_lvb + theme(panel.border = element_blank(),
                    axis.line = element_blank(),
                    axis.line.y = element_blank()) +
  draw_image("figures/protocol-horizontal.tiff", width = 3.25, height = 1, halign = 0)
fig_1a

# (B) Bulk tissue P100 m.5024C>T
df_1b <- bt_df %>% filter(day == "P100" & mutation == "m.5024")
ts_1b <- get_tissue_labs(df_1b$tissue)

fig_1b <- ggplot(df_1b, aes(x = tissue, y = het_bulk, shape = tissue, colour = tissue)) +
  theme_lvb + theme(legend.position = "bottom", panel.grid.major.x = element_blank(), legend.justification = 0.2) + 
  facet_grid(mouse~., labeller = master_lab) +
  geom_point(size = 2) +
  scale_x_discrete("", labels = NULL) +
  scale_y_continuous("Heteroplasmy %", limits = c(0, 100)) +
  scale_shape_manual("", values = tissue_shapes[ts_1b], labels = tissue_labs[ts_1b]) +
  scale_color_manual("", values = tissue_cols[ts_1b], labels = tissue_labs[ts_1b]) +
  guides(shape = guide_legend(ncol = 4))
fig_1b

# (C) Single cells P100 m.5024C>T
df_1c <- sc_df %>% filter(day == "P100" & mutation == "m.5024")

fig_1c <- ggplot(df_1c, aes(x = cell_type, y = het, shape = cell_type, colour = tissue)) +
  theme_lvb + theme(legend.position = "none") +
  facet_grid(mouse ~ tissue, labeller = master_lab, 
             scales = "free", switch = "x") +
  geom_hline(aes(yintercept = het_bulk)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = .1) +
  geom_quasirandom(alpha = 0.5) +
  scale_x_discrete("", labels = cell_labs) +
  scale_y_continuous("Heteroplasmy %", limits = c(0, 100), labels = seq(0, 100, 25)) +
  scale_shape_manual(values = cell_shapes, labels = cell_labs) +
  scale_colour_manual(values = tissue_cols)
fig_1c

# Compound figure
fig_1bc <- plot_grid(fig_1b, fig_1c,
                     labels = c("B", "C"), rel_widths = c(1, 2.25), 
                     ncol = 2, label_size = 10)
fig_1bc
fig1 <- plot_grid(fig_1a, fig_1bc,
                  labels = c("A", ""), rel_heights = c(1, 3),
                  ncol = 1, label_size = 10)
fig1

plot_save(fig1, "figures/04-fig1.jpg", size = 1, ar = 3.25/4)

#----- Main: Figure 2
# (A) Bulk tissue P0, P6, P365 m.5024C>T
df_2a <- bt_df %>% filter(day %in% c("P0", "P6", "P365") & mutation == "m.5024") 
ts_2a <- get_tissue_labs(df_2a$tissue)

fig_2a <- ggplot(df_2a, aes(x = tissue, y = het_bulk, shape = tissue, colour = tissue)) +
  theme_lvb + 
  theme(legend.position = "bottom", panel.grid.major.x = element_blank(), legend.justification = 0.2) +
  facet_grid(day ~ mouse, labeller = master_lab) +
  geom_point(size = 2) +
  scale_y_continuous("Heteroplasmy %", limits = c(0, 100)) +
  scale_x_discrete("", labels = NULL) +
  scale_shape_manual(values = tissue_shapes[ts_2a], labels = tissue_labs[ts_2a]) +
  scale_colour_manual(values = tissue_cols[ts_2a], labels = tissue_labs[ts_2a]) +
  guides(shape = guide_legend(ncol = 3))
fig_2a

# (B) Single cells P0, P6, P365 m.5024C>T
df_2b <- sc_df %>% filter(day %in% c("P0", "P6", "P365") & mutation == "m.5024")

fig_2b <- ggplot(df_2b, aes(x = tissue, y = het)) +
  theme_lvb + theme(legend.position = "none") +
  facet_grid(day ~ mouse, labeller = master_lab) +
  geom_violin(aes(fill = tissue)) +
  geom_boxplot(aes(fill = tissue), width = 0.1) +
  scale_y_continuous("Heteroplasmy %", limits = c(0, 105)) +
  scale_x_discrete("", labels = tissue_labs) +
  scale_fill_manual(values = tissue_cols)
fig_2b

# (C) Single cells E8.5 m.5024C>T
df_2c <- sc_df %>% filter(day  == "E8.5" & mutation == "m.5024")

fig_2c <- ggplot(df_2c, aes(x = mouse, y = het)) +
  theme_lvb + theme(legend.position = "none") +
  geom_violin(fill = mutation_cols["m.5024"]) +
  geom_boxplot(width = 0.1, fill = mutation_cols["m5024"]) +
  scale_y_continuous("Heteroplasmy %", limits = c(0, 105)) +
  scale_x_discrete("", labels = mouse_labs) 
fig_2c

# Compound figure
fig_2ac <- plot_grid(fig_2a, fig_2c,
                     labels = c("A", "C"), rel_heights = c(2, 1.1), 
                     ncol = 1, label_size = 10)
fig_2ac
fig2 <- plot_grid(fig_2ac, fig_2b, 
                  labels = c("", "B"), rel_widths = c(1, 1.5),
                  label_size = 10)
fig2

plot_save(fig2, "figures/04-fig2.jpg", size = 1, ar = 1/1.1)

#----- Main: Figure 3
# (A) Bulk tissue P100 m.5019A>G
df_3a <- bt_df %>% filter(day == "P100" & mutation == "m.5019")
ts_3a <- get_tissue_labs(df_3a$tissue)

fig_3a <- ggplot(df_3a, aes(x = tissue, y = het_bulk, shape = tissue, colour = tissue)) +
  theme_lvb + 
  theme(legend.position = "bottom", panel.grid.major.x = element_blank(), legend.justification = 0.2) +
  facet_grid(mouse~., labeller = master_lab) +
  geom_point(size = 2) +
  scale_x_discrete("", labels = NULL) +
  scale_y_continuous("Heteroplasmy %", limits = c(0, 100)) +
  scale_shape_manual("", values = tissue_shapes[ts_3a], labels = tissue_labs[ts_3a]) +
  scale_color_manual("", values = tissue_cols[ts_3a], labels = tissue_labs[ts_3a]) +
  guides(shape = guide_legend(ncol = 4))
fig_3a

# (B) Single cells P100 m.5019A>G
df_3b <- sc_df %>% filter(day == "P100" & mutation == "m.5019")

fig_3b <- ggplot(df_3b, aes(x = cell_type, y = het, shape = cell_type, colour = tissue)) +
  theme_lvb + theme(legend.position = "none") +
  facet_grid(mouse ~ tissue, labeller = master_lab, scales = "free", switch = "x") +
  geom_hline(aes(yintercept = het_bulk)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = .1) +
  geom_quasirandom(alpha = 0.5) +
  scale_x_discrete("", labels = cell_labs) +
  scale_y_continuous("Heteroplasmy %", limits = c(0, 100)) +
  scale_shape_manual(values = cell_shapes, labels = cell_labs) +
  scale_colour_manual(values = tissue_cols)
fig_3b

# Compound figure
fig_3 <- plot_grid(fig_3a, fig_3b, labels = c("A", "B"), rel_widths = c(1, 2), label_size = 10)
plot_save(fig_3, "figures/04-fig3.jpg", size = 1, ar = 4/3)

#----- Main: Figure 4
# (A) Bulk tissue  P0, P6, P365 m.5019A>G
df_4a <- bt_df %>% 
  filter(day %in% c("P0", "P6", "P365")) %>% 
  filter(mutation == "m.5019" & mouse != "mouse1.1")
ts_4a <- get_tissue_labs(df_4a$tissue)

fig_4a <- ggplot(df_4a, aes(x = tissue, y = het_bulk, shape = tissue, colour = tissue)) +
  theme_lvb + 
  theme(legend.position = "bottom", panel.grid.major.x = element_blank(), legend.justification = 0.2) +
  facet_grid(day ~ mouse, labeller = master_lab) +
  geom_point(size = 2) +
  scale_shape_manual(values = tissue_shapes[ts_4a], labels = tissue_labs[ts_4a]) +
  scale_y_continuous("Heteroplasmy %", limits = c(0, 100)) +
  scale_x_discrete("", labels = NULL) +
  scale_colour_manual(values = tissue_cols[ts_4a], labels = tissue_labs[ts_4a]) +
  guides(shape = guide_legend(ncol = 3))
fig_4a

# (B) Single cells P0, P6, P365 m.5019A>G
df_4b <- sc_df %>% filter(day %in% c("P0", "P6", "P365") & mutation == "m.5019")

fig_4b <- ggplot(df_4b, aes(x = tissue, y = het)) +
  theme_lvb + theme(legend.position = "none") +
  facet_grid(day ~ mouse, labeller = master_lab) +
  geom_violin(aes(fill = tissue)) +
  geom_boxplot(aes(fill = tissue), width = 0.1) +
  scale_y_continuous("Heteroplasmy %", limits = c(0, 105)) +
  scale_x_discrete("", labels = tissue_labs) +
  scale_fill_manual(values = tissue_cols)
fig_4b

# (C) Single cells E8.5 m.5019A>G
df_4c <- sc_df %>% filter(day  == "E8.5" & mutation == "m.5019")

fig_4c <- ggplot(df_4c, aes(x = mouse, y = het)) +
  theme_lvb + theme(legend.position = "none") +
  geom_violin(fill = mutation_cols["m.5019"]) +
  geom_boxplot(width = 0.1, fill = mutation_cols["m5019"]) +
  scale_y_continuous("Heteroplasmy %", limits = c(0, 105)) +
  scale_x_discrete("", labels = mouse_labs) 
fig_4c

# Compound figure
fig_4ac <- plot_grid(fig_4a, fig_4c,
                     labels = c("A", "C"), rel_heights = c(2, 1.1), 
                     ncol = 1, label_size = 10)
fig_4ac
fig4 <- plot_grid(fig_4ac, fig_4b, 
                  labels = c("", "B"), rel_widths = c(1, 1.5),
                  label_size = 10)
fig4

plot_save(fig4, "figures/04-fig4.jpg", size = 1, ar = 1/1.1)

#----- Main: Figure 5
# (A) Comparison between bulk-tissue and single-cell means
df_5a <- sc_df %>% 
  filter(day >= "P0" & mouse != "mouse1.1") %>% 
  group_by(mutation, day, mouse, tissue, mouse_id, het_bulk) %>% 
  summarise(het = mean(het))
df_5a_annot <- df_5a %>% filter(day == "P365" & mutation == "m.5024" & tissue == "spleen")

fig_5a <- ggplot(df_5a, aes(x = het, y = het_bulk, shape = tissue, colour = mutation)) +
  theme_lvb + 
  geom_point(size = 2, alpha = 0.75) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey30") +
  annotate("text", x = 35, y = 70, label = "P365", colour = "grey30") +
  annotate("curve", x = 35, y = 66, xend = df_5a_annot$het, yend = df_5a_annot$het_bulk, 
           curvature = .3, arrow = arrow(length = unit(2, "mm")), colour = "grey30") +
  scale_x_continuous("Mean single-cell heteroplasmy %", limits = c(20, 100)) +
  scale_y_continuous("Mean bulk tissue heteroplasmy %", limits = c(20, 100)) +
  scale_shape_manual(values = tissue_shapes[c("brain", "spleen")], labels = tissue_labs[c("brain", "spleen")]) +
  scale_colour_manual(values = mutation_cols, labels = mutation_labs) +
  guides(colour = "none")
fig_5a

# (B) Kimura distribution fits
# Idea: pick two representative distributions to plot: P100, Mouse 1, ACSA1+ve
# cells for m.5024C>T and m.5019A>G. Then calculate theoretical and sample
# quantiles and plot.
df_5b <- sc_df %>% 
  filter(day == "P100" & mouse == "mouse1" & cell_type == "acsa1_pos") %>% 
  left_join(sc_kimura_fits) %>% 
  select(mutation, het, kimura_p, kimura_b) %>% 
  as_tibble() %>% 
  group_split(mutation) %>% 
  lapply(function(df) {
    p <- df$kimura_p[1]
    b <- df$kimura_b[1]
    theo_cdf <- pkimura(seq(0, 1, .01), p, b)
    samp_het <- df$het
    samp_cdf <- ecdf(samp_het * .01)(seq(0, 1, .01))
    cbind(df[1, c(1, 3, 4)],
          prob = c(theo_cdf, samp_cdf),
          het = rep(seq(0, 1, .01), 2),
          type = rep(c("theoretical", "sample"), each = 101))
  }) %>% 
  bind_rows() %>% 
  mutate(type = case_when(type == "theoretical" ~ "theoretical",
                          mutation == "m.5024" ~ "m.5024",
                          T ~ "m.5019"))

type_labs <- c("theoretical" = "Theoretical", 
               "m.5024" = "m.5024C>T samples", 
               "m.5019" = "m.5019A>G samples")
type_cols <- c("grey40", mutation_cols)
names(type_cols)[1] <- "theoretical"

fig_5b <- ggplot(df_5b, aes(x = 100 * het, y = prob, colour = type)) +
  theme_lvb + theme(legend.position = "none") +
  facet_grid(~mutation, labeller = master_lab) +
  geom_line(size = 1.1, alpha = .7) +
  scale_x_continuous("Heteroplasmy %") +
  scale_y_continuous("CDF") +
  scale_colour_manual(values = type_cols, labels = type_labs)
fig_5b

# (C) Normalised variance
df_5c <- filter(sc_kimura_fits, mouse != "mouse1.1")
ts_5c <- get_tissue_labs(df_5c$tissue)

fig_5c <- ggplot(df_5c, aes(x = day, y = nvar)) +
  theme_lvb + theme(legend.position = "none") +
  facet_grid(.~mutation, labeller = labeller(mutation = mutation_labs)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(aes(colour = tissue)) +
  scale_x_discrete("Day") +
  scale_y_continuous("V'(h)", limits = c(0, 0.6), breaks = seq(0, 0.6, .1)) +
  scale_colour_manual(values = tissue_cols[ts_5c])
fig_5c

# (D) Variance estimates following Johnston & Jones (2016)
df_5d <- filter(jjsim_df, nc %in% c(200, 500, 1000))
sim_labs <- c("proliferating" =  "simulated 'spleen'",
              "quiescent" = "simulated 'brain'\nlower mtDNA turnover", 
              "quiescent_tau" = "simulated 'brain'\nsame mtDNA turnover")
sim_cols <- c("proliferating" =  "#e9851d", "quiescent" = "#788f33", "quiescent_tau" = "#788f33")
sim_ltys <- c("proliferating" =  "solid", "quiescent" = "solid", "quiescent_tau" = "dotted")

fig_5d <- ggplot(df_5d, aes(x = t, y = nvar, colour = cells, linetype = cells)) +
  theme_lvb +
  facet_grid(~nc, labeller = labeller(nc = c("200" = "Copy number 200",
                                             "500" = "Copy number 500",
                                             "1000" = "Copy number 1000"))) +
  geom_line(size = 1.2) +
  scale_x_continuous("Days") +
  scale_y_continuous("V'(h)", limits = c(0, 0.6)) +
  scale_color_manual(values = sim_cols, labels = sim_labs) +
  scale_linetype_manual(values = sim_ltys, labels = sim_labs)
fig_5d

# Compound figure
fig_5bc <- plot_grid(fig_5b, fig_5c, labels = c("B", "C"), ncol = 1, label_size = 10)
fig_5bc
fig_5abc <- plot_grid(fig_5a, fig_5bc, labels = c("A", ""), ncol = 2, label_size = 10)
fig_5abc
fig_5 <- plot_grid(fig_5abc, fig_5d, labels = c("", "D"), ncol = 1, label_size = 10, rel_heights = c(3, 2))
fig_5

plot_save(fig_5, "figures/04-fig5.jpg", size = 1, ar = 1)

#----- SI: Figure S3
# (A) Pyrosequencing repeats
fig_s3a <- ggplot(reps_df, aes(x = rep_1_1, y = rep_1_2)) +
  theme_lvb_nl +
  geom_point(size = 2, colour = mutation_cols["m.5024"]) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous("Pyrosequencing # 1", limits = c(40, 80), breaks = seq(40, 80, 10)) +
  scale_y_continuous("Pyrosequencing # 2", limits = c(40, 80), breaks = seq(40, 80, 10))
fig_s3a

# (B) PCR repeats
fig_s3b <- ggplot(reps_df, aes(x = rep_1_1, y = rep_2_1)) +
  theme_lvb_nl +
  geom_point(size = 2, colour = mutation_cols["m.5024"]) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous("PCR # 1", limits = c(40, 80), breaks = seq(40, 80, 10)) +
  scale_y_continuous("PCR # 2", limits = c(40, 80), breaks = seq(40, 80, 10))
fig_s3b

# (C) Correlation heatmap
df_s3c <- melt(cor(reps_df[,-(1:3)], use = "complete.obs"))
rep_labs <- paste0("PCR ", rep(1:3, each = 2), "\npyro ", rep(1:2, 3))
names(rep_labs) <- paste0("rep_", rep(1:3, each = 2), "_", rep(1:2, 3))
df_s3c$Var2 <- factor(df_s3c$Var2 , levels = rev(levels(df_s3c$Var1)))
df_s3c <- filter(df_s3c, Var1 != "mean_het" & Var2 != "mean_het")

fig_s3c <- ggplot(df_s3c, aes(x = Var1, y = Var2, fill = value)) +
  theme_lvb +
  geom_tile() +
  scale_fill_gradientn(colours = c("white", "#961f1f"), values = c(0, 1), breaks = seq(0, 1, .25), 
                       labels = seq(0, 1, .25), limits = c(0, 1)) +
  scale_x_discrete("", labels = rep_labs) +
  scale_y_discrete("", labels = rep_labs)
fig_s3c

# Compound figure
fig_s3ab <- plot_arrange(fig_s3a, fig_s3b, ncol = 1)
fig_s3ab
fig_s3 <- plot_grid(fig_s3ab, fig_s3c, 
                    labels = c("", "C"),  rel_widths = c(1, 2),
                    label_size = 10, ncol = 2)
fig_s3

plot_save(fig_s3, "figures/04-figS3.jpg", size = 1, ar = 1.75)

#----- SI: Figure S4
df_s4 <- cn_df %>% 
  mutate(lab = paste0(day, " ", mutation_labs[mutation], "\n", cell_labs[cell_type]))

fig_s4 <- ggplot(df_s4, aes(x = lab, y = n, fill = mutation)) +
  theme_lvb_nl +
  geom_violin(colour = "grey20") +
  geom_boxplot(width = .12, colour = "grey20") +
  scale_x_discrete("") +
  scale_y_continuous("Copy number") +
  scale_fill_manual(values = mutation_cols, labels = mutation_labs)
fig_s4

plot_save(fig_s4, "figures/04-figS4.jpg", ar = 3)

#----- SI: Tables S1 & S2
tab_s12 <- sc_kimura_fits[, c("mutation", "day", "mouse", "tissue", "cell_type",
                        "kimura_p", "kimura_b", "p.value")] %>% 
  mutate_if(is.numeric, function(x) round(x, 3)) %>% 
  mutate(mouse = mouse_labs[mouse],
         tissue = tissue_labs[tissue],
         cell_type = str_replace(cell_labs[cell_type], "\n", " "))

write_csv(tab_s12, "data/parsed/04-tabS12.csv")

#----- SI: Table S3
tab_s3 <- summary(model_best)$coefficients %>% 
  as.data.frame %>% 
  rownames_to_column("Name")

write_csv(tab_s3, "data/parsed/04-tabS3.csv")
