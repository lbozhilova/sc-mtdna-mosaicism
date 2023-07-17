#################################################
##### Summary statistics and other analysis #####
#################################################

# Last edited: 17/07/23 by LVB

# Description: Summary statistics and other analysis which appears ad-hoc in the
# paper and SI.

#----- Packages
library("tidyverse")

#----- Data
load("data/parsed/01-parsed-replicates.RData")
load("data/parsed/01-parsed-data.RData")
load("data/parsed/01-parsed-copy-numbers.RData")

#----- Reproducibility of pyrosequencing
# Mean absolute deviation of repeat pyrosequencing measurements
reps_df %>% 
  rowwise() %>% 
  mutate(mad = mean(abs(mean_het - c(rep_1_1, rep_1_2, 
                                     rep_2_1, rep_2_2, 
                                     rep_3_1, rep_3_2)), na.rm = T)) %>% 
  ungroup() %>% 
  summarise(mean_mad = mean(mad)) %>% 
  round(1)

# Corelations between PCR and pyrosequencing runs. Note that in rep_i_j, i is
# the index of the PCR run, and j is the index of the pyrosequencing run.
cor(reps_df[, 4:9], use = "pairwise.complete.obs") %>% round(3)

#----- Cell counts
# Total cell count
nrow(sc_df)

# m.5024C>T P100 cells per cell type
sc_df %>% 
  filter(mutation == "m.5024" & day == "P100") %>% 
  group_by(cell_type) %>% 
  summarise(n = n())

# m.5024C>T cells per day
sc_df %>% 
  filter(mutation == "m.5024") %>% 
  group_by(day) %>% 
  summarise(n = n())

# m.5024C>T number & proportion of homoplasmic cells 
sc_df %>% 
  filter(mutation == "m.5024") %>% 
  group_by(day) %>% 
  summarise(n = n(),
            n_homoplasmic = sum(het == 0 | het == 100)) %>% 
  mutate(prop_homoplasmic = round(100 * n_homoplasmic / n, 2))

# m.5019A>G cells per day
sc_df %>% 
  filter(mutation == "m.5019") %>% 
  group_by(day) %>% 
  summarise(n = n())

#----- Comparison between bulk-tissue and single-cell heteroplasmy
# m.5024C>T P100 single-cell mean vs bulk-tissue
sc_df %>% 
  filter(mutation == "m.5024" & day == "P100") %>% 
  group_by(mouse, tissue) %>% 
  summarise(mean_het = mean(het)) %>% 
  left_join((bt_df %>% filter(mutation == "m.5024" & day == "P100"))) %>% 
  mutate(diff = abs(mean_het - het_bulk)) %>%
  group_by(mutation) %>% 
  summarise(max_diff = round(max(diff), 2),
            mean_diff = round(mean(diff), 2),
            sd_diff = round(sd(diff), 2))

# >E8.5 correlation between single-cell and bulk-tissue heteroplasmy
sc_df %>% 
  filter(day > "E8.5") %>% 
  group_by(day, mouse, tissue, het_bulk) %>% 
  summarise(mean_het = mean(het)) %>%
  (function(x) cor(x$mean_het, x$het_bulk, use = "pairwise.complete.obs"))
  
#----- Single-cell heteroplasmy distributions
# m.5024C>T P100 single-cell IQRs
sc_df %>% 
  filter(mutation == "m.5024" & day == "P100") %>% 
  group_by(tissue) %>% 
  summarise(lower_q = quantile(het, .25),
            upper_q = quantile(het, .75))

# m.5019A>G P100 single-cell IQRs
sc_df %>% 
  filter(mutation == "m.5019" & day == "P100") %>% 
  group_by(tissue) %>% 
  summarise(lower_q = quantile(het, .25),
            upper_q = quantile(het, .75))

# m.5024C>T P100 test single-cell distributions within the same tissue
get_cell_type_pval <- function(het, cell_type) {
  x <- split(het, cell_type, drop = T)
  if (length(x) == 1) return(NA)
  ks.test(x[[1]], x[[2]])$p.value
}
sc_df %>% 
  filter(mutation == "m.5024" & day == "P100") %>% 
  group_by(mouse, tissue) %>% 
  summarise(pval = get_cell_type_pval(het, cell_type)) %>% 
  group_by(tissue) %>% 
  summarise(min_pval = min(pval))

# Table S2. All >E8.5 test single-cell distributions within the same tissue 
table_s2 <- sc_df %>% 
  group_by(mutation, day, mouse, tissue) %>% 
  summarise(n = n(),
            pval = get_cell_type_pval(het, cell_type)) %>% 
  group_by(tissue) %>% 
  mutate(is_significant = (pval < 0.05 / n())) %>% 
  mutate(pval = round(pval, 3))
write_csv(table_s2, "data/parsed/05-tabS2.csv")

#----- Normalised heteroplasmy variance
nvar <- function(H) {H <- .01 * H; var(H) / (mean(H) * (1 - mean(H)))}

# P100 bulk tissue normalised variance
bt_df %>% 
  filter(day == "P100") %>% 
  group_by(mutation, mouse) %>% 
  summarise(nvar = nvar(het_bulk)) %>% 
  group_by(mutation) %>% 
  summarise(mean = round(mean(nvar), 3),
            sd = round(sd(nvar), 3))

# >E8.5 bulk tissue normalised variance
bt_df %>% 
  filter(day > "E8.5") %>% 
  group_by(mutation, mouse) %>% 
  summarise(nvar = nvar(het_bulk)) %>% 
  group_by(mutation) %>% 
  summarise(mean = round(mean(nvar), 3),
            sd = round(sd(nvar), 3))

# P100 and E8.5 single cell normalised variance
sc_df %>% 
  filter(day %in% c("P100", "E8.5")) %>% 
  group_by(mutation, day, mouse, tissue, cell_type) %>% 
  summarise(nvar = nvar(het)) %>% 
  group_by(mutation, day) %>% 
  summarise(mean = round(mean(nvar), 3),
            sd = round(sd(nvar), 3))

# m.5024C>T P365 normalised variance compared to earlier time points
nvs <- sc_df %>% 
  filter(mutation == "m.5024") %>% 
  group_by(day, mouse, tissue, cell_type) %>% 
  summarise(nvar = nvar(het)) %>% 
  (function(x) split(x$nvar, x$day))
sapply(nvs[1:4], function(x) wilcox.test(x, nvs[[5]], alternative = "less")$p.value)

# m.5019A>G P365 normalised variance compared to earlier time points
nvs <- sc_df %>% 
  filter(mutation == "m.5019") %>% 
  group_by(day, mouse, tissue, cell_type) %>% 
  summarise(nvar = nvar(het)) %>% 
  (function(x) split(x$nvar, x$day))
sapply(nvs[1:4], function(x) wilcox.test(x, nvs[[5]], alternative = "less")$p.value)

# P6 to P365 increase in normalised variance
sc_df %>% 
  filter(day %in% c("P6", "P365")) %>% 
  group_by(mutation, day, mouse, tissue) %>% 
  summarise(nvar = nvar(het)) %>% 
  group_by(mutation, day) %>% 
  summarise(mean_nvar = round(mean(nvar), 3),
            sd_nvar = round(sd(nvar), 3))

# comparison of spleen and brain normalised variances
sc_df %>% 
  filter(day > "E8.5") %>% 
  group_by(mutation, day, mouse, tissue) %>% 
  summarise(nvar = nvar(het)) %>% 
  pivot_wider(names_from = tissue, values_from = nvar) %>% 
  (function(x) wilcox.test(x$brain, x$spleen, paired = T))

#----- Copy number
# E8.5 vs other days: copy number
cn_df %>% 
  group_by(day == "E8.5") %>% 
  summarise(mean_cn = mean(n),
            sd_cn = sd(n))

# m.5024C>T and m.5019A>G copy number comparison
cns <- filter(cn_df, day == "P0")
cns <- split(cns$n, cns$mutation)
wilcox.test(cns[[1]], cns[[2]])
