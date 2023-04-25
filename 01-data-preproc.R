##############################
##### Data preprocessing #####
##############################

# Last edited: 06/04/23 by LVB

# Description: Data pre-processing. The *only* purpose of this script is to
# clean the data stored in data/raw and save data files for analysis in
# data/parsed.

#----- Packages
library("tidyverse")

#----- Load data
sc_df <- read_csv("data/raw/single-cell-data-20220824.csv")
bt_df <- read_csv("data/raw/bulk-tissue-data-20220824.csv")
reps_df <- read_csv("data/raw/embryo-repeats-20230217.csv")
cn_df <- read_csv("data/raw/cn-20230328.csv")

# Pivot and tidy copy number data
cn_df <- cn_df %>% 
  pivot_longer(cols = 1:5, names_to = "name", values_to = "n") %>% 
  mutate(mutation = c("m.5024", "m.5019")[1 + (name != "ACSA/Prom 5024")]) %>% 
  mutate(cell_type = case_when(name %in% c("ACSA/Prom 5024", "ACSA/Prom 5019") ~ "acsa_prom",
                               name == "E8.5" ~ "unknown",
                               name == "B-cells 5019" ~ "cd19_pos",
                               T ~ "psd95_pos")) %>% 
  mutate(day = case_when(cell_type == "acsa_prom" ~ "P0",
                         cell_type == "unknown" ~ "E8.5",
                         T ~ "P100"))

#----- Fix columns
# Factor mutation
mutation_lvls <- c("m.5024", "m.5019")
sc_df$mutation <- factor(sc_df$mutation, levels = mutation_lvls)
bt_df$mutation <- factor(bt_df$mutation, levels = mutation_lvls)
reps_df$mutation <- factor(c("m.5024"), levels = mutation_lvls)
cn_df$mutation <- factor(cn_df$mutation, levels = mutation_lvls)

# Factor day
day_lvls <- c("E8.5", "P0", "P6", "P100", "P365")
sc_df$day <- factor(sc_df$day, levels = day_lvls, ordered = T)
bt_df$day <- factor(bt_df$day, levels = day_lvls, ordered = T)
reps_df$day <- factor(reps_df$day, levels = day_lvls, ordered = T)
cn_df$day <- factor(cn_df$day, levels = day_lvls, ordered = T)

# Factor tissue
tissue_lvls <- c("brain", "spleen", "liver", "gut", "gonads", 
                 "muscle", "skin", "heart", "ear", "unknown")
sc_df$tissue <- factor(sc_df$tissue, levels = tissue_lvls)
bt_df$tissue <- factor(bt_df$tissue, levels = tissue_lvls)

# Factor cell type
cell_lvls <- c("cd19_pos", "cd19_neg", "acsa1_pos", 
               "acsa_prom", "psd95_pos", "unknown")
sc_df$cell_type[is.na(sc_df$cell_type)] <- "unknown"
sc_df$cell_type <- factor(sc_df$cell_type, levels = cell_lvls)
cn_df$cell_type <- factor(cn_df$cell_type, levels = cell_lvls)

# Unique mouse and cell type identifiers
sc_df <- sc_df %>% 
  rowwise() %>% 
  mutate(mouse_id = paste(mutation, day, mouse, sep = "_")) %>% 
  mutate(cell_id = paste(mouse_id, cell_type, sep = "_"))
bt_df <- bt_df %>% 
  rowwise() %>% 
  mutate(mouse_id = paste(mutation, day, mouse, sep = "_")) 

#----- Add mean (bulk) heteroplasmy measurements
# Single cells
sc_df <- left_join(
  sc_df,
  (bt_df %>% group_by(mouse_id) %>% summarise(het_bulk = mean(het_bulk))))

# Repeats
reps_df <- reps_df %>%
  rowwise() %>% 
  mutate(mean_het = mean(c(rep_1_1, rep_1_2, rep_2_1, rep_2_2, rep_3_1, rep_3_2),
                         na.rm = T))

#----- Filter outliers + extra mouse
sc_df <- filter(sc_df, mouse != "mouse1.1")
bt_df <- filter(bt_df, mouse != "mouse1.1")

# Remove CN outliers as per Tukey's fences (k = 1.5)
cn_df <- inner_join(
  cn_df,
  (cn_df %>% 
     group_by(name) %>% 
     summarise(fence_lower = quantile(n, .25) - 1.5 * IQR(n),
               fence_higher = quantile(n, .75) + 1.5 * IQR(n)))) %>%
  filter(n >= fence_lower & n <= fence_higher) %>% 
  dplyr::select(mutation, day, cell_type, n)

#----- Save
save(sc_df, bt_df, file = "data/parsed/01-parsed-data.RData")
save(reps_df, file = "data/parsed/01-parsed-replicates.RData")
save(cn_df, file = "data/parsed/01-parsed-copy-numbers.RData")
