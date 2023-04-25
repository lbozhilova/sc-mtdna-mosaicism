######################
##### Plot setup #####
######################

# Last edited: 31/03/23 by LVB

# Description: Plot setup, including ggplot theme, labellers, colour palettes,
# and functions to combine and save figures with the appropriate dimensions.

#----- Packages
library("tidyverse")
library("MetBrewer")
library("cowplot")

#----- Theme
theme_lvb <- theme_minimal(base_size = 10) +
  theme(
    text = element_text(color = "gray20"),
    # Legend
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = 0.1,
    legend.title = element_blank(),
    # Axes
    axis.text = element_text(face = "italic"),
    axis.title.x = element_text(vjust = -1),        
    axis.title.y = element_text(vjust = 2),
    axis.ticks.x = element_line(color = "gray70", size = 0.2),
    axis.ticks.y = element_line(color = "gray70", size = 0.2),
    axis.line = element_line(color = "gray40", size = 0.3),
    axis.line.y = element_line(color = "gray40", size = 0.3),
    # Panel
    panel.grid.major = element_line(color = "gray70", size = 0.2),
    panel.grid.major.x = element_line(color = "gray70", size = 0.2),
    panel.border = element_rect(color = "gray70", fill = NA, size = 0.5))
theme_lvb_nl <- theme_lvb + 
  theme(legend.position = "none")

#----- Arranging and saving figures
plot_arrange <- function(...) plot_grid(..., label_size = 10, labels = "AUTO")

plot_save <- function(p, filename, size = 1, ar = 1, dev = "jpeg"){
  allowed_devs <- c("eps", "ps", "tex", "pdf", "jpeg", 
                    "tiff", "png", "bmp", "svg")
  if (!(dev %in% allowed_devs))
    stop("Invalid device.")
  if (dev != "jpeg" & !str_detect(filename, paste0("\\.", dev)))
    filename <- paste0(filename, ".", dev)
  if (dev == "jpeg" & !str_detect(filename, "\\.jpg|\\.jpeg"))
    filename <- paste0(filename, ".jpg")
  w <- round(180 * size)
  h <- w/ar
  ggsave(filename = filename,
         plot = p,
         width = w,
         height = h,
         units = "mm",
         device = dev)
}

#----- Colours and labels
# Mutations
mutation_names <- c("m.5024", "m.5019")
mutation_labs <- c("m.5024C>T", "m.5019A>G")
names(mutation_labs) <- mutation_names
mutation_cols <- met.brewer("Homer2", 4)[c(1, 3)]
names(mutation_cols) <- mutation_names

# Tissues
tissue_names <- c("brain", "spleen", "liver", "gut", "gonads", "muscle", "skin", "heart", "ear", "unknown")
tissue_labs <- c("Brain", "Spleen", "Liver", "Gut", "Gonads", "Muscle", "Skin", "Heart", "Ear biopsy", "Unknown")
names(tissue_labs) <- tissue_names
tissue_shapes <- c(19, 15, 17, 18, 5, 6, 7, 8, 10, 9)
names(tissue_shapes) <- tissue_names
tissue_cols <- c(met.brewer("Homer2", 2)[c(2, 1)], rep("grey40", 8))
names(tissue_cols) <- tissue_names

# Cell types
cell_names <- c("cd19_pos", "cd19_neg", "acsa1_pos", "acsa_prom", "psd95_pos", "unknown", "NA")
cell_labs <- c("CD19 +ve", "CD19 -ve", "ACSA-1 +ve", "ACSA-1 &\nPROM +ve", "PSD95 +ve", "Unknown", "Unknown")
names(cell_labs) <- cell_names
cell_shapes <- c(0, 15, 19, 1, 1, 9, 9)
names(cell_shapes) <- cell_names

# Mice
mouse_names <- c("mouse1", "mouse2", "mouse3", "mouse1.1", 
                 "embryo1", "embryo2", "embryo3")
mouse_labs <- c("Mouse 1", "Mouse 2", "Mouse 3", "Mouse 4",
                "Embryo 1", "Embryo 2", "Embryo 3")
names(mouse_labs) <- mouse_names

# Labeller
master_lab <- labeller(mutation = mutation_labs,
                       tissue = tissue_labs,
                       cell_type = cell_labs,
                       mouse = mouse_labs)