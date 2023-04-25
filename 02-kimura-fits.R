#######################
##### Kimura fits #####
#######################

# Last edited: 31/03/23 by LVB

# Description: Fitting the Kimura distribution to single-cell data.

#---- Packages
library("tidyverse")
library("kimura")
library("parallel")

#---- Data
load("data/parsed/01-parsed-data.RData")

#---- Kimura fits
summary_df <- function(df) {
  ht <- .01 * df$het
  cbind(df[1, setdiff(colnames(df), c("het"))],
        data.frame(no_cells = nrow(df),
                   mean = mean(ht),
                   median = median(ht),
                   nvar = var(ht) / (mean(ht) * (1 - mean(ht)))))
}

test_kim <- function(df) test_kimura(.01 * df$het)
reform_kim <- function(df) {
  tst <- test_kim(df)
  cbind(summary_df(df),
        kimura_p = tst$statistic["p"],
        kimura_b = tst$statistic["b"],
        p.value = tst$p.value)
}

sc_kimura_fits <- sc_df %>% 
  split(sc_df$cell_id) %>% 
  mclapply(reform_kim, mc.cores = 9) %>% 
  bind_rows()

save(sc_kimura_fits, file = "data/parsed/02-kimura-fits.RData")
