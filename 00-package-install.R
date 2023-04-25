############################
##### Install packages #####
############################

# Last edited: 11/04/23 by LVB

# Description: Installs and updates necessary packages. Note these scripts were
# written and executed in R v.4.2.3. While code should have reasonable
# back-compatibility, you may need to update R. The way you do this may depend on
# your OS. 

# Note script <02-kimura-fits.R> will fail on Windows as the OS does not
# support forking. It will/should work on MacOS and Linux.

#----- Main
reqs_cran <- matrix(c(
  "cowplot", "1.1.1",
  "devtools", "2.4.3",
  "ggbeeswarm", "0.6.0",
  "MASS", "7.3.57",
  "MetBrewer", "0.2.0",
  "parallel", "4.2.3",
  "reshape2", "1.4.4",
  "tidyverse", "1.3.1"),
  ncol = 2, byrow = T,
  dimnames = list(c(), c("Package", "Version")))

to_install <- which(!(reqs_cran[, 1] %in% installed.packages()))
install.packages(reqs_cran[to_install, 1])

to_update <- which(apply(reqs_cran, 1, function(x) packageVersion(x[1]) < x[2]))
install.packages(reqs_cran[to_update, 1])

devtools::install_github("lbozhilova/kimura")

# Check packages are installed. This will NOT check package versions and should
# return TRUE.
all(c(reqs_cran[, 1], "kimura") %in% installed.packages())
