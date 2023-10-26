# High-throughput single cell analysis reveals progressive mitochondrial DNA mosaicism throughout life

This repository contains all scripts necessary to reproduce the data analysis and figures presented in [Glynos _et al._ (2023)](https://www.science.org/doi/10.1126/sciadv.adi4038).

## Software requirements

### Operating system

This repository was developed and tested on MacOS, and should run in the same way on Linux distributions. If you are running code on Windows, which does not support forking via `mclapply()`, `L38` of `02-kimura-fits.R` will throw and error and should be replaced with a suitable `parLapply()` call.

### Packages

Scripts were written for and executed in `R v.4.2.3`. They employ a number of `R` packages available on CRAN:
- `cowplot v.1.1.1`,
- `ggbeeswarm v.0.6.0`,
- `MASS v.7.3.57`,
- `MetBrewer v.0.2.0`,
- `parallel v.4.2.3`,
- `reshape2 v.1.4.4`, and
- `tidyverse v.1.3.1`.

Fitting to the two-parameter Kimura distribution was done through [`kimura v.0.0.0.9000`](https://github.com/lbozhilova/kimura), which can be installed by running `devtools::install_github("lbozhilova/kimura")`. To install and update all necessary packages, run `source(00-package-install.R)`.

All code should be both forward- and back-compatible within reason, but please remember to run regular updates.

## Setup and script outline

Raw data can be found in `data/raw/`. Parsed data, including data tables and summaries, is in `data/parsed/`, and figures are in `figures/`.

Scripts `01-05` should principally be run in numerical order from the parent directory. Scripts starting with `00` are auxiliary. Whenever a file (e.g. a data file or a figure) is generated by a script, that file is given the prefix of the corresponding script. For example, the code for generating `04-fig1.jpg` can be found in `04-figures-tables.R `. Below is a brief breakdown of all scripts.

`00-package-install.R` installs and updates necessary packages.

`00-plot-setup.R` contains the `ggplot2` theme and colour palettes used for generating all figures.

`01-data-preproc.R` is used to parse the raw data, and do all necessary preprocessing.

`02-kimura-fits.R` is used to fit single-cell heteroplasmy samples to the two-parameter Kimura distribution, and to calculate normalised heteroplasmy shift.

`03-variance-models.R` is used to investigate normalised heteroplasmy variances, first using linear models, and then using the [Johnston & Jones (2016)](https://pubmed.ncbi.nlm.nih.gov/27843124/) model.

`04-figures-tables.R` generates main and supplementary figures and tables, as seen in the manuscript.

`05-summary-statistics.R` is used to calculate all summary statistics and hypothesis tests which appear _ad hoc_ in the manuscript. 

## References

If you use any of the code, methods or data presented here, please cite us:

A. Glynos, L. V. Bozhilova, M. Frison, S. Burr, J. B. Stewart, and P. F. Chinnery, "High-throughput single cell analysis reveals progressive mitochondrial DNA mosaicism throughout life," _Science Advances_ 9,eadi4038 (2023)
