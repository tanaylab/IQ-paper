CRAN_PACKAGES <- c(
  "tidyverse",   
  "tgutil",
  "tgstat",
  "glue",
  "here",
  "purrr",
  "Matrix",
  "patchwork",
  "princurve"
)

# install packages if not already installed
if (length(setdiff(CRAN_PACKAGES, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(CRAN_PACKAGES, rownames(installed.packages())))
}

GITHUB_PACKAGES <- c(
  "misha",
  "misha.ext",
  "prego",
  "iceqream",
  "mcATAC"
)

# install packages if not already installed
if (length(setdiff(GITHUB_PACKAGES, rownames(installed.packages()))) > 0) {
  remotes::install_github(paste0("tanaylab/", setdiff(GITHUB_PACKAGES, rownames(installed.packages()))))
}


