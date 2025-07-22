options(tidyverse.quiet = TRUE)
suppressPackageStartupMessages({
  library(misha)
  library(misha.ext)
  library(prego)
  library(iceqream)
  library(mcATAC)
  library(tgutil)
  library(tgstat)
  library(tidyverse)
  library(purrr)
  library(glue)
  library(patchwork)
  library(here)
  library(Matrix)
})

options(gmax.data.size = 1e+9)

source(here("code/themes.R"))

ggplot2::theme_set(theme_arial_classic(8))

scripts <- list.files(here::here("code"), pattern = "\\.(r|R)$", full.names = TRUE, recursive = TRUE)
exclude <- c(here("code/init.R"), here("code/themes.R"), here("code/download_data.R"), here("code/install.R"))
scripts <- scripts[!(scripts %in% exclude)]

purrr::walk(scripts, source)

gsetroot(here("data/mm10"))