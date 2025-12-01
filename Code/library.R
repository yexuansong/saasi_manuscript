required_packages <- c(
  "tidyr", "stringr", "treeio", "tidyverse", "ggalluvial",
  "ggpubr", "irr", "stats4", "ape", "deSolve", "strap", "ips",
  "ggplot2", "ggtree", "reshape2", "ggstance", "diversitree",
  "ggimage", "patchwork", "janitor", "phytools", "dplyr",
  "rstatix", "jsonlite", "rootSolve"
)

missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing_packages) > 0) {
  warning(paste0(
    "Missing packages: ", paste(missing_packages, collapse = ", ")
  ))
} else {
  invisible(lapply(required_packages, library, character.only = TRUE))
  cat("Load all required packages\n")
}