# Load required packages
if (!require("pacman", character.only = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(DiagrammeR, greta, dplyr, tidyr, binom, matrixStats, ggplot2, bayesplot,viridis, patchwork,ggridges,scoringRules,purrr)


