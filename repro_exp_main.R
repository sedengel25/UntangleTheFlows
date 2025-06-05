source("./src/exp_utils.R")
#################################################################################
# Config
################################################################################
# args <- commandArgs(trailingOnly = TRUE)
# int_row_index <- as.numeric(args[1])
int_row_index <- 1
df_pg_snn <- read.csv(here("input","exp_param.csv"))
param_row <- df_pg_snn[int_row_index, ]
int_k <- param_row$k
int_eps <- param_row$eps
int_minpts <- param_row$minpts


# PaCMAP-config
int_n_samples <- 10

# OD-flow config
# Options "euclid", "chebyshev", "network"
  # "euclid" corresponds to "manhattan_euclid"
  # "chebyshev" corresponds to "chebyshev_euclid"
  # "network" = corresponds to "manhattan_network"
char_dist_measure <- "chebyshev"
path_flows <- here("data", "synthetic_data_experiment")
char_dirs <- list.dirs(path_flows, recursive = FALSE)



lapply(char_dirs, function(x){
  char_dir <- basename(x)
  main(char_dir = char_dir,
       path_flows = path_flows,
       char_dist_measure = char_dist_measure,
       int_n_samples = int_n_samples,
       int_k = int_k,
       int_eps = int_eps,
       int_minpts = int_minpts)
})


