library(here)
library(readr)
library(dplyr)
library(data.table)
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

char_config <- paste0(int_k, "_",
                      int_eps, "_",
                      int_minpts)

# Options "euclid", "chebyshev", "network"
  # "euclid" corresponds to "manhattan_euclid"
  # "chebyshev" corresponds to "chebyshev_euclid"
  # "network" = corresponds to "manhattan_network"
char_dist_measure <- "chebyshev"
char_file <- paste0("cvi_h_dbscan_snn_", char_config, ".rds")
path_res <- here("res_exp_intermediate", char_dist_measure)
char_dirs <- list.dirs(path_res, recursive = FALSE)
char_results <- paste0("cvi_h_dbscan_snn_grouped_", char_config, ".rds")
path_res_final <- here("res_exp_processed", char_dist_measure)

if(!dir.exists(here(path_res_final))){
  dir.create(path_res_final, recursive = TRUE)
}

list_grouped_rds <- lapply(char_dirs, function(x){
  char_data <- basename(x)
  if(!file.exists(here(x, char_file))){
    print("Skip!")
    return(NULL)
  }
  dt <- read_rds(here(x, char_file))
  setDT(dt)
  dt_median <- dt[, .(
    median_cvi_val = median(cvi_value, na.rm = TRUE),
    k = unique(k),
    eps = unique(eps),
    minpts = unique(minpts)
  ), by = list(type, cvi_name, algo)]
  print(dt_median)
  dt_median$data <- rep(char_data, nrow(dt_median))
  dt_median
})



dt_grouped_rds <- rbindlist(list_grouped_rds)
write_rds(dt_grouped_rds, here(path_res_final, char_results))
