library(tidyverse)
library(here)
library(data.table)
#################################################################################
# Load data
################################################################################
# Options "euclid", "chebyshev", "network"
# "euclid" corresponds to "manhattan_euclid"
# "chebyshev" corresponds to "chebyshev_euclid"
# "network" = corresponds to "manhattan_network"
char.dist.measure <- "chebyshev"
path.import <- here("res_exp_processed", char.dist.measure)

list.files <- list.files(path.import, recursive = FALSE, pattern = "rds")
list.rds.files <- lapply(list.files, function(x){
  dt <- read_rds(here(path.import, x))
})
dt <- rbindlist(list.rds.files)
rm(list.rds.files)
gc()
#################################################################################
# Mutate
################################################################################
dt <- dt %>%
  mutate(
    algo_lab = toupper(algo), 
    cvi_name_lab = ifelse(cvi_name == "ari", "ARI",
                          ifelse(cvi_name == "folkes_mallows", "FMI",
                                 ifelse(cvi_name == "silhouette", "SWC", cvi_name))),
    type_lab = ifelse(type == "org", "Original",
                      ifelse(type == "pacmap", "PaCMAP", type))
    
  )
dt <- dt %>%
  group_by(data) %>%
  mutate(data_id = cur_group_id()) %>%
  ungroup() %>%
  as.data.table() %>%
  select(-data)
#################################################################################
# Check NAs
################################################################################
dt %>%
  group_by(type, algo, cvi_name) %>%
  summarise(
    na_count = sum(is.na(median_cvi_val))
  ) %>%
  as.data.frame()

dt <- dt %>%
  mutate(median_cvi_val = if_else(cvi_name == "silhouette" & is.na(median_cvi_val),
                                  -1, 
                                  median_cvi_val))

#################################################################################
# Check Duplicates
################################################################################
dt <- dt %>%
  distinct(type, data_id, cvi_name, algo, minpts, eps, k, .keep_all = TRUE) %>%
  ungroup()

if(!dir.exists(here("res_exp_final"))){
  dir.create(here("res_exp_final"), recursive = TRUE)
}

write_rds(dt, here("res_exp_final", paste0("dt_", char.dist.measure, ".rds")))