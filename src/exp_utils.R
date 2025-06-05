################################################################################
# Libraries
################################################################################
library(dplyr)
library(readr)
library(aricode)
library(here)
library(data.table)
library(dbscan)
library(reticulate)
library(clusterCrit)
library(cluster)
reticulate::use_virtualenv("r-reticulate")
np <- import("numpy")




cvi_snn <- function(c_cluster_true = NULL,
                    c_cluster_pred,
                    matrix_eval,
                    int_k=NULL,
                    int_eps=NULL,
                    int_minpts=NULL) {
  
  num_ari <- aricode::ARI(c1 = c_cluster_true, c2 = c_cluster_pred)
  list_cluster_crit_ext <- clusterCrit::extCriteria(
    part1 = c_cluster_true,
    part2 = c_cluster_pred,
    crit = c("Folkes_Mallows"))
  
  int_n_cl <- length(unique(c_cluster_pred)) - 1
  int_n_noise <- length(which(c_cluster_pred==0))

  sil_res <- NULL
  if (length(unique(c_cluster_pred)) < 2) {
    cat("Only one cluster detected. Returning -1 for silhouette score.\n")
    sil_res <- NA
  } else {
    sil_res <- silhouette(c_cluster_pred, dmatrix = as.matrix(matrix_eval))
    sil_res <- mean(sil_res[, 3])
  }
  

  c_cvi_char <- c(
    "ari",
    names(list_cluster_crit_ext),
    "silhouette"
  )
  c_cvi_num <- c(
    num_ari,
    as.numeric(list_cluster_crit_ext),
    sil_res
  )


  return(list("names" = c_cvi_char,
              "values" = c_cvi_num,
              "n_cl_pred" = int_n_cl,
              "n_noise" = int_n_noise))
}





run_all_get_cvi <- function(matrix_cluster,
                            matrix_eval,
                            int_minpts,
                            int_eps,
                            dist_eps,
                            int_k,
                            c_cluster_true) {
  
  
  algorithms <- list(
    "dbscan" = function(data, eps, minpts) dbscan::dbscan(data, eps = eps, minPts = minpts)$cluster,
    "hdbscan" = function(data, minpts) hdbscan(data, minPts = minpts)$cluster
  )
  
  if (int_k!=0) {
    algorithms[["snn"]] <- function(data, eps, minpts, k) sNNclust(data, k, eps, minpts)$cluster
  }
  results <- lapply(names(algorithms), function(algo) {
    print(algo)
    c_cluster_pred <- if (algo == "hdbscan") {
      algorithms[[algo]](matrix_cluster, int_minpts)
    } else if (algo == "snn") {
      algorithms[[algo]](matrix_cluster, int_eps, int_minpts, int_k)
    } else {
      algorithms[[algo]](matrix_cluster, dist_eps, int_minpts)
    }
    
    c_cluster_pred <- as.integer(c_cluster_pred)

    cvi_res <- cvi_snn(c_cluster_true = c_cluster_true, c_cluster_pred = c_cluster_pred,
                       matrix_eval,
                       int_k = int_k, int_eps = int_eps, int_minpts = int_minpts)
    
    k_value <- if (algo %in% c("hdbscan", "dbscan")) NA else int_k
    eps_value <- if (algo == "hdbscan") NA else int_eps
    print(cvi_res)
    num_metrics <- length(unlist(cvi_res$names))
    
    df_cvi <- data.frame(
      algo = rep(algo, num_metrics),
      k = rep(k_value, num_metrics),
      eps = rep(eps_value, num_metrics),
      minpts = rep(int_minpts, num_metrics),
      n_cl = rep(cvi_res$n_cl, num_metrics),
      n_noise = rep(cvi_res$n_noise, num_metrics),
      cvi_name = unlist(cvi_res$names),
      cvi_value = unlist(cvi_res$values)
    )
    df_cvi
  })
  df_results <- do.call(rbind, results)
  return(df_results)
}






run_pacmap_instances <- function(path_pacmap_input,
                                 char_dist_measure,
                                 int_n_samples,
                                 int_k,
                                 int_eps,
                                 int_minpts,
                                 c_cluster_true) {
  
  list_cvi <- lapply(1:int_n_samples, function(i){
    int_i <- i
    matrix_euclid_dist_org <- np$load(here(path_pacmap_input, char_dist_measure, "dist_mat.npy")) %>% as.matrix()
    print("Original OD flow distance matrix:")
    print(matrix_euclid_dist_org[1:5,1:5])
    
    matrix_euclid_dist_org <- as.dist(matrix_euclid_dist_org)


    char_embedding <- paste0("embedding_4d_", int_i)
    
    c_eps <- read_rds(here(path_pacmap_input, char_dist_measure, "eps_minpts_emb.rds"))
    dist_eps <- c_eps[int_minpts,int_eps] %>% as.numeric()
  
    matrix_euclid_embedding <- np$load(here(path_pacmap_input, char_dist_measure, 
                                            paste0(char_embedding, ".npy")))
    print("PaCMAP transformed embedding space:")
    print(head(matrix_euclid_embedding))
    df_sample_res <- run_all_get_cvi(
      matrix_cluster = matrix_euclid_embedding,
      matrix_eval = matrix_euclid_dist_org,
      int_minpts = int_minpts,
      dist_eps = dist_eps,
      int_eps = int_eps,
      int_k = int_k,
      c_cluster_true = c_cluster_true)
    df_sample_res
  })
  
  dt_res <- rbindlist(list_cvi)
  return(dt_res)
}






run_org_instance <- function(path_pacmap_input,
                             c_cluster_true,
                             char_dist_measure,
                             int_k,
                             int_eps,
                             int_minpts) {
  c_eps <- read_rds(here(path_pacmap_input, char_dist_measure, "eps_minpts_org.rds"))
  dist_eps <- c_eps[int_minpts,int_eps] %>% as.numeric()
  
  matrix_euclid_dist_org <- np$load(here(path_pacmap_input, char_dist_measure, "dist_mat.npy")) %>%
    as.matrix() %>% as.dist()
  df_res <- run_all_get_cvi(
    matrix_cluster = matrix_euclid_dist_org,
    matrix_eval = matrix_euclid_dist_org,
    int_minpts = int_minpts,
    dist_eps = dist_eps,
    int_eps = int_eps,
    int_k = int_k,
    c_cluster_true = c_cluster_true)
  return(df_res)
}


  
main <- function(char_dir,
                 path_flows,
                 char_dist_measure,
                 int_n_samples,
                 int_k,
                 int_eps,
                 int_minpts) {
  
  path_pacmap_input <- here(path_flows, char_dir)
  sf_trips <- here(path_pacmap_input, "synthetic_flows.rds") %>%
    read_rds()
  c_cluster_true <- sf_trips$cluster_id
  dt_res_pacmap <- run_pacmap_instances(path_pacmap_input = path_pacmap_input,
                                        char_dist_measure = char_dist_measure,
                                        int_n_samples = int_n_samples,
                                        int_k = int_k,
                                        int_eps = int_eps,
                                        int_minpts = int_minpts,
                                        c_cluster_true = c_cluster_true)
  
  dt_res_pacmap$type = rep("pacmap", nrow(dt_res_pacmap))
  path_file <- here("res_exp_intermediate", char_dist_measure, char_dir)
  
  if(!dir.exists(here(path_file))){
    dir.create(path_file, recursive = TRUE)
  }
  
  
  dt_res_org <- run_org_instance(path_pacmap_input = path_pacmap_input,
                                 c_cluster_true = c_cluster_true,
                                 char_dist_measure = char_dist_measure,
                                 int_k = int_k,
                                 int_eps = int_eps,
                                 int_minpts = int_minpts)
  
  
  dt_res_org$type = rep("org", nrow(dt_res_org))
  dt_res <- rbind(dt_res_pacmap, dt_res_org)
  char_file <- paste0("cvi_h_dbscan_snn_",
                      int_k,
                      "_",
                      int_eps,
                      "_",
                      int_minpts,
                      ".rds")
  write_rds(dt_res, here(path_file, char_file))
}


