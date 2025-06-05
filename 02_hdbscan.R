library(tidyverse)
library(shiny)
library(lwgeom)
library(leaflet)
library(sf)
library(RColorBrewer)
library(dbscan)
library(data.table)
library(reticulate)
library(here)
library(aricode)
library(cluster)
source("./src/utils_hdbscan.R")
library(Rcpp)
library(RcppParallel)
Rcpp::sourceCpp("./src/utils_hdbscan.cpp")

#set.seed(123)
reticulate::use_virtualenv("r-reticulate")
np <- import("numpy")


################################################################################
# Functions
################################################################################
sys_run_pacmap <- function(path.python,
                           path.pacmap.src.file,
                           path.data,
                           int.id) {
  
  cmd <- paste(path.python, 
               path.pacmap.src.file, 
               "--path_file_directory", 
               path.data, 
               "--i", int.id)
  
  cat("Command:\n", cmd, "\n")
  
  output <- system2(path.python, args = c(path.pacmap.src.file,
                                          "--path_file_directory ",
                                          path.data,
                                          "--i ", int.id),
                    stdout = stdout, stderr = stderr, wait = TRUE)
  cat("system2: ", output, "\n")
  if(output!=0){
    stop("ERROR")
  }
}

get_manhattan_euclidean_distmat <- function(df.coord) {
  num.dist.o <- as.matrix(dist(cbind(df.coord$o_x, df.coord$o_y))) 
  num.dist.d <- as.matrix(dist(cbind(df.coord$d_x, df.coord$d_y))) 
  num.dist.final <- num.dist.o + num.dist.d
  return(num.dist.final)
}

get_distmat_chebyshev_distances <- function(df.coord) {
  num.dist.o <- as.matrix(dist(cbind(df.coord$o_x, df.coord$o_y))) 
  num.dist.d <- as.matrix(dist(cbind(df.coord$d_x, df.coord$d_y))) 
  num.dist.final <- pmax(num.dist.o, num.dist.d) 
  return(num.dist.final)
}


# Description: Creates a (ram-heavy) matrix containing the ND between OD flows 
# depending on the given road network and the given local node distance matrix
get_distmat_network_distances <- function(sf_trips, sf_network, dt_dist_mat) {
  dt_trips <- add_dist_start_end(sf_trips, sf_network)
  dt_trips_origin <- dt_trips %>%
    select(cluster_id, origin_id, flow_id, dist_to_start_origin, dist_to_end_origin) %>%
    rename(dist_to_start = dist_to_start_origin,
           dist_to_end = dist_to_end_origin,
           id_edge = origin_id)
  
  dt_trips_dest <- dt_trips %>%
    select(cluster_id, dest_id, flow_id, dist_to_start_dest, dist_to_end_dest) %>%
    rename(dist_to_start = dist_to_start_dest,
           dist_to_end = dist_to_end_dest,
           id_edge = dest_id)
  dt_network <- sf_network %>%
    as.data.table()
  matrix_dist_origin <- parallel_process_networks_matrix(dt_trips_origin, 
                                                         dt_network, 
                                                         dt_dist_mat, 
                                                         1)
  
  matrix_dist_dest <- parallel_process_networks_matrix(dt_trips_dest, 
                                                       dt_network, 
                                                       dt_dist_mat,
                                                       1)
  
  matrix_dist <- matrix_dist_origin + matrix_dist_dest
  return(matrix_dist)
}

calculate_distmat <- function(dist.measure, 
                              df.coord,
                              sf.trips = NULL,
                              sf.network = NULL,
                              dt.dist.mat = NULL){
  if(dist.measure == "manhattan_euclid"){
    dist.mat <- get_manhattan_euclidean_distmat(df.coord)
  } else if (dist.measure == "chebyshev_euclid"){
    dist.mat <- get_distmat_chebyshev_distances(df.coord)
  } else if (dist.measure == "manhattan_network"){
    
    if (is.null(sf.trips) || is.null(sf.network) || is.null(dt.dist.mat)) {
      stop("To calculate the 'manhattan_network'-distance, 'sf.trips', 
           'sf.network' and 'dt.dist.mat' need to be provided.")
    }
    
    dist.mat <- get_distmat_network_distances(sf.trips,
                                              sf.network,
                                              dt.dist.mat)
  }
  return(dist.mat)
}
################################################################################
# Calculate OD flow distance matrix
################################################################################
sf.network <- read_rds(here("input", "sf_network.rds"))
dt.dist.mat <- read_rds(here("input", "local_node_dist_mat.rds"))
char.dist.measure <- "chebyshev_euclid" # further options: "chebyshev_euclid", "manhattan_network"
char.dirs <- list.dirs(here("data", "synthetic_data_generated"), full.names = FALSE, recursive = FALSE)
char.dirs
char.data <- char.dirs[2]
path.data <- here("data", "synthetic_data_generated", char.data, char.dist.measure)

if(!dir.exists(here(path.data))){
  dir.create(path.data, recursive = TRUE)
} 

path.python <- "/usr/bin/python3.8"
path.pacmap.src.file <- here("src", "pacmap.py")
sf.trips <- read_rds(here("data", "synthetic_data_generated", char.data, "synthetic_flows.rds"))
sf.trips$flow_id <- 1:nrow(sf.trips)
int.id <- 1 # if you run pacmap several times on the same dataset
df.coord <- data.frame(
  o_x = st_coordinates(lwgeom::st_startpoint(sf.trips$geometry))[,1],
  o_y = st_coordinates(lwgeom::st_startpoint(sf.trips$geometry))[,2],
  d_x = st_coordinates(lwgeom::st_endpoint(sf.trips$geometry))[,1],
  d_y = st_coordinates(lwgeom::st_endpoint(sf.trips$geometry))[,2])

# Write 4d dataset since it is used to initialize PaCMAP in the 4d embedding space
write_csv(df.coord, here(path.data, "df_coordinates.csv"))
# Calculate OD flow distance matrix based on specific dist-measure
od.flow.dist.mat <- calculate_distmat(dist.measure = char.dist.measure,
                                      df.coord = df.coord,
                                      sf.trips = sf.trips,
                                      sf.network = sf.network,
                                      dt.dist.mat = dt.dist.mat)

################################################################################
# PaCMAP
################################################################################
# Save it as numpy-array
np$save(here(path.data, "dist_mat.npy"), od.flow.dist.mat)
od.flow.dist.mat <- as.dist(od.flow.dist.mat)
# Run PaCMAP using the OD flow distance matrix
sys_run_pacmap(path.python = path.python,
               path.pacmap.src.file = path.pacmap.src.file,
               path.data = path.data,
               int.id = int.id)


matrix.embedding <- np$load(here(path.data, paste0("embedding_4d_", int.id,".npy")))


################################################################################
# Clustering (mini-experiment)
################################################################################
results <- lapply(2:50, function(minpts){
  print(minpts)
  hdbscan.org <- hdbscan(od.flow.dist.mat, minPts = minpts)
  hdbscan.emb <- hdbscan(matrix.embedding, minPts = minpts)
  ari.org <- aricode::ARI(sf.trips$cluster_id, hdbscan.org$cluster)
  ari.emb <- aricode::ARI(sf.trips$cluster_id, hdbscan.emb$cluster)
  c(minpts, ari.org, ari.emb)
})

df.res <- as.data.frame(do.call(rbind, results))
colnames(df.res) <- c("minpts", "ari.org", "ari.emb")
df.long <- df.res %>%
  pivot_longer(cols = c(ari.org, ari.emb), 
               names_to = "Data", 
               values_to = "ARI") %>%
  mutate(Data = recode(Data, "ari.org" = "Original", "ari.emb" = "PaCMAP")) %>%
  rename(`Minimum Cluster Size` = minpts)




################################################################################
# Plotting
################################################################################
viridis_pal <- c(viridis::viridis(n = 3)[c(2,3)])

ggplot(df.long, 
       aes(x = `Minimum Cluster Size`, y = ARI,
           color = Data
       )) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.5) +
  scale_color_manual(values = viridis_pal) +
  scale_x_continuous(breaks = seq(2, 50, 1),  
                     labels = seq(2, 50, 1)) +
  xlab("Minimum Cluster Size\n") +
  ylab("ARI") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14L),
    axis.text = element_text(size = 14L),
    axis.title = element_text(size = 16L),
    strip.text = element_text(size = 14L),
    axis.title.x = element_text(margin = margin(t = 20))
  )
