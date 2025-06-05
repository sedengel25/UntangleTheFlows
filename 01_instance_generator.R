source("./src/utils_instance_generator.R")
################################################################################
# Config
################################################################################
sf_network <- read_rds(here("input", "sf_network.rds"))
sf_convex_hull <- st_convex_hull(st_union(sf_network))
dt_dist_mat <- read_rds(here("input", "local_node_dist_mat.rds"))

num_avg_radius <- main_get_avg_radius_chull(sf_convex_hull)
num_min_length <- num_avg_radius/6
num_max_length <- num_avg_radius*0.8
int_min_n_cl <- 50
int_max_n_cl <- 100
int_n_clusters <- sample(int_min_n_cl:int_max_n_cl, 1)
int_min_cl_size <- 15
int_max_cl_size <- 20
################################################################################
# Main
################################################################################
sf_clusters <- main_create_synth_dataset(int_n_clusters = int_n_clusters, 
                                         num_avg_radius = num_avg_radius,
                                         int_min_cl_size = int_min_cl_size,
                                         int_max_cl_size = int_max_cl_size,
                                         num_var_angle_direction = num_var_angle_direction,
                                         num_min_length = num_min_length,
                                         num_max_length = num_max_length)

char_id <- uuid::UUIDgenerate(use.time = TRUE)
path <- here("data", "synthetic_data_generated", char_id) 
if(!dir.exists(here(path))){
  dir.create(path, recursive = TRUE)
} else {
  warning(paste0("Directory", path, " already exists."))
}


write_rds(sf_clusters, here(path, "synthetic_flows.rds"))



colors <- c(rgb(211/255, 211/255, 211/255, 0.5), rainbow(int_n_clusters))
names(colors) <- c("0", as.character(1:int_n_clusters))
ggplot(data = sf_clusters) +
  geom_sf(data=sf_network) +
  geom_sf(aes(color = as.character(cluster_id)), size = 1) +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 18L),
    axis.text = element_text(size = 16L),
    axis.title = element_text(size = 18L),
    strip.text = element_text(size = 18L)
  )



