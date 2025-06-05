################################################################################
# Libraries
################################################################################
library(sf)
library(dplyr)
library(lwgeom)
library(tidygraph)
library(here)
library(readr)
library(sfnetworks)
library(ggplot2)
library(units)
library(clue)
library(plotly)
library(igraph)
################################################################################
# Functions
################################################################################
# Documentation: helper_get_degree_furthest_direction
# Usage: helper_get_degree_furthest_direction(sf_convex_hull, sf_point)
# Description: Gets the direction of the furthest point on the convex hull and
# the distance to that point (from the given point 'sf_point')
# Args/Options: sf_convex_hull, sf_point
# Returns: list
# Output: ...
# Action: ...
helper_get_degree_furthest_direction <- function(sf_convex_hull, sf_point) {
  sfc_hull_vertices <- st_cast(sf_convex_hull, "POINT")
  num_distances <- st_distance(sf_point, sfc_hull_vertices)
  int_max_distance_index <- which.max(num_distances)
  num_max_dist <- max(num_distances) %>% as.numeric()
  sfc_farthest_point <- sfc_hull_vertices[int_max_distance_index, ]
  
  # Koordinaten extrahieren
  num_point_coords <- st_coordinates(sf_point)
  num_farthest_coords <- st_coordinates(sfc_farthest_point)
  
  # Winkel berechnen
  num_angle <- atan2(
    num_farthest_coords[2] - num_point_coords[2],
    num_farthest_coords[1] - num_point_coords[1]
  )
  
  # num_angle_degrees <- (num_angle * 180 / pi) %% 360
  # x_start <- num_point_coords[1]
  # y_start <- num_point_coords[2]
  # x_end <- x_start + num_max_dist * cos(num_angle)
  # y_end <- y_start + num_max_dist * sin(num_angle)
  # p <- ggplot() +
  #   geom_sf(data = sf_convex_hull, fill = NA, color = "black") +
  #   geom_sf(data = sf_point, color = "red", size = 3) +
  #   geom_segment(
  #     aes(x = x_start, y = y_start, 
  #         xend = x_end, yend = y_end),
  #     arrow = arrow(length = unit(0.3, "cm")), color = "green", linewidth = 1
  #   ) +
  #   annotate("text", 
  #            x = (x_start + x_end) / 2, 
  #            y = (y_start + y_end) / 2, 
  #            label = paste0(round(num_angle_degrees, 1), "°"), 
  #            color = "black",
  #            size = 5, fontface = "bold")
  # print(p)
  
  return(list("degree" = num_angle,
              "distance" = num_max_dist))
}


# Documentation: helper_get_line_ids_of_mapped_points
# Usage: helper_get_line_ids_of_mapped_points(sf_points, sf_network)
# Description: Joins buffered points with 'sf_network' and extracts id, source and target
# Args/Options: sf_points, sf_network
# Returns: sf-dataframe
# Output: ...
# Action: ...
helper_get_line_ids_of_mapped_points <- function(sf_points, sf_network) {
  # Create small buffers around the points...
  sf_points_buffered <- st_buffer(sf_points, dist = 1e-7)  
  #...to get the linestrings and its attributes...
  sf_points_with_ids <- st_join(sf_points_buffered, sf_network, join = st_intersects)
  # ...and then get the centroids of the buffer's polygons
  sf_centroids <- st_centroid(sf_points_with_ids)
  sf_centroids <- sf_centroids %>%
    select(point_id, id, source, target, x) %>%
    rename("geom" = "x")
  sf_centroids <- sf_centroids %>%
    group_by(point_id) %>%  
    slice(1) %>%               
    ungroup()
  return(sf_centroids)
}

# Documentation: helper_snap_single_point_to_network
# Usage: helper_snap_single_point_to_network(sf_point, sf_network)
# Description: Snap a single point to the road network
# Args/Options: sf_point, sf_network
# Returns: list
# Output: ...
# Action: ...
helper_snap_single_point_to_network <- function(sf_point, sf_network) {
  int_nearest_feature <- st_nearest_feature(sf_point, sf_network)
  sf_nearest_line <- sf_network[int_nearest_feature, ]
  sf_nearest_point <- st_nearest_points(sf_point, sf_nearest_line) %>% 
    st_cast("POINT") %>% 
    .[2]
  return(sf_nearest_point)
}



# Documentation: create_ellipse_with_points
# Usage: create_ellipse_with_points(sf_network, num_center_point, num_width_m,
# num_height_m, int_n)
# Description: 
# 1. Gets a random road segment
# 2. Draws an ellipse around the center-point of that road segment
# 3. Gets all road segments within that ellipse
# 4. Converts sf road segments into a sfnetwork to keep the ones that
# form the largest connected network
# 5. Samples n points on the corresponding road segments 
# Args/Options: sf_network, num_center_point, num_width_m, num_height_m, int_n, 
# Returns: sf-dataframe
# Output: ...
# Action: ...
create_ellipse_with_points <- function(sf_network, 
                                       num_center_point, 
                                       num_width_m, 
                                       num_height_m, 
                                       int_n) {
  
  
  # Draw the ellipse 
  num_angles <- seq(0, 2 * pi, length.out = 100)
  matrix_unit_circle <- cbind(cos(num_angles), sin(num_angles))
  matrix_ellipse_points <- matrix_unit_circle %*% diag(c(num_width_m / 2, num_height_m / 2))
  matrix_ellipse_points <- sweep(matrix_ellipse_points, 2, c(num_center_point[1], 
                                                             num_center_point[2]), "+")
  
  # Get sf polygon representing the ellipse
  sf_ellipse <- st_polygon(list(matrix_ellipse_points))
  sfc_ellipse <- st_sfc(sf_ellipse, crs = st_crs(sf_network))
  
  # Get road segments lying within the ellipse
  sf_roads <- st_intersection(sf_network, sfc_ellipse)
  
  
  # After performing the intersection, geometries such as MULTILINESTRING are removed,...
  #...which is totally fine since there are only few cases where MULTILINESTRINGS...
  #...
  
  sf_roads <- sf_roads %>%
    filter(st_geometry_type(geometry) == "LINESTRING")
  
  # Ensure that linestrings are connected so network distance between points is low
  sfs_net_small <- as_sfnetwork(sf_roads, directed = FALSE)
  # Get the components of the created network
  list_components <- tidygraph::to_components(sfs_net_small)
  list_components_lengths <- lapply(list_components, function(i){length(i)})
  # Choose the largest component...
  num_largest_comp <- which.max(unlist(list_components_lengths))
  sfs_net_small <- list_components[[num_largest_comp]]
  
  #...and use the corresponding road segments...
  sf_roads <- sfs_net_small %>%
    activate(edges) %>%
    as.data.frame() %>%
    st_as_sf()
  st_crs(sf_roads) <- st_crs(sf_network)
  
  
  
  
  #...to sample points
  sfc_points <- st_sample(sf_roads, size = int_n, type = "random")
  sfc_points <- sfc_points[(!st_is_empty(sfc_points)) %>% which]
  sfc_points <- st_cast(sfc_points, "POINT") 
  st_crs(sfc_points) <- st_crs(sfc_points)
  sf_points <- sfc_points %>% st_as_sf()
  sf_points$point_id <- 1:nrow(sf_points)
  sf_points <- helper_get_line_ids_of_mapped_points(sf_points, sf_roads)
  sf_points <- sf_points %>% st_as_sf()
  return(sf_points)
}


# Documentation: main_get_avg_radius_chull
# Usage: main_get_avg_radius_chull(sf_convex_hull)
# Description: Calculates distances between all edge-vertices of the networks 
# convex hull and returns the average distance to ensure properly scaled
# flow distances
# Args/Options: sf_convex_hull
# Returns: numeric
# Output: ...
# Action: ...
main_get_avg_radius_chull <- function(sf_convex_hull) {
  sf_centroid <- st_centroid(sf_convex_hull)
  sfc_boundary_points <- st_cast(sf_convex_hull, "POINT")  
  unit_distances <- st_distance(sf_centroid, sfc_boundary_points)
  num_average_radius <- mean(as.numeric(unit_distances))
  return(num_average_radius)
}





# Documentation: connect_od_points_to_linestrings
# Usage: connect_od_points_to_linestrings(sf_origin, sf_dest)
# Description: Creates sf-linestrings from sf-od-point-pairs
# Args/Options: sf_origin, sf_dest
# Returns: sf-dataframe
# Output: ...
# Action: ...
connect_od_points_to_linestrings <- function(sf_origin, sf_dest) {
  matrix_dist_mat <- as.matrix(st_distance(sf_origin, sf_dest))
  matrix_dist_mat <- drop_units(matrix_dist_mat)
  matrix_assignment <- solve_LSAP(matrix_dist_mat) %>% as.numeric()
  # Ergebnisse extrahieren
  df_results <- data.frame(
    start_id = 1:nrow(sf_origin),
    end_id = matrix_assignment,
    distance = matrix_dist_mat[cbind(1:nrow(sf_origin), matrix_assignment)]
  )
  # Liniengeometrien erstellen
  sfc_lines <- do.call(rbind, lapply(1:nrow(df_results), function(i) {
    st_sfc(st_linestring(rbind(
      st_coordinates(sf_origin[df_results$start_id[i], ]),
      st_coordinates(sf_dest[df_results$end_id[i], ])
    )), crs = 4326)
  }))
  
  sf_lines <- st_sf(df_results, geometry = sfc_lines)
  return(sf_lines)
}



get_endpoint_within_dist <- function(sf_network, sf_startpoint, dt_dist_mat, int_max_dist) {
  
  sf_network$dist_to_start <- st_distance(sf_network, st_sfc(sf_startpoint,
                                                             crs = st_crs(sf_network))) %>%
    as.numeric()
  
  int_source_node <- sf_network$source[which.min(sf_network$dist_to_start)]
  int_target_node <- dt_dist_mat[source==int_source_node][which.min(abs(m - int_max_dist))] %>%
    pull(target)
  
  sf_target_row <- sf_network %>% filter(target == int_target_node) %>% head(1)
  sf_source_row <- sf_network %>% filter(source == int_target_node) %>% head(1)
  
  sfc_endpoint <- st_sfc()  
  
  # Überprüfe, welche Bedingung zutrifft
  if (nrow(sf_target_row) > 0) {
    sfc_endpoint <- lwgeom::st_endpoint(sf_target_row$geometry)
  } else if (nrow(sf_source_row) > 0) {
    sfc_endpoint <- lwgeom::st_startpoint(sf_source_row$geometry)
  } 
  
  return(sfc_endpoint)
}


# Documentation: main_create_synth_dataset
# Usage: main_create_synth_dataset(int_n_clusters, num_avg_radius, int_min_cl_size,
# int_max_cl_size, num_var_angle_direction)
# Description: Main function that creates int_n clusters on a given road network
# Args/Options: int_n_clusters, num_avg_radius, int_min_cl_size,
# int_max_cl_size, num_var_angle_direction
# Returns: sf-dataframe
# Output: ...
# Action: ...
main_create_synth_dataset <- function(int_n_clusters, 
                                      num_avg_radius,
                                      int_min_cl_size,
                                      int_max_cl_size,
                                      num_var_angle_direction,
                                      num_min_length,
                                      num_max_length) {
  
  
  
  list_all_lines <- list()
  list_all_origin_pts <- list()
  list_all_dest_pts <- list()
  list_all_cluster_ids <- list()
  for (i in 1:int_n_clusters) {
    
    # Configuration for LINESTRINGS
    int_n_points <- sample(int_min_cl_size:int_max_cl_size, 1)
    # Configuration for ellipse
    num_width_min <- num_avg_radius/8
    num_width_m <- sample(num_width_min:num_width_min*3, 1)
    num_height_m <- num_width_min
    
    # Get random road within the network
    sf_sampled_road <- sample_n(sf_network, 1)
    num_center_point <-  st_coordinates(st_sample(sf_sampled_road, 1)) %>% as.numeric
    
    # Sample origin points
    sf_points_o_cl <- create_ellipse_with_points(sf_network = sf_network,
                                                 num_center_point = num_center_point,
                                                 num_width_m = num_width_m,
                                                 num_height_m = num_height_m,
                                                 int_n = int_n_points)
    
    sf_points_o_cl <- sf_points_o_cl %>%
      rename("origin_point_id" = "point_id",
             "origin_id" = "id",
             "origin_source" = "source",
             "origin_target"="target",
             "origin_geom"="geom")
    
    
    # Get one of the origin points...
    sf_startpoint <- sample(sf_points_o_cl$origin_geom, 1)
    num_length <- sample(num_min_length:num_max_length, 1)

    sf_endpoint <- get_endpoint_within_dist(sf_network, sf_startpoint, dt_dist_mat, num_length)
    # Sample destination points
    sf_points_d_cl <- create_ellipse_with_points(sf_network = sf_network,
                                                 num_center_point = st_coordinates(
                                                   sf_endpoint),
                                                 num_width_m = num_width_m,
                                                 num_height_m = num_height_m,
                                                 int_n = int_n_points)  
    sf_points_d_cl <- sf_points_d_cl %>%
      rename("dest_point_id" = "point_id",
             "dest_id" = "id",
             "dest_source" = "source",
             "dest_target"="target",
             "dest_geom"="geom")
    
    # Draw LINESTRINGS with little intersection between each other
    sf_linestrings_cl <- connect_od_points_to_linestrings(sf_origin = sf_points_o_cl,
                                                          sf_dest = sf_points_d_cl)
    
    # Save created OD points and corresponding OD flows in lists...
    list_all_lines[[i]] <- sf_linestrings_cl
    list_all_origin_pts[[i]] <- sf_points_o_cl
    list_all_dest_pts[[i]] <- sf_points_d_cl
    list_all_cluster_ids[[i]] <- rep(i, int_n_points)
    cat("Cluster", i, "created\n")
  }
  
  # ... and combine them into one matrix
  c_combined_lines <- do.call(rbind, list_all_lines)
  c_combined_origin_pts <- do.call(rbind, list_all_origin_pts)
  c_combined_dest_pts <- do.call(rbind, list_all_dest_pts)
  c_cluster_id <- unlist(list_all_cluster_ids) 
  matrix_cluster_data <- cbind(c_combined_lines, 
                               c_cluster_id, 
                               c_combined_origin_pts, 
                               c_combined_dest_pts)
  
  
  sf_cluster_data <- st_as_sf(data.frame(matrix_cluster_data), 
                              sf_column_name = "geometry")
  st_crs(sf_cluster_data) <- st_crs(sf_network)
  sf_cluster_data$cluster_id <- unlist(sf_cluster_data$c_cluster_id) 
  sf_cluster_data <- sf_cluster_data %>%
    select(cluster_id, origin_id, dest_id, geometry)
  return(sf_cluster_data)
 }