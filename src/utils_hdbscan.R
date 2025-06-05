# Documentation: helper_dist_to_start_end
# Usage: helper_dist_to_start_end(sf_trips, dt_network, char_road_id, char_point_geom,
# char_suffix)
# Description: Joins 'sf_trips' with 'dt_network' based on 'char_road_id' (indicating
# whether the distance from the origin or destination of the flow to start and
# end of the udnerlying road segment should be calculated)
# Args/Options: sf_trips, dt_network, char_road_id, char_point_geom, char_suffix
# Returns: sf-dataframe
# Output: ...
# Action: ...
helper_dist_to_start_end <- function(sf_trips, 
                                     dt_network,
                                     char_road_id, 
                                     char_point_geom, 
                                     char_suffix) {
  dt_network <- dt_network %>%
    rename(geometry_road = geometry) %>%
    mutate(length_road = as.numeric(st_length(geometry_road))) %>%
    select(id, all_of(c("geom_start", "geom_end")), geometry_road, length_road)
  
  sf_trips <- sf_trips %>%
    left_join(dt_network, by = setNames("id", char_road_id)) %>%
    mutate(
      dist_to_start = st_distance(sf_trips[[char_point_geom]], geom_start, by_element = TRUE) %>% as.numeric(),
      dist_to_end = st_distance(sf_trips[[char_point_geom]], geom_end, by_element = TRUE) %>% as.numeric()
    ) %>%
    rename_with(~ paste0(., "_", char_suffix), c(length_road, dist_to_start, dist_to_end)) %>%
    select(-geom_start, -geom_end, -geometry_road)
  
  return(sf_trips)
}

# Documentation: add_dist_start_end
# Usage: add_dist_start_end(sf_trips, sf_network)
# Description: Calculates distance of OD-points to start and end-node of the
# corresponding road segment
# Args/Options: sf_trips, sf_network
# Returns: data.table
# Output: ...
# Action: ...
add_dist_start_end <- function(sf_trips, sf_network) {
  
  # Preprocess network: Add start/end points and convert to data.table
  sf_network <- sf_network %>%
    mutate(
      geom_start = lwgeom::st_startpoint(geometry),
      geom_end = lwgeom::st_endpoint(geometry)
    )
  
  dt_network <- sf_network %>% 
    as.data.table() %>%
    select(id, source, target, geom_start, geom_end, geometry)
  
  # Precompute trip points
  sf_trips <- sf_trips %>%
    mutate(
      geom_origin = lwgeom::st_startpoint(geometry),
      geom_dest = lwgeom::st_endpoint(geometry)
    )
  
  
  
  # Process both sides
  sf_trips <- sf_trips %>%
    helper_dist_to_start_end(dt_network = dt_network,
                             char_road_id = "origin_id", 
                             char_point_geom = "geom_origin", 
                             char_suffix = "origin") %>%
    helper_dist_to_start_end(dt_network = dt_network,
                             char_road_id = "dest_id", 
                             char_point_geom = "geom_dest", 
                             char_suffix = "dest") 
  # Final cleanup
  dt_trips <- sf_trips %>%
    st_drop_geometry() %>%
    as.data.table() %>%
    select(
      cluster_id, origin_id, dest_id, flow_id,
      length_road_origin, dist_to_start_origin, dist_to_end_origin,
      length_road_dest, dist_to_start_dest, dist_to_end_dest
    )
  
  return(dt_trips)
}