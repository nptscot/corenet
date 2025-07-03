# Declare global variables to avoid 'no visible binding' notes
utils::globalVariables(c("edge_paths", "influence_network", "all_fastest_bicycle_go_dutch", 
                         "weight", "to_linegraph", "edges", "group", "mean_potential", "LAD23NM", 
                         "road_function",  "grid_id", "density", 
                         "max_value", "min_value", "arterialness", "road_score", "value", "key_attribute", "n_group",".data", "n_removeDangles", "path_type", "maxDistPts", "minDistPts","penalty_value","penalty", "use_stplanr", "group_column", "path_cache_env"))

#' Prepare a cohesive cycling network using NPT data
#'
#' This function prepares the base network for generating a cohesive cycling network 
#' using NPT data. It integrates two different road network data sources within a specified 
#' target zone to produce a cohesive cycling network that considers specific road attributes.
#'
#' @param base_network Base road network object from Open Road, class 'sf'.
#' @param influence_network NPT road network object, which contains influence factors 
#'        like all_fastest_bicycle_go_dutch, class 'sf'.
#' @param target_zone Target zone, a polygon of a study area, such as a 3km radius around 
#'        Edinburgh city centre, class 'sf'.
#' @param crs Coordinate reference system to use, default is "EPSG:27700".
#' @param key_attribute The attribute in the network data to filter by and influence the outcome.
#' @param attribute_values Values of the key_attribute to retain in the network.
#' @param use_stplanr Logical value indicating whether to use the stplanr package for merging
#' 
#' @return A list containing two elements: a cohesive network and zone data, both class 'sf'.
#' @export
#' @examples
#' library(sf)
#' library(zonebuilder)
#' library(dplyr)
#' library(tmaptools)
#'
#' # Load demo data 
#' os_edinburgh_demo_3km = sf::st_set_crs(os_edinburgh_demo_3km, 27700)
#' NPT_demo_3km = sf::st_set_crs(NPT_demo_3km, 27700)
#' base_network = sf::st_transform(os_edinburgh_demo_3km, crs = 27700)
#' influence_network = sf::st_transform(NPT_demo_3km, crs = 27700)
#' target_zone = zonebuilder::zb_zone("Edinburgh", n_circles = 3) |> sf::st_transform(crs = "EPSG:27700")        
#'
#' # Prepare the cohesive network
#' OS_NPT_demo = cohesive_network_prep( base_network = base_network, 
#'                                 influence_network = influence_network, 
#'                                 target_zone = target_zone, 
#'                                 key_attribute = "road_function", 
#'                                 crs = "EPSG:27700", 
#'                                 attribute_values = c("A Road", "B Road", "Minor Road"))
#'


cohesive_network_prep = function(base_network, influence_network, target_zone, crs = "EPSG:27700",  key_attribute = "road_function", attribute_values = c("A Road", "B Road", "Minor Road"), use_stplanr = FALSE) {
    base_network = sf::st_transform(base_network, crs)    
    influence_network = sf::st_transform(influence_network, crs)
    target_zone = sf::st_transform(target_zone, crs)

    # Intersect and transform NPT data with the calculated zones
    NPT_zones = influence_network[sf::st_union(target_zone), , op = sf::st_intersects] |> 
                sf::st_transform(crs)

    OS_zones = base_network[sf::st_union(target_zone), , op = sf::st_intersects]

    # Dynamically filter the data using provided column name and attributes
    filtered_OS_zones = OS_zones |> 
                        dplyr::filter(!!rlang::sym(key_attribute) %in% attribute_values ) |> 
                        sf::st_transform(crs) |> 
                        sf::st_zm()
    # if (use_stplanr) {
    #   # Assign functions for data aggregation based on network attribute values 
    #   name_list = names(NPT_zones)
      
    #   funs = list()
    #   for (name in name_list) {
    #     if (name == "geometry") {
    #       next  # Correctly skip the current iteration if the name is "geometry"
    #     } else if (name %in% c("gradient", "quietness")) {
    #       funs[[name]] = mean  # Assign mean function for specified fields
    #     } else {
    #       funs[[name]] = mean  # Assign sum function for all other fields
    #     }
    #   }
      
    #   # Merge road networks with specified parameters
    #   filtered_OS_NPT_zones = stplanr::rnet_merge(filtered_OS_zones, NPT_zones, dist = 10, funs = funs, segment_length = 20, max_angle_diff = 10)
      
    # } else {
    #   print("Using sf::st_join")
    #   filtered_OS_NPT_zones = sf::st_join(filtered_OS_zones, NPT_zones, join = sf::st_intersects)
    # }
    
    NPT_zones = sf::st_cast(NPT_zones, "LINESTRING")
    filtered_OS_zones = sf::st_cast(filtered_OS_zones, "LINESTRING")
    NPT_zones$id = 1:nrow(NPT_zones)
    filtered_OS_zones$id = 1:nrow(filtered_OS_zones)

    params = list(
                  list(
                    source = NPT_zones,
                    target = filtered_OS_zones,
                    attribute = "all_fastest_bicycle_go_dutch",
                    new_name = "all_fastest_bicycle_go_dutch",
                    agg_fun = sum,
                    weights = c("target_weighted")
                  )
                )

    results_list = purrr::map(params, function(p) {
      anime_join(
        source_data = p$source,
        target_data = p$target,
        attribute = p$attribute,
        new_name = p$new_name,
        agg_fun = p$agg_fun,
        weights = p$weights,
        angle_tolerance = 35,
        distance_tolerance = 15
      )
    })


    filtered_OS_NPT_zones = reduce(results_list, function(x, y) {
      left_join(x, y, by = "id")
    }, .init = filtered_OS_zones)

    filtered_OS_NPT_zones = filtered_OS_NPT_zones |>
                            dplyr::mutate(dplyr::across(dplyr::where(is.numeric), as.integer))

    print("Finished preparing the network data")
    
    return(filtered_OS_NPT_zones)
}


#' Generating cohesive cycling network
#'
#' This function processes the provided network data to identify and analyze 
#' significant routes based on a percentile threshold. It performs spatial 
#' operations, clusters data points, and calculates the largest network component 
#' without dangles, aiming to develop a cohesive cycling network infrastructure.
#'
#' @param influence_network The NPT network data, class 'sf'.
#' @param cohesive_base_network Spatial object representing the cohesive base network obtained via function cohesive_network_prep, class 'sf'.
#' @param target_zone Spatial object representing the study area or target zone, class 'sf'.
#' @param key_attribute Attribute used to determine significant network routes, default is "all_fastest_bicycle_go_dutch".
#' @param crs Coordinate reference system for transformation, default is "EPSG:27700".
#' @param maxDistPts Distance threshold used in path calculations, default is 1500 meters.
#' @param minDistPts Minimum distance threshold used in path calculations, default is 2 meters.
#' @param npt_threshold Threshold value for filtering the NPT network, default is 1500.
#' @param road_scores A list of road types and their corresponding scoring weights.
#' @param n_removeDangles Number of iterations to remove dangles from the network, default is 6.
#' @param penalty_value The penalty value for roads with low values, default is 1.
#' @param group_column The column name to group the network by edge betweenness, default is "name_1".
#' @return A spatial object representing the largest cohesive component of the network, free of dangles.
#' @export
#' @examples
#' library(sf)
#' library(dplyr)
#' library(dbscan)
#' library(zonebuilder)
#'
#' # Load demo data 
#' os_edinburgh_demo_3km = sf::st_set_crs(os_edinburgh_demo_3km, 27700)
#' NPT_demo_3km = sf::st_set_crs(NPT_demo_3km, 27700)
#' base_network = sf::st_transform(os_edinburgh_demo_3km, crs = 27700)
#' influence_network = sf::st_transform(NPT_demo_3km, crs = 27700)
#' target_zone = zonebuilder::zb_zone("Edinburgh", n_circles = 2) |>
#'                sf::st_transform(crs = "EPSG:27700")
#'
#' # Execute the function
#' OS_NPT_demo = cohesive_network_prep( base_network = base_network, 
#'                                 influence_network = influence_network, 
#'                                 target_zone = target_zone, 
#'                                 key_attribute = "road_function", 
#'                                 crs = "EPSG:27700", 
#'                                 attribute_values = c("A Road", "B Road", "Minor Road"))
#' OS_NPT_demo$geometry = OS_NPT_demo$geom
#' coherent_network = corenet(influence_network = OS_NPT_demo, 
#'                   cohesive_base_network = OS_NPT_demo, 
#'                   target_zone = target_zone, 
#'                   key_attribute = "all_fastest_bicycle_go_dutch", 
#'                   crs = "EPSG:27700",npt_threshold = 1500,maxDistPts = 1500,minDistPts = 2,
#'                   road_scores = list("A Road" = 1, "B Road" = 1, "Minor Road" = 1000),n_removeDangles = 6,penalty_value = 1)
#'

corenet = function(influence_network, cohesive_base_network, target_zone, key_attribute = "all_fastest_bicycle_go_dutch",  crs = "EPSG:27700", npt_threshold = 1500, maxDistPts = 1500,minDistPts = 2, road_scores = list("A Road" = 1, "B Road" = 1, "Minor Road" = 1000), n_removeDangles = 6 , penalty_value = 1, group_column = "name_1", max_path_weight =10) {

  if (key_attribute %in% names(influence_network)) {
    paste0("Using ", key_attribute, " as indicator for the network")
    influence_network = sf::st_transform(influence_network, crs)
    influence_network_split = stplanr::line_segment(influence_network, segment_length = 20, use_rsgeo = TRUE)

    # Assuming threshold is already defined somewhere in your code
    filtered_roads = influence_network_split[influence_network_split[[key_attribute]] > npt_threshold, ]

    # Calculate centroids of network segments
    centroids = sf::st_centroid(filtered_roads)

    # Perform DBSCAN clustering
    coordinates = sf::st_coordinates(centroids)
    coordinates_clean = coordinates[complete.cases(coordinates), ]
    clusters = dbscan::dbscan(coordinates_clean, eps = 18, minPts = 1)
    centroids_clean = centroids[complete.cases(coordinates), ]
    centroids_clean$cluster = clusters$cluster
    unique_centroids = centroids_clean[!duplicated(centroids_clean$cluster), ]   

    # create a buffer of 10 meters around the cohesive_base_network
    cohesive_base_network_buffer = sf::st_buffer(cohesive_base_network, dist = 10)

    # filter unique_centroids by the buffer
    unique_centroids = sf::st_intersection(unique_centroids, cohesive_base_network_buffer)

    if (!is.null(group_column)) {
      tryCatch({
        if (group_column %in% colnames(unique_centroids)) {
          # Ensure point_type column exists
          unique_centroids$point_type = "original_point"

          unique_centroids = unique_centroids |>
            group_by(!!sym(group_column)) |>
            group_modify(~ {
              # Compute the centroid of the group
              group_centroid = sf::st_union(.x$geometry) |> sf::st_centroid()
              centroid_row = .x[1, ]
              centroid_row$geometry = group_centroid
              centroid_row$point_type = "centroid"

              if (nrow(.x) < 2) {
                # For groups with less than 2 points, return the point and centroid
                # Use rbind to preserve sf class
                return(rbind(.x, centroid_row))
              } else {
                # Compute the distance matrix
                dist_matrix = st_distance(.x)
                # Exclude self-distances by setting diagonal to NA
                diag(dist_matrix) = NA
                # Find the maximum distance
                max_dist = max(dist_matrix, na.rm = TRUE)
                # Get the indices of the two points with maximum distance
                indices = which(dist_matrix == max_dist, arr.ind = TRUE)[1, ]
                # Extract the two points
                max_points = .x[indices, ]
                max_points$point_type = "max_distance_point"
                # Combine the two points with the centroid
                result = rbind(max_points, centroid_row)
                return(result)
              }
            }) |>
            ungroup()

          # Ensure the result is an sf object
          unique_centroids = sf::st_as_sf(unique_centroids) |> select(geometry)
        } else {
          message(sprintf("Column '%s' does not exist in unique_centroids.", group_column))
        }
      }, error = function(e) {
        message("An error occurred: ", e$message)
      })
    }
  } else {
    warning(paste0("Warning: ", key_attribute, " does not exist in the network"))
    # grid_sf = target_zone

    # grid_sf$grid_id = seq_len(nrow(grid_sf))

    # split_network = sf::st_intersection(influence_network, grid_sf)

    # select_by_density = function(density) {
    #   if (density < 10) {
    #     top_n = 1
    #   } else if (density >= 10 & density < 20) {
    #     top_n = 5
    #   } else {
    #     top_n = (density %/% 10) + 12
    #   }
      
    #   # Check if the value is NA (rest of the cases) and set top_n to 1
    #   if (is.na(density)) {
    #     top_n = 1
    #   }
    #   return(top_n)
    # }

    # distance_threshold = units::set_units(dist, "m") # Adjust the distance as needed

    # split_network = split_network |>
    #   dplyr::group_by(grid_id) |>
    #   dplyr::mutate(top_n = select_by_density(unique(density)))


    # split_network_max = split_network |>
    #   dplyr::group_by(grid_id) |>
    #   dplyr::group_modify(~ {
    #     # Step 1: Filter for the current grid_id
    #     split_network_filtered = .x
        
    #     # Step 2: Calculate Centroids for These Elements
    #     centroids = sf::st_centroid(split_network_filtered)
        
    #     # Step 3: Sort the Centroids Based on a Value Attribute
    #     sorted_centroids = centroids |>
    #       dplyr::arrange(desc(value))
        
    #     # Step 4: Identify the Centroid with the Highest Value
    #     highest_value_centroid = sorted_centroids[1, ]
        
    #     # Step 5: Calculate Distances from Each Centroid to the Highest Value Centroid
    #     distances_to_highest = sf::st_distance(centroids, highest_value_centroid)

    #     # Step 6: Select Line Strings Based on top_n and Distance Criteria
    #     top_n = unique(split_network_filtered$top_n)

    #     selected_lines = split_network_filtered |>
    #       dplyr::mutate(distance_to_highest = sf::st_distance(centroids, highest_value_centroid)[,1]) |>
    #       dplyr::filter(distance_to_highest >= distance_threshold) |>
    #       dplyr::slice_head(n = top_n)        
    #     return(selected_lines)
    #   }) |> 
    #   dplyr::bind_rows()

    # # Calculate centroids of these max segments
    # split_network_max = sf::st_as_sf(split_network_max)
    # unique_centroids = sf::st_centroid(split_network_max)
  }

  # Prepare network and calculate paths
  prepared_network = prepare_network(cohesive_base_network, key_attribute , road_scores, transform_crs = crs , penalty_value = penalty_value) 

  # Initialize cache environment for path calculations
  path_cache_env = new.env(parent = emptyenv())

  all_paths = purrr::map_dfr(
      seq_len(nrow(unique_centroids)),
      function(n) {
          calculate_paths_from_point_dist(prepared_network, point = unique_centroids[n,], centroids = unique_centroids, path_type = "shortest", maxDistPts = maxDistPts, minDistPts = minDistPts, max_path_weight = max_path_weight)
      }
  )

  # Remove duplicates and calculate the largest connected component
  all_paths_sf_percentile_duplicated_geoms = duplicated(all_paths$geometry)
  all_paths_sf_percentile = all_paths[!all_paths_sf_percentile_duplicated_geoms, ]
  largest_component_sf = calculate_largest_component(all_paths_sf_percentile) |> sf::st_as_sf('edges')
  largest_component_sf_without_dangles = removeDangles(largest_component_sf, tolerance = 0.001)

  # Remove dangles multiple times to ensure a clean network
  for (i in 1:n_removeDangles) {
    largest_component_sf_without_dangles = removeDangles(largest_component_sf_without_dangles)
  }

  return(largest_component_sf_without_dangles)
}


#' Generate a grouped network by edge betweenness from a preprocessed coherent network
#'
#' This function takes a preprocessed network and applies graph-based analysis to group
#' the network edges based on their betweenness centrality. The function assumes that
#' the input network has attributes relevant to cycling traffic dynamics, specifically
#' 'all_fastest_bicycle_go_dutch' and 'weight'. It outputs a transformed network where
#' edges are grouped and ranked according to their mean potential, facilitating further
#' analysis or visualization of critical network pathways.
#' @param coherent_network A preprocessed 'sf' object containing the network data,
#'        expected to have columns 'all_fastest_bicycle_go_dutch' and 'weight'.
#' @param key_attribute The attribute used to keep.
#' @param n_group The number of groups to divide the network into, based on edge
#'        betweenness centrality.
#' @return An 'sf' object with edges grouped and ranked based on their mean potential.
#' @export
#' @examples
#' # Assuming 'coherent_network' is obtained from a previous function corenet
#' # Generate the grouped network
#' # grouped_network = coherent_network_group(coherent_network)
coherent_network_group = function(coherent_network, key_attribute = "all_fastest_bicycle_go_dutch", n_group = 12) {

  rnet_coherent_selected = coherent_network |>
    dplyr::select({{ key_attribute }})
  
  # Group and process the network
  grouped_net = rnet_coherent_selected |>
    sfnetworks::as_sfnetwork(directed = FALSE) |>
    tidygraph::morph(tidygraph::to_linegraph) |>
    dplyr::mutate(group = tidygraph::group_edge_betweenness(n_groups = n_group)) |>
    tidygraph::unmorph() |>
    tidygraph::activate(edges) |>
    sf::st_as_sf() |>
    sf::st_transform("EPSG:4326") |>
    dplyr::group_by(group) |>
    dplyr::summarise(mean_potential  = sum(get(key_attribute), na.rm = TRUE)) |>
    dplyr::mutate(group  = rank(-mean_potential)) 

  return(grouped_net)
}



#' Prepare a network data structure by transforming, scoring, and weighting based on road types and conditions
#'
#' This function transforms a spatial network object into an 'sfnetwork', scoring it based on road conditions and classifications.
#' The transformation process includes casting the network to LINESTRING, converting it to 'sfnetwork', and adding attributes
#' like 'arterialness' and 'weight' which are calculated based on the given road scores and the importance of the roads in the network.
#'
#' @param network An 'sf' object representing a spatial network.
#' @param key_attribute The attribute name from the network used for normalization and scoring, default is "all_fastest_bicycle_go_dutch".
#' @param road_scores A list specifying scores for different road types. Example: list("A Road" = 1, "B Road" = 1, "Minor Road" = 10000000).
#' @param transform_crs The numeric CRS code for coordinate transformation, default is 27700.
#' @param penalty_value The penalty value for roads with low values, default is 1.
#' @return An 'sfnetwork' object with attributes 'arterialness' and 'weight' that account for road conditions and their relative importance.
#' @examples
#' library(sf)
#' library(sfnetworks)
#' library(dplyr)
#' library(purrr)
#' library(zonebuilder)
#'
#' # Load demo data 
#' os_edinburgh_demo_3km = sf::st_set_crs(os_edinburgh_demo_3km, 27700)
#' NPT_demo_3km = sf::st_set_crs(NPT_demo_3km, 27700)
#' base_network = sf::st_transform(os_edinburgh_demo_3km, crs = 27700)
#' influence_network = sf::st_transform(NPT_demo_3km, crs = 27700)
#' target_zone = zonebuilder::zb_zone("Edinburgh", n_circles = 2) |>
#'                sf::st_transform(crs = "EPSG:27700")
#'
#' # Prepare the cohesive network
#' network = cohesive_network_prep(base_network = base_network, 
#'                                 influence_network = influence_network, 
#'                                 target_zone = target_zone, 
#'                                 key_attribute = "road_function", 
#'                                 attribute_values = c("A Road", "B Road", "Minor Road"))
#'
#' # Define road scores
#' road_scores = list("A Road" = 1, "B Road" = 1, "Minor Road" = 100000)
#'
#' # Prepare the network
#' prepared_network = prepare_network(network, 
#'                                     key_attribute = "all_fastest_bicycle_go_dutch",
#'                                     road_scores = road_scores,
#'                                     transform_crs = 27700, penalty_value = 1)
#'
#' # Print the prepared network
#' print(prepared_network)
#' @export

prepare_network = function(network, key_attribute = "all_fastest_bicycle_go_dutch", 
                           road_scores = list("A Road" = 1, "B Road" = 2, "Minor Road" = 100000),penalty_value = 1,
                           transform_crs = 27700) {
    # Cast the network to LINESTRING and transform to the specified CRS
    network = network |>
        sf::st_cast("LINESTRING") |>
        sf::st_transform(transform_crs) |>
        sfnetworks::as_sfnetwork(directed = FALSE) |>
        sfnetworks::activate("edges") |>
        dplyr::mutate(
            # Handle NA values and normalize using key_attribute
            value = dplyr::if_else(is.na(!!rlang::sym(key_attribute)), 0, !!rlang::sym(key_attribute)),
            max_value = max(value, na.rm = TRUE),
            value = ifelse(value <= 50, 0, value),
            value = ifelse(value == 0, 0.01, value),
            # Calculate logarithmic arterialness and round it
            arterialness = ((max_value / value) ^ 3) / max(((max_value / value) ^ 3), na.rm = TRUE),
            # Dynamically apply road type scores from the road_scores list
            road_score = purrr::map_dbl(road_function, function(x) {
                if (x %in% names(road_scores)) {
                    road_scores[[x]]
                } else {
                    0  # Default score for unrecognized road types
                }
            }),
            # Calculate weight considering the road type influence
            weight = round(0.95 * arterialness + 0.05 * road_score, 6), 
            penalty = ifelse(value <= 200, penalty_value, 1),
            weight = weight * penalty    
        )
    # network = sfnetworks::activate(network, "nodes")
    return(network)
}


#' Calculate paths from a given point to centroids within a specified distance range
#'
#' This function determines the network paths from a specific point to multiple centroids based on distance thresholds.
#' It can compute the shortest, all shortest, or all simple paths depending on the path type specified.
#'
#' @param network An sfnetwork object representing the network.
#' @param point An sf or sfc object containing a single point feature from which paths will be calculated.
#' @param maxDistPts The maximum distance (in meters) to consider for path calculations.
#' @param minDistPts The minimum distance (in meters) to consider for path calculations.
#' @param centroids An sf object containing centroids to which paths are calculated.
#' @param path_type A character string indicating the type of path calculation: 'shortest', 'all_shortest', or 'all_simple'.
#' @param max_path_weight 
#' 
#'        - 'shortest': Calculates the shortest path considering weights.
#'        - 'all_shortest': Calculates all paths that tie for the shortest distance, considering weights.
#'        - 'all_simple': Calculates all simple paths, ignoring weights.
#' @return An sf object containing the paths that meet the criteria or NULL if no paths meet the criteria.
#' @export

calculate_paths_from_point_dist = function(
  network,
  point,
  minDistPts = 2,
  maxDistPts = 1500,
  centroids,
  path_type = "shortest",
  max_path_weight = 10,
  crs_transform = 27700
) {
  # 1) Ensure the network is in a meter-based CRS, transform if necessary
  net_crs = sf::st_crs(network)
  
  # Check if units are in "m"; if missing or different, transform
  if (is.na(net_crs) || (!is.null(net_crs$units) && net_crs$units != "m")) {
    message("Transforming network to CRS: ", crs_transform)
    network = sf::st_transform(network, crs = crs_transform)
  }
  
  # 2) Generate a unique key for caching based on the point's coordinates
  #    You can modify how this key is generated if your 'point' is already sf-compatible
  point_key = paste(sort(as.character(point)), collapse = "_")
  
  # 3) Check if result is already cached
  if (exists(point_key, envir = path_cache_env)) {
    return(get(point_key, envir = path_cache_env))
  }
  
  # 4) Convert the point and centroids to sfc if not already
  point_geom     = sf::st_as_sfc(point)
  centroids_geom = sf::st_as_sfc(centroids)
  
  # 5) Filter centroids by distance to the point
  distances = sf::st_distance(point_geom, centroids_geom)
  # distances is a 1 x n matrix; convert to a numeric vector (in meters)
  distances_m = units::set_units(as.vector(distances[1, ]), "m")
  
  valid_idx = which(
    distances_m >= units::set_units(minDistPts, "m") &
      distances_m <= units::set_units(maxDistPts, "m")
  )
  
  # If no centroids qualify, cache NULL and return
  if (!length(valid_idx)) {
    assign(point_key, NULL, envir = path_cache_env)
    return(NULL)
  }
  
  valid_centroids      = centroids[valid_idx, , drop = FALSE]
  valid_centroids_geom = centroids_geom[valid_idx]
  
  # 6) Determine which weight column to use
  #    If path_type == "all_simple", typically no weights are used.
  #    Adjust as needed if your network uses a different attribute for cost/distance.
  weights_to_use = if (path_type == "all_simple") NULL else "weight"
  
  # 7) Compute paths from the single point to all valid centroids
  #    This returns indices of nodes/edges for each path
  paths_from_point = sfnetworks::st_network_paths(
    network,
    from    = point_geom,
    to      = valid_centroids_geom,
    weights = weights_to_use, 
    type    = path_type
  )
  
  # 8) Summation of edge weights for each path
  #    Instead of re-activating edges for each path, extract them once
  net_edges        = network |> tidygraph::activate("edges")
  edges_weight_vec = net_edges |> dplyr::pull(weight)
  
  # Vectorized summation of weights for each path
  total_weights = sapply(paths_from_point$edge_paths, function(edge_ids) {
    sum(edges_weight_vec[edge_ids], na.rm = TRUE)
  })
  
  paths_from_point$total_weight = total_weights
  
  # 9) Filter out paths whose total weight exceeds threshold
  valid_paths = paths_from_point[total_weights <= max_path_weight, ]
  if (!nrow(valid_paths)) {
    assign(point_key, NULL, envir = path_cache_env)
    return(NULL)
  }
  
  # 10) Collect the edges used by these filtered paths
  edges_in_paths = unique(unlist(valid_paths$edge_paths))
  
  # 11) Slice the network edges to extract only those edges
  result = net_edges |> 
    dplyr::slice(edges_in_paths) |> 
    sf::st_as_sf()
  
  # 12) Store in cache and return
  assign(point_key, result, envir = path_cache_env)
  return(result)
}

#' Calculate the largest connected component of a network
#'
#' This function takes a spatial network represented as an `sf` object, converts it into a graph format,
#' identifies all connected components, and extracts the largest one.
#'
#' @param network An `sf` object representing the network to be analyzed.
#' @return An `sfnetwork` object representing the largest connected component of the network.
#' @export
#' 
calculate_largest_component = function(network) {
    # Ensure the network tile is in the correct format and convert to sfnetwork
    if (!inherits(network, "sfnetwork")) {
        network_sfn = sfnetworks::as_sfnetwork(network, directed = FALSE)
    } else {
        network_sfn = network
    }

    # Convert sfnetwork to tbl_graph for easier graph manipulations
    network_igraph = tidygraph::as_tbl_graph(network_sfn)

    # Compute connected components using igraph
    components = igraph::components(network_igraph)

    # Find the largest component by size
    largest_component_id = which.max(components$csize)
    largest_component_nodes = which(components$membership == largest_component_id)

    # Extract the subgraph corresponding to the largest component
    largest_component_subgraph = igraph::subgraph(network_igraph, largest_component_nodes)

    # Convert back to sfnetwork
    largest_component_sfn = sfnetworks::as_sfnetwork(largest_component_subgraph, directed = FALSE)

    return(largest_component_sfn)
}


#' Remove dangling segments from a road network
#'
#' This function identifies and removes dangling line segments (lines with one end not connected to another line)
#' within a specified tolerance.
#'
#' @param network An `sf` object of class LINESTRING representing the road network.
#' @param tolerance The distance tolerance for identifying isolated endpoints as dangling.
#' @return An `sf` object with dangling line segments removed.
#' @export
#' 
removeDangles = function(network, tolerance = 0.001) {
    # Convert to Spatial Lines if not already
    network_lines = sf::st_cast(network, "LINESTRING")

    # Extract and combine all end points of line segments
    end_points = do.call(rbind, lapply(network_lines$geometry, function(line) {
    endpoints = rbind(sf::st_coordinates(line)[1, ], sf::st_coordinates(line)[nrow(sf::st_coordinates(line)), ])
    sf::st_as_sf(data.frame(x = endpoints[, 1], y = endpoints[, 2]), coords = c("x", "y"), crs = sf::st_crs(network))
    }))

    # Identify unique end points (potential dangles) using a spatial join to find nearby points
    buffer_points = sf::st_buffer(end_points, dist = tolerance)
    overlaps = sf::st_intersects(buffer_points, buffer_points, sparse = FALSE)
    isolated_points = end_points[rowSums(overlaps) == 1,]

    # Filter out road segments that end in these isolated points
    segments_with_dangles = sapply(sf::st_geometry(network_lines), function(geom) {
    ends = sf::st_sfc(sf::st_point(sf::st_coordinates(geom)[1,]), sf::st_point(sf::st_coordinates(geom)[nrow(sf::st_coordinates(geom)),]), crs = sf::st_crs(network))
    any(sf::st_intersects(ends, isolated_points, sparse = FALSE))
    })

    network_without_dangles = network_lines[!segments_with_dangles, ]

    return(network_without_dangles)
}



#' Create coherent network PMtiles
#'
#' This function generates PMtiles for coherent network GeoJSON file using the Tippecanoe tool.
#' @param folder_path The directory path where the files will be saved.
#' @param city_filename The base name for the output files, in this case using city name.
#' @param cohesive_network An sf object representing the cohesive network.
#' @return The output from the Tippecanoe command as a character vector.
#' @export
#' 
create_coherent_network_PMtiles = function(folder_path, city_filename, cohesive_network) {
  # Check geometry types of the cohesive_network
  geometry_types = sf::st_geometry_type(cohesive_network)

  # Handle LINESTRING geometries
  cohesive_network_linestring = cohesive_network[geometry_types == "LINESTRING", ]

  # Handle MULTILINESTRING geometries
  cohesive_network_multilinestring = cohesive_network[geometry_types == "MULTILINESTRING", ]
  cohesive_network_multilinestring = sf::st_cast(cohesive_network_multilinestring, "LINESTRING")

  # Handle other geometries (e.g., POINT or POLYGON)
  cohesive_network_others = cohesive_network[geometry_types != "LINESTRING" & geometry_types != "MULTILINESTRING", ]

  # Combine LINESTRING and casted MULTILINESTRING geometries
  cohesive_network_combined = rbind(cohesive_network_linestring, cohesive_network_multilinestring)

  # Generate filenames for GeoJSON and pmtiles outputs
  coherent_geojson_filename = paste0(folder_path, city_filename, "_coherent_network", ".geojson")
  
  # Save cohesive_network_combined as GeoJSON
  sf::st_write(cohesive_network_combined, coherent_geojson_filename, delete_dsn = TRUE)
  
  coherent_pmtiles_filename = paste0(folder_path, city_filename, "_coherent_network", ".pmtiles")
  
  # Construct the Tippecanoe command
  command_tippecanoe = paste0(
    'tippecanoe -o ', coherent_pmtiles_filename,
    ' --name="', city_filename, '_coherent_network', '"',
    ' --layer=cohesivenetwork',
    ' --attribution="University of Leeds"',
    ' --minimum-zoom=6',
    ' --maximum-zoom=13',
    ' --maximum-tile-bytes=5000000',
    ' --simplification=10',
    ' --buffer=5',
    ' -rg',
    ' --force ',
    coherent_geojson_filename
  )
  
  # Execute the command and capture output
  system_output = system(command_tippecanoe, intern = TRUE)
  
  return(system_output)
}

