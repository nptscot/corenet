# Declare global variables to avoid 'no visible binding' notes
utils::globalVariables(c("edge_paths", "influence_network", "all_fastest_bicycle_go_dutch", 
                         "weight", "to_linegraph", "edges", "group", "mean_potential", "LAD23NM", 
                         "road_function",  "grid_id", "density", 
                         "max_value", "min_value", "arterialness", "road_score", "value", "key_attribute", "n_group",".data", "n_removeDangles"))

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
#' target_zone = zonebuilder::zb_zone("Edinburgh", n_circles = 2) |>
#'                sf::st_transform(crs = "EPSG:27700")
#'
#' # Prepare the cohesive network
#' OS_NPT_demo = cohesive_network_prep( base_network = base_network, 
#'                                 influence_network = influence_network, 
#'                                 target_zone = target_zone, 
#'                                 key_attribute = "road_function", 
#'                                 crs = "EPSG:27700", 
#'                                 attribute_values = c("A Road", "B Road", "Minor Road"))
#'


cohesive_network_prep = function(base_network, influence_network, target_zone, crs = "EPSG:27700",  key_attribute = "road_function", attribute_values = c("A Road", "B Road", "Minor Road")) {
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

    # Assign functions for data aggregation based on network attribute_values 
    name_list = names(NPT_zones)

    funs = list()
    for (name in name_list) {
      if (name == "geometry") {
        next  # Correctly skip the current iteration if the name is "geometry"
      } else if (name %in% c("gradient", "quietness")) {
        funs[[name]] = mean  # Assign mean function for specified fields
      } else {
        funs[[name]] = sum  # Assign sum function for all other fields
      }
    }

    # Merge road networks with specified parameters
    filtered_OS_NPT_zones = stplanr::rnet_merge(filtered_OS_zones, NPT_zones, dist = 10, funs = funs, segment_length = 20,max_angle_diff = 10)
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
#' @param dist Distance threshold used in path calculations, default is 10 meters.
#' @param threshold Value threshold for filtering significant routes, default is 1500.
#' @param road_scores A list of road types and their corresponding scoring weights.
#' @param n_removeDangles Number of iterations to remove dangles from the network, default is 6.
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
#'                   crs = "EPSG:27700", dist = 10, threshold = 1500,
#'                   road_scores = list("A Road" = 1, "B Road" = 1, "Minor Road" = 10000000))
#'

corenet = function(influence_network, cohesive_base_network, target_zone, key_attribute = "all_fastest_bicycle_go_dutch",  crs = "EPSG:27700", dist = 10, threshold = 1500, road_scores = list("A Road" = 1, "B Road" = 1, "Minor Road" = 10000000), n_removeDangles = 6) {

  if (key_attribute %in% names(influence_network)) {
    paste0("Using ", key_attribute, " as indicator for the network")
    influence_network = sf::st_transform(influence_network, crs)
    influence_network_split = stplanr::line_segment(influence_network, segment_length = 20, use_rsgeo = TRUE)

    # Assuming threshold is already defined somewhere in your code
    filtered_roads = influence_network_split[influence_network_split[[key_attribute]] > threshold, ]

    # Calculate centroids of network segments
    centroids = sf::st_centroid(filtered_roads)

    # Perform DBSCAN clustering
    coordinates = sf::st_coordinates(centroids)
    clusters = dbscan::dbscan(coordinates, eps = 18, minPts = 1)
    centroids$cluster = clusters$cluster
    unique_centroids = centroids[!duplicated(centroids$cluster), ]   

    # create a buffer of 10 meters around the cohesive_base_network
    cohesive_base_network_buffer = sf::st_buffer(cohesive_base_network, dist = 20)

    # filter unique_centroids by the buffer
    unique_centroids = sf::st_intersection(unique_centroids, cohesive_base_network_buffer)

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
  prepared_network = prepare_network(cohesive_base_network, key_attribute , road_scores, transform_crs = crs) 


  all_paths = purrr::map_dfr(
      seq_len(nrow(unique_centroids)),
      function(n) {
          calculate_paths_from_point_dist(prepared_network, point = unique_centroids[n,], centroids = unique_centroids, path_type = "shortest", dist = 1500)
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
#'
#' @param coherent_network A preprocessed 'sf' object containing the network data,
#'        expected to have columns 'all_fastest_bicycle_go_dutch' and 'weight'.
#' @param key_attribute The attribute used to keep.
#' @param n_group The number of groups to divide the network into, based on edge
#'        betweenness centrality.
#' 
#' @return An 'sf' object with edges grouped and ranked based on their mean potential.
#'
#' # Assuming 'coherent_network' is obtained from a previous function corenet
#' # Generate the grouped network
#' grouped_network = coherent_network_group(coherent_network)
#'
#' @export

coherent_network_group = function(coherent_network, key_attribute = "all_fastest_bicycle_go_dutch", n_group = 12) {

  rnet_coherent_selected = coherent_network |>
    dplyr::select({{ key_attribute }}, weight)
  
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
    dplyr::summarise(mean_potential = mean(weight, na.rm = TRUE)) |>
    dplyr::mutate(group = rank(-mean_potential))  

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
#'                                     transform_crs = 27700)
#'
#' # Print the prepared network
#' print(prepared_network)
#' @export

prepare_network = function(network, key_attribute = "all_fastest_bicycle_go_dutch", 
                           road_scores = list("A Road" = 1, "B Road" = 1, "Minor Road" = 10000000),
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
            min_value = min(value, na.rm = TRUE),
            max_value = max(value, na.rm = TRUE),
            arterialness = dplyr::if_else(
                max_value > min_value,
                (value - min_value) / (max_value - min_value),
                0  # Avoid division by zero
            ),
            # Dynamically apply road type scores from the road_scores list
            road_score = purrr::map_dbl(road_function, function(x) {
                if (x %in% names(road_scores)) {
                    road_scores[[x]]
                } else {
                    0  # Default score for unrecognized road types
                }
            }),
            # Calculate weight considering the road type influence
            weight = as.integer((1 - arterialness) * 100 * (1 + 0.1 * road_score))
        )

    return(network)
}


#' Calculate paths from a given point to centroids within a specified distance range
#'
#' This function determines the network paths from a specific point to multiple centroids based on distance thresholds.
#' It can compute the shortest, all shortest, or all simple paths depending on the path type specified.
#'
#' @param network An sfnetwork object representing the network.
#' @param point An sf or sfc object containing a single point feature from which paths will be calculated.
#' @param dist The maximum distance (in meters) to consider for path calculations.
#' @param centroids An sf object containing centroids to which paths are calculated.
#' @param path_type A character string indicating the type of path calculation: 'shortest', 'all_shortest', or 'all_simple'.
#'        - 'shortest': Calculates the shortest path considering weights.
#'        - 'all_shortest': Calculates all paths that tie for the shortest distance, considering weights.
#'        - 'all_simple': Calculates all simple paths, ignoring weights.
#' @return An sf object containing the paths that meet the criteria or NULL if no paths meet the criteria.
#' @export

calculate_paths_from_point_dist <- function(network, point, dist = 500, centroids, path_type = "shortest") {
    path_cache = list()
    
    # Ensure the network's CRS is correctly set for distance measurement in meters
    if (is.na(sf::st_crs(network)) || sf::st_crs(network)$units != "m") {
        network = sf::st_transform(network, crs = 27700)  # Example: UTM zone 32N for meters
    }

    # Generate a unique key for the cache based on the point's coordinates
    point_key = paste(sort(as.character(point)), collapse = "_")
    if (exists("path_cache") && point_key %in% names(path_cache)) {
        return(path_cache[[point_key]])
    }

    # Convert point and centroids to sfc if not already
    point_geom = sf::st_as_sfc(point)
    centroids_geom = sf::st_as_sfc(centroids)

    # Calculate distances between point and centroids
    distances = sf::st_distance(point_geom, centroids_geom)

    # Filter centroids based on the distance criteria
    valid_centroids = centroids[distances >= units::set_units(2, "m") & distances <= units::set_units(dist, "m"),]

    if (nrow(valid_centroids) > 0) {
        # Define weights based on path_type
        weights_to_use <- if (path_type == "all_simple") NA else "weight"
        
        # Calculate paths based on specified path_type
        paths_from_point = sfnetworks::st_network_paths(
            network, 
            from = point_geom, 
            to = sf::st_as_sfc(valid_centroids), 
            weights = weights_to_use,
            type = path_type
        )
        
        edges_in_paths = paths_from_point |>  
            dplyr::pull(edge_paths) |> 
            base::unlist() |> base::unique()

        result = network |> dplyr::slice(unique(edges_in_paths)) |> sf::st_as_sf()
    } else {
        result = NULL
    }

    path_cache[[point_key]] = result

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
  # Generate filenames for GeoJSON and pmtiles outputs
  coherent_geojson_filename = paste0(folder_path, city_filename, "_coherent_network",  ".geojson")
  # save coherent_geojson_filename as geojson
  sf::st_write(cohesive_network, coherent_geojson_filename, delete_dsn = TRUE)
  
  coherent_pmtiles_filename = paste0(folder_path, city_filename, "_coherent_network",  ".pmtiles")
  
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
