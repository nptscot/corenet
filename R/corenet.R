#' Generate a core network from road network data
#'
#' 
#' This function returns OSM keys that are relevant for active travel
#'

std_pkgs = c(
              "sf", 
              "dplyr", 
              "tidygraph", 
              "igraph", 
              "mapview", 
              "tidyverse", 
              "osmextract", 
              "geos", 
              "dbscan", 
              "sfnetworks", 
              "rsgeo",
              "targets",
              "zonebuilder"
            )

dev_pkgs = c(
              "stplanr", 
              "rsgeo"
            )

# Declare global variables to avoid 'no visible binding' notes
utils::globalVariables(c("edge_paths", "combined_network_tile", "all_fastest_bicycle_go_dutch", 
                         "weight", "to_linegraph", "edges", "group", "mean_potential", "LAD23NM", 
                         "road_function",  "grid_id", "density", 
                         "max_value", "min_value", "arterialness", "road_score"))

#' This function prepare the base netwrok for generating cohesive cycling network using NPT data.
#'
#' @param combined_network_tile NPT road network object, class sf.
#' @param crs Coordinate reference system to use, default is "EPSG:27700".
#' @param parameters A list containing specific parameters like 'coherent_area'.
#' @return A list of two elements: a cohesive network and zone data, both as sf objects.
#' @export


cohesive_network_prep = function(combined_network_tile, crs = "EPSG:27700", parameters) {
  # Initialize empty lists to store results and zones for each area
  cohesive_network_list = list()
  cohesive_zone_list = list()

  # Check for specified areas in the 'coherent_area' parameter
  if (!is.null(parameters$coherent_area) && length(parameters$coherent_area) > 0) {

    # Process each specified area
    for (area in parameters$coherent_area) {
      area_filename = gsub(" ", "_", area)
      print(paste("Preparing the network data for the", area))
      combined_network_tile = sf::st_transform(combined_network_tile, crs)

      # Read spatial data and filter for the specific area with transformation and buffering
      sg_intermediatezone_bdry_2011 = sf::st_read("inputdata/sg_intermediatezone_bdry_2011.gpkg") |> sf::st_transform(crs = crs)
      las_scotland_2023 = sf::st_read("inputdata/las_scotland_2023.geojson") |>
                           sf::st_transform(crs = crs) |>
                           dplyr::filter(LAD23NM == area) |>
                           sf::st_buffer(dist = 4000)

      # Calculate zones within the specified buffer area
      zones = sg_intermediatezone_bdry_2011[sf::st_union(las_scotland_2023), , op = sf::st_within]
      zones$density = zones$TotPop2011 / zones$StdAreaHa  

      # Intersect and transform NPT data with the calculated zones
      NPT_zones = combined_network_tile[sf::st_union(zones), , op = sf::st_intersects] |> 
                   sf::st_transform(crs)

      # Read and process road network data
      open_roads_national = sf::read_sf("inputdata/oproad_gb.gpkg", layer = "road_link") |> 
                   sf::st_transform(crs)

      OS_zones = open_roads_national[sf::st_union(zones), , op = sf::st_intersects]

      filtered_OS_zones = OS_zones |> 
                            dplyr::filter(road_function == 'A Road' | road_function == 'B Road' | road_function == 'Minor Road') |>
                            sf::st_transform(crs) |>
                            sf::st_zm()

      # Assign functions for data aggregation based on network attributes
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
      # OS_NPT_zones = stplanr::rnet_merge(OS_zones, NPT_zones, dist = 15, funs = funs, segment_length = 20,max_angle_diff = 10)
      filtered_OS_NPT_zones = stplanr::rnet_merge(filtered_OS_zones, NPT_zones, dist = 15, funs = funs, segment_length = 20,max_angle_diff = 10)
      # Store the processed data for the current area
      cohesive_network_list[[area]] = filtered_OS_NPT_zones
      cohesive_zone_list[[area]] = zones

      print(paste("Finished preparing the network data for the", area))
    }

    return(list(cohesive_zone = cohesive_zone_list, cohesive_network = cohesive_network_list))

  } else {
    print("No coherent area specified, proceeding with default settings.")
  }
}


#' Generate a base network for developing a cohesive cycling network
#'
#' This function processes the provided network data to extract and analyze
#' the most significant routes based on a percentile threshold, performs spatial operations,
#' clusters data points, and calculates the largest network component without dangles.
#' @param combined_network_tile The NPT network, class sf.
#' @param network_tile Spatial object ob.
#' @param combined_grid_buffer Additional spatial data, currently unused in the road_function
#' @param crs Coordinate reference system for transformation, default is "EPSG:27700".
#' @param dist The distance threshold used in path calculations, default is 10.
#' @return A spatial object representing the largest cohesive component of the network without dangles.
#' @export


corenet = function(combined_network_tile, network_tile, combined_grid_buffer, crs = "EPSG:27700", dist = 10) {

  if ("all_fastest_bicycle_go_dutch" %in% names(network_tile)) {

    combined_network_tile = sf::st_transform(combined_network_tile, crs)
    combined_network_tile_split = stplanr::line_segment(combined_network_tile, segment_length = 20, use_rsgeo = TRUE)

    filtered_roads = combined_network_tile_split[combined_network_tile_split$all_fastest_bicycle_go_dutch > 1500, ]

    # Calculate centroids of network segments
    centroids = sf::st_centroid(filtered_roads)

    # Perform DBSCAN clustering
    coordinates = sf::st_coordinates(centroids)
    clusters = dbscan::dbscan(coordinates, eps = 18, minPts = 1)
    centroids$cluster = clusters$cluster
    unique_centroids = centroids[!duplicated(centroids$cluster), ]    

  } else {
    grid_sf = combined_grid_buffer

    grid_sf$grid_id = seq_len(nrow(grid_sf))

    split_network = sf::st_intersection(combined_network_tile, grid_sf)

    select_by_density = function(density) {
      if (density < 10) {
        top_n = 1
      } else if (density >= 10 & density < 20) {
        top_n = 5
      } else {
        top_n = (density %/% 10) + 12
      }
      
      # Check if the value is NA (rest of the cases) and set top_n to 1
      if (is.na(density)) {
        top_n = 1
      }
      return(top_n)
    }

    distance_threshold = units::set_units(dist, "m") # Adjust the distance as needed

    split_network = split_network |>
      dplyr::group_by(grid_id) |>
      dplyr::mutate(top_n = select_by_density(unique(density)))


    split_network_max = split_network |>
      dplyr::group_by(grid_id) |>
      dplyr::group_modify(~ {
        # Step 1: Filter for the current grid_id
        split_network_filtered = .x
        
        # Step 2: Calculate Centroids for These Elements
        centroids = sf::st_centroid(split_network_filtered)
        
        # Step 3: Sort the Centroids Based on a Value Attribute
        sorted_centroids = centroids |>
          dplyr::arrange(desc(value))
        
        # Step 4: Identify the Centroid with the Highest Value
        highest_value_centroid = sorted_centroids[1, ]
        
        # Step 5: Calculate Distances from Each Centroid to the Highest Value Centroid
        distances_to_highest = sf::st_distance(centroids, highest_value_centroid)

        # Step 6: Select Line Strings Based on top_n and Distance Criteria
        top_n = unique(split_network_filtered$top_n)

        selected_lines = split_network_filtered |>
          dplyr::mutate(distance_to_highest = sf::st_distance(centroids, highest_value_centroid)[,1]) |>
          dplyr::filter(distance_to_highest >= distance_threshold) |>
          dplyr::slice_head(n = top_n)        
        return(selected_lines)
      }) |> 
      dplyr::bind_rows()

    # Calculate centroids of these max segments
    split_network_max = sf::st_as_sf(split_network_max)
    unique_centroids = sf::st_centroid(split_network_max)
  }

  # Prepare network and calculate paths
  prepared_network = prepare_network(network_tile, transform_crs = crs,  A_Road = 1, B_Road = 1, Minor_Road = 100000000)


  all_paths = purrr::map_dfr(
      seq_len(nrow(unique_centroids)),
      function(n) {
          calculate_paths_from_point_dist(prepared_network, point = unique_centroids[n,], centroids = unique_centroids, shortest = FALSE, dist = 1500)
      }
  )

  # Remove duplicates and calculate the largest connected component
  all_paths_sf_percentile_duplicated_geoms = duplicated(all_paths$geometry)
  all_paths_sf_percentile = all_paths[!all_paths_sf_percentile_duplicated_geoms, ]
  largest_component_sf = calculate_largest_component(all_paths_sf_percentile) |> sf::st_as_sf('edges')
  largest_component_sf_without_dangles = removeDangles(largest_component_sf, tolerance = 0.001)

  # Remove dangles multiple times to ensure a clean network
  for (i in 1:6) {
    largest_component_sf_without_dangles = removeDangles(largest_component_sf_without_dangles)
  }

  return(largest_component_sf_without_dangles)
}


#' This function generates a coherent network grouped by edge betweenness.
#'
#' @param CITY The road network for the city, expected to be an sf object.
#' @param ZONE Combined grid buffer, presumably an area of interest within the city, also an sf object.
#' @return A grouped sf network with ranked groups based on mean potential.
#' @export


coherent_network_group = function(CITY, ZONE) {
  # library(tidygraph)
  # Generate coherent network
  rnet_coherent = corenet(combined_network_tile, network_tile = CITY, combined_grid_buffer = ZONE)
  
  # Select relevant columns
  rnet_coherent_selected = rnet_coherent |>
    dplyr::select(all_fastest_bicycle_go_dutch, weight)
  
  # Group and process the network
  grouped_net = rnet_coherent_selected |>
    sfnetworks::as_sfnetwork(directed = FALSE) |>
    tidygraph::morph(tidygraph::to_linegraph) |>
    dplyr::mutate(group = tidygraph::group_edge_betweenness(n_groups = 12)) |>
    tidygraph::unmorph() |>
    tidygraph::activate(edges) |>
    sf::st_as_sf() |>
    sf::st_transform("EPSG:4326") |>
    dplyr::group_by(group) |>
    dplyr::summarise(mean_potential = mean(weight, na.rm = TRUE)) |>
    dplyr::mutate(group = rank(-mean_potential))
  
  # Return the processed network
  return(grouped_net)
}



#' Prepare a network data structure by transforming, scoring, and weighting based on road types and conditions
#'
#' This function takes a spatial network object, casts it to LINESTRING, converts it to an sfnetwork,
#' and computes scores based on road conditions and classifications. The transformation ensures that
#' the network is ready for further analytical processes.
#'
#' @param network An sf object representing a spatial network.
#' @param A_Road A numeric score for A Roads, default is 1.
#' @param B_Road A numeric score for B Roads, default is 1.
#' @param Minor_Road A numeric score for Minor Roads, default is 10000000.
#' Example usage
#' road_scores = list("A Road" = 1, "B Road" = 1, "Minor Road" = 100000)
#' @param transform_crs Numeric CRS code for coordinate transformation, default is 27700.
#' @return An sfnetwork object with additional attributes like 'arterialness' and 'weight' based on road conditions.
prepare_network = function(network, A_Road = 1, B_Road = 1, Minor_Road = 10000000, transform_crs = 27700) {
    # Cast the network to LINESTRING and transform to the specified CRS
    network = network |>
        sf::st_cast("LINESTRING") |>
        sf::st_transform(transform_crs) |>
        sfnetworks::as_sfnetwork(directed = FALSE) |>
        sfnetworks::activate("edges") |>
        dplyr::mutate(
            # Normalize 'all_fastest_bicycle_go_dutch' and handle NA values
            all_fastest_bicycle_go_dutch = dplyr::if_else(is.na(all_fastest_bicycle_go_dutch), 0, all_fastest_bicycle_go_dutch),
            min_value = min(all_fastest_bicycle_go_dutch, na.rm = TRUE),
            max_value = max(all_fastest_bicycle_go_dutch, na.rm = TRUE),
            arterialness = dplyr::if_else(
                max_value > min_value,
                (all_fastest_bicycle_go_dutch - min_value) / (max_value - min_value),
                0  # Avoid division by zero by setting arterialness to 0 when min_value equals max_value
            ),
            # Dynamically apply road type scores
            road_score = dplyr::case_when(
                road_function == "A Road" ~ A_Road,
                road_function == "B Road" ~ B_Road,
                road_function == "Minor Road" ~ Minor_Road,
                TRUE ~ 0  # Handle possible NA or unexpected values
            ),
            # Calculate weight considering the road type influence
            weight = (1 - arterialness) * 100 * (1 + 0.1 * road_score)
        )

    return(network)
}


#' Calculate paths from a given point to centroids within a specified distance range
#'
#' This function determines the network paths from a specific point to multiple centroids based on distance thresholds,
#' optionally returning the shortest paths.
#'
#' @param network An sfnetwork object representing the network.
#' @param point An sf point object from which paths will be calculated.
#' @param dist The maximum distance (in meters) to consider for path calculations.
#' @param centroids An sf object containing centroids to which paths are calculated.
#' @param shortest Logical indicating whether the shortest paths are calculated (TRUE) or weighted paths (FALSE).
#' @return An sf object containing the paths that meet the criteria or NULL if no paths meet the criteria.


calculate_paths_from_point_dist = function(network, point, dist = 500, centroids, shortest = FALSE) {
    
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
      if (shortest) {
          paths_from_point = sfnetworks::st_network_paths(network, from = point_geom, to = sf::st_as_sfc(valid_centroids), weights = NULL, type = "shortest") 
      } else {
          paths_from_point = sfnetworks::st_network_paths(network, from = point_geom, to = sf::st_as_sfc(valid_centroids), weights = "weight",type = "shortest") 
      }
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
#' @param network_tile An `sf` object representing the network to be analyzed.
#' @return An `sfnetwork` object representing the largest connected component of the network.
calculate_largest_component = function(network_tile) {
    # Ensure the network tile is in the correct format and convert to sfnetwork
    if (!inherits(network_tile, "sfnetwork")) {
        network_tile_sfn = sfnetworks::as_sfnetwork(network_tile, directed = FALSE)
    } else {
        network_tile_sfn = network_tile
    }

    # Convert sfnetwork to tbl_graph for easier graph manipulations
    network_tile_igraph = tidygraph::as_tbl_graph(network_tile_sfn)

    # Compute connected components using igraph
    components = igraph::components(network_tile_igraph)

    # Find the largest component by size
    largest_component_id = which.max(components$csize)
    largest_component_nodes = which(components$membership == largest_component_id)

    # Extract the subgraph corresponding to the largest component
    largest_component_subgraph = igraph::subgraph(network_tile_igraph, largest_component_nodes)

    # Convert back to sfnetwork
    largest_component_sfn = sfnetworks::as_sfnetwork(largest_component_subgraph, directed = FALSE)

    return(largest_component_sfn)
}


#' Remove dangling segments from a road network
#'
#' This function identifies and removes dangling line segments (lines with one end not connected to another line)
#' within a specified tolerance.
#'
#' @param road_network An `sf` object of class LINESTRING representing the road network.
#' @param tolerance The distance tolerance for identifying isolated endpoints as dangling.
#' @return An `sf` object with dangling line segments removed.
removeDangles = function(road_network, tolerance = 0.001) {
    # Convert to Spatial Lines if not already
    road_network_lines = sf::st_cast(road_network, "LINESTRING")

    # Extract and combine all end points of line segments
    end_points = do.call(rbind, lapply(road_network_lines$geometry, function(line) {
    endpoints = rbind(sf::st_coordinates(line)[1, ], sf::st_coordinates(line)[nrow(sf::st_coordinates(line)), ])
    sf::st_as_sf(data.frame(x = endpoints[, 1], y = endpoints[, 2]), coords = c("x", "y"), crs = sf::st_crs(road_network))
    }))

    # Identify unique end points (potential dangles) using a spatial join to find nearby points
    buffer_points = sf::st_buffer(end_points, dist = tolerance)
    overlaps = sf::st_intersects(buffer_points, buffer_points, sparse = FALSE)
    isolated_points = end_points[rowSums(overlaps) == 1,]

    # Filter out road segments that end in these isolated points
    segments_with_dangles = sapply(sf::st_geometry(road_network_lines), function(geom) {
    ends = sf::st_sfc(sf::st_point(sf::st_coordinates(geom)[1,]), sf::st_point(sf::st_coordinates(geom)[nrow(sf::st_coordinates(geom)),]), crs = sf::st_crs(road_network))
    any(sf::st_intersects(ends, isolated_points, sparse = FALSE))
    })

    road_network_without_dangles = road_network_lines[!segments_with_dangles, ]

    return(road_network_without_dangles)
}



