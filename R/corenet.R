#' Generate a core network from road network data
#'
#' 
#' This function returns OSM keys that are relevant for active travel
#'
#' @export

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

#' This function prepare the base netwrok for generating cohesive cycling network using NPT data.
#'
#' @param combined_network_tile NPT road network object, class sf.
#' @param crs Coordinate reference system to use, default is "EPSG:27700".
#' @param parameters A list containing specific parameters like 'coherent_area'.
#' @return A list of two elements: a cohesive network and zone data, both as sf objects.
#' @export
#' @examples
#' 

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
      sg_intermediatezone_bdry_2011 = sf::st_read("inputdata/sg_intermediatezone_bdry_2011.gpkg")
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
      open_roads_national = sf::read_sf("Data/oproad_gb.gpkg", layer = "road_link")
      OS_zones = open_roads_national[sf::st_union(zones), , op = sf::st_intersects]

      # Filter and transform OS zone data
      filtered_OS_zones = OS_zones |> 
                           dplyr::filter(function. == 'A Road' | function. == 'B Road' | function. == 'Minor Road') |>
                           sf::st_transform(crs) |>
                           sf::st_zm()

      # Assign functions for data aggregation based on network attributes
      name_list = names(NPT_zones)
      funs = list()
      for (name in name_list) {
        if (name == "geometry") {
          next  # Skip geometry field
        } elseif (name %in% c("gradient", "quietness")) {
          funs[[name]] = mean
        } else {
          funs[[name]] = sum
        }
      }

      # Merge road networks with specified parameters
      OS_NPT_zones = stplanr::rnet_merge(filtered_OS_zones, NPT_zones, dist = 15, segment_length = 10, funs = funs, max_angle_diff = 10)

      # Store the processed data for the current area
      cohesive_network_list[[area]] = OS_NPT_zones
      cohesive_zone_list[[area]] = zones

      print(paste("Finished preparing the network data for the", area))
    }

    return(list(cohesive_network = cohesive_network_list, cohesive_zone = cohesive_zone_list))

  } else {
    print("No coherent area specified, proceeding with default settings.")
  }
}


#' Generate a cohesive network based on the prepared network data
#'
#' This function processes the provided network data to extract and analyze
#' the most significant routes based on a percentile threshold, performs spatial operations,
#' clusters data points, and calculates the largest network component without dangles.
#'
#' @param network_tile Spatial object representing network data.
#' @param combined_grid_buffer Additional spatial data, currently unused in the function.
#' @param base_value Attribute column name to analyze in `network_tile`. Default is 'all_fastest_bicycle_go_dutch'.
#' @param crs Coordinate reference system for transformation, default is "EPSG:27700".
#' @param min_percentile Minimum percentile for filtering network data, default is 0.85.
#' @param dist The distance threshold used in path calculations, default is 10.
#' @return A spatial object representing the largest cohesive component of the network without dangles.
#' @export
#' @examples

corenet = function(network_tile, combined_grid_buffer, base_value = NULL, crs = "EPSG:27700", min_percentile = 0.85, dist = 10) {
  if (!is.null(base_value) && base_value %in% names(network_tile)) {
    # Calculate the minimum value at the specified percentile
    min_percentile_value = stats::quantile(network_tile[[base_value]], probs = min_percentile)
    network_tile_filtered = dplyr::filter(network_tile, !!rlang::sym(base_value) > min_percentile_value)

    # Split network data into manageable segments
    segment_lengths = sf::st_length(network_tile_filtered)
    n_segments = n_segments(segment_lengths, 20)
    network_tile_split = line_segment_rsgeo(network_tile_filtered, n_segments = n_segments)

    # Calculate centroids of network segments
    centroids = sf::st_centroid(network_tile_split)

    # Perform DBSCAN clustering
    coordinates = sf::st_coordinates(centroids)
    clusters = dbscan::dbscan(coordinates, eps = 18, minPts = 1)
    centroids$cluster = clusters$cluster
    unique_centroids = centroids[!duplicated(centroids$cluster), ]    

  } else {
    # New method for generating centroids
    grid_sf = combined_grid_buffer
    grid_sf$grid_id = seq_len(nrow(grid_sf))

    split_network = sf::st_intersection(network_tile, grid_sf)

    split_network = split_network |>
      group_by(grid_id) |>
      mutate(top_n = select_by_density(unique(density))) |>
      group_modify(~ {
        split_network_filtered = .x
        centroids = sf::st_centroid(split_network_filtered)
        sorted_centroids = centroids |>
          arrange(desc(value))
        highest_value_centroid = sorted_centroids[1, ]
        distances_to_highest = sf::st_distance(centroids, highest_value_centroid)
        selected_lines = split_network_filtered |>
          mutate(distance_to_highest = distances_to_highest[,1]) |>
          filter(distance_to_highest >= units::set_units(dist, "m")) |>
          slice_head(n = unique(.x$top_n))
        return(selected_lines)
      }) |>
      bind_rows()

    # Calculate centroids of these max segments
    unique_centroids = sf::st_centroid(sf::st_as_sf(split_network))
  }

  # Prepare network and calculate paths
  prepare_network = prepare_network(net, transform_crs = crs, attr1 = attr1_val, attr2 = attr2_val, attr3 = attr3_val)

  all_paths = purrr::map_dfr(
    seq_len(nrow(unique_centroids)),
    function(n) {
      calculate_paths_from_point_dist(prepare_network, point = unique_centroids[n, ], centroids = unique_centroids, shortest = FALSE, dist = 1500)
    }
  )

  # Remove duplicates and calculate the largest connected component
  all_paths_sf_percentile_duplicated_geoms = duplicated(all_paths$geometry)
  all_paths_sf_percentile = all_paths[!all_paths_sf_percentile_duplicated_geoms, ]
  largest_component_sf = calculate_largest_component(all_paths_sf_percentile) |> sf::st_as_sf('edges')
  largest_component_sf_without_dangles = removeDangles(largest_component_sf, tolerance = 0.001)

  # Remove dangles multiple times to ensure a clean network
  for (i in 1:4) {
    largest_component_sf_without_dangles = removeDangles(largest_component_sf_without_dangles)
  }

  return(largest_component_sf_without_dangles)
}



#' Get the Edinburgh road network, within radius of 6 km of the center
#'
#' @export
#' @examples
#' get_edinburgh_6km()
get_edinburgh_6km = function() {
  sf::read_sf("https://github.com/nptscot/corenet/releases/download/v0.0.1/open_roads_edinburgh_6km.gpkg")
}


#' Data from OS Open Roads, Edinburgh, 3 km radius
#'
#'
#' @docType data
#' @keywords datasets
#' @name os_edinburgh_demo
#' @format An sf data frame
#' @examples
#' library(sf)
#' names(os_edinburgh_demo)
#' plot(os_edinburgh_demo)
#' # Larger example:
#' # os_edinburgh = get_edinburgh_6km()
NULL


#' Prepare a network data structure by transforming, scoring, and weighting based on road types and conditions
#'
#' This function takes a spatial network object, casts it to LINESTRING, converts it to an sfnetwork,
#' and computes scores based on road conditions and classifications. The transformation ensures that
#' the network is ready for further analytical processes.
#'
#' @param network An sf object representing a spatial network.
#' @param road_scores A named list where names are road types and values are scores associated with each road type.
#' Example usage
#' road_scores = list("A Road" = 1, "B Road" = 1, "Minor Road" = 100000)
#' @param transform_crs Numeric CRS code for coordinate transformation, default is 27700.
#' @return An sfnetwork object with additional attributes like 'arterialness' and 'weight' based on road conditions.
prepare_network = function(network, road_scores, transform_crs = 27700) {
    # Cast the network to LINESTRING and transform to the specified CRS
    network = network |>
        sf::st_cast("LINESTRING") |>
        sf::st_transform(transform_crs) |>
        sfnetworks::as_sfnetwork(directed = FALSE) |>
        sfnetworks::activate("edges") |>
        dplyr::mutate(
            # Normalize 'all_fastest_bicycle_go_dutch' and handle NA values
            all_fastest_bicycle_go_dutch = if_else(is.na(all_fastest_bicycle_go_dutch), 0, all_fastest_bicycle_go_dutch),
            min_value = min(all_fastest_bicycle_go_dutch, na.rm = TRUE),
            max_value = max(all_fastest_bicycle_go_dutch, na.rm = TRUE),
            arterialness = if_else(
                max_value > min_value,
                (all_fastest_bicycle_go_dutch - min_value) / (max_value - min_value),
                0  # Avoid division by zero by setting arterialness to 0 when min_value equals max_value
            ),
            # Dynamically apply road type scores
            road_score = case_when(
                !!rlang::sym(names(road_scores)[1]) ~ road_scores[[1]],
                !!rlang::sym(names(road_scores)[2]) ~ road_scores[[2]],
                !!rlang::sym(names(road_scores)[3]) ~ road_scores[[3]],
                TRUE ~ 0  # Default score for undefined road types or missing data
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
        # Determine path type based on shortest parameter
        weight_column = if (shortest) NULL else "weight"
        paths_from_point = sfnetworks::st_network_paths(
            network, 
            from = point_geom, 
            to = sf::st_as_sfc(valid_centroids), 
            weights = weight_column,
            type = "shortest"
        )

        # Extract unique edges from the paths
        edges_in_paths = paths_from_point %>%
            dplyr::pull(edge_paths) %>%
            base::unlist() %>%
            base::unique()

        result = network %>%
            dplyr::slice(unique(edges_in_paths)) %>%
            sf::st_as_sf()
    } else {
        result = NULL
    }

    # Update the cache with the result
    if (exists("path_cache")) {
        path_cache[[point_key]] = result
    }

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
    # Ensure the network is cast to LINESTRING if not already
    if (sf::st_geometry_type(road_network) != "LINESTRING") {
        road_network = sf::st_cast(road_network, "LINESTRING")
    }

    # Extract and combine all end points of line segments
    end_points = do.call(rbind, lapply(road_network_lines$geometry, function(line) {
    endpoints = rbind(st_coordinates(line)[1, ], st_coordinates(line)[nrow(st_coordinates(line)), ])
    st_as_sf(data.frame(x = endpoints[, 1], y = endpoints[, 2]), coords = c("x", "y"), crs = st_crs(road_network))
    }))

    # Identify unique end points (potential dangles) using a spatial join to find nearby points
    buffer_points = st_buffer(end_points, dist = tolerance)
    overlaps = st_intersects(buffer_points, buffer_points, sparse = FALSE)
    isolated_points = end_points[rowSums(overlaps) == 1,]

    # Filter out road segments that end in these isolated points
    segments_with_dangles = sapply(st_geometry(road_network_lines), function(geom) {
    ends = st_sfc(st_point(st_coordinates(geom)[1,]), st_point(st_coordinates(geom)[nrow(st_coordinates(geom)),]), crs = st_crs(road_network))
    any(st_intersects(ends, isolated_points, sparse = FALSE))
    })

    road_network_without_dangles = road_network_lines[!segments_with_dangles, ]

    return(road_network_without_dangles)
}


#' Calculate the number of segments for line division
#'
#' This function computes the number of segments a line should be divided into,
#' ensuring that each segment does not exceed a specified maximum length.
#'
#' @param line_length The total length of the line.
#' @param max_segment_length The maximum allowable length of each segment.
#' @return An integer indicating the number of segments.
n_segments = function(line_length, max_segment_length) {
  pmax(ceiling(line_length / max_segment_length), 1)
}


#' Segmentize line geometries in an sf object
#'
#' This function divides each line in the input sf object into a specified number of segments,
#' handling coordinate systems and preserving original attributes.
#'
#' @param l An sf object containing LINESTRING geometries.
#' @param n_segments The number of segments to divide each line into.
#' @return An sf object with the original attributes and new geometries divided into segments.
line_segment_rsgeo = function(l, n_segments) {

  crs = sf::st_crs(l)
  # Test to see if the CRS is latlon or not and provide warning if so
  if (sf::st_is_longlat(l)) {
    warning(
      "The CRS of the input object is latlon.\n",
      "This may cause problems with the rsgeo implementation of line_segment()."
    )
  }

  # extract geometry and convert to rsgeo
  geo = rsgeo::as_rsgeo(sf::st_geometry(l))

  # segmentize the line strings
  res_rsgeo = rsgeo::line_segmentize(geo, n_segments)

  # make them into sfc_LINESTRING
  res = sf::st_cast(sf::st_as_sfc(res_rsgeo), "LINESTRING")

  # give them them CRS
  res = sf::st_set_crs(res, crs)

  # calculate the number of original geometries
  n_lines = length(geo)
  # create index ids to grab rows from
  ids = rep.int(seq_len(n_lines), n_segments)

  # index the original sf object
  res_tbl = sf::st_drop_geometry(l)[ids, , drop = FALSE]

  # assign the geometry column
  nrow(res_tbl)

  res_tbl[[attr(l, "sf_column")]] = res

  # convert to sf and return
  res_sf = sf::st_as_sf(res_tbl)
  res_sf
}

