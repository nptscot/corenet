### Function: Compute Spatial Coverage
compute_spatial_coverage = function(rnet_core, lads, city_name = "City of Edinburgh", buffer_distance = 500) {

  city_boundary  = lads |> dplyr::filter(LAD23NM == city_name)
  rnet_core_zone = sf::st_intersection(rnet_core, city_boundary)
  # Create a buffer around the network
  network_buffer = sf::st_buffer(rnet_core_zone, dist = buffer_distance)
  # Union all buffer polygons
  buffer_union = sf::st_union(network_buffer)
  # Intersect with city boundary
  buffer_intersection = sf::st_intersection(buffer_union, city_boundary)
  
  # Compute areas
  area_buffered = sf::st_area(buffer_intersection)
  area_city = sf::st_area(city_boundary)
  
  # Spatial coverage ratio
  spatial_coverage = as.numeric(area_buffered / area_city)
  return(list(
    coverage = spatial_coverage,
    buffered_area = buffer_intersection
  ))
}

### Function: Compute Zone Connectivity
compute_zone_connectivity = function(intermediate_zone, lads, city_name = "City of Edinburgh", rnet_core, buffer_distance = 500, density_quantile = 0.3) {
  # Filter zones by density threshold
  city_boundary  = lads |> dplyr::filter(LAD23NM == city_name)
  intermediate_zone = sf::st_intersection(intermediate_zone, city_boundary)
  intermediate_zone$density = intermediate_zone$ResPop2011 / intermediate_zone$StdAreaKm2
  density_threshold = stats::quantile(intermediate_zone$density, density_quantile, na.rm = TRUE)
  intermediate_zone = intermediate_zone |> dplyr::filter(density > density_threshold)
  
  zones = intermediate_zone |>
    dplyr::select(InterZone, geometry) |>
    sf::st_make_valid()
  
  # Calculate buffer intersection for THIS specific network
  rnet_core_zone = sf::st_intersection(rnet_core, city_boundary)
  network_buffer = st_buffer(rnet_core_zone, dist = buffer_distance)
  buffer_union = st_union(network_buffer)
  W = sf::st_intersection(buffer_union, city_boundary)
  
  # Compute intersections W_i
  zones = zones |>
    dplyr::rowwise() |>
    dplyr::mutate(W_i = list(sf::st_intersection(geometry, W))) |>
    dplyr::ungroup()
  
  # Check intersections
  zones = zones |> mutate(has_intersection = lengths(W_i) > 0)
  
  # Build adjacency based on W_i connectivity
  zone_ids = zones$InterZone
  num_zones = length(zone_ids)
  adj_matrix = matrix(0, nrow = num_zones, ncol = num_zones,
                      dimnames = list(zone_ids, zone_ids))
  
  for (i in 1:num_zones) {
    for (j in i:num_zones) {
      if (i == j) {
        adj_matrix[i, j] = 1
      } else {
        geom_i = zones$W_i[[i]]
        geom_j = zones$W_i[[j]]
        
        intersects = sf::st_intersects(geom_i, geom_j, sparse = FALSE)
        touches = sf::st_touches(geom_i, geom_j, sparse = FALSE)
        
        if (any(intersects) | any(touches)) {
          adj_matrix[i, j] = 1
          adj_matrix[j, i] = 1
        }
      }
    }
  }
  
  all_connected = all(adj_matrix == 1)
  cat("Are all zones inter-connected within W? ", all_connected, "\n")

  g = igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  comp = igraph::components(g)
  largest_comp_size = max(comp$csize)
  total_zones = length(igraph::V(g))
  fraction_connected = largest_comp_size / total_zones
  return(list(
    graph = g,
    all_connected = (comp$no == 1),
    fraction_connected = fraction_connected,
    adj_matrix = adj_matrix,
    buffered_area = W  # Return the buffered area for potential use by other functions
  ))
}

### Function: Compute Cycling Potential Coverage
compute_cycling_potential_coverage = function(rnet_npt, lads, city_name = "City of Edinburgh", rnet_core, crs_target, buffer_distance = 20) {
  
  # Filter city network to within the city boundary
  city_boundary  = lads |> dplyr::filter(LAD23NM == city_name)
  rnet_core_zone = sf::st_intersection(rnet_core, city_boundary)
  rnet_city = sf::st_intersection(rnet_npt, city_boundary)
  
  # Compute length of each segment
  rnet_city = rnet_city |>
    dplyr::mutate(length_m = as.numeric(sf::st_length(geometry))) 
  
  # Total city potential (sum of all_fastest_bicycle_go_dutch)
  P_total = sum(rnet_city$all_fastest_bicycle_go_dutch, na.rm = TRUE)
  
  # Total city length
  L_city = sum(rnet_city$length_m, na.rm = TRUE)
  
  # City-wide density of potential
  D_city = P_total / L_city
  
  # Create a buffer around the core network
  rnet_core_buffer = sf::st_buffer(rnet_core_zone, buffer_distance)
  
  # Extract segments within the buffer
  rnet_city_buffer = rnet_city[sf::st_union(rnet_core_buffer), , op = sf::st_within]
  
  # Buffered potential sum
  P_U_city = sum(rnet_city_buffer$all_fastest_bicycle_go_dutch, na.rm = TRUE)
  
  # Buffered length
  L_buffer = sum(rnet_city_buffer$length_m, na.rm = TRUE)
  
  # Density of potential within the buffered area
  D_buffer = P_U_city / L_buffer
  
  # Coverage ratio
  coverage_ratio = D_buffer / D_city
  
  cycling_potential_coverage = P_U_city / P_total
  
  return(list(
    coverage = cycling_potential_coverage,
    D_city = D_city,
    D_buffer = D_buffer,
    coverage_ratio = coverage_ratio
  ))
}

### Function: Compute Population Coverage
compute_population_coverage = function(intermediate_zone, lads, city_name = "City of Edinburgh", rnet_core, dist_threshold = 500) {
  city_boundary  = lads |> dplyr::filter(LAD23NM == city_name)
  intermediate_zone = sf::st_intersection(intermediate_zone, city_boundary)
  rnet_core_zone = sf::st_intersection(rnet_core, city_boundary)

  zones = intermediate_zone |>
    dplyr::select(InterZone, TotPop2011, StdAreaKm2, geometry) |>
    sf::st_make_valid() |>
    dplyr::mutate(pop_density = TotPop2011 / StdAreaKm2)
  
  rnet_core_buffer = sf::st_buffer(rnet_core_zone, dist_threshold)
  W = sf::st_intersection(sf::st_union(rnet_core_buffer), city_boundary)
  
  zones_coverage = sf::st_intersection(zones, W)
  zones_coverage$covered_area = sf::st_area(zones_coverage)
  zones_coverage$covered_area_km2 = units::set_units(zones_coverage$covered_area, km^2)
  
  zones_coverage = zones_coverage |>
    dplyr::mutate(covered_population = pop_density * as.numeric(covered_area_km2))
  
  P_covered = sum(zones_coverage$covered_population, na.rm = TRUE)
  P_total = sum(zones$TotPop2011, na.rm = TRUE)
  
  population_coverage = P_covered / P_total
  return(population_coverage)
}

### Function: Compute O-D Accessibility
compute_od_accessibility = function(od_data, rnet_core, lads, city_name = "City of Edinburgh", dist_threshold = 500) {
  city_boundary  = lads |> filter(LAD23NM == city_name)
  rnet_core_zone = sf::st_intersection(rnet_core, city_boundary)
  od_data = sf::st_intersection(od_data, city_boundary) |> sf::st_cast("LINESTRING")
  rnet_core_union = st_union(rnet_core_zone)
  od_midpoints = st_line_sample(od_data, sample = 0.5)
  od_data_mid = st_as_sf(data.frame(od_data), geometry = od_midpoints)
  
  od_data_mid = od_data_mid |>
    mutate(distance_to_network = as.numeric(st_distance(geometry, rnet_core_union)))
  
  avg_distance = mean(od_data_mid$distance_to_network, na.rm = TRUE)
  
  od_count_within = sum(od_data_mid$distance_to_network <= dist_threshold, na.rm = TRUE)
  od_count_total = nrow(od_data_mid)
  od_coverage_count = od_count_within / od_count_total
  
  total_trips = sum(od_data_mid$all, na.rm = TRUE)
  trips_within = sum(od_data_mid$all[od_data_mid$distance_to_network <= dist_threshold], na.rm = TRUE)
  od_coverage_demand = trips_within / total_trips
  
  return(list(
    avg_distance = avg_distance,
    od_coverage_count = od_coverage_count,
    od_coverage_demand = od_coverage_demand
  ))
}

### Function: Compute OD-Based Directness and Efficiency
compute_directness_efficiency = function(rnet_core, od_points, lads, city_name = "City of Edinburgh", use_all_od_points = TRUE) {
  city_boundary = lads |> filter(LAD23NM == city_name)
  rnet_core_zone = sf::st_intersection(rnet_core, city_boundary) |> st_cast("LINESTRING")
  od_points_zone = od_points[city_boundary, ]

  cat("Network edges:", nrow(rnet_core_zone), "Total OD points:", nrow(od_points_zone), "\n")
  
  # Check if we have data
  if (nrow(rnet_core_zone) == 0 || nrow(od_points_zone) == 0) {
    return(list(D = NA, E_glob = NA, E_loc = NA, total_od_count = nrow(od_points_zone), routable_pairs = 0, total_pairs = 0))
  }

  # Create network
  network_sfnet = as_sfnetwork(rnet_core_zone, directed = FALSE)
  
  # Calculate edge weights more robustly to avoid geometry issues
  edges_sf = network_sfnet |> activate("edges") |> st_as_sf()
  edge_lengths = as.numeric(st_length(edges_sf))
  
  network_sfnet = network_sfnet |>
    activate("edges") |>
    mutate(weight = edge_lengths)
  
  g = network_sfnet |> as.igraph()
  
  # Get network nodes
  nodes_sf = network_sfnet |>
    activate("nodes") |>
    st_as_sf()
  
  # Use ALL OD points for fair comparison across networks
  # This ensures same denominator for all network comparisons
  od_points_filtered = od_points_zone
  n_od = nrow(od_points_filtered)
  
  cat("Using ALL", n_od, "OD points for fair comparison across networks\n")
  
  # Find nearest network nodes for all OD points
  od_nearest_nodes = st_nearest_feature(od_points_filtered, nodes_sf)
  
  # Create OD pairs (all combinations of OD points)
  od_coords = st_coordinates(od_points_filtered)
  
  # Sample OD pairs to make computation manageable
  max_pairs = 10000  # Limit to prevent excessive computation
  if (n_od > sqrt(max_pairs)) {
    sample_size = min(n_od, as.integer(sqrt(max_pairs)))
    sampled_indices = sample(1:n_od, sample_size)
    od_coords_sample = od_coords[sampled_indices, ]
    od_nearest_nodes_sample = od_nearest_nodes[sampled_indices]
    od_weights_sample = od_points_filtered$all[sampled_indices]
  } else {
    od_coords_sample = od_coords
    od_nearest_nodes_sample = od_nearest_nodes
    od_weights_sample = od_points_filtered$all
  }
  
  n_sample = nrow(od_coords_sample)
  
  # Create all pairs from sampled OD points
  od_pairs = expand.grid(origin = 1:n_sample, destination = 1:n_sample)
  od_pairs = od_pairs[od_pairs$origin != od_pairs$destination, ]  # Remove same-point pairs
  
  total_pairs = nrow(od_pairs)
  cat("Processing", total_pairs, "OD pairs from", n_sample, "sampled OD points\n")
  
  # Calculate euclidean distances for OD pairs
  origin_coords = od_coords_sample[od_pairs$origin, ]
  dest_coords = od_coords_sample[od_pairs$destination, ]
  od_euclidean_dist = sqrt(rowSums((origin_coords - dest_coords)^2))
  
  # Calculate network distances for OD pairs
  origin_nodes = od_nearest_nodes_sample[od_pairs$origin]
  dest_nodes = od_nearest_nodes_sample[od_pairs$destination]
  
  od_network_dist = numeric(length(origin_nodes))
  
  for (i in seq_along(origin_nodes)) {
    if (origin_nodes[i] != dest_nodes[i]) {
      tryCatch({
        path_dist = distances(g, v = origin_nodes[i], to = dest_nodes[i], weights = E(g)$weight)
        od_network_dist[i] = path_dist[1, 1]
      }, error = function(e) {
        od_network_dist[i] = Inf
      })
    } else {
      od_network_dist[i] = 0
    }
  }
  
  # Count successful routing attempts
  routable_pairs = sum(is.finite(od_network_dist) & od_euclidean_dist > 0 & od_network_dist > 0)
  cat("Successfully routable pairs:", routable_pairs, "out of", total_pairs, 
      "(", round(100 * routable_pairs / total_pairs, 1), "%)\n")
  
  # Calculate metrics for ALL pairs (assign 0 efficiency to unsuccessful ones)
  # Get weights for all pairs (average of origin and destination weights)
  origin_weights = od_weights_sample[od_pairs$origin]
  dest_weights = od_weights_sample[od_pairs$destination]
  pair_weights = (origin_weights + dest_weights) / 2
  
  # Initialize efficiency vectors with zeros
  directness_values = numeric(total_pairs)
  global_efficiency_values = numeric(total_pairs)
  
  # Calculate values only for valid pairs
  valid_pairs = is.finite(od_network_dist) & od_euclidean_dist > 0 & od_network_dist > 0
  
  if (sum(valid_pairs) > 0) {
    # OD-based Directness
    directness_values[valid_pairs] = od_euclidean_dist[valid_pairs] / od_network_dist[valid_pairs]
    
    # OD-based Global Efficiency
    global_efficiency_values[valid_pairs] = (od_network_dist[valid_pairs]) / (od_euclidean_dist[valid_pairs])
  }
  
  # Calculate weighted averages (unsuccessful pairs contribute 0)
  D = weighted.mean(directness_values, pair_weights, na.rm = TRUE)
  E_glob = weighted.mean(global_efficiency_values, pair_weights, na.rm = TRUE)
  
  # OD-based Local Efficiency
  # For each OD point, calculate efficiency within its local neighborhood
  calc_local_eff_od = function(od_idx, od_coords_sample, od_nearest_nodes_sample, od_weights_sample, g, radius = 2000) {
    center_coord = od_coords_sample[od_idx, ]
    
    # Find OD points within radius
    distances_to_center = sqrt(rowSums((sweep(od_coords_sample, 2, center_coord))^2))
    local_indices = which(distances_to_center <= radius & distances_to_center > 0)
    
    if (length(local_indices) < 2) {
      return(0)  # Return 0 instead of NA for fair comparison
    }
    
    # Calculate efficiency within local area
    local_coords = od_coords_sample[local_indices, ]
    local_nodes = od_nearest_nodes_sample[local_indices]
    local_weights = od_weights_sample[local_indices]
    
    # Create pairs within local area
    n_local = length(local_indices)
    local_pairs = expand.grid(1:n_local, 1:n_local)
    local_pairs = local_pairs[local_pairs$Var1 != local_pairs$Var2, ]
    
    if (nrow(local_pairs) == 0) {
      return(0)  # Return 0 instead of NA
    }
    
    # Calculate distances for local pairs
    local_euc = sqrt(rowSums((local_coords[local_pairs$Var1, ] - local_coords[local_pairs$Var2, ])^2))
    
    local_net = numeric(nrow(local_pairs))
    for (i in 1:nrow(local_pairs)) {
      node1 = local_nodes[local_pairs$Var1[i]]
      node2 = local_nodes[local_pairs$Var2[i]]
      if (node1 != node2) {
        tryCatch({
          path_dist = distances(g, v = node1, to = node2, weights = E(g)$weight)
          local_net[i] = path_dist[1, 1]
        }, error = function(e) {
          local_net[i] = Inf
        })
      } else {
        local_net[i] = 0
      }
    }
    
    # Calculate efficiency for ALL local pairs (assign 0 to unsuccessful ones)
    local_efficiency_values = numeric(nrow(local_pairs))
    valid_local = is.finite(local_net) & local_euc > 0 & local_net > 0
    
    if (sum(valid_local) > 0) {
      local_efficiency_values[valid_local] = (1/local_net[valid_local]) / (1/local_euc[valid_local])
    }
    
    # Calculate weighted average (unsuccessful pairs contribute 0)
    local_weights_pairs = (local_weights[local_pairs$Var1] + local_weights[local_pairs$Var2]) / 2
    
    return(weighted.mean(local_efficiency_values, local_weights_pairs, na.rm = TRUE))
  }
  
  # Calculate local efficiency for each OD point
  E_loc_values = sapply(1:n_sample, function(i) {
    calc_local_eff_od(i, od_coords_sample, od_nearest_nodes_sample, od_weights_sample, g)
  })
  
  E_loc = mean(E_loc_values, na.rm = TRUE)
  
  return(list(
    D = D,
    E_glob = E_glob,
    E_loc = E_loc,
    total_od_count = nrow(od_points_zone),
    routable_pairs = routable_pairs,
    total_pairs = total_pairs,
    routing_success_rate = routable_pairs / total_pairs
  ))
}

generate_radar_chart = function(city_name, 
                                rnet_core, lads, intermediate_zone, 
                                rnet_npt, crs_target, od_data, od_points,
                                dist_threshold = 500, buffer_distance = 500, 
                                save_path = NULL) {
  
  if (!requireNamespace("fmsb", quietly = TRUE)) {
    stop("Package 'fmsb' is required for this function. Please install it with install.packages('fmsb')")
  }
  
  # 1. Spatial Coverage
  sp_cov_result = compute_spatial_coverage(rnet_core, lads, city_name = city_name, buffer_distance = buffer_distance)
  spatial_coverage = sp_cov_result$coverage * 100  # Convert to percentage
  print(paste("Spatial coverage: ", spatial_coverage, "%"))
  # 2. Zone Connectivity
  zone_conn_result = compute_zone_connectivity(intermediate_zone, lads, city_name = city_name, rnet_core, buffer_distance = buffer_distance)
  zone_connectivity = zone_conn_result$fraction_connected * 100
  print(paste("Zone connectivity: ", zone_connectivity, "%"))
  
  # 3. Cycling Potential Coverage
  cp_cov_result = compute_cycling_potential_coverage(rnet_npt, lads, city_name = city_name, rnet_core, crs_target)
  cycling_potential_coverage = cp_cov_result$coverage * 100
  citywide_density_ratio = cp_cov_result$coverage_ratio  # Might need different scaling if > 1
  print(paste("Cycling potential coverage: ", cycling_potential_coverage, "%"))
  
  # 4. Population Coverage
  pop_cov = compute_population_coverage(intermediate_zone, lads, city_name = city_name, rnet_core, dist_threshold = dist_threshold)
  population_coverage = pop_cov * 100
  print(paste("Population coverage: ", population_coverage, "%"))
  
  # 5. O-D Accessibility
  od_result = compute_od_accessibility(od_data, rnet_core, lads, city_name = city_name, dist_threshold = dist_threshold)
  od_coverage_count = od_result$od_coverage_count * 100
  od_coverage_demand = od_result$od_coverage_demand * 100
  print(paste("OD coverage (count-based): ", od_coverage_count, "%"))
  
  # 6. OD-Based Directness and Efficiency
  de_result = compute_directness_efficiency(rnet_core, od_points, lads, city_name = city_name)
  directness = de_result$D
  global_efficiency = de_result$E_glob
  local_efficiency = de_result$E_loc
  print(paste("OD-based Directness: ", directness))
  print(paste("OD-based Global Efficiency: ", global_efficiency))
  print(paste("OD-based Local Efficiency: ", local_efficiency))
  
  # Combine metrics
  df_metrics_numeric <- data.frame(
    SpatialCoverage          = spatial_coverage,
    ZoneConnectivity         = zone_connectivity,
    CyclingPotentialCoverage = cycling_potential_coverage,
    CitywideDensityRatio     = citywide_density_ratio,
    PopulationCoverage       = population_coverage,
    AvgODDistance            = od_result$avg_distance,
    ODCoverageCountBased     = od_coverage_count,
    ODCoverageDemandWeighted = od_coverage_demand,
    Directness               = directness,
    GlobalEfficiency         = global_efficiency,
    LocalEfficiency          = local_efficiency
  )
  
  # Scale metrics for radar chart
  df_metrics_scaled = data.frame(
    SpatialCoverage          = df_metrics_numeric$SpatialCoverage / 100,
    ZoneConnectivity         = df_metrics_numeric$ZoneConnectivity / 100,
    CyclingPotentialCoverage = df_metrics_numeric$CyclingPotentialCoverage / 100,
    CitywideDensityRatio     = df_metrics_numeric$CitywideDensityRatio / 5,
    PopulationCoverage       = df_metrics_numeric$PopulationCoverage / 100,
    AvgODDistance            = df_metrics_numeric$AvgODDistance / 1000,
    ODCoverageCountBased     = df_metrics_numeric$ODCoverageCountBased / 100,
    ODCoverageDemandWeighted = df_metrics_numeric$ODCoverageDemandWeighted / 100,
    Directness               = df_metrics_numeric$Directness * 10,
    GlobalEfficiency         = df_metrics_numeric$GlobalEfficiency * 10,
    LocalEfficiency          = df_metrics_numeric$LocalEfficiency * 10
  )
  
  # Prepare radar chart data
  max_vals = rep(1, ncol(df_metrics_scaled))
  min_vals = rep(0, ncol(df_metrics_scaled))
  
  df_radar = rbind(
    max_vals,
    min_vals,
    df_metrics_scaled
  )
  
  # Rename columns for radar chart
  colnames(df_radar) = c(
    "Spatial\nCoverage",
    "Zone\nConnectivity",
    "Cycling Potential\nCoverage",
    "Density\nRatio",
    "Pop.\nCoverage",
    "Avg.\nOD\nDistance",
    "OD\nCoverage\n(Count)",
    "OD\nCoverage\n(Demand)",
    "OD-based\nDirectness",
    "OD-based\nGlobal\nEfficiency",
    "OD-based\nLocal\nEfficiency"
  )
  
  # Generate radar chart
  if (!is.null(save_path)) {
    # Open a PNG device to save the plot
    grDevices::png(filename = save_path, width = 800, height = 800)
  }

  fmsb::radarchart(
    df_radar,
    axistype = 1,
    seg = 5,
    pcol  = grDevices::rgb(0.2, 0.5, 0.5, 0.9),
    pfcol = grDevices::rgb(0.2, 0.5, 0.5, 0.5),
    plwd  = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "grey",
    caxislabels = seq(0, max(max_vals), length.out = 6),
    title = paste(city_name),
    # Increase font sizes
    cex.axis = 2,    # Axis label size
    cex.lab = 2,     # Axis title size (if applicable)
    cex.main = 2,       # Main title size
    vlcex = 1.5
  )
  
  if (!is.null(save_path)) {
    grDevices::dev.off()  # Close the PNG device
  }
}

