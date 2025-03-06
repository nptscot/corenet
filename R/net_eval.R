### Function: Compute Spatial Coverage
compute_spatial_coverage = function(rnet_core, lads, city_name = "City of Edinburgh", buffer_distance = 500) {

  city_boundary  = lads |> filter(LAD23NM == city_name)
  rnet_core_zone = sf::st_intersection(rnet_core, city_boundary)
  # Create a buffer around the network
  network_buffer = st_buffer(rnet_core_zone, dist = buffer_distance)
  # Union all buffer polygons
  buffer_union = st_union(network_buffer)
  # Intersect with city boundary
  buffer_intersection = sf::st_intersection(buffer_union, city_boundary)
  
  # Compute areas
  area_buffered = st_area(buffer_intersection)
  area_city = st_area(city_boundary)
  
  # Spatial coverage ratio
  spatial_coverage = as.numeric(area_buffered / area_city)
  return(list(
    coverage = spatial_coverage,
    buffered_area = buffer_intersection
  ))
}

### Function: Compute Zone Connectivity
compute_zone_connectivity = function(intermediate_zone, lads, city_name = "City of Edinburgh", buffer_intersection, density_quantile = 0.3) {
  # Filter zones by density threshold
  city_boundary  = lads |> filter(LAD23NM == city_name)
  intermediate_zone = sf::st_intersection(intermediate_zone, city_boundary)
  intermediate_zone$density = intermediate_zone$ResPop2011 / intermediate_zone$StdAreaKm2
  density_threshold = quantile(intermediate_zone$density, density_quantile, na.rm = TRUE)
  intermediate_zone = intermediate_zone |> filter(density > density_threshold)
  
  zones = intermediate_zone |>
    select(InterZone, geometry) |>
    st_make_valid()
  
  # W = B_union ∩ A_city (already computed outside)
  W = buffer_intersection
  
  # Compute intersections W_i
  zones = zones |>
    rowwise() |>
    mutate(W_i = list(st_intersection(geometry, W))) |>
    ungroup()
  
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
        
        intersects = st_intersects(geom_i, geom_j, sparse = FALSE)
        touches = st_touches(geom_i, geom_j, sparse = FALSE)
        
        if (any(intersects) | any(touches)) {
          adj_matrix[i, j] = 1
          adj_matrix[j, i] = 1
        }
      }
    }
  }
  
  all_connected = all(adj_matrix == 1)
  cat("Are all zones inter-connected within W? ", all_connected, "\n")

  g = graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  comp = components(g)
  largest_comp_size = max(comp$csize)
  total_zones = length(V(g))
  fraction_connected = largest_comp_size / total_zones
  return(list(
    graph = g,
    all_connected = (comp$no == 1),
    fraction_connected = fraction_connected,
    adj_matrix = adj_matrix
  ))
}

### Function: Compute Cycling Potential Coverage
compute_cycling_potential_coverage = function(rnet_npt, lads, city_name = "City of Edinburgh", rnet_core, crs_target, buffer_distance = 20) {
  
  # Filter city network to within the city boundary
  city_boundary  = lads |> filter(LAD23NM == city_name)
  rnet_core_zone = sf::st_intersection(rnet_core, city_boundary)
  rnet_city = sf::st_intersection(rnet_npt, city_boundary)
  
  # Compute length of each segment
  rnet_city = rnet_city |>
    mutate(length_m = as.numeric(st_length(geometry))) 
  
  # Total city potential (sum of all_fastest_bicycle_go_dutch)
  P_total = sum(rnet_city$all_fastest_bicycle_go_dutch, na.rm = TRUE)
  
  # Total city length
  L_city = sum(rnet_city$length_m, na.rm = TRUE)
  
  # City-wide density of potential
  D_city = P_total / L_city
  
  # Create a buffer around the core network
  rnet_core_buffer = st_buffer(rnet_core_zone, buffer_distance)
  
  # Extract segments within the buffer
  rnet_city_buffer = rnet_city[st_union(rnet_core_buffer), , op = st_within]
  
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
  city_boundary  = lads |> filter(LAD23NM == city_name)
  intermediate_zone = sf::st_intersection(intermediate_zone, city_boundary)
  rnet_core_zone = sf::st_intersection(rnet_core, city_boundary)

  zones = intermediate_zone |>
    select(InterZone, TotPop2011, StdAreaKm2, geometry) |>
    st_make_valid() |>
    mutate(pop_density = TotPop2011 / StdAreaKm2)
  
  rnet_core_buffer = st_buffer(rnet_core_zone, dist_threshold)
  W = st_intersection(st_union(rnet_core_buffer), city_boundary)
  
  zones_coverage = st_intersection(zones, W)
  zones_coverage$covered_area = st_area(zones_coverage)
  zones_coverage$covered_area_km2 = units::set_units(zones_coverage$covered_area, km^2)
  
  zones_coverage = zones_coverage |>
    mutate(covered_population = pop_density * as.numeric(covered_area_km2))
  
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

### Function: Compute Directness and Efficiency
compute_directness_efficiency = function(rnet_core, lads, city_name = "City of Edinburgh") {
  
  # 1) Choose a planar coordinate system for distances in meters.
  #    For Great Britain, EPSG:27700 is common.
  target_crs = 27700
  
  # 2) Reproject boundary and network lines to the chosen CRS
  lads_proj       = st_transform(lads, target_crs)
  rnet_core_proj  = st_transform(rnet_core, target_crs)
  
  # 3) Clip the network to the chosen city
  city_boundary   = lads_proj |> filter(LAD23NM == city_name)
  rnet_core_zone  = st_intersection(rnet_core_proj, city_boundary)
  rnet_core_zone  = st_cast(rnet_core_zone, "LINESTRING")  # ensure LINESTRING geometry
  
  # 4) Build an initial sfnetwork (UNDIRECTED)
  network_sfnet = as_sfnetwork(rnet_core_zone, directed = FALSE)
  
  # 5) Subdivide lines at all intersections
  #    This ensures real intersections become shared nodes
  network_sfnet_subdiv = network_sfnet %>%
    convert(to_spatial_subdivision) %>%
    activate("edges") %>%
    mutate(weight = st_length(geometry))  # edge weights in meters
  
  # Convert to igraph for distance computations
  g = as.igraph(network_sfnet_subdiv)
  
  # 6) Compute distance matrices:
  #    dist_matrix_g = shortest-path distances on the network
  #    dist_matrix_e = straight-line (Euclidean) distances
  dist_matrix_g = distances(g, weights = E(g)$weight)
  
  # Extract node coordinates in the same planar CRS
  nodes_sf = network_sfnet_subdiv %>%
    activate("nodes") %>%
    st_as_sf()
  
  coords = st_coordinates(nodes_sf)
  dist_matrix_e = as.matrix(dist(coords))
  
  # ------------------------------------------------------------------
  # 7) Fraction-Weighted Directness & Global Efficiency
  #
  #    1) Among connected node pairs only
  #    2) Multiply by the fraction of all node pairs that are connected
  # ------------------------------------------------------------------
  all_pairs_mask = upper.tri(dist_matrix_g, diag = FALSE)
  
  connected_mask = all_pairs_mask &
                   is.finite(dist_matrix_g) &
                   (dist_matrix_g > 0)
  
  # Fraction of connected pairs
  n_all_pairs     = sum(all_pairs_mask)
  n_conn_pairs    = sum(connected_mask)
  frac_connected  = n_conn_pairs / n_all_pairs
  
  # Directness among connected pairs
  directness_ratios_conn = dist_matrix_e[connected_mask] / dist_matrix_g[connected_mask]
  D_conn = mean(directness_ratios_conn, na.rm = TRUE)
  
  # Fraction-weighted Directness
  D_fraction_weighted = D_conn * frac_connected
  
  # Global Efficiency among connected pairs
  E_glob_conn = sum(1 / dist_matrix_g[connected_mask], na.rm = TRUE) /
                sum(1 / dist_matrix_e[connected_mask], na.rm = TRUE)
  
  # Fraction-weighted Global Efficiency
  E_glob_fraction_weighted = E_glob_conn * frac_connected
  
  # ------------------------------------------------------------------
  # 8) Local Efficiency
  #
  #    Standard approach: for each node's neighbors, ignore pairs that
  #    aren't connected. If < 2 neighbors, local efficiency = NA.
  # ------------------------------------------------------------------
  calc_local_eff = function(node_id, distG, distE, igraph_obj) {
    nbrs = as.numeric(neighbors(igraph_obj, node_id))
    if (length(nbrs) < 2) {
      return(NA)
    }
    distG_sub = distG[nbrs, nbrs]
    distE_sub = distE[nbrs, nbrs]
    
    valid_sub = upper.tri(distG_sub, diag = FALSE) &
                is.finite(distG_sub) &
                (distG_sub > 0)
    if (sum(valid_sub) == 0) return(NA)
    
    sum_invG = sum(1 / distG_sub[valid_sub], na.rm = TRUE)
    sum_invE = sum(1 / distE_sub[valid_sub], na.rm = TRUE)
    if (sum_invE == 0) return(NA)
    
    return(sum_invG / sum_invE)
  }
  
  E_loc_values = sapply(seq_len(vcount(g)), function(i) {
    calc_local_eff(i, dist_matrix_g, dist_matrix_e, g)
  })
  E_loc = mean(E_loc_values, na.rm = TRUE)
  
  # ------------------------------------------------------------------
  # 9) Return all metrics
  # ------------------------------------------------------------------
  return(list(
    # Fraction of connected pairs
    frac_connected            = frac_connected,
    
    # Within-connected subgraphs only
    D_conn                    = D_conn,         # directness among connected pairs
    E_glob_conn               = E_glob_conn,    # global efficiency among connected pairs
    
    # Fraction-weighted directness & efficiency
    D_fraction_weighted       = D_fraction_weighted,
    E_glob_fraction_weighted  = E_glob_fraction_weighted,
    
    # Local Efficiency
    E_loc                     = E_loc
  ))
}

generate_radar_chart = function(city_name, 
                                rnet_core, lads, intermediate_zone, 
                                rnet_npt, crs_target, od_data, 
                                dist_threshold = 500, buffer_distance = 500, 
                                save_path = NULL) {
  
  library(fmsb)
  
  # 1. Spatial Coverage
  sp_cov_result = compute_spatial_coverage(rnet_core, lads, city_name = city_name, buffer_distance = buffer_distance)
  spatial_coverage = sp_cov_result$coverage * 100  # Convert to percentage
  print(paste("Spatial coverage: ", spatial_coverage, "%"))
  # 2. Zone Connectivity
  zone_conn_result = compute_zone_connectivity(intermediate_zone, lads, city_name = city_name, sp_cov_result$buffered_area)
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
  
  # 6. Directness and Efficiency
  de_result = compute_directness_efficiency(rnet_core, lads, city_name = city_name)
  directness = de_result$D_fraction_weighted
  global_efficiency = de_result$E_glob_fraction_weighted
  local_efficiency = de_result$E_loc
  print(paste("Directness: ", directness))
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
    Directness               = df_metrics_numeric$Directness,
    GlobalEfficiency         = df_metrics_numeric$GlobalEfficiency,
    LocalEfficiency          = df_metrics_numeric$LocalEfficiency
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
    "Directness",
    "Global\nEfficiency",
    "Local\nEfficiency"
  )
  
  # Generate radar chart
  if (!is.null(save_path)) {
    # Open a PNG device to save the plot
    png(filename = save_path, width = 800, height = 800)
  }

  radarchart(
    df_radar,
    axistype = 1,
    seg = 5,
    pcol  = rgb(0.2, 0.5, 0.5, 0.9),
    pfcol = rgb(0.2, 0.5, 0.5, 0.5),
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
    dev.off()  # Close the PNG device
  }
}

