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

corenet = function(net, tbc , zone , crs = "EPSG:27700", attr1, attr2, attr3, attr1_val, attr2_val, attr3_val) {

  prepare_network = prepare_network(net,  transform_crs = crs, attr1 = attr1_val, attr2 = attr2_val, attr3 = attr3_val)

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


pkgs = c(
  "sf", 
  "dplyr", 
  "tidygraph", 
  "igraph", 
  "mapview", 
  "tidyverse", 
  "osmextract", 
  "stplanr", 
  "geos", 
  "dbscan", 
  "sfnetworks", 
  "rsgeo",
  "targets"
)
pkgs_not_installed = pkgs[!pkgs %in% installed.packages()]
if (length(pkgs_not_installed) > 0) {
  install.packages(pkgs_not_installed)
}
# Load the packages with vapply:
vapply(pkgs, require, character.only = TRUE, logical(1))

# Load all functions from the package
devtools::load_all()

remotes::install_dev("rsgeo", force = TRUE)
install.packages("zonebuilder")
