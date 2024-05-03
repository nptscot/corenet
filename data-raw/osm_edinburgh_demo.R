## code to prepare `osm_edinburgh_demo` dataset goes here

# Load necessary libraries
library(tidyverse)
library(osmactive)
library(zonebuilder)
library(sf)

# Install required GitHub package if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("nptscot/osmactive", quietly = TRUE)

# Define the geographic zone for Edinburgh at different radii
edinburgh = zonebuilder::zb_zone("Edinburgh")
edinburgh_3km = edinburgh |>
  # Change number in next line to change zone size:
  filter(circle_id <= 2) |>
  sf::st_union() |>
  sf::st_transform(27700)
edinburgh_6km = edinburgh |>
  # Change number in next line to change zone size:
  filter(circle_id <= 3) |>
  sf::st_union() |>
  sf::st_transform(27700)

# Get Open Roads data
if (!file.exists("srn.gpkg")) {
  message("Missing the SRN dataset locally, download it from the DfT or releases")
  u = "https://api.os.uk/downloads/v1/products/OpenRoads/downloads?area=GB&format=GeoPackage&redirect"
  f = "oproad_gpkg_gb.zip"
  if (!file.exists(f)) download.file(u, f)
  unzip(f)
  list.files("Data")
  # See https://github.com/acteng/atip-data-prep/blob/d2a5d0058932a00e0130048cf407448f4b75a477/layers/srn.py#L6
  # for filter:             "SELECT name_1 as name, geometry FROM road_link WHERE trunk_road"
  # q = "SELECT name_1 as name, geometry FROM road_link WHERE trunk_road"
  # srn = sf::read_sf("Data/oproad_gb.gpkg", query = q) # failed:
  open_roads_national = sf::read_sf("Data/oproad_gb.gpkg", layer = "road_link")
  names(open_roads_national)
  #    [1] "id"                         "fictitious"
  #  [3] "road_classification"        "road_function"
  #  [5] "form_of_way"                "road_classification_number"
  #  [7] "name_1"                     "name_1_lang"
  #  [9] "name_2"                     "name_2_lang"
  # [11] "road_structure"             "length"
  # [13] "length_uom"                 "loop"
  # [15] "primary_route"              "trunk_road"
  # [17] "start_node"                 "end_node"
  # [19] "road_number_toid"           "road_name_toid"
  # [21] "geometry"
  table(open_roads_national$trunk_road)
  open_roads_cleaned = open_roads_national |>
    # filter(trunk_road) |>
    transmute(
      name = name_1,
      road_function = road_function,
      form_of_way = form_of_way,
      road_classification = road_classification
    )
  open_roads_edinburgh_3km = open_roads_national[edinburgh_3km, , op = sf::st_within] |>
    dplyr::select(road_function)
  # Check the size of the dataset:
  nrow(open_roads_edinburgh_3km) / nrow(open_roads_national)
  object.size(open_roads_edinburgh_3km) |>
    # Format in MB:
    format(units = "MB")
  sf::write_sf(open_roads_edinburgh_3km, "open_roads_edinburgh_3km.gpkg", delete_dsn = TRUE)
  
  open_roads_edinburgh_6km = open_roads_national[edinburgh_6km, , op = sf::st_within] |>
    dplyr::select(road_function)
  # Check the size of the dataset:
  nrow(open_roads_edinburgh_6km) / nrow(open_roads_national)
  object.size(open_roads_edinburgh_6km) |>
    # Format in MB:
    format(units = "MB")
  sf::write_sf(open_roads_edinburgh_6km, "open_roads_edinburgh_6km.gpkg", delete_dsn = TRUE)  
  
  # Release data:
  if (FALSE) {
    usethis::use_github_release()
    system("gh release upload v0.0.1 open_roads_edinburgh_3km.gpkg --clobber")
    system("gh release upload v0.0.1 open_roads_edinburgh_6km.gpkg --clobber")
  }
}

os_edinburgh_demo_3km = sf::read_sf("open_roads_edinburgh_3km.gpkg")|>
  sf::st_transform(27700)
os_edinburgh_demo_6km = sf::read_sf("open_roads_edinburgh_6km.gpkg")|>
  sf::st_transform(27700)

# Remove CRS info to avoid warning with non ascii characters:
sf::st_crs(os_edinburgh_demo_3km) = NA
usethis::use_data(os_edinburgh_demo_3km, overwrite = TRUE)

sf::st_crs(os_edinburgh_demo_6km) = NA
usethis::use_data(os_edinburgh_demo_6km, overwrite = TRUE)
# Let's see how big the file is:
fs::file_size("data/os_edinburgh_demo_3km.rda")
fs::file_size("data/os_edinburgh_demo_6km.rda")

# Get npt data
NPT = sf::st_read("https://github.com/nptscot/corenet/releases/download/testing_data/combined_network_tile.geojson") |>
  sf::st_transform(27700)
NPT_demo_3km = NPT[edinburgh_3km, , op = sf::st_within]|>
    dplyr::select(all_fastest_bicycle_go_dutch)
NPT_demo_6km = NPT[edinburgh_6km, , op = sf::st_within]|>
    dplyr::select(all_fastest_bicycle_go_dutch)

sf::write_sf(NPT_demo_3km, "NPT_demo_3km.gpkg", delete_dsn = TRUE)  
sf::write_sf(NPT_demo_6km, "NPT_demo_6km.gpkg", delete_dsn = TRUE)  

sf::st_crs(NPT_demo_3km) = NA
usethis::use_data(NPT_demo_3km, overwrite = TRUE)

sf::st_crs(NPT_demo_6km) = NA
usethis::use_data(NPT_demo_6km, overwrite = TRUE)

# Let's see how big the file is:
fs::file_size("data/NPT_demo_3km.rda")
fs::file_size("data/NPT_demo_6km.rda")

# Get OSM data (for future reference)
osm = get_travel_network("Scotland", boundary = edinburgh_3km, boundary_type = "clipsrc")
names(osm)
#  [1] "osm_id"                "name"                  "highway"              
#  [4] "man_made"              "maxspeed"              "oneway"               
#  [7] "bicycle"               "cycleway"              "cycleway_left"        
# [10] "cycleway_right"        "cycleway_both"         "lanes"                
# [13] "lanes_both_ways"       "lanes_forward"         "lanes_backward"       
# [16] "lanes_bus"             "lanes_bus_conditional" "width"                
# [19] "segregated"            "sidewalk"              "footway"              
# [22] "service"               "surface"               "tracktype"            
# [25] "smoothness"            "access"                "z_order"              
# [28] "other_tags"            "geometry" 
cycle_net = get_cycling_network(osm)
drive_net = get_driving_network_major(osm)
cycle_net = distance_to_road(cycle_net, drive_net)
cycle_net = classify_cycle_infrastructure(cycle_net) 
names(cycle_net)
cycle_net |>
  select(cycle_segregation) |>
  plot()

usethis::use_data(osm_edinburgh_demo, overwrite = TRUE)
