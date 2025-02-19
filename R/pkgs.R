get_pkgs = function() {
   c(
  "crew", # For targets with workers
  "collapse", # Needed for bind_sf
  "cyclestreets", # For routing
  "geojsonsf", # For converting geojson to sf
  "geos", # For geometric operations
  "gert", # For interactive with git
  "glue", # For string interpolation
  "lubridate", # For working with dates and times
  "lwgeom", # For working with spatial data
  "nngeo", # Nearest neighbour functions
  "osmextract", # For extracting OpenStreetMap data
  "pct", # PCT interface
  "remotes", # For installing packages from remote sources
  "sf", # For working with spatial data
  "simodels", # For spatial interaction models
  "snakecase", # For converting strings to snake case
  "stplanr", # For sustainable transport planning
  "targets", # For managing targets in a workflow
  "tidyverse", # Includes dplyr, ggplot2, tidyr, stringr etc.
  "zonebuilder", # For creating zones for spatial analysis
  "iterators", # For creating iterators
  "doParallel", # For parallel processing
  "httr"
  )
}
