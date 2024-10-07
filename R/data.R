#' Edinburgh OpenRoads Network - 3km Demo Data
#'
#' This data set contains network data for Edinburgh within a 3km radius, sourced from OpenRoads.
#' It is formatted as a spatial data frame (sf object), suitable for spatial analysis demonstrations.
#' 
#' @docType data
#' @keywords datasets
#' @name os_edinburgh_demo_3km
#' @format An object of class \code{sf} (inherits from \code{data.frame}).
#' @examples
#' data(os_edinburgh_demo_3km)
#' head(os_edinburgh_demo_3km)
"os_edinburgh_demo_3km"

#' Edinburgh OpenRoads Network - 6km Demo Data
#'
#' This data set contains network data for Edinburgh within a 6km radius, sourced from OpenRoads.
#' It is formatted as a spatial data frame (sf object), which is commonly used in spatial data analysis.
#'
#' @docType data
#' @keywords datasets
#' @name os_edinburgh_demo_6km
#' @format An object of class \code{sf} (inherits from \code{data.frame}).
#' @examples
#' data(os_edinburgh_demo_6km)
#' head(os_edinburgh_demo_6km)
"os_edinburgh_demo_6km"

#' Edinburgh NPT Network - 3km Demo Data
#'
#' This data set contains network data for Edinburgh within a 3km radius, sourced from NPT.
#' It highlights specific network characteristics important for transportation planning.
#'
#' @docType data
#' @keywords datasets
#' @name NPT_demo_3km
#' @format An object of class \code{sf} (inherits from \code{data.frame}).
#' @examples
#' data(NPT_demo_3km)
#' head(NPT_demo_3km)
"NPT_demo_3km"

#' Edinburgh NPT Network - 6km Demo Data
#'
#' This data set contains network data for Edinburgh within a 6km radius, sourced from NPT.
#' Useful for analyses that require detailed network data and planning metrics.
#'
#' @docType data
#' @keywords datasets
#' @name NPT_demo_6km
#' @format An object of class \code{sf} (inherits from \code{data.frame}).
#' @examples
#' data(NPT_demo_6km)
#' head(NPT_demo_6km)
"NPT_demo_6km"
NULL

#' Central Leeds OSM Network
#' 
#' See the `data-raw` folder for the code used to generate this data set.
#' 
#' @docType data
#' @keywords datasets
#' @name central_leeds_osm
#' @format An object of class \code{sf} (inherits from \code{data.frame}).
#' @examples
#' head(central_leeds_osm)
#' library(sf) # for plotting
#' plot(central_leeds_osm$geometry)
NULL