#' Generate a core network from road network data
#' 
#' Description to go here...
#' 
#' @param net A road network object, of class sf
#' @param tbc An object to be determined
#' @return A core network object, of class sf
#' @export
#' @examples
#' corenet()
corenet = function(net, tbc) {
    data.frame(x = 1, y = 2) |>
      utils::head()
}



#' Data from edinburgh's OSM network
#'
#'
#' @docType data
#' @keywords datasets
#' @name osm_edinburgh_demo
#' @format An sf data frame
#' @examples 
#' library(sf)
#' names(osm_edinburgh_demo)
#' head(osm_edinburgh_demo)
#' plot(osm_edinburgh_demo)
NULL