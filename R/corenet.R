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
#' # TODO: fix this
#' # data(os_edinburgh_demo)
#' # names(os_edinburgh_demo)
#' os_edinburgh = get_edinburgh_6km()
#' head(os_edinburgh)
NULL