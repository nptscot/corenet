#' Aggregate a Source Attribute onto Target Features Using 'anime'
#'
#' This function runs [anime::anime()] to match source features to target features
#' based on angle and distance tolerances. The matched source attribute(s) are then
#' aggregated (sum, max, mean, etc.) and joined back to the target data.
#'
#' @param source_data      An \code{sf} object representing the source data.
#' @param target_data      An \code{sf} object representing the target data.
#' @param attribute        A character string, the name of the source attribute to aggregate.
#' @param new_name         A character string for the output column name in the returned data frame.
#' @param agg_fun          A function used to aggregate the matched values (e.g., \code{sum}, \code{mean}, \code{max}).
#'                         Defaults to \code{sum}.
#' @param weights          A character vector of column names in the source data to be used as weights.
#'                         (For example, \code{"target_weighted"}).
#' @param angle_tolerance  A numeric value for angular tolerance used by \code{anime}. Defaults to 35.
#' @param distance_tolerance A numeric value for distance tolerance used by \code{anime}. Defaults to 15.
#'
#' @return A \code{data.frame} (or \code{tibble}) containing:
#' \itemize{
#'   \item \code{id}: The original \code{id} from \code{target_data}.
#'   \item \code{new_name}: The aggregated attribute (named according to \code{new_name}).
#' }
#'
#' @examples
#' \dontrun{
#'   result = anime_join(
#'     source_data      = my_source_sf,
#'     target_data      = my_target_sf,
#'     attribute        = "population",
#'     new_name         = "pop_sum",
#'     agg_fun          = sum,
#'     weights          = "target_weighted",
#'     angle_tolerance  = 35,
#'     distance_tolerance = 15
#'   )
#' }
#'
#' @export
#'
#' @importFrom dplyr mutate left_join group_by summarise row_number select any_of
#' @importFrom sf st_drop_geometry
#' @importFrom rlang sym syms expr
#' @importFrom purrr reduce
#' @importFrom anime anime get_matches
anime_join = function(
  source_data,
  target_data,
  attribute,
  new_name,
  agg_fun = sum,
  weights,
  angle_tolerance = 35,
  distance_tolerance = 15
) {
  # Helper to display the name of the aggregation function
  get_fun_name = function(fn) {
    # If the function is one of the standard ones, just return its name.
    # Otherwise, you can fallback to something generic, e.g. "user-defined function".
    if (identical(fn, sum)) return("sum")
    if (identical(fn, mean)) return("mean")
    if (identical(fn, max)) return("max")
    # fallback
    return("user-defined function")
  }

  # Let user know what's happening
  message(
    "Aggregating attribute: '", attribute,
    "' using function: '", get_fun_name(agg_fun), "'."
  )

  # 1. Run anime matching between source and target
  matches = anime::anime(
    source            = source_data,
    target            = target_data,
    angle_tolerance   = angle_tolerance,
    distance_tolerance = distance_tolerance
  )
  matches_df = anime::get_matches(matches)

  # 2. Add a source_id to source_data, drop geometry, and join with matches
  net_source_matches = dplyr::mutate(source_data, source_id = dplyr::row_number()) |>
    sf::st_drop_geometry() |>
    dplyr::left_join(matches_df, by = "source_id")

  # 3. Convert the attribute and weight names to symbols
  attr_sym    = rlang::sym(attribute)
  weight_syms = rlang::syms(weights)

  # 4. Create the expression for aggregation
  #    - If agg_fun is max, we don't multiply by weights.
  if (identical(agg_fun, max)) {
    mult_expr = purrr::reduce(
      weight_syms,
      function(acc, w) rlang::expr((!!acc)),  # no-op for max
      .init = attr_sym
    )
  } else {
    # For sum, mean, etc., multiply the attribute by the provided weights
    mult_expr = purrr::reduce(
      weight_syms,
      function(acc, w) rlang::expr((!!acc) * (!!w)),
      .init = attr_sym
    )
  }

  # 5. Group by target_id and summarize
  net_target_aggregated = net_source_matches |>
    dplyr::group_by(row_number = target_id) |>
    dplyr::summarise(
      !!rlang::sym(new_name) := agg_fun(!!mult_expr, na.rm = TRUE),
      .groups = "drop"
    )

  # 6. If aggregator is max, round the result
  if (identical(agg_fun, max)) {
    net_target_aggregated = dplyr::mutate(
      net_target_aggregated,
      !!rlang::sym(new_name) := round(!!rlang::sym(new_name))
    )
  }

  # 7. Join the aggregated values back to the target data
  #    - Drop geometry, remove any existing 'new_name' column, then join
  net_target_joined = target_data |>
    dplyr::mutate(row_number = dplyr::row_number()) |>
    sf::st_drop_geometry() |>
    dplyr::select(-dplyr::any_of(new_name)) |>
    dplyr::left_join(net_target_aggregated, by = "row_number") |>
    dplyr::select(id, !!rlang::sym(new_name))

  return(net_target_joined)
}
