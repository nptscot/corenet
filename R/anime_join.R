# Define a helper function that runs anime and aggregates an attribute.
anime_join = function(source_data,
                       target_data,
                       attribute,      # character name of the attribute to aggregate
                       new_name,       # character name for the output column
                       agg_fun = sum,  # aggregation function (e.g., sum, mean, max)
                       weights,        # character vector of weight columns (e.g., "target_weighted")
                       angle_tolerance = 35,
                       distance_tolerance = 15) {
                        library(dplyr)
  library(sf)
  library(rlang)
  library(purrr)
  library(anime)  
  
  get_fun_name = function(fn) {
    deparse(substitute(fn))
  }

  cat("Aggregating attribute:", attribute, "using function:", get_fun_name(sum), "\n")

  # Run anime matching between source and target
  matches = anime::anime(
    source = source_data,
    target = target_data,
    angle_tolerance = angle_tolerance,
    distance_tolerance = distance_tolerance
  )
  matches_df = anime::get_matches(matches)
  
  # Join the source data (with a source_id) to the matches.
  net_source_matches = source_data %>%
    mutate(source_id = row_number()) %>%
    st_drop_geometry() %>%
    left_join(matches_df, by = "source_id")
  
  # Convert the attribute and weight names to symbols.
  attr_sym = sym(attribute)
  weight_syms = syms(weights)
  
  # For max aggregatio
  if (identical(agg_fun, max)) {
    mult_expr = purrr::reduce(weight_syms, function(acc, w) expr((!!acc)), .init = attr_sym)
  } else {
    mult_expr = purrr::reduce(weight_syms, function(acc, w) expr((!!acc) * (!!w)), .init = attr_sym)
  }
  
  # Group by target id (assumed to be provided by anime as "target_id") and aggregate.
  net_target_aggregated = net_source_matches %>%
    group_by(row_number = target_id) %>%
    summarise(!!sym(new_name) := agg_fun(!!mult_expr, na.rm = TRUE))
  
  # For max, we also want to round the final result.
  if (identical(agg_fun, max)) {
    net_target_aggregated = net_target_aggregated %>%
      mutate(!!sym(new_name) := round(!!sym(new_name)))
  }
  
  # Join the aggregated values back to the target data.
  net_target_joined = target_data %>%
    mutate(row_number = row_number()) %>%
    st_drop_geometry() %>%
    select(-any_of(new_name)) %>%   # Remove the original column (if it exists)
    left_join(net_target_aggregated, by = "row_number") %>%
    select(id, !!sym(new_name))
  
  return(net_target_joined)
}
