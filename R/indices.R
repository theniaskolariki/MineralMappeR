#' Calculate Mineral Indices from Landsat 9
#'
#' @param image A SpatRaster object with 4 bands (B2, B4, B6, B7)
#' @return A SpatRaster with Iron and Clay indices
#' @export
calc_mineral_indices <- function(image) {

  # Iron Oxide Index (Red / Blue)
  iron_oxide <- image[[2]] / image[[1]]

  # Clay Minerals Index (SWIR1 / SWIR2)
  clay_minerals <- image[[3]] / image[[4]]

  # Ενώνουμε τα αποτελέσματα
  results <- c(iron_oxide, clay_minerals)
  names(results) <- c("Iron_Oxide", "Clay_Minerals")

  return(results)
}
