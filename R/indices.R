# ==========================================================================
# MineralMappeR: Core Processing & Visualization Engine
# Application: Landsat 9 Multi-spectral Analysis for Hydrothermal Alteration
# ==========================================================================



#' Calibrate Landsat 9 Digital Numbers (DN) to Surface Reflectance
#'
#' Applies the official USGS Collection 2 Level-2 scale factor (0.0000275) and
#' additive offset (-0.2). These coefficients are valid for Level-2 surface
#' reflectance bands only (SR_B2 ... SR_B7); Level-1 top-of-atmosphere DN use
#' different rescaling coefficients (see the scene MTL file).
#'
#' @param image A SpatRaster of Landsat 9 Collection 2 Level-2 surface
#'   reflectance bands (digital numbers).
#' @return A SpatRaster with surface reflectance values clamped to [0, 1].
#'   Small negative values, which occur in dark pixels of surface reflectance
#'   products, are set to 0.
#' @export
scale_landsat_data <- function(image) {
  # Apply official Collection 2 scaling factors and offsets
  scaled_image <- (image * 0.0000275) - 0.2

  # Clamp values to ensure physical validity within the [0, 1] range
  return(terra::clamp(scaled_image, lower = 0, upper = 1))
}

#' Integrated Atmospheric and Surface Noise Masking
#'
#' @param image Multispectral SpatRaster (surface reflectance).
#' @param qa_layer Landsat 9 QA_PIXEL layer for bitmasking (cloud, cloud shadow
#'   and water bits).
#' @param albedo_threshold Threshold on the mean reflectance across all bands
#'   of `image`, above which pixels are removed as bright surfaces such as salt
#'   flats (default 0.5).
#' @return A SpatRaster with noisy pixels set to NA.
#' @export
clean_image_noise <- function(image, qa_layer, albedo_threshold = 0.5) {
  # Bitmasking for Landsat Collection 2: Clouds (8), Shadows (16), Water (128). Sum = 152
  qa_mask <- terra::app(qa_layer, function(x) bitwAnd(x, 152) != 0)

  # Surface Masking: Identify high-reflectance anomalies (e.g., Salars)
  albedo_mask <- terra::mean(image) > albedo_threshold

  # Merge atmospheric and surface masks into a single noise layer
  noise_mask <- qa_mask | albedo_mask

  return(terra::mask(image, noise_mask, maskvalues = TRUE))
}

#' Environmental Interference Masking (Vegetation and Water)
#'
#' @param image A SpatRaster of surface reflectance containing layers named
#'   "B4" (Red) and "B5" (NIR).
#' @param threshold NDVI limit for masking vegetation (default 0.3). Pixels
#'   with NDVI above this value, or below 0 (water, deep shadow), are removed.
#' @return A SpatRaster with vegetated and water/shadow pixels set to NA.
#' @export
mask_env_noise <- function(image, threshold = 0.3) {
  # Calculate Normalized Difference Vegetation Index (NDVI)
  ndvi <- (image[["B5"]] - image[["B4"]]) / (image[["B5"]] + image[["B4"]])

  # Mask high NDVI (vegetation) and residual negative NDVI (water/shadows)
  env_mask <- (ndvi > threshold) | (ndvi < 0)

  return(terra::mask(image, env_mask, maskvalues = TRUE))
}

#' Calculate Diagnostic Spectral Mineral Indices
#'
#' @param image A SpatRaster of surface reflectance with layers named "B2",
#'   "B4", "B6" and "B7".
#' @return A SpatRaster with the layers `Iron_Oxide` (B4/B2) and
#'   `Clay_Minerals` (B6/B7). Non-finite values (division by zero, which occurs
#'   where reflectance was clamped to 0) are set to NA.
#' @export
calc_mineral_indices <- function(image) {
  # Mineral Ratios: Iron Oxide (Red/Blue) and Clay Minerals (SWIR1/SWIR2)
  iron_oxide <- image[["B4"]] / image[["B2"]]
  clay_minerals <- image[["B6"]] / image[["B7"]]

  results <- c(iron_oxide, clay_minerals)
  names(results) <- c("Iron_Oxide", "Clay_Minerals")

  # Division by zero (reflectance clamped to 0) gives Inf; set to NA so that
  # these pixels cannot be flagged as targets
  results[!is.finite(results)] <- NA

  return(results)
}

#' Identify High-Probability Mineral Exploration Targets
#'
#' @param mineral_raster SpatRaster with calculated mineral indices.
#' @param threshold Anomaly cut-off applied to both indices (default 1.5). This
#'   is an empirical starting value, not a universal standard, and should be
#'   tuned to the scene (sensor, atmosphere, lithology).
#' @return A SpatRaster named `Exploration_Targets` (1 = target, 0 = background).
#' @export
find_mineral_targets <- function(mineral_raster, threshold = 1.5) {
  # Logical intersection of Iron Oxide and Clay Mineral anomalies
  targets <- ((mineral_raster[["Iron_Oxide"]] > threshold) &
                (mineral_raster[["Clay_Minerals"]] > threshold)) * 1

  names(targets) <- "Exploration_Targets"
  return(targets)
}

#' Plot Mineral Diagnostic Maps
#'
#' @param mineral_raster A SpatRaster from calc_mineral_indices. Values above
#'   2.5 are displayed in the top colour of the scale.
#' @return A ggplot object.
#' @import ggplot2
#' @import tidyterra
#' @export
plot_mineral_indices <- function(mineral_raster) {
  ggplot() +
    geom_spatraster(data = mineral_raster) +
    facet_wrap(~lyr) +
    scale_fill_whitebox_c(palette = "viridi", limits = c(0, 2.5),
                          oob = scales::squish) +
    theme_minimal() +
    labs(
      title = "Mineral Alteration Mapping",
      subtitle = "Diagnostic Ratios: Iron Oxide and Clay Minerals",
      fill = "Index Value"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
}

#' Generate Professional Exploration Target Plot
#'
#' @param target_raster A SpatRaster containing identified targets
#'   (output of find_mineral_targets).
#' @return A ggplot object.
#' @import ggplot2
#' @import tidyterra
#' @import ggspatial
#' @export
plot_mineral_targets <- function(target_raster) {
  ggplot() +
    geom_spatraster(data = target_raster) +
    scale_fill_gradient(
      low = "#2F4F4F", high = "gold",
      name = "Target Status", breaks = c(0, 1), labels = c("Background", "Target")
    ) +
    annotation_north_arrow(
      location = "tr", height = unit(0.8, "cm"), width = unit(0.8, "cm"),
      style = north_arrow_fancy_orienteering(fill = c("grey20", "white"))
    ) +
    annotation_scale(location = "bl", width_hint = 0.2) +
    theme_minimal() +
    labs(
      title = "Final Exploration Targets",
      subtitle = "Hydrothermal Alteration zones: Iron Oxide & Clay Minerals",
      caption = "Data Source: Landsat 9 | Processed with MineralMappeR"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid = element_blank(),
      axis.title = element_blank()
    )
}

#' Batch Export Spectral Products to GeoTIFF
#'
#' @param mineral_raster SpatRaster to export.
#' @param output_dir Directory path for the output files (created if missing).
#' @return Invisibly returns NULL; called for its side effect of writing files.
#' @export
export_mineral_maps <- function(mineral_raster, output_dir = "Results") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  for (i in names(mineral_raster)) {
    terra::writeRaster(mineral_raster[[i]],
                       filename = file.path(output_dir, paste0(i, "_Map.tif")),
                       overwrite = TRUE)
  }
  message("Success: All maps exported to ", output_dir)
  invisible(NULL)
}
