# MineralMappeR

`MineralMappeR` is an R package designed for regional-scale mineral exploration targeting using Landsat 9 multispectral imagery. It automates the screening and detection of hydrothermal alteration zones associated with porphyry copper-gold systems by calculating diagnostic spectral indices while filtering out atmospheric, surface, and environmental noise.

## Features

-   **Radiometric Calibration:** Converts raw Landsat 9 DN values to Top-of-Atmosphere (TOA) Reflectance automatically.
-   **Multi-Tier Noise Masking:** Integrated QA bitmasking (clouds/shadows/water), Albedo filtering (to exclude highly reflective salt flats/Salars), and NDVI thresholding (to eliminate vegetation and topographic shadows).
-   **Diagnostic Mineral Ratios:** Computes Iron Oxide (B4/B2) and Clay Minerals (B6/B7) indices based on classic Rowan and Sabins remote sensing methodologies.
-   **Exploration Targeting:** Logically intersects anomalous pixel values to highlight high-probability alteration targets.
-   **Geospatial Visualization:** Built-in professional plotting capabilities optimized with `ggplot2` and `tidyterra`.

For a detailed, step-by-step scientific walkthrough and geological analysis, please refer to the package vignette: vignette("MineralMappeR_Tutorial", package = "MineralMappeR").

## Installation

You can install the development version of `MineralMappeR` directly from GitHub using the `devtools` package:

``` r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install MineralMappeR from GitHub
devtools::install_github("theniaskolariki/MineralMappeR")
```

## Quick Start

Here is a quick example of how to run the full processing pipeline using the built-in sample data from the Atacama Desert, Chile:

``` r
library(MineralMappeR)
library(terra)

# 1. Load built-in sample data using internal package paths
band_paths <- system.file("extdata", paste0("atacam_sample_B", c(2,4,5,6,7), ".tif"), package = "MineralMappeR")
image <- terra::rast(band_paths)
names(image) <- c("B2", "B4", "B5", "B6", "B7")

qa_path <- system.file("extdata", "atacam_sample_QA.tif", package = "MineralMappeR")
qa_layer <- terra::rast(qa_path)

# 2. Run the processing and masking pipeline
image_ref   <- scale_landsat_data(image)
image_clean <- clean_image_noise(image_ref, qa_layer, albedo_threshold = 0.5)
image_final <- mask_env_noise(image_clean, threshold = 0.3)

# 3. Calculate indices and identify targets
results     <- calc_mineral_indices(image_final)
targets     <- find_mineral_targets(results, threshold = 1.5)

# 4. Plot final exploration targets
plot_mineral_targets(targets)
```

## License

This project is licensed under the MIT License - see the DESCRIPTION file for details.
