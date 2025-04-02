# 1. LOAD LIBRARIES ------------------------------------------------------------
pacman::p_load(
  sf,
  terra,
  tidyverse
)

# 2. PREPARE TRACKING DATA -----------------------------------------------------
# Read and prepare tracking data
tracking_data <- read_delim("bobcat_coyotes_wa_gps.csv")

# Convert to SF object and project to Mollweide
animal_sf <- tracking_data %>%
  st_as_sf(
    coords = c("location-long", "location-lat"),
    crs = 4326  # WGS84
  ) %>%
  st_transform("ESRI:54009")  # Mollweide projection

# 3. CREATE STUDY AREA BUFFER --------------------------------------------------
# Generate 200km buffer around tracking points
tracking_buffer <- animal_sf %>%
  st_union() %>%
  st_convex_hull() %>%    # Create minimum convex polygon
  st_buffer(200000)       # 200km buffer in meters (Mollweide uses meters)

# 4. IDENTIFY RELEVANT HFP TILES -----------------------------------------------
# Function to check tile overlap with study area
check_overlap <- function(tif_file, boundary) {
  r <- rast(tif_file)
  r_extent <- ext(r)
  r_poly <- st_as_sfc(st_bbox(r_extent)) %>% st_set_crs(crs(r))
  return(st_intersects(r_poly, boundary, sparse = FALSE)[1,1])
}

# Get list of all HFP tiles
tif_files <- list.files("HFP-100m-2020/", pattern = "\\.tif$", full.names = TRUE)

# Find tiles intersecting our buffer
overlap_results <- sapply(tif_files, check_overlap, boundary = tracking_buffer)
relevant_files <- tif_files[overlap_results]

# 5. PROCESS HFP DATA ----------------------------------------------------------
# Load, crop, and merge relevant tiles
hfp_cropped <- relevant_files %>%
  map(~crop(rast(.x), tracking_buffer)) %>%
  do.call(merge, .)

# Save the final cropped HFP raster
writeRaster(
  hfp_cropped,
  "HFP_washington_mollweide.tif",
  overwrite = TRUE,
  datatype = "FLT4S",          # Maintain decimal precision for HFP values
  gdal = c("COMPRESS=LZW")     # Reduce file size
)

# 6. VISUAL CHECK --------------------------------------------------------------
# Quick verification plot
plot(hfp_cropped, main = "Human Footprint Index (Cropped)")
plot(st_geometry(tracking_buffer), add = TRUE, border = "red")
plot(st_geometry(animal_sf), add = TRUE, col = "blue", pch = 16, cex = 0.5)