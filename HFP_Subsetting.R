# 1. LOAD LIBRARIES ------------------------------------------------------------
pacman::p_load(
  sf,
  terra,
  tidyverse,
  rnaturalearth
)

# 2. LOAD WASHINGTON STATE BOUNDARY --------------------------------------------
# Load Washington boundary in Mollweide
wa_boundary <- ne_states(
  country = "United States of America",
  returnclass = "sf") %>%
  filter(name == "Washington") %>%
  st_transform("ESRI:54009")  # Convert to Mollweide

# 4. IDENTIFY RELEVANT HFP TILES -----------------------------------------------
# Function to check tile overlap with study area
check_overlap <- function(tif_file, boundary) {
  r <- rast(tif_file)
  r_extent <- ext(r)
  r_poly <- st_as_sfc(st_bbox(r_extent)) %>% st_set_crs(crs(r))
  return(st_intersects(r_poly, boundary, sparse = FALSE)[1,1])
}

# Get list of all HFP tiles
tif_files <- list.files("data/hfp-100m-2020/", pattern = "\\.tif$", full.names = TRUE)

# Find tiles intersecting our buffer
overlap_results <- sapply(tif_files, check_overlap, boundary = wa_boundary)
relevant_files <- tif_files[overlap_results]

# 5. PROCESS HFP DATA ----------------------------------------------------------
# Load, crop, and merge relevant tiles
hfp_combined <- relevant_files %>%
  map(rast) %>%             # Convert each file to a raster (no cropping)
  do.call(merge, .) %>%
  crop(vect(wa_boundary)) %>% 
  mask(vect(wa_boundary))

#Length of 1 degree longitude at latitude φ ≈ 111,320 × cos(φ) meters
#cos(48°) ≈ 0.6691
#So 1° longitude ≈ 111,320 × 0.6691 ≈ 74,500 meters
#1° latitude is still ≈ 111,000 meters
# Divide 100 by 111,000 and 74,500 to match resolution

hfp_wgs84 <- project(
  hfp_combined, "EPSG:4326",
  method = "near", #Using nearest neighbor method for data integrity
  res = c(0.00134, 0.0009)) # (x, y) = (long, lat)

# Quick verification plot
plot(hfp_wgs84)

# Save the final cropped HFP raster
writeRaster(
  hfp_wgs84,
  "data/HFP_washington.tif",
  overwrite = TRUE,
  datatype = "FLT4S",          # Maintain decimal precision for HFP values
  gdal = c("COMPRESS=LZW")     # Reduce file size
)
 