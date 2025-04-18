# 1. LOAD LIBRARIES ------------------------------------------------------------
pacman::p_load(
  rstac,       # For accessing STAC APIs
  sf,
  terra,
  tidyverse,
  tidyterra,
  rnaturalearth,
  ggplot2
)
# 2. LOAD WASHINGTON STATE BOUNDARY --------------------------------------------
# Load Washington boundary in WGS84
wa_boundary <- ne_states(
  country = "United States of America",
  returnclass = "sf") %>%
  filter(name == "Washington") %>%
  st_transform("EPSG:4326")  # Convert to WGS84

# 3. DOWNLOAD ESA WORLDCOVER DATA ----------------------------------------------
# Create 'data/esa-worldcover' folder
output_dir <- "data/esa-worldcover"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Connect to Microsoft Planetary Computer STAC
stac_api <- stac("https://planetarycomputer.microsoft.com/api/stac/v1")

# Search for ESA WorldCover 2021 data within study area
esa_items <- stac_search(
  q = stac_api,
  collections = "esa-worldcover",
  datetime = "2021-01-01/2021-12-31",  # 2021 version
  bbox = st_bbox(wa_boundary),           # WA bounding box
  limit = 100
) %>%
  get_request() %>%
  items_sign(sign_planetary_computer())  # Authenticate

# Download all assets
assets_download(
  items = esa_items,
  asset_names = "map",
  output_dir = "data/",
  overwrite = TRUE
)

# 4. PROCESS ESA DATA ----------------------------------------------------------
# Get all downloaded tiles
esa_tiles <- list.files(
  path = file.path("data/esa-worldcover/v200/2021/map"),
  pattern = "\\.tif$",
  full.names = TRUE
)

# Merge and crop tiles to study area
land_cover <- esa_tiles %>%
  lapply(rast) %>%
  do.call(merge, .) %>%    # Combine all tiles
  crop(vect(wa_boundary)) %>% 
  mask(vect(wa_boundary))

# Save the final cropped/merged land cover
writeRaster(
  land_cover,
  "data/ESA_washington.tif",   # Save inside 'data'
  overwrite = TRUE,            # Overwrite if exists
  datatype = "INT1U",          # Optimize for categorical data
  gdal = c("COMPRESS=LZW")     # Reduce file size
)

# 5. VISUALIZE RESULTS ---------------------------------------------------------

# Convert to factor for plotting
land_cover <- as.factor(land_cover)

# Class labels (ESA WorldCover 2021)
class_labels <- c(
  "10" = "Tree cover",
  "20" = "Shrubland",
  "30" = "Grassland",
  "40" = "Cropland",
  "50" = "Built-up",
  "60" = "Bare/sparse vegetation",
  "70" = "Snow and ice",
  "80" = "Permanent water bodies",
  "90" = "Herbaceous wetland",
  "95" = "Mangroves",
  "100" = "Moss and lichen"
)

# Create plot
ggplot() +
  tidyterra::geom_spatraster(
    data = land_cover,
    maxcell = 1e6
  ) +
  scale_fill_manual(
    values = c(
      "10" = "#006400", "20" = "#FFBB22", "30" = "#FFFF4C",
      "40" = "#F096FF", "50" = "#FA0000", "60" = "#B4B4B4",
      "70" = "#F0F0F0", "80" = "#0064C8", "90" = "#0096A0",
      "95" = "#00CF75", "100" = "#FAE6A0"
    ),
    labels = class_labels,
    na.value = NA
  ) +
  labs(
    title = "Washington, NA",
    fill = "Land cover class"
  ) +
  theme_minimal()
ggsave("img/esa_landcover_washington.png")
