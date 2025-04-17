# Load necessary packages
library(amt)
library(lubridate)
library(dplyr)
library(sf)
library(raster)
library(tidyverse)
library(leaflet)
library(geosphere)
library(terra)

options(scipen = 999) # turn off scientific notation
options(digits = 15) # set digits to 15 to ensure GPS coordinates aren't truncated
sf_use_s2(FALSE) # workaround for one of the amt functions 

# Read, rename, select, and filter tracking data
tracking_data <- read_delim("data/bobcat_coyotes_wa_gps.csv") %>%
  dplyr::rename(
    id = `event-id`,
    long = `location-long`, 
    lat = `location-lat`,
    animal_id = `individual-local-identifier`,
    species = `individual-taxon-canonical-name`
  ) %>%
  dplyr::arrange(animal_id, timestamp) %>%
  dplyr::select(animal_id, species, timestamp, lat, long) %>%
  dplyr::filter(
    !animal_id %in% c(
      "MVBOB87M", "NECOY12F", "MVCOY65M", "MVBOB88M", 
      "NEBOB8M", "MVBOB77M", "MVBOB91M", "NEBOB11M", 
      "NECOY42F", "MVCOY72F", "NECOY14M"
    ),
    !(animal_id == "MVBOB71M" & timestamp > as.POSIXct("2019-09-24 00:00:00"))
  )


# plot a random subset of the locations (to reduce processing time)
tracking_for_plot <- tracking_data %>% mutate(rand = round(runif(nrow(tracking_data), 0, 5)), mock = 1) %>% filter(rand == 1) 

cols25 <- function() {
  c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", 
    "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6",
    "#ffff99", "#b15928", "#8dd3c7", "#fb8072", "#80b1d3",
    "#bebada", "#ffed6f", "#fdb462", "#b3de69", "#fccde5",
    "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#7FC97F")
}

leaflet() %>% 
  addProviderTiles(providers$Esri.NatGeoWorldMap) %>%
  addCircles(data = tracking_data , 
             color = ~ colorFactor(cols25(), domain = factor(animal_id))(factor(animal_id)),
             label = paste(tracking_data$animal_id,
                           tracking_data$timestamp), opacity = 1) 


# Select variables, nest by individual animal, and make 'amt' track
track <- tracking_data %>%
  data.frame() %>%
  nest(data = c(-animal_id, -species)) %>%
  mutate(trk = lapply(data, function(d) {
    make_track(d, .x = long, .y = lat, .t = timestamp, crs = 4326)
  }))


# Summarize sampling rate
trackSummary <- track %>% mutate(sr = lapply(trk, summarize_sampling_rate, time_unit = "hour")) %>%
  dplyr::select(animal_id, sr) %>% unnest %>% 
  left_join(distinct(data.frame(track)[,c("animal_id", "species")])) %>%
  arrange(species, median)

print(trackSummary, n = 70)

# separate by species
coyote <- track %>% filter(grepl("COY", animal_id))
bobcat <- track %>% filter(grepl("BOB", animal_id))

coyote1 <- coyote %>%
  mutate(stp = map(trk, function(x)
    x %>% track_resample(rate = hours(4), tolerance = minutes(10)) %>%
      steps_by_burst() %>% 
      random_steps(n_control = 10) %>% 
      mutate(log_sl_ = log(sl_ + 1), cos_ta_ = cos(ta_)))) %>%  # Distribute 10 available points for each used location
  dplyr::select(-data, -trk) %>%
  unnest(cols = stp) %>%
  mutate(case_binary = ifelse(case_ == TRUE, 1, 0))

bobcat1 <- bobcat %>%
  mutate(stp = map(trk, function(x)
    x %>% track_resample(rate = hours(4), tolerance = minutes(10)) %>%
      steps_by_burst() %>% 
      random_steps(n_control = 10) %>% 
      mutate(log_sl_ = log(sl_ + 1), cos_ta_ = cos(ta_)))) %>%  # Distribute 10 available points for each used location
  dplyr::select(-data, -trk) %>%
  unnest(cols = stp) %>%
  mutate(case_binary = ifelse(case_ == TRUE, 1, 0))

write_csv(coyote1, "data/coyote_extracted.csv")
write_csv(bobcat1, "data/bobcat_extracted.csv")

coyote_extracted <- read_delim("coyote_extracted_20250415.csv") %>%
  dplyr::select(-log_sl_, -sl_) %>%
  mutate(# step lengths from amt are in degrees; overwrite to meters
    sl_ = distGeo(.[,c("x1_","y1_")],.[,c("x2_","y2_")]),
    log_sl = log(sl_),
    # instead of using calendar years, we start the year on 1 Dec of the prior year
    analysis_year = ifelse(month(t2_) == 12, year(t2_) + 1, year(t2_)),
    season = ifelse(month(t2_) %in% 4:11, "summer", "winter")) 

bobcat_extracted <- read_delim("bobcat_extracted_20250415.csv") %>%
  dplyr::select(-sl_) %>%
  mutate(# step lengths from amt are in degrees; overwrite to meters
    sl_ = distGeo(.[,c("x1_","y1_")],.[,c("x2_","y2_")]),
    log_sl = log(sl_),
    # instead of using calendar years, we start the year on 1 Dec of the prior year
    analysis_year = ifelse(month(t2_) == 12, year(t2_) + 1, year(t2_)),
    season = ifelse(month(t2_) %in% 4:11, "summer", "winter")) 


land_cover <- raster("ESA_land_cover_washington.tif") %>% project(hfp)

hfp <- raster("HFP_washington_mollweide.tif")
NAvalue(hfp) <- 64536

# Cap values > 50000 and scale
hfp_scaled <- calc(hfp, fun = function(x) {
  x[x > 50000] <- 50000
  x / 1000
})

