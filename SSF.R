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
library(glmmTMB)
library(DHARMa)

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

png("img/bobcat_coyote_mean_sampling_rates.png", width = 800, height = 600)
boxplot(mean ~ species, data = trackSummary, 
        xlab = "", 
        ylab = "Mean sampling interval in hours")
dev.off()

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

coyote_extracted <- read_delim("data/coyote_extracted.csv") %>%
  dplyr::select(-log_sl_, -sl_) %>%
  mutate(# step lengths from amt are in degrees; overwrite to meters
    sl_ = distGeo(.[,c("x1_","y1_")],.[,c("x2_","y2_")]),
    log_sl_ = log(sl_),
    # instead of using calendar years, we start the year on 1 Dec of the prior year
    analysis_year = ifelse(month(t2_) == 12, year(t2_) + 1, year(t2_)),
    season = ifelse(month(t2_) %in% 4:11, "summer", "winter")) 

bobcat_extracted <- read_delim("data/bobcat_extracted.csv") %>%
  dplyr::select(-sl_) %>%
  mutate(# step lengths from amt are in degrees; overwrite to meters
    sl_ = distGeo(.[,c("x1_","y1_")],.[,c("x2_","y2_")]),
    log_sl_ = log(sl_),
    # instead of using calendar years, we start the year on 1 Dec of the prior year
    analysis_year = ifelse(month(t2_) == 12, year(t2_) + 1, year(t2_)),
    season = ifelse(month(t2_) %in% 4:11, "summer", "winter")) 


# Load and preprocess rasters --------------------------------------------------
hfp <- rast("data/HFP_washington.tif")
NAflag(hfp) <- 64536  # Set no-data value
hfp_capped <- classify(hfp, matrix(c(50000, Inf, 50000), ncol = 3, byrow = TRUE)) # Cap at 50k
hfp_scaled <- app(hfp_capped, \(x) (x / 50000) - 0.5)  # Scale to -0.5–0.5
                       
land_use <- rast("data/ESA_washington.tif")  # Ensure CRS matches tracking data
                       
# Land use class labels (ESA WorldCover 2021)
class_labels <- c(
  "10" = "Tree cover", "20" = "Shrubland", "30" = "Grassland",
  "40" = "Cropland", "50" = "Built-up", "60" = "Bare/sparse vegetation",
  "70" = "Snow and Ice", "80" = "Permanent water bodies",
  "90" = "Herbaceous wetland", "95" = "Mangroves", "100" = "Moss and lichen"
)

# Extract covariates -----------------------------------------------------------
extract_covariates <- function(data) {
  data %>%
    mutate(
      human_footprint = terra::extract(hfp_scaled, cbind(x2_, y2_))[, 1],
      land_use_code = terra::extract(land_use, cbind(x2_, y2_))[, 1]
    ) %>%
    mutate(
      land_use = factor(land_use_code, 
                        levels = names(class_labels),
                        labels = class_labels)
    )
}
  
coyote_extracted1 <- extract_covariates(coyote_extracted)
bobcat_extracted1 <- extract_covariates(bobcat_extracted)

# Final formatting --------------------------------------------------------
format_for_ssf <- function(data) {
  data %>%
    mutate(
      step_id_ = paste(animal_id, step_id_, sep = "_")) %>% # Unique stratum ID
    group_by(animal_id) %>%
    mutate(n = n()/11) %>%  # 1 used + 10 available steps per stratum
    ungroup()
}

coyote_final <- format_for_ssf(coyote_extracted1)
bobcat_final <- format_for_ssf(bobcat_extracted1)

write_csv(coyote_final, "data/coyote_ssf_data.csv")
write_csv(bobcat_final, "data/bobcat_ssf_data.csv")

# Coyote Step Selection Function ---------------------------------------------

coyote_ssf <- glmmTMB(
  case_ ~ 
    human_footprint +          # Fixed effect of human footprint
    land_use +                 # Fixed effect of land use (categorical)
    (1 + human_footprint | animal_id) +  # Random slopes: individual variation in HFP response
    (1 | step_id_) +           # Stratum-specific intercept (matched used/available)
    offset(log_sl_),            # Control for step length
  family = poisson,            # Likelihood-equivalent to conditional logistic
  data = coyote_final
)

coyote_ssf <- glmmTMB(
  case_ ~ 
    # Habitat selection under human influence:
    land_use * human_footprint +          # Interaction: land use × human footprint
    I(human_footprint^2) +                # Non-linear human footprint effect
    # Random effects:
    (1 | animal_id) +                     # Individual baseline selection
    (1 | step_id_) +                      # Stratum (matched steps)
    offset(log_sl_),                      # Step length control
  family = poisson,
  data = coyote_final
)

summary(coyote_ssf)

sim_res <- simulateResiduals(coyote_ssf)
plot(sim_res)
