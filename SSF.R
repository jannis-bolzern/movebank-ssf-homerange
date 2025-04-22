# Load necessary packages -----------------------------------------------------
pacman::p_load(
  tidyverse,
  amt,
  sf,
  geosphere,
  terra,
  glmmTMB,
  DHARMa,
  emmeans
)

options(scipen = 999)
options(digits = 15)

# Load and prepare tracking data ---------------------------------------------
tracking_data <- read_delim("data/bobcat_coyotes_wa_gps.csv") |> 
  dplyr::rename(
    long = `location-long`, 
    lat = `location-lat`,
    id = `individual-local-identifier`,
    timestamp = `timestamp`,
    species = `individual-taxon-canonical-name`) |> 
  dplyr::arrange(id, timestamp) |> 
  dplyr::select(id, species, timestamp, lat, long) |> 
  dplyr::filter(
    !id %in% c(
      "MVBOB87M", "NECOY12F", "MVCOY65M", "MVBOB88M", 
      "NEBOB8M", "MVBOB77M", "MVBOB91M", "NEBOB11M", 
      "NECOY42F", "MVCOY72F", "NECOY14M"),
    !(id == "MVBOB71M" & timestamp > as.POSIXct("2019-09-24 00:00:00"))
  )

# Create and inspect tracks ---------------------------------------------------
track <- tracking_data |>
  nest(data = c(-id, -species)) |>
  mutate(trk = map(data, ~ make_track(.x, long, lat, timestamp, crs = 4326)))

# Summarize sampling rate
trackSummary <- track |> 
  mutate(sr = lapply(trk, summarize_sampling_rate, time_unit = "hour")) |> 
  dplyr::select(id, sr) |> 
  unnest(cols = sr) |>
  left_join(distinct(dplyr::select(tracking_data, id, species))) |>
  arrange(species, median)

# Get the individual sampling rates instead of the summary, for plotting.
trackSummarySamples <- track |>
  mutate(sr = lapply(trk, summarize_sampling_rate, time_unit = "hour", summarize = FALSE)) |>
  dplyr::select(id, sr) |>
  unnest(cols = sr) |>
  mutate(species = ifelse(grepl("BOB", id), "Bobcat", "Coyote"))

boxplot(sr ~ species, outline = FALSE, data = trackSummarySamples, 
        xlab = "", 
        ylab = "Sampling interval in hours, with outliers removed")

# Split into species ----------------------------------------------------------
coyote <- filter(track, grepl("COY", id))
bobcat <- filter(track, grepl("BOB", id))

# Step generation for SSF -----------------------------------------------------
coyote1 <- coyote |> 
  mutate(stp = map(trk, function(df)
    df |> 
      track_resample(rate = hours(4), tolerance = minutes(10)) |> 
      steps_by_burst() |> 
      random_steps(n_control = 10) %>% 
      mutate(log_sl_ = log(sl_ + 1), cos_ta_ = cos(ta_)))) |> 
  dplyr::select(-data, -trk) |> 
  unnest(cols = stp) |> 
  mutate(case_binary_ = ifelse(case_ == TRUE, 1, 0))

bobcat1 <- bobcat |> 
  mutate(stp = map(trk, function(df)
    df |> 
      track_resample(rate = hours(8), tolerance = minutes(10)) |> 
      steps_by_burst() |>  
      random_steps(n_control = 10) |> 
      mutate(log_sl_ = log(sl_ + 1), cos_ta_ = cos(ta_)))) |> 
  dplyr::select(-data, -trk) |> 
  unnest(cols = stp) |> 
  mutate(case_binary_ = ifelse(case_ == TRUE, 1, 0))

# Export resampled step data
write_csv(coyote1, "data/coyote_resampled.csv")
write_csv(bobcat1, "data/bobcat_resampled.csv")

# Read and recalculate step length (in meters) --------------------------------
recalc_steps <- function(file) {
  read_delim(file) |>
    dplyr::select(-log_sl_, -sl_) |>
    mutate(
      sl_ = distGeo(across(c(x1_, y1_)), across(c(x2_, y2_))),
      log_sl_ = log(sl_)
    )
}

coyote_resampled <- recalc_steps("data/coyote_resampled.csv")
bobcat_resampled <- recalc_steps("data/bobcat_resampled.csv")

# Load and prepare rasters ----------------------------------------------------
hfp <- rast("data/HFP_washington.tif")
NAflag(hfp) <- 64536  # Set no-data value
hfp_capped <- classify(hfp, matrix(c(50000, Inf, 50000), ncol = 3, byrow = TRUE)) # Cap at 50k
hfp_scaled <- app(hfp_capped, \(x) (x / 50000) - 0.5)  # Scale to -0.5–0.5

land_use <- rast("data/ESA_washington.tif")

# Land use class labels (ESA WorldCover 2021)
esa_labels <- c(
  "10" = "Tree cover", "20" = "Shrubland", "30" = "Grassland",
  "40" = "Cropland", "50" = "Built-up", "60" = "Bare or sparse vegetation",
  "70" = "Snow and ice", "80" = "Permanent water bodies",
  "90" = "Herbaceous wetland", "95" = "Mangroves", "100" = "Moss and lichen"
)

# Extract covariates from rasters ---------------------------------------------
extract_covariates <- function(df) {
  df |>
    mutate(
      human_footprint = terra::extract(hfp_scaled, cbind(x2_, y2_))[, 1],
      land_use_code = terra::extract(land_use, cbind(x2_, y2_))[, 1],
      land_use = factor(land_use_code, levels = names(esa_labels), labels = esa_labels)
    )
}

coyote_cov <- extract_covariates(coyote_resampled)
bobcat_cov <- extract_covariates(bobcat_resampled)

# Final formatting for SSF ----------------------------------------------------
prepare_ssf_data <- function(df) {
  df |>
    mutate(
      land_use = as.factor(land_use),
      land_use_grouped = fct_collapse(
        land_use,
        "TreeCover" = "Tree cover",
        "Open"      = c("Grassland", "Bare or sparse vegetation", "Moss and lichen"),
        "Cropland"  = "Cropland",
        "BuiltUp"   = "Built-up",
        "Water"     = c("Permanent water bodies", "Herbaceous wetland")
      ),
      step_id_ = paste(id, step_id_, sep = "_") # Unique stratum ID
    ) |>
    group_by(id) |>
    mutate(n = n() / 11) |> # 1 used + 10 available steps per stratum
    ungroup()
}

coyote_final <- prepare_ssf_data(coyote_cov)
bobcat_final <- prepare_ssf_data(bobcat_cov)

# Save final data
write_csv(coyote_final, "data/coyote_ssf_data.csv")
write_csv(bobcat_final, "data/bobcat_ssf_data.csv")

# Fit SSF: Coyote -------------------------------------------------------------
# Read SSF ready data
coyote_ssf_data <- read_delim("data/coyote_ssf_data.csv")

nt <- parallel::detectCores() - 2

ssf_coyote_struc <- glmmTMB(
  case_binary_ ~ -1 + 
    land_use_grouped * human_footprint + 
    offset(log_sl_) + 
    (1 | step_id_) + 
    (0 + land_use_grouped + human_footprint || id),
  family = poisson,
  data = coyote_ssf_data,
  control = glmmTMBControl(parallel = nt),
  doFit = FALSE
)

# Check theta values (to build map/start correctly)
length(ssf_coyote_struc$parameters$theta)

ssf_coyote_struc$parameters$theta[1] <- log(1e3)
ssf_coyote_struc$mapArg <- list(theta = factor(c(NA, 1:6)))
ssf_coyote <- glmmTMB:::fitTMB(ssf_coyote_struc)

#saveRDS(ssf_coyote, file = "models/ssf_coyote_model.rds")
#ssf_coyote <- readRDS("models/ssf_coyote_model.rds")

summary(ssf_coyote)
confint(ssf_coyote)

# Marginal effects of HFP by land cover
trends_coyote <- emtrends(ssf_coyote, ~ land_use_grouped, var = "human_footprint") |>
  summary(infer = c(TRUE, TRUE))

ggplot(trends_coyote, aes(x = land_use_grouped, y = human_footprint.trend)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Effect of Human Footprint") +
  xlab("Land Cover Type") +
  theme_minimal()

# Fit SSF: Bobcat -------------------------------------------------------------
# Read SSF ready data
bobcat_ssf_data <- read_delim("data/bobcat_ssf_data.csv")

# Dropping land use groups with too few used steps
# (virtually no variation in case status)
bobcat_ssf_filtered <- bobcat_ssf_data |> 
  filter(!(land_use_grouped %in% c("BuiltUp", "Water", "Cropland")))

nt <- parallel::detectCores() - 2

ssf_bobcat_struc <- glmmTMB(
  case_binary_ ~ -1 + 
    land_use_grouped * human_footprint + 
    offset(log_sl_) + 
    (1 | step_id_) + 
    (0 + land_use_grouped + human_footprint || id),
  family = poisson,
  data = bobcat_ssf_filtered,
  control = glmmTMBControl(parallel = nt),
  doFit = FALSE
)

# Check theta values (to build map/start correctly)
length(ssf_bobcat_struc$parameters$theta)

ssf_bobcat_struc$parameters$theta[1] <- log(1e3)
ssf_bobcat_struc$mapArg <- list(theta = factor(c(NA, 1:3)))
ssf_bobcat <- glmmTMB:::fitTMB(ssf_bobcat_struc)

#saveRDS(ssf_bobcat, file = "models/ssf_bobcat_model.rds")
#ssf_bobcat <- readRDS("models/ssf_coyote_model.rds")

summary(ssf_bobcat)
confint(ssf_bobcat)

# Marginal trends for bobcats
trends_bobcat <- emtrends(ssf_bobcat, ~ land_use_grouped, var = "human_footprint") |> 
  summary(infer = c(TRUE, TRUE))
