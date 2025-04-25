# Load necessary packages -----------------------------------------------------
pacman::p_load(
  tidyverse,
  amt,
  sf,
  geosphere,
  terra,
  glmmTMB,
  DHARMa,
  emmeans,
  gratia,
  ggridges
)

options(scipen = 999) # turn off scientific notation
options(digits = 15) # set digits to 15 to ensure GPS coordinates aren't truncated

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
  dplyr::filter(!(id == "MVBOB71M" & timestamp > as.POSIXct("2019-09-24 00:00:00")))

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

png("img/bobcat_coyote_sampling_rates.png")
boxplot(sr ~ species, outline = FALSE, data = trackSummarySamples, 
        xlab = "", 
        ylab = "Sampling interval in hours")
dev.off()

png("img/bobcat_coyote_sampling_rates_outliers.png")
boxplot(sr ~ species, outline = TRUE, data = trackSummarySamples, 
        xlab = "", 
        ylab = "Sampling interval in hours")
dev.off()

# Split into species ----------------------------------------------------------
coyote <- filter(track, grepl("COY", id))
bobcat <- filter(track, grepl("BOB", id))

# Step generation for SSF -----------------------------------------------------
coyote1 <- coyote[-c(8, 19), ] |> # omitting coyote in row 8 and 19; too few consecutive data points - causing function to fail
  mutate(stp = map(trk, function(df)
    df |> 
      track_resample(rate = hours(4), tolerance = minutes(10)) |> 
      steps_by_burst() |> 
      random_steps(n_control = 10) %>% 
      mutate(log_sl_ = log(sl_ + 1), cos_ta_ = cos(ta_)))) |> 
  dplyr::select(-data, -trk) |> 
  unnest(cols = stp) |> 
  mutate(case_binary_ = ifelse(case_ == TRUE, 1, 0))

bobcat1 <- bobcat[-c(15, 18), ] |> # omitting bobcat in row 15 and 18; too few consecutive data points - causing function to fail
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


# Settings for modeling ---------------------------------------------------

nt <- parallel::detectCores() - 2
options(scipen = 0)
options(digits = 7)

# Fit SSF: Coyote -------------------------------------------------------------
# Read SSF ready data
coyote_ssf_data <- read_delim("data/coyote_ssf_data.csv") |> 
  filter(n > 100) |>  # select animals with more than 100 fixes
  mutate(land_use_grouped = as.factor(land_use_grouped))

coyote_ssf_data$hfp_std <- scale(coyote_ssf_data$human_footprint)[, 1]

# fit SSF with glmmTMB following Muff et al (2019)
ssf_coyote <- glmmTMB(
  case_binary_ ~ -1 + 
    land_use_grouped * (hfp_std + I(hfp_std^2)) + 
    log_sl_ + 
    (0 + land_use_grouped + hfp_std + I(hfp_std^2) + log_sl_|| id) +
    (1 | step_id_),
  family = poisson,
  doFit = TRUE,
  data = coyote_ssf_data,
  map = list(theta = factor(c(1:8, NA))),
  start = list(theta = c(rep(0, times = 8),log(1e3))),
  control = glmmTMBControl(parallel = nt)
)

#saveRDS(ssf_coyote, file = "models/ssf_coyote_model.rds")
#ssf_coyote <- readRDS("models/ssf_coyote_model.rds")

summary(ssf_coyote)

coyote_ssf_pred <- coyote_ssf_data |> 
  filter(case_binary_ == 0) |> 
  mutate(id = NA)  # remove ID for population-level prediction

coy_pred <- predict(ssf_coyote, coyote_ssf_pred, re.form = NA, se.fit = TRUE)
coyote_ssf_pred$fit <- coy_pred$fit
coyote_ssf_pred$se <- coy_pred$se
coyote_ssf_pred <- coyote_ssf_pred |> ungroup()

#saveRDS(coyote_ssf_pred, "data/coyote_ssf_pred.rds")
#coyote_ssf_pred <- readRDS("data/coyote_ssf_pred.rds")

ggplot(coyote_ssf_pred, aes(x = exp(human_footprint), y = fit, color = land_use_grouped)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(ymin = fit - se, ymax = fit + se, fill = land_use_grouped), stat = "identity", alpha = 0.2) +
  theme_minimal()

# average effect plot function --------------------------------------------

avg_eff_plot_hfp_landuse <- function(fittedResponse, nsim = 100, k = 10, showPeakValue = TRUE) {
  
  set.seed(123)
  # Monte Carlo sampling of predictions based on fit ± se
  fit_sample <- sapply(rep(1, nsim), function(i) {
    rnorm(n = nrow(fittedResponse), mean = fittedResponse$fit, sd = fittedResponse$se)
  })
  
  fit_sample_matrix <- matrix(fit_sample, nrow = nrow(fittedResponse), ncol = nsim)
  
  smooth_list <- list()
  
  for (j in 1:nsim) {
    message("GAM smoothing sample ", j)
    smooth_list[[j]] <- mgcv::bam(
      fit_sample_matrix[, j] ~ s(human_footprint, by = land_use_grouped, bs = "ts", k = k) + land_use_grouped,
      data = fittedResponse,
      select = TRUE, discrete = TRUE, nthreads = nt
    ) |>
      gratia::smooth_estimates(unconditional = FALSE, overall_uncertainty = TRUE) |>
      gratia::add_confint() |>
      dplyr::rename(hfp = human_footprint)
  }
  
  avg_smooth <- dplyr::bind_rows(smooth_list) |>
    dplyr::group_by(.smooth, .by, land_use_grouped, hfp) |>
    dplyr::summarise(
      est = mean(.estimate),
      upper_ci = mean(.upper_ci),
      lower_ci = mean(.lower_ci),
      .groups = "drop"
    )
  
  # Faceted plot instead of overlaid lines
  p <- ggplot(avg_smooth, aes(x = hfp, y = est)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.25, fill = "steelblue", colour = NA) +
    geom_line(color = "steelblue", size = 0.9) +
    facet_wrap(~land_use_grouped, scales = "free_y") +
    theme_minimal() +
    labs(x = "Human Footprint", y = "Relative Use (log scale)") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
  
  if (showPeakValue) {
    peak_vals <- avg_smooth |>
      group_by(land_use_grouped) |>
      filter(est == max(est)) |>
      slice(rep(1:n(), each = 2)) |>
      mutate(est = ifelse(row_number() %% 2 == 1, est, -Inf))
    
    p <- p + geom_line(data = peak_vals, aes(x = hfp, y = est, group = land_use_grouped), 
                       linetype = "dashed", size = 0.6, alpha = 0.3, color = "black")
  }
  
  print(p)
  return(list(plot = p))
}


coyote_ssf_pred<- coyote_ssf_pred |> 
  mutate(land_use_grouped = as.factor(land_use_grouped))


avg_eff_output <- avg_eff_plot_hfp_landuse(coyote_ssf_pred, nsim = 100, showPeakValue = TRUE)

avg_eff_output[["plot"]][["plot_env"]][["peak_vals"]]



# Relative Selection Strength ---------------------------------------------

calc_rss_hfp_landcover <- function(model, data, land_use_col = "land_use_grouped", hfp_col = "hfp_std", n_points = 100) {
  
  levels <- unique(data[[land_use_col]])
  all_rss <- purrr::map_dfr(levels, function(lc) {
    
    # Subset data for this land use class
    dat_lc <- data |> filter(!!sym(land_use_col) == lc)
    
    # Range of human footprint within this class
    hfp_range <- range(dat_lc[[hfp_col]], na.rm = TRUE)
    hfp_seq <- seq(from = hfp_range[1], to = hfp_range[2], length.out = n_points)
    
    # Build newdata for prediction
    newdata <- expand.grid(
      hfp_std = hfp_seq,
      land_use_grouped = lc
    ) |> 
      mutate(
        `I(hfp_std^2)` = hfp_std^2,
        log_sl_ = mean(data$log_sl_, na.rm = TRUE),
        step_id_ = NA,
        id = NA
      )
    
    # Create a baseline: same data with lowest hfp
    baseline <- newdata |> filter(hfp_std == min(hfp_std))
    
    # Predict
    x1_pred <- predict(model, newdata = newdata, re.form = NA)
    x2_pred <- predict(model, newdata = baseline[rep(1, nrow(newdata)), ], re.form = NA)
    
    # Calculate logRSS
    rss_df <- newdata |>
      mutate(
        logRSS = x1_pred - x2_pred,
        RSS = exp(logRSS),
        land_use_grouped = lc
      )
    return(rss_df)
  })
  
  return(all_rss)
}

rss_data <- calc_rss_hfp_landcover(model = ssf_coyote, data = coyote_ssf_data)

ggplot(rss_data, aes(x = hfp_std, y = RSS, color = land_use_grouped)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  labs(x = "Human Footprint", y = "Relative Selection Strength (RSS)", color = "Land Cover") +
  theme_minimal()

# Some plot exploration ---------------------------------------------------
ggplot(coyote_ssf_data |> filter(case_binary_ == 1),
       aes(x = human_footprint, fill = land_use_grouped)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(x = "Human Footprint", y = "Density",
       title = "Distribution of Used HFP by Land Cover")


# Bivariate density hex plot (continuous y-variable vs human footprint)
hex_plot_fun <- function(data, yvar, ylab, title_) {
  ggplot(data, aes(x = human_footprint, y = .data[[yvar]])) +
    geom_bin2d(aes(fill = after_stat(log(density)))) +
    scale_fill_viridis_c(limits = c(-14, -1), option = "D") +
    stat_density_2d(
      aes(x = human_footprint, y = .data[[yvar]]),
      colour = "red", breaks = c(0.05), n = 15, size = 1
    ) +
    geom_vline(xintercept = 0, colour = "gray70", size = 1) +
    geom_hline(yintercept = 0, colour = "gray70", size = 1) +
    labs(
      x = "Human Footprint",
      y = ylab,
      title = title_
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}

hex_plot_fun(coyote_ssf_data, yvar = "log_sl_", ylab = "Log Step Length", title_ = "Coyote")
hex_plot_fun(coyote_ssf_data, yvar = "cos_ta_", ylab = "Cosine Turning Angle", title_ = "Coyote")

ggplot(coyote_ssf_data, aes(x = land_use_grouped, y = human_footprint)) +
  geom_boxplot(outlier.alpha = 0.1) +
  geom_jitter(width = 0.2, alpha = 0.05) +
  theme_minimal() +
  xlab("Land Use Type") +
  ylab("Human Footprint") +
  ggtitle("Human Footprint across Land Use Types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(coyote_ssf_data, aes(x = land_use_grouped, y = human_footprint)) +
  geom_violin(fill = "skyblue", alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  theme_minimal() +
  xlab("Land Use Type") +
  ylab("Human Footprint") +
  ggtitle("Distribution of HFP by Land Use")

ggplot(coyote_ssf_data, aes(x = human_footprint, y = land_use_grouped)) +
  geom_density_ridges(scale = 1.2, fill = "lightgreen", alpha = 0.7) +
  theme_minimal() +
  xlab("Human Footprint") +
  ylab("Land Use") +
  ggtitle("Distribution of Human Footprint across Land Use")


ggplot(coyote_ssf_data, aes(x = land_use_grouped, y = human_footprint, fill = factor(case_binary_))) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values = c("0" = "grey80", "1" = "dodgerblue"), name = "Used?") +
  theme_minimal() +
  xlab("Land Use Type") +
  ylab("Human Footprint") +
  ggtitle("HFP by Land Use (Used vs. Available)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Fit SSF: Bobcat -------------------------------------------------------------
# Read SSF ready data
bobcat_ssf_data <- read_delim("data/bobcat_ssf_data.csv") |> 
  filter(n > 100) # select animals with more than 100 fixes
table(bobcat_ssf_data$land_use_grouped)
# Dropping land use groups with too few used steps
# (virtually no variation in case status)
bobcat_ssf_filtered <- bobcat_ssf_data |> 
  filter(!(land_use_grouped %in% c("BuiltUp", "Snow and ice", "Water", "Cropland")))

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
