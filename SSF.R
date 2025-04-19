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
library(emmeans)

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

# Get the individual sampling rates instead of the summary, for plotting.
trackSummarySamples <- track |>
  mutate(sr = lapply(trk, summarize_sampling_rate, time_unit = "hour", summarize = FALSE)) |>
  dplyr::select(animal_id, sr) |>
  unnest(cols = c(sr)) |>
  mutate(species = ifelse(grepl("BOB", animal_id), "Bobcat", "Coyote"))
png("img/bobcat_coyote_sampling_rates.png", width = 800, height = 600)
boxplot(sr ~ species, outline = FALSE, data = trackSummarySamples, 
        xlab = "", 
        ylab = "Sampling interval in hours, with outliers removed")
dev.off()

# separate by species
coyote <- track %>% filter(grepl("COY", animal_id))
bobcat <- track %>% filter(grepl("BOB", animal_id))

coyote1 <- coyote %>%
  mutate(stp = map(trk, function(x)
    x %>% track_resample(rate = hours(4), tolerance = minutes(30)) %>%
      steps_by_burst() %>% 
      random_steps(n_control = 10) %>% 
      mutate(log_sl_ = log(sl_ + 1), cos_ta_ = cos(ta_)))) %>%  # Distribute 10 available points for each used location
  dplyr::select(-data, -trk) %>%
  unnest(cols = stp) %>%
  mutate(case_binary = ifelse(case_ == TRUE, 1, 0))

bobcat1 <- bobcat %>%
  mutate(stp = map(trk, function(x)
    x %>% track_resample(rate = hours(8), tolerance = minutes(30)) %>%
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
  "70" = "Snow and ice", "80" = "Permanent water bodies",
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

# Estimate the trends of human_footprint at different levels of land_use.
emtrends_result <- emtrends(coyote_ssf, ~ land_use, var = "human_footprint")
summary(emtrends_result)
plot(emtrends_result)

m2 <- coyote_final |> fit_clogit(case_ ~ human_footprint * land_use + strata(step_id_))
summary(m2)

#Call:
#  coxph(formula = Surv(rep(1, 666248L), case_) ~ human_footprint * 
#          land_use + strata(step_id_), data = data, method = "exact")
#
#n= 666248, number of events= 60568 

#                                                             coef           exp(coef)           se(coef)         z               Pr(>|z|)    
#human_footprint                                 3.116917755441647 22.576685425241724  0.130174030522241  23.94424 < 0.000000000000000222 ***
#land_useShrubland                                              NA                 NA  0.000000000000000        NA                     NA    
#land_useGrassland                              -0.495144762869329  0.609482670449573  0.053296513761235  -9.29038 < 0.000000000000000222 ***
#land_useCropland                               -1.602709580074205  0.201350203687871  0.126634261247812 -12.65621 < 0.000000000000000222 ***
#land_useBuilt-up                               -2.373815685806689  0.093124712847245  0.356767104268213  -6.65368      0.000000000028585 ***
#land_useBare/sparse vegetation                 -1.980960061528149  0.137936745958032  0.430800570678768  -4.59832      0.000004259058602 ***
#land_useSnow and ice                                           NA                 NA  0.000000000000000        NA                     NA    
#land_usePermanent water bodies                 -3.335856990469511  0.035584077924936  0.863941805008322  -3.86121             0.00011283 ***
#land_useHerbaceous wetland                                     NA                 NA  0.000000000000000        NA                     NA    
#land_useMangroves                                              NA                 NA  0.000000000000000        NA                     NA    
#land_useMoss and lichen                         0.671037302717489  1.956265508059556  2.744092165224100   0.24454             0.80681347    
#human_footprint:land_useShrubland                              NA                 NA  0.000000000000000        NA                     NA    
#human_footprint:land_useGrassland              -2.071207096254910  0.126033555215017  0.152963788015728 -13.54051 < 0.000000000000000222 ***
#human_footprint:land_useCropland               -6.000093388697127  0.002478520700039  0.469612663013891 -12.77669 < 0.000000000000000222 ***
#human_footprint:land_useBuilt-up               -6.604150121642969  0.001354734043679  1.242636630606389  -5.31463      0.000000106876053 ***
#human_footprint:land_useBare/sparse vegetation -4.158991940479760  0.015623299200072  1.311947071426306  -3.17009             0.00152391 ** 
#human_footprint:land_useSnow and ice                           NA                 NA  0.000000000000000        NA                     NA    
#human_footprint:land_usePermanent water bodies -7.141538211917665  0.000791533613127  2.845122742607975  -2.51010             0.01206975 *  
#human_footprint:land_useHerbaceous wetland                     NA                 NA  0.000000000000000        NA                     NA    
#human_footprint:land_useMangroves                              NA                 NA  0.000000000000000        NA                     NA    
#human_footprint:land_useMoss and lichen         1.839356554980204  6.292488088012208  6.392162143003827   0.28775             0.77353669    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Concordance= 0.539  (se = 0.001 )
#Likelihood ratio test= 1317.36  on 13 df,   p=< 0.000000000000000222045
#Wald test            = 1248.99  on 13 df,   p=< 0.000000000000000222045
#Score (logrank) test = 1272.25  on 13 df,   p=< 0.000000000000000222045

m3 <- bobcat_final |> fit_clogit(case_ ~ human_footprint * land_use + strata(step_id_))
summary(m3)

#Call:
#  coxph(formula = Surv(rep(1, 118206L), case_) ~ human_footprint * 
#          land_use + strata(step_id_), data = data, method = "exact")
#
#n= 118206, number of events= 10746 
#
#                                                                     coef                     exp(coef)                       se(coef)          z               Pr(>|z|)
#human_footprint                                   -0.0927934692695459167888     0.9113817091901329492387     0.2595776977215661918130    -0.35748               0.720734    
#land_useShrubland                                  0.0000000000000000000000     1.0000000000000000000000     0.0000000000000000000000         NaN                    NaN    
#land_useGrassland                                 -1.1343788010138180499808     0.3216218504596430327602     0.0313001356818177270402   -36.24198 < 0.000000000000000222 ***
#land_useCropland                                  -2.6955932992397100278481     0.0675023208157129045448     0.4151745357833080407595    -6.49267   0.000000000084325515 ***
#land_useBuilt-up                                  -9.3068582406613700186426     0.0000907993672801700994     1.0163126044938173286880    -9.15748 < 0.000000000000000222 ***
#land_useBare/sparse vegetation                    -0.9896770013548591427011     0.3716967291748484725211     0.1560035497221849698501    -6.34394   0.000000000223963554 ***
#land_useSnow and ice                               0.0000000000000000000000     1.0000000000000000000000     0.0000000000000000000000         NaN                    NaN    
#land_usePermanent water bodies                 -4088.0113792708016262622550     0.0000000000000000000000     1.2276406412183498151336 -3329.97397 < 0.000000000000000222 ***
#land_useHerbaceous wetland                         0.0000000000000000000000     1.0000000000000000000000     0.0000000000000000000000         NaN                    NaN    
#land_useMangroves                                  0.0000000000000000000000     1.0000000000000000000000     0.0000000000000000000000         NaN                    NaN    
#land_useMoss and lichen                            1.3645151147404965819021     3.9138248358502174539808     0.1042798555154512002430    13.08513 < 0.000000000000000222 ***
#human_footprint:land_useShrubland                  0.0000000000000000000000     1.0000000000000000000000     0.0000000000000000000000         NaN                    NaN    
#human_footprint:land_useGrassland                 -1.7990737280473276804571     0.1654520708790014060785     0.0776170970004177868118   -23.17883 < 0.000000000000000222 ***
#human_footprint:land_useCropland                  -3.9062198171585351680335     0.0201164011877258130934     1.6261773191233279689527    -2.40209               0.016302 *  
#human_footprint:land_useBuilt-up                 -23.5086809789049340224665     0.0000000000617034680146     3.0648141867883307298825    -7.67051   0.000000000000017132 ***
#human_footprint:land_useBare/sparse vegetation     0.2475950127030740188783     1.2809410622749177743884     0.3288371761102296120249     0.75294               0.451485    
#human_footprint:land_useSnow and ice               0.0000000000000000000000     1.0000000000000000000000     0.0000000000000000000000         NaN                    NaN    
#human_footprint:land_usePermanent water bodies -8335.7938215816575393546373     0.0000000000000000000000     2.5025300504979064442068 -3330.94654 < 0.000000000000000222 ***
#human_footprint:land_useHerbaceous wetland         0.0000000000000000000000     1.0000000000000000000000     0.0000000000000000000000         NaN                    NaN    
#human_footprint:land_useMangroves                  0.0000000000000000000000     1.0000000000000000000000     0.0000000000000000000000         NaN                    NaN    
#human_footprint:land_useMoss and lichen            3.8211372448783156308139    45.6561009749259554268974     0.2196364665097263746407    17.39755 < 0.000000000000000222 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Concordance= 0.53  (se = 0.003 )
#Likelihood ratio test= 358.45  on 21 df,   p=< 0.000000000000000222045
#Wald test            = 22186487.71  on 21 df,   p=< 0.000000000000000222045
#Score (logrank) test = 307.05  on 21 df,   p=< 0.000000000000000222045