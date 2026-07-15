library(tidyverse)
library(lubridate)
library(lme4)
library(car)
library(DHARMa)
library(MuMIn)
library(corrplot)
library(broom.mixed)
library(patchwork)
library(brms)
library(loo)
library(posterior)
library(bayesplot)
library(lmerTest)
library(overlap)
library(ggeffects)

#===============================================================================
# 1. DATA PROCESSING (Using your existing code structure)
#===============================================================================

# Load in files
cams <- read.csv('cropton_cams_behaviour.csv')
meta <- read_csv('crop_meta.csv')
locations <- read.csv('cam_locations.csv')
eDNA <- read.csv('crop_prc.csv', row.names = 1)
distances <- read.csv('distance_matrix.csv')

# Process meta data
meta_processed <- meta %>%
  mutate(
    sampling_datetime = dmy_hm(sampling_date),
    sampling_date_only = as.Date(sampling_datetime),
    camera = paste0("CROP", str_pad(str_sub(sample, 8, 9), width = 2, pad = "0")),
    month = month(sampling_date_only),
  )

# Process eDNA data
eDNA_processed <- eDNA %>%
  mutate(sample = row.names(eDNA)) %>%
  pivot_longer(cols = c("Castor_fiber", "Lutra_lutra"), names_to = "species", values_to = "prop_reads") %>%
  mutate(
    presence = as.integer(prop_reads > 0),
    camera_id_part = str_sub(sample, 8, 9),
    camera = paste0("CROP", str_pad(camera_id_part, width = 2, pad = "0")),
    month = as.integer(str_sub(sample, 6, 7))
  ) %>%
  filter(species %in% c("Castor_fiber", "Lutra_lutra")) %>%
  left_join(meta_processed %>% select(sample, sampling_datetime, sampling_date_only, temp, flow, pH, cond, vol), 
            by = "sample") %>%
  select(-camera_id_part)

# Merge with locations
eDNA_data <- eDNA_processed %>%
  left_join(locations %>% select(camera, lodge_distance), by = "camera") %>%
  mutate(
    lodge_distance = as.numeric(lodge_distance),
    lodge_distance = replace_na(lodge_distance, round(mean(as.numeric(locations$lodge_distance), na.rm = TRUE))),
    has_camera = !camera %in% paste0("CROP", str_pad(11:14, width = 2, pad = "0"))
  ) %>%
  filter(complete.cases(sampling_datetime))

# Process camera data
cams_clean <- cams %>%
  mutate(
    datetime = lubridate::dmy_hm(date),
    date_only = as.Date(datetime),
    camera = stringr::str_extract(camera, "CROP\\d+")
  ) %>%
  filter(species %in% c("Castor_fiber", "Lutra_lutra"), !is.na(datetime))

# Filter faulty deployments
dodge_cams <- c("CROP_1001", "CROP_1103", "CROP_0108", "CROP_0208", "CROP_1303", 
                "CROP_0901", "CROP_0903", "CROP_0908", "CROP_0909","CROP_0910", "CROP_0801")

# Create dataset for visits (30-day window)
cams_filtered <- cams_clean %>%
  inner_join(meta_processed %>% select(sample, camera, sampling_datetime), by = c("camera", "sample")) %>%
  mutate(days_before_sampling = as.numeric(difftime(sampling_datetime, datetime, units = "days"))) %>%
  filter(days_before_sampling >= 0, days_before_sampling <= 30, !sample %in% dodge_cams)

# Calculate visits per species
visits_per_species <- cams_filtered %>%
  arrange(camera, sample, species, datetime) %>%
  group_by(camera, sample, sampling_datetime, species) %>%
  mutate(
    time_diff = as.numeric(difftime(datetime, lag(datetime), units = "mins")),
    visit_id = cumsum(is.na(time_diff) | time_diff > 30)
  ) %>%
  summarise(n_visits = n_distinct(visit_id), .groups = "drop")

# Calculate days since last visit for all samples
calculate_days_since_last_visit <- function() {
  sampling_events <- eDNA_data %>%
    filter(species %in% c("Castor_fiber", "Lutra_lutra")) %>%
    distinct(sample, species, camera, sampling_datetime) %>%
    arrange(camera, sample, species)
  
  days_since_results <- tibble()
  
  for(i in 1:nrow(sampling_events)) {
    current_camera <- sampling_events$camera[i]
    current_sample <- sampling_events$sample[i]
    current_datetime <- sampling_events$sampling_datetime[i]
    current_species <- sampling_events$species[i]
    
    if (current_sample %in% dodge_cams) {
      days_result <- NA_real_
    } else {
      prior_detections <- cams_clean %>%
        filter(
          camera == current_camera,
          species == current_species,
          datetime < current_datetime
        ) %>%
        mutate(days_before = as.numeric(difftime(current_datetime, datetime, units = "days"))) %>%
        filter(days_before <= 60)
      
      if (nrow(prior_detections) > 0) {
        days_result <- min(prior_detections$days_before)
      } else {
        days_result <- 60
      }
    }
    
    result_row <- tibble(
      camera = current_camera,
      sample = current_sample,
      sampling_datetime = current_datetime,
      species = current_species,
      days_since_last_visit = days_result
    )
    
    days_since_results <- bind_rows(days_since_results, result_row)
  }
  
  return(days_since_results)
}

# Calculate all days since last visit
days_since_visit_data <- calculate_days_since_last_visit()

# Calculate maximum individuals per day
max0 <- function(x) if (length(x) == 0 || all(is.na(x))) 0 else max(x, na.rm = TRUE)

daily_max <- cams_filtered %>%
  group_by(camera, sample, sampling_datetime, date_only, species) %>%
  summarise(Nijt = max0(count), .groups = "drop") %>%
  group_by(camera, sample, sampling_datetime, species) %>%
  summarise(max_ind_day = max(Nijt, na.rm = TRUE), .groups = "drop")

# Calculate upstream detection
upstream_detection <- eDNA_data %>%
  select(sample, species, camera, sampling_date_only, presence, lodge_distance) %>%
  group_by(species, sampling_date_only) %>%
  arrange(lodge_distance) %>%
  mutate(upstream_presence = lag(presence, default = 0)) %>%
  ungroup() %>%
  select(sample, species, camera, sampling_date_only, upstream_presence)

# Calculate behaviour proportions AND absolute counts AND diversity measures
behaviour_metrics <- cams_filtered %>%
  mutate(
    behaviour = str_trim(behaviour),
    behaviour_group = case_when(
      behaviour %in% c("SW", "WA", "DV", "RU", "EW", "EX") ~ "locomotion",
      behaviour %in% c("FO", "EA", "SC", "LC", "HA", "HB", "EP") ~ "feeding",
      behaviour %in% c("PL", "SG", "AL", "KC", "NT", "VO") ~ "social",
      behaviour %in% c("SM", "FI", "FE", "AT", "SP") ~ "defence",
      TRUE ~ "other"
    )
  ) %>%
  group_by(camera, sample, species) %>%
  summarise(
    total_behaviours = n(),
    n_locomotion = sum(behaviour_group == "locomotion"),
    n_feeding = sum(behaviour_group == "feeding"),
    n_social = sum(behaviour_group == "social"),
    n_defence = sum(behaviour_group == "defence"),
    
    # Within-group diversity: diversity of specific behaviors within each group
    locomotion_diversity = {
      loco_behaviors <- behaviour[behaviour_group == "locomotion"]
      if(length(loco_behaviors) == 0) {
        0
      } else {
        loco_counts <- table(loco_behaviors)
        loco_props <- loco_counts / sum(loco_counts)
        -sum(loco_props * log(loco_props))
      }
    },
    
    feeding_diversity = {
      feed_behaviors <- behaviour[behaviour_group == "feeding"]
      if(length(feed_behaviors) == 0) {
        0
      } else {
        feed_counts <- table(feed_behaviors)
        feed_props <- feed_counts / sum(feed_counts)
        -sum(feed_props * log(feed_props))
      }
    },
    
    social_diversity = {
      social_behaviors <- behaviour[behaviour_group == "social"]
      if(length(social_behaviors) == 0) {
        0
      } else {
        social_counts <- table(social_behaviors)
        social_props <- social_counts / sum(social_counts)
        -sum(social_props * log(social_props))
      }
    },
    
    defence_diversity = {
      defence_behaviors <- behaviour[behaviour_group == "defence"]
      if(length(defence_behaviors) == 0) {
        0
      } else {
        defence_counts <- table(defence_behaviors)
        defence_props <- defence_counts / sum(defence_counts)
        -sum(defence_props * log(defence_props))
      }
    },
    
    .groups = "drop"
  ) %>%
  # Calculate proportions first
  mutate(
    prop_locomotion = ifelse(total_behaviours == 0, 0, n_locomotion / total_behaviours),
    prop_feeding = ifelse(total_behaviours == 0, 0, n_feeding / total_behaviours),
    prop_social = ifelse(total_behaviours == 0, 0, n_social / total_behaviours),
    prop_defence = ifelse(total_behaviours == 0, 0, n_defence / total_behaviours)
  ) %>%
  # Calculate hierarchical diversity using rowwise for proper evaluation
  rowwise() %>%
  mutate(
    # Level 1: Diversity across the 4 main behavioral groups
    group_diversity = {
      if(total_behaviours == 0) {
        0
      } else {
        H_groups <- 0
        if(prop_locomotion > 0) H_groups <- H_groups - (prop_locomotion * log(prop_locomotion))
        if(prop_feeding > 0) H_groups <- H_groups - (prop_feeding * log(prop_feeding))
        if(prop_social > 0) H_groups <- H_groups - (prop_social * log(prop_social))
        if(prop_defence > 0) H_groups <- H_groups - (prop_defence * log(prop_defence))
        H_groups
      }
    },
    
    # Level 2: Average within-group diversity
    within_group_diversity = {
      if(total_behaviours == 0) {
        0
      } else {
        n_active_groups <- (n_locomotion > 0) + (n_feeding > 0) + (n_social > 0) + (n_defence > 0)
        if(n_active_groups == 0) {
          0
        } else {
          total_within_diversity <- 0
          if(n_locomotion > 0) total_within_diversity <- total_within_diversity + locomotion_diversity
          if(n_feeding > 0) total_within_diversity <- total_within_diversity + feeding_diversity
          if(n_social > 0) total_within_diversity <- total_within_diversity + social_diversity
          if(n_defence > 0) total_within_diversity <- total_within_diversity + defence_diversity
          total_within_diversity / n_active_groups
        }
      }
    },
    
    # Combined hierarchical diversity
    hierarchical_diversity = group_diversity + within_group_diversity,
    
    # Simple behavior richness (number of different behavior groups used)
    behavior_richness = (n_locomotion > 0) + (n_feeding > 0) + (n_social > 0) + (n_defence > 0)
    
  ) %>%
  ungroup() %>%
  left_join(meta_processed %>% select(sample, camera, sampling_datetime), 
            by = c("sample", "camera"))

#===============================================================================
# 2. CALCULATE DETECTABILITY SCORES 
#===============================================================================

# Parameters for detectability calculation
spatial_param <- 200  # Half-weight distance in meters
temporal_param <- 14  # Half-weight time in days
max_days <- 60       # Maximum lookback period
max_distance <- 1000 # Maximum upstream distance

# Prepare locations data and distance matrix
locations_clean <- locations %>%
  select(camera, lodge_distance) %>%
  mutate(lodge_distance = as.numeric(lodge_distance)) %>%
  filter(!is.na(lodge_distance))

# Process distance matrix - handle triangular structure with upstream distances
distance_matrix_clean <- distances %>%
  column_to_rownames(var = "X")  # Use X column as row names

# Function to get upstream distance from detection_camera to sampling_camera
get_upstream_distance <- function(detection_camera, sampling_camera) {
  if(detection_camera %in% rownames(distance_matrix_clean) && 
     sampling_camera %in% colnames(distance_matrix_clean)) {
    
    distance <- distance_matrix_clean[detection_camera, sampling_camera]
    
    # If NA, try the reverse (since it's triangular matrix)
    if(is.na(distance) && sampling_camera %in% rownames(distance_matrix_clean) && 
       detection_camera %in% colnames(distance_matrix_clean)) {
      reverse_distance <- distance_matrix_clean[sampling_camera, detection_camera]
      # If reverse distance exists and is positive, detection is downstream (return NA)
      # If reverse distance is negative, detection is upstream (return absolute value)
      if(!is.na(reverse_distance)) {
        if(reverse_distance < 0) {
          return(abs(reverse_distance))  # Detection is upstream
        } else {
          return(NA_real_)  # Detection is downstream
        }
      }
    }
    
    # Return distance if positive (upstream), NA if negative (downstream)
    if(!is.na(distance)) {
      return(ifelse(distance >= 0, distance, NA_real_))
    }
  }
  
  return(NA_real_)
}

# Get unique eDNA samples to process
edna_samples <- eDNA_data %>%
  distinct(sample, species, camera, sampling_datetime, presence, lodge_distance) %>%
  filter(!is.na(sampling_datetime))

cat("Processing", nrow(edna_samples), "eDNA samples for detectability scores...\n")

# Calculate detectability scores - manual loop
detectability_results <- tibble()

for(i in 1:nrow(edna_samples)) {
  
  if(i %% 20 == 0) cat("Processing sample", i, "of", nrow(edna_samples), "\n")
  
  current_sample <- edna_samples[i, ]
  current_camera <- current_sample$camera
  current_species <- current_sample$species
  current_datetime <- current_sample$sampling_datetime
  current_lodge_distance <- current_sample$lodge_distance
  current_sample_id <- current_sample$sample
  
  # Check if current sampling site is faulty
  if (current_sample_id %in% dodge_cams) {
    result_row <- tibble(
      sample = current_sample_id,
      species = current_species,
      camera = current_camera,
      presence = current_sample$presence,
      detectability_score = NA_real_,
      n_contributing_detections = NA_real_
    )
    detectability_results <- bind_rows(detectability_results, result_row)
    next
  }
  
  # Find all camera detections for this species before sampling
  relevant_detections <- cams_clean %>%
    filter(
      species == current_species,
      datetime < current_datetime,
      as.numeric(difftime(current_datetime, datetime, units = "days")) <= max_days
    ) %>%
    mutate(
      days_since = as.numeric(difftime(current_datetime, datetime, units = "days")),
      # Get upstream distance from detection camera to current sampling camera
      upstream_distance = map_dbl(camera, ~get_upstream_distance(.x, current_camera))
    ) %>%
    # Only keep detections that are upstream (positive distance) and within max distance
    filter(!is.na(upstream_distance), upstream_distance <= max_distance)
  
  if(nrow(relevant_detections) == 0) {
    # No relevant detections found
    result_row <- tibble(
      sample = current_sample_id,
      species = current_species,
      camera = current_camera,
      presence = current_sample$presence,
      detectability_score = 0,
      n_contributing_detections = 0
    )
    detectability_results <- bind_rows(detectability_results, result_row)
    next
  }
  
  # Calculate squared inverse weights for each detection
  weighted_detections <- relevant_detections %>%
    mutate(
      spatial_weight = 1 / (1 + (upstream_distance / spatial_param)^2),
      temporal_weight = 1 / (1 + (days_since / temporal_param)^2),
      combined_weight = spatial_weight * temporal_weight,
      count_to_use = ifelse(is.na(count) | count == 0, 1, count),
      weighted_contribution = count_to_use * combined_weight
    )
  
  # Sum all weighted contributions
  total_score <- sum(weighted_detections$weighted_contribution, na.rm = TRUE)
  n_detections <- nrow(weighted_detections)
  
  result_row <- tibble(
    sample = current_sample_id,
    species = current_species,
    camera = current_camera,
    presence = current_sample$presence,
    detectability_score = total_score,
    n_contributing_detections = n_detections
  )
  
  detectability_results <- bind_rows(detectability_results, result_row)
}

cat("Completed detectability score calculation!\n")

#===============================================================================
# 3. CREATE MASTER DATASET
#===============================================================================

# Create the master dataset with ALL potential metrics
model_data_raw <- eDNA_data %>%
  left_join(visits_per_species, by = c("sample", "species", "camera", "sampling_datetime")) %>%
  left_join(days_since_visit_data, by = c("sample", "species", "camera", "sampling_datetime")) %>%
  left_join(daily_max, by = c("sample", "species", "camera", "sampling_datetime")) %>%
  left_join(behaviour_metrics, by = c("sample", "species", "camera", "sampling_datetime")) %>%
  left_join(detectability_results %>% select(sample, species, detectability_score, n_contributing_detections),
            by = c("sample", "species")) %>%
  left_join(upstream_detection, by = c("sample", "species", "camera", "sampling_date_only")) %>%
  mutate(
    # Handle missing data for cameras without deployments
    is_crop11 = str_detect(camera, "CROP11"),
    
    # n_visits handling
    n_visits = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(n_visits, 0)
    ),
    
    # days_since_last_visit handling
    days_since_last_visit = case_when(
      !has_camera & species == "Castor_fiber" ~ NA_real_,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ days_since_last_visit
    ),
    
    # max_ind_day handling
    max_ind_day = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(max_ind_day, 0)
    ),
    
    # ABSOLUTE BEHAVIOR COUNTS
    total_behaviours = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(total_behaviours, 0)
    ),
    n_locomotion = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(n_locomotion, 0)
    ),
    n_feeding = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(n_feeding, 0)
    ),
    n_social = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(n_social, 0)
    ),
    n_defence = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(n_defence, 0)
    ),
    
    # Behavior proportions handling
    prop_locomotion = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(prop_locomotion, 0)
    ),
    prop_feeding = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(prop_feeding, 0)
    ),
    prop_social = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(prop_social, 0)
    ),
    prop_defence = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(prop_defence, 0)
    ),
    
    # DIVERSITY MEASURES
    locomotion_diversity = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(locomotion_diversity, 0)
    ),
    feeding_diversity = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(feeding_diversity, 0)
    ),
    social_diversity = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(social_diversity, 0)
    ),
    defence_diversity = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(defence_diversity, 0)
    ),
    group_diversity = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(group_diversity, 0)
    ),
    within_group_diversity = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(within_group_diversity, 0)
    ),
    hierarchical_diversity = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(hierarchical_diversity, 0)
    ),
    behavior_richness = case_when(
      !has_camera & species == "Castor_fiber" ~ 0,
      !has_camera & species == "Lutra_lutra" ~ NA_real_,
      is_crop11 & species == "Castor_fiber" ~ NA_real_,
      sample %in% dodge_cams ~ NA_real_,
      TRUE ~ replace_na(behavior_richness, 0)
    ),
    
    upstream_presence = as.factor(upstream_presence)
  ) %>%
  # Remove the helper variable
  select(-is_crop11)

# Scale variables - ALL potential metrics for correlation analysis
model_data <- model_data_raw %>%
  mutate(
    # Environmental variables
    temp_z = as.numeric(scale(temp)),
    flow_z = as.numeric(scale(log(flow + 1))),
    pH_z = as.numeric(scale(pH)),
    cond_z = as.numeric(scale(sqrt(pmax(cond - min(cond, na.rm = TRUE) + 1, 0)))),
    vol_z = as.numeric(scale(log(vol + 1))),
    lodge_distance_z = as.numeric(scale(lodge_distance)),
    
    # Camera metrics
    n_visits_z = as.numeric(scale(log(n_visits + 1))),
    detectability_score_z = as.numeric(scale(log(detectability_score + 0.001))),
    days_since_last_visit_z = as.numeric(scale(log(days_since_last_visit + 1))),
    max_ind_day_z = as.numeric(scale(log(max_ind_day + 1))),
    
    # Behavior proportions (arcsine transformation)
    prop_locomotion_z = as.numeric(scale(asin(sqrt(prop_locomotion)))),
    prop_feeding_z = as.numeric(scale(asin(sqrt(prop_feeding)))),
    prop_social_z = as.numeric(scale(asin(sqrt(prop_social)))),
    prop_defence_z = as.numeric(scale(asin(sqrt(prop_defence)))),
    
    # Absolute behavior counts (log transformation)
    total_behaviours_z = as.numeric(scale(log(total_behaviours + 1))),
    n_locomotion_z = as.numeric(scale(log(n_locomotion + 1))),
    n_feeding_z = as.numeric(scale(log(n_feeding + 1))),
    n_social_z = as.numeric(scale(log(n_social + 1))),
    n_defence_z = as.numeric(scale(log(n_defence + 1))),
    
    # Diversity measures
    group_diversity_z = as.numeric(scale(group_diversity)),
    within_group_diversity_z = as.numeric(scale(within_group_diversity)),
    hierarchical_diversity_z = as.numeric(scale(hierarchical_diversity)),
    behavior_richness_z = as.numeric(scale(behavior_richness)),
    locomotion_diversity_z = as.numeric(scale(locomotion_diversity)),
    feeding_diversity_z = as.numeric(scale(feeding_diversity)),
    social_diversity_z = as.numeric(scale(social_diversity)),
    defence_diversity_z = as.numeric(scale(defence_diversity)),
    
    # Seasonality variables
    day_of_year = yday(sampling_date_only),
    sin_doy = sin(2 * pi * day_of_year / 365.25),
    cos_doy = cos(2 * pi * day_of_year / 365.25),
    sin_doy_z = as.numeric(scale(sin_doy)),
    cos_doy_z = as.numeric(scale(cos_doy)),
    
    # Month for random effects
    month_factor = as.factor(month),
    site = camera
  )

# Remove first 83 columns
model_data <- model_data[, -(1:92)]

#===============================================================================
# 4. OTTER DATA
#===============================================================================

# Filter to otter data only and remove rows with missing key variables
otter_data <- model_data %>%
  filter(species == "Lutra_lutra") %>%
  filter(complete.cases(presence, temp_z, flow_z, pH_z, vol_z, lodge_distance_z, month_factor))

cat("Otter dataset:", nrow(otter_data), "samples\n")
cat("Otter detection rate:", round(mean(otter_data$presence, na.rm = TRUE) * 100, 1), "%\n")

# Check data availability for both models
otter_visits_complete <- otter_data %>%
  filter(complete.cases(n_visits_z, prop_locomotion_z, prop_feeding_z, prop_social_z, prop_defence_z))

otter_detectability_complete <- otter_data %>%
  filter(complete.cases(detectability_score_z, prop_locomotion_z, prop_feeding_z, prop_social_z, prop_defence_z))

cat("Otter n_visits model complete cases:", nrow(otter_visits_complete), "\n")
cat("Otter detectability model complete cases:", nrow(otter_detectability_complete), "\n")

#===============================================================================
# 5. CORRELATION AND MULTICOLLINEARITY ANALYSIS
#===============================================================================

# First check which variables have zero variance
all_predictors <- otter_data %>%
  select(# Camera metrics
    n_visits_z, detectability_score_z, days_since_last_visit_z, max_ind_day_z,
    # Behavior proportions
    prop_social_z, prop_locomotion_z, prop_feeding_z, prop_defence_z,
    # # Absolute behavior counts
    # total_behaviours_z, n_locomotion_z, n_feeding_z, n_social_z, n_defence_z,
    # Diversity measures
    group_diversity_z, hierarchical_diversity_z, behavior_richness_z,
    # Environmental variables
    temp_z, flow_z, pH, cond_z, vol_z, lodge_distance_z,
    # Seasonality
    sin_doy_z, cos_doy_z) %>%
  filter(complete.cases(.))

all_cor_matrix <- cor(all_predictors, use = "complete.obs")

cat("\nFull Correlation Matrix Plot:\n")
corrplot(all_cor_matrix, method = "color", type = "upper", tl.col = "black", 
         tl.srt = 45, addCoef.col = "black", number.cex = 0.5)

# N_VISITS MODEL - VIF analysis
visits_vif_model <- lm(presence ~ prop_social_z + prop_locomotion_z + prop_feeding_z + prop_defence_z +
                         temp_z + flow_z + pH + vol_z + lodge_distance_z + n_visits_z,
                       data = otter_visits_complete)

visits_vif_results <- vif(visits_vif_model)
print(round(visits_vif_results, 2))

# DETECTABILITY MODEL VIF
detectability_vif_model <- lm(presence ~ prop_social_z + prop_locomotion_z + prop_feeding_z + prop_defence_z +
                                hierarchical_diversity_z + temp_z + flow_z + pH + vol_z + lodge_distance_z + detectability_score_z,
                              data = otter_detectability_complete)

detectability_vif_results <- vif(detectability_vif_model)
print(round(detectability_vif_results, 2))

#===============================================================================
# 7. MODEL 1: OTTER PRESENCE WITH N_VISITS
#===============================================================================

# Full model with n_visits
otter_visits_full <- glmer(
  presence ~ prop_locomotion_z + prop_feeding_z + prop_social_z + prop_defence_z + 
    n_visits_z + temp_z + flow_z + pH + vol_z + lodge_distance_z + (1|month_factor),
  data = otter_visits_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

summary(otter_visits_full)

# Manual stepwise removal - you can modify this as needed
otter_visits_no_ph <- glmer(
  presence ~ prop_locomotion_z + prop_feeding_z + prop_social_z + prop_defence_z + 
    n_visits_z + temp_z + flow_z + vol_z + lodge_distance_z + (1|month_factor),
  data = otter_visits_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

# Compare models
AIC(otter_visits_full, otter_visits_no_ph) # Model improves with pH exclusion
summary(otter_visits_no_ph)

# Remove feeding
otter_visits_no_feeding <- glmer(
  presence ~ prop_locomotion_z + prop_defence_z + prop_social_z + 
    n_visits_z + temp_z + flow_z + vol_z + lodge_distance_z + (1|month_factor),
  data = otter_visits_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

AIC(otter_visits_no_ph, otter_visits_no_feeding)
summary(otter_visits_no_feeding)

# Remove volume
otter_visits_no_vol <- glmer(
  presence ~ prop_locomotion_z + prop_defence_z + prop_social_z + 
    n_visits_z + temp_z + flow_z + lodge_distance_z + (1|month_factor),
  data = otter_visits_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

AIC(otter_visits_no_vol, otter_visits_no_feeding)
summary(otter_visits_no_vol)

# Remove defence from visits model (CRUCIAL FIX)
otter_visits_no_defence <- glmer(
  presence ~ prop_locomotion_z + prop_social_z + 
    n_visits_z + temp_z + flow_z + lodge_distance_z + (1|month_factor),
  data = otter_visits_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

AIC(otter_visits_no_vol, otter_visits_no_defence)
summary(otter_visits_no_defence)

# Test removing lodge distance
otter_visits_no_distance <- glmer(
  presence ~ prop_locomotion_z + prop_social_z + 
    n_visits_z + temp_z + flow_z + (1|month_factor),
  data = otter_visits_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

AIC(otter_visits_no_defence, otter_visits_no_distance) # Check which is better

# Choose the best visits model (without defence)
otter_visits_final <- otter_visits_no_defence

summary(otter_visits_final)
r2_visits <- r.squaredGLMM(otter_visits_final)
cat("R-squared for visits model - Marginal:", round(r2_visits[1], 3), "Conditional:", round(r2_visits[2], 3), "\n")

#===============================================================================
# 8. MODEL 2: OTTER PRESENCE WITH DETECTABILITY SCORE
#===============================================================================

cat("\n=== MODEL 2: OTTER PRESENCE WITH DETECTABILITY SCORE ===\n")

# Full model with detectability score
otter_detectability_full <- glmer(
  presence ~ prop_locomotion_z + prop_feeding_z + prop_social_z + prop_defence_z + hierarchical_diversity +
    detectability_score_z + temp_z + flow_z + pH_z + vol_z + lodge_distance_z + (1|month_factor),
  data = otter_detectability_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

summary(otter_detectability_full)

# Manual stepwise removal
otter_detectability_no_ph <- glmer(
  presence ~ prop_locomotion_z + prop_feeding_z + prop_social_z + prop_defence_z + hierarchical_diversity +
    detectability_score_z + temp_z + flow_z + vol_z + lodge_distance_z + (1|month_factor),
  data = otter_detectability_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

AIC(otter_detectability_full, otter_detectability_no_ph)
summary(otter_detectability_no_ph)

# Remove behavioural diversity
otter_detectability_no_div <- glmer(
  presence ~ prop_locomotion_z + prop_feeding_z + prop_social_z + prop_defence_z +
    detectability_score_z + temp_z + flow_z + vol_z + lodge_distance_z + (1|month_factor),
  data = otter_detectability_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

AIC(otter_detectability_no_div, otter_detectability_no_ph)
summary(otter_detectability_no_div)

# Remove feeding
otter_detectability_no_feeding <- glmer(
  presence ~ prop_locomotion_z + prop_social_z + prop_defence_z +
    detectability_score_z + temp_z + flow_z + vol_z + lodge_distance_z + (1|month_factor),
  data = otter_detectability_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

AIC(otter_detectability_no_div, otter_detectability_no_feeding)
summary(otter_detectability_no_feeding)

# Remove volume
otter_detectability_no_vol <- glmer(
  presence ~ prop_locomotion_z + prop_social_z + prop_defence_z +
    detectability_score_z + temp_z + flow_z + lodge_distance_z + (1|month_factor),
  data = otter_detectability_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

AIC(otter_detectability_no_feeding, otter_detectability_no_vol)
summary(otter_detectability_no_vol)

# Remove defence from detectability model (CRUCIAL FIX)
otter_detectability_no_defence <- glmer(
  presence ~ prop_locomotion_z + prop_social_z +
    detectability_score_z + temp_z + flow_z + lodge_distance_z + (1|month_factor),
  data = otter_detectability_complete,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

AIC(otter_detectability_no_defence, otter_detectability_no_vol)
summary(otter_detectability_no_defence) # This should be better

# Choose the final detectability model (without defence)
otter_detectability_final <- otter_detectability_no_defence

summary(otter_detectability_final)
r2_detectability <- r.squaredGLMM(otter_detectability_final)
cat("R-squared for detectability model - Marginal:", round(r2_detectability[1], 3), "Conditional:", round(r2_detectability[2], 3), "\n")

#===============================================================================
# 8B. OTTER MODEL DIAGNOSTICS
#===============================================================================

# Residual plots for visits model
cat("Visits model diagnostics:\n")
res_visits <- simulateResiduals(fittedModel = otter_visits_final, n = 1000)
plot(res_visits)
testUniformity(res_visits)
testDispersion(res_visits)

# Residual plots for detectability model
cat("\nDetectability model diagnostics:\n")
res_detectability <- simulateResiduals(fittedModel = otter_detectability_final, n = 1000)
plot(res_detectability)
testUniformity(res_detectability)
testDispersion(res_detectability)

#===============================================================================
# 8C. OTTER MODEL COMPARISON - FIXED TO USE CORRECT MODELS
#===============================================================================

# R-squared comparison using the corrected models (both without defence)
cat("R-squared comparison (both models without defence):\n")
cat("Visits model - Marginal:", round(r2_visits[1], 3), "Conditional:", round(r2_visits[2], 3), "\n")
cat("Detectability model - Marginal:", round(r2_detectability[1], 3), "Conditional:", round(r2_detectability[2], 3), "\n")

# Choose best otter model (using corrected models)
if(r2_detectability[2] > r2_visits[2]) {
  cat("→ Detectability model is better for otters\n")
  otter_best_model <- otter_detectability_final
  otter_best_type <- "detectability"
} else {
  cat("→ Visits model is better for otters\n")
  otter_best_model <- otter_visits_final
  otter_best_type <- "visits"
}

# VERIFICATION: Check that defence is removed from final model
cat("\nFINAL OTTER MODEL VERIFICATION:\n")
cat("Selected model type:", otter_best_type, "\n")
cat("Model predictors:", paste(names(fixef(otter_best_model)), collapse = ", "), "\n")

has_defence_final <- any(grepl("defence", names(fixef(otter_best_model))))
if(has_defence_final) {
  cat("❌ ERROR: Defence is still in the final otter model!\n")
} else {
  cat("✅ SUCCESS: Defence has been correctly removed from the final otter model\n")
}

#===============================================================================
# BEAVER DATA
#===============================================================================

# Filter to beaver data only and remove rows with missing key variables
beaver_data <- model_data %>%
  filter(species == "Castor_fiber") %>%
  filter(complete.cases(presence, temp_z, flow_z, vol_z, lodge_distance_z, month_factor))

cat("Beaver dataset:", nrow(beaver_data), "samples\n")
cat("Beaver detection rate:", round(mean(beaver_data$presence, na.rm = TRUE) * 100, 1), "%\n")

# Check data availability for both models
beaver_visits_complete <- beaver_data %>%
  filter(complete.cases(n_visits_z, prop_locomotion_z, prop_social_z, prop_defence_z))

beaver_detectability_complete <- beaver_data %>%
  filter(complete.cases(detectability_score_z, prop_locomotion_z, prop_social_z, prop_defence_z))

cat("Beaver n_visits model complete cases:", nrow(beaver_visits_complete), "\n")
cat("Beaver detectability model complete cases:", nrow(beaver_detectability_complete), "\n")

#===============================================================================
# BEAVER CORRELATION AND MULTICOLLINEARITY ANALYSIS
#===============================================================================

# Create correlation matrix for beaver predictors
beaver_predictors <- beaver_data %>%
  select(# Camera metrics
    n_visits_z, detectability_score_z, days_since_last_visit_z, max_ind_day_z,
    # Behavior proportions
    prop_social_z, prop_locomotion_z, prop_feeding_z, prop_defence_z,
    # Diversity measures
    group_diversity_z, hierarchical_diversity_z, behavior_richness_z,
    # Environmental variables
    temp_z, flow_z, pH, cond_z, vol_z, lodge_distance_z,
    # Seasonality
    sin_doy_z, cos_doy_z) %>%
  filter(complete.cases(.))

# Create correlation matrix
beaver_cor_matrix <- cor(beaver_predictors, use = "complete.obs")

corrplot(beaver_cor_matrix, method = "color", type = "upper", tl.col = "black", 
         tl.srt = 45, addCoef.col = "black", number.cex = 0.5,
         title = "Beaver Predictor Correlations", mar = c(0,0,1,0))

# BEAVER N_VISITS MODEL - VIF analysis
beaver_visits_vif_model <- lm(presence ~ prop_social_z + prop_locomotion_z + prop_defence_z +
                                temp_z + flow_z + vol_z + sin_doy + pH + lodge_distance_z + n_visits_z,
                              data = beaver_visits_complete)

beaver_visits_vif_results <- vif(beaver_visits_vif_model)
print(round(beaver_visits_vif_results, 2))

# BEAVER DETECTABILITY MODEL VIF
beaver_detectability_vif_model <- lm(presence ~ prop_social_z + prop_locomotion_z + prop_defence_z + hierarchical_diversity +
                                       temp_z + flow_z + vol_z + sin_doy + pH + lodge_distance_z + detectability_score_z,
                                     data = beaver_detectability_complete)

beaver_detectability_vif_results <- vif(beaver_detectability_vif_model)
print(round(beaver_detectability_vif_results, 2))

#===============================================================================
# 9. BEAVER PRESENCE ANALYSIS (BAYESIAN GLMM) - WITH VOLUME/DIVERSITY TESTS
#===============================================================================

cat("\n=== BEAVER PRESENCE ANALYSIS (BAYESIAN GLMM) ===\n")

beaver_detection_rate <- mean(beaver_data$presence, na.rm = TRUE)
cat("Beaver detection rate:", round(beaver_detection_rate * 100, 1), "%\n")
cat("Using Bayesian GLMM due to high detection rate\n")

# Prepare beaver datasets (matching your existing structure)
beaver_visits_complete <- beaver_data %>%
  filter(complete.cases(presence, n_visits_z, prop_locomotion_z, 
                        prop_social_z, prop_defence_z,
                        temp_z, flow_z, vol_z, pH, lodge_distance_z, sin_doy))

beaver_detectability_complete <- beaver_data %>%
  filter(complete.cases(presence, detectability_score_z, prop_locomotion_z, 
                        prop_social_z, prop_defence_z, hierarchical_diversity,
                        temp_z, flow_z, vol_z, lodge_distance_z))

cat("Beaver complete cases - Visits:", nrow(beaver_visits_complete), 
    "Detectability:", nrow(beaver_detectability_complete), "\n")

# Set informative priors for high detection rate scenario
beaver_priors <- c(
  set_prior("normal(0, 1.5)", class = "b"),           # Regularizing priors for coefficients
  set_prior("normal(2, 1)", class = "Intercept"),     # Informed prior for high detection
  set_prior("student_t(3, 0, 2.5)", class = "sd")     # Random effects variance
)

#===============================================================================
# BEAVER MODEL 1: N_VISITS (BAYESIAN)
#===============================================================================

cat("\n--- Beaver N_Visits Model (Bayesian) ---\n")

# Full model
beaver_visits_full_brm <- brm(
  presence ~ prop_locomotion_z + prop_social_z + prop_defence_z + prop_feeding_z +
    hierarchical_diversity + n_visits_z + temp_z + flow_z + vol_z + sin_doy + pH + lodge_distance_z +
    (1 | month_factor),
  data = beaver_visits_complete,
  family = bernoulli(link = "logit"),
  prior = beaver_priors,
  chains = 4, iter = 4000, warmup = 2000, cores = 4,
  control = list(adapt_delta = 0.99),  # Better for high detection rates
  save_pars = save_pars(all = TRUE),
  silent = 2, refresh = 0  # Reduce output clutter
)

# Enhanced removal function that tests coefficient stability
test_predictor_removal <- function(model, predictor, model_name = "model") {
  cat("Testing removal of", predictor, "from", model_name, "...\n")
  
  tryCatch({
    # Update model
    updated <- update(model, formula. = as.formula(paste0(". ~ . -", predictor)),
                      silent = 2, refresh = 0)
    
    # Calculate log-likelihoods (more robust than LOO for problematic cases)
    loglik_orig <- sum(colMeans(log_lik(model)))
    loglik_new <- sum(colMeans(log_lik(updated)))
    loglik_diff <- loglik_new - loglik_orig
    
    # Check coefficient stability for key predictors
    orig_summary <- summary(model)$fixed
    new_summary <- summary(updated)$fixed
    
    # Identify key predictors that should remain stable
    key_predictors <- c("detectability_score_z", "n_visits_z", "lodge_distance_z")
    key_predictors <- key_predictors[key_predictors %in% rownames(orig_summary)]
    
    if(length(key_predictors) > 0) {
      max_change <- 0
      for(key_pred in key_predictors) {
        if(key_pred %in% rownames(new_summary)) {
          change <- abs(new_summary[key_pred, "Estimate"] - orig_summary[key_pred, "Estimate"])
          max_change <- max(max_change, change)
        }
      }
      
      cat("   Log-likelihood change:", round(loglik_diff, 2), "\n")
      cat("   Max coefficient change:", round(max_change, 3), "\n")
      
      # Decision criteria: remove if log-likelihood change is small AND coefficients are stable
      if(abs(loglik_diff) < 1.0 && max_change < 0.2) {
        cat("→ Safe to remove", predictor, "(minimal impact on fit and coefficients)\n")
        return(updated)
      } else if(abs(loglik_diff) < 2.0 && max_change < 0.15) {
        cat("→ Acceptable to remove", predictor, "(small impact)\n")
        return(updated)
      } else {
        cat("→ Keep", predictor, "(substantial impact detected)\n")
        return(model)
      }
    } else {
      # Fallback to log-likelihood only
      if(abs(loglik_diff) < 1.5) {
        cat("→ Remove", predictor, "(log-likelihood change:", round(loglik_diff, 2), ")\n")
        return(updated)
      } else {
        cat("→ Keep", predictor, "(log-likelihood change:", round(loglik_diff, 2), ")\n")
        return(model)
      }
    }
    
  }, error = function(e) {
    cat("Error with", predictor, "- keeping original model:", e$message, "\n")
    return(model)
  })
}

# Stepwise removal for visits model - prioritizing vol_z and hierarchical_diversity
cat("\nStepwise selection for beaver visits model:\n")
beaver_visits_step1 <- test_predictor_removal(beaver_visits_full_brm, "vol_z", "visits")
beaver_visits_step2 <- test_predictor_removal(beaver_visits_step1, "hierarchical_diversity", "visits")
beaver_visits_step3 <- test_predictor_removal(beaver_visits_step2, "prop_feeding_z", "visits")
beaver_visits_step4 <- test_predictor_removal(beaver_visits_step3, "sin_doy", "visits")
beaver_visits_step5 <- test_predictor_removal(beaver_visits_step4, "pH", "visits")
beaver_visits_final_brm <- beaver_visits_step5

cat("Final beaver visits model formula:", as.character(formula(beaver_visits_final_brm))[3], "\n")

#===============================================================================
# BEAVER MODEL 2: DETECTABILITY (BAYESIAN)
#===============================================================================

cat("\n--- Beaver Detectability Model (Bayesian) ---\n")

# Full model
beaver_detectability_full_brm <- brm(
  presence ~ prop_locomotion_z + prop_social_z + prop_defence_z + 
    hierarchical_diversity + detectability_score_z + temp_z + flow_z + vol_z + lodge_distance_z +
    (1 | month_factor),
  data = beaver_detectability_complete,
  family = bernoulli(link = "logit"),
  prior = beaver_priors,
  chains = 4, iter = 4000, warmup = 2000, cores = 4,
  control = list(adapt_delta = 0.99),
  save_pars = save_pars(all = TRUE),
  silent = 2, refresh = 0
)

# Stepwise removal for detectability model - prioritizing vol_z and hierarchical_diversity
cat("\nStepwise selection for beaver detectability model:\n")
beaver_detect_step1 <- test_predictor_removal(beaver_detectability_full_brm, "vol_z", "detectability")
beaver_detect_step2 <- test_predictor_removal(beaver_detect_step1, "hierarchical_diversity", "detectability")
beaver_detect_step3 <- test_predictor_removal(beaver_detect_step2, "prop_defence_z", "detectability")
beaver_detectability_final_brm <- beaver_detect_step3

cat("Final beaver detectability model formula:", as.character(formula(beaver_detectability_final_brm))[3], "\n")

#===============================================================================
# BEAVER MODEL COMPARISON (BAYESIAN) - FIXED FOR DIFFERENT SAMPLE SIZES
#===============================================================================

cat("\n=== BEAVER MODEL COMPARISON (BAYESIAN) ===\n")

# Check sample sizes
cat("Sample sizes:\n")
cat("Visits model:", nrow(beaver_visits_complete), "observations\n")
cat("Detectability model:", nrow(beaver_detectability_complete), "observations\n")

# Create common dataset for fair comparison (without vol_z since it should be removed)
beaver_common_complete <- beaver_data %>%
  filter(complete.cases(presence, n_visits_z, detectability_score_z, prop_locomotion_z, 
                        prop_social_z, prop_defence_z, temp_z, flow_z, 
                        lodge_distance_z, month_factor))

cat("Common complete cases:", nrow(beaver_common_complete), "observations\n")

if(nrow(beaver_common_complete) < 50) {
  cat("⚠️  Small sample size for comparison. Using separate model evaluation.\n")
  
  # Use log-likelihood for comparison (most robust approach)
  visits_loglik <- log_lik(beaver_visits_final_brm)
  detect_loglik <- log_lik(beaver_detectability_final_brm)
  
  visits_total_loglik <- sum(colMeans(visits_loglik))
  detect_total_loglik <- sum(colMeans(detect_loglik))
  
  cat("Visits model total log-likelihood:", round(visits_total_loglik, 2), "\n")
  cat("Detectability model total log-likelihood:", round(detect_total_loglik, 2), "\n")
  
  # Decide based on total log-likelihood
  if(detect_total_loglik > visits_total_loglik) {
    cat("→ Detectability model is better for beavers\n")
    beaver_best_model_brm <- beaver_detectability_final_brm
    beaver_best_type <- "detectability"
  } else {
    cat("→ Visits model is better for beavers\n")
    beaver_best_model_brm <- beaver_visits_final_brm
    beaver_best_type <- "visits"
  }
  
} else {
  cat("✅ Sufficient common data for LOO comparison. Refitting models.\n")
  
  # Refit both models on common dataset for fair comparison
  beaver_visits_common_brm <- update(beaver_visits_final_brm, 
                                     newdata = beaver_common_complete,
                                     silent = 2, refresh = 0)
  
  beaver_detect_common_brm <- update(beaver_detectability_final_brm, 
                                     newdata = beaver_common_complete,
                                     silent = 2, refresh = 0)
  
  # Now compare with LOO
  visits_loo_common <- loo(beaver_visits_common_brm, moment_match = TRUE)
  detect_loo_common <- loo(beaver_detect_common_brm, moment_match = TRUE)
  
  final_comparison <- loo_compare(visits_loo_common, detect_loo_common)
  
  cat("Beaver model comparison (LOO on common dataset):\n")
  print(final_comparison)
  
  # Determine best model
  best_model_name <- rownames(final_comparison)[1]
  if(grepl("detect", best_model_name)) {
    cat("→ Detectability model is better for beavers\n")
    beaver_best_model_brm <- beaver_detectability_final_brm  # Use original model
    beaver_best_type <- "detectability"
  } else {
    cat("→ Visits model is better for beavers\n")
    beaver_best_model_brm <- beaver_visits_final_brm  # Use original model
    beaver_best_type <- "visits"
  }
  
  cat("Best beaver model ELPD difference:", round(final_comparison[1, "elpd_diff"], 2), "±", round(final_comparison[1, "se_diff"], 2), "\n")
}

cat("\nFinal selection: Using", beaver_best_type, "model for beaver analysis\n")

#===============================================================================
# VERIFY MODEL STRUCTURE CONSISTENCY
#===============================================================================

cat("\n=== VERIFYING MODEL STRUCTURE CONSISTENCY ===\n")

# Extract predictor names (excluding intercept)
beaver_predictors <- names(fixef(beaver_best_model_brm))[-1]
otter_predictors <- names(fixef(otter_best_model))[-1]

cat("Beaver final predictors:", paste(beaver_predictors, collapse = ", "), "\n")
cat("Otter final predictors:", paste(otter_predictors, collapse = ", "), "\n")

# Check if structures match (accounting for different naming)
beaver_clean <- gsub("b_", "", beaver_predictors)
otter_clean <- otter_predictors

common_predictors <- intersect(beaver_clean, otter_clean)
cat("Common predictors:", paste(common_predictors, collapse = ", "), "\n")

if(length(common_predictors) >= 5) {
  cat("✅ Models have consistent structure for meaningful comparison\n")
} else {
  cat("⚠️  Models have different structures - comparison may be limited\n")
}

#===============================================================================
# FINAL MODEL SUMMARIES
#===============================================================================

cat("\n=== FINAL UNIFIED MODEL SUMMARIES ===\n")

cat("Beaver final model (", beaver_best_type, "):\n")
print(summary(beaver_best_model_brm))

cat("\nOtter final model (", otter_best_type, "):\n")
print(summary(otter_best_model))

#===============================================================================
# PART E: FOREST PLOT COMPARISON (2x2 LAYOUT) ----
#===============================================================================

cat("\n=== CREATING 2x2 FOREST PLOT COMPARISON ===\n")

# Function to standardize predictor names
standardize_predictor_names <- function(term) {
  case_when(
    str_detect(term, "Intercept") ~ "Intercept",
    str_detect(term, "prop_locomotion_z") ~ "Locomotion Behavior",
    str_detect(term, "prop_social_z") ~ "Social Behavior",
    str_detect(term, "prop_defence_z") ~ "Defence Behavior",
    str_detect(term, "n_visits_z") ~ "Number of Visits",
    str_detect(term, "detectability_score_z") ~ "Detectability Score",
    str_detect(term, "temp_z") ~ "Temperature",
    str_detect(term, "flow_z") ~ "Flow Rate",
    str_detect(term, "vol_z") ~ "Sample Volume",
    str_detect(term, "lodge_distance_z") ~ "Lodge Distance",
    TRUE ~ term
  )
}

# Extract coefficients for best models
extract_model_coefs_updated <- function(model, model_type, species, analysis_type) {
  if(model_type == "bayesian") {
    # Extract Bayesian results
    posterior_summary <- posterior_summary(model)
    fixed_effects <- posterior_summary[grepl("^b_", rownames(posterior_summary)), ]
    
    coef_df <- data.frame(
      term = gsub("^b_", "", rownames(fixed_effects)),
      estimate = fixed_effects[, "Estimate"],
      se = fixed_effects[, "Est.Error"],
      ci_lower = fixed_effects[, "Q2.5"],
      ci_upper = fixed_effects[, "Q97.5"],
      stringsAsFactors = FALSE
    )
    
    # Calculate probability of direction (convert to p-value-like)
    draws <- as_draws_df(model)
    fixed_cols <- grep("^b_", names(draws), value = TRUE)
    prob_direction <- sapply(fixed_cols, function(col) {
      samples <- draws[[col]]
      max(mean(samples > 0), mean(samples < 0))
    })
    
    coef_df$p_value <- 1 - prob_direction[paste0("b_", coef_df$term)]
    
  } else if(model_type == "glmer") {
    # Extract frequentist results - FIX: properly rename p.value to p_value
    coef_df <- tidy(model, conf.int = TRUE) %>%
      filter(effect == "fixed") %>%
      select(term, estimate, std.error, p.value, conf.low, conf.high) %>%
      rename(se = std.error, ci_lower = conf.low, ci_upper = conf.high, p_value = p.value)  # THIS IS THE FIX!
  }
  
  coef_df$species <- species
  coef_df$analysis_type <- analysis_type
  coef_df$model_type <- model_type
  
  return(coef_df)
}

# Extract data for best models with updated function
beaver_presence_data <- extract_model_coefs_updated(beaver_best_model_brm, "bayesian", "Beaver", "Presence")
otter_presence_data <- extract_model_coefs_updated(otter_best_model, "glmer", "Otter", "Presence")

# Rest of your forest plot code stays the same...
plot_data_presence <- bind_rows(beaver_presence_data, otter_presence_data) %>%
  filter(!str_detect(term, "Intercept")) %>%
  mutate(
    predictor = standardize_predictor_names(term),
    significant = p_value < 0.05
  )

create_forest_plot <- function(data, title, color, species_name) {
  # Filter out intercept and prepare data
  plot_data <- data %>%
    filter(!str_detect(term, "Intercept")) %>%
    mutate(
      # Calculate odds ratios
      OR = exp(estimate),
      OR_lower = exp(ci_lower),
      OR_upper = exp(ci_upper),
      # Create significance indicator
      significant = p_value < 0.05,
      # Clean predictor names if not already done
      predictor = if_else(is.na(predictor), standardize_predictor_names(term), predictor)
    ) %>%
    # Reorder by effect size for better visualization
    arrange(desc(abs(estimate)))
  
  # Create the plot
  ggplot(plot_data, aes(x = estimate, y = reorder(predictor, estimate))) +
    # Add vertical reference line at 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
    
    # Add confidence intervals
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), 
                   height = 0.2, color = color, alpha = 0.7, size = 1) +
    
    # Add points for estimates
    geom_point(aes(fill = significant), 
               size = 3, shape = 21, color = color, stroke = 1.2) +
    
    # Customize fill for significance
    scale_fill_manual(values = c("TRUE" = color, "FALSE" = "white"),
                      name = "Significant\n(p < 0.05)",
                      labels = c("FALSE" = "No", "TRUE" = "Yes")) +
    
    # Labels and theme
    labs(
      title = title,
      x = "Effect Size (Log Odds)",
      y = ""
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 9),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", size = 0.5),
      panel.grid.major.x = element_line(color = "gray95", size = 0.3)
    )
}

# Create the plots
beaver_presence_plot <- create_forest_plot(
  filter(plot_data_presence, species == "Beaver"),
  paste("Beaver Presence", paste0("(", beaver_best_type, ")")),
  "#8B4513",
  "Beaver"
)

otter_presence_plot <- create_forest_plot(
  filter(plot_data_presence, species == "Otter"),
  paste("Otter Presence", paste0("(", otter_best_type, ")")),
  "#2E8B57",
  "Otter"
)

# Display the comparison
presence_comparison <- beaver_presence_plot + otter_presence_plot +
  plot_annotation(title = "Species Comparison: Final Models (Defence Removed from Both)")


#===============================================================================
# OTTER PREY MODELLING ----
#===============================================================================

otter_prey <- read.csv('otter_prey_eDNA.csv')

# Remove non-prey species and calculate richness
otter_prey_clean <- otter_prey %>%
  mutate(
    sample = otter_prey$X,  # Use the X column which should contain sample IDs
    # Extract site and month from the actual sample ID
    site = paste0("CROP", str_pad(str_sub(sample, 8, 9), width = 2, pad = "0")),
    month = as.integer(str_sub(sample, 6, 7))
  ) %>%
  select(-X, -Castor_fiber, -Lutra_lutra) %>%  
  # Calculate Shannon diversity
  rowwise() %>%
  mutate(
    # Get total reads for this sample
    total_reads = sum(c_across(Cobitis_taenia:Anas), na.rm = TRUE),
    # Calculate Shannon diversity
    prey_diversity = if_else(total_reads == 0, 0, {
      props <- c_across(Cobitis_taenia:Anas) / total_reads
      props <- props[props > 0]  # Remove zeros
      -sum(props * log(props), na.rm = TRUE)
    })
  ) %>%
  ungroup() %>%
  select(sample, site, month, prey_diversity)

# Monthly otter behavior from cameras
monthly_otter_behavior <- cams_filtered %>%
  filter(species == "Lutra_lutra") %>%
  mutate(
    month = month(date_only),
    behaviour_group = case_when(
      behaviour %in% c("SW", "WA", "DV", "RU", "EW", "EX") ~ "locomotion",
      behaviour %in% c("FO", "EA", "SC", "LC", "HA", "HB", "EP") ~ "feeding",
      behaviour %in% c("PL", "SG", "AL", "KC", "NT", "VO") ~ "social",
      behaviour %in% c("SM", "FI", "FE", "AT", "SP") ~ "defence",
      TRUE ~ "other"
    ),
    feeding_behavior = behaviour_group == "feeding"
  ) %>%
  arrange(camera, datetime) %>%
  group_by(camera, month) %>%
  mutate(
    time_diff = as.numeric(difftime(datetime, lag(datetime), units = "mins")),
    visit_id = cumsum(is.na(time_diff) | time_diff > 30)
  ) %>%
  summarise(
    total_detections = n(),
    n_visits = n_distinct(visit_id),
    n_feeding_detections = sum(feeding_behavior, na.rm = TRUE),
    prop_feeding_detections = n_feeding_detections / total_detections,
    feeding_intensity = n_feeding_detections / n_visits,  # Feeding events per visit
    .groups = "drop"
  ) %>%
  mutate(site = paste0("CROP", str_pad(str_extract(camera, "\\d+"), width = 2, pad = "0"))) %>%
  select(-camera)

# Combine with prey data
otter_prey_behavior <- monthly_otter_behavior %>%
  left_join(otter_prey_clean, by = c("site", "month")) %>%
  filter(!is.na(prey_diversity), total_detections > 0)

# Add environmental variables
env_data <- model_data_raw %>%  # Use the raw version before column removal
  filter(species == "Lutra_lutra") %>%
  select(sample, camera, month, temp, flow, pH, vol) %>%
  mutate(
    site = paste0("CROP", str_pad(str_extract(camera, "\\d+"), width = 2, pad = "0")),
    # Create scaled versions
    temp_z = as.numeric(scale(temp)),
    flow_z = as.numeric(scale(log(flow + 1))),
    pH_z = as.numeric(scale(pH)),
    vol_z = as.numeric(scale(log(vol + 1)))
  ) %>%
  select(site, month, temp_z, flow_z, pH, vol_z)

# Join with your prey-behavior data
otter_prey_behavior_env <- otter_prey_behavior %>%
  left_join(env_data, by = c("site", "month")) %>%
  filter(!is.na(temp_z))  # Remove rows without environmental data

### Models ----

hist(otter_prey_behavior$n_visits)
hist(otter_prey_behavior$prop_feeding_detections)


activity_model <- glmer(
  n_visits ~ prey_diversity + (1|site) + (1|month),
  data = otter_prey_behavior,
  family = poisson  # Count data
)

summary(activity_model)

feeding_prop_model <- glmer(
  prop_feeding_detections ~ prey_diversity + (1|site),
  data = otter_prey_behavior,
  family = binomial,
  weights = total_detections  # Weight by sample size
)

summary(feeding_prop_model)


activity_model_env <- glmer(
  n_visits ~ prey_diversity + flow_z + pH + (1|site) + (1|month),
  data = otter_prey_behavior_env,
  family = poisson
)
summary(activity_model_env)

feeding_model_env <- glmer(
  prop_feeding_detections ~ prey_diversity + pH + (1|site),
  data = otter_prey_behavior_env,
  family = binomial,
  weights = total_detections
)

summary(feeding_model_env)

# Compare the activity models
AIC(activity_model)
AIC(activity_model_env)

# Compare the feeding models
AIC(feeding_prop_model)
AIC(feeding_model_env)

r.squaredGLMM(activity_model)
r.squaredGLMM(activity_model_env)

r.squaredGLMM(feeding_prop_model)
r.squaredGLMM(feeding_model_env)


library(ggplot2)
library(ggeffects)
library(patchwork)

# ------------------------------------
# 1. Plotting the Activity Model
# ------------------------------------
# Use the best-fitting model (with environmental variables)
# Generate predictions for the relationship between prey_diversity and n_visits,
# while holding other variables (temp_z, flow_z, pH) at their mean values.
plot_data_activity_env <- ggpredict(activity_model_env, terms = "prey_diversity", type = "fixed")

# Create the plot for otter activity
plot_activity_env <- ggplot() +
  # Add the raw data points
  geom_point(data = otter_prey_behavior_env, aes(x = prey_diversity, y = n_visits), alpha = 0.5) +
  # Add the model's predicted line and confidence interval
  geom_line(data = plot_data_activity_env, aes(x = x, y = predicted), size = 1.2, color = "darkblue") +
  geom_ribbon(data = plot_data_activity_env, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightblue") +
  labs(
    title = "Otter Site Use vs. Prey Diversity",
    x = "Prey Diversity (Shannon Index)",
    y = "Predicted Number of Otter Visits"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# ------------------------------------
# 2. Plotting the Feeding Model
# ------------------------------------
# Use the best-fitting model (with environmental variables)
# Generate predictions for the relationship between prey_diversity and
# prop_feeding_detections, holding other variables at their means.
plot_data_feeding_env <- ggpredict(feeding_model_env, terms = "prey_diversity", type = "fixed")

# Create the plot for feeding proportion
plot_feeding_env <- ggplot() +
  # Add the raw data points
  geom_point(data = otter_prey_behavior_env, aes(x = prey_diversity, y = prop_feeding_detections), alpha = 0.5) +
  # Add the model's predicted line and confidence interval
  geom_line(data = plot_data_feeding_env, aes(x = x, y = predicted), size = 1.2, color = "darkgreen") +
  geom_ribbon(data = plot_data_feeding_env, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "lightgreen") +
  labs(
    title = "Feeding Proportion vs. Prey Diversity",
    x = "Prey Diversity (Shannon Index)",
    y = "Predicted Proportion of Feeding Detections"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# ------------------------------------
# 3. Combine the Plots
# ------------------------------------
# Use the patchwork package to arrange the two plots side-by-side with labels.
combined_plot <- (plot_activity_env + plot_feeding_env) +
  plot_annotation(tag_levels = 'A')

# Display the final combined plot
combined_plot

#===============================================================================
# STREAMLINED TEMPORAL AVOIDANCE ANALYSIS WITH RANDOMIZATION TEST
#===============================================================================

# Data preparation (keep your original)
activity_data <- cams_clean %>%
  filter(species %in% c("Castor_fiber", "Lutra_lutra")) %>%
  mutate(
    season = case_when(
      month(datetime) %in% c(12,1,2) ~ "winter",
      month(datetime) %in% c(3,4,5) ~ "spring", 
      month(datetime) %in% c(6,7,8) ~ "summer",
      month(datetime) %in% c(9,10,11) ~ "autumn"
    ),
    site = camera
  )

# Calculate succession intervals
calculate_succession_intervals <- function(data, time_window_days = 7) {
  beaver_data <- data %>% filter(species == "Castor_fiber") %>% arrange(site, datetime)
  otter_data <- data %>% filter(species == "Lutra_lutra") %>% arrange(site, datetime)
  
  succession_intervals <- tibble()
  
  for(i in 1:nrow(otter_data)) {
    current_otter <- otter_data[i, ]
    
    recent_beaver <- beaver_data %>%
      filter(
        site == current_otter$site,
        datetime < current_otter$datetime,
        datetime >= (current_otter$datetime - days(time_window_days))
      ) %>%
      slice_max(datetime, n = 1)
    
    if(nrow(recent_beaver) > 0) {
      hours_between <- as.numeric(difftime(current_otter$datetime, 
                                           recent_beaver$datetime, 
                                           units = "hours"))
      
      succession_intervals <- bind_rows(succession_intervals,
                                        tibble(
                                          site = current_otter$site,
                                          hours_between = hours_between,
                                          otter_time = current_otter$datetime,
                                          season = current_otter$season
                                        ))
    }
  }
  return(succession_intervals)
}

# Simple randomization function
randomize_otter_times <- function(data) {
  otter_data <- data %>% filter(species == "Lutra_lutra")
  beaver_data <- data %>% filter(species == "Castor_fiber")
  
  # Completely shuffle otter times
  otter_shuffled <- otter_data %>%
    mutate(datetime = sample(datetime, size = n(), replace = FALSE))
  
  # Combine and return
  bind_rows(beaver_data, otter_shuffled) %>%
    arrange(site, datetime)
}

# Main randomization test
run_temporal_avoidance_test <- function(data, n_permutations = 1000, 
                                        time_windows = c(3, 7, 14)) {
  
  results_list <- list()
  
  for(window in time_windows) {
    cat("\n=== TESTING", window, "DAY WINDOW ===\n")
    
    # Calculate observed statistics
    observed_data <- calculate_succession_intervals(data, time_window_days = window)
    
    if(nrow(observed_data) == 0) {
      cat("No succession events found for", window, "day window\n")
      next
    }
    
    observed_stats <- observed_data %>%
      summarise(
        n_events = n(),
        median_hours = median(hours_between, na.rm = TRUE),
        mean_hours = mean(hours_between, na.rm = TRUE),
        sd_hours = sd(hours_between, na.rm = TRUE),
        prop_within_24h = mean(hours_between < 24, na.rm = TRUE),
        prop_within_6h = mean(hours_between < 6, na.rm = TRUE),
        prop_within_1h = mean(hours_between < 1, na.rm = TRUE)
      )
    
    cat("Observed succession events:", observed_stats$n_events, "\n")
    cat("Observed median hours:", round(observed_stats$median_hours, 1), "\n")
    
    # Run permutations
    cat("Running", n_permutations, "permutations...\n")
    
    null_stats <- tibble()
    
    for(i in 1:n_permutations) {
      if(i %% 100 == 0) cat("  ", i, "/", n_permutations, "\n")
      
      # Randomize and calculate
      randomized_data <- randomize_otter_times(data)
      random_succession <- calculate_succession_intervals(randomized_data, time_window_days = window)
      
      if(nrow(random_succession) > 0) {
        random_stats <- random_succession %>%
          summarise(
            median_hours = median(hours_between, na.rm = TRUE),
            mean_hours = mean(hours_between, na.rm = TRUE),
            prop_within_24h = mean(hours_between < 24, na.rm = TRUE),
            prop_within_6h = mean(hours_between < 6, na.rm = TRUE),
            prop_within_1h = mean(hours_between < 1, na.rm = TRUE),
            permutation = i
          )
        
        null_stats <- bind_rows(null_stats, random_stats)
      }
    }
    
    # Calculate p-values and effect sizes
    p_values <- tibble(
      metric = c("median_hours", "mean_hours", "prop_within_24h", "prop_within_6h", "prop_within_1h"),
      observed = c(observed_stats$median_hours, observed_stats$mean_hours, 
                   observed_stats$prop_within_24h, observed_stats$prop_within_6h, 
                   observed_stats$prop_within_1h),
      expected_null = c(mean(null_stats$median_hours, na.rm = TRUE),
                        mean(null_stats$mean_hours, na.rm = TRUE),
                        mean(null_stats$prop_within_24h, na.rm = TRUE),
                        mean(null_stats$prop_within_6h, na.rm = TRUE),
                        mean(null_stats$prop_within_1h, na.rm = TRUE)),
      # Two-sided p-values
      p_value = c(
        2 * min(mean(null_stats$median_hours <= observed_stats$median_hours, na.rm = TRUE),
                mean(null_stats$median_hours >= observed_stats$median_hours, na.rm = TRUE)),
        2 * min(mean(null_stats$mean_hours <= observed_stats$mean_hours, na.rm = TRUE),
                mean(null_stats$mean_hours >= observed_stats$mean_hours, na.rm = TRUE)),
        2 * min(mean(null_stats$prop_within_24h <= observed_stats$prop_within_24h, na.rm = TRUE),
                mean(null_stats$prop_within_24h >= observed_stats$prop_within_24h, na.rm = TRUE)),
        2 * min(mean(null_stats$prop_within_6h <= observed_stats$prop_within_6h, na.rm = TRUE),
                mean(null_stats$prop_within_6h >= observed_stats$prop_within_6h, na.rm = TRUE)),
        2 * min(mean(null_stats$prop_within_1h <= observed_stats$prop_within_1h, na.rm = TRUE),
                mean(null_stats$prop_within_1h >= observed_stats$prop_within_1h, na.rm = TRUE))
      ),
      # One-sided for avoidance hypothesis (longer intervals = avoidance)
      p_value_avoidance = c(
        mean(null_stats$median_hours <= observed_stats$median_hours, na.rm = TRUE),
        mean(null_stats$mean_hours <= observed_stats$mean_hours, na.rm = TRUE),
        mean(null_stats$prop_within_24h >= observed_stats$prop_within_24h, na.rm = TRUE),
        mean(null_stats$prop_within_6h >= observed_stats$prop_within_6h, na.rm = TRUE),
        mean(null_stats$prop_within_1h >= observed_stats$prop_within_1h, na.rm = TRUE)
      )
    )
    
    # Print results
    cat("\nRESULTS:\n")
    cat("Expected median hours:", round(p_values$expected_null[1], 1), "\n")
    cat("Difference from expected:", round(observed_stats$median_hours - p_values$expected_null[1], 1), "hours\n")
    
    # Focus on median results first
    median_p <- p_values$p_value[1]  # median_hours p-value
    median_avoidance_p <- p_values$p_value_avoidance[1]  # median avoidance p-value
    
    cat("Median hours p-value (two-sided):", round(median_p, 4), "\n")
    cat("Median hours p-value (avoidance):", round(median_avoidance_p, 4), "\n")
    
    # Interpret median results
    if(median_p < 0.05) {
      if(observed_stats$median_hours > p_values$expected_null[1]) {
        cat("\n*** SIGNIFICANT TEMPORAL AVOIDANCE DETECTED (MEDIAN) ***\n")
        cat("Otters wait significantly LONGER than expected by chance\n")
      } else {
        cat("\n*** SIGNIFICANT TEMPORAL ATTRACTION DETECTED (MEDIAN) ***\n") 
        cat("Otters arrive significantly SOONER than expected by chance\n")
      }
    } else {
      cat("\n*** NO SIGNIFICANT MEDIAN DIFFERENCE ***\n")
      cat("Median succession time is consistent with random expectation\n")
    }
    
    # Then check other metrics
    other_sig_results <- p_values %>% 
      filter(metric != "median_hours", p_value < 0.05)
    
    if(nrow(other_sig_results) > 0) {
      cat("\nOther significant patterns detected:\n")
      for(j in 1:nrow(other_sig_results)) {
        cat("  ", other_sig_results$metric[j], ": p =", round(other_sig_results$p_value[j], 4), "\n")
      }
      cat("(These may be driven by outliers/distribution shape rather than central tendency)\n")
    }
    
    # Store results
    results_list[[paste0("window_", window, "d")]] <- list(
      window_days = window,
      observed_stats = observed_stats,
      null_distribution = null_stats,
      p_values = p_values,
      observed_data = observed_data
    )
  }
  
  return(results_list)
}

# Plotting function
plot_randomization_results <- function(results, window_days = 7) {
  
  result <- results[[paste0("window_", window_days, "d")]]
  
  if(is.null(result)) {
    cat("No results found for", window_days, "day window\n")
    return(NULL)
  }
  
  # Create plots
  par(mfrow = c(2, 2))
  
  # 1. Null distribution for median hours
  hist(result$null_distribution$median_hours, 
       main = paste("Null Distribution: Median Hours\n(", window_days, "day window)"),
       xlab = "Median Hours Between Detections", 
       col = "lightblue", 
       breaks = 30)
  abline(v = result$observed_stats$median_hours, col = "red", lwd = 2)
  abline(v = mean(result$null_distribution$median_hours), col = "blue", lwd = 2, lty = 2)
  legend("topright", 
         c(paste("Observed:", round(result$observed_stats$median_hours, 1)), 
           paste("Expected:", round(mean(result$null_distribution$median_hours), 1))),
         col = c("red", "blue"), lwd = 2, lty = c(1, 2))
  
  # 2. Proportion within 24h
  hist(result$null_distribution$prop_within_24h, 
       main = "Null Distribution: Prop < 24h",
       xlab = "Proportion Within 24 Hours", 
       col = "lightgreen", 
       breaks = 20)
  abline(v = result$observed_stats$prop_within_24h, col = "red", lwd = 2)
  abline(v = mean(result$null_distribution$prop_within_24h), col = "blue", lwd = 2, lty = 2)
  
  # 3. Observed succession intervals
  hist(result$observed_data$hours_between, 
       main = paste("Observed Succession Intervals\n(n =", nrow(result$observed_data), ")"),
       xlab = "Hours Between Beaver and Otter", 
       col = "orange", 
       breaks = 20)
  abline(v = result$observed_stats$median_hours, col = "blue", lwd = 2, lty = 2)
  
  # 4. P-values
  barplot(result$p_values$p_value, 
          names.arg = result$p_values$metric,
          main = "P-values (Two-sided)",
          ylab = "P-value",
          col = ifelse(result$p_values$p_value < 0.05, "red", "gray"),
          las = 2)
  abline(h = 0.05, col = "red", lty = 2)
  
  par(mfrow = c(1, 1))
}

#===============================================================================
# RUN THE ANALYSIS
#===============================================================================

# Run the streamlined test
cat("TEMPORAL AVOIDANCE ANALYSIS WITH RANDOMIZATION TEST\n")
cat("==================================================\n")

results <- run_temporal_avoidance_test(activity_data, 
                                       n_permutations = 1000,
                                       time_windows = c(7))  # Start with just 7-day window

# Create plots
plot_randomization_results(results, window_days = 7)