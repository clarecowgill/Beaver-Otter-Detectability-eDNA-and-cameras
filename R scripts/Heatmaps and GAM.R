#_______________________________________________________________________________
# Complete GAM Analysis Script with Camera Type Random Effects
# Includes: Data wrangling, Initial Heatmap Plots, GAM fitting, model comparison, plotting, and diagnostics
#_______________________________________________________________________________
# Load required libraries
library(mgcv)
library(tidyverse)
library(lubridate)
library(patchwork)
library(emmeans)
library(lme4)

#_______________________________________________________________________________
# INITIAL DATA WRANGLING AND PREPARATION ----
#_______________________________________________________________________________

# Load and clean camera trap data ----

cams <- read.csv('cropton_cams_behaviour.csv')

# Clean camera data
cams <- cams %>%
  mutate(
    datetime = dmy_hm(date),
    date_only = as.Date(datetime),
    camera = str_extract(camera, "CROP\\d+")
  ) %>%
  filter(
    !is.na(datetime),                  # Remove parsing failures
    date_only >= as.Date("2023-09-29") # Keep only dates from 29 Sept 2023 onward
  )
#_______________________________________________________________________________
# Load and clean metadata with camera type preservation ----

# Read the original metadata CSV fresh and preserve camera type
meta_original <- read.csv('crop_meta.csv')

# Process meta correctly, preserving camera type
meta <- meta_original %>%
  mutate(
    sampling_datetime = dmy_hm(sampling_date),
    sampling_date_only = as.Date(sampling_datetime),
    # PRESERVE the original camera column as camera_type 
    camera_type = camera,  # This should contain "me" or "fe"
    # Extract CROP## for site matching
    camera_site = str_extract(sample, "CROP\\d+")
  ) %>%
  filter(!is.na(sampling_datetime)) %>%
  # Remove old camera column and rename for clarity
  select(-camera) %>%
  rename(camera = camera_site)

# Create camera type lookup
camera_type_lookup <- meta %>%
  distinct(camera, camera_type) %>%
  filter(!is.na(camera), !is.na(camera_type))

#_______________________________________________________________________________
# Proper Faulty Deployment Handling ----

# Faulty camera deployments (specific deployment periods that failed)
dodge_cams <- c("CROP_1001", "CROP_1103", "CROP_0108", "CROP_0208", "CROP_1303", 
                "CROP_0901", "CROP_0903", "CROP_0908", "CROP_0909","CROP_0910", "CROP_0801")

# Parse faulty deployments correctly
faulty_deployments <- data.frame(
  deployment = dodge_cams
) %>%
  mutate(
    camera = paste0("CROP", str_sub(deployment, 8, 9))
  )

# Get deployment periods with camera type preserved
meta_deployments <- meta %>%
  mutate(
    camera_extracted = paste0("CROP", str_sub(sample, 8, 9))
  ) %>%
  arrange(camera_extracted, sampling_date_only) %>%
  group_by(camera_extracted) %>%
  mutate(
    period_start = lag(sampling_date_only),
    period_end = sampling_date_only,
    period_start = if_else(is.na(period_start), 
                           min(cams$date_only[cams$camera == first(camera_extracted)], na.rm = TRUE), 
                           period_start + days(1))
  ) %>%
  ungroup() %>%
  mutate(
    is_faulty = sample %in% dodge_cams,
    camera = camera_extracted  # Use extracted camera name
  ) %>%
  select(sample, camera, period_start, period_end, is_faulty, camera_type)

# Filter camera data with camera type preserved
cams_with_deployment_info <- cams %>%
  mutate(original_row = row_number()) %>%
  crossing(
    meta_deployments %>% 
      select(camera_meta = camera, sample_meta = sample, period_start, period_end, is_faulty, camera_type)
  ) %>%
  filter(
    camera == camera_meta,
    date_only >= period_start,
    date_only <= period_end
  ) %>%
  group_by(original_row) %>%
  slice_head(n = 1) %>%
  ungroup()

# Final filtered data with camera type
cams_filtered <- cams_with_deployment_info %>%
  filter(!is_faulty) %>%
  select(-camera_meta, -sample_meta, -period_start, -period_end, -is_faulty, -original_row) %>%
  mutate(camera_type = factor(camera_type))

#_______________________________________________________________________________
# Load and Process eDNA Data ----
edna <- read.csv("crop_filtered.csv")

# Convert eDNA data from wide to long format
edna_long <- edna %>%
  # Convert to long format
  pivot_longer(cols = -X, names_to = "sample", values_to = "reads") %>%
  # Rename species column
  rename(species = X) %>%
  # Convert reads to detection (presence/absence)
  mutate(detected = ifelse(reads > 0, 1, 0)) %>%
  # Filter out control samples (NEG, POS, EB, FB)
  filter(!str_detect(sample, "NEG|POS|EB|FB")) %>%
  # JOIN WITH METADATA TO GET ACTUAL SAMPLING DATES
  left_join(meta %>% select(sample, sampling_date_only), by = "sample") %>%
  # Remove any samples without metadata
  filter(!is.na(sampling_date_only)) %>%
  # Rename and process the date column
  rename(sampling_date = sampling_date_only) %>%
  # Ensure sampling_date is in Date format
  mutate(
    sampling_date = as.Date(sampling_date),
    # Extract site information from sample name
    site = paste0("CROP", str_sub(sample, 8, 9))  # Get SS part for site
  ) %>%
  # Remove sites 11–14
  filter(!site %in% c("CROP11", "CROP12", "CROP13", "CROP14")) %>%
  # Handle multiple samples per site-date combination
  group_by(species, site, sampling_date) %>%
  summarise(
    detected = max(detected),  # If detected in any sample, count as detected
    reads = sum(reads),        # Sum reads if combining samples
    .groups = "drop"
  ) %>%
  arrange(species, site, sampling_date)

#_______________________________________________________________________________
# DETECTION HEATMAPS FUNCTION ----
#_______________________________________________________________________________
create_detection_heatmaps <- function(camera_data, edna_data) {
  
  # Site order: CROP02, CROP01, CROP03-10
  site_order <- paste0("CROP", sprintf("%02d", c(2, 1, 3:10)))
  
  # Process camera data
  camera_processed <- camera_data %>%
    filter(species %in% c("Castor_fiber", "Lutra_lutra")) %>%
    mutate(
      year_month = floor_date(date_only, "month"),
      month_year = paste(month.abb[month(year_month)], year(year_month)),
      site = camera
    ) %>%
    # Count detections per species, site, month
    group_by(species, site, year_month, month_year) %>%
    summarise(detections = n(), .groups = "drop") %>%
    # Complete grid for all sites and months
    complete(
      species,
      site = site_order,
      year_month = seq(from = as.Date("2023-09-01"), to = as.Date("2024-09-01"), by = "month")
    ) %>%
    mutate(
      detections = replace_na(detections, 0),
      month_year = paste(month.abb[month(year_month)], year(year_month)),
      month_year = factor(month_year, 
                          levels = paste(month.abb[c(9:12, 1:9)], 
                                         c(rep(2023, 4), rep(2024, 9))),
                          ordered = TRUE),
      species_label = case_when(
        species == "Castor_fiber" ~ "Beaver",
        species == "Lutra_lutra" ~ "Otter"
      ),
      site = factor(site, levels = site_order),
      sqrt_detections = sqrt(detections)
    )
  
  # Process eDNA data
  edna_processed <- edna_data %>%
    mutate(
      year_month = floor_date(sampling_date, "month"),
      reads = pmax(0, reads)
    ) %>%
    # Calculate proportional reads within each species-month
    group_by(site, sampling_date) %>%
    mutate(
      total_reads_sample = sum(reads, na.rm = TRUE),
      prop_reads = ifelse(total_reads_sample > 0, reads / total_reads_sample, 0)
    ) %>%
    ungroup() %>%
    # Complete grid for all sites and months
    complete(
      species,
      site = site_order,
      year_month = seq(from = as.Date("2023-09-01"), to = as.Date("2024-09-01"), by = "month")
    ) %>%
    mutate(
      prop_reads = replace_na(prop_reads, 0),
      month_year = paste(month.abb[month(year_month)], year(year_month)),
      month_year = factor(month_year, 
                          levels = paste(month.abb[c(9:12, 1:9)], 
                                         c(rep(2023, 4), rep(2024, 9))),
                          ordered = TRUE),
      species_label = case_when(
        species == "Castor_fiber" ~ "Beaver",
        species == "Lutra_lutra" ~ "Otter"
      ),
      site = factor(site, levels = site_order)
    )
  
  edna_processed <- edna_processed %>% 
    filter(species %in% c("Castor_fiber", "Lutra_lutra"))
  
  # Apply square root transformation for better scale visibility
  edna_processed <- edna_processed %>%
    mutate(
      sqrt_prop_reads = sqrt(prop_reads)
    )
  
  # Calculate shared scales
  camera_max <- max(camera_processed$sqrt_detections)
  edna_max <- max(edna_processed$sqrt_prop_reads, na.rm = TRUE)
  
  # Beaver Camera plot (shared camera scale)
  p1 <- camera_processed %>%
    filter(species_label == "Beaver") %>%
    ggplot(aes(x = month_year, y = site, fill = sqrt_detections)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(low = "white", high = "#2c7a56", name = "√Det",
                        limits = c(0, camera_max)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # Remove x labels (shared below)
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      plot.margin = margin(5, 2, 2, 5),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  # Beaver eDNA plot (shared eDNA scale)
  p2 <- edna_processed %>%
    filter(species_label == "Beaver") %>%
    ggplot(aes(x = month_year, y = site, fill = sqrt_prop_reads)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(low = "white", high = "#323b8a", name = "Prop",
                        limits = c(0, edna_max)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # Remove x labels (shared below)
      axis.text.y = element_blank(),  # Remove y labels (shared with left)
      axis.title = element_blank(),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      plot.margin = margin(5, 5, 2, 2),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  # Otter Camera plot (shared camera scale)
  p3 <- camera_processed %>%
    filter(species_label == "Otter") %>%
    ggplot(aes(x = month_year, y = site, fill = sqrt_detections)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(low = "white", high = "#2c7a56", name = "√Det",
                        limits = c(0, camera_max)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      plot.margin = margin(2, 2, 5, 5),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  # Otter eDNA plot (shared eDNA scale)
  p4 <- edna_processed %>%
    filter(species_label == "Otter") %>%
    ggplot(aes(x = month_year, y = site, fill = sqrt_prop_reads)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient(low = "white", high = "#323b8a", name = "Prop",
                        limits = c(0, edna_max)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_blank(),  # Remove y labels (shared with left)
      axis.title = element_blank(),
      plot.title = element_blank(),
      panel.grid = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      plot.margin = margin(2, 5, 5, 2),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  return(list(
    p1 = p1,
    p2 = p2,
    p3 = p3,
    p4 = p4,
    camera_data = camera_processed,
    edna_data = edna_processed
  ))
}

heatmaps <- create_detection_heatmaps(cams_filtered, edna_long)
library(patchwork)

p_final <- (heatmaps$p1 + heatmaps$p2) /
  (heatmaps$p3 + heatmaps$p4)
    

ggsave(filename = "Plots/heatmaps.png",
       plot = p_final +
         theme(
           panel.background = element_rect(fill = "white"),
           panel.border = element_blank(),
           plot.background = element_rect(fill = "white", color = NA)
         ),
       width = 6, height = 4, dpi = 300)

#_______________________________________________________________________________
# UNIFIED MODELLING WITH CAMERA TYPE RANDOM EFFECTS — MODEL-BASED SITE EFFECTS ----
#_______________________________________________________________________________

# --- helper: model-based site effects (works for eDNA & cameras) ----------------
# Returns one row per species x site with model-predicted prob and 95% CI,
# converted to (log) relative effect vs the species mean (matching your plots).
extract_site_effects_from_model <- function(unified_model, raw_data, method_name = "Camera") {
  
  stopifnot(!is.null(unified_model$best_model))
  best_mod <- unified_model$best_model
  
  # reference time: mean of observed time (on same units as time_numeric)
  ref_time <- mean(raw_data$time_numeric, na.rm = TRUE)
  
  # prediction grid: species x site at reference time
  pred_grid <- expand_grid(
    species = levels(raw_data$species),
    camera  = levels(raw_data$camera)
  ) %>%
    mutate(
      species = factor(species, levels = levels(raw_data$species)),
      camera  = factor(camera,  levels = levels(raw_data$camera)),
      time_numeric = ref_time
    )
  
  # if the model/data contain camera_type (weekly camera data), include a "typical" level
  if ("camera_type" %in% names(raw_data)) {
    pred_grid <- pred_grid %>%
      mutate(camera_type = fct_explicit_na(raw_data$camera_type)[1]) %>%   # first level (fe/me)
      mutate(camera_type = factor(camera_type, levels = levels(raw_data$camera_type)))
  }
  
  # predict on link scale with SE, then to response + CIs
  pr <- predict(best_mod, newdata = pred_grid, type = "link", se.fit = TRUE)
  
  eff <- pred_grid %>%
    mutate(
      fit_link   = pr$fit,
      se_link    = pr$se.fit,
      lower_link = fit_link - 1.96 * se_link,
      upper_link = fit_link + 1.96 * se_link,
      prob       = plogis(fit_link),
      lower.CL   = plogis(lower_link),
      upper.CL   = plogis(upper_link)
    )
  
  # species mean probability (at ref time; average across sites)
  species_means <- eff %>%
    group_by(species) %>%
    summarise(mean_prob = mean(prob, na.rm = TRUE), .groups = "drop")
  
  effects_df <- eff %>%
    left_join(species_means, by = "species") %>%
    mutate(
      relative_response      = prob / mean_prob,
      log_relative_response  = log(relative_response),
      log_lower_rel          = log(pmax(0.001, lower.CL / mean_prob)),
      log_upper_rel          = log(pmin(1000, upper.CL / mean_prob)),
      # caps for plotting stability (same as your current plots)
      log_relative_response  = pmax(-6, pmin(6, log_relative_response)),
      log_lower_rel          = pmax(-6, pmin(6, log_lower_rel)),
      log_upper_rel          = pmax(-6, pmin(6, log_upper_rel)),
      # significance vs the within-species mean on the response scale
      sig_vs_average         = !(lower.CL <= mean_prob & upper.CL >= mean_prob)
    )
  
  # for eDNA downstream you expect a `site` column (you rename back to camera when plotting)
  if (method_name == "eDNA") {
    effects_df <- effects_df %>% rename(site = camera)
  }
  
  list(effects_df = effects_df)
}
#_______________________________________________________________________________
# --- model comparison ----------------------------------------------
compare_unified_models_with_camera_type <- function(data_long, method_name = "Camera", 
                                                    k_camera = 10, k_edna = 8) {
  data_long <- data_long %>% filter(species %in% c("Castor_fiber", "Lutra_lutra"))
  
  if (method_name != "eDNA" && !"camera_type" %in% names(data_long)) {
    data_long <- data_long %>% left_join(camera_type_lookup, by = "camera")
    if ("camera_type" %in% names(data_long)) {
      data_long <- data_long %>% mutate(camera_type = factor(camera_type))
    }
  }
  
  if (method_name == "eDNA") {
    cat("eDNA analysis - camera_type not applicable\n")
    cat(sprintf("Using k = %d for eDNA models\n", k_edna))
  } else {
    cat(sprintf("Using k = %d for %s models\n", k_camera, method_name))
  }
  k_value <- ifelse(method_name == "eDNA", k_edna, k_camera)
  
  # prepare data
  if (method_name == "Camera" || method_name == "Camera Trap" || str_detect(method_name, "Weekly")) {
    individual_data <- data_long %>%
      mutate(week_start = floor_date(date_only, "week", week_start = 1)) %>%
      group_by(species, camera, camera_type, week_start) %>%
      summarise(detected = ifelse(any(!is.na(camera)), 1, 0), .groups = "drop") %>%
      complete(
        nesting(species),
        nesting(camera, camera_type),
        week_start = seq(from = min(week_start, na.rm = TRUE), 
                         to = max(week_start, na.rm = TRUE), by = "week"),
        fill = list(detected = 0)
      ) %>%
      filter(!is.na(week_start))
  } else if (method_name == "Camera Monthly") {
    individual_data <- data_long %>% rename(week_start = date_only) %>% filter(!is.na(week_start))
  } else {
    individual_data <- data_long %>%
      filter(!is.na(sampling_date)) %>%
      rename(week_start = sampling_date, camera = site)
  }
  
  model_data <- individual_data %>%
    arrange(species, if ("site" %in% names(individual_data)) site else camera, week_start) %>%
    mutate(
      time_numeric = as.numeric(week_start - min(week_start)),
      species = factor(species),
      camera  = factor(if ("site" %in% names(individual_data)) site else camera),
      camera_type = if ("camera_type" %in% names(individual_data)) factor(camera_type) else NULL
    )
  
  models_list <- list(); model_names <- character()
  
  if (method_name == "eDNA") {
    tryCatch({
      mod1 <- gam(detected ~ species, data = model_data, family = binomial(), method = "REML",
                  optimizer = c("outer","bfgs"), control = gam.control(maxit = 200))
      models_list$species_only <- mod1; model_names <- c(model_names, "Species Only")
    }, error=function(e){})
    
    tryCatch({
      mod2 <- gam(detected ~ species + s(camera, bs = "re"),
                  data = model_data, family = binomial(), method = "REML",
                  optimizer = c("outer","bfgs"), control = gam.control(maxit = 200))
      models_list$species_site_re <- mod2; model_names <- c(model_names, "Species + Site RE")
    }, error=function(e){})
    
    tryCatch({
      mod3 <- gam(detected ~ species + s(time_numeric, k = k_value) + s(camera, bs = "re"),
                  data = model_data, family = binomial(), method = "REML",
                  optimizer = c("outer","bfgs"), control = gam.control(maxit = 200))
      models_list$species_time_site_re <- mod3; model_names <- c(model_names, "Species + Time + Site RE")
    }, error=function(e){})
    
    if (length(unique(model_data$time_numeric)) >= 6) {
      tryCatch({
        mod4 <- gam(detected ~ species + s(time_numeric, by = species, k = k_value) + s(camera, bs = "re"),
                    data = model_data, family = binomial(), method = "REML",
                    optimizer = c("outer","bfgs"), control = gam.control(maxit = 200))
        models_list$species_species_time_site_re <- mod4
        model_names <- c(model_names, "Species + Species Time + Site RE")
      }, error=function(e){})
    }
    
    tryCatch({
      mod5 <- gam(detected ~ species + camera + s(time_numeric, k = k_value),
                  data = model_data, family = binomial(), method = "REML",
                  optimizer = c("outer","bfgs"), control = gam.control(maxit = 200))
      models_list$species_camera_simple <- mod5
      model_names <- c(model_names, "Species + Site Simple + Time")
    }, error=function(e){})
    
    tryCatch({
      mod6 <- gam(detected ~ species + s(time_numeric, by = species, k = k_value) +
                    s(camera, by = species, bs = "re"),
                  data = model_data, family = binomial(), method = "REML",
                  optimizer = c("outer","bfgs"), control = gam.control(maxit = 200))
      models_list$species_time_site_speciesRE <- mod6
      model_names <- c(model_names, "Species + Species Time + Species-specific Site RE")
    }, error=function(e){})
    
  } else {
    tryCatch({
      mod1 <- gam(detected ~ species * camera, data = model_data, family = binomial(), method = "REML",
                  optimizer = c("outer","bfgs"), control = gam.control(maxit = 200))
      models_list$species_camera_fixed <- mod1; model_names <- c(model_names, "Species * Camera Fixed")
    }, error=function(e){})
    
    tryCatch({
      mod2 <- gam(detected ~ species * camera + s(time_numeric, k = k_value),
                  data = model_data, family = binomial(), method = "REML",
                  optimizer = c("outer","bfgs"), control = gam.control(maxit = 200))
      models_list$species_camera_time <- mod2; model_names <- c(model_names, "Species * Camera + Time")
    }, error=function(e){})
    
    tryCatch({
      mod3 <- gam(detected ~ species * camera + s(time_numeric, by = species, k = k_value),
                  data = model_data, family = binomial(), method = "REML",
                  optimizer = c("outer","bfgs"), control = gam.control(maxit = 200))
      models_list$species_camera_species_time <- mod3; model_names <- c(model_names, "Species * Camera + Species Time")
    }, error=function(e){})
    
    tryCatch({
      mod4 <- gam(detected ~ species + s(time_numeric, by = species, k = k_value) + s(camera, bs = "re"),
                  data = model_data, family = binomial(), method = "REML",
                  optimizer = c("outer","bfgs"), control = gam.control(maxit = 200))
      models_list$species_time_camera_re <- mod4; model_names <- c(model_names, "Species + Species Time + Camera RE")
    }, error=function(e){})
  }
  
  if (length(models_list) == 0) return(NULL)
  
  aic_values <- sapply(models_list, function(m) tryCatch(AIC(m), error=function(e) 1e6))
  dev_expl   <- sapply(models_list, function(m) round(summary(m)$dev.expl * 100, 1))
  conv_stat  <- sapply(models_list, function(m) tryCatch({ if (!m$converged) "WARNING" else "OK" }, error=function(e) "ERROR"))
  
  model_comparison <- tibble(
    Model = model_names,
    AIC = round(aic_values, 2),
    Deviance_Explained = dev_expl,
    df_residual = sapply(models_list, function(m) round(m$df.residual, 1)),
    Convergence = conv_stat
  ) %>%
    filter(!is.na(AIC)) %>%
    mutate(Delta_AIC = AIC - min(AIC, na.rm = TRUE),
           AIC_Weight = exp(-0.5 * Delta_AIC) / sum(exp(-0.5 * Delta_AIC), na.rm = TRUE)) %>%
    arrange(AIC)
  
  best_candidates <- model_comparison %>% filter(Convergence != "ERROR") %>% arrange(AIC)
  if (nrow(best_candidates) == 0) return(NULL)
  
  best_model_name <- best_candidates$Model[1]
  best_model_key  <- names(models_list)[which(model_names == best_model_name)]
  best_model      <- models_list[[best_model_key]]
  
  list(
    method = method_name,
    model_comparison = model_comparison,
    raw_data = model_data,
    models = models_list,
    best_model = best_model,
    best_model_name = best_model_name,
    convergence_status = best_candidates$Convergence[1],
    k_used = k_value
  )
}

#_______________________________________________________________________________
# --- predictions over time ----
# --- NEW: temporal predictions for camera models (fixed-effect site model) ----
create_temporal_predictions_with_camera_type <- function(model_result, raw_data,
                                                         extend_to = NULL, ci_level = 0.95) {
  best_mod <- model_result$best_model
  z <- qnorm((1 + ci_level) / 2)
  
  # build 100-point time grid
  if (is.null(extend_to)) {
    max_date <- max(raw_data$week_start, na.rm = TRUE)
  } else {
    max_date <- extend_to
  }
  
  t_seq <- seq(min(raw_data$week_start, na.rm = TRUE),
               max_date, length.out = 100)
  
  grid <- expand_grid(
    species      = levels(raw_data$species),
    week_start   = t_seq
  ) %>%
    mutate(
      time_numeric = as.numeric(week_start - min(raw_data$week_start, na.rm = TRUE)),
      species      = factor(species, levels = levels(raw_data$species))
    )
  
  # predict per site first, then average
  site_grid <- grid %>%
    crossing(camera = levels(raw_data$camera)) %>%
    mutate(camera = factor(camera, levels = levels(raw_data$camera)))
  
  pr <- predict(best_mod, newdata = site_grid, type = "link", se.fit = TRUE)
  
  site_grid <- site_grid %>%
    mutate(
      fit_link   = pr$fit,
      se_link    = pr$se.fit,
      lower_link = fit_link - z * se_link,
      upper_link = fit_link + z * se_link,
      prob       = plogis(fit_link),
      lower_prob = plogis(lower_link),
      upper_prob = plogis(upper_link)
    )
  
  # mean across cameras (simple unweighted mean per species × time)
  out <- site_grid %>%
    group_by(species, time_numeric, week_start) %>%
    summarise(
      detection_prob = mean(prob, na.rm = TRUE),
      lower_ci       = mean(lower_prob, na.rm = TRUE),
      upper_ci       = mean(upper_prob, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      species_label = ifelse(species == "Castor_fiber", "Beaver", "Otter")
    )
  
  return(out)
}

# Smooth-only CI for temporal trend, compatible with your plotting columns
predict_temporal_partial_compat <- function(best_model, raw_data, method_label,
                                            drop_site = NULL, n_grid = 200) {
  # Ensure factors are consistent
  rd <- raw_data %>%
    mutate(
      species = factor(species, levels = levels(raw_data$species)),
      camera  = factor(camera,  levels = levels(raw_data$camera))
    )
  
  # Build prediction grid: species x time x site
  t_seq <- seq(min(rd$time_numeric, na.rm = TRUE),
               max(rd$time_numeric, na.rm = TRUE), length.out = n_grid)
  
  grid <- tidyr::expand_grid(
    species      = levels(rd$species),
    time_numeric = t_seq,
    camera       = levels(rd$camera)
  )
  
  # Drop CROP02 only from averaging grid (if present)
  if (!is.null(drop_site) && drop_site %in% levels(rd$camera)) {
    grid <- dplyr::filter(grid, camera != drop_site)
  }
  
  # Back to dates
  t0 <- min(rd$week_start, na.rm = TRUE)
  grid <- dplyr::mutate(grid, week_start = t0 + time_numeric)
  
  # Linear predictor pieces
  X   <- predict(best_model, newdata = grid, type = "lpmatrix")
  Vp  <- vcov(best_model)
  bet <- coef(best_model)
  
  # Keep only columns for s(time_numeric) smooths (works for species-specific too)
  keep_cols <- grepl("s\\(time_numeric\\)", colnames(X))
  if (!any(keep_cols)) stop("No s(time_numeric) smooth found in model; cannot build smooth-only CI.")
  
  X_smooth <- X
  X_smooth[, !keep_cols] <- 0
  
  fixed_part  <- as.numeric(X %*% bet) - as.numeric(X_smooth %*% bet)
  smooth_mean <- as.numeric(X_smooth %*% bet)
  smooth_var  <- rowSums((X_smooth %*% Vp) * X_smooth)
  smooth_se   <- sqrt(pmax(smooth_var, 0))
  
  fit_link  <- fixed_part + smooth_mean
  low_link  <- fixed_part + (smooth_mean - 1.96 * smooth_se)
  high_link <- fixed_part + (smooth_mean + 1.96 * smooth_se)
  
  grid %>%
    mutate(
      fit   = plogis(fit_link),
      lower = plogis(low_link),
      upper = plogis(high_link)
    ) %>%
    group_by(species, time_numeric, week_start) %>%
    summarise(
      detection_prob = mean(fit,   na.rm = TRUE),
      lower_ci       = mean(lower, na.rm = TRUE),
      upper_ci       = mean(upper, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      method = method_label,
      species_label = ifelse(species == "Castor_fiber", "Beaver", "Otter")
    )
}

#_______________________________________________________________________________
# --- main runner: now uses MODEL-BASED site effects for ALL methods ----
run_unified_analysis_with_camera_type <- function(camera_data, edna_data, 
                                                  k_weekly = 8, k_monthly = 6, k_edna = 6) {
  
  cat("=== RUNNING ANALYSIS WITH CUSTOM K VALUES ===\n")
  cat(sprintf("Camera Weekly k = %d\nCamera Monthly k = %d\neDNA k = %d\n", k_weekly, k_monthly, k_edna))
  cat("===============================================\n\n")
  
  unified_camera_weekly <- compare_unified_models_with_camera_type(camera_data, "Camera Trap", k_camera = k_weekly)
  
  monthly_camera_data <- camera_data %>%
    filter(species %in% c("Castor_fiber", "Lutra_lutra")) %>%
    mutate(year_month = floor_date(date_only, "month"),
           year = year(date_only), month = month(date_only)) %>%
    filter(!(month == 8 & year == 2023), !(month == 9 & year == 2024)) %>%
    group_by(species, camera, camera_type, year_month) %>%
    summarise(detected = 1, .groups = "drop") %>%
    left_join(camera_data %>% distinct(camera, camera_type) %>% filter(!is.na(camera_type)), by = "camera") %>%
    mutate(camera_type = coalesce(camera_type.x, camera_type.y)) %>%
    select(-camera_type.x, -camera_type.y) %>%
    complete(
      nesting(species),
      nesting(camera, camera_type),
      year_month = seq(from = as.Date("2023-09-01"), to = as.Date("2024-08-01"), by = "month"),
      fill = list(detected = 0)
    ) %>%
    filter(!is.na(year_month), !is.na(camera_type)) %>%
    rename(date_only = year_month)
  
  unified_camera_monthly <- compare_unified_models_with_camera_type(monthly_camera_data, "Camera Monthly", k_camera = k_monthly)
  unified_edna          <- compare_unified_models_with_camera_type(edna_data, "eDNA", k_edna = k_edna)
  
  cat("\n=== MODEL CONVERGENCE STATUS ===\n")
  cat(sprintf("Camera Weekly: %s\n", unified_camera_weekly$convergence_status))
  cat(sprintf("Camera Monthly: %s\n", unified_camera_monthly$convergence_status))
  cat(sprintf("eDNA: %s\n", unified_edna$convergence_status))
  
  # separation check (unchanged)
  cat("\n=== CHECKING FOR SEPARATION ISSUES ===\n")
  edna_separation_check <- unified_edna$raw_data %>%
    group_by(species, camera) %>%
    summarise(n_obs = n(), n_detected = sum(detected), detection_rate = mean(detected), .groups = "drop") %>%
    mutate(separation_issue = detection_rate %in% c(0, 1))
  cat(sprintf("Separation found in %d out of %d site-species combinations\n",
              sum(edna_separation_check$separation_issue), nrow(edna_separation_check)))

  # edna_separation_check <- unified_edna$raw_data %>%
  #   group_by(species, camera) %>%
  #   summarise(
  #     n_obs  = n(),
  #     n_zero = sum(prop == 0, na.rm = TRUE),
  #     n_one  = sum(prop == 1, na.rm = TRUE),
  #     prop_min = min(prop, na.rm = TRUE),
  #     prop_max = max(prop, na.rm = TRUE),
  #     .groups = "drop"
  #   ) %>%
  #   mutate(
  #     # Separation-like issue if a site has ONLY zeros or ONLY ones
  #     separation_issue = (prop_min == 0 & prop_max == 0) | (prop_min == 1 & prop_max == 1)
  #   )
  
  # === MODEL-BASED SITE EFFECTS here ===
  camera_weekly_effects <- extract_site_effects_from_model(unified_camera_weekly, unified_camera_weekly$raw_data, "Camera Trap")
  camera_monthly_effects <- extract_site_effects_from_model(unified_camera_monthly, unified_camera_monthly$raw_data, "Camera Monthly")
  edna_effects <- extract_site_effects_from_model(unified_edna, unified_edna$raw_data, "eDNA")
  
  # temporal predictions (unchanged)
  camera_weekly_temporal  <- create_temporal_predictions_with_camera_type(unified_camera_weekly, unified_camera_weekly$raw_data, as.Date("2024-09-02"))
  camera_monthly_temporal <- create_temporal_predictions_with_camera_type(unified_camera_monthly, unified_camera_monthly$raw_data, as.Date("2024-08-01"))
  # Build the drop_site robustly in case of level naming differences
  edna_levels <- levels(unified_edna$raw_data$camera)
  drop_site_val <- if ("CROP02" %in% edna_levels) "CROP02" else edna_levels[2]
  
  edna_temporal <- predict_temporal_partial_compat(
    best_model = unified_edna$best_model,
    raw_data   = unified_edna$raw_data,
    method_label = "eDNA (Sampling Events)"  )
  
  # --- your existing plotting code below will now use model-based effects ---
  # Weekly CT vs eDNA
  if (!is.null(camera_weekly_temporal) && !is.null(edna_temporal)) {
    camera_weekly_temporal <- camera_weekly_temporal %>% mutate(species_label = as.character(species_label), method = "Camera Trap (Weekly)")
    edna_temporal <- edna_temporal %>% mutate(species_label = as.character(ifelse(species == "Castor_fiber", "Beaver", "Otter")), method = "eDNA (Sampling Events)")
    weekly_plot_data <- bind_rows(camera_weekly_temporal, edna_temporal) %>%
      filter(week_start >= as.Date("2023-09-01")) %>%
      mutate(method = factor(method, levels = c("Camera Trap (Weekly)", "eDNA (Sampling Events)")),
             species_label = factor(species_label, levels = c("Beaver", "Otter")))
    p_weekly_vs_edna <- ggplot(weekly_plot_data, aes(x = week_start, y = detection_prob, color = species_label, fill = species_label)) +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, color = NA) +
      geom_line(size = 1.2) +
      facet_wrap(~method, scales = "free_x") +
      scale_color_manual(values = c("Beaver" = "#1e3a8a", "Otter" = "#2a9d8f")) +
      scale_fill_manual(values = c("Beaver" = "#1e3a8a", "Otter" = "#2a9d8f")) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
      labs(title = "Temporal Detection Patterns: Camera Traps vs eDNA",
           subtitle = "GAM models with improved convergence and optimization",
           x = NULL, y = "Detection Probability", color = NULL, fill = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            legend.position = "bottom",
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 10, color = "gray60"))
  } else {
    p_weekly_vs_edna <- NULL
  }
  
  # Monthly CT vs eDNA
  if (!is.null(camera_monthly_temporal) && !is.null(edna_temporal)) {
    camera_monthly_temporal <- camera_monthly_temporal %>% mutate(species_label = as.character(species_label), method = "Camera Trap (Monthly)")
    edna_temporal_monthly   <- edna_temporal %>% mutate(species_label = as.character(ifelse(species == "Castor_fiber","Beaver","Otter")), method = "eDNA (Sampling Events)")
    monthly_plot_data <- bind_rows(camera_monthly_temporal, edna_temporal_monthly) %>%
      filter(week_start >= as.Date("2023-09-01"), week_start <= as.Date("2024-08-01")) %>%
      mutate(method = factor(method, levels = c("Camera Trap (Monthly)", "eDNA (Sampling Events)")),
             species_label = factor(species_label, levels = c("Beaver", "Otter")))
    p_monthly_vs_edna <- ggplot(monthly_plot_data, aes(x = week_start, y = detection_prob, color = species_label, fill = species_label)) +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.3, color = NA) +
      geom_line(size = 1.2) +
      facet_wrap(~method, scales = "free_x") +
      scale_color_manual(values = c("Beaver" = "#1e3a8a", "Otter" = "#2a9d8f")) +
      scale_fill_manual(values = c("Beaver" = "#1e3a8a", "Otter" = "#2a9d8f")) +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
      labs(title = "Monthly Detection Patterns: Camera Traps vs eDNA",
           subtitle = "Note: Monthly models may have convergence issues",
           x = NULL, y = "Detection Probability", color = NULL, fill = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            legend.position = "bottom",
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 10, color = "gray60"))
  } else {
    p_monthly_vs_edna <- NULL
  }
  
  #_______________________________________________________________________________
  #--- Site effects plots ---------------------------
  site_order <- rev(c("CROP02", "CROP01", paste0("CROP", sprintf("%02d", 3:10))))
  weekly_data    <- camera_weekly_effects$effects_df %>% mutate(method = "Camera Trap (Weekly)")
  edna_data_plot <- edna_effects$effects_df %>% rename(camera = site) %>% mutate(method = "eDNA")
  weekly_data <- weekly_data %>%
    mutate(camera = factor(camera, levels = site_order))
  edna_data_plot <- edna_data_plot %>%
    mutate(camera = factor(camera, levels = site_order))
  combined_data_week <- bind_rows(weekly_data, edna_data_plot) %>%
    mutate(
      method = factor(method, levels = c("Camera Trap (Weekly)", "eDNA")),
      camera = factor(camera, levels = site_order)
    )
  species_differences_week <- combined_data_week %>%
    filter(!is.infinite(log_relative_response)) %>%
    select(camera, method, species, log_relative_response, log_lower_rel, log_upper_rel) %>%
    pivot_wider(names_from = species, values_from = c(log_relative_response, log_lower_rel, log_upper_rel), names_sep = "_") %>%
    filter(!is.na(log_relative_response_Castor_fiber) & !is.na(log_relative_response_Lutra_lutra)) %>%
    mutate(
      log_diff  = log_relative_response_Castor_fiber - log_relative_response_Lutra_lutra,
      se_beaver = (log_upper_rel_Castor_fiber - log_lower_rel_Castor_fiber)/(2*1.96),
      se_otter  = (log_upper_rel_Lutra_lutra - log_lower_rel_Lutra_lutra)/(2*1.96),
      se_diff   = sqrt(se_beaver^2 + se_otter^2),
      diff_lower = log_diff - 1.96*se_diff,
      diff_upper = log_diff + 1.96*se_diff,
      sig_between_species = !(diff_lower <= 0 & diff_upper >= 0)
    ) %>%
    select(camera, method, sig_between_species)

  combined_data_week <- combined_data_week %>%
    left_join(species_differences_week, by = c("camera","method")) %>%
    mutate(sig_between_species = replace_na(sig_between_species, FALSE))

  combined_data_week_fixed <- combined_data_week %>%
    mutate(
      log_relative_response = ifelse(is.infinite(log_relative_response), ifelse(log_relative_response < 0, -6, 6), log_relative_response),
      log_lower_rel = ifelse(is.infinite(log_lower_rel), ifelse(log_lower_rel < 0, -6, 6), log_lower_rel),
      log_upper_rel = ifelse(is.infinite(log_upper_rel), ifelse(log_upper_rel < 0, -6, 6), log_upper_rel),
      log_lower_rel = pmax(log_lower_rel, -6),
      log_upper_rel = pmin(log_upper_rel,  6),
      alpha_level   = ifelse(sig_vs_average, 1.0, 0.4)
    )

  axis_limits_week <- combined_data_week_fixed %>%
    group_by(method) %>%
    summarise(min_x = min(log_lower_rel, na.rm = TRUE), max_x = max(log_upper_rel, na.rm = TRUE), .groups = "drop") %>%
    mutate(x_range = max_x - min_x, asterisk_x = max_x + x_range*0.08, min_x_padded = min_x - x_range*0.05, max_x_padded = max_x + x_range*0.15)

  asterisk_positions_week <- combined_data_week_fixed %>%
    group_by(camera, method) %>%
    filter(any(sig_between_species, na.rm = TRUE)) %>%
    summarise(sig_between_species = first(sig_between_species[!is.na(sig_between_species)]), .groups = "drop") %>%
    filter(sig_between_species == TRUE) %>%
    left_join(axis_limits_week %>% select(method, asterisk_x), by = "method")

  p_site_effects_weekly <- ggplot(combined_data_week_fixed, aes(y = camera, x = log_relative_response, color = species)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    geom_point(aes(alpha = I(alpha_level)), size = 3, position = position_dodge(width = 0.4)) +
    geom_errorbarh(aes(xmin = log_lower_rel, xmax = log_upper_rel, alpha = I(alpha_level)), height = 0.2, position = position_dodge(width = 0.4)) +
    { if (nrow(asterisk_positions_week) > 0)
      geom_text(data = asterisk_positions_week, aes(x = asterisk_x, y = camera, label = "*"),
                color = "black", size = 6, fontface = "bold", inherit.aes = FALSE) } +
    geom_blank(data = axis_limits_week, aes(x = min_x_padded, y = "CROP01"), inherit.aes = FALSE) +
    geom_blank(data = axis_limits_week, aes(x = max_x_padded, y = "CROP01"), inherit.aes = FALSE) +
    facet_wrap(~method, ncol = 1, scales = "free_y") + 
    coord_flip() + 
    scale_color_manual(values = c("Castor_fiber" = "#1e3a8a", "Lutra_lutra" = "#2a9d8f"),
                       labels = c("Castor_fiber" = "Beaver", "Lutra_lutra" = "Otter")) +
    labs(title = "Site-Specific Detection Effects: Camera Traps vs eDNA",
         subtitle = "Faded colors = not significantly different from species average | * = significant difference between species",
         x = "Log Relative Likelihood (Compared to Species Average)", y = "Site", color = "Species") +
    theme_minimal() +
    theme(legend.position = "bottom", strip.text = element_text(size = 12),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "gray60"))
  #_______________________________________________________________________________
  
  # show plots
  if (!is.null(p_weekly_vs_edna)) print(p_weekly_vs_edna)
  if (!is.null(p_monthly_vs_edna)) print(p_monthly_vs_edna)
  print(p_site_effects_weekly)
  
  list(
    unified_models = list(camera_weekly = unified_camera_weekly,
                          camera_monthly = unified_camera_monthly,
                          edna = unified_edna),
    site_effects = list(camera_weekly  = camera_weekly_effects,
                        camera_monthly = camera_monthly_effects,
                        edna = edna_effects),
    predictions = list(camera_weekly_temporal  = camera_weekly_temporal,
                       camera_monthly_temporal = camera_monthly_temporal,
                       edna_temporal           = edna_temporal),
    plots = list(weekly_vs_edna = p_weekly_vs_edna,
                 monthly_vs_edna = p_monthly_vs_edna,
                 site_effects_weekly = p_site_effects_weekly),
    k_values_used = list(weekly = k_weekly, monthly = k_monthly, edna = k_edna),
    separation_check = edna_separation_check
  )
}
#_______________________________________________________________________________
# RUN THE COMPLETE UPDATED ANALYSIS ----
#_______________________________________________________________________________

unified_results <- run_unified_analysis_with_camera_type(
  cams_filtered, edna_long,
  k_weekly = 8,
  k_monthly = 8,
  k_edna = 6
)

summary(unified_results$unified_models$camera_weekly$best_model)
summary(unified_results$unified_models$camera_monthly$best_model)
summary(unified_results$unified_models$edna$best_model)

gam.check(unified_results$unified_models$camera_weekly$best_model)
gam.check(unified_results$unified_models$camera_monthly$best_model)
gam.check(unified_results$unified_models$edna$best_model)

print(unified_results$unified_models$camera_weekly$model_comparison)
print(unified_results$unified_models$camera_monthly$model_comparison)
print(unified_results$unified_models$edna$model_comparison)


# Save all plots
if (!is.null(unified_results$plots$weekly_vs_edna)) {
  ggsave("Plots/final_temporal_weekly.png", 
         plot = unified_results$plots$weekly_vs_edna,
         width = 7.5, height = 4.8, dpi = 300, bg = "white")
}

if (!is.null(unified_results$plots$monthly_vs_edna)) {
  ggsave("Plots/final_temporal_monthly.png", 
         plot = unified_results$plots$monthly_vs_edna,
         width = 10, height = 6, dpi = 300, bg = "white")
}

if (!is.null(unified_results$plots$site_effects_weekly)) {
  ggsave("Plots/final_site_effects_weekly.png", 
         plot = unified_results$plots$site_effects_weekly,
         width = 8, height = 6, dpi = 300, bg = "white")
}

if (!is.null(unified_results$plots$site_effects_monthly)) {
  ggsave("Plots/final_site_effects_monthly.png", 
         plot = unified_results$plots$site_effects_monthly,
         width = 10, height = 8, dpi = 300, bg = "white")
}

# Print final summary
cat(sprintf("Camera Weekly: %s - %s\n", 
            unified_results$unified_models$camera_weekly$convergence_status,
            unified_results$unified_models$camera_weekly$best_model_name))
cat(sprintf("Camera Monthly: %s - %s\n", 
            unified_results$unified_models$camera_monthly$convergence_status,
            unified_results$unified_models$camera_monthly$best_model_name))
cat(sprintf("eDNA: %s - %s\n", 
            unified_results$unified_models$edna$convergence_status,
            unified_results$unified_models$edna$best_model_name))

# Show GAM carpet plots----
plot(unified_results$unified_models$camera_weekly$best_model, pages = 1)
plot(unified_results$unified_models$edna$best_model, pages = 1)

