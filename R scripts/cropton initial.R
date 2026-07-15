library(tidyverse)
library(tibble)
crop = read.csv('Cropton_v2_blast98_denoise.tsv', sep = '\t', row.names = 1, header = T)

# to take out dundreggan data
dun = crop[, c(172:257)]
dun = dun[rowSums(dun[, !colnames(dun) %in% 'taxonomy']) != 0, ]

# remove non-cropton data
crop = crop[, -c(172:234)]
write.csv(crop, 'crop_dat.csv')

# LOD threshold
blank_association <- read.csv("blank_association.csv")
crop_long <- as.data.frame(as.table(as.matrix(crop)))
colnames(crop_long) <- c("species", "sample", "reads")
crop_long$reads <- as.numeric(crop_long$reads)
blank_reads <- crop_long %>% 
  filter(sample %in% blank_association$blank)
blank_stats <- blank_reads %>% 
  group_by(species) %>% 
  summarise(
    mean_blank_reads = mean(reads, na.rm = TRUE),
    sd_blank_reads   = sd(reads, na.rm = TRUE)
  ) %>% 
  mutate(
    LOD = mean_blank_reads + 3 * sd_blank_reads
  )
non_blank <- crop_long %>% 
  filter(!sample %in% blank_association$blank)
lod_data <- non_blank %>%
  left_join(blank_stats, by = "species") %>%
  mutate(
    reads = ifelse(is.na(LOD) | reads >= LOD, reads, 0)
  )
lod_wide <- lod_data %>%
  select(species, sample, reads) %>%
  pivot_wider(names_from = sample, values_from = reads, values_fill = 0) %>%
  as.data.frame()
rownames(lod_wide) <- lod_wide$species
lod_wide$species <- NULL

crop <- lod_wide
crop = crop[, sapply(crop, is.numeric)]

# 0.1% threshold of total reads in a sample:
crop[t(t(crop) / colSums(crop)) < 0.001] = 0

# drop OTUs with no reads:
crop = crop[rowSums(crop[, !colnames(crop) %in% 'taxonomy']) != 0, ]

# drop unassigned:
crop = crop[!rownames(crop) %in% 'unassigned',]

# write csv to be used for initial analysis
write.csv(crop, 'crop_filtered.csv')

# species composition (proportion reads)
sapply(crop, class)
crop_prc <- t(crop) / colSums(crop)

rows_to_remove <- grepl("POS|NEG|EB|FB", rownames(crop_prc))
prc_filtered <- crop_prc[!rows_to_remove, ]

write.csv(prc_filtered, 'crop_prc.csv')

# ______________________________________________________________________________
library(ggplot2)
library(dplyr)
library(tidyverse)

crop_tax <- read.csv("Cropton_v2_blast98_denoise.tsv",
                     sep = "\t", row.names = 1, header = TRUE)
crop_joined <- crop %>%
  mutate(species = rownames(crop)) %>%
  left_join(
    crop_tax %>% 
      select(taxonomy) %>% 
      mutate(species = rownames(crop_tax)),
    by = "species"
  )

# Restore rownames and drop species column if desired
rownames(crop_joined) <- crop_joined$species
crop_joined$species <- NULL

crop_joined <- crop
# Convert row names to a column named 'species'
crop$species = rownames(crop)

# use reshape to convert from wide to long format
crop_long <- crop %>%
  pivot_longer(
    cols = -c(species, taxonomy),
    names_to = "sample",
    values_to = "reads"
  ) %>%
  filter(!species %in% c(
    "Homo_sapiens", "Bos_taurus", "Capra_hircus",
    "Ovis_aries", "Sus_scrofa", "Canis_lupus", "Equus_caballus"
  ))


rownames(crop_long) = NULL

# ______________________________________________________________________________

# Function to generate exploration plots for a different taxonomic groups

exploration_plots = function(data, taxonomic_group, title_prefix) {
  # Filter data based on taxonomic group using grep to include only the specified group
  filtered_data = data[grep(taxonomic_group, data$taxonomy), ]
  
  # Check if filtered_data is empty
  if (nrow(filtered_data) == 0) {
    cat("No data found for taxonomic group:", taxonomic_group, "\n")
    return(NULL)
  }
  
  # Stacked bar plot
  stacked_plot = ggplot(filtered_data, aes(x = sample, y = reads, fill = species)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = paste0(title_prefix, " (a)")) +
    scale_fill_manual(values = rainbow(length(unique(filtered_data$species)))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  print(stacked_plot)
  
  # Filled bar plot
  filled_plot = ggplot(filtered_data, aes(x = sample, fill = species, y = reads)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = rainbow(length(unique(filtered_data$species)))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  print(filled_plot)
  
  # Bubble plot
  filtered_data = filtered_data %>%
    replace(., . == 0, NA)
  
  bubble_plot = ggplot(filtered_data, aes(x = sample, y = species, size = reads)) +
    geom_point(alpha = 0.4, color = "royalblue4") +
    labs(title = paste0(title_prefix, " (b)")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  print(bubble_plot)
  
  # Filter samples with any reads
  filtered_samples = filtered_data %>%
    group_by(sample) %>%
    filter(any(!is.na(reads))) %>%
    ungroup() %>%
    pull(sample)
  
  filtered_data = filtered_data %>%
    filter(sample %in% filtered_samples)
  
  # Bubble plot with filtered data
  filtered_bubble_plot = ggplot(filtered_data, aes(x = sample, y = species, size = reads)) +
    geom_point(alpha = 0.4, color = "royalblue4") +
    labs(title = paste0(title_prefix, " (c)")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  print(filtered_bubble_plot)
  
  ggsave(filtered_bubble_plot, filename = paste0("Exploration plots/filtered_bubble_", title_prefix, ".png"))
  
}

# 
# # Taxonomic groups
exploration_plots(crop_long, 'Mammalia', 'Mammals')
exploration_plots(crop_long, 'Aves', 'Birds')
exploration_plots(crop_long, 'Amphibia', 'Amphibians')
exploration_plots(crop_long, 'Actinopteri', 'Fish')

# ______________________________________________________________________________
# extract month and sample location metadata

# Extract month and location from the sample names
crop_long = crop_long %>%
  mutate(month = substr(sample, 6, 7),
         location = substr(sample, 8, 9))

# Define function for bubble plot with species as an input
species_bubble_plot = function(data, target_species) {
  # Filter for the specified species
  species_data = data %>% 
    filter(species == target_species, reads > 0) %>%
    mutate(month = substr(sample, 6, 7),    # Extract month
           location = substr(sample, 8, 9)) # Extract location
  
  fb_eb_data = species_data %>%
    filter(location %in% c("FB", "EB")) %>%
    dplyr::select(location, species, reads)
  
  # Is the target species in the field or extraction blanks? 
  if (nrow(fb_eb_data) > 0) {
    print(fb_eb_data)
  } else {
    cat("Blanks are clear for", target_species, "\n")
  }
  
  # Create the bubble plot
  spec_plot = ggplot(species_data, aes(x = month, y = location, size = reads)) +
    geom_point(alpha = 0.6, color = "royalblue4") +
    labs(title = paste("Reads by Month and Location for", target_species),
         x = "Month", y = "Sampling Point", size = "Reads") +
    scale_x_discrete(limits = c(sprintf("%02d", 9:12), sprintf("%02d", 1:8), "13")) +
    scale_y_discrete(limits = c(sprintf("%02d", 1:14), "FB", "EB")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(spec_plot, filename = paste0("Exploration plots/by_month_", 
                                      target_species, ".png"))
  
  print(spec_plot)
}

species_bubble_plot(crop_long, "Castor_fiber")
species_bubble_plot(crop_long, "Lutra_lutra")


#_____________________________________________________________________________
# proportional read counts

# Calculate the total reads per sample
total_reads_per_sample = crop_long %>% 
  group_by(sample) %>% 
  summarise(total_reads = sum(reads))

# Join total reads back to the long format data
crop_long_prc = crop_long %>% 
  left_join(total_reads_per_sample, by = "sample")

# Calculate the proportion of reads for each species in each sample
crop_long_prc$proportion = crop_long_prc$reads / crop_long_prc$total_reads

# View the resulting data
head(crop_long_prc)

# Define function for bubble plot with species as an input
species_bubble_plot_prc = function(data, target_species) {
  # Filter for the specified species
  species_data = data %>% 
    filter(species == target_species, proportion > 0) %>%
    mutate(month = substr(sample, 6, 7),    # Extract month
           location = substr(sample, 8, 9)) # Extract location
  
  fb_eb_data = species_data %>%
    filter(location %in% c("FB", "EB")) %>%
    dplyr::select(location, species, proportion)
  
  # Is the target species in the field or extraction blanks? 
  if (nrow(fb_eb_data) > 0) {
    print(fb_eb_data)
  } else {
    cat("Blanks are clear for", target_species, "\n")
  }
  
  # Create the bubble plot
  spec_plot = ggplot(species_data, aes(x = month, y = location, size = proportion)) +
    geom_point(alpha = 0.6, color = "royalblue4") +
    labs(title = paste("Reads by Month and Location for", target_species),
         x = "Month", y = "Sampling Point", size = "Reads") +
    scale_x_discrete(limits = c(sprintf("%02d", 9:12), sprintf("%02d", 1:8), "13")) +
    scale_y_discrete(limits = c("02", "01", sprintf("%02d", 3:14), "FB", "EB")) +
  theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(spec_plot, filename = paste0("Exploration plots/prc_by_month_", 
                                      target_species, ".png"))
  
  print(spec_plot)
}

species_bubble_plot_prc(crop_long_prc, "Castor_fiber")
species_bubble_plot_prc(crop_long_prc, "Lutra_lutra")
species_bubble_plot_prc(crop_long_prc, "Sciurus_carolinensis")
species_bubble_plot_prc(crop_long_prc, "Ardea_cinerea")
species_bubble_plot_prc(crop_long_prc, "Capreolus_capreolus")
species_bubble_plot_prc(crop_long_prc, "Phasianus_colchicus")
species_bubble_plot_prc(crop_long_prc, "Meles_meles")
species_bubble_plot_prc(crop_long_prc, "Rana_temporaria")
species_bubble_plot_prc(crop_long_prc, "Lampetra_fluviatilis")
species_bubble_plot_prc(crop_long_prc, "Cottus_gobio")
species_bubble_plot_prc(crop_long_prc, "Anguilla_anguilla")


# Filter to mammals for mammal-only prc dataframe
crop_long_prc_noblanks <- crop_long_prc %>%
  filter(!grepl("EB|FB|POS|NEG", sample))

mammal_long <- crop_long_prc %>%
  filter(grepl("Mammalia", taxonomy, ignore.case = TRUE))

# Recalculate mammal-only total reads per sample
mammal_totals <- mammal_long %>%
  group_by(sample) %>%
  summarise(total_mammal_reads = sum(reads))

# Join totals back in
mammal_long <- mammal_long %>%
  left_join(mammal_totals, by = "sample")

# Calculate mammal-only proportion
mammal_long$prop_mammal <- mammal_long$reads / mammal_long$total_mammal_reads

# Clamp to (0,1) range for beta regression (optional but recommended)
mammal_long$prop_mammal_beta <- pmin(pmax(mammal_long$prop_mammal, 1e-6), 1 - 1e-6)

write.csv(mammal_long, "crop_prc_mams.csv", row.names = FALSE)
