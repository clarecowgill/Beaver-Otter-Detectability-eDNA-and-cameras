library(ggplot2)
library(dplyr)
library(stringr)
library(vegan)
library(lubridate)
library(tidyr)
library(RColorBrewer)
library(VennDiagram)
library(grid)
library(gridExtra)
library(tibble)
library(pairwiseAdonis)
library(iNEXT)
library(ggrepel)

## Clean eDNA data --------------------------
edna_raw = read.csv('crop_dat.csv', row.names = 1, header = T)

edna_raw$species = rownames(edna_raw)
edna_long = reshape(edna_raw,
                    varying = list(names(edna_raw)[!names(edna_raw) %in% c("taxonomy", "species")]),
                    v.names = "reads", timevar = "sample",
                    times = names(edna_raw)[!names(edna_raw) %in% c("taxonomy", "species")],
                    idvar = c("species", "taxonomy"), direction = "long")
rownames(edna_long) = NULL

edna_long <- edna_long %>%
  filter(!species %in% c('Homo_sapiens', 'Bos_taurus', 'Capra_hircus',
                         'Ovis_aries', 'Sus_scrofa', 'Canis_lupus', 'Equus_caballus',
                         'Astatotilapia_calliptera'))

edna_long$class <- str_extract(edna_long$taxonomy, "c__([A-Za-z_]+)") %>%
  str_remove("c__")
edna_long$genus_clean <- str_extract(edna_long$taxonomy, "g__([A-Za-z_]+)") %>%
  str_remove("g__")

edna_long_filtered <- edna_long %>%
  filter(!is.na(genus_clean))

edna_long_filtered$clean_species <- edna_long_filtered$species
edna_long_filtered$clean_species[!grepl("_", edna_long_filtered$species)] <-
  edna_long_filtered$genus_clean[!grepl("_", edna_long_filtered$species)]

# keep counts, then normalise per sample
edna <- edna_long_filtered %>%
  select(sample, reads, class, species = clean_species, genus = genus_clean) %>%
  rename(count = reads) %>%
  mutate(method = "eDNA",
         is_full_species = grepl("_", species)) %>%
  group_by(sample) %>%
  mutate(prop_count = count / sum(count, na.rm = TRUE)) %>%
  ungroup()

## Clean camera data and join--------------------------
cams_raw <- read.csv('cropton_cams_behaviour.csv')

cams_long_raw <- cams_raw %>%
  select(sample, species, date) %>%
  mutate(date = as.POSIXct(date, format = "%d/%m/%Y %H:%M"))

cams_final <- cams_long_raw %>%
  group_by(sample, species) %>%
  arrange(date) %>%
  mutate(is_new_visit = c(TRUE, diff(date, units = "mins") > 30)) %>%
  mutate(visit_id = cumsum(is_new_visit)) %>%
  summarise(count = n_distinct(visit_id), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(prop_count = count / sum(count, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(method = "camera",
         is_full_species = grepl("_", species),
         class = NA,
         genus = word(species, 1, sep = "_"))

## Combine dataframes --------------------------
combined_data <- bind_rows(edna, cams_final) %>%
  filter(count > 0)

combined_data_aligned <- combined_data %>%
  group_by(sample, genus) %>%
  mutate(
    aligned_species = case_when(
      genus %in% c("Turdus", "Apodemus", "Anas", "Fringilla", "Columba", "Myotis") ~ genus,
      any(is_full_species) ~ first(species[is_full_species]),
      TRUE ~ first(species)
    )
  ) %>%
  ungroup() %>%
  select(sample, method, aligned_species, class, count, prop_count)

# write csv and manually add in classes for camera entries
write.csv(combined_data_aligned, 'combined_with_props.csv', row.names = FALSE)

## Diversity estimates --------------------------
combined <- read.csv('combined_with_props_update.csv', header = T)

combined <- combined %>%
  filter(!grepl("POS|NEG|EB|FB", sample, ignore.case = TRUE)) %>%
  filter(!grepl("11$|12$|13$|14$", sample)) %>%
  mutate(
    class2 = ifelse(class %in% c("Actinopteri","Hyperoartia", "Amphibia"), "Fish_Amphib", class),
    month = substr(sample, 6, 7),
    site  = substr(sample, 8, 9)
  )

#-------------------------------
# Alpha diversity (richness, Shannon, Simpson)
#-------------------------------
excluded_species <- c("Aegithalos_caudatus", "Anas", "Apodemus", "Columba", 
                      "Erithacus_rubecula", "Fringilla", "Motacilla_cinerea", 
                      "Phasianus_colchicus", "Rattus_norvegicus", "Sylvia_atricapilla", 
                      "Troglodytes_troglodytes", "Turdus")

# Venn diagrams of camera and eDNA data ----
#-------------------------------
# 1) Overall Venn for all species
#-------------------------------
all_species <- combined %>%
  group_by(method) %>%
  summarise(species = list(unique(aligned_species)), .groups = "drop") %>%
  deframe()

edna_species   <- all_species$eDNA
camera_species <- all_species$camera

grid.newpage()
draw.pairwise.venn(
  area1 = length(edna_species),
  area2 = length(camera_species),
  cross.area = length(intersect(edna_species, camera_species)),
  category = c("eDNA", "Camera"),
  fill = c("#00BFC4", "#F8766D"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.2,
  main = "Overall species: eDNA vs Camera"
)

#-------------------------------
# 2) Overall Venn for birds + mammals
#-------------------------------
bm_data <- combined %>% filter(class2 %in% c("Aves", "Mammalia"))

edna_species   <- unique(bm_data$aligned_species[bm_data$method == "eDNA"])
camera_species <- unique(bm_data$aligned_species[bm_data$method == "camera"])

grid.newpage()
draw.pairwise.venn(
  area1 = length(edna_species),
  area2 = length(camera_species),
  cross.area = length(intersect(edna_species, camera_species)),
  category = c("eDNA", "Camera"),
  fill = c("#00BFC4", "#F8766D"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.2,
  main = "Overall Birds + Mammals: eDNA vs Camera"
)

# Create presence/absence (1 = detected, 0 = not detected) for each method
pa_df <- combined %>%
  select(aligned_species, class2, method) %>%
  distinct() %>%  # only unique species-method combinations
  mutate(presence = 1) %>%
  pivot_wider(names_from = method, values_from = presence, values_fill = 0) %>%
  rename(class = class2) %>%
  arrange(aligned_species)

# Inspect
head(pa_df)


#---------------------------------
# Detection probabilities

# Define constants
faulty_samples <- c("CROP_1001", "CROP_1103", "CROP_0108", "CROP_0208", "CROP_1303",
                    "CROP_0901", "CROP_0903", "CROP_0908", "CROP_0909", 
                    "CROP_0910", "CROP_0801")

problem_sites <- c("02", "04", "05", "06", "07")
problem_taxa <- c("Aegithalos_caudatus", "Anas", "Apodemus", 
                  "Columba", "Erithacus_rubecula", "Fringilla", 
                  "Motacilla_cinerea", "Phasianus_colchicus", 
                  "Rattus_norvegicus", "Sylvia_atricapilla",
                  "Troglodytes_troglodytes", "Turdus")

# Filter valid samples: exclude blanks, faulty, controls etc.
valid_samples <- combined %>%
  filter(
    !grepl("POS|NEG|EB|FB", sample, ignore.case = TRUE),
    !grepl("11$|12$|13$|14$", sample),
    !sample %in% faulty_samples
  ) %>%
  distinct(sample, method) %>%
  mutate(site = str_sub(sample, -2, -1))  # Extract site from sample ID

# Define species list
species_list <- unique(combined$aligned_species)

# Calculate eDNA detection probs (no site exclusions needed)
edna_samples <- valid_samples %>% filter(method == "eDNA") %>% pull(sample)
edna_detections <- combined %>%
  filter(method == "eDNA", sample %in% edna_samples) %>%
  group_by(species = aligned_species, sample) %>%
  summarise(detected = any(count > 0), .groups = "drop")
edna_grid <- expand_grid(species = unique(edna_detections$species), sample = edna_samples)
edna_full <- edna_grid %>%
  left_join(edna_detections, by = c("species", "sample")) %>%
  mutate(detected = replace_na(detected, FALSE))
edna_prob <- edna_full %>%
  group_by(species) %>%
  summarise(
    total_samples = n(),
    samples_detected = sum(detected),
    eDNA_detection_prob = samples_detected / total_samples,
    .groups = "drop"
  )

# Calculate camera detection probs with problem taxa-site exclusions
camera_samples <- valid_samples %>% filter(method == "camera")
camera_grid <- expand_grid(species = species_list, sample = camera_samples$sample) %>%
  left_join(camera_samples, by = "sample")  # to add site info

# Exclude samples from problem_sites for problem_taxa in camera_grid
camera_grid_filtered <- camera_grid %>%
  filter(!(species %in% problem_taxa & site %in% problem_sites))

camera_detections <- combined %>%
  filter(method == "camera", sample %in% camera_samples$sample) %>%
  group_by(species = aligned_species, sample) %>%
  summarise(detected = any(count > 0), .groups = "drop")

camera_full <- camera_grid_filtered %>%
  left_join(camera_detections, by = c("species", "sample")) %>%
  mutate(detected = replace_na(detected, FALSE))

camera_prob <- camera_full %>%
  group_by(species) %>%
  summarise(
    total_samples = n(),
    samples_detected = sum(detected),
    camera_detection_prob = samples_detected / total_samples,
    .groups = "drop"
  )

# Combine detection probabilities and join taxon class
detection_comparison <- full_join(
  edna_prob %>% select(species, eDNA_detection_prob),
  camera_prob %>% select(species, camera_detection_prob),
  by = "species"
) %>%
  replace_na(list(eDNA_detection_prob = 0, camera_detection_prob = 0)) %>%
  filter(eDNA_detection_prob > 0 & camera_detection_prob > 0) %>%
  left_join(
    combined %>% distinct(aligned_species, class) %>% rename(species = aligned_species),
    by = "species"
  ) %>%
  mutate(
    species_label = case_when(
      species == "Apodemus" ~ "A. sylvaticus",
      str_detect(species, "_") ~ paste0(
        str_sub(species, 1, 1), ". ",
        str_replace(species, "^[^_]+_", "")
      ),
      TRUE ~ paste0(species, " spp.")
    )
  )

# Example manual color values; adjust class names and colors to your data
class_colors <- c(
  "Mammalia" = "#f46d43",
  "Aves" = "#66c2a5",
  "Amphibia" = "#053061"
)

# Plot
comparison_plot <- ggplot(detection_comparison, aes(x = eDNA_detection_prob, y = camera_detection_prob, color = class)) +
  geom_point(size = 1, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_text_repel(
    aes(label = species_label),
    size = 3,
    fontface = "bold.italic",
    max.overlaps = 20,
    box.padding = 0.35,
    point.padding = 0.3,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_color_manual(values = class_colors) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none"
  ) +
  labs(
    x = "eDNA Detection Probability",
    y = "Camera Detection Probability",
    color = "Taxonomic Class"
  )

print(comparison_plot)

ggsave("comparison_plot.png", comparison_plot, dpi = 300)
