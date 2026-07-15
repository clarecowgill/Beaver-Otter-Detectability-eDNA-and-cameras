# Environmental DNA and Camera Trap Monitoring of Beaver and Otter Activity

## Overview

This repository contains all data and code required to reproduce analyses from:

**"Environmental DNA and camera traps capture different ecological signals of beaver and otter spatial and temporal activity"**

This study integrates **riverine environmental DNA (12S metabarcoding)** with **year-round camera trapping** to investigate the spatial and temporal ecology of **Eurasian beaver (*Castor fiber*)** and **Eurasian otter (*Lutra lutra*)**.

The repository contains the complete analytical workflow used to:

* process raw eDNA metabarcoding data
* compare biodiversity detected by eDNA and camera traps
* quantify temporal and spatial patterns in species detections
* characterise behavioural activity from camera trap observations
* develop spatially and temporally weighted detectability metrics
* investigate environmental and behavioural drivers of eDNA detection using mixed-effects modelling.

---

## FAIR eDNA Metadata

- **`12S_Cropton_FAIRe_checklist.xlsx`**

  FAIR eDNA metadata record for the Cropton environmental DNA dataset, completed using the FAIR eDNA reporting standard. The checklist documents sample collection, laboratory protocols, bioinformatic processing, quality assurance and dataset provenance to facilitate transparency, reproducibility and long-term data reuse.

---

# 1. Data Files

## Raw Data

### eDNA

* **`Cropton_blast98_denoise.tsv`**
  Raw denoised 12S metabarcoding read table containing taxonomic assignments for all samples.

* **`blank_association.csv`**
  Reference file linking blank samples to field samples for species-specific limit of detection (LOD) calculations.

---

### Camera Trap Data

* **`cropton_cams_behaviour.csv`**
  Camera trap detections including species identifications, timestamps, behavioural observations, counts and deployment information.

---

### Metadata

* **`crop_meta.csv`**
  Sample metadata including sampling dates, environmental variables (temperature, flow, pH, conductivity and filtered volume), deployment information and camera metadata.

* **`cam_locations.csv`**
  Spatial information for each camera location, including distance from the primary beaver lodge.

* **`distance_matrix.csv`**
  Along-stream distance matrix describing distances between camera locations and eDNA sampling sites, used to calculate spatially weighted detectability scores.

---

## Processed Data

These datasets are generated during the processing workflow and are provided to facilitate reproducibility without requiring users to rerun the complete filtering pipeline.

* **`cropton_edna_filtered.csv`**
  Filtered eDNA read table following contamination removal, species-specific limit of detection filtering, proportional read thresholding and removal of unassigned taxa.

* **`cropton_edna_proportions.csv`**
  Sample-by-species proportional read table derived from the filtered eDNA dataset.

* **`combined_edna_camera_data.csv`**
  Combined eDNA and camera trap dataset containing harmonised species names, proportional abundances and taxonomic classifications used for diversity and method comparison analyses.

---

# 2. Reference Database

Database used for taxonomic identification of UK vertebrates from 12S metabarcoding.

* **`12S_verts.fasta`**
  Curated reference database of 12S sequences for UK vertebrate taxa.

* **`12S_verts_tax_map.txt`**
  Mapping file linking reference sequences to their corresponding GenBank taxonomic IDs.

These files were used during BLAST-based taxonomic assignment of vertebrate sequences following bioinformatic processing.

---

# 3. R Scripts

### `cropton_initial.R`

Processes the raw 12S metabarcoding dataset.

This script:

* separates Cropton samples from other study datasets
* calculates species-specific limits of detection from blank samples
* applies contamination filtering
* applies a 0.1% proportional read threshold
* removes unassigned taxa and non-target species
* generates filtered and proportional read datasets
* produces exploratory visualisations of species detections.

This script generates the processed eDNA datasets required for all downstream analyses.

---

### `overall_diversity.R`

Compares biodiversity detected by environmental DNA and camera traps.

Analyses include:

* species richness
* Shannon and Simpson diversity
* community composition
* taxonomic overlap
* detection probabilities
* species accumulation curves
* rarefaction analyses
* community dissimilarity
* comparisons between monitoring methods.

---

### `Heatmaps and GAM.R`

Investigates spatial and temporal detection dynamics of beavers and otters.

This script:

* generates monthly heatmaps of camera detections and eDNA proportional reads
* fits Generalised Additive Models (GAMs) describing temporal detection probability
* compares alternative model structures
* estimates site-level effects
* visualises temporal trends
* performs model diagnostics.

---

### `cams_behaviour.R`

Processes behavioural observations recorded from camera traps.

Calculates:

* monthly behavioural budgets
* weekly behavioural budgets
* site-level behavioural composition
* water-associated behaviour
* behavioural diversity metrics
* temporal persistence
* spatial autocorrelation.

This script also generates behavioural predictor variables used in subsequent detectability modelling.

---

### `detectability_models.R`

Constructs camera-derived detectability metrics and investigates factors influencing eDNA detection.

Analyses include:

* spatially and temporally weighted detectability scores
* camera activity metrics
* behavioural predictor calculation
* integration of environmental variables
* correlation and multicollinearity analyses
* Generalised Linear Mixed Models (GLMMs) for otter detection
* Bayesian Generalised Linear Mixed Models for beaver detection
* model comparison and validation
* forest plot visualisation of predictor effects
* analyses relating prey diversity to otter activity and feeding behaviour.

---

# 4. Bioinformatic Processing

Raw 12S metabarcoding data were processed using the **Tapirs** metabarcoding pipeline, including:

* quality filtering
* paired-end read merging
* dereplication
* chimera removal
* VSEARCH clustering
* BLAST-based taxonomic assignment.

Taxonomic identification was performed using the curated UK vertebrate reference database included within this repository.

Within this repository, additional filtering includes:

* species-specific limits of detection calculated from blank samples
* removal of human and domestic taxa
* removal of unassigned sequences
* proportional read filtering (0.1% of total reads per sample).

---

# 5. Reproducibility

To reproduce all analyses:

1. Run:

   * `cropton_initial.R`

   This script generates the processed eDNA datasets required for downstream analyses.

2. Run the remaining scripts independently:

   * `overall_diversity.R`
   * `Heatmaps and GAM.R`
   * `cams_behaviour.R`
   * `detectability_models.R`

---

# Repository Structure

```text
12S_Cropton_FAIRe_checklist.xlsx

Data Files/
├── Raw Data/
│   ├── Cropton_v2_blast98_denoise.tsv
│   ├── blank_association.csv
│   ├── cropton_cams_behaviour.csv
│   ├── crop_meta.csv
│   ├── cam_locations.csv
│   ├── distance_matrix.csv
│   └── otter_prey_eDNA.csv
│
├── Processed Data/
│   ├── cropton_edna_filtered.csv
│   ├── cropton_edna_proportions.csv
│   └── combined_edna_camera_data.csv
│
Reference Database/
├── 12S_verts.fasta
└── 12S_verts_tax_map.txt

R Scripts/
├── cropton_initial.R
├── overall_diversity.R
├── Heatmaps and GAM.R
├── cams_behaviour.R
└── detectability_models.R
```
