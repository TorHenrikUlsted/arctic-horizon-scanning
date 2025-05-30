# Arctic Horizon Scanning: Potential Alien Vascular Plant Species

## Overview

This repository contains a generalizable analysis pipeline for horizon scanning of potential alien vascular plant species based on climatic niche modeling through hypervolume analysis. The framework identifies species from a candidate pool that could potentially establish in a reference region under current climate conditions.

While the provided example implementation focuses on Arctic plant invasions (analyzing non-Arctic species that could establish in the Arctic), the pipeline is designed to be adapted for other regions and invasion scenarios through customizable components.

**⚠️ Important Note**: The validation sequence and some visualization elements contain Arctic-specific assumptions that would need modification for other regions.

## Project Structure

### Core Pipeline (Region-Agnostic)
```
R/
├── initiate/          # Main workflow orchestration
├── utils/             # General utility functions
├── hypervolume/       # Core hypervolume analysis
├── filter/            # Filtering framework
├── setup/             # Setup framework
└── visualize/         # Visualization framework
```

### Customizable Components
```
R/
├── setup/custom_setup/
│   ├── region/              # YOUR REGION DEFINITION HERE
│   └── wrangle/
│       ├── known_species/   # YOUR REFERENCE SPECIES DATA
│       └── unknown_species/ # YOUR CANDIDATE SPECIES DATA
│
└── filter/custom_filters/
├── known_species/       # YOUR REFERENCE FILTERING LOGIC
└── unknown_species/     # YOUR CANDIDATE FILTERING LOGIC
```

### Example Implementation (Arctic)
```
example/
└── R/
├── setup/custom_setup/  # Arctic-specific setup
└── filter/              # Arctic-specific filters
```

## Installation

1. Clone the repository:
```bash
git clone https://github.com/TorHenrikUlsted/arctic-horizon-scanning.git
cd arctic-horizon-scanning
```

## Quick Start

### Option 1: Use Arctic Example
```r
# Uses the provided Arctic implementation
source("src/initiate/run.R")
```

This executes the complete pipeline:

- Setup sequence (data preparation)
- Filter sequence (species filtering)
- Hypervolume sequence (niche modeling)
- Visualization sequence (generating figures)
- Validation sequence (model validation)

### Option 2: Create Your Own Implementation

1. Copy templates to create your custom functions:

```bash
# Create your wrangling functions
cp src/setup/custom_setup/wrangle/unknown_species/wrangle_template.R \
   src/setup/custom_setup/wrangle/unknown_species/wrangle_yourdata.R

# Create your filter functions  
cp src/filter/custom_filters/unknown_species/unknown_template.R \
   src/filter/custom_filters/unknown_species/filter_yourdata.R
```

2. Edit ```config.yaml``` to use your functions:

```yaml
dataset:
  known: "yourknown"    # matches filter_yourknown.R
  unknown: "yourdata"   # matches filter_yourdata.R and wrangle_yourdata.R
```

3. Run your analysis:

```r
source("src/initiate/run.R")
```

### Adapting for Your Study

1. Define Your Study System
Decide on:

* **Reference region:** Where species might establish (e.g., Arctic, Mediterranean, specific country)
* **Candidate pool:** Source of potential invaders (e.g., global flora, neighboring regions)

2. Prepare Your Data
Region Definition
Edit ```src/setup/custom_setup/region/handle_region.R```:

```r
setup_region <- function() {
  # Load your shapefile
  region <- load_region("path/to/your/region.shp")
  
  # Process as needed
  # ...
  
  return(region_path)
}
```

### Species Data
Create wrangling functions following the template structure:

* Place reference species wranglers in wrangle/known_species/
* Place candidate species wranglers in wrangle/unknown_species/
* Each function should output species lists with scientificName column

### Filtering Logic
Create filter functions to:

* Combine multiple data sources
* Apply inclusion/exclusion criteria
* Separate present/absent species

3. Configure Analysis Parameters
Edit ```config.yaml```:

```yaml
simulation:
  projection: "laea"  # or appropriate for your region
  
dataset:
  known: "your_reference_dataset"
  unknown: "your_candidate_dataset"

visualization:
  shape.name: "your-region-shapefile"
  region.name: "Your Region Name"
```
  
4. Select Climate Variables
Choose appropriate bioclimatic variables for your region:

```yaml
hypervolume:
  dimensions: [1, 2, 3, 12]  # BioClim indices 1-19
```

## Requirements

### System Requirements
- R (≥ 4.0.0)
- Minimum 16GB RAM (32GB+ recommended for full analysis)
- ~50GB free disk space for data processing
- Multi-core processor recommended for parallel processing

### Required R Packages
The project automatically installs and loads required packages including:

- Core: `data.table`, `terra`, `hypervolume`, `parallel`
- Data access: `rgbif`, `geodata`, `WorldFlora`
- Visualization: `ggplot2`, `tidyterra`, `gridExtra`
- Analysis: `gamlss`, `mgcv`, `vegan`

See `src/utils/utils.R` for the complete package list.

### Running Individual Components

You can also run specific sequences:

```r
# Data setup only
source("src/setup/setup_sequence.R")
setup_sequence(coord.uncertainty = 1000, hv.method = "box", ...)

# Species filtering
source("src/filter/filter_sequence.R")
filter_sequence(spec.known = "arctic", spec.unknown = "glonaf", ...)

# Hypervolume analysis
source("src/hypervolume/parallel_hypervolume.R")
hypervolume_sequence(spec.list = species_files, ...)

# Visualization
source("src/visualize/visualize.R")
visualize_sequence(res.unknown = "glonaf", res.known = "arctic", ...)
```

### Using Example Data

The repository includes example data for testing:

```r
# Run with test data
main(
  spec.known = "test_known",
  spec.unknown = "test_small",  # or "test_big"
  # ... other parameters
)
```

## Outputs

The analysis generates several types of outputs in the `outputs/` directory:

### Data Outputs
- **filter/**: Filtered species lists (present/absent)
- **hypervolume/**: Hypervolume analysis results and projections
- **setup/**: Processed climate and species data

### Visualizations
- **Figure 1**: Species richness frequency distributions
- **Figure 2**: Potential invasion hotspots
- **Figure 3**: Area of occupancy and climatic suitability
- **Figure 4**: Taxonomic composition patterns
- **Figure 5**: Geographic connections and origin regions
- **Figure 6-8**: Latitudinal distribution patterns
- **Figure 9**: Sankey diagram of species flows

### Validation Results
- Model performance metrics
- Comparison between known and predicted Arctic species

### Custom Visualizations

Modify visualization parameters in `src/visualize/components/visualize_figures.R`

## Data Sources

- **Species occurrences**: GBIF (Global Biodiversity Information Facility)
- **Climate data**: WorldClim 2.1 or CHELSA
- **Taxonomic backbone**: World Flora Online (WFO)
- **Arctic regions**: Circumpolar Arctic Vegetation Map (CAVM)

## Limitations

Validation sequence: Currently Arctic-specific, comparing against known Arctic flora
Some visualizations: Assume Northern Hemisphere and invasion "into" reference region
Geographic projections: Default to Arctic-centric projections

## Citation

If you use this code, please cite using the "cite this repository" button in the project's about menu.

## Contact

For questions or issues, please open an issue on GitHub or contact:

- Tor Henrik Ulsted (ORCID: [0000-0001-8854-2696](https://orcid.org/0000-0001-8854-2696))

## Acknowledgments

This project was developed as part of a master thesis at NTNU, investigating biological invasions in the Arctic under current and future climate scenarios.