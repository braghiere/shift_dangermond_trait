# SHIFT Dangermond Traits vs. PFT Project

This repository contains scripts and Jupyter notebooks for data processing, analysis, and visualization comparing trait-based vs. PFT-based approaches in the SHIFT Dangermond campaign (February-May 2022).

## Repository Structure

### 📁 `notebooks/`
Jupyter notebooks organized by function:

#### `final/` - Publication Figures
- `figure_3_combined_traits_temporal.ipynb` - Figure 3: Temporal dynamics of leaf traits
- `figure_5_spatial_trait_maps.ipynb` - Figure 5: Spatial distribution of traits
- `figure_6_updated_clean.ipynb` - Figure 6: GPP and SIF comparisons with density insets
- `figures_S1_S2_final.ipynb` - Supplementary Figures S1/S2: Chl and LMA temporal comparisons
- `topography_trait_correlation_full_updated.ipynb` - Topography-trait relationships

#### `validation/` - Model Validation
- `compare_sif_all_lr_week_prescribed_lai_ci.ipynb` - SIF validation against TROPOMI
- `plot_gpp_sif_fit_prescribed_lai_ci.ipynb` - GPP/SIF temporal validation

#### `trait_analysis/` - Trait Analysis
- `pft_trait_map_all_fit_prescribed_lai_ci.ipynb` - PFT-based trait mapping
- `plot_trait.ipynb` - Individual trait visualization
- `plsr_trait_histograms.ipynb` - PLSR trait distribution analysis

#### `processing/` - Data Processing
- `prepare_fitted_traits_prescribed_lai_ci.ipynb` - Prepare fitted trait datasets
- `plot_cumulative_gpp.ipynb` - Cumulative GPP calculations
- `plot_gpp_difference_clima_fit_prescribed_lai_ci.ipynb` - GPP difference maps

#### `exploratory/` - Exploratory Analysis
- `plot_hyperspectral_curves.ipynb` - Hyperspectral reflectance curves
- `plt_*_gifs*.ipynb` - Animated GIFs of trait evolution

#### `archive/` - Deprecated/Old Versions
- Historical notebooks and superseded analyses

### 📁 `figures/`
Generated figures organized into:
- `publication/` - Final publication-ready figures (300 & 600 DPI)
- `working/` - Intermediate figures from active analysis
- `archive/` - Old figure versions

### 📁 `src/`
Shared Python utilities (NEW):
- `data_loading.py` - Standard functions for loading NetCDF datasets, PFT maps, and flux data
- `plotting_utils.py` - Publication-quality plotting functions with consistent styling
- `__init__.py` - Package initialization

### 📁 `data_processing/`
Data preparation scripts:
- `clip_dangermond.py`, `clip_dangermond_updated.py` - Clip data to study region
- `regrid_netcdf.py`, `regrid_script.py` - Regrid to consistent resolution
- `spatial_ref_dangermond.py` - Set spatial references

### 📁 `analysis/`
Analysis scripts:
- `compare_*.py` - Statistical comparisons (SIF, GPP, traits)
- `extract_*.py` - Extract specific variables
- `pft_trait_map*.py` - PFT-trait mapping

### 📁 `plotting/`
Plotting scripts:
- `plot_trait_pdf.py` - Generate trait PDFs

### 📁 `julia_scripts/`
Julia scripts for CliMA Land model:
- `1_gpp.jl`, `shift_spectra.jl` - GPP calculations and spectral analysis

### 📁 `shell_scripts/`
Automation scripts:
- `clip_bash.sh`, `regrid_files.sh` - Batch processing

## Key Datasets

**13 AVIRIS-NG Flight Dates** (2022):
- Day 0: Feb 24 (campaign start)
- Days 1-11: Mar 4 - May 13
- Day 12: May 29 (campaign end)

**Plant Functional Types (PFTs)**:
- PFT 2: Scrub (blue)
- PFT 3: Chaparral (green)
- PFT 4: Woodland (orange)

**Analyzed Traits**:
- Chlorophyll content (Chl) - µg/cm²
- Leaf mass per area (LMA) - g/m²
- Leaf water content (LWC) - g/cm²
- Carotenoid content (Car) - µg/cm²

**Model Outputs**:
- Gross Primary Productivity (GPP) - µmol/m²/s
- Solar-Induced Fluorescence (SIF) - mW/m²/sr/nm

## Usage

Clone the repository and navigate to the desired directory to run specific scripts or notebooks:

```bash
git clone https://github.com/braghiere/shift_dangermond_trait.git
cd shift_dangermond_trait
```
## Data Availability 

The datasets used in this study are accessible at the Caltech Data Library under https://doi.org/10.22002/z15fj-44h89, including the SHIFT AVIRIS-NG dataset for the Dangermond Preserve at 30m resolution, which provides detailed spectral reflectance information across 399 wavelengths. This also includes TROPOMI sun-induced fluorescence data at 740 nm, CliMA Land model outputs of ecosystem fluxes, and trait variations for different plant functional types. These datasets cover the period from February to May 2022.

The Level 2A (L2A) unrectified surface reflectance images from NASA's Airborne Visible / Infrared Imaging Spectrometer-Next Generation (AVIRIS-NG) instrument, collected as part of the Surface Biology and Geology High-Frequency Time Series (SHIFT) campaign during February to May 2022, including reflected radiance at 5-nm intervals in the Visible to Shortwave Infrared (VSWIR) spectral range from 380-2510 nm, can be found at https://doi.org/10.3334/ORNLDAAC/2183 (Brodrick et al., 2023).

The California Multi-Source Vegetation Layer, specifically depicting Wildlife Habitat Relationship classes (WHRTYPE) (California Department of Forestry and Fire Protection, 2014), is available at https://gis.data.ca.gov/maps/CALFIRE-Forestry::california-vegetation-whrtype/about.

The CliMA Land model used to simulate surface hyperspectral reflectance and transmittance, energy, and carbon fluxes is available at https://github.com/CliMA/Land.

The codes used in the analysis presented in this paper are available on GitHub and can be accessed at https://github.com/braghiere/shift_dangermond_trait/.


## References

Brodrick, P., R. Pavlick, M. Bernas, J.W. Chapman, R. Eckert, M. Helmlinger, M. Hess-Flores, L.M. Rios, F.D. Schneider, M.M. Smyth, M. Eastwood, R.O. Green, D.R. Thompson, K.D. Chadwick, & D.S. Schimel. (2023). SHIFT: AVIRIS-NG L2A Unrectified Reflectance. ORNL DAAC. https://doi.org/10.3334/ORNLDAAC/2183

Chadwick, K. D., Davis, F., Miner, K. R., Pavlick, R., Reynolds, M., Townsend, P. A., Brodrick, P. G., Ade, C., Allen, J., Anderegg, L., Angel, Y., Boving, I., Byrd, K. B., Campbell, P., Carberry, L., Cavanaugh, K. C., Cavanaugh, K. C., Easterday, K., Eckert, R., … Schimel, D. (2024). Unlocking Ecological Insights from Subseasonal Visible-to-Shortwave Infrared Imaging Spectroscopy: The SHIFT Campaign. Ecosphere.

Queally, N., Davis, F. W., Chadwick, K. D., Ade, C., Anderegg, L., Angel, Y., Baker, B., Boving, I., Braghiere, R. K., Brodrick, P., Campbell, P., Cryer, J., Cushman, K. C., Dao, P. D., Dibartolo, A., Eckert, R., Grant, K., Heberlein, B., Johnson, M., … Schimel, D. S. (2024). SHIFT: Vegetation Plot Characterization, Santa Barbara County, CA, 2022. ORNL Distributed Active Archive Center. https://doi.org/10.3334/ORNLDAAC/2295

Nature Conservancy (2022). Lidar Survey of Dangermond Preserve, CA. OpenTopography. https://doi.org/10.5069/G9T43R8K

