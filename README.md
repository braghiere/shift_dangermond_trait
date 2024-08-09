# SHIFT Dangermond Traits vs. PFT Project

This repository contains scripts and Jupyter notebooks utilized for data processing, analysis, and visualization in the SHIFT Dangermond Traits vs. PFT project. The repository is structured to help navigate through various phases of data manipulation, from raw data processing to detailed visualizations and analysis.

## Repository Structure

### `data_processing/`
Scripts dedicated to initial data manipulation, including clipping and regridding operations. Essential for preparing the data for further analysis.

- `clip_dangermond.py`: Clips data to a specified region.
- `regrid_netcdf.py`: Regrids data from one spatial resolution to another.
- `clip_dangermond_updated.py`: Updated version of the initial clipping script.
- `regrid_script.py`: Another script for regridding data with additional functionality.
- `spatial_ref_dangermond.py`: Sets spatial references for datasets.

### `analysis/`
Contains scripts that perform complex data analysis tasks, including statistical comparisons and trait extraction.

- Various `compare_*.py` scripts: Perform statistical comparisons across different data metrics like solar-induced fluorescence and plant functional traits.
- `extract_*.py` scripts: Extract specific traits or variables from datasets.
- `pft_trait_map_all.py` and `pft_trait_map.py`: Map and analyze plant functional traits.

### `plotting/`
Scripts for generating plots and visual representations of data.

- `plot_trait_pdf.py`: Generates PDFs of plant traits plots.

### `julia_scripts/`
Julia language scripts for specific analytical tasks, providing efficient computation and processing capabilities.

- Scripts like `1_gpp.jl` and `shift_spectra.jl`: Focus on gross primary productivity calculations and spectral data analysis.

### `shell_scripts/`
Shell scripts to automate routine tasks such as data merging and regridding.

- Scripts such as `clip_bash.sh` and `regrid_files.sh`: Help in automating the workflow and data management.

### `notebooks/`
Jupyter notebooks that provide interactive environments to explore and visualize data dynamically.

- Notebooks like `plot_hyperspectral_curves.ipynb` and `plt_chl_gifs.ipynb`: Offer detailed and interactive visual analysis of the data.

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

