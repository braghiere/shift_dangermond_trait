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

