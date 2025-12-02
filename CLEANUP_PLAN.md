# Repository Cleanup & Improvement Plan
## shift_dangermond_trait

**Date:** December 2, 2025  
**Current Status:** 29 notebooks, multiple analysis scripts, ~33MB of figures

---

## 🎯 Priority Actions

### 1. **Notebook Organization & Archival**

#### A. **Keep - Publication-Ready Notebooks** (Move to `/notebooks/final/`)
```
✅ figure_3_combined_traits_temporal.ipynb (8.3MB) - Figure 3
✅ figure_5_spatial_trait_maps.ipynb (4.4MB) - Figure 5
✅ figure_6_updated_clean.ipynb (2.2MB) - Figure 6 (clean version)
✅ figures_S1_S2_final.ipynb - Supplementary figures
✅ topography_trait_correlation_full_updated.ipynb (2.5MB) - Topography analysis
```

#### B. **Archive - Deprecated/Older Versions** (Move to `/notebooks/archive/`)
```
📦 figure_6_updated.ipynb - Superseded by _clean version
📦 compare_sif_all_lr_week.ipynb (2.1MB) - Old version
📦 compare_sif_all_lr_week_prescribed_lai_ci_res_005.ipynb (2.3MB) - Experiment
📦 plot_gpp_difference.ipynb - Superseded by clima_fit version
📦 pft_trait_map_all.ipynb - Old version without prescribed lai_ci
📦 prepare_fitted_traits_prescribed_lai.ipynb - Old version
📦 plt_chl_gifs.ipynb - Old version
📦 plt_lma_gifs.ipynb - Old version
```

#### C. **Keep - Active Analysis Notebooks** (Organize by category)
```
/notebooks/validation/
  - compare_sif_all_hr_week.ipynb
  - compare_sif_all_lr_week_prescribed_lai_ci.ipynb
  - plot_gpp_sif_fit_prescribed_lai_ci.ipynb (1.2MB)

/notebooks/trait_analysis/
  - pft_trait_map_all_fit_prescribed_lai_ci.ipynb
  - plt_chl_gifs_fit_lai_ci.ipynb
  - plt_lma_gifs_fit_lai_ci.ipynb
  - plt_lwc_gifs_fit_lai_ci.ipynb (1.4MB)
  - plsr_trait_histograms.ipynb
  - plot_trait.ipynb

/notebooks/processing/
  - prepare_fitted_traits_prescribed_lai_ci.ipynb
  - plot_cumulative_gpp.ipynb
  - plot_gpp_difference_clima_fit_prescribed_lai_ci.ipynb

/notebooks/exploratory/
  - plot_hyperspectral_curves.ipynb (4.0MB)
  - plt_lai_gifs.ipynb
  - compare_sif_all_hr.ipynb
```

---

### 2. **Figure Management**

#### A. **Organize Final Figures** (Create `/figures/publication/`)
```
✅ figure_S1_chl_temporal_revised.png (2.3MB)
✅ figure_S1_chl_temporal_600dpi.png (3.1MB)
✅ figure_S2_lma_temporal_revised.png (2.1MB)
✅ figure_S2_lma_temporal_600dpi.png (2.9MB)
✅ figure_6_gpp_sif_spatial_temporal_revised.png (1.5MB)
✅ figure_6_gpp_sif_spatial_temporal_600dpi.png (2.2MB)
✅ topography_trait_correlation_heatmap.png
✅ topography_trait_scatterplots.png (6.5MB)
```

#### B. **Archive Old/Intermediate Figures** (Create `/figures/archive/`)
```
📦 figure_6_gpp_sif_spatial_revised.png - Superseded
📦 figure_6_gpp_sif_spatial_600dpi.png - Superseded
📦 figure_6_complete_with_tropomi.png - Experimental
📦 figure_6_with_tropomi_*.png - Experimental
📦 *_hist.png files (18 files, histograms)
```

#### C. **Keep Working Figures** (`/figures/working/`)
```
- topography_and_traits_spatial.png
- topography_vs_traits_comprehensive.png
```

---

### 3. **Code Organization**

#### A. **Create Module Structure**
```
/src/
  __init__.py
  data_loading.py      # Consolidate data loading functions
  plotting_utils.py    # Common plotting functions (density insets, etc.)
  spatial_analysis.py  # Spatial analysis utilities
  trait_processing.py  # Trait calculation and processing
```

#### B. **Consolidate Analysis Scripts**
```
/analysis/
  sif_validation/
    - compare_sif_*.py scripts
  trait_extraction/
    - extract_*.py scripts
  pft_analysis/
    - pft_trait_map_all.py
```

---

### 4. **Documentation Improvements**

#### A. **Update README.md**
- Add section on final manuscript figures
- Document notebook organization
- Add requirements.txt or environment.yml
- Include data provenance and sources
- Add citation information

#### B. **Create Additional Documentation**
```
CONTRIBUTING.md - Guidelines for adding new analyses
FIGURES.md - Description of all final figures
DATA_PROCESSING.md - Data pipeline documentation
REPRODUCIBILITY.md - Steps to reproduce results
```

#### C. **Add Notebook Metadata**
Each publication notebook should have:
- Clear title and purpose
- Author and date
- Data sources
- Output files generated
- Reviewer comments addressed

---

### 5. **Data Management**

#### A. **Document Data Dependencies**
Create `DATA_SOURCES.md`:
```markdown
- Trait maps: `/traits/datasets/clima_fit_prescribed_lai_ci/`
- PFT classifications: `/California_Vegetation_WHRTYPE_Dangermond/`
- Flux data: `/fitting/fitted_prescribed_lai_ci/`
- Topography: [document source]
```

#### B. **Create Data Manifest**
Track file sizes, dates, and checksums for key datasets

---

### 6. **Version Control Hygiene**

#### A. **Update .gitignore**
```
# Jupyter checkpoints
.ipynb_checkpoints/
**/.ipynb_checkpoints/

# Python cache
__pycache__/
*.pyc
*.pyo

# Large data files (if not already tracked)
*.nc
*.tif
*.hdf

# OS files
.DS_Store
Thumbs.db

# Temporary files
*~
*.swp
*.swo

# Output that regenerates
/figures/archive/
/notebooks/archive/
```

#### B. **Clean Git History**
- Consider git-lfs for large notebooks if needed
- Document branching strategy (main, dev, review-responses)

---

### 7. **Code Quality Improvements**

#### A. **Standardize Imports**
Create common import blocks across notebooks:
```python
# Standard imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path

# Project imports
from src.plotting_utils import add_density_inset_pft
from src.data_loading import load_trait_data, load_pft_map
```

#### B. **Extract Repeated Code**
Common patterns across notebooks:
- Data loading boilerplate
- Figure styling (plt.rcParams)
- Path definitions
- PFT masking operations
- Statistical calculations

#### C. **Add Error Handling**
```python
# Example
try:
    ds = xr.open_dataset(file_path)
except FileNotFoundError:
    print(f"Error: {file_path} not found")
    print("Expected location: ...")
```

---

### 8. **Testing & Validation**

#### A. **Create Test Suite**
```
/tests/
  test_data_loading.py
  test_plotting.py
  test_spatial_analysis.py
```

#### B. **Add Data Validation**
- Check coordinate ranges
- Verify PFT classifications
- Validate time dimensions
- Check for NaN patterns

---

## 📋 Implementation Checklist

### Phase 1: Organization (Week 1)
- [ ] Create new directory structure
- [ ] Move publication notebooks to `/notebooks/final/`
- [ ] Move deprecated notebooks to `/notebooks/archive/`
- [ ] Organize active notebooks by category
- [ ] Move final figures to `/figures/publication/`
- [ ] Archive old figures

### Phase 2: Code Cleanup (Week 2)
- [ ] Create `/src/` module structure
- [ ] Extract common functions to modules
- [ ] Update notebooks to use shared functions
- [ ] Standardize imports across notebooks
- [ ] Add error handling

### Phase 3: Documentation (Week 2-3)
- [ ] Update README.md
- [ ] Create FIGURES.md
- [ ] Create DATA_SOURCES.md
- [ ] Add notebook headers
- [ ] Document reviewer responses

### Phase 4: Version Control (Week 3)
- [ ] Update .gitignore
- [ ] Clean up git history if needed
- [ ] Tag final publication version
- [ ] Create release notes

### Phase 5: Testing (Week 4)
- [ ] Add data validation checks
- [ ] Create reproducibility test
- [ ] Verify all publication figures regenerate correctly

---

## 🗂️ Proposed Directory Structure

```
shift_dangermond_trait/
├── README.md (updated)
├── CONTRIBUTING.md (new)
├── FIGURES.md (new)
├── DATA_SOURCES.md (new)
├── REPRODUCIBILITY.md (new)
├── LICENSE
├── requirements.txt (new)
├── .gitignore (updated)
│
├── src/                           # Shared Python modules
│   ├── __init__.py
│   ├── data_loading.py
│   ├── plotting_utils.py
│   ├── spatial_analysis.py
│   └── trait_processing.py
│
├── notebooks/
│   ├── final/                     # Publication-ready notebooks
│   │   ├── figure_3_combined_traits_temporal.ipynb
│   │   ├── figure_5_spatial_trait_maps.ipynb
│   │   ├── figure_6_updated_clean.ipynb
│   │   ├── figures_S1_S2_final.ipynb
│   │   └── topography_trait_correlation_full_updated.ipynb
│   │
│   ├── validation/                # SIF/GPP validation analyses
│   │   ├── compare_sif_all_hr_week.ipynb
│   │   ├── compare_sif_all_lr_week_prescribed_lai_ci.ipynb
│   │   └── plot_gpp_sif_fit_prescribed_lai_ci.ipynb
│   │
│   ├── trait_analysis/            # Trait mapping and GIFs
│   │   ├── pft_trait_map_all_fit_prescribed_lai_ci.ipynb
│   │   ├── plt_chl_gifs_fit_lai_ci.ipynb
│   │   ├── plt_lma_gifs_fit_lai_ci.ipynb
│   │   ├── plt_lwc_gifs_fit_lai_ci.ipynb
│   │   └── plsr_trait_histograms.ipynb
│   │
│   ├── processing/                # Data preparation
│   │   ├── prepare_fitted_traits_prescribed_lai_ci.ipynb
│   │   ├── plot_cumulative_gpp.ipynb
│   │   └── plot_gpp_difference_clima_fit_prescribed_lai_ci.ipynb
│   │
│   ├── exploratory/               # Exploratory analyses
│   │   ├── plot_hyperspectral_curves.ipynb
│   │   ├── plt_lai_gifs.ipynb
│   │   └── compare_sif_all_hr.ipynb
│   │
│   └── archive/                   # Old/superseded notebooks
│       ├── figure_6_updated.ipynb
│       ├── compare_sif_all_lr_week.ipynb
│       └── [other deprecated notebooks]
│
├── figures/
│   ├── publication/               # Final manuscript figures
│   │   ├── figure_S1_chl_temporal_revised.png
│   │   ├── figure_S2_lma_temporal_revised.png
│   │   ├── figure_6_gpp_sif_spatial_temporal_revised.png
│   │   └── [high-res versions]
│   │
│   ├── working/                   # Work-in-progress figures
│   │   └── topography_and_traits_spatial.png
│   │
│   └── archive/                   # Old/superseded figures
│       └── [deprecated figures]
│
├── analysis/                      # Analysis scripts
│   ├── sif_validation/
│   │   ├── compare_sif_mean_all_trait_hr_week.py
│   │   └── compare_sif_all_trait_hr_week.py
│   │
│   ├── trait_extraction/
│   │   ├── extract_lai.py
│   │   ├── extract_lma.py
│   │   └── extract_all_vars_v1.py
│   │
│   └── pft_analysis/
│       └── pft_trait_map_all.py
│
├── data_processing/               # Data preparation scripts
│   ├── clip_dangermond_updated.py
│   ├── regrid_netcdf.py
│   └── spatial_ref_dangermond.py
│
├── plotting/                      # Plotting utilities
│   └── plot_trait_pdf.py
│
├── shell_scripts/                 # Automation scripts
│   └── [bash scripts]
│
├── julia_scripts/                 # Julia analyses
│   └── [julia scripts]
│
├── tests/                         # Unit tests (new)
│   ├── test_data_loading.py
│   └── test_plotting.py
│
└── docs/                          # Additional documentation (new)
    ├── notebook_guide.md
    └── data_pipeline.md
```

---

## 💡 Best Practices Going Forward

### For New Notebooks:
1. **Clear naming convention**: `[figure|analysis]_[topic]_[version].ipynb`
2. **Header cell** with title, author, date, purpose
3. **Import shared functions** from `/src/`
4. **Document data sources** and paths
5. **Include reproducibility info** (environment, versions)

### For Figures:
1. **Semantic naming**: `figure_[N]_[description]_[resolution].png`
2. **Keep both 300dpi and 600dpi** versions
3. **Include PDF versions** for vector graphics
4. **Archive superseded versions** rather than deleting

### For Code:
1. **Extract repeated code** to shared modules
2. **Add docstrings** to all functions
3. **Use type hints** where helpful
4. **Add error handling** for file operations

### For Version Control:
1. **Meaningful commit messages** referencing reviewer comments
2. **Branch per major change** (e.g., `rb/rv-02` for review round 2)
3. **Tag releases** for publication milestones
4. **Keep clean history** by archiving rather than deleting

---

## 📊 Current Statistics

- **Total notebooks**: 29
- **Publication notebooks**: 5 (17%)
- **Total figure size**: ~33MB
- **Largest notebook**: 8.3MB (figure_3)
- **Deprecated notebooks**: ~10 (35%)
- **Analysis scripts**: ~15 Python files

## 🎯 Target Statistics

- **Active notebooks**: ~20 (organized)
- **Archived notebooks**: ~10
- **Shared modules**: 4-5 files
- **Documentation files**: 5-6
- **Test coverage**: >70% of core functions

---

## 🚀 Quick Start Commands

```bash
# Create new directory structure
mkdir -p notebooks/{final,validation,trait_analysis,processing,exploratory,archive}
mkdir -p figures/{publication,working,archive}
mkdir -p analysis/{sif_validation,trait_extraction,pft_analysis}
mkdir -p src tests docs

# Move publication notebooks
mv notebooks/figure_*.ipynb notebooks/final/
mv notebooks/figures_S*.ipynb notebooks/final/
mv notebooks/topography_trait_correlation_full_updated.ipynb notebooks/final/

# Move final figures
mv figures/figure_S*_revised.png figures/publication/
mv figures/figure_S*_600dpi.png figures/publication/
mv figures/figure_6_gpp_sif_spatial_temporal*.png figures/publication/

# Archive old versions
mv figures/figure_6_gpp_sif_spatial_revised.png figures/archive/
mv figures/*_hist*.png figures/archive/

# Create initialization files
touch src/__init__.py
touch tests/__init__.py

# Update git
git add .gitignore
git status
```

---

## 📝 Notes

- Keep `analysis_summary.txt` and `topography_trait_correlations.csv` in root for quick reference
- Consider creating a Jupyter Book for comprehensive documentation
- May want to add pre-commit hooks for code quality
- Consider adding a CITATION.cff file for proper citation

**Priority Level**: HIGH - Many duplicates and unclear organization  
**Estimated Time**: 2-3 weeks for full cleanup  
**Impact**: Significant improvement in maintainability and reproducibility
