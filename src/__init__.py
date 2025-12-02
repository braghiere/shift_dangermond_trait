"""
AVIRIS Dangermond Trait Analysis - Shared Utilities

This package provides common data loading and plotting functions
used across multiple analysis notebooks.
"""

from .data_loading import (
    load_trait_dataset,
    load_pft_map,
    load_flux_data,
    load_pft_averaged_trait,
    get_date_info,
    extract_coordinates,
    BASE_DATA_DIR,
    TRAIT_DATA_DIR,
    PFT_MAP_PATH,
    FLUX_DATA_DIR
)

from .plotting_utils import (
    setup_publication_style,
    add_density_inset,
    add_pft_density_inset,
    add_statistics_box,
    add_panel_label,
    add_coastlines,
    create_trait_colorbar,
    get_trait_config,
    PFT_COLORS,
    PFT_NAMES
)

__version__ = '0.1.0'

__all__ = [
    # Data loading
    'load_trait_dataset',
    'load_pft_map',
    'load_flux_data',
    'load_pft_averaged_trait',
    'get_date_info',
    'extract_coordinates',
    'BASE_DATA_DIR',
    'TRAIT_DATA_DIR',
    'PFT_MAP_PATH',
    'FLUX_DATA_DIR',
    # Plotting
    'setup_publication_style',
    'add_density_inset',
    'add_pft_density_inset',
    'add_statistics_box',
    'add_panel_label',
    'add_coastlines',
    'create_trait_colorbar',
    'get_trait_config',
    'PFT_COLORS',
    'PFT_NAMES'
]
