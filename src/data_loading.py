"""
Data loading utilities for AVIRIS Dangermond trait analysis.

This module provides standardized functions for loading NetCDF datasets,
PFT maps, and flux data used across multiple notebooks.
"""

import xarray as xr
import numpy as np
from pathlib import Path
from typing import Tuple, Dict, Optional


# Standard base paths
BASE_DATA_DIR = Path('/home/renatob/data/FluoData1/aviris_dangermond')
TRAIT_DATA_DIR = BASE_DATA_DIR / 'traits/datasets/clima_fit_prescribed_lai_ci'
PFT_MAP_PATH = BASE_DATA_DIR / 'California_Vegetation_WHRTYPE_Dangermond/output_latlon.nc'
FLUX_DATA_DIR = BASE_DATA_DIR / 'fitting'


def load_trait_dataset(trait_name: str, data_dir: Optional[Path] = None) -> xr.Dataset:
    """
    Load a trait dataset from NetCDF file.
    
    Parameters
    ----------
    trait_name : str
        Name of the trait ('chl', 'lma', 'lwc', 'cab', 'car')
    data_dir : Path, optional
        Directory containing trait datasets. Defaults to TRAIT_DATA_DIR.
    
    Returns
    -------
    xr.Dataset
        Loaded trait dataset with coordinates and time dimension
    """
    if data_dir is None:
        data_dir = TRAIT_DATA_DIR
    
    file_path = data_dir / f'{trait_name}_aviris_dangermond_clima_fit_reg.nc'
    
    if not file_path.exists():
        raise FileNotFoundError(f"Trait file not found: {file_path}")
    
    return xr.open_dataset(file_path)


def load_pft_map(pft_path: Optional[Path] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load PFT (Plant Functional Type) map and create valid PFT mask.
    
    Parameters
    ----------
    pft_path : Path, optional
        Path to PFT NetCDF file. Defaults to PFT_MAP_PATH.
    
    Returns
    -------
    pft_map : np.ndarray
        2D array of PFT classifications
    pft_mask : np.ndarray
        Boolean mask for valid PFTs (2, 3, 4)
    """
    if pft_path is None:
        pft_path = PFT_MAP_PATH
    
    if not pft_path.exists():
        raise FileNotFoundError(f"PFT map not found: {pft_path}")
    
    pft_ds = xr.open_dataset(pft_path)
    pft_map = pft_ds['Band1'].values
    
    # Valid PFTs for this study: 2, 3, 4
    pft_mask = np.isin(pft_map, [2, 3, 4])
    
    return pft_map, pft_mask


def load_flux_data(day_idx: int, flux_type: str = 'trait', 
                   flux_dir: Optional[Path] = None) -> xr.Dataset:
    """
    Load flux data (GPP, SIF) for a specific day.
    
    Parameters
    ----------
    day_idx : int
        Day index (0-12 for the 13 AVIRIS flights)
    flux_type : str
        Type of flux data: 'trait' for trait-based or 'pft' for PFT-based
    flux_dir : Path, optional
        Directory containing flux files. Defaults to FLUX_DATA_DIR.
    
    Returns
    -------
    xr.Dataset
        Flux dataset with GPP and SIF variables
    """
    if flux_dir is None:
        flux_dir = FLUX_DATA_DIR
    
    if flux_type == 'trait':
        subdir = 'fitted_prescribed_lai_ci'
        filename = f'shift_fluxes_day_{day_idx:02d}_clima_fit_reg_jmax.nc'
    elif flux_type == 'pft':
        subdir = ''
        filename = f'pft_shift_fluxes_day_{day_idx:02d}_clima_fit_reg_jmax.nc'
    else:
        raise ValueError(f"Invalid flux_type: {flux_type}. Must be 'trait' or 'pft'")
    
    file_path = flux_dir / subdir / filename if subdir else flux_dir / filename
    
    if not file_path.exists():
        raise FileNotFoundError(f"Flux file not found: {file_path}")
    
    return xr.open_dataset(file_path, decode_times=False)


def load_pft_averaged_trait(trait_name: str, data_dir: Optional[Path] = None) -> xr.Dataset:
    """
    Load PFT-averaged (masked) trait dataset.
    
    Parameters
    ----------
    trait_name : str
        Name of the trait ('chl', 'lma', 'lwc')
    data_dir : Path, optional
        Directory containing trait datasets. Defaults to TRAIT_DATA_DIR.
    
    Returns
    -------
    xr.Dataset
        PFT-averaged trait dataset
    """
    if data_dir is None:
        data_dir = TRAIT_DATA_DIR
    
    file_path = data_dir / f'mean_masked_{trait_name}_aviris_dangermond_clima_fit.nc'
    
    if not file_path.exists():
        raise FileNotFoundError(f"PFT-averaged trait file not found: {file_path}")
    
    return xr.open_dataset(file_path)


def get_date_info() -> Dict[int, str]:
    """
    Get mapping of day indices to flight dates.
    
    Returns
    -------
    dict
        Dictionary mapping day index (0-12) to date string (YYYY-MM-DD)
    """
    return {
        0: '2022-02-24',
        1: '2022-03-04',
        2: '2022-03-12',
        3: '2022-03-18',
        4: '2022-03-26',
        5: '2022-04-01',
        6: '2022-04-09',
        7: '2022-04-15',
        8: '2022-04-23',
        9: '2022-04-29',
        10: '2022-05-07',
        11: '2022-05-13',
        12: '2022-05-29'
    }


def extract_coordinates(dataset: xr.Dataset) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract latitude and longitude coordinates from a dataset.
    
    Parameters
    ----------
    dataset : xr.Dataset
        Dataset containing 'lat' and 'lon' coordinates
    
    Returns
    -------
    lats : np.ndarray
        Latitude values
    lons : np.ndarray
        Longitude values
    """
    lats = dataset['lat'].values
    lons = dataset['lon'].values
    return lats, lons
