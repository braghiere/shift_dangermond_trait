"""
Plotting utilities for AVIRIS Dangermond trait analysis.

This module provides standardized plotting functions for creating
publication-quality figures with consistent styling.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colorbar import Colorbar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from typing import Dict, Tuple, Optional


# PFT color scheme (consistent across all figures)
PFT_COLORS = {
    2: '#1f77b4',  # Blue - Scrub
    3: '#2ca02c',  # Green - Chaparral
    4: '#ff7f0e'   # Orange - Woodland
}

PFT_NAMES = {
    2: 'Scrub',
    3: 'Chaparral',
    4: 'Woodland'
}


def setup_publication_style(scale: float = 1.2):
    """
    Set matplotlib rcParams for publication-quality figures.
    
    Parameters
    ----------
    scale : float
        Scaling factor for font sizes (default: 1.2)
    """
    plt.rcParams['font.size'] = 14 * scale
    plt.rcParams['axes.labelsize'] = 16 * scale
    plt.rcParams['axes.titlesize'] = 18 * scale
    plt.rcParams['xtick.labelsize'] = 14 * scale
    plt.rcParams['ytick.labelsize'] = 14 * scale
    plt.rcParams['legend.fontsize'] = 12 * scale


def add_density_inset(ax, data: np.ndarray, config: Dict, 
                     position: str = 'lower right',
                     width: str = "28%", height: str = "28%") -> plt.Axes:
    """
    Add a density distribution inset to an axes.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to add inset to
    data : np.ndarray
        Array of difference values
    config : dict
        Configuration dictionary with 'unit' key
    position : str
        Position of inset ('lower right', 'upper right', etc.)
    width : str
        Width of inset as percentage
    height : str
        Height of inset as percentage
    
    Returns
    -------
    matplotlib.axes.Axes
        The created inset axes
    """
    # Create inset axes
    if position == 'lower right':
        inset_ax = inset_axes(ax, width=width, height=height, loc='lower right',
                             bbox_to_anchor=(0.02, 0.05, 1, 1), 
                             bbox_transform=ax.transAxes,
                             borderpad=1.0)
    else:
        inset_ax = inset_axes(ax, width=width, height=height, 
                             loc=position, borderpad=1.0)
    
    # Remove NaNs and zeros
    valid_data = data[np.isfinite(data) & (data != 0)]
    
    if len(valid_data) > 0:
        # Plot histogram
        inset_ax.hist(valid_data, bins=50, color='gray', alpha=0.7, 
                     edgecolor='black', linewidth=0.5, density=True)
        
        # Add vertical line at zero (red dashed)
        inset_ax.axvline(0, color='red', linestyle='--', 
                        linewidth=1.5, alpha=0.8)
        
        # Add mean line (blue dashed)
        mean_val = np.nanmean(valid_data)
        inset_ax.axvline(mean_val, color='blue', linestyle='--', 
                        linewidth=1.5, alpha=0.8)
        
        # Style the inset
        unit_label = config.get('unit', '')
        inset_ax.set_xlabel(f'Δ {unit_label}', fontsize=7)
        inset_ax.set_ylabel('', fontsize=7)
        inset_ax.tick_params(labelsize=5)
        inset_ax.grid(True, alpha=0.3, linewidth=0.5)
    
    return inset_ax


def add_pft_density_inset(ax, data: np.ndarray, pft_map: np.ndarray,
                          config: Dict, position: str = 'lower right',
                          width: str = "28%", height: str = "28%") -> plt.Axes:
    """
    Add a PFT-separated density distribution inset.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to add inset to
    data : np.ndarray
        Array of values to plot
    pft_map : np.ndarray
        Array of PFT classifications
    config : dict
        Configuration dictionary with 'unit' key
    position : str
        Position of inset
    width : str
        Width of inset as percentage
    height : str
        Height of inset as percentage
    
    Returns
    -------
    matplotlib.axes.Axes
        The created inset axes
    """
    # Create inset axes
    if position == 'lower right':
        inset_ax = inset_axes(ax, width=width, height=height, loc='lower right',
                             bbox_to_anchor=(0.02, 0.05, 1, 1), 
                             bbox_transform=ax.transAxes,
                             borderpad=1.0)
    else:
        inset_ax = inset_axes(ax, width=width, height=height, 
                             loc=position, borderpad=1.0)
    
    # Plot separate histograms for each PFT
    for pft_id in [2, 3, 4]:
        pft_data = data[pft_map == pft_id]
        valid_data = pft_data[np.isfinite(pft_data) & (pft_data != 0)]
        
        if len(valid_data) > 0:
            inset_ax.hist(valid_data, bins=30, color=PFT_COLORS[pft_id],
                         alpha=0.5, edgecolor='black', linewidth=0.5,
                         density=True, label=PFT_NAMES[pft_id])
    
    # Add reference lines
    inset_ax.axvline(0, color='red', linestyle='--', linewidth=1.5, alpha=0.8)
    
    all_valid = data[np.isfinite(data) & (data != 0)]
    if len(all_valid) > 0:
        mean_val = np.nanmean(all_valid)
        inset_ax.axvline(mean_val, color='black', linestyle='--', 
                        linewidth=1.5, alpha=0.8)
    
    # Style
    unit_label = config.get('unit', '')
    inset_ax.set_xlabel(f'Δ {unit_label}', fontsize=7)
    inset_ax.set_ylabel('', fontsize=7)
    inset_ax.tick_params(labelsize=5)
    inset_ax.legend(fontsize=5, framealpha=0.9)
    inset_ax.grid(True, alpha=0.3, linewidth=0.5)
    
    return inset_ax


def add_statistics_box(ax, data: np.ndarray, position: Tuple[float, float] = (0.98, 0.98),
                       ha: str = 'right', va: str = 'top',
                       fontsize: int = 16) -> None:
    """
    Add a statistics box showing mean and standard deviation.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to add statistics box to
    data : np.ndarray
        Data array for computing statistics
    position : tuple
        (x, y) position in axes coordinates
    ha : str
        Horizontal alignment
    va : str
        Vertical alignment
    fontsize : int
        Font size for text
    """
    mean_val = np.nanmean(data)
    std_val = np.nanstd(data)
    
    stats_text = f'Mean: {mean_val:.2f}\nSTD: {std_val:.2f}'
    
    ax.text(position[0], position[1], stats_text, 
           transform=ax.transAxes,
           fontsize=fontsize, va=va, ha=ha,
           bbox=dict(boxstyle='round', facecolor='white', 
                    alpha=0.9, edgecolor='black'))


def add_panel_label(ax, label: str, position: Tuple[float, float] = (-0.05, 1.05),
                   fontsize: int = 22) -> None:
    """
    Add a panel label (e.g., '(a)', '(b)') to an axes.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to add label to
    label : str
        Label text (e.g., 'a', 'b', 'c')
    position : tuple
        (x, y) position in axes coordinates
    fontsize : int
        Font size for label
    """
    ax.text(position[0], position[1], f'({label})', 
           transform=ax.transAxes,
           fontsize=fontsize, fontweight='bold', 
           va='bottom', ha='right')


def add_coastlines(ax, lons: np.ndarray, lats: np.ndarray,
                  color: str = 'k', linewidth: float = 1.5, 
                  alpha: float = 0.7) -> None:
    """
    Add coastline boundary to a map plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to add coastlines to
    lons : np.ndarray
        Longitude values (1D array)
    lats : np.ndarray
        Latitude values (1D array)
    color : str
        Line color
    linewidth : float
        Line width
    alpha : float
        Line transparency
    """
    # Draw rectangle around study area
    ax.plot(lons[[0, -1, -1, 0, 0]], 
           lats[[0, 0, -1, -1, 0]], 
           color=color, linewidth=linewidth, alpha=alpha)


def create_trait_colorbar(im, ax, label: str, unit: str,
                          orientation: str = 'vertical',
                          fontsize_label: int = 17,
                          fontsize_ticks: int = 14) -> Colorbar:
    """
    Create a standardized colorbar for trait maps.
    
    Parameters
    ----------
    im : matplotlib.image.AxesImage
        Image object from pcolormesh
    ax : matplotlib.axes.Axes
        Axes associated with the colorbar
    label : str
        Colorbar label (trait name)
    unit : str
        Unit of measurement
    orientation : str
        'vertical' or 'horizontal'
    fontsize_label : int
        Font size for label
    fontsize_ticks : int
        Font size for tick labels
    
    Returns
    -------
    matplotlib.colorbar.Colorbar
        The created colorbar
    """
    cbar = plt.colorbar(im, ax=ax, orientation=orientation, 
                       pad=0.02, fraction=0.046)
    cbar.set_label(f"{label} ({unit})", fontsize=fontsize_label, 
                   fontweight='bold')
    cbar.ax.tick_params(labelsize=fontsize_ticks)
    
    return cbar


def get_trait_config(trait_name: str) -> Dict:
    """
    Get standard configuration for a trait.
    
    Parameters
    ----------
    trait_name : str
        Name of trait ('chl', 'lma', 'lwc', 'gpp', 'sif')
    
    Returns
    -------
    dict
        Configuration dictionary with label, unit, cmap, vmin, vmax
    """
    configs = {
        'chl': {
            'label': 'Chlorophyll Content',
            'unit': 'µg/cm²',
            'cmap': 'YlGn',
            'cmap_diff': 'RdBu_r',
            'vmin': 0,
            'vmax': 80,
            'vmin_diff': -30,
            'vmax_diff': 30
        },
        'lma': {
            'label': 'Leaf Mass per Area',
            'unit': 'g/m²',
            'cmap': 'YlOrBr',
            'cmap_diff': 'RdBu_r',
            'vmin': 0,
            'vmax': 200,
            'vmin_diff': -80,
            'vmax_diff': 80
        },
        'lwc': {
            'label': 'Leaf Water Content',
            'unit': 'g/cm²',
            'cmap': 'Blues',
            'cmap_diff': 'RdBu_r',
            'vmin': 0,
            'vmax': 0.03,
            'vmin_diff': -0.01,
            'vmax_diff': 0.01
        },
        'gpp': {
            'label': 'GPP',
            'unit': 'µmol/m²/s',
            'cmap': 'YlGn',
            'cmap_diff': 'RdBu_r',
            'vmin': 0,
            'vmax': 30,
            'vmin_diff': -10,
            'vmax_diff': 10
        },
        'sif': {
            'label': 'SIF',
            'unit': 'mW/m²/sr/nm',
            'cmap': 'plasma',
            'cmap_diff': 'RdBu_r',
            'vmin': 0,
            'vmax': 2,
            'vmin_diff': -1,
            'vmax_diff': 1
        }
    }
    
    if trait_name.lower() not in configs:
        raise ValueError(f"Unknown trait: {trait_name}. "
                        f"Available: {list(configs.keys())}")
    
    return configs[trait_name.lower()]
