"""
Figure 5 (spatial trait maps) rebuilt for print legibility and unit consistency.
Editor item 9: "scale numbers in Figure 5" too small at 18 cm print width.

Fixes vs. original notebook (figure_5_spatial_trait_maps.ipynb):
  - lat/lon tick + colorbar "scale numbers": 12/14 -> 22 pt (~8-9 pt at 18 cm)
  - colorbar labels 17 -> 20 pt; density-inset label/ticks 7/5 -> 14/12 pt
  - panel letters (a)-(i) moved clear of the column titles (no overlap)
  - LWC shown as g/cm2 (x10^-3), matching the Figure 4 histograms and Sec 3.2
    (trait-mean ~4.1 x10^-3 g/cm2); previously x10^-2.
Output: figures/figure5_remapcon_spatial_trait_maps_legible.{png,pdf}
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pathlib import Path

BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')
base_path = BASE / 'traits' / 'datasets' / 'clima_fit_prescribed_lai_ci'
pft_path = BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'
FIG_DIR = BASE / 'shift_dangermond_trait' / 'figures'

# ---- Load trait-based (spatially explicit) maps ----
chl_trait_ds = xr.open_dataset(base_path / 'chl_aviris_dangermond_clima_fit_reg.nc')
lma_trait_ds = xr.open_dataset(base_path / 'lma_aviris_dangermond_clima_fit_reg.nc')
lwc_trait_ds = xr.open_dataset(base_path / 'lwc_aviris_dangermond_clima_fit_reg.nc')

chl_trait = chl_trait_ds['chl'].mean(dim='time').values          # ug/cm2
lma_trait = lma_trait_ds['lma'].mean(dim='time').values * 1e4    # g/cm2 -> g/m2
lwc_trait_raw = lwc_trait_ds['lwc'].mean(dim='time').values
lwc_trait = lwc_trait_raw * 0.0018 if np.nanmax(lwc_trait_raw) > 1 else lwc_trait_raw

lats = chl_trait_ds['lat'].values
lons = chl_trait_ds['lon'].values

# ---- Load PFT-based (averaged) maps ----
chl_pft = xr.open_dataset(base_path / 'mean_masked_chl_aviris_dangermond_clima_fit.nc')['chl'].values.squeeze()
lma_pft = xr.open_dataset(base_path / 'mean_masked_lma_aviris_dangermond_clima_fit.nc')['lma'].values.squeeze()
lwc_pft = xr.open_dataset(base_path / 'mean_masked_lwc_aviris_dangermond_clima_fit.nc')['lwc'].values.squeeze() * 0.0018

pft_map = xr.open_dataset(pft_path)['Band1'].values
pft_mask = np.isin(pft_map, [2, 3, 4])

chl_diff = chl_trait - chl_pft
lma_diff = lma_trait - lma_pft
lwc_diff = lwc_trait - lwc_pft

msk = lambda a: np.where(pft_mask, a, np.nan)
traits_data = [
    ('chl', msk(chl_trait), msk(chl_pft), msk(chl_diff)),
    ('lma', msk(lma_trait), msk(lma_pft), msk(lma_diff)),
    ('lwc', msk(lwc_trait), msk(lwc_pft), msk(lwc_diff)),
]

trait_configs = {
    'chl': {'label': 'Chlorophyll Content', 'unit': 'µg/cm²', 'cmap': 'YlGn',
            'cmap_diff': 'RdBu_r', 'vmin': 0, 'vmax': 80, 'vmin_diff': -30, 'vmax_diff': 30},
    'lma': {'label': 'Leaf Mass per Area', 'unit': 'g/m²', 'cmap': 'YlOrBr',
            'cmap_diff': 'RdBu_r', 'vmin': 0, 'vmax': 200, 'vmin_diff': -80, 'vmax_diff': 80},
    'lwc': {'label': 'Leaf Water Content', 'unit': 'g/cm²', 'unit_raw': 'g/cm²',
            'unit_multiplier': 1000, 'cmap': 'Blues', 'cmap_diff': 'RdBu_r',
            'vmin': 0, 'vmax': 15, 'vmin_diff': -6, 'vmax_diff': 6},
}

# ---- Enlarged font sizes (the fix) ----
# Design canvas is 18 in wide; at the 18 cm print width fonts scale ~x0.39,
# so these sizes land at ~8-9 pt on the page (editor item 9: scale numbers legible).
TICK = 22        # lat/lon "scale numbers"  (was 12)
CBAR_TICK = 22   # colorbar scale numbers    (was 14)
CBAR_LBL = 20    # colorbar labels           (was 17)
STAT = 16
PANEL = 24
TITLE = 23
INSET_LBL = 14   # was 7
INSET_TICK = 12  # was 5

def add_density_inset(ax, data, config):
    inset_ax = inset_axes(ax, width="30%", height="30%", loc='lower right',
                          bbox_to_anchor=(0.02, 0.05, 1, 1), bbox_transform=ax.transAxes,
                          borderpad=1.0)
    v = data[np.isfinite(data) & (data != 0)]
    if len(v) > 0:
        inset_ax.hist(v, bins=50, color='gray', alpha=0.7, edgecolor='black',
                      linewidth=0.5, density=True)
        inset_ax.axvline(0, color='red', linestyle='--', linewidth=1.5, alpha=0.8)
        inset_ax.axvline(np.nanmean(v), color='blue', linestyle='--', linewidth=1.5, alpha=0.8)
        unit_label = config.get('unit_raw', config['unit']) if 'unit_multiplier' in config else config['unit']
        inset_ax.set_xlabel('Δ ' + unit_label, fontsize=INSET_LBL)
        inset_ax.tick_params(labelsize=INSET_TICK)
        inset_ax.grid(True, alpha=0.3, linewidth=0.5)

fig, axes = plt.subplots(3, 3, figsize=(18, 16), dpi=150)
column_titles = ['Trait-based Approach', 'PFT-based Approach', 'Difference (Trait - PFT)']
LON, LAT = np.meshgrid(lons, lats)
panel_labels = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']

for row_idx, (trait_name, trait_data, pft_data, diff_data) in enumerate(traits_data):
    config = trait_configs[trait_name]
    pidx = row_idx * 3
    cols = [
        (trait_data, config['cmap'], config['vmin'], config['vmax'], config['label'], False),
        (pft_data, config['cmap'], config['vmin'], config['vmax'], config['label'], False),
        (diff_data, config['cmap_diff'], config['vmin_diff'], config['vmax_diff'], config['label'], True),
    ]
    for col_idx, (data, cmap, vmin, vmax, label, is_diff) in enumerate(cols):
        ax = axes[row_idx, col_idx]
        disp = data * config.get('unit_multiplier', 1) if trait_name == 'lwc' else data
        im = ax.pcolormesh(LON, LAT, disp, cmap=cmap, vmin=vmin, vmax=vmax,
                           shading='auto', rasterized=True)
        ax.set_aspect('equal')
        ax.annotate(f'({panel_labels[pidx + col_idx]})', xy=(0, 1), xycoords='axes fraction',
                    xytext=(-4, 4), textcoords='offset points',
                    fontsize=PANEL, fontweight='bold', va='bottom', ha='right')
        stats_text = f'Mean: {np.nanmean(disp):.2f}\nSTD: {np.nanstd(disp):.2f}'
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=STAT,
                va='top', ha='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'))
        ax.tick_params(labelsize=TICK)
        ax.locator_params(axis='x', nbins=4)
        ax.locator_params(axis='y', nbins=4)
        cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02, fraction=0.046)
        pref = 'Δ ' if is_diff else ''
        cbar.set_label(f"{pref}{label} ({config['unit']})", fontsize=CBAR_LBL, fontweight='bold')
        if trait_name == 'lwc':
            cbar.ax.set_title('(×10⁻³)', fontsize=18, pad=10)
        cbar.ax.tick_params(labelsize=CBAR_TICK)
        if is_diff:
            add_density_inset(ax, disp, config)

for col_idx, col_title in enumerate(column_titles):
    axes[0, col_idx].set_title(col_title, fontsize=TITLE, fontweight='bold', pad=30)

plt.tight_layout(rect=[0.02, 0, 1, 0.99])
out_png = FIG_DIR / 'figure5_remapcon_spatial_trait_maps_legible.png'
out_pdf = FIG_DIR / 'figure5_remapcon_spatial_trait_maps_legible.pdf'
fig.savefig(out_png, dpi=300, facecolor='white', bbox_inches='tight')
fig.savefig(out_pdf, facecolor='white', bbox_inches='tight')
plt.close(fig)
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
