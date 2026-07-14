"""
Figure 5 (spatial trait maps) rebuilt for print legibility, a tighter layout, and
unit consistency. Editor item 9: "scale numbers in Figure 5" too small at 18 cm.

Layout / fixes vs. the original notebook (figure_5_spatial_trait_maps.ipynb):
  - lat/lon tick + colorbar "scale numbers": 12/14 -> 22 pt (~8-9 pt at 18 cm)
  - shared axes: latitude labels only on column 1, longitude only on the bottom row
  - ONE shared colorbar for the Trait + PFT columns (identical scale); the
    Difference column keeps its own colorbar
  - constrained_layout removes the inter-row white space and enlarges the maps
  - density insets on the Difference maps sit fully inside (their ticks preserved)
  - panel letters (a)-(i) inside the top-left corner, lowercase, matching Fig 4
  - LWC shown as g/cm2 (x10^-3), matching the Fig 4 histograms and Sec 3.2
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
CBAR_LBL = 18    # colorbar labels           (was 17)
STAT = 16
PANEL = 24
TITLE = 23
INSET_LBL = 14   # was 7
INSET_TICK = 12  # was 5

def add_density_inset(ax, data, config):
    # fully inside the lower-right corner, opaque white background, own ticks kept
    inset_ax = inset_axes(ax, width="40%", height="32%", loc='lower right', borderpad=1.1)
    v = data[np.isfinite(data) & (data != 0)]
    if len(v) == 0:
        return
    inset_ax.hist(v, bins=50, color='0.55', alpha=0.95, edgecolor='black',
                  linewidth=0.4, density=True)
    inset_ax.axvline(0, color='red', linestyle='--', linewidth=1.4)
    inset_ax.axvline(np.nanmean(v), color='blue', linestyle='--', linewidth=1.4)
    unit_label = config.get('unit_raw', config['unit'])
    if 'unit_multiplier' in config:
        unit_label += ' (×10⁻³)'
    inset_ax.set_xlabel('Δ ' + unit_label, fontsize=INSET_LBL, labelpad=1.5)
    inset_ax.tick_params(labelsize=INSET_TICK, length=2.5, pad=1.5)
    inset_ax.locator_params(axis='x', nbins=3)
    inset_ax.locator_params(axis='y', nbins=3)
    inset_ax.patch.set_facecolor('white')
    inset_ax.patch.set_alpha(0.9)
    for s in inset_ax.spines.values():
        s.set_linewidth(0.6)

LON, LAT = np.meshgrid(lons, lats)
panel_labels = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
column_titles = ['Trait-based Approach', 'PFT-based Approach', 'Difference (Trait - PFT)']

fig, axes = plt.subplots(3, 3, figsize=(18.0, 12.0), dpi=150, constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=0.10, h_pad=0.02, wspace=0.02, hspace=0.015)

for row_idx, (trait_name, trait_data, pft_data, diff_data) in enumerate(traits_data):
    config = trait_configs[trait_name]
    mult = config.get('unit_multiplier', 1) if trait_name == 'lwc' else 1
    panels = [
        (trait_data, config['cmap'], config['vmin'], config['vmax'], False),
        (pft_data,   config['cmap'], config['vmin'], config['vmax'], False),
        (diff_data,  config['cmap_diff'], config['vmin_diff'], config['vmax_diff'], True),
    ]
    ims = []
    for col_idx, (data, cmap, vmin, vmax, is_diff) in enumerate(panels):
        ax = axes[row_idx, col_idx]
        disp = data * mult
        im = ax.pcolormesh(LON, LAT, disp, cmap=cmap, vmin=vmin, vmax=vmax,
                           shading='auto', rasterized=True)
        ims.append(im)
        ax.set_aspect('equal')
        # panel letter inside the top-left corner (lowercase, consistent with Fig 4)
        ax.text(0.035, 0.965, f'({panel_labels[row_idx * 3 + col_idx]})', transform=ax.transAxes,
                fontsize=PANEL, fontweight='bold', va='top', ha='left',
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.85, edgecolor='none'))
        ax.text(0.965, 0.965, f'Mean: {np.nanmean(disp):.2f}\nSTD: {np.nanstd(disp):.2f}',
                transform=ax.transAxes, fontsize=STAT, va='top', ha='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'))
        ax.tick_params(labelsize=TICK)
        ax.locator_params(axis='x', nbins=3)
        ax.locator_params(axis='y', nbins=4)
        if col_idx != 0:                       # latitude labels only on the first column
            ax.tick_params(labelleft=False)
        if row_idx != 2:                       # longitude labels only on the bottom row
            ax.tick_params(labelbottom=False)
        if row_idx == 0:
            ax.set_title(column_titles[col_idx], fontsize=TITLE, fontweight='bold', pad=8)
        if is_diff:
            add_density_inset(ax, disp, config)

    # one shared colorbar for the Trait + PFT columns (identical colour scale)
    cb = fig.colorbar(ims[1], ax=[axes[row_idx, 0], axes[row_idx, 1]], fraction=0.046, pad=0.02)
    cb.set_label(f"{config['label']} ({config['unit']})", fontsize=CBAR_LBL, fontweight='bold')
    cb.ax.tick_params(labelsize=CBAR_TICK)
    # separate colorbar for the Difference column
    cbd = fig.colorbar(ims[2], ax=axes[row_idx, 2], fraction=0.046, pad=0.02)
    cbd.set_label(f"Δ {config['label']} ({config['unit']})", fontsize=CBAR_LBL, fontweight='bold')
    cbd.ax.tick_params(labelsize=CBAR_TICK)
    if trait_name == 'lwc':
        cb.ax.set_title('(×10⁻³)', fontsize=16, pad=6)
        cbd.ax.set_title('(×10⁻³)', fontsize=16, pad=6)

out_png = FIG_DIR / 'figure5_remapcon_spatial_trait_maps_legible.png'
out_pdf = FIG_DIR / 'figure5_remapcon_spatial_trait_maps_legible.pdf'
fig.savefig(out_png, dpi=300, facecolor='white')
fig.savefig(out_pdf, facecolor='white')
plt.close(fig)
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
