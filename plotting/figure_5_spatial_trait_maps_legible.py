"""
Figure 5 (spatial trait maps) SPLIT into two print-legible figures to avoid the
crowding of a single 3x3 grid (editor item 9: scale numbers legible at 18 cm).

  figure5a_trait_vs_pft : Trait vs PFT, 2 columns x 3 traits (CHL/LMA/LWC rows),
                          one shared colorbar per trait row.
  figure5b_difference   : the three Trait-minus-PFT maps as a 1x3 ROW, each with
                          its own colorbar and a Δ-distribution inset.

Design choices:
  - DATA: conservative-remapping (remapcon) retrievals — the same dataset as Fig 4
    and Sec 3.2 — reproducing the manuscript Fig 5 values exactly. The earlier
    bilinear *_reg.nc files preserved the domain means but gave ~10% higher
    per-pixel STD and did NOT match the paper.
  - lat/lon "scale numbers" 22 pt; only 2 ticks each ([34.45, 34.55] / [-120.50,
    -120.40]) to declutter; latitude shown on the left-most map of each block.
  - colorbar titles split over two lines so the rotated labels never overlap.
  - panel letters (a).. OUTSIDE the maps (top-left), lowercase, matching Fig 4.
  - constrained_layout tightens columns/rows and enlarges the maps.
  - density insets on the Difference maps sit fully inside the lower-right corner.
  - LWC shown as g/cm2 (x10^-3), matching the Fig 4 histograms and Sec 3.2
    (trait-mean ~4.1 x10^-3 g/cm2); previously x10^-2.
Output: figures/figure5a_trait_vs_pft.{png,pdf}, figures/figure5b_difference.{png,pdf}
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pathlib import Path

BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')
# Use the SAME conservative-remapping (remapcon) dataset as Figure 4 and Results
# Sec 3.2 — this is what the manuscript Figure 5 was built from. (The old
# `clima_fit_prescribed_lai_ci` *_reg.nc files are bilinear-resampled: same domain
# means but ~10% higher per-pixel STD, which did NOT match the paper.)
RC = BASE / 'traits' / 'datasets' / 'clima_fit_prescribed_lai_ci_remapcon'
pft_path = BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'
FIG_DIR = BASE / 'shift_dangermond_trait' / 'figures'

CONV = {'chl': 1.0, 'lma': 1e4, 'lwc': 0.0018}   # -> µg/cm², g/m², g/cm²
TIMES = range(13)

# ---- Trait maps: time-average of the per-date remapcon retrievals ----
trait = {}
for _t in ['chl', 'lma', 'lwc']:
    stack = np.stack([
        xr.open_dataset(RC / f'shift_traits_time_{i:02d}_remapcon.nc')[_t]
          .transpose('lat', 'lon').values
        for i in TIMES])
    trait[_t] = np.nanmean(stack, axis=0) * CONV[_t]

# ---- PFT maps: pre-computed per-PFT means (remapcon) ----
pftm = {}
for _t in ['chl', 'lma', 'lwc']:
    da = xr.open_dataset(RC / f'mean_masked_{_t}_aviris_dangermond_clima_fit_remapcon.nc')[_t].squeeze()
    pftm[_t] = da.transpose('lat', 'lon').values * CONV[_t]

_grid = xr.open_dataset(RC / 'shift_traits_time_00_remapcon.nc')
lats = _grid['lat'].values
lons = _grid['lon'].values

pft_map = xr.open_dataset(pft_path)['Band1'].transpose('lat', 'lon').values
pft_mask = np.isin(pft_map, [2, 3, 4])

chl_trait, lma_trait, lwc_trait = trait['chl'], trait['lma'], trait['lwc']
chl_pft, lma_pft, lwc_pft = pftm['chl'], pftm['lma'], pftm['lwc']
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
    # Δ-distribution box tucked in the lower-right corner. Opaque white so it reads
    # cleanly over the map; x-units are omitted (given on the colorbar) so no text
    # runs off the map edge. Dashed lines mark zero (red) and the mean (blue).
    inset_ax = inset_axes(ax, width="33%", height="28%", loc='lower right', borderpad=0.9)
    v = data[np.isfinite(data) & (data != 0)]
    if len(v) == 0:
        return
    inset_ax.hist(v, bins=45, color='0.55', edgecolor='black', linewidth=0.4, density=True)
    inset_ax.axvline(0, color='red', linestyle='--', linewidth=1.4)
    inset_ax.axvline(np.nanmean(v), color='blue', linestyle='--', linewidth=1.4)
    inset_ax.tick_params(labelsize=INSET_TICK, length=2.5, pad=1.0)
    inset_ax.locator_params(axis='x', nbins=3)
    inset_ax.locator_params(axis='y', nbins=3)
    inset_ax.set_facecolor('white')
    inset_ax.patch.set_alpha(1.0)
    inset_ax.set_zorder(5)
    for s in inset_ax.spines.values():
        s.set_linewidth(0.7)

LON, LAT = np.meshgrid(lons, lats)

# two-line colorbar labels so the (rotated) titles do not overrun / overlap
CLABEL = {
    'chl': 'Chlorophyll\nContent (µg/cm²)',
    'lma': 'Leaf Mass per\nArea (g/m²)',
    'lwc': 'Leaf Water\nContent (g/cm²)',
}
LAT_TICKS = [34.45, 34.55]      # only 2 ticks each, to declutter
LON_TICKS = [-120.50, -120.40]


def style_map(ax, show_lat, show_lon):
    ax.set_aspect('equal')
    ax.set_xticks(LON_TICKS)
    ax.set_yticks(LAT_TICKS)
    ax.tick_params(labelsize=TICK)
    ax.tick_params(labelleft=show_lat, labelbottom=show_lon)


def panel_letter(ax, letter):
    # outside the map, above the top-left corner
    ax.annotate(f'({letter})', xy=(0, 1), xycoords='axes fraction',
                xytext=(-2, 3), textcoords='offset points',
                fontsize=PANEL, fontweight='bold', va='bottom', ha='right')


def stats_box(ax, disp):
    ax.text(0.965, 0.96, f'Mean: {np.nanmean(disp):.2f}\nSTD: {np.nanstd(disp):.2f}',
            transform=ax.transAxes, fontsize=STAT, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'))


# ===================== Figure 5a: Trait vs PFT (2 columns) =====================
figA, axesA = plt.subplots(3, 2, figsize=(11.0, 14.5), dpi=150, constrained_layout=True)
figA.set_constrained_layout_pads(w_pad=0.08, h_pad=0.02, wspace=0.01, hspace=0.02)
col_titles = ['Trait-based Approach', 'PFT-based Approach']
lettersA = [['a', 'b'], ['c', 'd'], ['e', 'f']]

for r, (trait_name, trait_data, pft_data, _diff) in enumerate(traits_data):
    config = trait_configs[trait_name]
    mult = config.get('unit_multiplier', 1) if trait_name == 'lwc' else 1
    im = None
    for c, data in enumerate([trait_data, pft_data]):
        ax = axesA[r, c]
        disp = data * mult
        im = ax.pcolormesh(LON, LAT, disp, cmap=config['cmap'], vmin=config['vmin'],
                           vmax=config['vmax'], shading='auto', rasterized=True)
        style_map(ax, show_lat=(c == 0), show_lon=(r == 2))
        panel_letter(ax, lettersA[r][c])
        stats_box(ax, disp)
        if r == 0:
            ax.set_title(col_titles[c], fontsize=TITLE, fontweight='bold', pad=22)
    cb = figA.colorbar(im, ax=[axesA[r, 0], axesA[r, 1]], fraction=0.05, pad=0.02)
    cb.set_label(CLABEL[trait_name], fontsize=CBAR_LBL, fontweight='bold')
    cb.ax.tick_params(labelsize=CBAR_TICK)
    if trait_name == 'lwc':
        cb.ax.set_title('(×10⁻³)', fontsize=16, pad=6)

out_a = FIG_DIR / 'figure5a_trait_vs_pft'
figA.savefig(f'{out_a}.png', dpi=300, facecolor='white')
figA.savefig(f'{out_a}.pdf', facecolor='white')
plt.close(figA)

# ================= Figure 5b: Difference (Trait - PFT) as a 1x3 row =================
figB, axesB = plt.subplots(1, 3, figsize=(18.0, 6.0), dpi=150, constrained_layout=True)
figB.set_constrained_layout_pads(w_pad=0.06, h_pad=0.04, wspace=0.04, hspace=0.02)
lettersB = ['a', 'b', 'c']
map_titles = {'chl': 'Chlorophyll Content', 'lma': 'Leaf Mass per Area', 'lwc': 'Leaf Water Content'}

for c, (trait_name, _t, _p, diff_data) in enumerate(traits_data):
    config = trait_configs[trait_name]
    mult = config.get('unit_multiplier', 1) if trait_name == 'lwc' else 1
    ax = axesB[c]
    disp = diff_data * mult
    im = ax.pcolormesh(LON, LAT, disp, cmap=config['cmap_diff'], vmin=config['vmin_diff'],
                       vmax=config['vmax_diff'], shading='auto', rasterized=True)
    style_map(ax, show_lat=(c == 0), show_lon=True)
    panel_letter(ax, lettersB[c])
    stats_box(ax, disp)
    ax.set_title(map_titles[trait_name], fontsize=TITLE, fontweight='bold', pad=22)
    cb = figB.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.set_label('Δ ' + CLABEL[trait_name], fontsize=CBAR_LBL, fontweight='bold')
    cb.ax.tick_params(labelsize=CBAR_TICK)
    if trait_name == 'lwc':
        cb.ax.set_title('(×10⁻³)', fontsize=16, pad=6)
    add_density_inset(ax, disp, config)

figB.suptitle('Difference (Trait − PFT)', fontsize=TITLE + 3, fontweight='bold')
out_b = FIG_DIR / 'figure5b_difference'
figB.savefig(f'{out_b}.png', dpi=300, facecolor='white')
figB.savefig(f'{out_b}.pdf', facecolor='white')
plt.close(figB)

print(f'Saved: {out_a}.png / .pdf')
print(f'Saved: {out_b}.png / .pdf')
