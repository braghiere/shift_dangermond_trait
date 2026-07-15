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
  - 5b: the map extent is padded right/below the preserve so each frame's lower-right
    is blank white; the Δ-distribution inset (with axis names + units) lives there,
    clear of the data.
  - LWC shown as g/cm2 (x10^-3), matching the Fig 4 histograms and Sec 3.2
    (trait-mean ~4.1 x10^-3 g/cm2); previously x10^-2.
Output: figures/figure5a_trait_vs_pft.{png,pdf}, figures/figure5b_difference.{png,pdf}
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
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

# ---- Font sizes (editor item 9: legible scale numbers at the 18 cm print width) ----
TICK = 22        # lat/lon scale numbers
CBAR_TICK = 22   # colorbar scale numbers
CBAR_LBL = 18    # colorbar / histogram axis labels
STAT = 16        # in-map Mean/STD box
PANEL = 24       # (a).. panel letters
TITLE = 23       # column / map titles
HIST_TICK = 16   # Fig 5b histogram tick numbers

LON, LAT = np.meshgrid(lons, lats)
ASP = float((lons.max() - lons.min()) / (lats.max() - lats.min()))   # map width:height, equal-degree

# two-line colorbar labels so the (rotated) titles do not overrun / overlap
CLABEL = {
    'chl': 'Chlorophyll\nContent (µg/cm²)',
    'lma': 'Leaf Mass per\nArea (g/m²)',
    'lwc': 'Leaf Water\nContent (g/cm²)',
}
# Δ-distribution histogram x-axis labels (units omitted — they are on the colorbar)
DLABEL = {'chl': 'Δ CHL', 'lma': 'Δ LMA', 'lwc': 'Δ LWC'}
DRANGE = {'chl': 50, 'lma': 200, 'lwc': 20}        # symmetric x-limits for the Δ insets
LAT_TICKS = [34.45, 34.55]      # only 2 ticks each, to declutter
LON_TICKS = [-120.50, -120.40]


def style_map(ax, show_lat, show_lon, xlim=None, ylim=None):
    # maps fill their (aspect-matched) axes box, so a same-box colorbar == map height
    ax.set_xlim(xlim if xlim else (lons.min(), lons.max()))
    ax.set_ylim(ylim if ylim else (lats.min(), lats.max()))
    ax.set_xticks(LON_TICKS)
    ax.set_yticks(LAT_TICKS)
    ax.tick_params(labelsize=TICK)
    ax.tick_params(labelleft=show_lat, labelbottom=show_lon)


def panel_letter(ax, letter):
    ax.annotate(f'({letter})', xy=(0, 1), xycoords='axes fraction',
                xytext=(-2, 4), textcoords='offset points',
                fontsize=PANEL, fontweight='bold', va='bottom', ha='right')


def stats_box(ax, disp):
    ax.text(0.965, 0.96, f'Mean: {np.nanmean(disp):.2f}\nSTD: {np.nanstd(disp):.2f}',
            transform=ax.transAxes, fontsize=STAT, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'))


# ============ Figure 5a: Trait vs PFT (manual layout: colorbar == map height) ============
Wf, Hf = 11.0, 13.3
L, G_COL, G_CB, CW, R_LAB = 0.100, 0.020, 0.014, 0.022, 0.140
MW = (1 - L - G_COL - G_CB - CW - R_LAB) / 2      # map width (fig fraction)
MH = MW * Wf / (ASP * Hf)                          # map height: box aspect == ASP
G_ROW, B = 0.038, 0.050                            # row gap just fits the outside letters (#3)
TOP = 1 - 3 * MH - 2 * G_ROW - B                   # region above row-0 maps (titles)
COL_X = [L, L + MW + G_COL]
CB_X = L + 2 * MW + G_COL + G_CB
col_titles = ['Trait-based Approach', 'PFT-based Approach']
lettersA = [['a', 'b'], ['c', 'd'], ['e', 'f']]

figA = plt.figure(figsize=(Wf, Hf), dpi=150)
for r, (trait_name, trait_data, pft_data, _diff) in enumerate(traits_data):
    config = trait_configs[trait_name]
    mult = config.get('unit_multiplier', 1) if trait_name == 'lwc' else 1
    y0 = 1 - TOP - r * (MH + G_ROW) - MH
    im = None
    for c, data in enumerate([trait_data, pft_data]):
        ax = figA.add_axes([COL_X[c], y0, MW, MH])
        disp = data * mult
        im = ax.pcolormesh(LON, LAT, disp, cmap=config['cmap'], vmin=config['vmin'],
                           vmax=config['vmax'], shading='auto', rasterized=True)
        style_map(ax, show_lat=(c == 0), show_lon=(r == 2))
        panel_letter(ax, lettersA[r][c])
        stats_box(ax, disp)
        if r == 0:                              # extra title->map space (#2)
            ax.set_title(col_titles[c], fontsize=TITLE, fontweight='bold', pad=34)
    cax = figA.add_axes([CB_X, y0, CW, MH])     # colorbar exactly the map height (#1)
    cb = figA.colorbar(im, cax=cax)
    cb.set_label(CLABEL[trait_name], fontsize=CBAR_LBL, fontweight='bold')
    cb.ax.tick_params(labelsize=CBAR_TICK)
    if trait_name == 'lwc':
        cb.ax.set_title('(×10⁻³)', fontsize=16, pad=6)

out_a = FIG_DIR / 'figure5a_trait_vs_pft'
figA.savefig(f'{out_a}.png', dpi=300, facecolor='white')
figA.savefig(f'{out_a}.pdf', facecolor='white')
plt.close(figA)

# ===== Figure 5b: Difference (Trait − PFT), 1x3 row with Δ-distribution insets =====
# The map extent is padded to the right and below the preserve so the lower-right
# of each frame is blank white; the Δ-distribution inset lives there with its full
# axis names + units, clear of the data (user request). Aspect stays consistent so
# the colorbars remain exactly the map height.
lon_r = lons.max() - lons.min()
lat_r = lats.max() - lats.min()
PAD_R, PAD_B = 0.35, 0.28                          # modest white pad (the map is the star)
XLIM = (lons.min(), lons.max() + PAD_R * lon_r)
YLIM = (lats.min() - PAD_B * lat_r, lats.max())
ASP_B = (XLIM[1] - XLIM[0]) / (YLIM[1] - YLIM[0])  # padded-frame aspect

Wb, Hb = 18.0, 4.9
Lb, GCB, CWb, CLAB, GUNIT, RB = 0.055, 0.006, 0.013, 0.052, 0.035, 0.006  # GUNIT = space between panels
MWb = (1 - Lb - 3 * (GCB + CWb + CLAB) - 2 * GUNIT - RB) / 3
MHb = MWb * Wb / (ASP_B * Hb)                      # map height: box aspect == ASP_B
TOPb = 0.17                                         # top margin (suptitle + map title)
map_y = 1 - TOPb - MHb
lettersB = ['a', 'b', 'c']
map_titles = {'chl': 'Chlorophyll Content', 'lma': 'Leaf Mass per Area', 'lwc': 'Leaf Water Content'}
INSET_LBL, INSET_TICK = 15, 13

figB = plt.figure(figsize=(Wb, Hb), dpi=150)
for c, (trait_name, _t, _p, diff_data) in enumerate(traits_data):
    config = trait_configs[trait_name]
    mult = config.get('unit_multiplier', 1) if trait_name == 'lwc' else 1
    x = Lb + c * (MWb + GCB + CWb + CLAB + GUNIT)
    disp = diff_data * mult
    ax = figB.add_axes([x, map_y, MWb, MHb])
    im = ax.pcolormesh(LON, LAT, disp, cmap=config['cmap_diff'], vmin=config['vmin_diff'],
                       vmax=config['vmax_diff'], shading='auto', rasterized=True)
    style_map(ax, show_lat=(c == 0), show_lon=True, xlim=XLIM, ylim=YLIM)
    panel_letter(ax, lettersB[c])
    stats_box(ax, disp)
    ax.set_title(map_titles[trait_name], fontsize=TITLE, fontweight='bold', pad=8)
    cax = figB.add_axes([x + MWb + GCB, map_y, CWb, MHb])   # colorbar == map height (#1)
    cb = figB.colorbar(im, cax=cax)
    cb.set_label('Δ ' + CLABEL[trait_name], fontsize=CBAR_LBL, fontweight='bold')
    cb.ax.tick_params(labelsize=CBAR_TICK)
    if trait_name == 'lwc':
        cb.ax.set_title('(×10⁻³)', fontsize=16, pad=6)
    # Δ-distribution histogram INSET tucked into the blank lower-right corner
    # (secondary to the map); fixed symmetric x-range per trait, no units on label.
    iw, ih = 0.36 * MWb, 0.30 * MHb
    hax = figB.add_axes([x + MWb - iw - 0.02 * MWb, map_y + 0.16 * MHb, iw, ih])
    v = disp[np.isfinite(disp) & (disp != 0)]
    hax.hist(v, bins=45, color='0.55', edgecolor='black', linewidth=0.4, density=True)
    hax.axvline(0, color='red', linestyle='--', linewidth=1.5)
    hax.axvline(np.nanmean(v), color='blue', linestyle='--', linewidth=1.5)
    hax.set_xlim(-DRANGE[trait_name], DRANGE[trait_name])
    hax.set_xticks([-DRANGE[trait_name], 0, DRANGE[trait_name]])
    hax.set_xlabel(DLABEL[trait_name], fontsize=INSET_LBL, fontweight='bold', labelpad=1.5)
    hax.tick_params(labelsize=INSET_TICK, length=2.5, pad=1.5)
    hax.locator_params(axis='y', nbins=2)

figB.text(0.5, 0.975, 'Difference (Trait − PFT)', ha='center', va='top',
          fontsize=TITLE + 3, fontweight='bold')
out_b = FIG_DIR / 'figure5b_difference'
figB.savefig(f'{out_b}.png', dpi=300, facecolor='white')
figB.savefig(f'{out_b}.pdf', facecolor='white')
plt.close(figB)

print(f'Saved: {out_a}.png / .pdf')
print(f'Saved: {out_b}.png / .pdf')
