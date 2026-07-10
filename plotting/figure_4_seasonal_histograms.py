"""Regenerate figure4_remapcon_seasonal_histograms.png with PFT 1/2/3 labels."""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import warnings
from pathlib import Path
from datetime import datetime
from scipy import ndimage as _ndi

warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.size': 14, 'axes.labelsize': 14, 'axes.titlesize': 14,
    'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'legend.fontsize': 12, 'axes.labelweight': 'bold',
    'font.family': 'sans-serif',
})

BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')
OUT_DIR = BASE / 'traits' / 'datasets' / 'clima_fit_prescribed_lai_ci_remapcon'
PFT_FILE = BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'
FIG_DIR = BASE / 'shift_dangermond_trait' / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)

TRAIT_REMAPCON_PAT = str(OUT_DIR / 'shift_traits_time_{t}_remapcon.nc')
TIMES = [f'{i:02d}' for i in range(13)]
PFT_VALUES = [2, 3, 4]
PFT_COLORS = {2: '#2E7D32', 3: '#1976D2', 4: '#D32F2F'}
SEASONS = {
    'Early': ['00', '01', '02', '03'],
    'Mid':   ['04', '05', '06', '07', '08'],
    'Late':  ['09', '10', '11', '12'],
}
TRAIT_VARS  = ['chl', 'lma', 'lwc']
TRAIT_UNITS = {'chl': 'CHL [µg cm⁻²]', 'lma': 'LMA [g m⁻²]', 'lwc': 'LWC [g cm⁻² ×10⁻³]'}
CONVERSIONS = {'chl': 1.0, 'lma': 1e4, 'lwc': 18.0 / 10000.0}
_LWC_DISPLAY = 1e3

_pft_names = {2: 'PFT 1', 3: 'PFT 2', 4: 'PFT 3'}


def _load_latlon(ds, var):
    da  = ds[var].squeeze()
    arr = da.values.astype(float)
    arr[arr <= -9000] = np.nan
    _lat = next((c for c in ds.coords if c.lower() in ('lat', 'latitude', 'y')), None)
    _lon = next((c for c in ds.coords if c.lower() in ('lon', 'longitude', 'x')), None)
    if _lat and _lon:
        nlon = len(ds[_lon])
        if arr.shape[0] == nlon:
            arr = arr.T
        if ds[_lat].values[0] > ds[_lat].values[-1]:
            arr = arr[::-1, :]
    return arr


# Load PFT map
pft_ds  = xr.open_dataset(PFT_FILE)
pft_map = _load_latlon(pft_ds, 'Band1')
pft_ds.close()

_inside = np.isin(pft_map, [2, 3, 4])
_labeled, _nc = _ndi.label(_inside)
_sizes = _ndi.sum(_inside, _labeled, range(1, _nc + 1))
_main  = np.argmax(_sizes) + 1
pft_map[_inside & (_labeled != _main)] = np.nan
pft_flat = pft_map.flatten()
print(f"PFT map: {pft_map.shape}")


def collect_season_data(var, times, conv=1.0):
    acc = {p: [] for p in PFT_VALUES}
    for t in times:
        fpath = Path(TRAIT_REMAPCON_PAT.format(t=t))
        if not fpath.exists():
            continue
        with xr.open_dataset(fpath) as ds:
            arr = _load_latlon(ds, var) * conv
        for p in PFT_VALUES:
            mask = pft_flat == p
            vals = arr.ravel()[mask]
            vals = vals[np.isfinite(vals)]
            acc[p].append(vals)
    return {p: np.concatenate(acc[p]) if acc[p] else np.array([]) for p in PFT_VALUES}


season_acc = {}
for sname, stimes in SEASONS.items():
    season_acc[sname] = {}
    for var in TRAIT_VARS:
        season_acc[sname][var] = collect_season_data(var, stimes, CONVERSIONS[var])
    print(f"  {sname} done")

# Global bin edges per trait
N_BINS = 60
bin_edges = {}
for var in TRAIT_VARS:
    all_vals = np.concatenate([season_acc[sn][var][p]
                               for sn in SEASONS for p in PFT_VALUES
                               if len(season_acc[sn][var][p]) > 0])
    if var == 'lwc':
        all_vals = all_vals * _LWC_DISPLAY
    p2, p98 = np.nanpercentile(all_vals, [1, 99])
    bin_edges[var] = np.linspace(p2, p98, N_BINS + 1)

_col_titles = {'chl': 'CHL [µg cm⁻²]', 'lma': 'LMA [g m⁻²]', 'lwc': 'LWC [g cm⁻² ×10⁻³]'}
_row_titles = {
    'Early': 'Early Season\n(late Feb – Mar)',
    'Mid':   'Mid Season\n(late Mar – Apr)',
    'Late':  'Late Season\n(May)',
}
_panel_letters = [['a.', 'b.', 'c.'], ['d.', 'e.', 'f.'], ['g.', 'h.', 'i.']]

fig, axes = plt.subplots(3, 3, figsize=(16, 13), dpi=150)
season_list = list(SEASONS.keys())

for ri, sname in enumerate(season_list):
    for ci, var in enumerate(TRAIT_VARS):
        ax   = axes[ri, ci]
        bins = bin_edges[var]

        for p in PFT_VALUES:
            raw_vals = season_acc[sname][var][p]
            if len(raw_vals) == 0:
                continue
            vals = raw_vals * _LWC_DISPLAY if var == 'lwc' else raw_vals
            col  = PFT_COLORS[p]
            mu   = np.nanmean(vals)
            sd   = np.nanstd(vals)
            lbl  = f'{_pft_names[p]}: {mu:.2f} ± {sd:.2f}'

            ax.hist(vals, bins=bins, density=True, histtype='stepfilled',
                    color=col, alpha=0.35, linewidth=0)
            ax.hist(vals, bins=bins, density=True, histtype='step',
                    color=col, linewidth=1.5, label=lbl)
            ax.axvline(mu, color=col, linestyle='--', linewidth=2.0, alpha=0.90)

        if ri == 0:
            ax.set_title(_col_titles[var], fontsize=14, fontweight='bold')
        if ci == 0:
            ax.set_ylabel(_row_titles[sname], fontsize=12, fontweight='bold')
        else:
            ax.set_ylabel('Density', fontsize=11)
        if ri == 2:
            ax.set_xlabel(TRAIT_UNITS[var], fontsize=12, fontweight='bold')

        ax.annotate(_panel_letters[ri][ci], xy=(0, 1), xycoords='axes fraction',
                    xytext=(-2, 6), textcoords='offset points',
                    fontsize=14, fontweight='bold', va='bottom', ha='right')

        if var == 'lma':
            ax.set_yscale('log')

        ax.tick_params(axis='both', labelsize=10)
        ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, zorder=0)

        leg = ax.legend(loc='upper right', frameon=True, fancybox=True,
                        shadow=True, fontsize=9)
        leg.set_zorder(4)

plt.tight_layout()
out_png = FIG_DIR / 'figure4_remapcon_seasonal_histograms.png'
out_pdf = FIG_DIR / 'figure4_remapcon_seasonal_histograms.pdf'
fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(out_pdf, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
