"""
Seasonal trait histograms, split by season into three separately-named,
print-legible figures (Ecosphere editor item 8). Kept as sub-figures of
Figure 4 (4a/4b/4c) for now, so the downstream figure numbers are undisturbed;
full renumbering will follow once Figures 5 and 6 are also revised.

    Figure 4a  ->  Early season (late Feb - Mar)   figure4a.{png,pdf}
    Figure 4b  ->  Mid season   (late Mar - Apr)   figure4b.{png,pdf}
    Figure 4c  ->  Late season  (May)              figure4c.{png,pdf}

Each figure is a single wide row of three panels (a) CHL, (b) LMA, (c) LWC, so
at the 18 cm page width the lettering is large and readable.
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import warnings
from pathlib import Path
from scipy import ndimage as _ndi

warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'Nimbus Sans',
    'font.size': 14, 'axes.labelsize': 15,
    'xtick.labelsize': 13, 'ytick.labelsize': 13,
    'legend.fontsize': 11, 'axes.labelweight': 'bold',
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Nimbus Sans',
    'mathtext.it': 'Nimbus Sans:italic', 'mathtext.bf': 'Nimbus Sans:bold',
})

BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')
OUT_DIR = BASE / 'traits' / 'datasets' / 'clima_fit_prescribed_lai_ci_remapcon'
PFT_FILE = BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'
FIG_DIR = BASE / 'shift_dangermond_trait' / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)

TRAIT_PAT = str(OUT_DIR / 'shift_traits_time_{t}_remapcon.nc')
PFT_VALUES = [2, 3, 4]
PFT_COLORS = {2: '#2E7D32', 3: '#1976D2', 4: '#D32F2F'}
PFT_NAMES = {2: 'PFT 1', 3: 'PFT 2', 4: 'PFT 3'}
SEASONS = {
    'Early': (['00', '01', '02', '03'],       'Early season (late Feb–Mar)', 'early'),
    'Mid':   (['04', '05', '06', '07', '08'],  'Mid season (late Mar–Apr)',   'mid'),
    'Late':  (['09', '10', '11', '12'],        'Late season (May)',           'late'),
}
# split kept as sub-figures 4a/4b/4c for now (downstream numbering undisturbed)
FIG_PANEL = {'Early': 'a', 'Mid': 'b', 'Late': 'c'}
FIG_LABEL = {'Early': '4a', 'Mid': '4b', 'Late': '4c'}
TRAIT_VARS = ['chl', 'lma', 'lwc']
TRAIT_NAME = {'chl': 'Chlorophyll Content', 'lma': 'Leaf Mass per Area', 'lwc': 'Leaf Water Content'}
TRAIT_UNIT = {'chl': r'CHL ($\mu$g cm$^{-2}$)', 'lma': r'LMA (g m$^{-2}$)',
              'lwc': r'LWC (g cm$^{-2}$, $\times$10$^{-3}$)'}
CONVERSIONS = {'chl': 1.0, 'lma': 1e4, 'lwc': 18.0 / 10000.0}
LWC_DISPLAY = 1e3
PANEL = ['a', 'b', 'c']


def _load_latlon(ds, var):
    da = ds[var].squeeze(); arr = da.values.astype(float); arr[arr <= -9000] = np.nan
    _lat = next((c for c in ds.coords if c.lower() in ('lat', 'latitude', 'y')), None)
    _lon = next((c for c in ds.coords if c.lower() in ('lon', 'longitude', 'x')), None)
    if _lat and _lon:
        if arr.shape[0] == len(ds[_lon]):
            arr = arr.T
        if ds[_lat].values[0] > ds[_lat].values[-1]:
            arr = arr[::-1, :]
    return arr

# PFT map (keep largest contiguous vegetated region, as in the original figure)
pft_ds = xr.open_dataset(PFT_FILE); pft_map = _load_latlon(pft_ds, 'Band1'); pft_ds.close()
_inside = np.isin(pft_map, PFT_VALUES)
_lab, _nc = _ndi.label(_inside)
_main = np.argmax(_ndi.sum(_inside, _lab, range(1, _nc + 1))) + 1
pft_map[_inside & (_lab != _main)] = np.nan
pft_flat = pft_map.flatten()


def collect(var, times, conv):
    acc = {p: [] for p in PFT_VALUES}
    for t in times:
        fp = Path(TRAIT_PAT.format(t=t))
        if not fp.exists():
            continue
        with xr.open_dataset(fp) as ds:
            arr = _load_latlon(ds, var) * conv
        for p in PFT_VALUES:
            vals = arr.ravel()[pft_flat == p]
            acc[p].append(vals[np.isfinite(vals)])
    return {p: (np.concatenate(acc[p]) if acc[p] else np.array([])) for p in PFT_VALUES}


# gather all seasons, and shared per-trait bin edges (so seasons are comparable)
data = {s: {v: collect(v, SEASONS[s][0], CONVERSIONS[v]) for v in TRAIT_VARS} for s in SEASONS}
bins = {}
for v in TRAIT_VARS:
    allv = np.concatenate([data[s][v][p] for s in SEASONS for p in PFT_VALUES if len(data[s][v][p])])
    if v == 'lwc':
        allv = allv * LWC_DISPLAY
    p2, p98 = np.nanpercentile(allv, [1, 99])
    bins[v] = np.linspace(p2, p98, 61)


def make_figure(season):
    times, subtitle, tag = SEASONS[season]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), dpi=150)
    for ci, v in enumerate(TRAIT_VARS):
        ax = axes[ci]
        for p in PFT_VALUES:
            raw = data[season][v][p]
            if len(raw) == 0:
                continue
            vals = raw * LWC_DISPLAY if v == 'lwc' else raw
            mu, sd = np.nanmean(vals), np.nanstd(vals)
            col = PFT_COLORS[p]
            ax.hist(vals, bins=bins[v], density=True, histtype='stepfilled', color=col, alpha=0.35, lw=0)
            ax.hist(vals, bins=bins[v], density=True, histtype='step', color=col, lw=1.8,
                    label=f'{PFT_NAMES[p]}: {mu:.2f} ± {sd:.2f}')
            ax.axvline(mu, color=col, ls='--', lw=2.0, alpha=0.9)
        ax.set_xlabel(TRAIT_UNIT[v], fontweight='bold')
        ax.set_ylabel('Density', fontweight='bold')
        if v == 'lma':
            ax.set_yscale('log')
        if v == 'chl':
            # only CHL needs headroom: its rising right tail runs under the upper-right legend.
            # LMA/LWC peak at the left, so the upper-right legend sits over empty space already.
            ax.set_ylim(top=ax.get_ylim()[1] * 1.4)
        ax.tick_params(length=3, width=0.8)
        ax.grid(True, alpha=0.3, ls='--', lw=0.5)
        # panel letter over the y-axis tick numbers
        ax.annotate(f'({PANEL[ci]})', xy=(0, 1), xycoords='axes fraction', xytext=(-22, 5),
                    textcoords='offset points', fontsize=18, fontweight='bold', va='bottom', ha='center')
        # legend INSIDE, upper-right (clear of the left-side mean lines and the peaks)
        ax.legend(loc='upper right', fontsize=11, frameon=True, framealpha=0.9, fancybox=True,
                  handlelength=1.3, borderpad=0.4, labelspacing=0.3)
    fig.suptitle(subtitle, fontsize=20, fontweight='bold', y=0.99)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.82, bottom=0.17, wspace=0.26)
    stem = FIG_DIR / f'figure4{FIG_PANEL[season]}'
    fig.savefig(f'{stem}.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(f'{stem}.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'saved Figure {FIG_LABEL[season]} ({season}) -> {stem}.png / .pdf')


for s in SEASONS:
    make_figure(s)
print('DONE')
