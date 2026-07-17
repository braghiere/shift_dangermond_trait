"""
Figure 3 — Temporal dynamics of LAI/CHL/LMA/LWC by PFT (remapcon 30 m data).

Ported from review/figures_5m_remapcon_traits_fixed.ipynb (cells 2/3/7/8/9) into a
standalone, reproducible script so Figure 3 matches the rest of the paper's
publication styling:
  * Nimbus Sans everywhere + custom mathtext (superscripts as mathtext, so no
    missing-glyph boxes — Nimbus Sans lacks the literal superscript-minus)
  * pdf/ps.fonttype = 42 (embed TrueType; no Type-3 fonts)
  * savefig at 600 dpi (PNG + PDF)
  * LWC standardized to ×10⁻³ g cm⁻² (consistent with Figs 2/4/5)

2×2 grid: (a) LAI + ERA5 precip/temp, (b) CHL, (c) LMA, (d) LWC.
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator
from pathlib import Path
from datetime import datetime
import pandas as pd
from scipy import ndimage as _ndi
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'Nimbus Sans',
    'font.size': 14, 'axes.labelsize': 14, 'axes.titlesize': 14,
    'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'legend.fontsize': 12, 'axes.labelweight': 'bold',
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Nimbus Sans',
    'mathtext.it': 'Nimbus Sans:italic', 'mathtext.bf': 'Nimbus Sans:bold',
    'pdf.fonttype': 42, 'ps.fonttype': 42,   # embed TrueType, no Type-3
})

# ---------------- config ----------------
BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')
RC = BASE / 'traits' / 'datasets' / 'clima_fit_prescribed_lai_ci_remapcon'
PFT_FILE = BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'
ERA5_TEMP = '/home/renatob/era5/era5_monthly_ILAMB/output_cdo_ilamb_forcing_cdo/ERA5/tas/era5_tas_mon_1979-2024.nc'
ERA5_PRECIP = '/home/renatob/era5/era5_monthly_ILAMB/output_cdo_ilamb_forcing_cdo/ERA5/pr/era5_pr_mon_1979-2024.nc'
FIG_DIR = BASE / 'shift_dangermond_trait' / 'figures'

TRAIT_PAT = str(RC / 'shift_traits_time_{t}_remapcon.nc')
LAI_PAT = str(RC / 'lai_aviris_dangermond_time_{t}_remapcon.nc')

DATES = [datetime(2022, 2, 24), datetime(2022, 2, 28), datetime(2022, 3, 7),
         datetime(2022, 3, 14), datetime(2022, 3, 21), datetime(2022, 3, 28),
         datetime(2022, 4, 4), datetime(2022, 4, 11), datetime(2022, 4, 18),
         datetime(2022, 4, 25), datetime(2022, 5, 2), datetime(2022, 5, 9),
         datetime(2022, 5, 23)]
TIMES = [f'{i:02d}' for i in range(13)]

PFT_VALUES = [2, 3, 4]
PFT_LABELS = {2: 'PFT 1', 3: 'PFT 2', 4: 'PFT 3'}
PFT_COLORS = {2: '#2E7D32', 3: '#1976D2', 4: '#D32F2F'}

# CHL: µg cm⁻² (×1); LMA: g cm⁻² ×1e4 -> g m⁻²; LWC: mol m⁻² ×(18/1e4) -> g cm⁻²,
# then ×1000 -> ×10⁻³ g cm⁻² (paper-standard display); LAI: m² m⁻² (×1)
CONVERSIONS = {'chl': 1.0, 'lma': 1e4, 'lwc': (18. / 10000.) * 1000., 'lai': 1.0}

_panel_cfg = [
    ('lai', '(a) Temporal Dynamics of LAI by PFT',         r'Leaf Area Index (LAI) [m$^2$ m$^{-2}$]'),
    ('chl', '(b) Temporal Dynamics of Chlorophyll by PFT', r'Chlorophyll [$\mu$g cm$^{-2}$]'),
    ('lma', '(c) Temporal Dynamics of LMA by PFT',         r'LMA [g m$^{-2}$]'),
    ('lwc', '(d) Temporal Dynamics of LWC by PFT',         r'LWC ($\times$10$^{-3}$ g cm$^{-2}$)'),
]


# ---------------- data ----------------
def _load_latlon(ds, var):
    """Return 2-D (lat_ascending, lon) float array — fill-masked, transposed, lat-flipped."""
    arr = ds[var].squeeze().values.astype(float)
    arr[arr <= -9000] = np.nan
    _lat = next((c for c in ds.coords if c.lower() in ('lat', 'latitude', 'y')), None)
    _lon = next((c for c in ds.coords if c.lower() in ('lon', 'longitude', 'x')), None)
    if _lat and _lon:
        if arr.shape[0] == len(ds[_lon]):          # first dim is lon -> transpose
            arr = arr.T
        if ds[_lat].values[0] > ds[_lat].values[-1]:   # north-first -> flip to south-first
            arr = arr[::-1, :]
    return arr


pft_ds = xr.open_dataset(PFT_FILE)
pft_map = _load_latlon(pft_ds, 'Band1')
pft_ds.close()
# remove stray ocean/coastal PFT fragments (keep the main connected preserve)
_inside = np.isin(pft_map, [2, 3, 4])
_lab, _n = _ndi.label(_inside)
_sizes = _ndi.sum(_inside, _lab, range(1, _n + 1))
_main = int(np.argmax(_sizes)) + 1
pft_map[_inside & (_lab != _main)] = np.nan
pft_flat = pft_map.flatten()


def load_var(var, t):
    path = LAI_PAT.format(t=t) if var == 'lai' else TRAIT_PAT.format(t=t)
    nc_var = 'lai' if var == 'lai' else var
    if not Path(path).exists():
        return None
    ds = xr.open_dataset(path)
    if nc_var not in ds:
        ds.close()
        return None
    data = _load_latlon(ds, nc_var).ravel()
    ds.close()
    return data


stats = {v: {k: {p: [] for p in PFT_VALUES} for k in ('mean', 'p25', 'p75')}
         for v in ['lai', 'chl', 'lma', 'lwc']}
for var in ['lai', 'chl', 'lma', 'lwc']:
    conv = CONVERSIONS[var]
    for t in TIMES:
        raw = load_var(var, t)
        for p in PFT_VALUES:
            if raw is None:
                vals = np.array([])
            else:
                data = raw * conv
                vals = data[(pft_flat == p) & np.isfinite(data)]
            stats[var]['mean'][p].append(np.nanmean(vals) if vals.size else np.nan)
            stats[var]['p25'][p].append(np.nanpercentile(vals, 25) if vals.size else np.nan)
            stats[var]['p75'][p].append(np.nanpercentile(vals, 75) if vals.size else np.nan)
for var in stats:
    for k in stats[var]:
        for p in PFT_VALUES:
            stats[var][k][p] = np.array(stats[var][k][p])

# ERA5 climate for panel (a)
DANG_LAT, DANG_LON = 34.51, -120.50 + 360.
era5_dates = era5_precip = era5_temp = None
try:
    ds_t = xr.open_dataset(ERA5_TEMP)
    ds_p = xr.open_dataset(ERA5_PRECIP)
    t22 = ds_t.sel(time=slice('2022-01-01', '2022-12-31'))
    p22 = ds_p.sel(time=slice('2022-01-01', '2022-12-31'))
    temp_ts = t22['tas'].sel(lat=DANG_LAT, lon=DANG_LON, method='nearest')
    precip_ts = p22['pr'].sel(lat=DANG_LAT, lon=DANG_LON, method='nearest')
    era5_dates = [datetime(int(str(d)[:4]), int(str(d)[5:7]), 15) for d in temp_ts.time.values]
    era5_temp = temp_ts.values - 273.15
    era5_precip = precip_ts.values * 86400. * 30.
    ds_t.close(); ds_p.close()
except Exception as e:
    print(f'ERA5 load failed: {e} — panel (a) climate overlay omitted')

# ---------------- figure ----------------
fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=150)
axes = axes.flatten()
for idx, (tkey, panel_title, ylabel) in enumerate(_panel_cfg):
    ax = axes[idx]
    for p in PFT_VALUES:
        col = PFT_COLORS[p]
        mean_v, p25_v, p75_v = stats[tkey]['mean'][p], stats[tkey]['p25'][p], stats[tkey]['p75'][p]
        ax.fill_between(DATES, p25_v, p75_v, color=col, alpha=0.2, linewidth=0, zorder=2)
        ax.plot(DATES, mean_v, '-o', color=col, lw=2.5, ms=6, markeredgecolor='white',
                markeredgewidth=1, label=PFT_LABELS[p], zorder=3)
        ax.axhline(np.nanmean(mean_v), color=col, lw=1.8, ls='-', alpha=0.75, zorder=3)
    ax.set_xlabel('Date', fontsize=14, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')
    ax.set_title(panel_title, fontsize=14, fontweight='bold', loc='left', pad=10)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=7))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, zorder=1)
    ax.set_xlim(DATES[0] - pd.Timedelta(days=3), DATES[-1] + pd.Timedelta(days=3))
    ax.legend(loc='upper right', frameon=True, fancybox=True, shadow=True, fontsize=11).set_zorder(4)

    if tkey == 'lai' and era5_dates is not None:
        ax_p = ax.twinx()
        ax_p.bar(era5_dates, era5_precip, width=20, color='#1976D2', alpha=0.25,
                 label='Precip.', zorder=1, align='center')
        ax_p.set_ylabel('Precipitation [mm/month]', fontsize=12, color='#1976D2', fontweight='bold')
        ax_p.tick_params(axis='y', labelcolor='#1976D2', labelsize=10)
        ax_p.set_ylim(0, max(era5_precip) * 2.5)
        ax_t = ax.twinx()
        ax_t.spines['right'].set_position(('outward', 60))
        ax_t.plot(era5_dates, era5_temp, color='#FBC02D', lw=2.0, marker='s', ms=4,
                  markeredgecolor='white', markeredgewidth=0.5, label='Temp.', alpha=0.9, zorder=2)
        ax_t.set_ylabel('Temperature [°C]', fontsize=12, color='#FBC02D', fontweight='bold')
        ax_t.tick_params(axis='y', labelcolor='#FBC02D', labelsize=10)
        lp, llp = ax_p.get_legend_handles_labels()
        lt, llt = ax_t.get_legend_handles_labels()
        ax_p.legend(lp + lt, llp + llt, loc='upper left', frameon=True, fancybox=True,
                    shadow=True, fontsize=10).set_zorder(4)

plt.tight_layout()
out = FIG_DIR / 'figure3_remapcon_traits_temporal'
fig.savefig(f'{out}.png', dpi=600, bbox_inches='tight', facecolor='white')
fig.savefig(f'{out}.pdf', dpi=600, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {out}.png / .pdf')
