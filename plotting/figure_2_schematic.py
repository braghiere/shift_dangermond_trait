"""
Figure 2 — CliMA Land workflow schematic (publication build).

  LEFT  : Traits  (Chlorophyll, Leaf Mass per Area, Leaf Water Content)
  CENTER: AVIRIS-NG reflectance spectrum + LiDAR point cloud -> CliMA Land hub;
          Inverse modeling (left) / Forward modeling (right);
          ERA5 reanalysis climograph (bottom)
  RIGHT : Fluxes  (GPP, SIF740)

Geographic maps (lat/lon graticule, scale bar + north arrow on the top map of
each column), cropped to the preserve polygon. Vector PDF (map layers rasterized
at 600 dpi) + 600-dpi PNG.  Units match the manuscript; panel means/STD go in the
caption.
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyBboxPatch, FancyArrowPatch
from matplotlib.ticker import FuncFormatter, FixedLocator
from pathlib import Path

plt.rcParams.update({
    'font.family': 'Nimbus Sans', 'font.size': 7,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Nimbus Sans',
    'mathtext.it': 'Nimbus Sans:italic', 'mathtext.bf': 'Nimbus Sans:bold',
    'axes.linewidth': 0.5, 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
TEAL = '#2b6c8f'; LABEL = '#cfe2f3'; LABEL_TX = '#123449'
PBLUE = '#1976D2'; TRED = '#D84315'
FS_HUB, FS_BOX, FS_TAG, FS_TITLE, FS_CB, FS_TICK = 15, 10, 9, 8.5, 7, 5.5

BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')
TRAIT = BASE / 'traits' / 'datasets' / 'clima_fit_prescribed_lai_ci'
FLUX = BASE / 'fitting' / 'fitted_prescribed_lai_ci'
PFT_FILE = BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'
REFL_FILE = BASE / 'aviris_dangermond_time_00' / 'output_clipped_reduced.nc'
ERA5_T = Path('/home/renatob/era5/era5_monthly_ILAMB/output_cdo_ilamb_forcing_cdo/ERA5/tas/era5_tas_mon_1979-2024.nc')
ERA5_P = Path('/home/renatob/era5/era5_monthly_ILAMB/output_cdo_ilamb_forcing_cdo/ERA5/pr/era5_pr_mon_1979-2024.nc')
OUT = BASE / 'shift_dangermond_trait' / 'figures'

def orient(a, nlat, nlon):
    a = np.asarray(a, float); a[a <= -9000] = np.nan
    return a.T if a.shape == (nlon, nlat) else a

# ---------------- data ----------------
chl_ds = xr.open_dataset(TRAIT / 'chl_aviris_dangermond_clima_fit_reg.nc')
lon = chl_ds['lon'].values; lat = chl_ds['lat'].values
nlat, nlon = len(lat), len(lon)
LON, LAT = np.meshgrid(lon, lat)
pft = orient(xr.open_dataset(PFT_FILE)['Band1'].values, nlat, nlon)
mask = np.isin(pft, [2, 3, 4])
mp = lambda a: np.where(mask, a, np.nan)

chl = mp(chl_ds['chl'].mean('time').values)
lma = mp(xr.open_dataset(TRAIT / 'lma_aviris_dangermond_clima_fit_reg.nc')['lma'].mean('time').values * 1e4)
lwc = mp(xr.open_dataset(TRAIT / 'lwc_aviris_dangermond_clima_fit_reg.nc')['lwc'].mean('time').values * 0.0018 * 100)
g, s = [], []
for t in range(13):
    f = FLUX / f'shift_fluxes_day_{t:02d}_clima_fit_reg_jmax.nc'
    if f.exists():
        d = xr.open_dataset(f, decode_times=False)
        g.append(orient(d['gpp'].squeeze().values, nlat, nlon))
        s.append(orient(d['sif740'].squeeze().values, nlat, nlon))
gpp = mp(np.nanmean(np.dstack(g), 2)); sif = mp(np.nanmean(np.dstack(s), 2))

ERA_MONTHS = ['Feb', 'Mar', 'Apr', 'May']
try:
    _lon = -120.44 + 360
    _t = (xr.open_dataset(ERA5_T)['tas'].sel(lat=34.5, lon=_lon, method='nearest') - 273.15).sel(time='2022')
    _p = (xr.open_dataset(ERA5_P)['pr'].sel(lat=34.5, lon=_lon, method='nearest') * 86400 * 30).sel(time='2022')
    mo = _t['time.month'].values
    ERA_T = [float(_t.values[mo == m][0]) for m in (2, 3, 4, 5)]
    ERA_P = [float(_p.values[mo == m][0]) for m in (2, 3, 4, 5)]
except Exception as e:
    print('ERA5 load failed, using stored campaign values:', e)
    ERA_T = [13.0, 13.5, 14.6, 16.0]; ERA_P = [204.3, 104.1, 37.3, 59.8]

vlon = LON[mask]; vlat = LAT[mask]
XLIM = (vlon.min() - 0.006, vlon.max() + 0.006)
YLIM = (vlat.min() - 0.006, vlat.max() + 0.006)
WX, HY = XLIM[1] - XLIM[0], YLIM[1] - YLIM[0]
MEANLAT = float(np.mean(vlat)); GEO_ASPECT = 1 / np.cos(np.deg2rad(MEANLAT))

# ---------------- figure & layout ----------------
fig = plt.figure(figsize=(7.09, 4.72), dpi=200)
FASP = 7.09 / 4.72
ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, 1); ov.set_ylim(0, 1); ov.axis('off')

MW, MH = 0.138, 0.215
LX, RX = 0.040, 0.782
CB_W, CB_GAP = 0.012, 0.004
TOP, BOT = 0.860, 0.061   # shared top/bottom of both map columns
lon_fmt = FuncFormatter(lambda v, _: f"{abs(v):.2f}°W")
lat_fmt = FuncFormatter(lambda v, _: f"{v:.2f}°N")

def _scalebar(ax, km=3):
    dlon = km / (111.320 * np.cos(np.deg2rad(MEANLAT)))
    x0 = XLIM[1] - 0.06 * WX - dlon; y0 = YLIM[0] + 0.07 * HY   # lower-right (white)
    ax.plot([x0, x0 + dlon], [y0, y0], color='black', lw=1.4, solid_capstyle='butt', zorder=6)
    ax.text(x0 + dlon / 2, y0 + 0.02 * HY, f'{km} km', ha='center', va='bottom', fontsize=5, zorder=6)

def _north(ax):
    x = XLIM[1] - 0.09 * WX; y0 = YLIM[1] - 0.30 * HY; dy = 0.14 * HY   # upper-right (scale-bar side)
    ax.annotate('', xy=(x, y0 + dy), xytext=(x, y0),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=1.0), zorder=6)
    ax.text(x, y0 + dy + 0.008 * HY, 'N', ha='center', va='bottom', fontsize=5, fontweight='bold', zorder=6)

def add_map(x, b, data, cmap, vmin, vmax, name, unit, show_lat=False, show_lon=False, deco=False):
    ax = fig.add_axes([x, b, MW, MH])
    cm = plt.get_cmap(cmap).copy(); cm.set_bad(alpha=0)
    im = ax.pcolormesh(LON, LAT, data, cmap=cm, vmin=vmin, vmax=vmax, shading='auto', rasterized=True)
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM); ax.set_aspect(GEO_ASPECT)
    ax.grid(True, color='0.7', lw=0.3, alpha=0.6)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.xaxis.set_major_locator(FixedLocator([-120.50, -120.45, -120.40]))
    ax.yaxis.set_major_locator(FixedLocator([34.45, 34.50, 34.55]))
    ax.tick_params(length=1.8, width=0.4, labelsize=FS_TICK, color='0.5', pad=1.5)
    ax.tick_params(labelbottom=show_lon, labelleft=show_lat)
    if show_lon:
        ax.xaxis.set_major_formatter(lon_fmt)
        for lab in ax.get_xticklabels():
            lab.set_rotation(30); lab.set_ha('right')
    if show_lat:
        ax.yaxis.set_major_formatter(lat_fmt)
    ax.set_title(name, fontsize=FS_TITLE, fontweight='bold', pad=3)
    cax = fig.add_axes([x + MW + CB_GAP, b + 0.10 * MH, CB_W, 0.80 * MH])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(unit, fontsize=FS_CB, labelpad=1); cb.ax.tick_params(labelsize=FS_TICK, width=0.4, length=1.8)
    cb.outline.set_linewidth(0.5)
    if deco:
        _scalebar(ax); _north(ax)
    return ax

# Trait column (left)
add_map(LX, 0.645, chl, 'YlGn',   0, 80,  'Chlorophyll Content', r'$\mu$g cm$^{-2}$',       show_lat=True, deco=True)
add_map(LX, 0.353, lma, 'YlOrBr', 0, 250, 'Leaf Mass per Area',  r'g m$^{-2}$',             show_lat=True)
add_map(LX, 0.061, lwc, 'Blues',  0, 1.5, 'Leaf Water Content',  r'g cm$^{-2}$ ($\times$10$^{-2}$)', show_lat=True, show_lon=True)
# Flux column (right) — mirrors the trait column's full height
add_map(RX, 0.553, gpp, 'viridis', 0, 12,  'GPP',        r'$\mu$mol CO$_2$ m$^{-2}$ s$^{-1}$', deco=True)
add_map(RX, 0.153, sif, 'plasma',  0, 1.5, 'SIF$_{740}$', r'mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$',  show_lon=True)

# ---------------- reflectance spectrum (symmetric top-left of centre) ----------------
rds = xr.open_dataset(REFL_FILE, decode_times=False)
wl = rds['wavelength'].values
refl = rds['reflectance'].isel(time=0).transpose('lat', 'lon', 'wavelength').values
refl = np.where(refl <= -9000, np.nan, refl)
ATM = [(1340, 1450), (1800, 1950)]
gap = np.zeros_like(wl, bool)
for lo, hi in ATM:
    gap |= (wl >= lo) & (wl <= hi)
sp = fig.add_axes([0.3125, 0.725, 0.175, 0.150])   # centre 0.400
pcol = {2: '#2E7D32', 3: '#1565C0', 4: '#C62828'}; pname = {2: 'PFT 1', 3: 'PFT 2', 4: 'PFT 3'}
for lo, hi in ATM:
    sp.axvspan(lo, hi, color='0.88', zorder=0)
sp.axvspan(636, 673, color='#C62828', alpha=0.10, zorder=0)
sp.axvspan(851, 879, color='#7B4A12', alpha=0.10, zorder=0)
mus = {p: np.nanmean(refl[pft == p], 0) for p in [2, 3, 4]}
sds = {p: np.nanstd(refl[pft == p], 0) for p in [2, 3, 4]}
grand = np.nanmean(np.vstack([mus[p] for p in [2, 3, 4]]), 0)
K = 4.0  # schematic: exaggerate inter-PFT differences (real deviations x4) to highlight
         # the distinct spectral features that make trait retrieval possible
for p in [2, 3, 4]:
    ex = np.clip(grand + K * (mus[p] - grand), 0.005, 0.6)
    band = 0.25 * sds[p]
    sp.fill_between(wl, ex - band, ex + band, where=~gap, color=pcol[p], alpha=0.12, lw=0)
    sp.plot(wl, np.where(gap, np.nan, ex), color=pcol[p], lw=1.6, label=pname[p])
sp.set_xlim(wl.min(), wl.max()); sp.set_ylim(0, 0.45)
sp.set_xlabel('Wavelength (nm)', fontsize=FS_CB, fontweight='bold', labelpad=1)
sp.set_ylabel('Reflectance', fontsize=FS_CB, fontweight='bold', labelpad=1)
sp.tick_params(labelsize=FS_TICK, width=0.4, length=1.8)
sp.grid(True, alpha=0.25, lw=0.3)
sp.legend(fontsize=5.5, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.0),
          frameon=True, facecolor='white', edgecolor='0.8', framealpha=0.85,
          columnspacing=1.0, handlelength=1.1, handletextpad=0.4, borderpad=0.3)

# ---------------- LiDAR point cloud (symmetric top-right of centre) ----------------
rng = np.random.default_rng(0)
n = 1100; cx = rng.uniform(0, 1, n); cy = rng.uniform(0, 1, n)
ht = (0.5 + 0.5 * np.sin(6 * cx) * np.cos(6 * cy)) * rng.uniform(0.35, 1.0, n) * 12.0
lax = fig.add_axes([0.515, 0.710, 0.170, 0.180], projection='3d')   # centre 0.600
scat = lax.scatter(cx, cy, ht, c=ht, cmap='viridis', vmin=0, vmax=12, s=3.0, depthshade=True, edgecolors='none')
lax.view_init(elev=24, azim=-60)
lax.set_xlim(0, 1); lax.set_ylim(0, 1); lax.set_zlim(0, 12)
lax.set_xticks([0, 1]); lax.set_yticks([0, 1]); lax.set_zticks([])
lax.set_xticklabels([]); lax.set_yticklabels([])
lax.set_xlabel('x', fontsize=6, labelpad=-24); lax.set_ylabel('y', fontsize=6, labelpad=-24)
for axis in (lax.xaxis, lax.yaxis, lax.zaxis):
    axis.pane.set_facecolor('white'); axis.pane.set_alpha(0.25); axis.pane.set_edgecolor('0.75')
lax.grid(True)
try:
    lax.set_box_aspect((1, 1, 0.55), zoom=1.25)
except TypeError:
    lax.set_box_aspect((1, 1, 0.55))
lcax = fig.add_axes([0.690, 0.732, 0.011, 0.095])
lcb = fig.colorbar(scat, cax=lcax)
lcb.set_label('Canopy height (m)', fontsize=FS_CB, labelpad=1)
lcb.ax.tick_params(labelsize=FS_TICK, width=0.4, length=1.8); lcb.outline.set_linewidth(0.5)

# ---------------- ERA5 climograph (bottom) ----------------
ea = fig.add_axes([0.4225, 0.020, 0.155, 0.110])
xpos = np.arange(len(ERA_MONTHS))
ea.bar(xpos, ERA_P, width=0.62, color=PBLUE, alpha=0.65)
ea.set_ylabel('P (mm)', color=PBLUE, fontsize=FS_CB, labelpad=1); ea.set_ylim(0, max(ERA_P) * 1.25)
ea.set_xticks(xpos); ea.set_xticklabels(ERA_MONTHS, fontsize=FS_TICK)
ea.tick_params(labelsize=FS_TICK, width=0.4, length=1.8, colors=PBLUE); ea.tick_params(axis='x', colors='black')
ea2 = ea.twinx()
ea2.plot(xpos, ERA_T, color=TRED, marker='o', ms=2.6, lw=1.4)
ea2.set_ylabel('T (°C)', color=TRED, fontsize=FS_CB, labelpad=1)
ea2.tick_params(labelsize=FS_TICK, width=0.4, length=1.8, colors=TRED)
ea.spines['top'].set_visible(False); ea.spines['right'].set_visible(False)
ea.spines['left'].set_color(PBLUE); ea.spines['left'].set_linewidth(0.8)
ea2.spines['top'].set_visible(False); ea2.spines['left'].set_visible(False)
ea2.spines['right'].set_color(TRED); ea2.spines['right'].set_linewidth(0.8)

# ---------------- schematic overlay ----------------
def tag(x, y, t, fs=FS_TAG):
    ov.text(x, y, t, ha='center', va='center', fontsize=fs, fontweight='bold', color=LABEL_TX,
            zorder=6, bbox=dict(boxstyle='round,pad=0.4', fc=LABEL, ec='none'))

def hub(x, y, t, rx=0.066):
    ov.add_patch(Ellipse((x, y), 2 * rx, 2 * rx * FASP, fc=TEAL, ec='none', zorder=4))
    ov.text(x, y, t, ha='center', va='center', color='white', fontsize=FS_HUB, fontweight='bold', zorder=5)

def pbox(x, y, t, w=0.110, h=0.115):
    ov.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0.006',
                                fc=TEAL, ec='none', zorder=4))
    ov.text(x, y, t, ha='center', va='center', color='white', fontsize=FS_BOX, fontweight='bold', zorder=5)

def badge(x, y, n):
    ov.add_patch(Ellipse((x, y), 0.026, 0.026 * FASP, fc=TEAL, ec='white', lw=1.3, zorder=10))
    ov.text(x, y, str(n), ha='center', va='center', color='white', fontsize=8, fontweight='bold', zorder=11)

def arr(p0, p1):
    ov.add_patch(FancyArrowPatch(p0, p1, arrowstyle='-|>', mutation_scale=12, lw=2.2, color=TEAL, zorder=3))

def ln(xs, ys):
    ov.plot(xs, ys, color=TEAL, lw=2.2, zorder=3, solid_capstyle='round')

def brace(x, y0, y1, tick):
    ov.plot([x + tick, x, x, x + tick], [y1, y1, y0, y0], color=TEAL, lw=2.2, zorder=3,
            solid_capstyle='round', solid_joinstyle='round')

CY = 0.45; HT = CY + 0.066 * FASP; HB = CY - 0.066 * FASP
hub(0.5, CY, 'CliMA\nLand')
pbox(0.335, CY, 'Inverse\nmodeling')
pbox(0.665, CY, 'Forward\nmodeling')
tag(0.400, 0.945, 'AVIRIS-NG')
tag(0.600, 0.945, 'LiDAR clumping index', fs=8)
tag(0.108, 0.945, 'Traits'); tag(0.850, 0.945, 'Fluxes')
ov.text(0.500, 0.780, '+', ha='center', va='center', fontsize=17, fontweight='bold', color=TEAL, zorder=6)

# inverse <-> hub : two arrows, gapped from both the box and the hub
arr((0.397, CY + 0.019), (0.430, CY + 0.019))
arr((0.430, CY - 0.019), (0.397, CY - 0.019))
# hub -> forward (gapped)
arr((0.572, CY), (0.604, CY))
# (1) AVIRIS spectrum + LiDAR -> hub ; symmetric about x=0.5, rounded corners
ov.plot([0.455, 0.455, 0.545, 0.545], [0.665, 0.595, 0.595, 0.665], color=TEAL, lw=2.2,
        zorder=3, solid_capstyle='round', solid_joinstyle='round')
arr((0.5, 0.595), (0.5, HT)); badge(0.5, 0.595, 1)
# (2) traits -> inverse
brace(0.248, BOT, TOP, -0.010); arr((0.248, CY), (0.270, CY)); badge(0.335, CY + 0.0575, 2)
# (3) ERA5 -> hub
arr((0.5, 0.200), (0.5, HB)); badge(0.5, 0.290, 3)
# (4) forward -> fluxes
brace(0.752, 0.153, 0.768, 0.010); arr((0.752, CY), (0.730, CY)); badge(0.665, CY + 0.0575, 4)
tag(0.5, 0.170, 'ERA5 reanalysis')

for ext in ('pdf', 'png'):
    fig.savefig(OUT / f'figure2_schematic.{ext}', dpi=600, facecolor='white', bbox_inches='tight', pad_inches=0.03)
plt.close(fig)
print('saved', OUT / 'figure2_schematic.pdf', '+ .png (600 dpi)')
