"""
Figure 2 — CliMA Land workflow schematic (publication build).

Single, grid-based, full-page (180 mm) figure:
  LEFT  : Traits  (Chlorophyll, Leaf Mass per Area, Leaf Water Content)
  CENTER: AVIRIS-NG spectrum + LiDAR point cloud -> CliMA Land hub;
          Inverse modeling (left) / Forward modeling (right); ERA5 (bottom)
  RIGHT : Fluxes  (GPP, SIF740)

Maps are true geographic axes (lat/lon graticule, scale bar + north arrow on the
top map of each column), cropped to the preserve polygon. Exports vector PDF
(map layers rasterized at 600 dpi, everything else vector) + 600-dpi PNG.
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
TEAL = '#2b6c8f'; LABEL = '#cfe2f3'; LABEL_TX = '#123449'; SUB = '#333333'
FS_HUB, FS_BOX, FS_TAG, FS_TITLE, FS_SUB, FS_CB, FS_TICK = 15, 10, 9, 8.5, 7, 7.5, 6.5

BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')
TRAIT = BASE / 'traits' / 'datasets' / 'clima_fit_prescribed_lai_ci'
FLUX = BASE / 'fitting' / 'fitted_prescribed_lai_ci'
PFT_FILE = BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'
REFL_FILE = BASE / 'aviris_dangermond_time_00' / 'output_clipped_reduced.nc'
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
lwc = mp(xr.open_dataset(TRAIT / 'lwc_aviris_dangermond_clima_fit_reg.nc')['lwc'].mean('time').values * 18.0)  # mol m-2 -> g m-2
g, s = [], []
for t in range(13):
    f = FLUX / f'shift_fluxes_day_{t:02d}_clima_fit_reg_jmax.nc'
    if f.exists():
        d = xr.open_dataset(f, decode_times=False)
        g.append(orient(d['gpp'].squeeze().values, nlat, nlon))
        s.append(orient(d['sif740'].squeeze().values, nlat, nlon))
gpp = mp(np.nanmean(np.dstack(g), 2)); sif = mp(np.nanmean(np.dstack(s), 2))

vlon = LON[mask]; vlat = LAT[mask]
XLIM = (vlon.min() - 0.006, vlon.max() + 0.006)
YLIM = (vlat.min() - 0.006, vlat.max() + 0.006)
WX, HY = XLIM[1] - XLIM[0], YLIM[1] - YLIM[0]
MEANLAT = float(np.mean(vlat)); GEO_ASPECT = 1 / np.cos(np.deg2rad(MEANLAT))

# ---------------- figure & layout ----------------
fig = plt.figure(figsize=(7.09, 4.72), dpi=200)
FASP = 7.09 / 4.72
ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, 1); ov.set_ylim(0, 1); ov.axis('off')

MW, MH = 0.160, 0.205
LX, RX = 0.050, 0.775
CB_W, CB_GAP = 0.012, 0.006
lon_fmt = FuncFormatter(lambda v, _: f"{abs(v):.2f}°W")
lat_fmt = FuncFormatter(lambda v, _: f"{v:.2f}°N")

def _scalebar(ax, km=3):
    dlon = km / (111.320 * np.cos(np.deg2rad(MEANLAT)))
    x0 = XLIM[0] + 0.09 * WX; y0 = YLIM[0] + 0.07 * HY
    ax.plot([x0, x0 + dlon], [y0, y0], color='black', lw=1.4, solid_capstyle='butt', zorder=6)
    ax.text(x0 + dlon / 2, y0 + 0.02 * HY, f'{km} km', ha='center', va='bottom', fontsize=5.5, zorder=6)

def _north(ax):
    x = XLIM[1] - 0.09 * WX; y0 = YLIM[1] - 0.34 * HY; dy = 0.15 * HY
    ax.annotate('', xy=(x, y0 + dy), xytext=(x, y0),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=1.0), zorder=6)
    ax.text(x, y0 + dy + 0.008 * HY, 'N', ha='center', va='bottom', fontsize=5.5, fontweight='bold', zorder=6)

def add_map(x, b, data, cmap, vmin, vmax, name, unit, letter,
            show_lat=False, show_lon=False, deco=False):
    ax = fig.add_axes([x, b, MW, MH])
    cm = plt.get_cmap(cmap).copy(); cm.set_bad(alpha=0)
    im = ax.pcolormesh(LON, LAT, data, cmap=cm, vmin=vmin, vmax=vmax, shading='auto', rasterized=True)
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM); ax.set_aspect(GEO_ASPECT)
    ax.grid(True, color='0.7', lw=0.3, alpha=0.6)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.xaxis.set_major_locator(FixedLocator([-120.50, -120.45, -120.40]))
    ax.yaxis.set_major_locator(FixedLocator([34.45, 34.50, 34.55]))
    ax.tick_params(length=2, width=0.5, labelsize=FS_TICK, color='0.5')
    ax.tick_params(labelbottom=show_lon, labelleft=show_lat)
    if show_lon:
        ax.xaxis.set_major_formatter(lon_fmt)
        for lab in ax.get_xticklabels():
            lab.set_rotation(30); lab.set_ha('right')
    if show_lat:
        ax.yaxis.set_major_formatter(lat_fmt)
    m_, s_ = np.nanmean(data), np.nanstd(data)
    ax.set_title(name, fontsize=FS_TITLE, fontweight='bold', pad=13)
    ax.text(0.5, 1.015, f'{m_:.1f} ± {s_:.1f}', transform=ax.transAxes,
            ha='center', va='bottom', fontsize=FS_SUB, color=SUB)
    ax.text(0.035, 0.965, f'({letter})', transform=ax.transAxes, ha='left', va='top',
            fontsize=FS_TITLE, fontweight='bold', zorder=7,
            bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none', alpha=0.75))
    cax = fig.add_axes([x + MW + CB_GAP, b + 0.10 * MH, CB_W, 0.80 * MH])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(unit, fontsize=FS_CB); cb.ax.tick_params(labelsize=FS_TICK, width=0.5, length=2)
    cb.outline.set_linewidth(0.5)
    if deco:
        _scalebar(ax); _north(ax)
    return ax

# Trait column (left)
add_map(LX, 0.665, chl, 'YlGn',   0, 80,  'Chlorophyll Content', r'$\mu$g cm$^{-2}$', 'a', show_lat=True, deco=True)
add_map(LX, 0.395, lma, 'YlOrBr', 0, 250, 'Leaf Mass per Area',  r'g m$^{-2}$',       'b', show_lat=True)
add_map(LX, 0.125, lwc, 'Blues',  0, 90,  'Leaf Water Content',  r'g m$^{-2}$',       'c', show_lat=True, show_lon=True)
# Flux column (right; centered on trait column)
add_map(RX, 0.5275, gpp, 'viridis', 0, 12, 'GPP',        r'$\mu$mol CO$_2$ m$^{-2}$ s$^{-1}$', 'd', deco=True)
add_map(RX, 0.2625, sif, 'plasma',  0, 2,  'SIF$_{740}$', r'mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$',  'e', show_lon=True)

# ---------------- reflectance spectrum ----------------
rds = xr.open_dataset(REFL_FILE, decode_times=False)
wl = rds['wavelength'].values
refl = rds['reflectance'].isel(time=0).transpose('lat', 'lon', 'wavelength').values
refl = np.where(refl <= -9000, np.nan, refl)
ATM = [(1340, 1450), (1800, 1950)]
gap = np.zeros_like(wl, bool)
for lo, hi in ATM:
    gap |= (wl >= lo) & (wl <= hi)
sp = fig.add_axes([0.3925, 0.700, 0.215, 0.150])
pcol = {2: '#2E7D32', 3: '#1565C0', 4: '#C62828'}; pname = {2: 'PFT 1', 3: 'PFT 2', 4: 'PFT 3'}
for lo, hi in ATM:
    sp.axvspan(lo, hi, color='0.88', zorder=0)
sp.axvspan(636, 673, color='#C62828', alpha=0.10, zorder=0)
sp.axvspan(851, 879, color='#7B4A12', alpha=0.10, zorder=0)
for p in [2, 3, 4]:
    vals = refl[pft == p]
    mu = np.nanmean(vals, 0); sd = np.nanstd(vals, 0)
    sp.fill_between(wl, mu - sd, mu + sd, where=~gap, color=pcol[p], alpha=0.15, lw=0)
    sp.plot(wl, np.where(gap, np.nan, mu), color=pcol[p], lw=1.5, label=pname[p])
sp.set_xlim(wl.min(), wl.max()); sp.set_ylim(0, 0.45)
sp.set_xlabel('Wavelength (nm)', fontsize=FS_CB, fontweight='bold', labelpad=1)
sp.set_ylabel('Reflectance', fontsize=FS_CB, fontweight='bold', labelpad=1)
sp.set_title('Mean Reflectance Spectrum by PFT', fontsize=FS_SUB + 0.5, fontweight='bold', pad=8)
sp.text(0.02, 0.965, '(f)', transform=sp.transAxes, ha='left', va='top', fontsize=FS_TITLE, fontweight='bold', zorder=7, bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none', alpha=0.75))
sp.tick_params(labelsize=FS_TICK, width=0.5, length=2)
sp.grid(True, alpha=0.25, lw=0.3)
sp.legend(fontsize=6, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.0),
          frameon=True, facecolor='white', edgecolor='0.8', framealpha=0.85,
          columnspacing=1.1, handlelength=1.2, handletextpad=0.4, borderpad=0.3)

# ---------------- LiDAR point-cloud icon (3D) ----------------
rng = np.random.default_rng(0)
n = 500; cx = rng.uniform(0, 1, n); cy = rng.uniform(0, 1, n)
cz = (0.5 + 0.5 * np.sin(6 * cx) * np.cos(6 * cy)) * rng.uniform(0.4, 1.0, n)
lax = fig.add_axes([0.628, 0.720, 0.110, 0.120], projection='3d')
lax.scatter(cx, cy, cz, c=cz, cmap='viridis', s=1.6, depthshade=True, edgecolors='none')
lax.view_init(elev=22, azim=-60); lax.set_axis_off()
try:
    lax.set_box_aspect((1, 1, 0.6))
except Exception:
    pass

# ---------------- schematic overlay ----------------
def tag(x, y, t, fs=FS_TAG):
    ov.text(x, y, t, ha='center', va='center', fontsize=fs, fontweight='bold', color=LABEL_TX,
            zorder=6, bbox=dict(boxstyle='round,pad=0.4', fc=LABEL, ec='none'))

def hub(x, y, t, rx=0.072):
    ov.add_patch(Ellipse((x, y), 2 * rx, 2 * rx * FASP, fc=TEAL, ec='none', zorder=4))
    ov.text(x, y, t, ha='center', va='center', color='white', fontsize=FS_HUB, fontweight='bold', zorder=5)

def pbox(x, y, t, w=0.115, h=0.115):
    ov.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0.006',
                                fc=TEAL, ec='none', zorder=4))
    ov.text(x, y, t, ha='center', va='center', color='white', fontsize=FS_BOX, fontweight='bold', zorder=5)

def badge(x, y, n):
    ov.add_patch(Ellipse((x, y), 0.026, 0.026 * FASP, fc=TEAL, ec='white', lw=1.2, zorder=8))
    ov.text(x, y, str(n), ha='center', va='center', color='white', fontsize=8, fontweight='bold', zorder=9)

def arr(p0, p1, dbl=False):
    ov.add_patch(FancyArrowPatch(p0, p1, arrowstyle='<|-|>' if dbl else '-|>',
                                 mutation_scale=13, lw=2.0, color=TEAL, zorder=3))

def ln(xs, ys):
    ov.plot(xs, ys, color=TEAL, lw=2.0, zorder=3, solid_capstyle='round')

def brace(x, y0, y1, tick, cx, cy):
    ln([x, x], [y0, y1]); ln([x, x + tick], [y0, y0]); ln([x, x + tick], [y1, y1]); ln([x, cx], [cy, cy])

CY = 0.45; HT = CY + 0.072 * FASP; HB = CY - 0.072 * FASP
hub(0.5, CY, 'CliMA\nLand')
pbox(0.335, CY, 'Inverse\nmodeling')
pbox(0.665, CY, 'Forward\nmodeling')
tag(0.5, 0.945, 'AVIRIS-NG')
tag(0.130, 0.945, 'Traits'); tag(0.855, 0.945, 'Fluxes')
tag(0.688, 0.700, 'LiDAR clumping index', fs=6.5)
tag(0.5, 0.055, 'ERA5 reanalysis')
ov.text(0.612, 0.790, '+', ha='center', va='center', fontsize=16, fontweight='bold', color=TEAL, zorder=6)

arr((0.3925, CY), (0.428, CY), dbl=True)
arr((0.572, CY), (0.6075, CY))
# (1) AVIRIS spectrum + LiDAR -> hub
ln([0.5, 0.5], [0.644, 0.605]); ln([0.683, 0.683], [0.678, 0.605]); ln([0.5, 0.683], [0.605, 0.605])
arr((0.5, 0.605), (0.5, HT)); badge(0.5, 0.605, 1)
# (2) traits -> inverse
brace(0.268, 0.128, 0.868, -0.010, 0.2775, CY); badge(0.268, 0.640, 2)
# (3) ERA5 -> hub
arr((0.5, 0.112), (0.5, HB)); badge(0.5, 0.250, 3)
# (4) forward -> fluxes
brace(0.750, 0.263, 0.732, 0.010, 0.7225, CY); badge(0.750, 0.640, 4)

for ext in ('pdf', 'png'):
    fig.savefig(OUT / f'figure2_schematic.{ext}', dpi=600, facecolor='white', bbox_inches='tight', pad_inches=0.03)
plt.close(fig)
print('saved', OUT / 'figure2_schematic.pdf', '+ .png (600 dpi)')
