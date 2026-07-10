"""
Regenerate the six panels embedded in Figure 2 (the CliMA Land workflow
schematic) as publication-quality thumbnails.

Design choices (so the panels read as a figure, not auto-generated plots):
  - Title = variable name (not "Mean"); Mean/STD shown as a subtitle line.
  - Colorbar label = units only (identity is in the title).
  - Map spines off, cropped tightly to the preserve polygon (no white frame).
  - Unified font (Nimbus Sans, Helvetica-metric) and font sizes across panels.
  - All five maps share one footprint/crop -> identical on-page size.
  - LWC in g cm-2 (x10-2), matching Figures 4 and 5.
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.image as _mpl  # noqa: F401 (ensures backend init)
from pathlib import Path

plt.rcParams.update({
    'font.family': 'Nimbus Sans',
    'font.size': 20,
    'axes.titlesize': 25,
    'axes.titleweight': 'bold',
    'figure.dpi': 150,
    # render super/subscripts in the same sans font (Nimbus lacks unicode ⁻² etc.)
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'Nimbus Sans',
    'mathtext.it': 'Nimbus Sans:italic',
    'mathtext.bf': 'Nimbus Sans:bold',
})

BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')
TRAIT = BASE / 'traits' / 'datasets' / 'clima_fit_prescribed_lai_ci'
FLUX = BASE / 'fitting' / 'fitted_prescribed_lai_ci'
PFT_FILE = BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'
REFL_FILE = BASE / 'aviris_dangermond_time_00' / 'output_clipped_reduced.nc'
OUT = BASE / 'shift_dangermond_trait' / 'figures'
OUT.mkdir(parents=True, exist_ok=True)

FIGSIZE = (6.6, 5.0)          # identical for all five maps
STAT_FS, CBLBL_FS, CBTICK_FS = 16, 20, 20

def orient(arr, nlat, nlon):
    a = np.asarray(arr, float)
    a[a <= -9000] = np.nan
    return a.T if a.shape == (nlon, nlat) else a

# ---- geometry & shared crop footprint ----
chl_ds = xr.open_dataset(TRAIT / 'chl_aviris_dangermond_clima_fit_reg.nc')
nlat, nlon = len(chl_ds['lat']), len(chl_ds['lon'])
pft_map = orient(xr.open_dataset(PFT_FILE)['Band1'].values, nlat, nlon)
pft_mask = np.isin(pft_map, [2, 3, 4])
_rows, _cols = np.any(pft_mask, axis=1), np.any(pft_mask, axis=0)
r0, r1 = np.where(_rows)[0][[0, -1]]
c0, c1 = np.where(_cols)[0][[0, -1]]
PAD = 6
CROP = dict(xlim=(c0 - PAD, c1 + PAD), ylim=(r0 - PAD, r1 + PAD))

def mask_pft(a):
    return np.where(pft_mask, a, np.nan)

def save_map(data, cmap, vmin, vmax, name, unit, fname):
    m, s = np.nanmean(data), np.nanstd(data)
    cm = plt.get_cmap(cmap).copy(); cm.set_bad(alpha=0.0)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    im = ax.imshow(data, cmap=cm, vmin=vmin, vmax=vmax, origin='lower', aspect='equal')
    ax.set_xlim(*CROP['xlim']); ax.set_ylim(*CROP['ylim'])
    ax.set_axis_off()
    ax.set_title(name, pad=26)
    ax.text(0.5, 1.01, f'Mean {m:.1f} · STD {s:.1f}', transform=ax.transAxes,
            ha='center', va='bottom', fontsize=STAT_FS, color='0.30')
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label(unit, fontsize=CBLBL_FS)
    cbar.ax.tick_params(labelsize=CBTICK_FS)
    cbar.outline.set_linewidth(0.6)
    fig.savefig(OUT / fname, dpi=300, facecolor='white', bbox_inches='tight', pad_inches=0.15)
    plt.close(fig)
    print(f'  saved {fname}  (mean={m:.3f} std={s:.3f})')

# ---------------- TRAIT MAPS ----------------
print('Trait maps:')
chl = mask_pft(chl_ds['chl'].mean(dim='time').values)
lma = mask_pft(xr.open_dataset(TRAIT / 'lma_aviris_dangermond_clima_fit_reg.nc')['lma'].mean(dim='time').values * 1e4)
lwc_raw = xr.open_dataset(TRAIT / 'lwc_aviris_dangermond_clima_fit_reg.nc')['lwc'].mean(dim='time').values
lwc_gcm2 = lwc_raw * 0.0018            # mol m-2 -> g cm-2
lwc_disp = mask_pft(lwc_gcm2 * 100.0)  # display as x10-2 (matches Figs 4/5)
save_map(chl,      'YlGn',   0, 80,  'Chlorophyll Content', r'$\mu$g cm$^{-2}$', 'fig2_trait_CHL.png')
save_map(lma,      'YlOrBr', 0, 250, 'Leaf Mass per Area',  r'g m$^{-2}$',   'fig2_trait_LMA.png')
save_map(lwc_disp, 'Blues',  0, 1.5, 'Leaf Water Content',  r'g cm$^{-2}$ ($\times$10$^{-2}$)', 'fig2_trait_LWC.png')

# ---------------- FLUX MAPS (season mean over 13 days) ----------------
print('Flux maps:')
gpp_stack, sif_stack = [], []
for t in range(13):
    f = FLUX / f'shift_fluxes_day_{t:02d}_clima_fit_reg_jmax.nc'
    if not f.exists():
        continue
    ds = xr.open_dataset(f, decode_times=False)
    gpp_stack.append(orient(ds['gpp'].squeeze().values, nlat, nlon))
    sif_stack.append(orient(ds['sif740'].squeeze().values, nlat, nlon))
gpp = mask_pft(np.nanmean(np.dstack(gpp_stack), axis=2))
sif = mask_pft(np.nanmean(np.dstack(sif_stack), axis=2))
save_map(gpp, 'YlGn', 0, 12, 'GPP',        r'$\mu$mol CO$_2$ m$^{-2}$ s$^{-1}$',   'fig2_flux_GPP.png')
save_map(sif, 'YlGn', 0, 2,  'SIF$_{740}$', r'mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$',  'fig2_flux_SIF.png')

# ---------------- REFLECTANCE SPECTRUM ----------------
print('Reflectance spectrum:')
rds = xr.open_dataset(REFL_FILE, decode_times=False)
wl = rds['wavelength'].values
refl = rds['reflectance'].isel(time=0).transpose('lat', 'lon', 'wavelength').values
refl = np.where(refl <= -9000, np.nan, refl)
pft_colors = {2: '#2E7D32', 3: '#1976D2', 4: '#D32F2F'}
pft_names = {2: 'PFT 1', 3: 'PFT 2', 4: 'PFT 3'}
fig, ax = plt.subplots(figsize=(8.3, 5.0))
ax.axvspan(636, 673, color='red', alpha=0.12, label='LAI red range')
ax.axvspan(851, 879, color='saddlebrown', alpha=0.12, label='LAI NIR range')
for lo, hi in [(1340, 1440), (1800, 1950)]:
    ax.axvspan(lo, hi, color='0.85', alpha=0.6)
for p in [2, 3, 4]:
    vals = refl[pft_map == p]
    mu, sd = np.nanmean(vals, axis=0), np.nanstd(vals, axis=0)
    ax.plot(wl, mu, color=pft_colors[p], lw=2.2, label=pft_names[p])
    ax.fill_between(wl, mu - sd, mu + sd, color=pft_colors[p], alpha=0.18)
ax.set_xlabel('Wavelength (nm)', fontsize=CBLBL_FS, fontweight='bold')
ax.set_ylabel('Reflectance', fontsize=CBLBL_FS, fontweight='bold')
ax.set_title('Mean Reflectance Spectrum for PFTs', fontsize=22)
ax.tick_params(axis='both', labelsize=CBTICK_FS)
ax.set_xlim(wl.min(), wl.max())
ax.legend(fontsize=15, framealpha=0.9, loc='upper right')
ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
fig.savefig(OUT / 'fig2_spectrum.png', dpi=300, facecolor='white', bbox_inches='tight', pad_inches=0.15)
plt.close(fig)
print('  saved fig2_spectrum.png')
print('DONE ->', OUT)
