"""
Regenerate the six raster panels embedded in Figure 2 (the CliMA Land workflow
schematic) with print-legible lettering. Editor item 7: the Trait and Flux panel
lettering is too small.

Panels (each saved as its own PNG to drop back into the PowerPoint frames):
  Traits column : CHL, LMA, LWC time-mean maps          (aspect ~1.33:1)
  Fluxes column : GPP, SIF740 season-mean maps           (aspect ~1.33:1)
  Top           : mean AVIRIS-NG reflectance spectrum    (aspect ~1.66:1)

Fonts are enlarged; lat/lon axis ticks are dropped (not load-bearing in a
schematic thumbnail) so the colorbar numbers + title carry the information.
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')
TRAIT = BASE / 'traits' / 'datasets' / 'clima_fit_prescribed_lai_ci'
FLUX = BASE / 'fitting' / 'fitted_prescribed_lai_ci'
PFT_FILE = BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'
REFL_FILE = BASE / 'aviris_dangermond_time_00' / 'output_clipped_reduced.nc'
OUT = BASE / 'shift_dangermond_trait' / 'figures'
OUT.mkdir(parents=True, exist_ok=True)

# ---- enlarged fonts ----
TITLE, CBLBL, CBTICK, STAT = 30, 20, 22, 22

def orient(arr, nlat, nlon):
    """Return arr as (nlat, nlon)."""
    a = np.asarray(arr, float)
    a[a <= -9000] = np.nan
    if a.shape == (nlon, nlat):
        a = a.T
    return a

def save_map(data, cmap, vmin, vmax, label, unit, fname, scale_title=None):
    m, s = np.nanmean(data), np.nanstd(data)
    fig, ax = plt.subplots(figsize=(7, 5.25), dpi=150)
    im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower', aspect='equal')
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(1.2)
    ax.set_title('Mean', fontsize=TITLE, fontweight='bold', pad=10)
    ax.text(0.97, 0.97, f'Mean: {m:.2f}\nSTD: {s:.2f}', transform=ax.transAxes,
            fontsize=STAT, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'))
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label(f'{label} ({unit})', fontsize=CBLBL, fontweight='bold')
    cbar.ax.tick_params(labelsize=CBTICK)
    if scale_title:
        cbar.ax.set_title(scale_title, fontsize=18, pad=8)
    plt.tight_layout()
    fig.savefig(OUT / fname, dpi=300, facecolor='white', bbox_inches='tight', pad_inches=0.4)
    plt.close(fig)
    print(f'  saved {fname}  (mean={m:.3f} std={s:.3f})')

# ---------- geometry from a trait file ----------
chl_ds = xr.open_dataset(TRAIT / 'chl_aviris_dangermond_clima_fit_reg.nc')
lats = chl_ds['lat'].values; lons = chl_ds['lon'].values
nlat, nlon = len(lats), len(lons)
pft_map = orient(xr.open_dataset(PFT_FILE)['Band1'].values, nlat, nlon)
pft_mask = np.isin(pft_map, [2, 3, 4])

def mask_pft(a):
    return np.where(pft_mask, a, np.nan)

# ---------- TRAIT MAPS ----------
print('Trait maps:')
chl = mask_pft(chl_ds['chl'].mean(dim='time').values)
lma = mask_pft(xr.open_dataset(TRAIT / 'lma_aviris_dangermond_clima_fit_reg.nc')['lma'].mean(dim='time').values * 1e4)
lwc_raw = xr.open_dataset(TRAIT / 'lwc_aviris_dangermond_clima_fit_reg.nc')['lwc'].mean(dim='time').values
lwc = mask_pft(lwc_raw)  # mol/m^2 as in the original schematic panel
save_map(chl, 'YlGn', 0, 80, 'Chlorophyll Content', 'µg cm⁻²', 'fig2_trait_CHL.png')
save_map(lma, 'YlOrBr', 0, 250, 'Leaf Mass per Area', 'g m⁻²', 'fig2_trait_LMA.png')
save_map(lwc, 'Blues', 0, 10, 'Leaf Water Content', 'mol m⁻²', 'fig2_trait_LWC.png')

# ---------- FLUX MAPS (season mean over 13 days) ----------
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
save_map(gpp, 'YlGn', 0, 12, 'GPP', 'µmol CO₂ m⁻² s⁻¹', 'fig2_flux_GPP.png')
save_map(sif, 'YlGn', 0, 2, 'SIF$_{740nm}$', 'mW m⁻² sr⁻¹ nm⁻¹', 'fig2_flux_SIF.png')

# ---------- REFLECTANCE SPECTRUM ----------
print('Reflectance spectrum:')
rds = xr.open_dataset(REFL_FILE, decode_times=False)
wl = rds['wavelength'].values
refl = rds['reflectance'].isel(time=0).transpose('lat', 'lon', 'wavelength').values  # (lat, lon, wl)
refl = np.where(refl <= -9000, np.nan, refl)

pft_colors = {2: '#2E7D32', 3: '#1976D2', 4: '#D32F2F'}
pft_names = {2: 'PFT 1', 3: 'PFT 2', 4: 'PFT 3'}
fig, ax = plt.subplots(figsize=(8.3, 5.0), dpi=150)
# LAI simple-ratio bands
ax.axvspan(636, 673, color='red', alpha=0.12, label='LAI red range')
ax.axvspan(851, 879, color='saddlebrown', alpha=0.12, label='LAI NIR range')
# atmospheric water-vapour gap
for lo, hi in [(1340, 1440), (1800, 1950)]:
    ax.axvspan(lo, hi, color='0.85', alpha=0.6)
for p in [2, 3, 4]:
    m = pft_map == p
    vals = refl[m]                                  # (Npix, nwl)
    mu = np.nanmean(vals, axis=0)
    sd = np.nanstd(vals, axis=0)
    ax.plot(wl, mu, color=pft_colors[p], lw=2.2, label=pft_names[p])
    ax.fill_between(wl, mu - sd, mu + sd, color=pft_colors[p], alpha=0.18)
ax.set_xlabel('Wavelength (nm)', fontsize=CBLBL, fontweight='bold')
ax.set_ylabel('Reflectance', fontsize=CBLBL, fontweight='bold')
ax.set_title('Mean Reflectance Spectrum for PFTs', fontsize=24, fontweight='bold', pad=10)
ax.tick_params(axis='both', labelsize=CBTICK)
ax.set_xlim(wl.min(), wl.max())
ax.legend(fontsize=16, framealpha=0.9, loc='upper right')
ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
plt.tight_layout()
fig.savefig(OUT / 'fig2_spectrum.png', dpi=300, facecolor='white', bbox_inches='tight')
plt.close(fig)
print('  saved fig2_spectrum.png')
print('DONE ->', OUT)
