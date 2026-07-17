"""
Figure 6 — Trait- vs PFT-based flux comparison (remapcon, RELATIVE SIF).

Faithful port of the accepted-manuscript Figure 6 generator
(review/figures_5m_remapcon_traits_fixed.ipynb, cell "Figure 6 remake … RELATIVE
SIF"), SPLIT into two print-legible figures (Ecosphere editor item 9) so each
panel — and the panel (a)/(b) Δ-distribution insets — is large enough to read at
the 18 cm print width. Same data/values as the single accepted figure.

  Figure 6a (figure6a_difference_maps): the two difference maps
    (a) GPP Trait-PFT spatial difference (30 m) + per-PFT ΔGPP inset
    (b) SIF740 Trait-PFT spatial difference (30 m, mW m-2 sr-1 nm-1) + per-PFT ΔSIF inset
  Figure 6b (figure6b_time_series): the two time series
    (c) GPP time series (Trait vs PFT, 95% CI)
    (d) RELATIVE SIF740 (%), twin axis: model (left) vs TROPOMI (right);
        TROPOMI ±4-day scatter + Gaussian-smoothed ±SEM. R² on relative SIF (0.69/0.68).

Output: figures/figure6a_difference_maps.{png,pdf}, figures/figure6b_time_series.{png,pdf}
"""
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.ndimage import gaussian_filter1d as _gf1d
from pathlib import Path
from scipy import stats as scipy_stats
from scipy.stats import linregress

# Publication style: match the other figures (Nimbus Sans) and embed TrueType
# fonts in vector output (pdf.fonttype 42) so the PDF/EPS carry no Type-3 fonts.
plt.rcParams.update({
    'font.family': 'Nimbus Sans',
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Nimbus Sans',
    'mathtext.it': 'Nimbus Sans:italic', 'mathtext.bf': 'Nimbus Sans:bold',
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

_BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')

# OLD runs (GPP + absolute SIF, panels a, b, c)
_TR_OLD  = str(_BASE / 'results_remapcon' / 'trait_based' /
               'shift_fluxes_day_{:02d}_clima_fit_reg_jmax_fixed_005.nc')
_PFT_OLD = str(_BASE / 'results_remapcon' / 'pft_based' /
               'pft_shift_fluxes_day_{:02d}_clima_fit_reg_jmax_fixed_005.nc')
# NEW APAR runs (relative SIF, panel d)
_TR_APAR  = str(_BASE / 'results_remapcon' / 'trait_based_apar' /
                'shift_fluxes_day_{:02d}_clima_fit_reg_jmax_apar.nc')
_PFT_APAR = str(_BASE / 'results_remapcon' / 'pft_based_apar' /
                'pft_shift_fluxes_day_{:02d}_clima_fit_reg_jmax_apar.nc')

_PFT_NC  = _BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'
_TROP_NC = _BASE / 'TROPOMI_dangermond' / 'TROPOMI_SIF740nm-v1.005deg_regrid_Dangermond_tll_clipped.nc'
_FIG_DIR = _BASE / 'shift_dangermond_trait' / 'figures'
_FIG_DIR.mkdir(parents=True, exist_ok=True)

_DATES = pd.to_datetime([
    "2022-02-24", "2022-02-28", "2022-03-08", "2022-03-16", "2022-03-22",
    "2022-04-05", "2022-04-12", "2022-04-20", "2022-04-29",
    "2022-05-03", "2022-05-11", "2022-05-17", "2022-05-29",
])
_ND = len(_DATES)
_DX = _DY = 0.04999984
_WIN = 4
_WIN14 = 14

# ── grid & PFT ────────────────────────────────────────────────────────────────
with xr.open_dataset(_PFT_NC) as _d:
    _pft_map = _d['Band1'].values
with xr.open_dataset(_TR_OLD.format(0), decode_times=False) as _d:
    _slats = _d['lat'].values; _slons = _d['lon'].values
_LON, _LAT = np.meshgrid(_slons, _slats)

# ── Step 1: GPP & absolute SIF difference maps + GPP ts (OLD runs) ────────────
print("Step 1: GPP & SIF difference …")
_gpp_diff_list, _sif_diff_list = [], []
_gpp_tr_ts, _gpp_pft_ts = [], []
for _i in range(_ND):
    with xr.open_dataset(_TR_OLD.format(_i), decode_times=False) as _dt, \
         xr.open_dataset(_PFT_OLD.format(_i), decode_times=False) as _dp:
        _gt = _dt['gpp'].values.squeeze(); _gp = _dp['gpp'].values.squeeze()
        _st = _dt['sif740'].values.squeeze(); _sp = _dp['sif740'].values.squeeze()
    _gpp_diff_list.append(_gt - _gp); _sif_diff_list.append(_st - _sp)
    _gpp_tr_ts.append(float(np.nanmean(_gt))); _gpp_pft_ts.append(float(np.nanmean(_gp)))
_gpp_diff = np.nanmean(_gpp_diff_list, axis=0)
_sif_diff = np.nanmean(_sif_diff_list, axis=0)

# ── Step 2: GPP CI95 ─────────────────────────────────────────────────────────
print("Step 2: GPP CI95 …")
_gpp_tr_ci, _gpp_pft_ci = [], []
for _i in range(_ND):
    with xr.open_dataset(_TR_OLD.format(_i), decode_times=False) as _dt, \
         xr.open_dataset(_PFT_OLD.format(_i), decode_times=False) as _dp:
        _tv = _dt['gpp'].values.flatten(); _tv = _tv[~np.isnan(_tv)]
        _pv = _dp['gpp'].values.flatten(); _pv = _pv[~np.isnan(_pv)]
    _gpp_tr_ci.append(np.std(_tv) / np.sqrt(len(_tv)) * scipy_stats.t.ppf(0.975, len(_tv) - 1))
    _gpp_pft_ci.append(np.std(_pv) / np.sqrt(len(_pv)) * scipy_stats.t.ppf(0.975, len(_pv) - 1))
_gpp_tr_ci = np.array(_gpp_tr_ci); _gpp_pft_ci = np.array(_gpp_pft_ci)

# ── Step 4: spatial-average upscaling helper ─────────────────────────────────
def _to_trop(src_2d, slats, slons, tlats, tlons):
    out = np.full((len(tlats), len(tlons)), np.nan)
    for _i2, _tlat in enumerate(tlats):
        _lm = (slats >= _tlat - _DY / 2) & (slats < _tlat + _DY / 2)
        if not _lm.any():
            continue
        for _j2, _tlon in enumerate(tlons):
            _cm = (slons >= _tlon - _DX / 2) & (slons < _tlon + _DX / 2)
            _p = src_2d[np.ix_(_lm, _cm)].ravel(); _v = _p[np.isfinite(_p)]
            if _v.size > 0:
                out[_i2, _j2] = np.mean(_v)
    return out

print("Step 4: Dangermond TROPOMI pixels + model relative SIF …")
with xr.open_dataset(_TROP_NC) as _td:
    _tlats = _td.lat.values; _tlons = _td.lon.values
with xr.open_dataset(_TR_APAR.format(0), decode_times=False) as _d0:
    _sr0 = _d0['sif_rel_740'].values.squeeze() * 100
_sa0 = _to_trop(_sr0, _slats, _slons, _tlats, _tlons)
_ri, _ci = np.where(np.isfinite(_sa0))

_sr_tr_ts, _sr_pft_ts = [], []
_sr_lo_tr, _sr_hi_tr = [], []
for _i in range(_ND):
    for (_pat, _ts, _lo, _hi) in [(_TR_APAR, _sr_tr_ts, _sr_lo_tr, _sr_hi_tr),
                                   (_PFT_APAR, _sr_pft_ts, [], [])]:
        with xr.open_dataset(_pat.format(_i), decode_times=False) as _ds:
            _sr2d = _ds['sif_rel_740'].values.squeeze() * 100
        _sa = _to_trop(_sr2d, _slats, _slons, _tlats, _tlons)
        _vv = _sa[_ri, _ci]; _vv = _vv[np.isfinite(_vv)]
        if len(_vv) > 1:
            _h = np.std(_vv) / np.sqrt(len(_vv)) * scipy_stats.t.ppf(0.975, len(_vv) - 1)
            _m = float(np.nanmean(_vv)); _ts.append(_m); _lo.append(_m - _h); _hi.append(_m + _h)
        elif len(_vv) == 1:
            _ts.append(float(_vv[0])); _lo.append(np.nan); _hi.append(np.nan)
        else:
            _ts.append(np.nan); _lo.append(np.nan); _hi.append(np.nan)

# ── Step 5: TROPOMI relative SIF (negatives included; ±4d and ±14d) ──────────
print("Step 5: TROPOMI sif_relative …")
with xr.open_dataset(_TROP_NC) as _td:
    _ttime = pd.date_range('2022-02-01', periods=len(_td.time), freq='D')
    _sifrel_raw = _td['sif_relative'].where(_td['sif_relative'] != -999).values

_trop_4d_vals = []
for _tdate in _DATES:
    _msk = ((pd.DatetimeIndex(_ttime) >= _tdate - pd.Timedelta(days=_WIN)) &
            (pd.DatetimeIndex(_ttime) <= _tdate + pd.Timedelta(days=_WIN)))
    _vals = _sifrel_raw[_msk][:, _ri, _ci]
    _trop_4d_vals.append(float(np.nanmean(_vals)) if np.any(np.isfinite(_vals)) else np.nan)
_trop_4d_vals = np.array(_trop_4d_vals)

_trop_14d_mean = []
for _tdate in _DATES:
    _msk = ((pd.DatetimeIndex(_ttime) >= _tdate - pd.Timedelta(days=_WIN14)) &
            (pd.DatetimeIndex(_ttime) <= _tdate + pd.Timedelta(days=_WIN14)))
    _dvals = np.nanmean(_sifrel_raw[_msk][:, _ri, _ci], axis=1); _dvals = _dvals[np.isfinite(_dvals)]
    _trop_14d_mean.append(float(np.mean(_dvals)) if len(_dvals) >= 2 else np.nan)
_trop_14d_mean = np.array(_trop_14d_mean)

_arr_tr = np.array(_sr_tr_ts); _arr_pft = np.array(_sr_pft_ts); _arr_ob = _trop_4d_vals
_ok_tr = np.isfinite(_arr_tr) & np.isfinite(_arr_ob)
_ok_pft = np.isfinite(_arr_pft) & np.isfinite(_arr_ob)
_r2_tr = linregress(_arr_tr[_ok_tr], _arr_ob[_ok_tr]).rvalue ** 2 if _ok_tr.sum() >= 3 else np.nan
_r2_pft = linregress(_arr_pft[_ok_pft], _arr_ob[_ok_pft]).rvalue ** 2 if _ok_pft.sum() >= 3 else np.nan
print(f"Relative SIF R2:  Trait={_r2_tr:.3f}  PFT={_r2_pft:.3f}")

# ── per-PFT density inset (ENLARGED for print legibility, editor item 9) ──────
# Figure 6 is split into 6a (maps) and 6b (time series) so each map is a full
# half-width panel; the inset can therefore be large with big fonts, legible at
# the 18 cm print width.
INSET_SIZE = "46%"
INSET_XLBL, INSET_YLBL, INSET_TICK, INSET_LEG = 19, 18, 17, 17


def _density_inset(ax, data, pft_map, xlabel, xlim=None):
    axins = inset_axes(ax, width=INSET_SIZE, height=INSET_SIZE, loc='lower right',
                       bbox_to_anchor=(0.05, 0.08, 0.93, 0.93),
                       bbox_transform=ax.transAxes, borderpad=0)
    _cols = {2: 'green', 3: 'blue', 4: 'red'}
    _nms = {2: 'PFT 1', 3: 'PFT 2', 4: 'PFT 3'}
    _fd = data.flatten(); _fp = pft_map.flatten()
    _ok = ~np.isnan(_fd) & ~np.isnan(_fp)
    _fd = _fd[_ok]; _fp = _fp[_ok]
    for _pid in [2, 3, 4]:
        _d = _fd[_fp == _pid]
        if len(_d) > 0:
            axins.hist(_d, bins=40, alpha=0.6, density=True, color=_cols[_pid],
                       label=_nms[_pid], edgecolor='black', linewidth=0.5)
    axins.axvline(0, color='red', ls='--', lw=1.5, alpha=0.8)
    axins.axvline(np.mean(_fd), color='black', ls='--', lw=1.5, alpha=0.8)
    axins.set_xlabel(xlabel, fontsize=INSET_XLBL, fontweight='bold')
    if xlim is not None:
        axins.set_xlim(*xlim)
    else:
        _p2, _p98 = np.nanpercentile(_fd, 2), np.nanpercentile(_fd, 98)
        _pad = (_p98 - _p2) * 0.15
        axins.set_xlim(_p2 - _pad, _p98 + _pad)
    axins.set_ylabel('Density', fontsize=INSET_YLBL, fontweight='bold')
    axins.tick_params(labelsize=INSET_TICK)
    axins.locator_params(axis='x', nbins=3)
    axins.legend(fontsize=INSET_LEG, loc='upper right', framealpha=0.9,
                 handlelength=1.1, handletextpad=0.4, borderpad=0.3, labelspacing=0.25)
    axins.grid(True, alpha=0.3, lw=0.5)
    axins.set_facecolor('white'); axins.patch.set_alpha(0.95)
    for sp in axins.spines.values():
        sp.set_edgecolor('black'); sp.set_linewidth(1.5)
    return axins

# Fonts for Figure 6b time series — matched to Figure 6a (panel 24, ticks 22, labels 18)
TS_PANEL, TS_AXLBL, TS_TICK, TS_LEG = 24, 18, 22, 15

# ============ Figure 6a: difference maps — matched to Figure 5b style ============
# Manual add_axes layout so each colorbar is exactly the map height; padded map
# extent so the lower-right is blank white for a transparent Δ-distribution inset;
# 2 lat/lon ticks and NO "Latitude"/"Longitude" words; panel letters outside;
# per-map title + "Difference (Trait − PFT)" suptitle. (Same look as Fig 5b.)
print("Building Figure 6a (maps, Fig-5b style) …")
A_TICK = A_CBTICK = 22; A_CBLBL = 18; A_STAT = 16; A_PANEL = 24; A_TITLE = 23
A_INSET_LBL, A_INSET_TICK, A_INSET_LEG = 13, 11, 8
LON_TICKS, LAT_TICKS = [-120.50, -120.40], [34.45, 34.55]

_lon_r = _slons.max() - _slons.min(); _lat_r = _slats.max() - _slats.min()
A_PAD_R, A_PAD_B = 0.14, 0.22
A_XLIM = (_slons.min(), _slons.max() + A_PAD_R * _lon_r)
A_YLIM = (_slats.min() - A_PAD_B * _lat_r, _slats.max())
A_ASP = (A_XLIM[1] - A_XLIM[0]) / (A_YLIM[1] - A_YLIM[0])

# Element sizes matched to Figure 5b (in INCHES) so per-panel proportions are
# identical. Fig 5b is 18 in wide with 3 maps: left .99, map-cbar gap .108,
# cbar .234, cbar-label .936, panel gap .90, right .63, each map 3.582 in wide.
# Reuse those exact inch values with 2 maps → figure ~12.2 in wide.
_L_in, _GCB_in, _CW_in, _CLAB_in, _GUNIT_in, _RB_in = 0.99, 0.108, 0.234, 0.936, 0.90, 0.63
_MAPW_in = 3.582
A_TOP_IN, A_BOT_IN = 1.05, 0.55
A_W = _L_in + 2 * (_MAPW_in + _GCB_in + _CW_in + _CLAB_in) + _GUNIT_in + _RB_in
_a_map_h_in = _MAPW_in / A_ASP
A_H = A_TOP_IN + _a_map_h_in + A_BOT_IN
A_L, A_GCB, A_CW, A_CLAB, A_GUNIT = (_L_in / A_W, _GCB_in / A_W, _CW_in / A_W,
                                     _CLAB_in / A_W, _GUNIT_in / A_W)
A_MW = _MAPW_in / A_W
A_MH = _a_map_h_in / A_H
A_MAPY = A_BOT_IN / A_H
A_INSET_POS = (0.645, 0.12, 0.325, 0.31)   # anchored to the lower-right white corner

_maps6a = [
    ('gpp', _gpp_diff, -5.0, 5.0, r'$\Delta$ GPP' + '\n' + r'($\mu$mol CO$_2$ m$^{-2}$ s$^{-1}$)',
     r'$\Delta$ GPP', 'GPP', (-10, 10), 2),
    ('sif', _sif_diff, -0.06, 0.06, r'$\Delta$ SIF$_{740}$' + '\n' + r'(mW m$^{-2}$ sr$^{-1}$ nm$^{-1}$)',
     r'$\Delta$ SIF$_{740}$', r'SIF$_{740}$', None, 3),
]
_letters6a = ['a', 'b']
figA = plt.figure(figsize=(A_W, A_H), dpi=150)
for _c, (_nm, _disp, _vmn, _vmx, _cblab, _dlab, _ttl, _ixlim, _dec) in enumerate(_maps6a):
    _x = A_L + _c * (A_MW + A_GCB + A_CW + A_CLAB + A_GUNIT)
    _ax = figA.add_axes([_x, A_MAPY, A_MW, A_MH])
    _im = _ax.pcolormesh(_LON, _LAT, _disp, cmap='RdBu_r', vmin=_vmn, vmax=_vmx, shading='auto', rasterized=True)
    _ax.set_xlim(*A_XLIM); _ax.set_ylim(*A_YLIM)
    _ax.set_xticks(LON_TICKS); _ax.set_yticks(LAT_TICKS)
    _ax.tick_params(labelsize=A_TICK, labelleft=(_c == 0), labelbottom=True)
    _ax.annotate(f'({_letters6a[_c]})', xy=(0, 1), xycoords='axes fraction', xytext=(-2, 4),
                 textcoords='offset points', fontsize=A_PANEL, fontweight='bold', va='bottom', ha='right')
    _ax.text(0.965, 0.96, f"Mean: {np.nanmean(_disp):.{_dec}f}\nSTD: {np.nanstd(_disp):.{_dec}f}",
             transform=_ax.transAxes, fontsize=A_STAT, va='top', ha='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'))
    _ax.set_title(_ttl, fontsize=A_TITLE, fontweight='bold', pad=8)
    _cax = figA.add_axes([_x + A_MW + A_GCB, A_MAPY, A_CW, A_MH])
    _cb = figA.colorbar(_im, cax=_cax)
    _cb.set_label(_cblab, fontsize=A_CBLBL, fontweight='bold')
    _cb.ax.tick_params(labelsize=A_CBTICK)
    if _nm == 'sif':
        _cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    # transparent per-PFT Δ-distribution inset in the lower-right blank white area
    _ix, _iy, _iw, _ih = A_INSET_POS
    _hax = figA.add_axes([_x + _ix * A_MW, A_MAPY + _iy * A_MH, _iw * A_MW, _ih * A_MH])
    _hax.patch.set_alpha(0)
    _fd = _disp.flatten(); _fp = _pft_map.flatten()
    _okf = ~np.isnan(_fd) & ~np.isnan(_fp); _fd = _fd[_okf]; _fp = _fp[_okf]
    for _pid, _col, _lab in [(2, 'green', 'PFT 1'), (3, 'blue', 'PFT 2'), (4, 'red', 'PFT 3')]:
        _dd = _fd[_fp == _pid]
        if len(_dd):
            _hax.hist(_dd, bins=40, alpha=0.6, density=True, color=_col, label=_lab,
                      edgecolor='black', linewidth=0.4)
    _hax.axvline(0, color='red', ls='--', lw=1.5, alpha=0.8)
    _hax.axvline(np.mean(_fd), color='black', ls='--', lw=1.5, alpha=0.8)
    _hax.set_ylim(top=_hax.get_ylim()[1] * 1.5)   # headroom so the legend clears the bars
    if _ixlim:
        _hax.set_xlim(*_ixlim)
    _hax.set_xlabel(_dlab, fontsize=A_INSET_LBL, fontweight='bold', labelpad=1.5)
    _hax.tick_params(labelsize=A_INSET_TICK)
    _hax.locator_params(axis='x', nbins=3); _hax.locator_params(axis='y', nbins=2)
    _hax.legend(fontsize=A_INSET_LEG, loc='upper right', framealpha=0.9, handlelength=1.0,
                handletextpad=0.4, borderpad=0.3, labelspacing=0.2)
figA.text(0.5, 0.975, 'Difference (Trait − PFT)', ha='center', va='top',
          fontsize=A_TITLE + 3, fontweight='bold')
_outA = _FIG_DIR / 'figure6a_difference_maps.png'
_outA_pdf = _FIG_DIR / 'figure6a_difference_maps.pdf'
figA.savefig(_outA, dpi=600, facecolor='white')
figA.savefig(_outA_pdf, format='pdf', dpi=600, facecolor='white')
plt.close(figA)

# ============ Figure 6b: time series (panels c, d) ============
print("Building Figure 6b (time series) …")
figB, (axB1, axB2) = plt.subplots(1, 2, figsize=(A_W, A_H), dpi=150)  # same dimensions as Fig 6a

# Panel (c): GPP time series ± CI95
axB1.plot(_DATES, _gpp_tr_ts, label='Trait', marker='o', color='blue', lw=2.4, ms=8)
axB1.fill_between(_DATES, np.array(_gpp_tr_ts) - _gpp_tr_ci, np.array(_gpp_tr_ts) + _gpp_tr_ci, color='blue', alpha=0.18)
axB1.plot(_DATES, _gpp_pft_ts, label='PFT', marker='s', color='orange', lw=2.4, ms=8)
axB1.fill_between(_DATES, np.array(_gpp_pft_ts) - _gpp_pft_ci, np.array(_gpp_pft_ts) + _gpp_pft_ci, color='orange', alpha=0.18)
axB1.text(-0.13, 1.03, '(a)', transform=axB1.transAxes, fontsize=TS_PANEL, fontweight='bold', va='bottom', ha='right')
axB1.set_xlabel('2022', fontsize=TS_TICK, fontweight='normal')  # year in the month-tick style
axB1.set_ylabel(r'GPP ($\mu$mol CO$_2$ m$^{-2}$ s$^{-1}$)', fontsize=TS_AXLBL, fontweight='bold')
axB1.tick_params(axis='both', labelsize=TS_TICK)
axB1.legend(loc='best', fontsize=TS_LEG, framealpha=0.9)
axB1.grid(axis='y', ls='--', alpha=0.55)
_gpp_max = max(np.nanmax(np.array(_gpp_tr_ts) + _gpp_tr_ci), np.nanmax(np.array(_gpp_pft_ts) + _gpp_pft_ci))
_gpp_min = min(np.nanmin(np.array(_gpp_tr_ts) - _gpp_tr_ci), np.nanmin(np.array(_gpp_pft_ts) - _gpp_pft_ci))
_gpp_span = _gpp_max - _gpp_min
axB1.set_ylim(_gpp_min - _gpp_span * 0.15, _gpp_max + _gpp_span * 0.15)
axB1.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
axB1.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
plt.setp(axB1.xaxis.get_majorticklabels(), rotation=35, ha='right')

# Panel (d): relative SIF — twin y-axis (model left, TROPOMI right)
_ax4_mod = axB2
_ax4_trop = axB2.twinx()
_ok_align = np.isfinite(_arr_tr) & np.isfinite(_trop_14d_mean)
_mod_mean_val = float(np.nanmean(_arr_tr[_ok_align])); _mod_std_val = float(np.nanstd(_arr_tr[_ok_align]))
_trop_mean_val = float(np.nanmean(_trop_14d_mean[_ok_align])); _trop_std_val = float(np.nanstd(_trop_14d_mean[_ok_align]))
_scale = _trop_std_val / _mod_std_val if _mod_std_val > 0 else 1.0
_offset = _trop_mean_val - _scale * _mod_mean_val
_left_lo, _left_hi = 0.625, 1.05
_right_lo = _scale * _left_lo + _offset
_right_hi = _scale * _left_hi + _offset

_ax4_mod.plot(_DATES, _sr_tr_ts, label=f'Trait (R²={_r2_tr:.2f})', marker='o', color='blue', lw=2.4, ms=8, zorder=4)
_ax4_mod.plot(_DATES, _sr_pft_ts, label=f'PFT (R²={_r2_pft:.2f})', marker='s', color='orange', lw=2.4, ms=8, zorder=4)
_ax4_mod.set_ylabel(r'Model Relative SIF$_{740}$ (%)', fontsize=TS_AXLBL, fontweight='bold', color='black')
_ax4_mod.tick_params(axis='y', labelcolor='black', labelsize=TS_TICK)
_ax4_mod.set_ylim(_left_lo, _left_hi)
_ax4_mod.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
_ax4_mod.yaxis.set_major_locator(plt.MaxNLocator(nbins=5, prune='both'))

# TROPOMI ±4d smoothed line + Gaussian-gradient SEM band
_ok4 = np.isfinite(_trop_4d_vals)
_trop_4d_smooth = np.full_like(_trop_4d_vals, np.nan)
if _ok4.sum() > 3:
    _trop_4d_smooth[_ok4] = _gf1d(_trop_4d_vals[_ok4], sigma=2.0)
_trop_4d_lo = np.full_like(_trop_4d_vals, np.nan)
_trop_4d_hi = np.full_like(_trop_4d_vals, np.nan)
for _ii, _tdate in enumerate(_DATES):
    _msk = ((pd.DatetimeIndex(_ttime) >= _tdate - pd.Timedelta(days=_WIN)) &
            (pd.DatetimeIndex(_ttime) <= _tdate + pd.Timedelta(days=_WIN)))
    _dv = np.nanmean(_sifrel_raw[_msk][:, _ri, _ci], axis=1); _dv = _dv[np.isfinite(_dv)]
    if len(_dv) >= 2:
        _sem = float(np.std(_dv, ddof=1) / np.sqrt(len(_dv)))
        _trop_4d_lo[_ii] = _trop_4d_smooth[_ii] - _sem
        _trop_4d_hi[_ii] = _trop_4d_smooth[_ii] + _sem

_ax4_trop.scatter(_DATES, _trop_4d_vals, marker='^', color='dimgray', s=70, alpha=0.45,
                  zorder=3, label='TROPOMI ±4d', clip_on=True)
# clean single ±1 SEM band around the smoothed TROPOMI line (no gradient cloud)
_ax4_trop.fill_between(_DATES, _trop_4d_lo, _trop_4d_hi, color='0.6',
                       alpha=0.15, zorder=2, linewidth=0)
_ax4_trop.plot(_DATES, _trop_4d_smooth, color='black', lw=2.4, ls='-', marker='D', ms=7,
               zorder=5, label='TROPOMI\nsmoothed ±SEM')
_ax4_trop.set_ylabel(r'TROPOMI Relative SIF$_{740}$ (%)', fontsize=TS_AXLBL, fontweight='bold',
                     color='black', rotation=270, labelpad=24)
_ax4_trop.tick_params(axis='y', labelcolor='black', labelsize=TS_TICK)
_ax4_trop.set_ylim(_right_lo, _right_hi)

for _spine in _ax4_mod.spines.values():
    _spine.set_edgecolor('black'); _spine.set_linewidth(1.2)
for _spine in _ax4_trop.spines.values():
    _spine.set_edgecolor('black'); _spine.set_linewidth(1.2)
axB2.text(-0.13, 1.03, '(b)', transform=axB2.transAxes, fontsize=TS_PANEL, fontweight='bold', va='bottom', ha='right')
_ax4_mod.set_xlabel('2022', fontsize=TS_TICK, fontweight='normal')  # year in the month-tick style
_ax4_mod.tick_params(axis='x', colors='black', labelsize=TS_TICK)
_ax4_mod.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
_ax4_mod.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
plt.setp(_ax4_mod.xaxis.get_majorticklabels(), rotation=35, ha='right')
_ax4_mod.set_xlim(_DATES[0] - pd.Timedelta(days=5), _DATES[-1] + pd.Timedelta(days=5))
_ax4_mod.grid(axis='y', ls='--', alpha=0.55)

_lines_mod, _labs_mod = _ax4_mod.get_legend_handles_labels()
_lines_trop, _labs_trop = _ax4_trop.get_legend_handles_labels()
_ax4_mod.legend(_lines_mod + _lines_trop, _labs_mod + _labs_trop, loc='lower left',
                fontsize=12, framealpha=0.95, edgecolor='gray', borderpad=0.4,
                labelspacing=0.3, handlelength=1.4, handletextpad=0.5)

figB.tight_layout()
figB.subplots_adjust(wspace=0.46)
_outB = _FIG_DIR / 'figure6b_time_series.png'
_outB_pdf = _FIG_DIR / 'figure6b_time_series.pdf'
figB.savefig(_outB, dpi=600, facecolor='white')
figB.savefig(_outB_pdf, format='pdf', dpi=600, facecolor='white')
plt.close(figB)

print(f"Saved: {_outA}")
print(f"Saved: {_outB}")
print(f"GPP diff mean={np.nanmean(_gpp_diff):.3f} std={np.nanstd(_gpp_diff):.3f}")
print(f"SIF diff mean={np.nanmean(_sif_diff):.4f} std={np.nanstd(_sif_diff):.4f}")
print(f"Relative SIF R2:  Trait={_r2_tr:.3f}  PFT={_r2_pft:.3f}")
