"""
Interactive layout helper for Figure 5b.

Pre-renders the real difference maps (remapcon data) at a grid of white-padding
values and writes a single self-contained HTML page where you:
  * move the PAD_R / PAD_B sliders to change the white area around the preserve, and
  * drag / resize the Δ-histogram inset box (shown over map a).
The page live-reports PAD_R, PAD_B and the inset position/size as fractions of the
map frame, plus the exact lines to paste into figure_5_spatial_trait_maps_legible.py.

Run:  <venv>/bin/python plotting/fig5b_layout_tool.py
Open: serve plotting/ and open fig5b_layout_tool.html
"""
import io
import base64
import json
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')
RC = BASE / 'traits' / 'datasets' / 'clima_fit_prescribed_lai_ci_remapcon'
PFT = BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'
OUT = BASE / 'shift_dangermond_trait' / 'plotting' / 'fig5b_layout_tool.html'

CONV = {'chl': 1.0, 'lma': 1e4, 'lwc': 0.0018}
cfg = {
    'chl': dict(cmap='RdBu_r', vmin=-30, vmax=30, mult=1,    title='Chlorophyll Content',
                clab='Δ Chlorophyll\nContent (µg/cm²)', drange=50,  dlab='Δ CHL'),
    'lma': dict(cmap='RdBu_r', vmin=-80, vmax=80, mult=1,    title='Leaf Mass per Area',
                clab='Δ Leaf Mass per\nArea (g/m²)',    drange=200, dlab='Δ LMA'),
    'lwc': dict(cmap='RdBu_r', vmin=-6,  vmax=6,  mult=1000, title='Leaf Water Content',
                clab='Δ Leaf Water\nContent (g/cm²)',   drange=20,  dlab='Δ LWC'),
}
order = ['chl', 'lma', 'lwc']
letters = ['a', 'b', 'c']

# ---- data (remapcon, identical to the figure) ----
trait = {}
for t in order:
    st = np.stack([xr.open_dataset(RC / f'shift_traits_time_{i:02d}_remapcon.nc')[t]
                   .transpose('lat', 'lon').values for i in range(13)])
    trait[t] = np.nanmean(st, 0) * CONV[t]
pftm = {t: xr.open_dataset(RC / f'mean_masked_{t}_aviris_dangermond_clima_fit_remapcon.nc')[t]
        .squeeze().transpose('lat', 'lon').values * CONV[t] for t in order}
g = xr.open_dataset(RC / 'shift_traits_time_00_remapcon.nc')
lats, lons = g['lat'].values, g['lon'].values
pft = xr.open_dataset(PFT)['Band1'].transpose('lat', 'lon').values
mask = np.isin(pft, [2, 3, 4])
LON, LAT = np.meshgrid(lons, lats)
diff = {t: np.where(mask, trait[t] - pftm[t], np.nan) for t in order}

lon_r, lat_r = lons.max() - lons.min(), lats.max() - lats.min()
LATT, LONT = [34.45, 34.55], [-120.50, -120.40]
TICK = CBAR_TICK = 22; CBAR_LBL = 18; STAT = 16; PANEL = 24; TITLE = 23
Wb = 18.0
Lb, GCB, CWb, CLAB, GUNIT, RB = 0.055, 0.006, 0.013, 0.052, 0.035, 0.006
MWb = (1 - Lb - 3 * (GCB + CWb + CLAB) - 2 * GUNIT - RB) / 3
TOP_IN, BOT_IN = 0.80, 0.55
DPI = 85
Wpx = int(round(Wb * DPI))


def render_bg(pad_r, pad_b):
    xlim = (lons.min(), lons.max() + pad_r * lon_r)
    ylim = (lats.min() - pad_b * lat_r, lats.max())
    asp = (xlim[1] - xlim[0]) / (ylim[1] - ylim[0])
    map_h_in = MWb * Wb / asp
    Hb = TOP_IN + map_h_in + BOT_IN
    MHb = map_h_in / Hb
    map_y = BOT_IN / Hb
    fig = plt.figure(figsize=(Wb, Hb), dpi=DPI)
    for c, t in enumerate(order):
        x = Lb + c * (MWb + GCB + CWb + CLAB + GUNIT)
        disp = diff[t] * cfg[t]['mult']
        ax = fig.add_axes([x, map_y, MWb, MHb])
        im = ax.pcolormesh(LON, LAT, disp, cmap=cfg[t]['cmap'], vmin=cfg[t]['vmin'],
                           vmax=cfg[t]['vmax'], shading='auto', rasterized=True)
        ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xticks(LONT); ax.set_yticks(LATT)
        ax.tick_params(labelsize=TICK, labelleft=(c == 0), labelbottom=True)
        ax.annotate(f'({letters[c]})', xy=(0, 1), xycoords='axes fraction', xytext=(-2, 4),
                    textcoords='offset points', fontsize=PANEL, fontweight='bold',
                    va='bottom', ha='right')
        ax.text(0.965, 0.96, f'Mean: {np.nanmean(disp):.2f}\nSTD: {np.nanstd(disp):.2f}',
                transform=ax.transAxes, fontsize=STAT, va='top', ha='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'))
        ax.set_title(cfg[t]['title'], fontsize=TITLE, fontweight='bold', pad=8)
        cax = fig.add_axes([x + MWb + GCB, map_y, CWb, MHb])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label(cfg[t]['clab'], fontsize=CBAR_LBL, fontweight='bold')
        cb.ax.tick_params(labelsize=CBAR_TICK)
        if t == 'lwc':
            cb.ax.set_title('(×10⁻³)', fontsize=16, pad=6)
    fig.text(0.5, 0.975, 'Difference (Trait − PFT)', ha='center', va='top',
             fontsize=TITLE + 3, fontweight='bold')
    buf = io.BytesIO(); fig.savefig(buf, format='png', dpi=DPI, facecolor='white'); plt.close(fig)
    Hpx = int(round(Hb * DPI))
    ra = [Lb, map_y, MWb, MHb]                       # map a, figure fraction
    mapA = [ra[0] * Wpx, (1 - ra[1] - ra[3]) * Hpx, ra[2] * Wpx, ra[3] * Hpx]
    return dict(bg=base64.b64encode(buf.getvalue()).decode(), Hpx=Hpx,
                ax=round(mapA[0], 1), ay=round(mapA[1], 1), aw=round(mapA[2], 1), ah=round(mapA[3], 1))


PADR = [0.15, 0.25, 0.35, 0.45]
PADB = [0.15, 0.25, 0.35, 0.45]
COMBO = {}
for ri, pr in enumerate(PADR):
    for bi, pb in enumerate(PADB):
        COMBO[f'{ri}_{bi}'] = render_bg(pr, pb)
        print(f'rendered PAD_R={pr} PAD_B={pb}')

# histogram preview (map a)
t0 = order[0]
hf = plt.figure(figsize=(2.6, 2.1), dpi=110)
hax = hf.add_axes([0.26, 0.30, 0.70, 0.64])
v = diff[t0] * cfg[t0]['mult']; v = v[np.isfinite(v) & (v != 0)]
hax.hist(v, bins=45, color='0.55', edgecolor='black', linewidth=0.4, density=True)
hax.axvline(0, color='red', ls='--', lw=1.5); hax.axvline(np.nanmean(v), color='blue', ls='--', lw=1.5)
hax.set_xlim(-cfg[t0]['drange'], cfg[t0]['drange']); hax.set_xticks([-cfg[t0]['drange'], 0, cfg[t0]['drange']])
hax.set_xlabel(cfg[t0]['dlab'], fontsize=13, fontweight='bold', labelpad=1.5)
hax.tick_params(labelsize=11); hax.locator_params(axis='y', nbins=2)
hbuf = io.BytesIO(); hf.savefig(hbuf, format='png', dpi=110, transparent=True); plt.close(hf)
hist_b64 = base64.b64encode(hbuf.getvalue()).decode()

html = """<!doctype html><html><head><meta charset="utf-8"><title>Fig 5b layout</title>
<style>
 body{font-family:sans-serif;margin:16px;background:#fafafa}
 #ctl{margin:8px 0 14px} label{font-weight:bold;margin-right:6px}
 #wrap{position:relative;display:inline-block;border:1px solid #ccc}
 #bg{display:block;width:1200px;height:auto}
 #box{position:absolute;border:2px solid #d00;box-sizing:border-box;cursor:move;background:transparent}
 #box img{width:100%;height:100%;object-fit:fill;pointer-events:none}
 #h{position:absolute;right:-7px;bottom:-7px;width:14px;height:14px;background:#d00;cursor:se-resize;border-radius:2px}
 #out{margin-top:12px;font-family:monospace;font-size:14px;background:#111;color:#0f0;padding:12px;border-radius:6px;white-space:pre;display:inline-block}
 .hint{color:#555;font-size:13px;margin:4px 0}
</style></head><body>
<h3>Figure 5b — padding sliders + draggable inset</h3>
<div class="hint">Move sliders to change the white area; drag/resize the red inset over map (a). Paste the block below to me.</div>
<div id="ctl">
 <label>PAD_R (right white)</label><input id="sr" type="range" min="0" max="__NR__" step="1" value="__DR__"><span id="lr"></span>
 &nbsp;&nbsp;<label>PAD_B (bottom white)</label><input id="sb" type="range" min="0" max="__NB__" step="1" value="__DB__"><span id="lb"></span>
</div>
<div id="wrap"><img id="bg"><div id="box"><img src="data:image/png;base64,__HIST__"><div id="h"></div></div></div>
<div id="out"></div>
<script>
const DPI=__DPI__, Wpx=__WPX__, DISPW=1200;
const PADR=__PADR__, PADB=__PADB__, COMBO=__COMBO__;
let ri=__DR__, bi=__DB__;
let frac={ix:0.55, iy:0.12, iw:0.34, ih:0.30};
const bg=document.getElementById('bg'), box=document.getElementById('box'), out=document.getElementById('out');
const sr=document.getElementById('sr'), sb=document.getElementById('sb'), lr=document.getElementById('lr'), lb=document.getElementById('lb');
function C(){return COMBO[ri+'_'+bi];}
function S(){return DISPW/Wpx;}
function placeBox(){ const c=C(), s=S();
  box.style.left=((c.ax+frac.ix*c.aw)*s)+'px';
  box.style.top=((c.ay+(1-frac.iy-frac.ih)*c.ah)*s)+'px';
  box.style.width=(frac.iw*c.aw*s)+'px';
  box.style.height=(frac.ih*c.ah*s)+'px'; }
function swap(){ const c=C(); bg.src='data:image/png;base64,'+c.bg; bg.style.width=DISPW+'px';
  lr.textContent=PADR[ri].toFixed(2); lb.textContent=PADB[bi].toFixed(2); placeBox(); report(); }
function fracFromBox(){ const c=C(), s=S();
  const bx=parseFloat(box.style.left)/s, by=parseFloat(box.style.top)/s;
  const bw=parseFloat(box.style.width)/s, bh=parseFloat(box.style.height)/s;
  frac.ix=(bx-c.ax)/c.aw; frac.iw=bw/c.aw; frac.ih=bh/c.ah; frac.iy=((c.ay+c.ah)-(by+bh))/c.ah; }
function report(){ out.textContent =
  `PAD_R, PAD_B = ${PADR[ri].toFixed(2)}, ${PADB[bi].toFixed(2)}\\n`+
  `inset: ix=${frac.ix.toFixed(3)} iy=${frac.iy.toFixed(3)} iw=${frac.iw.toFixed(3)} ih=${frac.ih.toFixed(3)}\\n\\n`+
  `PAD_R, PAD_B = ${PADR[ri].toFixed(2)}, ${PADB[bi].toFixed(2)}\\n`+
  `iw, ih = ${frac.iw.toFixed(3)} * MWb, ${frac.ih.toFixed(3)} * MHb\\n`+
  `hax = figB.add_axes([x + ${frac.ix.toFixed(3)} * MWb, map_y + ${frac.iy.toFixed(3)} * MHb, iw, ih])`; }
sr.addEventListener('input',()=>{ri=+sr.value; swap();});
sb.addEventListener('input',()=>{bi=+sb.value; swap();});
let drag=null;
box.addEventListener('mousedown',e=>{ if(e.target.id==='h')return; drag={m:'move',x:e.clientX,y:e.clientY,l:parseFloat(box.style.left),t:parseFloat(box.style.top)}; e.preventDefault();});
document.getElementById('h').addEventListener('mousedown',e=>{ const w=parseFloat(box.style.width),h=parseFloat(box.style.height); drag={m:'size',x:e.clientX,y:e.clientY,w:w,h:h,r:w/h}; e.stopPropagation(); e.preventDefault();});
window.addEventListener('mousemove',e=>{ if(!drag)return;
  if(drag.m==='move'){ box.style.left=(drag.l+e.clientX-drag.x)+'px'; box.style.top=(drag.t+e.clientY-drag.y)+'px'; }
  else { const nw=Math.max(40,drag.w+e.clientX-drag.x); box.style.width=nw+'px'; box.style.height=(nw/drag.r)+'px'; }  /* locked aspect ratio */
  fracFromBox(); report(); });
window.addEventListener('mouseup',()=>drag=null);
bg.addEventListener('load',()=>{placeBox();report();});
swap();
</script></body></html>"""

html = (html
        .replace('__DPI__', str(DPI)).replace('__WPX__', str(Wpx))
        .replace('__PADR__', json.dumps(PADR)).replace('__PADB__', json.dumps(PADB))
        .replace('__COMBO__', json.dumps(COMBO))
        .replace('__NR__', str(len(PADR) - 1)).replace('__NB__', str(len(PADB) - 1))
        .replace('__DR__', '2').replace('__DB__', '1')     # default PAD_R=0.35, PAD_B=0.25
        .replace('__HIST__', hist_b64))
OUT.write_text(html)
print(f'wrote {OUT}  ({OUT.stat().st_size // 1024} KB)')
