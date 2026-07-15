"""
LIVE interactive layout helper for Figure 5b (WYSIWYG).

Unlike a static preview, this starts a tiny local web server that RE-RENDERS the
real Figure 5b (identical code, transparent Δ-distribution insets) every time you
move the box or a padding slider. What you see is exactly the figure, so the
coordinates it reports plug straight into figure_5_spatial_trait_maps_legible.py
(PAD_R/PAD_B and INSET_POS) and reproduce the same result.

Run:   <venv>/bin/python plotting/fig5b_layout_tool.py
Open:  http://localhost:8000/    (Ctrl-C in the terminal to stop)
"""
import io
import json
import base64
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PORT = 8000
DPI = 80
DISPW = 1200

BASE = Path('/home/renatob/data/FluoData1/aviris_dangermond')
RC = BASE / 'traits' / 'datasets' / 'clima_fit_prescribed_lai_ci_remapcon'
PFT = BASE / 'California_Vegetation_WHRTYPE_Dangermond' / 'output_latlon.nc'

# ---- data (remapcon; identical to the figure) ----
CONV = {'chl': 1.0, 'lma': 1e4, 'lwc': 0.0018}
trait = {}
for _t in ['chl', 'lma', 'lwc']:
    st = np.stack([xr.open_dataset(RC / f'shift_traits_time_{i:02d}_remapcon.nc')[_t]
                   .transpose('lat', 'lon').values for i in range(13)])
    trait[_t] = np.nanmean(st, 0) * CONV[_t]
pftm = {_t: xr.open_dataset(RC / f'mean_masked_{_t}_aviris_dangermond_clima_fit_remapcon.nc')[_t]
        .squeeze().transpose('lat', 'lon').values * CONV[_t] for _t in ['chl', 'lma', 'lwc']}
_g = xr.open_dataset(RC / 'shift_traits_time_00_remapcon.nc')
lats, lons = _g['lat'].values, _g['lon'].values
pft_map = xr.open_dataset(PFT)['Band1'].transpose('lat', 'lon').values
pft_mask = np.isin(pft_map, [2, 3, 4])
LON, LAT = np.meshgrid(lons, lats)

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
diff = {t: np.where(pft_mask, trait[t] - pftm[t], np.nan) for t in order}

lon_r, lat_r = lons.max() - lons.min(), lats.max() - lats.min()
LATT, LONT = [34.45, 34.55], [-120.50, -120.40]
TICK = CBAR_TICK = 22; CBAR_LBL = 18; STAT = 16; PANEL = 24; TITLE = 23
INSET_LBL, INSET_TICK = 13, 11
Wb = 18.0
Lb, GCB, CWb, CLAB, GUNIT, RB = 0.055, 0.006, 0.013, 0.052, 0.050, 0.035
MWb = (1 - Lb - 3 * (GCB + CWb + CLAB) - 2 * GUNIT - RB) / 3
TOP_IN, BOT_IN = 1.05, 0.55


def render(padr, padb, ix, iy, iw, ih):
    """Render the real Figure 5b (transparent insets) and return (png, mapA_px, W, H)."""
    xlim = (lons.min(), lons.max() + padr * lon_r)
    ylim = (lats.min() - padb * lat_r, lats.max())
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
        # transparent Δ-distribution inset at the requested coordinates
        hax = fig.add_axes([x + ix * MWb, map_y + iy * MHb, iw * MWb, ih * MHb])
        hax.patch.set_alpha(0)
        v = disp[np.isfinite(disp) & (disp != 0)]
        hax.hist(v, bins=45, color='0.55', edgecolor='black', linewidth=0.4, density=True)
        hax.axvline(0, color='red', ls='--', lw=1.5); hax.axvline(np.nanmean(v), color='blue', ls='--', lw=1.5)
        hax.set_xlim(-cfg[t]['drange'], cfg[t]['drange']); hax.set_xticks([-cfg[t]['drange'], 0, cfg[t]['drange']])
        hax.set_xlabel(cfg[t]['dlab'], fontsize=INSET_LBL, fontweight='bold', labelpad=1.5)
        hax.tick_params(labelsize=INSET_TICK); hax.locator_params(axis='y', nbins=2)
    fig.text(0.5, 0.975, 'Difference (Trait − PFT)', ha='center', va='top',
             fontsize=TITLE + 3, fontweight='bold')
    buf = io.BytesIO(); fig.savefig(buf, format='png', dpi=DPI, facecolor='white'); plt.close(fig)
    Wpx, Hpx = int(round(Wb * DPI)), int(round(Hb * DPI))
    # panel-a map frame in image pixels (top-left origin)
    mapA = dict(ax=round(Lb * Wpx, 1), ay=round((1 - map_y - MHb) * Hpx, 1),
                aw=round(MWb * Wpx, 1), ah=round(MHb * Hpx, 1))
    return buf.getvalue(), mapA, Wpx, Hpx


HTML = """<!doctype html><html><head><meta charset="utf-8"><title>Fig 5b live layout</title>
<style>
 body{font-family:sans-serif;margin:16px;background:#fafafa}
 #ctl{margin:8px 0 12px} label{font-weight:bold;margin-right:4px}
 #wrap{position:relative;display:inline-block;border:1px solid #ccc}
 #fig{display:block;width:__DISPW__px;height:auto}
 #box{position:absolute;border:2px solid #d00;box-sizing:border-box;cursor:move;background:transparent}
 #h{position:absolute;right:-7px;bottom:-7px;width:14px;height:14px;background:#d00;cursor:se-resize;border-radius:2px}
 #busy{position:absolute;top:8px;left:8px;background:#d00;color:#fff;padding:2px 8px;border-radius:4px;display:none}
 #out{margin-top:12px;font-family:monospace;font-size:14px;background:#111;color:#0f0;padding:12px;border-radius:6px;white-space:pre;display:inline-block}
 .hint{color:#555;font-size:13px}
</style></head><body>
<h3>Figure 5b — live layout (this IS the real figure)</h3>
<div class="hint">Drag / resize the red box (= the histogram plot axes). It re-renders on release. Move sliders for white padding. Paste the block below to me.</div>
<div id="ctl">
 <label>PAD_R</label><input id="sr" type="range" min="0.05" max="0.60" step="0.01" value="0.15"><span id="lr">0.15</span>
 &nbsp;&nbsp;<label>PAD_B</label><input id="sb" type="range" min="0.05" max="0.60" step="0.01" value="0.15"><span id="lb">0.15</span>
 &nbsp;&nbsp;<label><input id="lock" type="checkbox"> lock aspect</label>
</div>
<div id="wrap"><img id="fig"><div id="busy">rendering…</div><div id="box"><div id="h"></div></div></div>
<div id="out"></div>
<script>
const DISPW=__DISPW__;
let padr=0.15, padb=0.15;
let ix=0.462, iy=-0.034, iw=0.534, ih=0.471;
let mapA={ax:0,ay:0,aw:1,ah:1}, Wpx=1, Hpx=1;
const fig=document.getElementById('fig'), box=document.getElementById('box'), out=document.getElementById('out');
const busy=document.getElementById('busy');
const sr=document.getElementById('sr'), sb=document.getElementById('sb'), lr=document.getElementById('lr'), lb=document.getElementById('lb');
function S(){return DISPW/Wpx;}
function placeBox(){ const s=S();
  box.style.left =((mapA.ax+ix*mapA.aw)*s)+'px';
  box.style.width =(iw*mapA.aw*s)+'px';
  box.style.height=(ih*mapA.ah*s)+'px';
  box.style.top  =(((mapA.ay+mapA.ah)-(iy+ih)*mapA.ah)*s)+'px'; }
function coordsFromBox(){ const s=S();
  const l=parseFloat(box.style.left)/s, t=parseFloat(box.style.top)/s;
  const w=parseFloat(box.style.width)/s, h=parseFloat(box.style.height)/s;
  ix=(l-mapA.ax)/mapA.aw; iw=w/mapA.aw; ih=h/mapA.ah;
  iy=((mapA.ay+mapA.ah)-(t+h))/mapA.ah; }
function report(){ out.textContent =
  `ix=${ix.toFixed(3)} iy=${iy.toFixed(3)} iw=${iw.toFixed(3)} ih=${ih.toFixed(3)}\\n\\n`+
  `PAD_R, PAD_B = ${padr.toFixed(2)}, ${padb.toFixed(2)}\\n`+
  `INSET_POS = (${ix.toFixed(3)}, ${iy.toFixed(3)}, ${iw.toFixed(3)}, ${ih.toFixed(3)})`; }
async function rerender(){ busy.style.display='block';
  const u=`/render?padr=${padr}&padb=${padb}&ix=${ix}&iy=${iy}&iw=${iw}&ih=${ih}`;
  const r=await fetch(u); const j=await r.json();
  fig.src='data:image/png;base64,'+j.png; fig.style.width=DISPW+'px';
  mapA={ax:j.ax,ay:j.ay,aw:j.aw,ah:j.ah}; Wpx=j.Wpx; Hpx=j.Hpx;
  busy.style.display='none'; placeBox(); report(); }
let deb; function debounced(){ clearTimeout(deb); deb=setTimeout(rerender,180); }
sr.addEventListener('input',()=>{padr=+sr.value; lr.textContent=padr.toFixed(2); debounced();});
sb.addEventListener('input',()=>{padb=+sb.value; lb.textContent=padb.toFixed(2); debounced();});
let drag=null;
box.addEventListener('mousedown',e=>{ if(e.target.id==='h')return; drag={m:'move',x:e.clientX,y:e.clientY,l:parseFloat(box.style.left),t:parseFloat(box.style.top)}; e.preventDefault();});
document.getElementById('h').addEventListener('mousedown',e=>{ const w=parseFloat(box.style.width),h=parseFloat(box.style.height); drag={m:'size',x:e.clientX,y:e.clientY,w:w,h:h,r:w/h}; e.stopPropagation(); e.preventDefault();});
window.addEventListener('mousemove',e=>{ if(!drag)return;
  if(drag.m==='move'){ box.style.left=(drag.l+e.clientX-drag.x)+'px'; box.style.top=(drag.t+e.clientY-drag.y)+'px'; }
  else { const nw=Math.max(30,drag.w+e.clientX-drag.x); box.style.width=nw+'px';
         box.style.height=(document.getElementById('lock').checked?(nw/drag.r):Math.max(24,drag.h+e.clientY-drag.y))+'px'; }
  coordsFromBox(); report(); });
window.addEventListener('mouseup',()=>{ if(drag){drag=null; coordsFromBox(); rerender();} });
fig.addEventListener('load',()=>{placeBox();});
rerender();
</script></body></html>"""
HTML = HTML.replace('__DISPW__', str(DISPW))


class H(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        if u.path == '/':
            body = HTML.encode()
            self.send_response(200); self.send_header('Content-Type', 'text/html')
            self.send_header('Content-Length', str(len(body))); self.end_headers()
            self.wfile.write(body)
        elif u.path == '/render':
            q = urllib.parse.parse_qs(u.query)
            g = lambda k, d: float(q.get(k, [d])[0])
            png, mapA, Wpx, Hpx = render(g('padr', 0.15), g('padb', 0.15),
                                         g('ix', 0.462), g('iy', -0.034), g('iw', 0.534), g('ih', 0.471))
            payload = json.dumps({'png': base64.b64encode(png).decode(), 'Wpx': Wpx, 'Hpx': Hpx, **mapA}).encode()
            self.send_response(200); self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Length', str(len(payload))); self.end_headers()
            self.wfile.write(payload)
        else:
            self.send_response(404); self.end_headers()


if __name__ == '__main__':
    ThreadingHTTPServer.allow_reuse_address = True
    srv = ThreadingHTTPServer(('127.0.0.1', PORT), H)
    print(f'Fig 5b live layout tool  ->  http://localhost:{PORT}/   (Ctrl-C to stop)')
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print('\nstopped')
