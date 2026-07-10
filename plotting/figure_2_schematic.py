"""
Assemble the complete Figure 2 workflow schematic from the six regenerated
panels (figure_2_schematic_panels.py). Vector schematic + embedded raster maps,
exported as PDF (vector) and 600-dpi PNG.

Layout:
  AVIRIS-NG spectrum + LiDAR clumping index --(1)--> CliMA Land
  Traits (CHL/LMA/LWC) --(2)--> Inverse modeling <-> CliMA Land
  ERA5 reanalysis --(3)--> CliMA Land --> Forward modeling
  Forward modeling --(4)--> Fluxes (GPP/SIF)
"""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Ellipse, FancyBboxPatch, FancyArrowPatch
from pathlib import Path

plt.rcParams.update({'font.family': 'Nimbus Sans'})

FIG = Path('/home/renatob/data/FluoData1/aviris_dangermond/shift_dangermond_trait/figures')
TEAL, LABEL, LABEL_TX = '#2b6c8f', '#bdd7ee', '#14364d'
ASPECT = 12 / 7
LW = 3.2

fig = plt.figure(figsize=(12, 7), dpi=150)
ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, 1); ov.set_ylim(0, 1); ov.axis('off')

def panel(path, l, b, w, h):
    ax = fig.add_axes([l, b, w, h]); ax.imshow(mpimg.imread(path)); ax.axis('off')

def label(x, y, text, fs=16):
    ov.text(x, y, text, ha='center', va='center', fontsize=fs, fontweight='bold',
            color=LABEL_TX, zorder=6, bbox=dict(boxstyle='round,pad=0.35', fc=LABEL, ec='none'))

def node_circle(x, y, text, rx=0.078):
    ov.add_patch(Ellipse((x, y), 2*rx, 2*rx*ASPECT, fc=TEAL, ec='none', zorder=4))
    ov.text(x, y, text, ha='center', va='center', color='white', fontsize=25, fontweight='bold', zorder=5)

def box(x, y, text, w=0.115, h=0.135):
    ov.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle='round,pad=0.006', fc=TEAL, ec='none', zorder=4))
    ov.text(x, y, text, ha='center', va='center', color='white', fontsize=16, fontweight='bold', zorder=5)

def marker(x, y, n):
    ov.add_patch(Ellipse((x, y), 0.034, 0.034*ASPECT, fc=TEAL, ec='white', lw=1.8, zorder=8))
    ov.text(x, y, str(n), ha='center', va='center', color='white', fontsize=15, fontweight='bold', zorder=9)

def arrow(p0, p1, double=False):
    ov.add_patch(FancyArrowPatch(p0, p1, arrowstyle='<|-|>' if double else '-|>',
                                 mutation_scale=26, lw=LW, color=TEAL, zorder=3))

def line(xs, ys):
    ov.plot(xs, ys, color=TEAL, lw=LW, zorder=3, solid_capstyle='round')

def bracket(x, y0, y1, tick, connect_to, cy):
    line([x, x], [y0, y1]); line([x, x+tick], [y0, y0]); line([x, x+tick], [y1, y1]); line([x, connect_to], [cy, cy])

# ---- panels (columns share top=0.86 / bottom=0.10; flux centered on middle trait) ----
TL, TW, TH = 0.030, 0.185, 0.235
panel(FIG/'fig2_trait_CHL.png', TL, 0.625, TW, TH)
panel(FIG/'fig2_trait_LMA.png', TL, 0.3625, TW, TH)
panel(FIG/'fig2_trait_LWC.png', TL, 0.100, TW, TH)
label(0.123, 0.905, 'Traits')
FL = 0.795
panel(FIG/'fig2_flux_GPP.png', FL, 0.495, TW, TH)
panel(FIG/'fig2_flux_SIF.png', FL, 0.230, TW, TH)
label(0.888, 0.905, 'Fluxes')
panel(FIG/'fig2_spectrum.png', 0.400, 0.660, 0.200, 0.200)
label(0.500, 0.955, 'AVIRIS-NG')
label(0.775, 0.790, 'LiDAR Clumping Index', fs=14)
ov.text(0.632, 0.790, '+', ha='center', va='center', fontsize=30, fontweight='bold', color=TEAL, zorder=6)

# ---- nodes (vertical spine at x=0.5) ----
node_circle(0.500, 0.470, 'CliMA\nLand')
box(0.300, 0.470, 'Inverse\nmodeling')
box(0.700, 0.470, 'Forward\nmodeling')
label(0.500, 0.070, 'ERA5 reanalysis')

# ---- connectors ----
arrow((0.352, 0.470), (0.420, 0.470), double=True)   # inverse <-> CliMA
arrow((0.580, 0.470), (0.648, 0.470))                # CliMA -> forward
# (1) AVIRIS-NG + LiDAR merge -> CliMA (top)
line([0.500, 0.500], [0.655, 0.600])                 # spectrum down
line([0.775, 0.775], [0.755, 0.630])                 # lidar down
line([0.500, 0.775], [0.630, 0.630])                 # merge bar
arrow((0.500, 0.630), (0.500, 0.600))
marker(0.500, 0.630, 1)
# (2) Traits -> Inverse
bracket(0.225, 0.100, 0.860, 0.013, 0.245, 0.470); marker(0.245, 0.560, 2)
# (3) ERA5 -> CliMA (bottom)
arrow((0.500, 0.140), (0.500, 0.338)); marker(0.500, 0.235, 3)
# (4) Forward -> Fluxes
bracket(0.780, 0.230, 0.730, -0.013, 0.755, 0.470); marker(0.755, 0.560, 4)

for ext in ('pdf', 'png'):
    fig.savefig(FIG/f'figure2_schematic.{ext}', dpi=600, facecolor='white', bbox_inches='tight', pad_inches=0.05)
plt.close(fig)
print('saved', FIG/'figure2_schematic.pdf', '+ .png (600 dpi)')
