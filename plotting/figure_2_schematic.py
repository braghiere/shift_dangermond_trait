"""
Assemble the complete Figure 2 workflow schematic from the six regenerated
panels (see figure_2_schematic_panels.py). Produces a single high-resolution
figure file — no PowerPoint step required.

Layout mirrors the accepted manuscript:
  AVIRIS-NG spectrum (+ LiDAR clumping index)  --(1)-->  CliMA Land
  Traits column (CHL/LMA/LWC)  --(2)-->  Inverse modeling  <-> CliMA Land
  ERA5 reanalysis  --(3)-->  CliMA Land  --> Forward modeling
  Forward modeling  --(4)-->  Fluxes column (GPP/SIF)
"""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Ellipse, FancyBboxPatch, FancyArrowPatch
from pathlib import Path

FIG = Path('/home/renatob/data/FluoData1/aviris_dangermond/shift_dangermond_trait/figures')

TEAL = '#2b6c8f'
LABEL = '#bdd7ee'
LABEL_TX = '#14364d'
ASPECT = 12 / 7  # fig width/height, for drawing true circles

fig = plt.figure(figsize=(12, 7), dpi=150)
ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, 1); ov.set_ylim(0, 1); ov.axis('off')

def panel(path, l, b, w, h):
    ax = fig.add_axes([l, b, w, h]); ax.imshow(mpimg.imread(path)); ax.axis('off')

def label(x, y, text, fs=15, weight='bold'):
    ov.text(x, y, text, ha='center', va='center', fontsize=fs, fontweight=weight,
            color=LABEL_TX, zorder=6,
            bbox=dict(boxstyle='round,pad=0.35', fc=LABEL, ec='none'))

def node_circle(x, y, text, rx=0.075):
    ov.add_patch(Ellipse((x, y), width=2*rx, height=2*rx*ASPECT, fc=TEAL, ec='none', zorder=4))
    ov.text(x, y, text, ha='center', va='center', color='white', fontsize=24,
            fontweight='bold', zorder=5)

def box(x, y, text, w=0.11, h=0.13):
    ov.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle='round,pad=0.005',
                                fc=TEAL, ec='none', zorder=4))
    ov.text(x, y, text, ha='center', va='center', color='white', fontsize=15,
            fontweight='bold', zorder=5)

def marker(x, y, n):
    ov.add_patch(Ellipse((x, y), 0.032, 0.032*ASPECT, fc=TEAL, ec='white', lw=1.5, zorder=8))
    ov.text(x, y, str(n), ha='center', va='center', color='white', fontsize=14,
            fontweight='bold', zorder=9)

def arrow(p0, p1, double=False):
    style = '<|-|>' if double else '-|>'
    ov.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=22,
                                 lw=2.2, color=TEAL, zorder=3))

def bracket(x, y0, y1, tick, connect_to, cy):
    """Vertical bracket at x spanning [y0,y1]; ticks of length `tick` (signed);
    connector from mid to point connect_to at height cy."""
    ov.plot([x, x], [y0, y1], color=TEAL, lw=2.2, zorder=3)
    ov.plot([x, x + tick], [y0, y0], color=TEAL, lw=2.2, zorder=3)
    ov.plot([x, x + tick], [y1, y1], color=TEAL, lw=2.2, zorder=3)
    ov.plot([x, connect_to], [cy, cy], color=TEAL, lw=2.2, zorder=3)

# ---------------- panels ----------------
# Traits column (left)
panel(FIG / 'fig2_trait_CHL.png', 0.035, 0.605, 0.175, 0.235)
panel(FIG / 'fig2_trait_LMA.png', 0.035, 0.360, 0.175, 0.235)
panel(FIG / 'fig2_trait_LWC.png', 0.035, 0.115, 0.175, 0.235)
label(0.122, 0.885, 'Traits')
# Fluxes column (right)
panel(FIG / 'fig2_flux_GPP.png', 0.790, 0.470, 0.175, 0.235)
panel(FIG / 'fig2_flux_SIF.png', 0.790, 0.210, 0.175, 0.235)
label(0.878, 0.885, 'Fluxes')
# AVIRIS-NG spectrum (top center)
panel(FIG / 'fig2_spectrum.png', 0.360, 0.660, 0.230, 0.235)
label(0.475, 0.945, 'AVIRIS-NG')
label(0.700, 0.820, '+ LiDAR Clumping Index', fs=13)

# ---------------- nodes ----------------
node_circle(0.500, 0.500, 'CliMA\nLand')
box(0.300, 0.500, 'Inverse\nmodeling')
box(0.700, 0.500, 'Forward\nmodeling')
label(0.500, 0.075, 'ERA5 reanalysis')

# ---------------- connectors ----------------
arrow((0.352, 0.500), (0.423, 0.500), double=True)     # inverse <-> CliMA
arrow((0.577, 0.500), (0.648, 0.500))                  # CliMA -> forward
# (1) inputs -> CliMA (top)
ov.plot([0.475, 0.475], [0.655, 0.615], color=TEAL, lw=2.2, zorder=3)
ov.plot([0.700, 0.700], [0.795, 0.700], color=TEAL, lw=2.2, zorder=3)
ov.plot([0.475, 0.700], [0.615, 0.615], color=TEAL, lw=2.2, zorder=3)
arrow((0.5875, 0.615), (0.5875, 0.575))
marker(0.5875, 0.615, 1)
# (2) Traits -> Inverse modeling
bracket(0.220, 0.130, 0.825, 0.012, 0.245, 0.478)
marker(0.240, 0.560, 2)
# (3) ERA5 -> CliMA (bottom)
arrow((0.500, 0.145), (0.500, 0.358))
marker(0.500, 0.250, 3)
# (4) Forward modeling -> Fluxes
bracket(0.780, 0.220, 0.690, -0.012, 0.755, 0.455)
marker(0.760, 0.520, 4)

for ext in ('png', 'pdf'):
    fig.savefig(FIG / f'figure2_schematic.{ext}', dpi=300, facecolor='white', bbox_inches='tight')
plt.close(fig)
print('saved', FIG / 'figure2_schematic.png')
