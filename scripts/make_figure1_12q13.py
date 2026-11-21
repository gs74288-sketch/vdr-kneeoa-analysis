# -*- coding: utf-8 -*-
"""Figure 1 v2 - rebuilt to the review checklist:
   A mirrored dual track (OA above, eQTL below, broken axis), inset, sentinel diamonds
   B full range, IVW line + CI band, one colour system
   C 100% stacked H0-H4 per window
   No in-panel titles; drawn at 170 mm print width.
"""
import sys, io, os, gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7,
                     'axes.linewidth': 0.7, 'xtick.major.width': 0.7,
                     'ytick.major.width': 0.7, 'xtick.major.size': 2.2,
                     'ytick.major.size': 2.2, 'axes.labelsize': 7.5})

DATA = r"C:\Users\liuha\OpenClawVDR\data"
TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"
OUT = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\Figures"

C_VDR, C_COL, C_OTH = '#2166AC', '#B2182B', '#BDBDBD'
C_PUR = '#6A3D9A'

COLS = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'odds_ratio',
        'standard_error', 'effect_allele_frequency', 'p_value', 'variant_id', 'ci_upper',
        'ci_lower', 'NCASES', 'NCONTROLS', 'hm', 'hm_code', 'rsid']
oa = pd.read_csv(os.path.join(DATA, 'GCST90566800_chr12_46.8-48.9Mb.tsv'),
                 sep='\t', names=COLS, header=None, dtype={'chromosome': str})
for c in ('base_pair_location', 'p_value', 'NCASES'):
    oa[c] = pd.to_numeric(oa[c], errors='coerce')
oa = oa.sort_values('NCASES', ascending=False).drop_duplicates('rsid')
oa = oa[oa.p_value.notna() & (oa.p_value > 0)].copy()
oa['mlp'] = -np.log10(oa.p_value)

ci = pd.read_csv(os.path.join(TAB, 'coloc_harmonised_input_GCST90566800.csv'))
ci['mlp_eqtl'] = -np.log10(ci.p1.clip(lower=1e-300))

VDR = (47_841_537, 47_943_048)
COL = (47_972_967, 48_004_554)
EQ_POS, EQ_RS = 47_845_054, 'rs7975232'
OA_POS, OA_RS = 48_009_982, 'rs11168351'
LO, HI = 46_800_000, 48_900_000

# ---- geometry ----------------------------------------------------------------
FIG_W, FIG_H = 6.70, 5.90          # 170 mm x 150 mm at 300 dpi
fig = plt.figure(figsize=(FIG_W, FIG_H))

ax1 = fig.add_axes([0.088, 0.735, 0.880, 0.210])   # OA track
ax2 = fig.add_axes([0.088, 0.480, 0.880, 0.200], sharex=ax1)  # eQTL track
ax3 = fig.add_axes([0.088, 0.410, 0.880, 0.052], sharex=ax1)  # gene models
axB = fig.add_axes([0.088, 0.068, 0.380, 0.268])
axC = fig.add_axes([0.615, 0.068, 0.350, 0.268])

x_oa = oa.base_pair_location.values / 1e6
x_eq = ci.pos.values / 1e6

# ============================ panel A: OA track ================================
for lo_, hi_, c in ((VDR[0], VDR[1], C_VDR), (COL[0], COL[1], C_COL)):
    ax1.axvspan(lo_ / 1e6, hi_ / 1e6, color=c, alpha=0.10, lw=0)

in_vdr = (oa.base_pair_location >= VDR[0]) & (oa.base_pair_location <= VDR[1])
in_col = (oa.base_pair_location >= COL[0]) & (oa.base_pair_location <= COL[1])
oth = ~(in_vdr | in_col)
ax1.scatter(x_oa[oth], oa.mlp[oth], s=1.6, c=C_OTH, alpha=0.55, linewidths=0, zorder=2)
ax1.scatter(x_oa[in_col], oa.mlp[in_col], s=2.6, c=C_COL, alpha=0.9, linewidths=0, zorder=4)
ax1.scatter(x_oa[in_vdr], oa.mlp[in_vdr], s=2.6, c=C_VDR, alpha=0.9, linewidths=0, zorder=4)

ax1.axhline(-np.log10(5e-8), color='0.35', ls=(0, (5, 3)), lw=0.6, zorder=1)
ax1.text(48.88, -np.log10(5e-8) + 0.35, r'$P = 5\times10^{-8}$', fontsize=5.8,
         color='0.35', ha='right', va='bottom')

# sentinels on the OA track
oa_peak = oa.loc[(oa.base_pair_location - OA_POS).abs().idxmin()]
eq_snp = oa.loc[(oa.base_pair_location - EQ_POS).abs().idxmin()]
ax1.scatter([EQ_POS / 1e6], [eq_snp.mlp], s=44, marker='D', facecolor=C_PUR,
            edgecolor='white', linewidths=0.5, zorder=8)
ax1.scatter([OA_POS / 1e6], [oa_peak.mlp], s=44, marker='D', facecolor='white',
            edgecolor=C_PUR, linewidths=0.9, zorder=8)
ax1.annotate('', xy=(OA_POS / 1e6, oa_peak.mlp + 0.5), xytext=(OA_POS / 1e6, 12.6),
             arrowprops=dict(arrowstyle='-', lw=0.5, color='0.45'))
ax1.text(OA_POS / 1e6 + 0.035, 12.7, '%s\nOA peak' % OA_RS, fontsize=6.0,
         color='0.25', ha='left', va='top')
ax1.annotate('', xy=(EQ_POS / 1e6, eq_snp.mlp + 0.4), xytext=(EQ_POS / 1e6, 8.4),
             arrowprops=dict(arrowstyle='-', lw=0.5, color='0.45'))
ax1.text(EQ_POS / 1e6 - 0.04, 8.5, '%s\nno OA signal' % EQ_RS, fontsize=6.0,
         color='0.25', ha='right', va='top')

ax1.set_ylim(0, 14.0)
ax1.set_xlim(LO / 1e6, HI / 1e6)
ax1.set_ylabel(r'Knee OA  $-\log_{10}$ P', fontsize=7.5)
ax1.tick_params(labelbottom=False)
for sp in ('top', 'right'):
    ax1.spines[sp].set_visible(False)

# inset: VDR +/-100 kb zoom
axI = fig.add_axes([0.108, 0.782, 0.150, 0.125])
z = oa[(oa.base_pair_location >= EQ_POS - 100_000) & (oa.base_pair_location <= EQ_POS + 100_000)]
axI.axvspan(VDR[0] / 1e6, VDR[1] / 1e6, color=C_VDR, alpha=0.12, lw=0)
axI.scatter(z.base_pair_location / 1e6, z.mlp, s=1.8, c=C_VDR, alpha=0.85, linewidths=0)
axI.axhline(-np.log10(5e-8), color='0.35', ls=(0, (5, 3)), lw=0.5)
axI.set_ylim(0, 8.2)
axI.set_xlim((EQ_POS - 100_000) / 1e6, (EQ_POS + 100_000) / 1e6)
axI.set_yticks([0, 2, 4, 6, 8])
axI.tick_params(labelsize=5.2, length=1.6, pad=1.5)
axI.set_title('VDR \u00b1100 kb', fontsize=5.6, pad=2)
axI.text(0.5, 0.56, 'no variant\n$P<5\\times10^{-8}$', transform=axI.transAxes,
         fontsize=5.0, ha='center', va='center', color='0.3')
for sp in ('top', 'right'):
    axI.spines[sp].set_visible(False)
axI.patch.set_alpha(0.92)

# ============================ panel A: eQTL track ==============================
CAP = 14.6
y_eq = np.minimum(ci.mlp_eqtl.values, CAP)
ax2.axvspan(VDR[0] / 1e6, VDR[1] / 1e6, color=C_VDR, alpha=0.10, lw=0)
ax2.axvspan(COL[0] / 1e6, COL[1] / 1e6, color=C_COL, alpha=0.10, lw=0)
ax2.scatter(x_eq, y_eq, s=1.6, c='#4393C3', alpha=0.75, linewidths=0, zorder=3)
ax2.axhline(0, color='0.6', lw=0.5)
ax2.set_ylim(0, 15.2)
ax2.set_ylabel(r'VDR eQTL  $-\log_{10}$ P', fontsize=7.5)
ax2.tick_params(labelbottom=False)
for sp in ('top', 'right'):
    ax2.spines[sp].set_visible(False)

# break marks on the left spine
for yy in (0.955, 0.905):
    ax2.plot([-0.010, 0.010], [yy - 0.018, yy + 0.018], transform=ax2.transAxes,
             color='k', lw=0.8, clip_on=False)
ax2.scatter([EQ_POS / 1e6], [CAP], s=26, marker='^', facecolor=C_PUR,
            edgecolor='white', linewidths=0.4, zorder=8, clip_on=False)
ax2.text(EQ_POS / 1e6 - 0.03, CAP + 0.15, '55.2', fontsize=5.9, color=C_PUR,
         ha='right', va='bottom')
ax2.scatter([OA_POS / 1e6], [CAP], s=26, marker='^', facecolor=C_PUR,
            edgecolor='white', linewidths=0.4, zorder=8, clip_on=False)
ax2.text(OA_POS / 1e6 + 0.03, CAP + 0.15, '34.3', fontsize=5.9, color=C_PUR,
         ha='left', va='bottom')
ax2.text(EQ_POS / 1e6, CAP + 1.55, 'both sentinel variants are\nstrong VDR eQTLs',
         fontsize=5.4, color='0.35', ha='center', va='bottom')

# ============================ panel A: gene models =============================
ax3.set_ylim(0, 1)
ax3.set_yticks([])
for sp in ('top', 'right', 'left'):
    ax3.spines[sp].set_visible(False)
for (g, yy, c) in (('VDR', 0.68, C_VDR), ('COL2A1', 0.20, C_COL)):
    s0, e0 = VDR if g == 'VDR' else COL
    ax3.add_patch(Rectangle((s0 / 1e6, yy - 0.10), (e0 - s0) / 1e6, 0.20,
                            facecolor=c, edgecolor='none', alpha=0.85))
    ax3.annotate('', xy=(s0 / 1e6 - 0.028, yy), xytext=(s0 / 1e6 - 0.006, yy),
                 arrowprops=dict(arrowstyle='-|>', lw=0.6, color=c, mutation_scale=5))
    ax3.text((s0 + e0) / 2e6, yy + 0.17, g, fontsize=6.2, color=c, ha='center',
             va='bottom')
ax3.set_xlabel('Chromosome 12 position (Mb, GRCh38)', fontsize=7.5)
ax3.tick_params(labelsize=6.5, pad=2)

# ============================ panel B ==========================================
ci['mlp'] = ci.mlp_eqtl
b_vdr = ci[(ci.pos >= VDR[0]) & (ci.pos <= VDR[1])]
b_col = ci[(ci.pos >= COL[0]) & (ci.pos <= COL[1])]
b_oth = ci[~ci.index.isin(b_vdr.index) & ~ci.index.isin(b_col.index)]
axB.scatter(b_oth.b1, b_oth.b2h, s=1.2, c=C_OTH, alpha=0.25, linewidths=0, zorder=2)
axB.scatter(b_vdr.b1, b_vdr.b2h, s=2.4, c=C_VDR, alpha=0.75, linewidths=0, zorder=3)
axB.scatter(b_col.b1, b_col.b2h, s=2.4, c=C_COL, alpha=0.75, linewidths=0, zorder=3)

theta = np.log(0.954)
se_theta = (np.log(1.061) - np.log(0.858)) / (2 * 1.959964)
xs = np.array([-0.36, 0.36])
axB.fill_between(xs, (theta - 1.96 * se_theta) * xs, (theta + 1.96 * se_theta) * xs,
                 color='0.25', alpha=0.13, lw=0, zorder=1)
axB.plot(xs, theta * xs, color='0.15', lw=1.0, zorder=5)
axB.text(0.035, 0.955, 'cis-IVW  0.954 (0.858\u20131.061)', transform=axB.transAxes,
         fontsize=6.2, ha='left', va='top')
axB.plot([0.035, 0.075], [0.895, 0.895], transform=axB.transAxes, color='0.15', lw=1.0)
axB.text(0.082, 0.895, 'IVW fit (3 instruments)', transform=axB.transAxes,
         fontsize=5.8, ha='left', va='center', color='0.3')

for snp, filled in ((EQ_RS, True), (OA_RS, False)):
    r = ci[ci.snp == snp]
    if r.empty:
        continue
    bx, by = float(r.b1.iloc[0]), float(r.b2h.iloc[0])
    axB.scatter([bx], [by], s=40, marker='D',
                facecolor=C_PUR if filled else 'white',
                edgecolor='white' if filled else C_PUR,
                linewidths=0.5 if filled else 0.9, zorder=9)
    dy = 0.016 if by < 0.05 else -0.020
    axB.text(bx, by + dy, snp, fontsize=5.6, ha='center',
             va='bottom' if dy > 0 else 'top', color=C_PUR, zorder=10)

axB.axhline(0, color='0.65', ls=(0, (3, 3)), lw=0.5, zorder=1)
axB.axvline(0, color='0.65', ls=(0, (3, 3)), lw=0.5, zorder=1)
axB.set_xlim(-0.355, 0.355)
axB.set_ylim(-0.215, 0.285)
axB.set_xlabel('VDR eQTL effect (SD per allele)', fontsize=7.2)
axB.set_ylabel('Knee OA log OR', fontsize=7.2)
axB.tick_params(labelsize=6.4)
for sp in ('top', 'right'):
    axB.spines[sp].set_visible(False)

handles = [Line2D([], [], marker='o', ls='', color=C_VDR, markersize=3.4, label='within VDR'),
           Line2D([], [], marker='o', ls='', color=C_COL, markersize=3.4, label='within COL2A1'),
           Line2D([], [], marker='o', ls='', color=C_OTH, markersize=3.4, label='other variants'),
           Line2D([], [], marker='D', ls='', markerfacecolor=C_PUR, markeredgecolor='white',
                  markersize=4.0, label='%s (expression peak)' % EQ_RS),
           Line2D([], [], marker='D', ls='', markerfacecolor='white', markeredgecolor=C_PUR,
                  markersize=4.0, label='%s (OA peak)' % OA_RS)]
axB.legend(handles=handles, loc='upper right', fontsize=5.6, frameon=False,
           handletextpad=0.35, borderpad=0.15, labelspacing=0.28)

# ============================ panel C ==========================================
w = pd.read_csv(os.path.join(TAB, 'coloc_GCST90566800_windows.csv')).set_index('window')
sel = ['eQTL peak +/- 10 kb', 'eQTL peak +/- 100 kb', 'OA peak +/- 25 kb',
       'full region 46.8-48.9 Mb']
lab = ['\u00b110 kb\nprecision', '\u00b1100 kb\nLD scale', 'OA peak\n\u00b125 kb',
       'full region\nsensitivity']
S = w.loc[sel]
SEG = [('H0', '#F0F0F0'), ('H1', '#C9C9C9'), ('H2', '#A0A0A0'),
       ('PP.H3', '#5C5C5C'), ('PP.H4', C_PUR)]
show = ['H0', 'H1', 'H2', 'PP.H3', 'PP.H4']
name = {'H0': 'H0 no association', 'H1': 'H1 expression only',
        'H2': 'H2 OA only', 'PP.H3': 'H3 distinct variants',
        'PP.H4': 'H4 shared variant'}

ys = np.arange(len(sel))
for i in range(len(sel)):
    left = 0.0
    for key, c in SEG:
        v = float(S.loc[sel[i], key])
        if not np.isfinite(v) or v <= 0:
            continue
        axC.barh(ys[i], v, left=left, height=0.62, color=c, edgecolor='white',
                 linewidth=0.35, zorder=3)
        left += v

axC.set_yticks(ys)
axC.set_yticklabels(lab, fontsize=6.2)
axC.invert_yaxis()
axC.set_xlim(0, 1.20)
axC.set_xticks([0, 0.5, 1.0])
axC.set_xlabel('Posterior probability', fontsize=7.2)
axC.set_ylabel('Colocalisation window', fontsize=7.2)
axC.tick_params(labelsize=6.4)
for sp in ('top', 'right'):
    axC.spines[sp].set_visible(False)

# dominant value inside the bar, H4 at the right edge
for i in range(len(sel)):
    row = S.loc[sel[i]]
    dom_key = max(show, key=lambda k: row[k])
    dom_v = float(row[dom_key])
    left = 0.0
    for key, _ in SEG:
        if key == dom_key:
            break
        left += float(row[key])
    if dom_v > 0.10:
        axC.text(left + dom_v / 2, ys[i], ('%.3f' % dom_v) if dom_v < 0.9995 else '1.000',
                 fontsize=5.6, ha='center', va='center',
                 color='white' if dom_key in ('PP.H3',) else '0.2', zorder=5)
    h4 = float(row['PP.H4'])
    txt = '%.1e' % h4 if h4 < 0.001 else '%.3f' % h4
    axC.text(1.015, ys[i], txt, fontsize=5.6, ha='left', va='center', color=C_PUR)

axC.text(1.015, -0.86, 'H4', fontsize=5.8, ha='left', va='center', color=C_PUR,
         fontweight='bold', clip_on=False)

# ============================ panel letters ===================================
for letter, xx, yy in (('A', 0.088, 0.958), ('B', 0.088, 0.345), ('C', 0.615, 0.345)):
    fig.text(xx, yy, letter, fontsize=11, fontweight='bold', ha='left', va='bottom')

os.makedirs(OUT, exist_ok=True)
png = os.path.join(OUT, 'Figure_1_v2.png')
tif = os.path.join(OUT, 'Figure_1_v2.tif')
fig.savefig(png, dpi=300)
fig.savefig(tif, dpi=300)
plt.close(fig)

from PIL import Image
im = Image.open(tif)
if im.mode != 'RGB':
    im = im.convert('RGB')
im.save(tif, dpi=(300, 300), compression='tiff_lzw')
im2 = Image.open(tif)
print('png :', png)
print('tif :', tif, im2.size, im2.mode, round(os.path.getsize(tif) / 1e6, 2), 'MB')
