# -*- coding: utf-8 -*-
"""Figure: the 12q13 knee OA signal maps to COL2A1, not VDR."""
import sys, io, gzip, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                     'axes.linewidth': 0.8, 'xtick.major.width': 0.8,
                     'ytick.major.width': 0.8})

DATA = r"C:\Users\liuha\OpenClawVDR\data"
TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"
OUT = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\Figures"

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

VDR = (47_841_537, 47_943_048)
COL = (47_972_967, 48_004_554)
EQTL_PEAK, EQTL_RS = 47_845_054, 'rs7975232'
OA_PEAK, OA_RS = 48_009_982, 'rs11168351'
LO, HI = 46_800_000, 48_900_000

col = {
    'VDR': '#2166AC',
    'COL2A1': '#B2182B',
    'other': '#BDBDBD',
}


def region(pos):
    if VDR[0] <= pos <= VDR[1]:
        return 'VDR'
    if COL[0] <= pos <= COL[1]:
        return 'COL2A1'
    return 'other'


oa['reg'] = [region(p) for p in oa.base_pair_location]

# ---- exon models from the GTF ----
GTF = os.path.join(DATA, 'Homo_sapiens.GRCh38.112.gtf.gz')
want = {'VDR', 'COL2A1'}
exons = {g: [] for g in want}
gene_span = {}
with gzip.open(GTF, 'rt', encoding='utf-8', errors='replace') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != '12':
            continue
        s, e = int(f[3]), int(f[4])
        if e < LO or s > HI:
            continue
        name = None
        for part in f[8].split(';'):
            part = part.strip()
            if part.startswith('gene_name '):
                name = part.split(' ', 1)[1].strip('"')
                break
        if name not in want:
            continue
        if f[2] == 'exon':
            exons[name].append((s, e))
        elif f[2] == 'gene':
            gene_span[name] = (s, e, f[6])

for g in exons:
    exons[g] = sorted(exons[g])
print('exons:', {g: len(v) for g, v in exons.items()})
print('spans:', gene_span)

# ---- windowed posteriors ----
w = pd.read_csv(os.path.join(TAB, 'coloc_GCST90566800_windows.csv'))
print(w[['window', 'H1', 'PP.H3', 'PP.H4']].to_string(index=False))

# ---- figure ----
fig = plt.figure(figsize=(7.5, 8.3))
gs = GridSpec(3, 2, height_ratios=[2.0, 0.50, 1.45], width_ratios=[1.15, 1.0],
              hspace=0.10, wspace=0.30,
              left=0.105, right=0.975, top=0.965, bottom=0.115)

axA = fig.add_subplot(gs[0, :])
axG = fig.add_subplot(gs[1, :], sharex=axA)
axB = fig.add_subplot(gs[2, 0])
axC = fig.add_subplot(gs[2, 1])

x = oa.base_pair_location.values / 1e6

for ax in (axA,):
    ax.axvspan(VDR[0] / 1e6, VDR[1] / 1e6, color=col['VDR'], alpha=0.10, lw=0)
    ax.axvspan(COL[0] / 1e6, COL[1] / 1e6, color=col['COL2A1'], alpha=0.10, lw=0)

for key, z in (('other', 2), ('VDR', 3), ('COL2A1', 4)):
    m = oa.reg == key
    axA.scatter(oa.base_pair_location[m] / 1e6, oa.mlp[m], s=4, c=col[key],
                alpha=0.75 if key == 'other' else 0.95, linewidths=0, zorder=z,
                label={'other': 'other variants',
                       'VDR': 'within VDR',
                       'COL2A1': 'within COL2A1'}[key])

axA.axhline(-np.log10(5e-8), color='k', ls='--', lw=0.8, zorder=1)
axA.text(46.85, -np.log10(5e-8) + 0.45, 'genome-wide significance', fontsize=6.5, va='bottom')

axA.axvline(EQTL_PEAK / 1e6, color=col['VDR'], ls=(0, (4, 2)), lw=1.0, zorder=1)
axA.annotate('VDR eQTL peak\n%s  P = 6.5 × 10⁻⁵⁶\nno OA association (P = 0.58)' % EQTL_RS,
             xy=(EQTL_PEAK / 1e6, 4.2), xytext=(47.02, 9.6), fontsize=6.8, color=col['VDR'],
             ha='left', va='top',
             arrowprops=dict(arrowstyle='-', lw=0.7, color=col['VDR']))
axA.annotate('OA peak\n%s  P = 1.5 × 10⁻¹²' % OA_RS,
             xy=(OA_PEAK / 1e6, 11.8), xytext=(48.24, 10.2), fontsize=6.8, color=col['COL2A1'],
             ha='left', va='top',
             arrowprops=dict(arrowstyle='-', lw=0.7, color=col['COL2A1']))

axA.set_ylabel('Knee OA  −log$_{10}$(P)', fontsize=8.5)
axA.set_ylim(0, 14.2)
axA.set_xlim(LO / 1e6, HI / 1e6)
axA.tick_params(labelbottom=False, labelsize=7)
axA.legend(loc='upper left', fontsize=6.8, frameon=False, markerscale=2.2,
           handletextpad=0.4, borderpad=0.2)
for sp in ('top', 'right'):
    axA.spines[sp].set_visible(False)
axA.set_title('(A)  Knee OA association across 12q13 (161,228 cases)',
              fontsize=8.8, loc='left', pad=4)

# gene models
axG.set_ylim(0.3, 2.7)
axG.set_yticks([])
for sp in ('top', 'right', 'left'):
    axG.spines[sp].set_visible(False)
for i, (g, yy) in enumerate((('VDR', 2.0), ('COL2A1', 1.0))):
    s0, e0, strand = gene_span[g]
    c = col[g]
    axG.plot([s0 / 1e6, e0 / 1e6], [yy, yy], color=c, lw=0.8, zorder=2)
    for (es, ee) in exons[g]:
        axG.add_patch(Rectangle((es / 1e6, yy - 0.16), max((ee - es) / 1e6, 0.0008), 0.32,
                                facecolor=c, edgecolor='none', zorder=3))
    axG.text(e0 / 1e6 + 0.012, yy, '%s (%s)' % (g, strand), fontsize=6.8, color=c,
             va='center', ha='left')
    # strand arrow
    xa = s0 / 1e6 + (e0 - s0) / 1e6 * 0.5
    dx = 0.012 if strand == '+' else -0.012
    axG.annotate('', xy=(xa + dx, yy), xytext=(xa, yy),
                 arrowprops=dict(arrowstyle='-|>', lw=0.6, color=c, mutation_scale=6))
axG.set_xlabel('Position on chromosome 12 (Mb, GRCh38)', fontsize=8.5)
axG.tick_params(labelsize=7)

# ---- panel B: eQTL vs OA effect ----
ci = pd.read_csv(os.path.join(TAB, 'coloc_harmonised_input_GCST90566800.csv'))
ci['reg'] = [region(p) for p in ci.pos]
for key, z in (('other', 2), ('VDR', 3), ('COL2A1', 4)):
    m = ci.reg == key
    axB.scatter(ci.b1[m], ci.b2h[m], s=4, c=col[key], alpha=0.55 if key == 'other' else 0.9,
                linewidths=0, zorder=z)
axB.axhline(0, color='k', lw=0.6)
axB.axvline(0, color='k', lw=0.6)

for snp, c, mk, dy in (('rs7975232', col['VDR'], '*', 0.006),
                       ('rs4506714', col['COL2A1'], '*', -0.008),
                       ('rs61919782', col['VDR'], '^', 0.006)):
    row = ci[ci.snp == snp]
    if row.empty:
        continue
    b1, b2 = float(row.b1.iloc[0]), float(row.b2h.iloc[0])
    axB.scatter([b1], [b2], s=64, marker=mk, facecolor=c, edgecolor='k', linewidths=0.5, zorder=6)
    axB.annotate(snp, xy=(b1, b2), xytext=(b1, b2 + dy), fontsize=6.4, ha='center',
                 color=c, zorder=7)

axB.set_xlabel('VDR eQTL effect (SD per allele)', fontsize=7.8)
axB.set_ylabel('Knee OA log OR', fontsize=7.8)
axB.set_xlim(-0.32, 0.32)
axB.set_ylim(-0.105, 0.105)
axB.set_yticks([-0.10, -0.05, 0.0, 0.05, 0.10])
axB.set_yticklabels(['-0.10', '-0.05', '0', '0.05', '0.10'], fontsize=6.8)
axB.text(0.98, 0.03, 'variants beyond the axis\nlimits are not shown',
         transform=axB.transAxes, fontsize=6.0, ha='right', va='bottom', color='#666666')
axB.tick_params(labelsize=6.8)
for sp in ('top', 'right'):
    axB.spines[sp].set_visible(False)
axB.set_title('(B)  Expression and disease effects dissociate', fontsize=8.4, loc='left', pad=4)

# ---- panel C: windowed posteriors ----
sel = ['eQTL peak +/- 10 kb', 'eQTL peak +/- 100 kb', 'OA peak +/- 25 kb',
       'full region 46.8-48.9 Mb']
lab = ['VDR\n±10 kb', 'VDR\n±100 kb', 'OA peak\n±25 kb', 'full\nregion']
sub = w.set_index('window').loc[sel]
xs = np.arange(len(sel))
wd = 0.26
axC.bar(xs - wd, sub.H1.values, wd, color='#4575B4', label='H1  expression only')
axC.bar(xs, sub['PP.H3'].values, wd, color='#F4A582', label='H3  distinct variants')
axC.bar(xs + wd, sub['PP.H4'].values, wd, color='#B2182B', label='H4  shared variant')
axC.set_xticks(xs)
axC.set_xticklabels(lab, fontsize=6.6)
axC.set_ylim(0, 1.06)
axC.set_ylabel('Posterior probability', fontsize=7.8)
axC.tick_params(labelsize=6.8)
for sp in ('top', 'right'):
    axC.spines[sp].set_visible(False)
axC.set_title('(C)  Colocalisation by window', fontsize=8.4, loc='left', pad=4)

handles, labels = axC.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.008),
           ncol=3, frameon=False, fontsize=7.2, handletextpad=0.5, columnspacing=1.6)

os.makedirs(OUT, exist_ok=True)
png = os.path.join(OUT, 'Figure_1_new.png')
fig.savefig(png, dpi=300)
tif = os.path.join(OUT, 'Figure_1_new.tif')
fig.savefig(tif, dpi=300, pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig)

from PIL import Image
im = Image.open(tif)
if im.mode != 'RGB':
    im = im.convert('RGB')
    im.save(tif, dpi=(300, 300), compression='tiff_lzw')
print('\nwritten:', png)
print('written:', tif, Image.open(tif).size, Image.open(tif).mode,
      round(os.path.getsize(tif) / 1e6, 2), 'MB')
