# -*- coding: utf-8 -*-
"""Recompute cis-MR IVW for VDR expression -> knee OA using the 2025 GWAS outcome."""
import sys, io, os
import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

INST = r"C:\Users\liuha\OpenClawVDR\repo\tables\cisMR_LDclumped_discovery_harmonised.csv"
NEW = r"C:\Users\liuha\OpenClawVDR\data\GCST90566800_chr12_46.8-48.9Mb.tsv"
OUT = r"C:\Users\liuha\OpenClawVDR\repo\tables"

COLS = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'odds_ratio',
        'standard_error', 'effect_allele_frequency', 'p_value', 'variant_id', 'ci_upper',
        'ci_lower', 'NCASES', 'NCONTROLS', 'hm_coordinate_conversion', 'hm_code', 'rsid']

ins = pd.read_csv(INST)
oa = pd.read_csv(NEW, sep='\t', names=COLS, header=None, dtype={'chromosome': str})
for c in ('odds_ratio', 'p_value', 'ci_upper', 'ci_lower', 'NCASES', 'NCONTROLS',
          'effect_allele_frequency'):
    oa[c] = pd.to_numeric(oa[c], errors='coerce')
oa['b'] = np.log(oa.odds_ratio)
oa['se'] = (np.log(oa.ci_upper) - np.log(oa.ci_lower)) / (2 * 1.959964)
oa = oa[np.isfinite(oa.se) & (oa.se > 0)]
oa = oa.sort_values('NCASES', ascending=False).drop_duplicates('rsid')

COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
rows = []
for _, r in ins.iterrows():
    snp = r['SNP']
    hit = oa[oa.rsid == snp]
    if hit.empty:
        print('MISSING in new GWAS:', snp)
        continue
    h = hit.iloc[0]
    ea_x = str(r['effect_allele.exposure']).upper()
    oa_x = str(r['other_allele.exposure']).upper()
    ea_o = str(h.effect_allele).upper()
    oa_o = str(h.other_allele).upper()

    if ea_x == ea_o and oa_x == oa_o:
        bo, seo = h.b, h.se
    elif ea_x == oa_o and oa_x == ea_o:
        bo, seo = -h.b, h.se
    elif ea_x == COMP.get(ea_o) and oa_x == COMP.get(oa_o):
        bo, seo = h.b, h.se
    elif ea_x == COMP.get(oa_o) and oa_x == COMP.get(ea_o):
        bo, seo = -h.b, h.se
    else:
        print('ALLELE MISMATCH', snp, ea_x, oa_x, '|', ea_o, oa_o)
        continue

    bx = r['beta.exposure']
    sex = r['se.exposure']
    wald = bo / bx
    wald_se = abs(seo / bx)
    rows.append(dict(snp=snp, bx=bx, sex=sex, bo=bo, seo=seo, p_out=h.p_value,
                     wald=wald, wald_se=wald_se,
                     F=(bx / sex) ** 2))

d = pd.DataFrame(rows)
d['w'] = 1 / d.wald_se ** 2
print('\n=== instruments ===')
print(d[['snp', 'bx', 'sex', 'bo', 'seo', 'p_out', 'wald', 'wald_se', 'F']].to_string(index=False))

def ivw(dd, random=True):
    w = 1 / dd.wald_se ** 2
    b = (w * dd.wald).sum() / w.sum()
    se_fixed = np.sqrt(1 / w.sum())
    q = (w * (dd.wald - b) ** 2).sum()
    dfree = len(dd) - 1
    p_q = 1 - stats.chi2.cdf(q, dfree) if dfree > 0 else np.nan
    if random and dfree > 0:
        phi = max(1.0, q / dfree)
        se = se_fixed * np.sqrt(phi)
    else:
        se = se_fixed
    z = b / se
    return dict(b=b, se=se, OR=np.exp(b), lo=np.exp(b - 1.959964 * se),
                hi=np.exp(b + 1.959964 * se), p=2 * (1 - stats.norm.cdf(abs(z))),
                Q=q, Qdf=dfree, Qp=p_q)

all4 = ivw(d)
gw3 = ivw(d[d.p_out < 5e-8]) if (d.p_out < 5e-8).sum() >= 2 else None

print('\n=== IVW, all 4 instruments (new outcome) ===')
print('  OR = %.3f (95%% CI %.3f-%.3f), P = %.4f' % (all4['OR'], all4['lo'], all4['hi'], all4['p']))
print('  Cochran Q = %.2f on %d df, P = %.3f' % (all4['Q'], all4['Qdf'], all4['Qp']))
b_fixed = (d.w * d.wald).sum() / d.w.sum()
se_f = np.sqrt(1 / d.w.sum())
print('  fixed-effect: OR = %.3f (%.3f-%.3f), P = %.4f' % (
    np.exp(b_fixed), np.exp(b_fixed - 1.96 * se_f), np.exp(b_fixed + 1.96 * se_f),
    2 * (1 - stats.norm.cdf(abs(b_fixed / se_f)))))

if gw3:
    print('\n=== IVW, genome-wide significant subset ===')
    print('  OR = %.3f (95%% CI %.3f-%.3f), P = %.4f' % (gw3['OR'], gw3['lo'], gw3['hi'], gw3['p']))
    print('  Cochran Q = %.2f on %d df, P = %.3f' % (gw3['Q'], gw3['Qdf'], gw3['Qp']))

top = d.loc[d.p_out.idxmin()]
print('\n=== single strongest cis variant (Wald) ===')
print('  %s: OR = %.3f (%.3f-%.3f), P = %.4f' % (
    top.snp, np.exp(top.wald), np.exp(top.wald - 1.96 * top.wald_se),
    np.exp(top.wald + 1.96 * top.wald_se),
    2 * (1 - stats.norm.cdf(abs(top.wald / top.wald_se)))))

d.to_csv(os.path.join(OUT, 'cisMR_GCST90566800_instruments.csv'), index=False)
print('\nwritten cisMR_GCST90566800_instruments.csv')
