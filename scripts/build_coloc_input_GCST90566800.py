# -*- coding: utf-8 -*-
"""Build a new coloc input: eQTLGen VDR cis (dataset 1) + GCST90566800 knee OA (dataset 2)."""
import sys, io, os
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

OLD = r"C:\Users\liuha\OpenClawVDR\repo\tables\coloc_harmonised_input.csv"
NEW = r"C:\Users\liuha\OpenClawVDR\data\GCST90566800_chr12_46.8-48.9Mb.tsv"
OUT = r"C:\Users\liuha\OpenClawVDR\repo\tables\coloc_harmonised_input_GCST90566800.csv"

COLS = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'odds_ratio',
        'standard_error', 'effect_allele_frequency', 'p_value', 'variant_id', 'ci_upper',
        'ci_lower', 'NCASES', 'NCONTROLS', 'hm_coordinate_conversion', 'hm_code', 'rsid']

eq = pd.read_csv(OLD)
eq = eq[['snp', 'ea1', 'oa1', 'eaf1', 'b1', 'se1', 'p1']].copy()
eq.columns = ['snp', 'ea', 'oa', 'eaf', 'b_eqtl', 'se_eqtl', 'p_eqtl']
print('eQTL variants:', len(eq))

oa = pd.read_csv(NEW, sep='\t', names=COLS, header=None, dtype={'chromosome': str})
for c in ('odds_ratio', 'p_value', 'effect_allele_frequency', 'ci_upper', 'ci_lower',
          'NCASES', 'NCONTROLS'):
    oa[c] = pd.to_numeric(oa[c], errors='coerce')

# per-SNP SE from the 95% CI (standard_error column is NA)
oa['b_oa'] = np.log(oa.odds_ratio)
oa['se_oa'] = (np.log(oa.ci_upper) - np.log(oa.ci_lower)) / (2 * 1.959964)
bad = ~np.isfinite(oa.se_oa) | (oa.se_oa <= 0)
print('rows without usable CI:', int(bad.sum()))
oa = oa[~bad].copy()
oa['eaf_oa'] = oa.effect_allele_frequency

# collapse duplicate rsids: keep the row with the largest case count
oa = oa.sort_values('NCASES', ascending=False).drop_duplicates('rsid')
print('OA unique rsids:', len(oa))

m = eq.merge(oa[['rsid', 'effect_allele', 'other_allele', 'eaf_oa', 'b_oa', 'se_oa',
                 'p_value', 'NCASES', 'NCONTROLS', 'base_pair_location']],
             left_on='snp', right_on='rsid', how='inner')
print('after rsID merge:', len(m))

COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
PALINDO = {('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')}
AMBIG = np.zeros(len(m), dtype=bool)


def norm(s):
    return s.astype(str).str.upper()


ea_e, oa_e = norm(m.ea), norm(m.oa)
ea_o, oa_o = norm(m.effect_allele), norm(m.other_allele)

direct = (ea_e == ea_o) & (oa_e == oa_o)
swapped = (ea_e == oa_o) & (oa_e == ea_o)
comp_direct = (ea_e == ea_o.map(COMP)) & (oa_e == oa_o.map(COMP))
comp_swapped = (ea_e == oa_o.map(COMP)) & (oa_e == ea_o.map(COMP))

b2 = np.full(len(m), np.nan)
flip = np.zeros(len(m), dtype=bool)
b2[direct.values] = m.b_oa.values[direct.values]
b2[swapped.values] = -m.b_oa.values[swapped.values]
flip[swapped.values] = True
cd = comp_direct.fillna(False).values
cs = comp_swapped.fillna(False).values
b2[cd] = m.b_oa.values[cd]
b2[cs] = -m.b_oa.values[cs]
flip[cs] = True

m['b_oa_h'] = b2
m['flip'] = flip
m['aligned'] = np.isfinite(b2)

pal = np.zeros(len(m), dtype=bool)
for i in range(len(m)):
    if (ea_e.iloc[i], oa_e.iloc[i]) in PALINDO:
        pal[i] = True
m['palindromic'] = pal

print('  aligned:', int(m.aligned.sum()))
print('  not alignable:', int((~m.aligned).sum()))
print('  palindromic:', int(pal.sum()))

keep = m[m.aligned].copy()
# drop palindromic SNPs with intermediate frequency (strand-ambiguous)
keep = keep[~((keep.palindromic) & (keep.eaf.between(0.4, 0.6)))]
print('  after palindromic filter:', len(keep))

keep['maf'] = np.minimum(keep.eaf, 1 - keep.eaf)
keep = keep[keep.maf > 0.01].copy()
print('  after MAF>0.01:', len(keep))

out = pd.DataFrame({
    'snp': keep.snp,
    'ea1': keep.ea, 'oa1': keep.oa, 'eaf1': keep.eaf,
    'b1': keep.b_eqtl, 'se1': keep.se_eqtl, 'p1': keep.p_eqtl,
    'ea2': keep.effect_allele, 'oa2': keep.other_allele, 'eaf2': keep.eaf_oa,
    'b2': keep.b_oa, 'se2': keep.se_oa, 'p2': keep.p_value,
    'flip': keep.flip.values, 'b2h': keep.b_oa_h, 'maf': keep.maf,
    'pos': keep.base_pair_location, 'ncase': keep.NCASES, 'ncontrol': keep.NCONTROLS,
})
out.to_csv(OUT, index=False)
print('\nwritten', OUT, os.path.getsize(OUT), 'bytes', out.shape)

print('\nN in region: cases %.0f-%.0f (median %.0f), controls %.0f-%.0f' % (
    out.ncase.min(), out.ncase.max(), out.ncase.median(),
    out.ncontrol.min(), out.ncontrol.max()))
print('OA p in merged set: min %.3g  median %.3g' % (out.p2.min(), out.p2.median()))
print('eQTL p: min %.3g  median %.3g' % (out.p1.min(), out.p1.median()))
print('GW-significant OA SNPs in merged set:', int((out.p2 < 5e-8).sum()))
