# -*- coding: utf-8 -*-
"""Quantify: VDR body vs COL2A1 body - where does the knee OA signal actually sit?"""
import sys, io
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

D = r"C:\Users\liuha\OpenClawVDR\data\GCST90566800_chr12_46.8-48.9Mb.tsv"
COLS = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'odds_ratio',
        'standard_error', 'effect_allele_frequency', 'p_value', 'variant_id', 'ci_upper',
        'ci_lower', 'NCASES', 'NCONTROLS', 'hm', 'hm_code', 'rsid']
d = pd.read_csv(D, sep='\t', names=COLS, header=None, dtype={'chromosome': str})
for c in ('base_pair_location', 'odds_ratio', 'p_value'):
    d[c] = pd.to_numeric(d[c], errors='coerce')
d = d.sort_values('NCASES', ascending=False).drop_duplicates('rsid')
d = d[d.p_value.notna()]

VDR = (47_841_537, 47_943_048)
COL2A1 = (47_972_967, 48_004_554)


def summarise(lo, hi, label, pad=0):
    lo_, hi_ = lo - pad, hi + pad
    s = d[(d.base_pair_location >= lo_) & (d.base_pair_location <= hi_)]
    if s.empty:
        print("  %-28s no variants" % label)
        return
    gw = int((s.p_value < 5e-8).sum())
    top = s.loc[s.p_value.idxmin()]
    print("  %-28s n=%5d  GW-sig=%3d  best P=%.2e  (%s @ %s)" % (
        label, len(s), gw, top.p_value, top.rsid, format(int(top.base_pair_location), ',')))
    return s


print("=== knee OA signal by gene body (GCST90566800, 163,447 cases) ===")
summarise(*VDR, 'VDR body')
summarise(*VDR, 'VDR body +/-10 kb', pad=10_000)
summarise(*VDR, 'VDR body +/-50 kb', pad=50_000)
print()
summarise(*COL2A1, 'COL2A1 body')
summarise(*COL2A1, 'COL2A1 body +/-10 kb', pad=10_000)
summarise(*COL2A1, 'COL2A1 body +/-25 kb', pad=25_000)
print()
summarise(46_800_000, 48_900_000, 'whole interval')

# how many of the 77 GW-significant SNPs sit closer to COL2A1 than to VDR?
gw = d[d.p_value < 5e-8].copy()
vdr_mid = (VDR[0] + VDR[1]) / 2
col_mid = (COL2A1[0] + COL2A1[1]) / 2
gw['d_vdr'] = (gw.base_pair_location - vdr_mid).abs()
gw['d_col'] = (gw.base_pair_location - col_mid).abs()
print("\nGW-significant SNPs: %d" % len(gw))
print("  closer to COL2A1 than VDR: %d" % int((gw.d_col < gw.d_vdr).sum()))
print("  span: %s - %s" % (format(int(gw.base_pair_location.min()), ','),
                           format(int(gw.base_pair_location.max()), ',')))
print("  median distance to COL2A1 midpoint: %s bp" % format(int(gw.d_col.median()), ','))

# distance of each GW SNP to the nearest of the two genes' bodies
def dist_to_body(pos, lo, hi):
    if pos < lo:
        return lo - pos
    if pos > hi:
        return pos - hi
    return 0

gw['in_vdr'] = [dist_to_body(p, *VDR) for p in gw.base_pair_location]
gw['in_col'] = [dist_to_body(p, *COL2A1) for p in gw.base_pair_location]
print("  within VDR body: %d" % int((gw.in_vdr == 0).sum()))
print("  within COL2A1 body: %d" % int((gw.in_col == 0).sum()))
print("  within 25 kb of COL2A1: %d" % int((gw.in_col <= 25_000).sum()))
print("  within 25 kb of VDR: %d" % int((gw.in_vdr <= 25_000).sum()))
