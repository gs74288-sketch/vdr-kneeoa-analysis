# -*- coding: utf-8 -*-
"""Rebuild the susie input with strict allele matching (no strand flips,
palindromic sites dropped) and verify that the LD matrix is consistent
with the summary statistics before running susie.
"""
import sys, io, os, gzip, struct, zlib, urllib.request
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

BASE = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
VCF = BASE + 'ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'
LD = r"C:\Users\liuha\OpenClawVDR\data\ld1000g"
TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"
MAP = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\data\oa_region.tsv"
OUTDIR = os.path.join(TAB, 'susie_OAwindow')

L38, H38 = 47_960_000, 48_060_000


def get(url, start=None, end=None):
    h = {'User-Agent': 'openclaw'}
    if start is not None:
        h['Range'] = 'bytes=%d-%d' % (start, end)
    return urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=300).read()


def bgzf_all(blob):
    out, rest = [], blob
    while rest:
        d = zlib.decompressobj(31)
        try:
            out.append(d.decompress(rest))
        except Exception:
            break
        rest = d.unused_data
        if not rest:
            break
    return b''.join(out)


ci = pd.read_csv(os.path.join(TAB, 'coloc_harmonised_input_GCST90566800.csv'))
win = ci[(ci.pos >= L38) & (ci.pos <= H38)].copy()
m = pd.read_csv(MAP, sep='\t').drop_duplicates('rsid').set_index('rsid')
win['pos37'] = win.snp.map(m['position'].to_dict())
win = win[win.pos37.notna()].copy()
win['pos37'] = win.pos37.astype(int)
print('variants in window:', len(win))

L37, H37 = int(win.pos37.min()), int(win.pos37.max())

tbi = gzip.decompress(open(os.path.join(LD, 'chr12.vcf.gz.tbi'), 'rb').read())
off = 4
vals = struct.unpack_from('<8i', tbi, off)
n_ref, l_nm = vals[0], vals[7]
off += 32
names = tbi[off:off + l_nm].rstrip(b'\x00').decode().split('\x00')
off += l_nm


def reg2bins(beg, end):
    end -= 1
    out = [0]
    for sh, base in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        for k in range(1 + (beg >> sh), 1 + (end >> sh) + 1):
            out.append(base + k)
    return out


for i in range(n_ref):
    n_bin = struct.unpack_from('<i', tbi, off)[0]
    off += 4
    bins = {}
    for _ in range(n_bin):
        b, nc = struct.unpack_from('<Ii', tbi, off)
        off += 8
        bins[b] = [struct.unpack_from('<QQ', tbi, off + 16 * k) for k in range(nc)]
        off += 16 * nc
    n_intv = struct.unpack_from('<i', tbi, off)[0]
    off += 4 + 8 * n_intv
    if names[i] == '12':
        my = bins
        break

chunks = []
for b in reg2bins(L37, H37):
    if b in my:
        chunks.extend(my[b])
chunks.sort()
merged = []
for c in chunks:
    if merged and c[0] <= merged[-1][1]:
        merged[-1][1] = max(merged[-1][1], c[1])
    else:
        merged.append(list(c))
rng = [(c[0] >> 16, (c[1] >> 16) + 65536) for c in merged]
txt = bgzf_all(b''.join(get(VCF, a, b - 1) for a, b in rng)).decode('utf-8', 'replace')

panel = pd.read_csv(os.path.join(LD, 'panel.txt'), sep='\t')
eur = panel.loc[panel['super_pop'] == 'EUR', 'sample'].tolist()
samples = open(os.path.join(LD, 'samples.txt')).read().split('\n')
idx = {s: i for i, s in enumerate(samples)}
eur_idx = [idx[s] for s in eur if s in idx]
print('EUR columns:', len(eur_idx))

recs = {}
for l in txt.split('\n'):
    if not l or l.startswith('#'):
        continue
    f = l.split('\t')
    if len(f) < 10 or f[0] != '12':
        continue
    try:
        p = int(f[1])
    except ValueError:
        continue
    if L37 <= p <= H37:
        recs[p] = (f[3], f[4], f[9:])

PAL = {('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')}
counts = {'direct': 0, 'swapped': 0, 'palindromic': 0, 'absent': 0, 'allele_mismatch': 0}
rows, ids, meta = [], [], []

for _, r in win.iterrows():
    p = int(r.pos37)
    if p not in recs:
        counts['absent'] += 1
        continue
    ref, alt, gts = recs[p]
    if len(ref) != 1 or len(alt) != 1:
        counts['allele_mismatch'] += 1
        continue
    if (ref.upper(), alt.upper()) in PAL:
        counts['palindromic'] += 1
        continue
    ea = str(r.ea2).upper()
    oa = str(r.oa2).upper()
    if ref.upper() == oa and alt.upper() == ea:
        flip = False
        counts['direct'] += 1
    elif ref.upper() == ea and alt.upper() == oa:
        flip = True
        counts['swapped'] += 1
    else:
        counts['allele_mismatch'] += 1
        continue

    vec = np.full(len(eur_idx), np.nan)
    for j, c_ in enumerate(eur_idx):
        g = gts[c_].split(':')[0]
        if g in ('0|0', '0/0'):
            v = 0
        elif g in ('0|1', '1|0', '0/1', '1/0'):
            v = 1
        elif g in ('1|1', '1/1'):
            v = 2
        else:
            continue
        vec[j] = (2 - v) if flip else v
    rows.append(vec)
    ids.append(r.snp)
    meta.append((r.snp, p, r.b1, r.se1, r.b2h, r.se2, r.p1, r.p2, ea, oa))

print('\nmatch modes:', counts)

G = np.array(rows, dtype=float)
good = (np.isnan(G).mean(1) < 0.02) & (np.nanstd(G, axis=1) > 0.05)
G = G[good]
ids = [i for i, k in zip(ids, good) if k]
meta = [x for x, k in zip(meta, good) if k]
print('variants kept:', G.shape[0])

mu = np.nanmean(G, axis=1, keepdims=True)
G = np.where(np.isnan(G), mu, G)
Z = G - G.mean(1, keepdims=True)
sd = Z.std(1, keepdims=True)
sd[sd == 0] = 1
Z = Z / sd
R = (Z @ Z.T) / Z.shape[1]
np.fill_diagonal(R, 1.0)

w, V = np.linalg.eigh(R)
if w.min() < 1e-6:
    w = np.clip(w, 1e-6, None)
    R = V @ np.diag(w) @ V.T
    dd = np.sqrt(np.diag(R))
    R = R / np.outer(dd, dd)
print('LD eigen min (corrected): %.2e' % np.linalg.eigvalsh(R).min())

df = pd.DataFrame(meta, columns=['snp', 'pos37', 'b1', 'se1', 'b2h', 'se2', 'p1', 'p2', 'ea', 'oa'])
df.to_csv(os.path.join(OUTDIR, 'susie_input.csv'), index=False)
np.save(os.path.join(OUTDIR, 'susie_ld.npy'), R)
open(os.path.join(OUTDIR, 'susie_ld_snps.txt'), 'w').write('\n'.join(ids) + '\n')
pd.DataFrame(R, index=ids, columns=ids).to_csv(os.path.join(OUTDIR, 'susie_ld.csv'))

# ---- consistency diagnostic ----
z = (df.b2h / df.se2).values
lead = int(np.argmax(np.abs(z)))
print('\nlead OA SNP: %s  z = %.2f' % (df.snp[lead], z[lead]))
nb = [i for i in np.argsort(-np.abs(R[lead])) if i != lead][:10]
agree = sum(1 for i in nb if abs(R[lead, i]) > 0.15
            and np.sign(R[lead, i]) == np.sign(z[i] / z[lead]))
tot = sum(1 for i in nb if abs(R[lead, i]) > 0.15)
print('  neighbour sign agreement: %d/%d' % (agree, tot))

top = np.argsort(-np.abs(z))[:200]
A = R[np.ix_(top, top)]
zz = np.outer(z[top], z[top])
iu = np.triu_indices(len(top), 1)
print('  corr(LD r, z_i*z_j) top-200 pairs: %.3f' % np.corrcoef(A[iu], zz[iu])[0, 1])

print('\nvariants: %d  |  OA p<5e-8: %d' % (len(df), (df.p2 < 5e-8).sum()))
print('written to', OUTDIR)
