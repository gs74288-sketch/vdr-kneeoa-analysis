# -*- coding: utf-8 -*-
"""Minimal tabix/BGZF remote range reader (no pysam needed).

Reads the TBI index, finds chunks overlapping a region, fetches only those
byte ranges over HTTP, and decompresses the BGZF blocks.
"""
import sys, io, os, gzip, zlib, struct, urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

BASE = ('https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
        'GCST90566001-GCST90567000/GCST90566800/harmonised/GCST90566800.h.tsv.gz')
TBI = BASE + '.tbi'


def get_range(url, start, end):
    req = urllib.request.Request(url, headers={'Range': 'bytes=%d-%d' % (start, end),
                                              'User-Agent': 'openclaw'})
    return urllib.request.urlopen(req, timeout=180).read()


# ---- fetch and parse the index -------------------------------------------------
raw = get_range(TBI, 0, 1 << 62) if False else get_range(TBI, 0, 2000000)
tbi = gzip.decompress(raw)
assert tbi[:4] == b'TBI\x01', 'bad magic %r' % tbi[:4]

off = 4
n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack_from('<8i', tbi, off)
off += 32
names_blob = tbi[off:off + l_nm]
names = names_blob.rstrip(b'\x00').decode().split('\x00')
off += l_nm
print('n_ref=%d format=%d col_seq=%d col_beg=%d col_end=%d' % (n_ref, fmt, col_seq, col_beg, col_end))
print('refs:', names[:25])

TARGET = '12'
ri = names.index(TARGET)

for i in range(n_ref):
    n_bin = struct.unpack_from('<i', tbi, off)[0]
    off += 4
    bins = {}
    for _ in range(n_bin):
        b, n_chunk = struct.unpack_from('<Ii', tbi, off)
        off += 8
        ch = [struct.unpack_from('<QQ', tbi, off + 16 * k) for k in range(n_chunk)]
        off += 16 * n_chunk
        bins[b] = ch
    n_intv = struct.unpack_from('<i', tbi, off)[0]
    off += 4 + 8 * n_intv
    if i == ri:
        my_bins = bins
        break

print('chr12 bins:', len(my_bins))

# ---- reg2bins ------------------------------------------------------------------
BEG, END = 46800000, 48900000


def reg2bins(beg, end):
    end -= 1
    out = [0]
    for k in range(1 + (beg >> 26), 1 + (end >> 26) + 1):
        out.append(k)
    for k in range(1 + (beg >> 23), 1 + (end >> 23) + 1):
        out.append(9 + k)
    for k in range(1 + (beg >> 20), 1 + (end >> 20) + 1):
        out.append(73 + k)
    for k in range(1 + (beg >> 17), 1 + (end >> 17) + 1):
        out.append(585 + k)
    for k in range(1 + (beg >> 14), 1 + (end >> 14) + 1):
        out.append(4681 + k)
    return out


want = reg2bins(BEG, END)
chunks = []
for b in want:
    if b in my_bins:
        chunks.extend(my_bins[b])

# merge overlapping / adjacent chunks
chunks.sort()
merged = []
for c in chunks:
    if merged and c[0] <= merged[-1][1]:
        merged[-1] = (merged[-1][0], max(merged[-1][1], c[1]))
    else:
        merged.append(list(c))
print('chunks:', len(chunks), '-> merged', len(merged))

# virtual offset -> byte offset (block start)
starts = [(c[0] >> 16, (c[1] >> 16) + 65536) for c in merged]
total = sum(b - a for a, b in starts)
print('bytes to fetch: %s' % format(total, ','))

blobs = []
for a, b in starts:
    blobs.append(get_range(BASE, a, b - 1))
data = b''.join(blobs)

# BGZF is concatenated gzip members
try:
    txt = gzip.decompress(data).decode('utf-8', 'replace')
except Exception as e:
    print('gzip.decompress failed (%s); streaming members' % e)
    out = []
    rest = data
    while rest:
        d = zlib.decompressobj(31)
        try:
            out.append(d.decompress(rest))
        except Exception:
            break
        rest = d.unused_data
        if not rest:
            break
    txt = b''.join(out).decode('utf-8', 'replace')

lines = txt.split('\n')
print('lines:', len(lines))

keep = []
for ln in lines:
    f = ln.split('\t')
    if len(f) < 9:
        continue
    if f[0] != TARGET:
        continue
    try:
        pos = int(f[1])
    except ValueError:
        continue
    if BEG <= pos <= END:
        keep.append(ln)

print('kept rows:', len(keep))
out_path = r"C:\Users\liuha\OpenClawVDR\data\GCST90566800_chr12_46.8-48.9Mb.tsv"
hdr = lines[0] if lines and lines[0].startswith('chromosome') else None
with open(out_path, 'w', encoding='utf-8', newline='') as fh:
    if hdr:
        fh.write(hdr + '\n')
    for ln in keep:
        fh.write(ln + '\n')
print('written ->', out_path, os.path.getsize(out_path), 'bytes')
