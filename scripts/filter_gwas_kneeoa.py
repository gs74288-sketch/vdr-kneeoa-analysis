import sys, gzip, io, urllib.request, csv, time

URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007090/harmonised/30664745-GCST007090-EFO_0004616.h.tsv.gz"
OUT = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\data\GCST007090_chr12_VDR_region.tsv"
LOOP = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\logs\gwas_filter.log"
CHR_WINDOW = ("12", 47500000, 49000000)

def log(m):
    with open(LOOP, "a", encoding="utf-8") as f:
        f.write(m + "\n")

log(f"start {time.ctime()}")
t0 = time.time()
req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
resp = urllib.request.urlopen(req, timeout=120)
gz = gzip.GzipFile(fileobj=resp)
text = io.TextIOWrapper(gz, encoding="utf-8", newline="")
r = csv.reader(text, delimiter="\t")
header = next(r)
log("HEADER: " + "|".join(header))
# locate columns
def idx(name, alt=None):
    for i, h in enumerate(header):
        if h == name:
            return i
    return None
ci_chr = idx("hm_chrom")
ci_pos = idx("hm_pos")
if ci_chr is None or ci_pos is None:
    log("ERROR: harmonised columns not found; header=" + "|".join(header))
    sys.exit(2)
chrom, lo, hi = CHR_WINDOW
kept = 0
with open(OUT, "w", encoding="utf-8", newline="") as fo:
    w = csv.writer(fo, delimiter="\t")
    w.writerow(header)
    n = 0
    for row in r:
        n += 1
        if n % 1000000 == 0:
            log(f"scanned {n} rows, kept {kept}, {time.time()-t0:.0f}s")
        if len(row) <= max(ci_chr, ci_pos):
            continue
        if row[ci_chr] == chrom:
            try:
                p = int(row[ci_pos])
            except ValueError:
                continue
            if lo <= p <= hi:
                w.writerow(row)
                kept += 1
    log(f"done scanned={n} kept={kept} elapsed={time.time()-t0:.0f}s")
log("END " + time.ctime())
