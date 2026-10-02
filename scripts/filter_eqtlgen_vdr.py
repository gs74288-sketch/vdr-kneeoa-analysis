import sys, gzip, io, urllib.request, csv, time

URL = "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
OUT = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\data\eQTLGen_VDR_cis.tsv"
LOOP = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\logs\eqtlgen_filter.log"
GENE = "ENSG00000111424"  # VDR

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
gi = header.index("Gene")
kept = 0
with open(OUT, "w", encoding="utf-8", newline="") as fo:
    w = csv.writer(fo, delimiter="\t")
    w.writerow(header)
    n = 0
    for row in r:
        n += 1
        if n % 1000000 == 0:
            log(f"scanned {n} rows, kept {kept}, {time.time()-t0:.0f}s")
        if len(row) > gi and row[gi] == GENE:
            w.writerow(row)
            kept += 1
    log(f"done scanned={n} kept={kept} elapsed={time.time()-t0:.0f}s")
log("END " + time.ctime())
