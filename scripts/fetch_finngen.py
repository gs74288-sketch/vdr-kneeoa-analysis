# -*- coding: utf-8 -*-
"""Fast stream-filter of FinnGen R11 knee-OA summary stats for the VDR instrument SNPs.
Uses a single compiled regex per line (much faster than 11 substring scans)."""
import io, gzip, re, urllib.request, time, os

SNPS = ["rs10783222","rs12943633","rs149110519","rs2594836","rs28498283","rs35979828",
        "rs372883","rs4142441","rs58067800","rs7485057","rs7975232"]
PAT = re.compile("|".join(SNPS))
OUTDIR = r"C:\Users\liuha\OpenClawVDR\finngen"
os.makedirs(OUTDIR, exist_ok=True)

TARGETS = {
  "M13_ARTHROSIS_KNEE": "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_M13_ARTHROSIS_KNEE.gz",
  "M13_ARTHROSIS_KNEE_PRIM_KNEESURG": "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_M13_ARTHROSIS_KNEE_PRIM_KNEESURG.gz",
}

for tag, url in TARGETS.items():
    t0 = time.time()
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    resp = urllib.request.urlopen(req, timeout=600)
    gz = gzip.GzipFile(fileobj=resp)
    txt = io.TextIOWrapper(gz, encoding="utf-8", errors="replace", newline="")
    header = txt.readline().rstrip("\n").split("\t")
    ri = None
    for i, h in enumerate(header):
        if h.strip().lower() in ("rsids", "rsid", "snp"):
            ri = i
    if ri is None:
        ri = 4
    print(tag, "| header:", header, "| rsid col:", ri, flush=True)
    out = []
    n = 0
    for line in txt:
        n += 1
        if PAT.search(line) is None:
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > ri and (set(f[ri].split(",")) & set(SNPS)):
            out.append(f)
    with open(os.path.join(OUTDIR, tag + "_instrument_snps.tsv"), "w", encoding="utf-8", newline="") as fo:
        fo.write("\t".join(header) + "\n")
        for f in out:
            fo.write("\t".join(f) + "\n")
    print(f"{tag}: scanned {n} rows, matched {len(out)}, {time.time()-t0:.0f}s", flush=True)

print("DONE_FINNGEN")

