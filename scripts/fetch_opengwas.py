# -*- coding: utf-8 -*-
"""Fetch VDR cis-eQTL and knee-OA associations via the OpenGWAS API in large batches."""
import json, urllib.request, time, csv, sys

BASE = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002"
ENV = r"C:\Users\liuha\Documents\.Renviron"
LOG = BASE + r"\logs\fetch_opengwas.log"

def log(m):
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(m + "\n")

# token
tok = None
for line in open(ENV, encoding="utf-8", errors="replace"):
    if "OPENGWAS_JWT" in line:
        tok = line.split("=", 1)[1].strip().strip("'").strip('"')
log(f"token present={bool(tok)} len={len(tok or '')} start={time.ctime()}")

def post(path, payload, tries=6):
    url = "https://api.opengwas.io/api/" + path
    data = json.dumps(payload).encode()
    for t in range(tries):
        try:
            req = urllib.request.Request(url, data=data, method="POST",
                                         headers={"Authorization": "Bearer " + tok,
                                                  "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            log(f"  retry {t+1}/{tries} {path}: {e}")
            time.sleep(3 * (t + 1))
    return None

rs = [l.strip() for l in open(BASE + r"\data\ens_vdr_region_rsids.txt") if l.strip()]
log(f"region rsIDs: {len(rs)}")

BATCH = 500
rows = []
for i in range(0, len(rs), BATCH):
    batch = rs[i:i + BATCH]
    res = post("associations", {"id": ["eqtl-a-ENSG00000111424"], "variants": batch})
    if res is None:
        log(f"  batch {i//BATCH+1} FAILED"); continue
    for r in res:
        rows.append(r)
    if (i // BATCH) % 10 == 0:
        log(f"  eQTL batch {i//BATCH+1}/{(len(rs)+BATCH-1)//BATCH} total_rows={len(rows)}")

log(f"eQTL rows fetched: {len(rows)}")
if rows:
    cols = list(rows[0].keys())
    with open(BASE + r"\data\eqtl_opengwas.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)

snp_ids = sorted({r["rsid"] for r in rows if r.get("rsid")})
log(f"eQTL unique rsids: {len(snp_ids)}")

rows2 = []
for i in range(0, len(snp_ids), BATCH):
    batch = snp_ids[i:i + BATCH]
    res = post("associations", {"id": ["ebi-a-GCST007090"], "variants": batch})
    if res is None:
        log(f"  OA batch {i//BATCH+1} FAILED"); continue
    rows2 += res
log(f"OA rows fetched: {len(rows2)}")
if rows2:
    cols = list(rows2[0].keys())
    with open(BASE + r"\data\oa_opengwas.tsv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows2)
log("DONE_FETCH")
