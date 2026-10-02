# -*- coding: utf-8 -*-
import json, urllib.request, urllib.error, csv
BASE = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002"
tok = None
for line in open(r"C:\Users\liuha\Documents\.Renviron", encoding="utf-8", errors="replace"):
    if "OPENGWAS_JWT" in line:
        tok = line.split("=", 1)[1].strip().strip("'").strip('"')

def call(ids, region):
    req = urllib.request.Request("https://api.opengwas.io/api/associations",
        data=json.dumps({"id": ids, "variant": [region]}).encode(), method="POST",
        headers={"Authorization": "Bearer " + tok, "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.loads(r.read().decode())

REGION = "12:46800000-48900000"
for ids, out in [(["eqtl-a-ENSG00000111424"], "eqtl_region.tsv"),
                 (["ebi-a-GCST007090"], "oa_region.tsv")]:
    d = call(ids, REGION)
    cols = list(d[0].keys())
    with open(BASE + r"\data\\" + out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(d)
    print(out, len(d), "rows, cols:", cols)
