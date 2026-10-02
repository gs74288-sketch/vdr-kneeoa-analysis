import json, urllib.request, urllib.error
BASE = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002"
tok = None
for line in open(r"C:\Users\liuha\Documents\.Renviron", encoding="utf-8", errors="replace"):
    if "OPENGWAS_JWT" in line:
        tok = line.split("=", 1)[1].strip().strip("'").strip('"')
rs = [l.strip() for l in open(BASE + r"\data\ens_vdr_region_rsids.txt") if l.strip()]
for n in (1000, 500, 100, 10):
    payload = {"id": ["eqtl-a-ENSG00000111424"], "variants": rs[:n]}
    req = urllib.request.Request("https://api.opengwas.io/api/associations",
                                 data=json.dumps(payload).encode(), method="POST",
                                 headers={"Authorization": "Bearer " + tok,
                                          "Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            d = json.loads(r.read().decode())
            print(f"batch {n}: OK rows={len(d)}")
    except urllib.error.HTTPError as e:
        body = e.read().decode()[:300]
        print(f"batch {n}: HTTP {e.code} -> {body}")
    except Exception as e:
        print(f"batch {n}: ERR {e}")
