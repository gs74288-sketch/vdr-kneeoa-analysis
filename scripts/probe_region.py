import json, urllib.request, urllib.error
tok=None
for line in open(r"C:\Users\liuha\Documents\.Renviron",encoding="utf-8",errors="replace"):
    if "OPENGWAS_JWT" in line: tok=line.split("=",1)[1].strip().strip("'").strip('"')
def call(payload):
    req=urllib.request.Request("https://api.opengwas.io/api/associations",data=json.dumps(payload).encode(),method="POST",
        headers={"Authorization":"Bearer "+tok,"Content-Type":"application/json"})
    try:
        with urllib.request.urlopen(req,timeout=180) as r: return json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        return {"__err__":e.code,"body":e.read().decode()[:300]}
r1=call({"id":["eqtl-a-ENSG00000111424"],"variant":["12:46800000-48900000"]})
print("eQTL region rows:", len(r1) if isinstance(r1,list) else r1)
if isinstance(r1,list) and r1: print("cols:", list(r1[0].keys())); print("first:", r1[0])
r2=call({"id":["ebi-a-GCST007090"],"variant":["12:46800000-48900000"]})
print("OA region rows:", len(r2) if isinstance(r2,list) else r2)
