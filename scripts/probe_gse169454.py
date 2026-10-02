import gzip, io, urllib.request, itertools, time

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169454/suppl/GSE169454_OA_gene_expression_matrix.tsv.gz"
req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
resp = urllib.request.urlopen(req, timeout=120)
gz = gzip.GzipFile(fileobj=resp)
txt = io.TextIOWrapper(gz, encoding="utf-8", newline="")
t0 = time.time()
for k, line in enumerate(itertools.islice(txt, 12)):
    print(f"L{k}: {line[:300]}")
print("elapsed", round(time.time()-t0, 1))
