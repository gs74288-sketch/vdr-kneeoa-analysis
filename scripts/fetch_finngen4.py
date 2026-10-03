import io,gzip,re,urllib.request,time,os
SNPS=["rs61919782","rs7975232","rs4506714","rs191597566"]
PAT=re.compile("|".join(SNPS))
OUT=r"C:\Users\liuha\OpenClawVDR\finngen"
T={"M13_ARTHROSIS_KNEE":"https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_M13_ARTHROSIS_KNEE.gz",
   "M13_ARTHROSIS_KNEE_PRIM_KNEESURG":"https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_M13_ARTHROSIS_KNEE_PRIM_KNEESURG.gz"}
for tag,url in T.items():
    t0=time.time(); req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0"})
    gz=gzip.GzipFile(fileobj=urllib.request.urlopen(req,timeout=600))
    txt=io.TextIOWrapper(gz,encoding="utf-8",errors="replace",newline="")
    header=txt.readline().rstrip("\n").split("\t"); ri=4
    out=[]
    for line in txt:
        if PAT.search(line) is None: continue
        f=line.rstrip("\n").split("\t")
        if len(f)>ri and (set(f[ri].split(",")) & set(SNPS)): out.append(f)
    open(os.path.join(OUT,tag+"_cisld_snps.tsv"),"w",encoding="utf-8",newline="").write("\t".join(header)+"\n"+"\n".join("\t".join(f) for f in out)+"\n")
    print(tag,"matched",len(out),f"{time.time()-t0:.0f}s",flush=True)
print("DONE4")
