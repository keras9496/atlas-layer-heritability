"""실험 3e w01: 실측 주석 원자료 (PROTOCOL_exp3e_2026-10-08). 계통 5개 × 모달리티 3개(CAGE, DNASE, H3K27ac)의 Atlas 트랙 용어 → FANTOM5 라이브러리 / ENCODE 기본 피크 파일.
 FANTOM5: 라이브러리 FF:xxxxx-xxxxx 의 is_a 조상을 따라 올라가 가장 가까운 단계의 derives_from UBERON/CL 용어가 Atlas 용어와 같으면 배정한다.
 ENCODE: biosample_ontology.term_id = 용어, 분류 = Atlas biosample_type, GRCh38, released, preferred_default, narrowPeak (DNase는 output_type = peaks).
출력: exp3e/fantom5_library_terms.csv, exp3e/encode_files.csv, exp3e/raw/encode/*.bed.gz, exp3e/term_sources.csv, results/w01_sources.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, time, urllib.parse
from collections import deque
from pathlib import Path
import pandas as pd, requests
EXP = Path(__file__).resolve().parents[1]; C3 = Path(f"{DATA}/exp3c"); W = Path(f"{DATA}/exp3e"); F5 = W/"raw"/"fantom5"; EN = W/"raw"/"encode"
EN.mkdir(parents=True, exist_ok=True); (EXP/"results").mkdir(exist_ok=True)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LING = ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]; MODS = ["CAGE", "DNASE", "H3K27ac"]
CLS = {"tissue": "tissue", "primary_cell": "primary cell", "in_vitro_differentiated_cells": "in vitro differentiated cells"}
TL = pd.read_csv(C3/"track_lineage.csv")
TL = TL[TL.lineage.isin(LING) & TL["mod"].isin(MODS)].copy()
TL["term"] = [k.split("|")[0].split()[1] if m == "CAGE" else k.split()[0] for k, m in zip(TL.key, TL["mod"])]
TERMS = TL[["lineage", "mod", "term", "biosample_name", "biosample_type"]].drop_duplicates(["lineage", "mod", "term", "biosample_type"]).reset_index(drop=True)
P(f"용어 {len(TERMS)}개 (계통 × 모달리티)\n" + TERMS.groupby(["lineage", "mod"]).size().unstack().to_string())

# ---------- FANTOM5 ----------
isa, dfrom, cur = {}, {}, None
for line in open(F5/"ff-phase2-170801.obo.txt", encoding="utf-8"):
    line = line.strip()
    if line.startswith("id: "): cur = line[4:]; isa[cur] = []; dfrom[cur] = []
    elif cur and line.startswith("is_a: "): isa[cur].append(line[6:].split(" !")[0].strip())
    elif cur and re.match(r"(relationship|intersection_of): derives_from ", line):
        t = line.split("derives_from ")[1].split(" !")[0].strip()
        if re.match(r"(UBERON|CL):", t) and t not in dfrom[cur]: dfrom[cur].append(t)
    elif line == "[Term]" or line == "[Typedef]": cur = None if line == "[Typedef]" else cur
def closest_dfrom(ff):
    seen, q = {ff}, deque([(ff, 0)]); best, lvl = [], None
    while q:
        x, d = q.popleft()
        if lvl is not None and d > lvl: break
        if dfrom.get(x): best += [t for t in dfrom[x] if t not in best]; lvl = d
        for y in isa.get(x, []):
            if y.startswith("FF:") and y not in seen: seen.add(y); q.append((y, d + 1))
    return best
import gzip
with gzip.open(F5/"hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz", "rt") as fh:
    for line in fh:
        if line.startswith("00Annotation"): hdr = line.rstrip("\n").split("\t"); break
libs = []
for i, c in enumerate(hdr[1:], start=1):
    m = re.search(r"\.(CNhs\d+)\.(\d+-\d+[A-Z]\d+)\.hg38", c)
    if not m: continue
    d = urllib.parse.unquote(c.split(".CNhs")[0][4:])
    libs.append(dict(col=i, cnhs=m.group(1), ff=f"FF:{m.group(2)}", desc=d, derives_from=";".join(closest_dfrom(f"FF:{m.group(2)}"))))
L = pd.DataFrame(libs); P(f"FANTOM5 라이브러리 {len(L):,}개, derives_from 있음 {int((L.derives_from != '').sum()):,}")
cage_terms = set(TERMS[TERMS["mod"] == "CAGE"].term)
L["atlas_terms"] = [";".join(t for t in s.split(";") if t in cage_terms) for s in L.derives_from]
L.to_csv(W/"fantom5_library_terms.csv", index=False)

# ---------- ENCODE ----------
S = requests.Session(); S.headers["Accept"] = "application/json"
def search(params):
    for k in range(5):
        try:
            r = S.get("https://www.encodeproject.org/search/", params=params, timeout=180)
            return r.json().get("@graph", []) if r.status_code in (200, 404) else []
        except Exception as e: time.sleep(5 * (k + 1)); err = e
    raise SystemExit(f"ENCODE 검색 실패: {err}")
rows = []
for r in TERMS[TERMS["mod"] != "CAGE"].itertuples():
    p = dict(type="File", status="released", assembly="GRCh38", file_format="bed", format="json", limit="all", preferred_default="true")
    p["biosample_ontology.term_id"] = r.term; p["biosample_ontology.classification"] = CLS[r.biosample_type]
    if r.mod == "DNASE": p["assay_title"] = "DNase-seq"; p["output_type"] = "peaks"
    else: p["assay_title"] = "Histone ChIP-seq"; p["target.label"] = "H3K27ac"
    for x in search(p):
        if x.get("file_type") != "bed narrowPeak": continue
        rows.append(dict(lineage=r.lineage, mod=r.mod, term=r.term, biosample_type=r.biosample_type, experiment=x["dataset"].strip("/").split("/")[-1],
                         file=x["accession"], output_type=x.get("output_type"), size=x.get("file_size"), href=x.get("href")))
    time.sleep(0.2)
E = pd.DataFrame(rows); E.to_csv(W/"encode_files.csv", index=False)
P(f"ENCODE 파일 {len(E):,}개 (실험 {E.experiment.nunique():,}개, 고유 파일 {E.file.nunique():,}개, {E.drop_duplicates('file')['size'].sum() / 1e9:.2f} GB)")
for f, h in E.drop_duplicates("file")[["file", "href"]].itertuples(index=False):
    out = EN/f"{f}.bed.gz"
    if out.exists() and out.stat().st_size > 0: continue
    for k in range(5):
        try:
            b = S.get("https://www.encodeproject.org" + h, timeout=300, headers={"Accept": "*/*"}); b.raise_for_status(); out.write_bytes(b.content); break
        except Exception as e:
            if k == 4: raise SystemExit(f"{f} 받기 실패: {e}")
            time.sleep(5 * (k + 1))
P("ENCODE 파일 받기 완료")

# ---------- 용어별 자료 수 ----------
nl = {t: int(L.atlas_terms.str.split(";").apply(lambda s: t in s).sum()) for t in cage_terms}
ne = E.groupby(["lineage", "mod", "term", "biosample_type"]).experiment.nunique()
TERMS["n_sources"] = [nl.get(r.term, 0) if r.mod == "CAGE" else int(ne.get((r.lineage, r.mod, r.term, r.biosample_type), 0)) for r in TERMS.itertuples()]
TERMS.to_csv(W/"term_sources.csv", index=False)
P("\n용어별 자료 수 (CAGE = FANTOM5 라이브러리, DNASE·H3K27ac = ENCODE 실험)")
for (g, m), d in TERMS.groupby(["lineage", "mod"]):
    z = d[d.n_sources == 0]
    P(f"  {g:11s} {m:8s} 용어 {len(d):3d}, 자료 있음 {int((d.n_sources > 0).sum()):3d}, 자료 합계 {int(d.n_sources.sum()):5d}" + (f" | 자료 없음: {', '.join(z.term + ' ' + z.biosample_name)}" if len(z) else ""))
(EXP/"results"/"w01_sources.txt").write_text("\n".join(log), encoding="utf-8")
