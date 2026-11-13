"""실험 3u u03: 같은 형질의 독립 GWAS 두 개 사이의 CAGE–DNase 균형 차이. 새 적합 없음. u02가 끝난 뒤 돌린다.
 원고에서의 자리: Results "Skeletal, connective-tissue and other traits enriched in the stromal group share a shift toward chromatin accessibility"
 (쌍 하나의 균형이 얼마나 부정확한지를 적은 문장), Supplementary Methods "Enriched pairs of the panel and the mixed model", Supplementary Table S16의 마지막 두 줄.
 형질 표(traits_exp3u.csv)의 pair 열이 같은 형질의 독립 GWAS를 묶는다. GWAS마다 τ*가 가장 큰 조직 묶음(u02의 top_lineage; 단측 P < 0.05/12일 때만)에서
 전체 트랙 세 층위 결합 모형(u01 fit의 T3all)을 읽고, 두 GWAS가 모두 조건을 채운 쌍의 차이를 낸다. 조건 미달 = 그 묶음의 단측 P ≥ 0.05/12이거나 균형이 정의되지 않음.
 두 GWAS의 묶음이 다른 쌍은 차이를 적되 표시해 둔다(same_lineage = False). 원고가 인용하는 것은 두 GWAS의 묶음이 같은 쌍이다.
 τ* = 계수 × sd(공통 SNP) × M_5_50 / h²(로그). 균형 φD = atan2(DNase τ*, CAGE τ*) (도). delete 값의 각도는 추정값 둘레 ±180° 안에 놓는다.
 차이의 표준오차는 두 GWAS의 정렬 블록 delete 값 차이의 jackknife.
 출력: results/u03_pairs.csv, u03_pairs.txt
 실행: python -I -X utf8 u03_pairs.py
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re
import numpy as np, pandas as pd
from pathlib import Path
DATA = Path(DATA); EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"
X = DATA/"exp3"; B = X/"GRCh38"; Y = DATA/"exp3d"; F = DATA/"exp3f"; RUN = DATA/"exp3u"/"runs"
LAY = ["CAGE", "DNASE", "H3K27ac"]
A = pd.read_csv(RES/"u02_analysis_lineage.csv").set_index("id"); PAIRS = pd.read_csv(EXP/"traits_exp3u.csv")[["id", "pair"]].dropna()
SD = pd.concat([pd.read_csv(p/"annot_sd_common.csv", index_col=0).sd_common for p in [Y, F]]); SD = SD[~SD.index.duplicated()]
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
strip = lambda n: re.sub(r"L2_\d+$", "", n)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def jse(v): v = np.asarray(v, float); n = len(v); return float(np.sqrt((n - 1) / n * ((v - v.mean()) ** 2).sum()))
def load(path):
    """층위 → (τ*, delete 값)"""
    z = np.load(path); h2 = float(z["h2_log"]); names, coef, D = [strip(n) for n in z["names"]], z["coef"], z["part_delete"]; res = {}
    assert h2 > 0, (path, h2)
    for k, n in enumerate(names):
        if n in SD.index: f = SD[n] * M_TOT / h2; res[[m for m in LAY if n.endswith("_" + m)][0]] = (float(coef[k] * f), D[:, k] * f)
    assert sorted(res) == sorted(LAY), (path, list(res)); return res
def balance(Lr):
    """균형 φD(도)와 그 delete 값"""
    e = float(np.degrees(np.arctan2(Lr["DNASE"][0], Lr["CAGE"][0]))); d = np.degrees(np.arctan2(Lr["DNASE"][1], Lr["CAGE"][1])); return e, e + (d - e + 180) % 360 - 180
def name(f): return f.split(".", 1)[1] if "." in f else f
SIDE = ["id", "file", "lineage", "lineage_p", "enriched", "cage", "dnase", "phiD", "phiD_se", "phiD_defined"]
COLS = ["pair"] + [f"{c}_a" for c in SIDE] + [f"{c}_b" for c in SIDE] + ["same_lineage", "diff", "diff_se"]
P("=== 같은 형질의 독립 GWAS 두 개: CAGE–DNase 균형(전체 트랙 세 층위 결합 모형)과 그 차이 ===")
P("GWAS마다 τ*가 가장 큰 조직 묶음에서 읽는다(단측 P < 0.05/12일 때만). 차이의 표준오차는 두 GWAS의 delete 값 차이의 jackknife.")
rows = []
for pn, g in PAIRS[PAIRS.id.isin(A.index)].groupby("pair"):
    assert len(g) == 2, (pn, list(g.id))
    rec = dict(pair=pn); PH = {}; ids = dict(zip("ab", list(g.id)))
    for k, i in ids.items():
        r = A.loc[i]; ok = bool(r.top_eligible); rec.update({f"id_{k}": i, f"file_{k}": r.file, f"lineage_{k}": r.top_lineage, f"lineage_p_{k}": r.top_p1, f"enriched_{k}": ok})
        if not ok: continue
        Lr = load(RUN/f"{i}__T3all_{r.top_lineage}.npz"); ph = balance(Lr); ok = bool(Lr["CAGE"][0] + Lr["DNASE"][0] > 0)   # 균형은 CAGE + DNase > 0일 때만 정의한다
        rec.update({f"cage_{k}": Lr["CAGE"][0], f"dnase_{k}": Lr["DNASE"][0], f"phiD_{k}": ph[0], f"phiD_se_{k}": jse(ph[1]), f"phiD_defined_{k}": ok})
        if ok: PH[k] = ph
    fa, fb = (name(A.loc[ids[k]].file) for k in "ab")
    if len(PH) == 2:
        d = PH["a"][0] - PH["b"][0]; s = jse(PH["a"][1] - PH["b"][1]); same = bool(rec["lineage_a"] == rec["lineage_b"]); rec.update(same_lineage=same, diff=d, diff_se=s)
        P(f"  {pn}: {fa} ({rec['lineage_a']}) {PH['a'][0]:+.1f}° vs {fb} ({rec['lineage_b']}) {PH['b'][0]:+.1f}° — 차이 {d:+.1f}° ± {s:.1f}" + ("" if same else " [두 GWAS의 조직 묶음이 다르다]"))
    else: P(f"  {pn}: {fa} {'있음' if 'a' in PH else '조건 미달'}, {fb} {'있음' if 'b' in PH else '조건 미달'}")
    rows.append(rec)
pd.DataFrame(rows, columns=COLS).to_csv(RES/"u03_pairs.csv", index=False)
(RES/"u03_pairs.txt").write_text("\n".join(log), encoding="utf-8")
