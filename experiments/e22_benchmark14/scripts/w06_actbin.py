"""실험 3w w06: 활성도를 네 구간(상위 0–1, 1–5, 5–10, 10–20%)의 공변량으로 넣은 모형 (PROTOCOL_exp3w 4.3).
 모형 EAb = baselineLD + 효과 주석 셋(L5_<g>_{CAGE, DNASE, H3K27ac}) + 활성도 구간 주석 12개(AB1–AB4 × 세 층위). 정렬 블록.
 보는 값: 형질의 층위(키·BMD는 DNase, 나머지는 CAGE)의 효과 주석 τ* (단측, 14형질의 BH). 상위 5% 하나를 넣은 모형(w04_activity.csv)과 나란히 적는다.
실행: python -I -X utf8 w06_actbin.py      출력: results/w06_actbin.csv, w06_summary.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; DATA = Path(DATA); X = DATA/"exp3"; Y = DATA/"exp3d"; W = DATA/"exp3w"
OWN = {"EA": "CNS_tissue", "N_SCZ": "CNS_tissue", "N_BIP": "CNS_tissue", "NEURO": "CNS_tissue", "N_BMI": "CNS_tissue", "N_AD": "myeloid", "MONO": "myeloid", "NEUT": "myeloid",
       "IBD": "lymphoid", "RA": "lymphoid", "N_LDL": "liver", "TC": "liver", "N_HEIGHT": "stromal", "BMD": "stromal"}
EN = {"EA": "Educational attainment", "N_SCZ": "Schizophrenia", "N_BIP": "Bipolar disorder", "NEURO": "Neuroticism", "N_BMI": "BMI", "N_AD": "Alzheimer's disease", "MONO": "Monocyte count",
      "NEUT": "Neutrophil count", "IBD": "IBD", "RA": "Rheumatoid arthritis", "N_LDL": "LDL cholesterol", "TC": "Total cholesterol", "N_HEIGHT": "Height", "BMD": "Bone mineral density"}
OWNL = {t: ("DNASE" if t in ("N_HEIGHT", "BMD") else "CAGE") for t in OWN}
M_TOT = sum(float(open(X/"GRCh38"/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
SD = pd.concat([pd.read_csv(p/"annot_sd_common.csv", index_col=0).sd_common for p in (Y, W)]); assert not SD.index.duplicated().any()
def jse(v): v = np.asarray(v, float); B = len(v); return float(np.sqrt((B - 1) / B * np.sum((v - v.mean()) ** 2)))
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); m = 1.0
    for r in range(n - 1, -1, -1): m = min(m, p[o[r]] * n / (r + 1)); q[o[r]] = m
    return q
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
A = pd.read_csv(RES/"w04_activity.csv").set_index("trait"); rows = []
for t, g in OWN.items():
    z = np.load(W/"runs"/f"{t}__EAb_{g}.npz"); h2 = float(z["h2_log"]); names = [re.sub(r"L2_\d+$", "", n) for n in z["names"]]; assert len(names) == 97 + 15, len(names)
    d = {n: (float(z["coef"][k] * SD[n] * M_TOT / h2), z["part_delete"][:, k] * SD[n] * M_TOT / h2) for k, n in enumerate(names) if n in SD.index}; assert len(d) == 15
    m = OWNL[t]; e = d[f"L5_{g}_{m}"]; rec = dict(trait=t, group=g, layer=m, effect_only=A.effect_only[t], binary=A.effect[t], binary_se=A.effect_se[t], binary_q=A.q[t], bins=e[0], bins_se=jse(e[1]), retained=e[0] / A.effect_only[t])
    for b in ("AB1", "AB2", "AB3", "AB4"): x = d[f"{b}_{g}_{m}"]; rec[b] = x[0]; rec[b + "_se"] = jse(x[1])
    rows.append(rec)
R = pd.DataFrame(rows); R["p1"] = stats.norm.sf(R.bins / R.bins_se); R["q"] = bh(R.p1); R.to_csv(RES/"w06_actbin.csv", index=False)
P("효과 주석의 τ* (형질의 층위): 효과만 → 활성도 상위 5%와 함께 → 활성도 네 구간과 함께 | 같은 층위 활성도 구간의 τ* (0–1%, 1–5%, 5–10%, 10–20%)")
for r in R.itertuples():
    P(f"  {EN[r.trait]:24s} {r.effect_only:.2f} → {r.binary:.2f} ± {r.binary_se:.2f} (q {r.binary_q:.3f}) → {r.bins:.2f} ± {r.bins_se:.2f} ({r.retained * 100:.0f}%; q {r.q:.3f}{' *' if r.q < .05 else ''}) | {r.AB1:.2f} ± {r.AB1_se:.2f}, {r.AB2:.2f} ± {r.AB2_se:.2f}, {r.AB3:.2f} ± {r.AB3_se:.2f}, {r.AB4:.2f} ± {r.AB4_se:.2f}")
P(f"효과 주석이 남는 형질(FDR < 0.05): 상위 5% 하나 {int((R.binary_q < .05).sum())}/14, 네 구간 {int((R.q < .05).sum())}/14 — " + ", ".join(EN[t] for t in R.trait[R.q < .05]))
(RES/"w06_summary.txt").write_text("\n".join(log), encoding="utf-8")
