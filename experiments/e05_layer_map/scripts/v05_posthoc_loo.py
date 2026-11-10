"""실험 3f v05 (사후 기술, 판정 아님): L-H1 다변량 Q의 형질 하나 빼기. 1단계 결과를 본 뒤 정한 분석이다.
 질문: 형질 간 층위 이질성이 키(유일한 DNase형) 하나에 달려 있는가? 그 밖의 형질끼리도 다른가?
 v03과 같은 시드 평균 jackknife(삭제값 시드 평균)를 쓴다.
출력: results/v05_posthoc_loo.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); F = Path(f"{DATA}/exp3f"); SL = F/"sldsc"
SD = pd.concat([pd.read_csv(F/"annot_sd_common.csv", index_col=0).sd_common, pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common])
B0 = X/"GRCh38"/"baselineLD_v2.2"; M_TOT = sum(float(open(B0/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
TR = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]; LAY = ["CAGE", "DNASE", "H3K27ac"]; SEEDS = [f"{s:02d}" for s in range(1, 21)]
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def load(t, model):
    out = SL/f"{t}__{model}"; R = pd.read_csv(f"{out}.results", sep="\t"); R["name"] = R.Category.str.replace(r"L2_\d+$", "", regex=True)
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    D = pd.read_csv(f"{out}.part_delete", sep=r"\s+", header=None).to_numpy(); res = {}
    for n in R.name:
        if n.startswith("L5r"):
            k = int(R.index[R.name == n][0]); f = SD[n] * M_TOT / h2; res[[m for m in LAY if n.endswith("_" + m)][0]] = (R.Coefficient[k] * f, D[:, k] * f)
    return res
d, S = {}, {}
for t in TR:
    rs = [load(t, f"T3_{s}") for s in SEEDS]; est = {m: np.mean([r[m][0] for r in rs]) for m in LAY}; dl = {m: np.mean([r[m][1] for r in rs], axis=0) for m in LAY}
    d[t] = np.array([est["CAGE"] - est["DNASE"], est["CAGE"] - est["H3K27ac"]]); Z = np.column_stack([dl["CAGE"] - dl["DNASE"], dl["CAGE"] - dl["H3K27ac"]]); B = Z.shape[0]
    Z = Z - Z.mean(0); S[t] = (B - 1) / B * Z.T @ Z
def Q(ts):
    W = {t: np.linalg.inv(S[t]) for t in ts}; dbar = np.linalg.solve(sum(W.values()), sum(W[t] @ d[t] for t in ts))
    q = float(sum((d[t] - dbar) @ W[t] @ (d[t] - dbar) for t in ts)); df = 2 * (len(ts) - 1); return q, df, stats.chi2.sf(q, df)
q, df, p = Q(TR); P(f"전체 11개: Q {q:.1f}, 자유도 {df}, p {p:.2g}")
P("형질 하나 빼기:")
for t in TR: q, df, p = Q([x for x in TR if x != t]); P(f"  {t:7s} 제외: Q {q:5.1f}, 자유도 {df}, p {p:.2g}")
q, df, p = Q([x for x in TR if x != "HEIGHT"]); P(f"\n키를 뺀 10개(모두 CAGE형 또는 CAGE 최빈): Q {q:.1f}, 자유도 {df}, p {p:.2g}")
(EXP/"results"/"v05_posthoc_loo.txt").write_text("\n".join(log), encoding="utf-8")
