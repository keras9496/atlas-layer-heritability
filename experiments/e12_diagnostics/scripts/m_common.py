"""실험 3m 공용: 3j j03_tests.py의 τ*·대비·공동 jackknife 계산을 그대로 옮긴 것(정렬 블록 판 (다)만).
 τ* = 계수 × sd(공통 SNP) × M_5_50 / h²(로그의 반올림 값). delete 값도 같은 배율. 재추출이 여럿이면 블록마다 평균.
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); F = Path(f"{DATA}/exp3f"); G = Path(f"{DATA}/exp3g")
I = Path(f"{DATA}/exp3i"); J = Path(f"{DATA}/exp3j"); RUNS = J/"runs"
LAY = ["CAGE", "DNASE", "H3K27ac"]
SD = pd.concat([pd.read_csv(p/"annot_sd_common.csv", index_col=0).sd_common for p in [Y, F, G, I]]); SD = SD[~SD.index.duplicated()]
B0 = X/"GRCh38"/"baselineLD_v2.2"; M_TOT = sum(float(open(B0/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
MAP = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]; VAL = ["BMD", "NEUT", "TC"]; ALL = MAP + VAL
KO = {"AD": "AD", "MONO": "단핵구 수", "EA": "교육 연수", "SCZ": "조현병", "BIP": "양극성장애", "NEURO": "신경증", "BMI": "BMI", "IBD": "IBD", "RA": "RA", "LDL": "LDL", "HEIGHT": "키",
      "BMD": "골밀도", "NEUT": "호중구 수", "TC": "총콜레스테롤"}
EN = {"AD": "Alzheimer's disease", "MONO": "Monocyte count", "EA": "Educational attainment", "SCZ": "Schizophrenia", "BIP": "Bipolar disorder", "NEURO": "Neuroticism", "BMI": "BMI",
      "IBD": "IBD", "RA": "Rheumatoid arthritis", "LDL": "LDL cholesterol", "HEIGHT": "Height", "BMD": "Bone mineral density", "NEUT": "Neutrophil count", "TC": "Total cholesterol"}
SEEDS = [f"T3_{s:02d}" for s in range(1, 21)]
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def lay_of(n): return [m for m in LAY if n.endswith("_" + m)][0]
def load(t, model, runs=RUNS):
    """층위 → (τ*, delete 값 B개), 정렬 블록 실행(npz)에서"""
    z = np.load(runs/f"{t}__{model}.npz"); h2 = float(z["h2_log"]); names, coef, D = [strip(n) for n in z["names"]], z["coef"], z["part_delete"]
    res = {}
    for k, n in enumerate(names):
        if n in SD.index: f = SD[n] * M_TOT / h2; res[lay_of(n)] = (coef[k] * f, D[:, k] * f)
    assert sorted(res) == sorted(LAY), (t, model, list(res)); return res
def layers(t, models, runs=RUNS):
    """층위별 (τ*, delete 값) — 모형이 여럿이면 블록마다 평균"""
    rs = [load(t, m, runs) for m in models]
    return {m: (float(np.mean([r[m][0] for r in rs])), np.mean([r[m][1] for r in rs], axis=0)) for m in LAY}
def prof(t, models, runs=RUNS):
    L = layers(t, models, runs)
    return (np.array([L["CAGE"][0] - L["DNASE"][0], L["CAGE"][0] - L["H3K27ac"][0]]),
            np.column_stack([L["CAGE"][1] - L["DNASE"][1], L["CAGE"][1] - L["H3K27ac"][1]]))
def jcov(Z): B = Z.shape[0]; z = Z - Z.mean(0); return (B - 1) / B * z.T @ z
def jse(v): return float(np.sqrt(jcov(np.asarray(v).reshape(-1, 1))[0, 0]))
def Qtest(PR, ts):
    """공동 jackknife 다변량 Q (판 다): 형질 K개의 (c1, c2), 공통 평균 GLS, 자유도 2(K−1), F 근사"""
    d = np.concatenate([PR[t][0] for t in ts]); S = jcov(np.hstack([PR[t][1] for t in ts])); K = len(ts)
    Si = np.linalg.inv(S); Xm = np.tile(np.eye(2), (K, 1)); mu = np.linalg.solve(Xm.T @ Si @ Xm, Xm.T @ Si @ d); r = d - Xm @ mu
    Q = float(r @ Si @ r); df = 2 * (K - 1); nu = PR[ts[0]][1].shape[0] - 1; Fs = (nu - df + 1) / (nu * df) * Q
    return Q, df, float(stats.chi2.sf(Q, df)), float(stats.f.sf(Fs, df, nu - df + 1))
def comb_weights(PR, groups, k):
    """groups: [(+1, [형질]), (−1, [형질])] → 형질별 가중치 a (묶음마다 IVW, 단일 형질 jackknife SE)"""
    a = {}
    for sgn, ts in groups:
        w = np.array([1 / jse(PR[t][1][:, k]) ** 2 for t in ts])
        for wi, t in zip(w, ts): a[t] = a.get(t, 0) + sgn * wi / w.sum()
    return a
def comb(PR, groups, k, a=None):
    """가중 결합 대비: (추정치, 공동 jackknife SE, 양측 P, 결합 delete 값). a를 주면 그 가중치를 쓴다."""
    a = a or comb_weights(PR, groups, k)
    est = float(sum(ai * PR[t][0][k] for t, ai in a.items())); dv = sum(ai * PR[t][1][:, k] for t, ai in a.items()); se = jse(dv)
    return est, se, float(2 * stats.norm.sf(abs(est / se))), dv
CON = [("AD vs height (c1)", [(1, ["AD"]), (-1, ["HEIGHT"])], 0), ("Brain traits vs height (c1)", [(1, ["EA", "SCZ", "BIP", "NEURO"]), (-1, ["HEIGHT"])], 0),
       ("IBD, RA vs AD, MONO (c2)", [(1, ["IBD", "RA"]), (-1, ["AD", "MONO"])], 1)]
