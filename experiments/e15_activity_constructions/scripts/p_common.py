"""실험 3p 공용: τ* 계산과 공동 jackknife 검정 (3l l04와 같은 식). p03_tests.py와 p04_fig_export.py가 함께 쓴다.
 τ* = 계수 × sd(공통 SNP) × M_5_50 / h²(로그의 반올림 값), delete 값도 같은 배율. 정렬 블록 200개(3j).
 주석 종류: 효과 L5_·S5_ → "E", 활성도 A5_·Ac_ → "A", 실측 M_ → "M".
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); W = Path(f"{DATA}/exp3e"); J = Path(f"{DATA}/exp3j")
L = Path(f"{DATA}/exp3l"); N3 = Path(f"{DATA}/exp3n"); O = Path(f"{DATA}/exp3p")
LAY = ["CAGE", "DNASE", "H3K27ac"]
SD = pd.concat([pd.read_csv(p, index_col=0).sd_common for p in [Y/"annot_sd_common.csv", W/"measured_sd_common.csv", L/"annot_sd_common.csv", O/"annot_sd_common.csv"]]); SD = SD[~SD.index.duplicated()]
B0 = X/"GRCh38"/"baselineLD_v2.2"; M_TOT = sum(float(open(B0/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
MAP = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]; VAL = ["BMD", "NEUT", "TC"]; ALL = MAP + VAL
NEW = ["HEIGHT", "BMI", "SCZ", "BIP", "LDL"]
OWNL = {t: ("DNASE" if t in ("HEIGHT", "BMD") else "CAGE") for t in ALL}
KO = {"AD": "AD", "MONO": "단핵구 수", "EA": "교육 연수", "SCZ": "조현병", "BIP": "양극성장애", "NEURO": "신경증", "BMI": "BMI", "IBD": "IBD", "RA": "RA", "LDL": "LDL", "HEIGHT": "키",
      "BMD": "골밀도", "NEUT": "호중구 수", "TC": "총콜레스테롤", "FAKE1": "가짜"}
KIND = {"L5_": "E", "S5_": "E", "A5_": "A", "Ac_": "A", "M_": "M"}
MODELS = ["E", "EA", "EAc", "S", "SA", "EAM"]
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def path(t, model, gwas="original", runs=None):
    """모형 실행 파일. 원래 GWAS: E = 3j T3all, EA·EAM = 3l, 새 모형 = 3p. 새 GWAS: E·EA = 3n, 새 모형 = 3p(N_ 접두어)"""
    if runs is not None: return runs/f"{t}__{model}.npz"
    if gwas == "original":
        return {"E": J/"runs"/f"{t}__T3all.npz", "EA": L/"runs"/f"{t}__EA.npz", "EAM": L/"runs"/f"{t}__EAM.npz"}.get(model, O/"runs"/f"{t}__{model}.npz")
    return {"E": N3/"runs"/f"{t}__T3all.npz", "EA": N3/"runs"/f"{t}__EA.npz"}.get(model, O/"runs"/f"N_{t}__{model}.npz")
def load(t, model, gwas="original", runs=None):
    """종류(E/A/M) → 층위 → (τ*, delete 값)"""
    z = np.load(path(t, model, gwas, runs)); h2 = float(z["h2_log"]); res = {}
    for k, n in enumerate(strip(x) for x in z["names"]):
        if n not in SD.index: continue
        kind = [v for p, v in KIND.items() if n.startswith(p)][0]; m = [x for x in LAY if n.endswith("_" + x)][0]; f = SD[n] * M_TOT / h2
        assert m not in res.get(kind, {}), f"{t} {model}: {kind} {m} 중복"
        res.setdefault(kind, {})[m] = (z["coef"][k] * f, z["part_delete"][:, k] * f)
    return res
def jcov(Z): B = Z.shape[0]; z = Z - Z.mean(0); return (B - 1) / B * z.T @ z
def jse(v): return float(np.sqrt(jcov(np.asarray(v).reshape(-1, 1))[0, 0]))
def prof(L3):
    return (np.array([L3["CAGE"][0] - L3["DNASE"][0], L3["CAGE"][0] - L3["H3K27ac"][0]]), np.column_stack([L3["CAGE"][1] - L3["DNASE"][1], L3["CAGE"][1] - L3["H3K27ac"][1]]))
def Qtest(PR, ts):
    d = np.concatenate([PR[t][0] for t in ts]); S = jcov(np.hstack([PR[t][1] for t in ts])); K = len(ts)
    Si = np.linalg.inv(S); Xm = np.tile(np.eye(2), (K, 1)); mu = np.linalg.solve(Xm.T @ Si @ Xm, Xm.T @ Si @ d); r = d - Xm @ mu
    Q = float(r @ Si @ r); df = 2 * (K - 1); nu = PR[ts[0]][1].shape[0] - 1; Fs = (nu - df + 1) / (nu * df) * Q
    return Q, df, float(stats.chi2.sf(Q, df)), float(stats.f.sf(Fs, df, nu - df + 1))
def comb(PR, groups, k, side):
    a = {}
    for sgn, ts in groups:
        w = np.array([1 / jse(PR[t][1][:, k]) ** 2 for t in ts])
        for wi, t in zip(w, ts): a[t] = a.get(t, 0) + sgn * wi / w.sum()
    est = float(sum(ai * PR[t][0][k] for t, ai in a.items())); se = jse(sum(ai * PR[t][1][:, k] for t, ai in a.items()))
    return est, se, float(stats.norm.sf(est / se) if side == 1 else 2 * stats.norm.sf(abs(est / se)))
def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]; return np.minimum(q, 1)
CON = [("AD 대 키 (c1)", [(1, ["AD"]), (-1, ["HEIGHT"])], 0, 2), ("뇌 형질 4개 대 키 (c1)", [(1, ["EA", "SCZ", "BIP", "NEURO"]), (-1, ["HEIGHT"])], 0, 2),
       ("IBD·RA 대 AD·단핵구 수 (c2)", [(1, ["IBD", "RA"]), (-1, ["AD", "MONO"])], 1, 2), ("검증 대비: 호중구 수·총콜레스테롤 대 골밀도 (c1)", [(1, ["NEUT", "TC"]), (-1, ["BMD"])], 0, 1)]
