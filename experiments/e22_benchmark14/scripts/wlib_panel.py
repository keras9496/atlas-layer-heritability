"""실험 3w 공용: 106개 GWAS의 농축된 쌍 288개와 쌍마다의 균형(추정값, delete 값 200개)을 읽는다. w05·w07과 같은 정의다.
 방식: all = 전체 트랙(3u T3all), res = 층위당 세 단위 재추출 20회 평균(3u·3w T3r), mat = 세포형 일치(3w M3), dist·prox = 원위·근위(3w PD; 적합이 있을 때만).
 cage/dnase의 τ*와 delete 값도 함께 둔다(대비 계산용).
 쓰는 법: from wlib_panel import *  →  PAIRS, D[(id, 묶음, 방식)] = dict(e, dv, ok, c, d, cdv, ddv), TR, CATN, valid(), gmeans(), ivw(), jse(), Qtest()
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; DATA = Path(DATA); X = DATA/"exp3"; Y = DATA/"exp3d"; F3 = DATA/"exp3f"; I3 = DATA/"exp3i"; K3 = DATA/"exp3k"; U = DATA/"exp3u"; W = DATA/"exp3w"
U3 = EXP.parent/"e20_gwas_panel"
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
FIVE = ["CNS_tissue", "liver", "lymphoid", "myeloid", "stromal"]; NON = [g for g in FIVE if g != "stromal"]; LAY = ["CAGE", "DNASE", "H3K27ac"]
KO = {"CNS_tissue": "뇌 조직", "neural_cells": "신경 세포", "myeloid": "myeloid", "lymphoid": "lymphoid", "liver": "간", "GI": "위장관", "cardiovascular": "심혈관", "skeletal_muscle": "골격근", "adipose": "지방", "stromal": "stromal",
      "reproductive_endocrine": "생식·내분비", "epithelial_barrier": "상피·장벽"}
NAME = {"all": "전체 트랙", "res": "세 단위", "mat": "세포형 일치", "dist": "원위", "prox": "근위"}
M_TOT = sum(float(open(X/"GRCh38"/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
SD = pd.concat([pd.read_csv(p/"annot_sd_common.csv", index_col=0).sd_common for p in (Y, F3, I3, K3)]); assert not SD.index.duplicated().any()
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def jse(v): v = np.asarray(v, float); B = len(v); return float(np.sqrt((B - 1) / B * np.sum((v - v.mean()) ** 2)))
def jcov(Z): Z = np.asarray(Z, float); B = Z.shape[0]; z = Z - Z.mean(0); return (B - 1) / B * z.T @ z
def _load(path):
    z = np.load(path); h2 = float(z["h2_log"]); d = {}
    for k, n in enumerate([strip(n) for n in z["names"]]):
        if n in SD.index: f = SD[n] * M_TOT / h2; d[n] = (float(z["coef"][k] * f), z["part_delete"][:, k] * f)
    return d
def _lay(d, prefix=""):
    out = {}
    for m in LAY:
        ks = [n for n in d if n.endswith("_" + m) and n.startswith(prefix)]; assert len(ks) == 1, (m, prefix, list(d)); out[m] = d[ks[0]]
    return out
def _avg(ds): return {m: (float(np.mean([d[m][0] for d in ds])), np.mean([d[m][1] for d in ds], axis=0)) for m in LAY}
def _phi(Lr):
    c, d = Lr["CAGE"], Lr["DNASE"]; e = float(np.degrees(np.arctan2(d[0], c[0]))); dv = np.degrees(np.arctan2(d[1], c[1])); dv = e + (dv - e + 180) % 360 - 180
    return dict(e=e, dv=dv, ok=bool(c[0] + d[0] > 0), c=c[0], d=d[0], cdv=c[1], ddv=d[1])
def ivw(items):
    w = np.array([1 / jse(d) ** 2 for _, d in items]); w = w / w.sum(); return float(sum(wi * e for wi, (e, _) in zip(w, items))), sum(wi * d for wi, (_, d) in zip(w, items))
TR = pd.read_csv(U3/"traits_exp3u.csv"); TR = TR[TR.include].set_index("id"); LG = pd.read_csv(U3/"results"/"u02_lineage.csv", encoding="utf-8-sig")
E = LG[(LG.p1 < .05 / 12) & (LG.tau > 0)].copy(); assert len(E) == 288
CATN = {"S": "skeletal", "BI": "blood and immune", "LL": "lipid", "BR": "brain and behaviour", "U": "other"}
def cat(i): return CATN[TR.group[i]]
def _p_t3r(i, g, s):
    a = U/"runs"/f"{i}__T3r{s:02d}_{g}.npz"; return a if a.exists() else W/"runs"/f"{i}__T3r{s:02d}_{g}.npz"
PAIRS = [(r.id, r.lineage) for r in E.itertuples()]; D = {}
PD_READY = all((W/"runs"/f"{i}__PD_{g}.npz").exists() for i, g in PAIRS if g in FIVE)          # 근위·원위 적합이 모두 끝났을 때만 읽는다
for i, g in PAIRS:
    D[(i, g, "all")] = _phi(_lay(_load(U/"runs"/f"{i}__T3all_{g}.npz")))
    if g in FIVE:
        D[(i, g, "res")] = _phi(_avg([_lay(_load(_p_t3r(i, g, s))) for s in range(1, 21)])); D[(i, g, "mat")] = _phi(_lay(_load(W/"runs"/f"{i}__M3_{g}.npz")))
        pd_ = W/"runs"/f"{i}__PD_{g}.npz"
        if PD_READY: d = _load(pd_); D[(i, g, "dist")] = _phi(_lay(d, "D_")); D[(i, g, "prox")] = _phi(_lay(d, "P_"))
_ref = pd.read_csv(RES/"w05_pairs.csv"); assert max(abs(D[(r.id, r.group, "all")]["e"] - r.all_phiD) for r in _ref.itertuples()) < 1e-9
def valid(lab, groups): return [(i, g) for i, g in PAIRS if g in groups and (i, g, lab) in D and D[(i, g, lab)]["ok"]]
def gmeans(lab, groups, pairs=None, key=("e", "dv")):
    pairs = valid(lab, groups) if pairs is None else pairs; out = {}
    for g in groups:
        ks = [(i, g) for i, gg in pairs if gg == g]
        if ks: m, dv = ivw([(D[(i, g, lab)][key[0]], D[(i, g, lab)][key[1]]) for i, g in ks]); out[g] = (m, dv, len(ks))
    return out
def Qtest(means, dvs):
    G = len(means); C = np.hstack([-np.ones((G - 1, 1)), np.eye(G - 1)]); d = C @ np.asarray(means); V = jcov(np.asarray(dvs) @ C.T); B = len(dvs)
    Q = float(d @ np.linalg.solve(V, d)); p = G - 1; nu = B - 1; Fs = Q * (nu - p + 1) / (nu * p); return Q, p, float(stats.chi2.sf(Q, p)), float(stats.f.sf(Fs, p, nu - p + 1))
