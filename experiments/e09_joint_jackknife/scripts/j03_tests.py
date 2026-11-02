"""실험 3j j03: 형질 간 공분산을 반영한 재검정 (PROTOCOL_exp3j_2026-10-14).
 세 판: (가) 원래 값 = 형질별 등분 블록, 형질 독립 (기존 보고값) / (나) 정렬 블록, 형질 독립 / (다) 정렬 블록, 형질 간 공분산 포함 (주).
 다변량 Q: 형질 K개의 (c1 = CAGE − DNase, c2 = CAGE − H3K27ac)를 이은 2K 벡터, 공동 jackknife 공분산(200 블록), 공통 평균 GLS, 자유도 2(K−1).
   (나)는 형질 사이 칸을 0으로 둔 같은 식(= 기존 식). F 근사: F = (ν − p + 1)/(ν p) · Q, p = 2(K−1), ν = 199, 자유도 (p, ν − p + 1).
 대비: 점 추정치는 기존 IVW 가중치(원래 SE)로 고정. (가)·(나) SE = √Σ aᵢ² seᵢ², (다) SE = 가중 결합한 delete 값의 jackknife SE.
 3f는 재추출 20회의 추정치와 delete 값을 블록마다 평균한 뒤 계산한다(기존과 같다).
 사용: python j03_tests.py dry | main
출력: results/j03_<단계>_tests.txt, j03_<단계>_summary.csv
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys, glob
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; R0 = EXP.parent
X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); F = Path(f"{DATA}/exp3f"); G = Path(f"{DATA}/exp3g")
H = Path(f"{DATA}/exp3h"); I = Path(f"{DATA}/exp3i"); J = Path(f"{DATA}/exp3j")
log, SUM = [], []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LAY = ["CAGE", "DNASE", "H3K27ac"]
SD = pd.concat([pd.read_csv(p/"annot_sd_common.csv", index_col=0).sd_common for p in [Y, F, G, I]]); SD = SD[~SD.index.duplicated()]
B0 = X/"GRCh38"/"baselineLD_v2.2"; M_TOT = sum(float(open(B0/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
MAP = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]; VAL = ["BMD", "NEUT", "TC"]
KO = {"AD": "AD", "MONO": "단핵구 수", "EA": "교육 연수", "SCZ": "조현병", "BIP": "양극성장애", "NEURO": "신경증", "BMI": "BMI", "IBD": "IBD", "RA": "RA", "LDL": "LDL", "HEIGHT": "키",
      "BMD": "골밀도", "NEUT": "호중구 수", "TC": "총콜레스테롤"}
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def lay_of(n): return [m for m in LAY if n.endswith("_" + m)][0]
def orig_path(t, model):
    if model.startswith("T3all_"): g = model[6:]; return F/"sldsc"/f"{t}__T3all" if g == OWN[t] else H/"sldsc"/f"{t}__{model}"
    if model.startswith(("T3_", "T3all")): return F/"sldsc"/f"{t}__{model}"
    if model.startswith("G3_"): return G/"sldsc"/f"{t}__{model}"
    if model.startswith("M3_"): return I/"sldsc"/f"{t}__{model}"
def load(t, model, src):
    """층위 → (τ*, delete 값 200개). src: orig = 기존 출력, aligned = 3j 정렬 블록"""
    if src == "orig":
        out = orig_path(t, model); R = pd.read_csv(f"{out}.results", sep="\t")
        h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
        names, coef, D = [strip(n) for n in R.Category], R.Coefficient.to_numpy(), np.loadtxt(f"{out}.part_delete")
    else:
        if model.startswith("T3all_") and model[6:] == OWN.get(t): model = "T3all"
        z = np.load((J/"test" if STAGE == "dry" else J/"runs")/f"{t}__{model}.npz"); h2 = float(z["h2_log"])
        names, coef, D = [strip(n) for n in z["names"]], z["coef"], z["part_delete"]
    res = {}
    for k, n in enumerate(names):
        if n in SD.index: f = SD[n] * M_TOT / h2; res[lay_of(n)] = (coef[k] * f, D[:, k] * f)
    assert sorted(res) == sorted(LAY), (t, model, list(res)); return res
def prof(t, models, src):
    """대비 (c1, c2)와 delete 값 (200 × 2). 모형이 여럿이면(재추출) 블록마다 평균."""
    rs = [load(t, m, src) for m in models]; est = {m: np.mean([r[m][0] for r in rs]) for m in LAY}; dl = {m: np.mean([r[m][1] for r in rs], axis=0) for m in LAY}
    return np.array([est["CAGE"] - est["DNASE"], est["CAGE"] - est["H3K27ac"]]), np.column_stack([dl["CAGE"] - dl["DNASE"], dl["CAGE"] - dl["H3K27ac"]])
def jcov(Z): B = Z.shape[0]; z = Z - Z.mean(0); return (B - 1) / B * z.T @ z
def jse(v): return float(np.sqrt(jcov(v.reshape(-1, 1))[0, 0]))
def Qtest(PR, ts, joint):
    d = np.concatenate([PR[t][0] for t in ts]); S = jcov(np.hstack([PR[t][1] for t in ts])); K = len(ts)
    if not joint:
        S0 = np.zeros_like(S)
        for i in range(K): S0[2 * i:2 * i + 2, 2 * i:2 * i + 2] = S[2 * i:2 * i + 2, 2 * i:2 * i + 2]
        S = S0
    Si = np.linalg.inv(S); Xm = np.tile(np.eye(2), (K, 1)); mu = np.linalg.solve(Xm.T @ Si @ Xm, Xm.T @ Si @ d); r = d - Xm @ mu
    Q = float(r @ Si @ r); df = 2 * (K - 1); nu = PR[ts[0]][1].shape[0] - 1; Fs = (nu - df + 1) / (nu * df) * Q
    return Q, df, stats.chi2.sf(Q, df), stats.f.sf(Fs, df, nu - df + 1)
def Qall(lab, PRo, PRa, ts, ref=None):
    qa, df, pa, _ = Qtest(PRo, ts, False) if PRo else (np.nan, 2 * (len(ts) - 1), np.nan, np.nan); qb, _, pb, _ = Qtest(PRa, ts, False); qc, _, pc, pf = Qtest(PRa, ts, True)
    if ref is not None: assert abs(qa - ref) < .06, (lab, qa, ref)   # 기존 보고값 재현
    P(f"  {lab}: (가) Q {qa:.1f}, p {pa:.2g} | (나) Q {qb:.1f}, p {pb:.2g} | (다) Q {qc:.1f}, p {pc:.2g} (F 근사 p {pf:.2g}) [자유도 {df}]")
    SUM.append(dict(test=lab, kind="Q", est=np.nan, p_a=pa, p_b=pb, p_c=pc, p_c_F=pf, stat_a=qa, stat_b=qb, stat_c=qc)); return pc, pf
def comb(lab, PRo, PRa, groups, k, side, ref=None):
    """groups: [(+1, [형질들]), (−1, [형질들])], 묶음마다 IVW(원래 SE) 평균의 차이. k: 0 = c1, 1 = c2"""
    src = PRo if PRo else PRa; a, keys = [], []
    for sgn, ts in groups:
        w = np.array([1 / jse(src[t][1][:, k]) ** 2 for t in ts]); a += list(sgn * w / w.sum()); keys += ts
    a = np.array(a); est = float(sum(ai * src[t][0][k] for ai, t in zip(a, keys)))
    se_a = np.sqrt(sum(ai ** 2 * jse(PRo[t][1][:, k]) ** 2 for ai, t in zip(a, keys))) if PRo else np.nan
    se_b = np.sqrt(sum(ai ** 2 * jse(PRa[t][1][:, k]) ** 2 for ai, t in zip(a, keys))); se_c = jse(sum(ai * PRa[t][1][:, k] for ai, t in zip(a, keys)))
    pv = (lambda s: stats.norm.sf(est / s)) if side == 1 else (lambda s: 2 * stats.norm.sf(abs(est / s)))
    if ref is not None: assert abs(est - ref) < .006, (lab, est, ref)
    P(f"  {lab}: 차이 {est:+.2f} | (가) SE {se_a:.2f}, p {pv(se_a):.2g} | (나) SE {se_b:.2f}, p {pv(se_b):.2g} | (다) SE {se_c:.2f}, p {pv(se_c):.2g}" + (" [단측]" if side == 1 else " [양측]"))
    SUM.append(dict(test=lab, kind="contrast", est=est, p_a=pv(se_a), p_b=pv(se_b), p_c=pv(se_c), p_c_F=np.nan, stat_a=se_a, stat_b=se_b, stat_c=se_c)); return pv(se_c)
CON = [("AD 대 키 (c1)", [(1, ["AD"]), (-1, ["HEIGHT"])], 0, 2), ("뇌 형질 4개 대 키 (c1)", [(1, ["EA", "SCZ", "BIP", "NEURO"]), (-1, ["HEIGHT"])], 0, 2),
       ("IBD·RA 대 AD·단핵구 수 (c2)", [(1, ["IBD", "RA"]), (-1, ["AD", "MONO"])], 1, 2)]
VCON = ("검증 대비: 호중구 수·총콜레스테롤 대 골밀도 (c1)", [(1, ["NEUT", "TC"]), (-1, ["BMD"])], 0, 1)
if STAGE == "dry":
    PRa = {t: prof(t, ["T3all"], "aligned") for t in ["FAKE1", "FAKE2", "FAKE3"]}
    for t in PRa: OWN[t] = "myeloid"
    P("=== 시험 (가짜 형질 3개, 결과 아님) ==="); Qall("가짜 3개 Q", None, PRa, list(PRa)); comb("FAKE1 대 FAKE2 (c1)", None, PRa, [(1, ["FAKE1"]), (-1, ["FAKE2"])], 0, 2)
    (EXP/"results"/"j03_dry_tests.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()
ALL = MAP + VAL; SEEDS = [f"T3_{s:02d}" for s in range(1, 21)]
# ---------------- 3f (주) ----------------
PRo = {t: prof(t, SEEDS, "orig") for t in ALL}; PRa = {t: prof(t, SEEDS, "aligned") for t in ALL}
P("=== 3f 층위 지도 (재추출 20회 평균) ===")
pc, pf = Qall("J-H1 (주) 다변량 Q, 지도 형질 11개", PRo, PRa, MAP, ref=66.2)
v = "유지" if (pc < .05 and pf < .05) else ("경계" if (pc < .05 or pf < .05) else "유지 안 됨"); P(f"  → J-H1 판정: {v}")
P("J-H2 사전 대비 (양측, α 0.0167):")
ok_old = {}
for lab, gr, k, side in CON:
    comb(lab, PRo, PRa, gr, k, side); r = SUM[-1]; ok_old[lab] = (r["p_a"] < .0167, r["p_c"] < .0167)
P("  결론 변화: " + ", ".join(f"{lab} {'유의' if a else '미달'} → {'유의' if c else '미달'}" for lab, (a, c) in ok_old.items()))
P("J-H3 검증 대비 (단측, α 0.05):")
pv = comb(VCON[0], PRo, PRa, VCON[1], VCON[2], VCON[3]); P(f"  → 검증 대비 {'유의 유지' if pv < .05 else '유의하지 않음'}")
P("\n형질 하나 빼기 (기술):")
for t in MAP: Qall(f"{KO[t]} 제외", PRo, PRa, [x for x in MAP if x != t])
C = np.corrcoef(np.column_stack([PRa[t][1][:, 0] for t in ALL]).T); ii = np.triu_indices(len(ALL), 1); o = np.argsort(-np.abs(C[ii]))
P("\n형질 간 c1 delete 값 상관 (정렬 블록, 절댓값 상위 8쌍; 기술): " + ", ".join(f"{KO[ALL[ii[0][j]]]}–{KO[ALL[ii[1][j]]]} {C[ii][j]:+.2f}" for j in o[:8]))
P(f"  전체 {len(o)}쌍 중앙값 {np.median(C[ii]):+.2f}, 범위 {C[ii].min():+.2f}–{C[ii].max():+.2f}")
# ---------------- 3g ----------------
P("\n=== 3g 역치 (전체 트랙, 단일 실행; 기술) ===")
for q, lab, ref in [("G3_L0p1", "0.1 %", 41.4), ("G3_L0p5", "0.5 %", 33.2), ("G3_L1p0", "1 %", 30.4), ("T3all", "5 %", 61.9)]:
    Po = {t: prof(t, [q], "orig") for t in ALL}; Pa = {t: prof(t, [q], "aligned") for t in ALL}
    Qall(f"[{lab}] Q", Po, Pa, MAP, ref=ref)
    for l2, gr, k, side in CON: comb(f"[{lab}] {l2}", Po, Pa, gr, k, side)
    comb(f"[{lab}] {VCON[0]}", Po, Pa, VCON[1], VCON[2], VCON[3])
# ---------------- 3i ----------------
P("\n=== 3i 세포형 일치 (단일 실행; 기술) ===")
Po = {t: prof(t, [f"M3_{OWN[t]}"], "orig") for t in ALL}; Pa = {t: prof(t, [f"M3_{OWN[t]}"], "aligned") for t in ALL}
Qall("[3i] Q", Po, Pa, MAP, ref=56.8)
for l2, gr, k, side in CON: comb(f"[3i] {l2}", Po, Pa, gr, k, side)
comb(f"[3i] {VCON[0]}", Po, Pa, VCON[1], VCON[2], VCON[3])
# ---------------- 3h ----------------
P("\n=== 3h 공유 계통 비교 (단일 실행; 기술) ===")
Hp = pd.read_csv(R0/"e07_shared_lineages"/"results"/"h02_main_profiles.csv")
DT, CT = ["HEIGHT", "BMD"], ["AD", "MONO", "EA", "SCZ", "BIP", "BMI", "LDL", "TC"]
E = {t: sorted(g.lineage) for t, g in Hp.groupby("trait")}
LG = sorted({g for g in Hp.lineage if any(g in E.get(t, []) for t in DT) and any(g in E.get(t, []) for t in CT)}); assert len(LG) == 6, LG
key = lambda t, g: f"{t}@{g}"
Po, Pa = {}, {}
for t, gs in E.items():
    for g in gs: Po[key(t, g)] = prof(t, [f"T3all_{g}"], "orig"); Pa[key(t, g)] = prof(t, [f"T3all_{g}"], "aligned")
rows, A = [], {}
for g in LG:
    grp = [(1, [key(t, g) for t in CT if g in E.get(t, [])]), (-1, [key(t, g) for t in DT if g in E.get(t, [])])]
    comb(f"[3h] {g}", Po, Pa, grp, 0, 1); r = SUM[-1]; rows.append((g, r["est"], r["stat_a"]))
    a = {}
    for sgn, ks in grp:
        w = np.array([1 / jse(Po[x][1][:, 0]) ** 2 for x in ks])
        for wi, x in zip(w, ks): a[x] = a.get(x, 0) + sgn * wi / w.sum()
    A[g] = a
wl = np.array([1 / s ** 2 for _, _, s in rows]); wl = wl / wl.sum(); Db = float(sum(wi * d for wi, (_, d, _) in zip(wl, rows)))
se_con = float(sum(wi * s for wi, (_, _, s) in zip(wl, rows))); se_ind = float(np.sqrt(sum(wi ** 2 * s ** 2 for wi, (_, _, s) in zip(wl, rows))))
tot = {}
for wi, g in zip(wl, LG):
    for x, ai in A[g].items(): tot[x] = tot.get(x, 0) + wi * ai
se_c = jse(sum(ai * Pa[x][1][:, 0] for x, ai in tot.items())); assert abs(Db - 2.29) < .006, Db
P(f"  [3h] 결합 Δ̄ {Db:+.2f}: (가) 보수적 SE {se_con:.2f}, 단측 p {stats.norm.sf(Db / se_con):.2g} (기존 주 판정) | (가) 독립 SE {se_ind:.2f}, p {stats.norm.sf(Db / se_ind):.2g}"
  f" | (다) 공동 SE {se_c:.2f}, p {stats.norm.sf(Db / se_c):.2g}")
SUM.append(dict(test="[3h] 결합 Δ̄ (가 = 보수적 SE, 나 = 독립 가정 원래 SE, 다 = 공동)", kind="contrast", est=Db, p_a=stats.norm.sf(Db / se_con), p_b=stats.norm.sf(Db / se_ind), p_c=stats.norm.sf(Db / se_c), p_c_F=np.nan,
                stat_a=se_con, stat_b=se_ind, stat_c=se_c))
pd.DataFrame(SUM).to_csv(EXP/"results"/"j03_main_summary.csv", index=False)
(EXP/"results"/"j03_main_tests.txt").write_text("\n".join(log), encoding="utf-8")
