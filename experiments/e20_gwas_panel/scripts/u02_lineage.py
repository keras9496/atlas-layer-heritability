"""실험 3u u02: 계통 지도와 분석 계통 (PROTOCOL_exp3u 4.1). 새 적합 없음. u01 fit이 끝난 뒤 돌린다.
 계통 모형 M1_<계통>의 τ* = 계수 × sd(공통 SNP) × M_5_50 / h²(로그). 표준오차는 정렬 블록 delete 값의 jackknife. 단측 P는 τ* > 0.
 분석 계통: S 묶음은 stromal(조건: 단측 P < 0.05), 그 밖은 τ*가 가장 큰 계통(조건: 단측 P < 0.05/12). 조건을 못 채우면 eligible = False.
 출력: results/u02_lineage.csv (형질 × 계통), u02_analysis_lineage.csv (형질), u02_lineage.txt
 실행: python -I -X utf8 u02_lineage.py
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
DATA = Path(DATA); EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"
X = DATA/"exp3"; B = X/"GRCh38"; C = DATA/"exp3c"; N3 = DATA/"exp3n"; U = DATA/"exp3u"; RUN = U/"runs"
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
TR = pd.read_csv(EXP/"traits_exp3u.csv"); TR = TR[TR.include].reset_index(drop=True); TR["discovery_trait"] = TR.discovery_trait.fillna("")
ONLY = sys.argv[1].split(",") if len(sys.argv) > 1 else None   # 시험용: 일부 형질만 (본 실행은 인자 없이)
if ONLY: TR = TR[TR.id.isin(ONLY)].reset_index(drop=True); RES = Path(sys.argv[2]); RES.mkdir(parents=True, exist_ok=True)
SDC = pd.read_csv(C/"annot_sd_common.csv", index_col=0).sd_common
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
strip = lambda n: re.sub(r"L2_\d+$", "", n)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def jse(v): v = np.asarray(v, float); n = len(v); return float(np.sqrt((n - 1) / n * ((v - v.mean()) ** 2).sum()))
def tau1(path, g):
    z = np.load(path); h2 = float(z["h2_log"]); names = [strip(n) for n in z["names"]]; k = names.index(f"L5_{g}"); f = SDC[f"L5_{g}"] * M_TOT / h2
    return float(z["coef"][k] * f), jse(z["part_delete"][:, k] * f), h2, int(z["n_snp"])
rows = []; nfail = 0; nmiss = 0
for r in TR.itertuples():
    for g in LIN:
        p = RUN/f"{r.id}__M1_{g}.npz"
        if not p.exists():
            failed = (RUN/f"{r.id}__M1_{g}.FAILED.txt").exists(); nfail += failed; nmiss += not failed
            rows.append(dict(id=r.id, file=r.file, group=r.group, status=r.status, lineage=g, tau=np.nan, se=np.nan, h2=np.nan, n_snp=np.nan, fit="failed" if failed else "missing")); continue
        t, s, h2, n = tau1(p, g); rows.append(dict(id=r.id, file=r.file, group=r.group, status=r.status, lineage=g, tau=t, se=s, h2=h2, n_snp=n, fit="ok"))
M = pd.DataFrame(rows); M["z"] = M.tau / M.se; M["p1"] = stats.norm.sf(M.z)
ok = M.fit == "ok"; pv = M.loc[ok, "p1"].to_numpy(); o = np.argsort(pv); q = np.empty(len(pv)); q[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]; M.loc[ok, "q"] = np.minimum(q, 1)
M.to_csv(RES/"u02_lineage.csv", index=False)
P(f"=== 실험 3u 계통 지도: 형질 {len(TR)}개 × 계통 12개 ===")
P(f"적합 {int(ok.sum())}개, 실패(빈 블록) {nfail}개, 없음 {nmiss}개")
assert nmiss == 0, "아직 끝나지 않은 적합이 있다"
P(f"h² ≤ 0인 형질: {sorted(set(M[ok & (M.h2 <= 0)].id))}")
P(f"FDR < 0.05인 칸: {int((M.q < .05).sum())} / {int(ok.sum())}")
# 확인(기술): 같은 GWAS(키, Yengo 2022)를 3n에서 따로 정리해 돌린 값과 나란히
if (TR.file == "PASS.Height.Yengo2022").any():
    t3n = tau1(N3/"runs"/"HEIGHT__M1_stromal.npz", "stromal"); hid = TR[TR.file == "PASS.Height.Yengo2022"].id.iloc[0]; t3u = M[(M.id == hid) & (M.lineage == "stromal")].iloc[0]
    P(f"확인(기술): 키(Yengo 2022) stromal τ* — 3n 정리본 {t3n[0]:.2f} ± {t3n[1]:.2f}, 이 묶음의 파일 {t3u.tau:.2f} ± {t3u.se:.2f}")
A = []
for r in TR.itertuples():
    m = M[(M.id == r.id) & ok]
    if len(m) == 0: A.append(dict(id=r.id, file=r.file, trait=r.trait, group=r.group, status=r.status, discovery_trait=r.discovery_trait, eligible=False, analysis_lineage="", note="적합 실패")); continue
    top = m.sort_values("tau", ascending=False).iloc[0]; st = m[m.lineage == "stromal"]; st = st.iloc[0] if len(st) else None
    enr = sorted(m[m.p1 < .05 / 12].lineage)
    if r.group == "S": al = "stromal"; el = bool(st is not None and st.p1 < .05); a = st
    else: al = top.lineage; el = bool(top.p1 < .05 / 12); a = top
    A.append(dict(id=r.id, file=r.file, trait=r.trait, group=r.group, status=r.status, discovery_trait=r.discovery_trait, disease=r.disease, h2=float(m.h2.iloc[0]),
                  top_lineage=top.lineage, top_tau=top.tau, top_se=top.se, top_p1=top.p1, top_eligible=bool(top.p1 < .05 / 12),
                  stromal_tau=(st.tau if st is not None else np.nan), stromal_se=(st.se if st is not None else np.nan), stromal_p1=(st.p1 if st is not None else np.nan),
                  n_enriched=len(enr), enriched=";".join(enr), analysis_lineage=al, analysis_tau=(a.tau if a is not None else np.nan), analysis_p1=(a.p1 if a is not None else np.nan), eligible=el, note=""))
A = pd.DataFrame(A); A.to_csv(RES/"u02_analysis_lineage.csv", index=False)
P("\n=== 분석 계통 (4.1) ===")
P("조건을 채운 형질 수 (묶음 × 새 형질/14형질과 같은 형질)")
P(A.groupby(["group", "status"]).eligible.agg(["sum", "count"]).rename(columns={"sum": "조건 충족", "count": "전체"}).to_string())
P("\n조건을 채운 형질의 분석 계통 분포"); P(A[A.eligible].groupby(["analysis_lineage", "group"]).size().unstack(fill_value=0).to_string())
P("\nS 묶음 (stromal에서; 참고로 τ*가 가장 큰 계통)")
for r in A[A.group == "S"].itertuples(): P(f"  {r.file:42s} stromal τ* {r.stromal_tau:5.2f} ± {r.stromal_se:.2f} (P {r.stromal_p1:.2g}) | 1위 {r.top_lineage} {r.top_tau:.2f} (P {r.top_p1:.2g}) | {'충족' if r.eligible else '미달'}")
P("\n조건을 못 채운 형질 (층위 분석에서 빠짐)")
for r in A[~A.eligible].itertuples(): P(f"  [{r.group}/{r.status}] {r.file:60s} 1위 {r.top_lineage} τ* {r.top_tau:.2f} (P {r.top_p1:.2g})")
P("\nstromal 계통에 농축된(단측 P < 0.05) 새 형질 가운데 S가 아닌 것 (E2의 비교군)")
e2 = A[(A.group != "S") & (A.status == "new") & (A.stromal_p1 < .05)]; P(f"  {len(e2)}개: " + ", ".join(f"{r.file.split('.', 1)[1]}({r.stromal_p1:.2g})" for r in e2.sort_values("stromal_p1").itertuples()))
(RES/"u02_lineage.txt").write_text("\n".join(log), encoding="utf-8")
