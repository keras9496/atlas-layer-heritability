"""실험 3m M4a: 블록 하나 빼기와 블록 영향력 (PROTOCOL_exp3m_2026-10-17, 기술; 새 실행 없음).
 모형: 3f 전체 트랙 T3all(형질마다 자기 계통 3층위), 3j 정렬 블록 200개의 delete 값(= 그 블록을 뺀 추정치).
 대비: (i) AD 대 키 c1, (ii) 뇌 형질 4개 대 키 c1(역분산 가중, 단일 형질 jackknife SE로 가중치), (iii) 지도 형질 11개 다변량 Q(공동 jackknife).
 영향력: (i)·(ii) 전체 추정치 − 블록을 뺀 추정치(양수 = 대비를 지지), (iii) Q − 블록을 뺀 Q(블록을 뺀 Q는 delete 값을 점 추정치로, 공분산은 전체 값 — 근사).
 4b 제외 집합: 대비마다 영향력 상위 k = 5, 10 블록 → results/m04_exclusion_sets.json
 블록 안 연관 신호: 형질별 최대 χ²와 χ² > 29.72(P < 5e-8) SNP 수, AD는 75 좌위 표(1e e13)의 이름.
출력: results/m04a_lobo.csv, m04a_top_blocks.csv, m04_exclusion_sets.json, m04a_tests.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, json
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent)); from m_common import *
EXP = Path(__file__).resolve().parents[1]; R0 = EXP.parent
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
PR = {t: prof(t, ["T3all"]) for t in ALL}; NB = PR["AD"][1].shape[0]; assert NB == 200
BND = np.load(J/"blocks.npz")["bounds"]
lo = np.concatenate([[0], BND]); hi = np.concatenate([BND, [23 * 10**10]])
def span(b): c0, p0, c1, p1 = lo[b] // 10**10, lo[b] % 10**10, (hi[b] - 1) // 10**10, (hi[b] - 1) % 10**10; return f"chr{c0}:{p0 / 1e6:.1f}" + (f"–{p1 / 1e6:.1f} Mb" if c0 == c1 else f"–chr{c1}:{p1 / 1e6:.1f} Mb")
res, sets = [], {}
for lab, groups, k in CON[:2]:
    a = comb_weights(PR, groups, k); est, se, p, dv = comb(PR, groups, k, a); inf = est - dv
    res.append(pd.DataFrame(dict(contrast=lab, block=np.arange(NB), lobo=dv, influence=inf)))
    o = np.argsort(-inf); sets[lab] = {f"k{kk}": [int(b) for b in o[:kk]] for kk in (5, 10)}
    P(f"[{lab}] 전체 {est:+.2f} ± {se:.2f} (P {p:.2g}) | 블록 하나 빼기 범위 {dv.min():+.2f} ~ {dv.max():+.2f}, 부호 바뀜 {int((np.sign(dv) != np.sign(est)).sum())}/200"
      f" | 영향력 최대 {inf.max():+.2f} (= SE의 {inf.max() / se:.2f}배), 상위 5개 합 {inf[o[:5]].sum():+.2f}")
Q, df, pq, pf = Qtest(PR, MAP); S = jcov(np.hstack([PR[t][1] for t in MAP])); Si = np.linalg.inv(S); Xm = np.tile(np.eye(2), (len(MAP), 1))
def Qd(d): mu = np.linalg.solve(Xm.T @ Si @ Xm, Xm.T @ Si @ d); r = d - Xm @ mu; return float(r @ Si @ r)
Qb = np.array([Qd(np.concatenate([PR[t][1][b] for t in MAP])) for b in range(NB)]); infQ = Q - Qb; oQ = np.argsort(-infQ)
sets["Q (11 map traits)"] = {f"k{kk}": [int(b) for b in oQ[:kk]] for kk in (5, 10)}
res.append(pd.DataFrame(dict(contrast="Q (11 map traits)", block=np.arange(NB), lobo=Qb, influence=infQ)))
P(f"[Q] 전체 {Q:.1f} (P {pq:.2g}) | 블록 하나 빼기 Q 범위 {Qb.min():.1f} ~ {Qb.max():.1f} | 감소 최대 {infQ.max():.1f}, 상위 5개 합 {infQ[oQ[:5]].sum():.1f} (근사)")
L = pd.concat(res); L.to_csv(EXP/"results"/"m04a_lobo.csv", index=False)
# 영향력 상위 블록의 위치와 연관 신호
KEY = pd.read_parquet(J/"snp_key.parquet").drop_duplicates("SNP").set_index("SNP").key
from m_common import X as X3
SS = {"AD": X3/"bellenguez_noAPOE.sumstats.gz", "HEIGHT": Path(f"{DATA}/exp3b/sumstats/PASS_Height1.sumstats"), "EA": Path(f"{DATA}/exp3b/sumstats/PASS_Years_of_Education2.sumstats"),
      "SCZ": Path(f"{DATA}/exp3b/sumstats/PASS_Schizophrenia.sumstats"), "BIP": Path(f"{DATA}/exp3b/sumstats/PASS_Bipolar_Disorder.sumstats"),
      "NEURO": Path(f"{DATA}/exp3b/sumstats/PASS_Neuroticism.sumstats")}
CH = {}
for t, f in SS.items():
    d = pd.read_csv(f, sep="\t", usecols=["SNP", "Z"]); d["key"] = KEY.reindex(d.SNP).to_numpy(); d = d.dropna(subset=["key", "Z"])
    d["b"] = np.searchsorted(BND, d.key.astype(np.int64).to_numpy(), side="right"); d["c2"] = d.Z ** 2; CH[t] = d.groupby("b").c2.agg(["max", lambda x: int((x > 29.72).sum())])
LOC = pd.read_csv(Path(__file__).resolve().parents[1]/"inputs"/"loci75_table_S3.csv", encoding="utf-8-sig")
LOC["b"] = np.searchsorted(BND, LOC.chr.astype(np.int64).to_numpy() * 10**10 + LOC.pos.astype(np.int64).to_numpy(), side="right")
top = []
for lab, traits in [(CON[0][0], ["AD", "HEIGHT"]), (CON[1][0], ["EA", "SCZ", "BIP", "NEURO", "HEIGHT"]), ("Q (11 map traits)", ["AD", "HEIGHT"])]:
    s = L[L.contrast == lab].sort_values("influence", ascending=False).head(10)
    for rk, r in enumerate(s.itertuples(), 1):
        b = int(r.block); sig = {t: (float(CH[t].loc[b, "max"]) if b in CH[t].index else np.nan, int(CH[t].loc[b].iloc[1]) if b in CH[t].index else 0) for t in traits}
        loci = ", ".join(f"{g}" for g in LOC[LOC.b == b].nearest_gene_TSS.astype(str))
        top.append(dict(contrast=lab, rank=rk, block=b, region=span(b), influence=r.influence, lobo=r.lobo, AD_loci=loci, **{f"maxchi2_{t}": v[0] for t, v in sig.items()}, **{f"ngws_{t}": v[1] for t, v in sig.items()}))
        P(f"  [{lab}] {rk:2d}. 블록 {b:3d} {span(b):28s} 영향력 {r.influence:+.2f} | " + ", ".join(f"{KO[t]} 최대 χ² {v[0]:.0f} (GWS {v[1]})" for t, v in sig.items()) + (f" | AD 좌위 {loci}" if loci else ""))
pd.DataFrame(top).to_csv(EXP/"results"/"m04a_top_blocks.csv", index=False)
(EXP/"results"/"m04_exclusion_sets.json").write_text(json.dumps(sets, ensure_ascii=False, indent=1), encoding="utf-8")
(EXP/"results"/"m04a_tests.txt").write_text("\n".join(log), encoding="utf-8")
