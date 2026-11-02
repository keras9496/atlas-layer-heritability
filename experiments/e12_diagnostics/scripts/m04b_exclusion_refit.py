"""실험 3m M4b: 영향력 상위 블록 제외 재적합 (PROTOCOL_exp3m_2026-10-17, 규칙·판정 사전 고정).
 3j j02와 같은 ldsc 호출·패치(정렬 블록)를 쓰고, 제외 블록의 회귀 SNP를 요약통계에서 뺀 뒤 남은 블록(200 − k)으로 jackknife한다.
 제외 집합: m04a의 results/m04_exclusion_sets.json (대비마다 영향력 상위 k = 5(주), 10(부)).
   (i) AD 대 키 c1 → AD·키, (ii) 뇌 형질 4개 대 키 c1 → 교육 연수·조현병·양극성장애·신경증·키, (iii) Q → 지도 형질 11개. k마다 18회, 합계 36회.
 모형: T3all (baselineLD v2.2 + 자기 계통 3층위, 전체 트랙).
 단계: check = 제외 없이 키 T3all 1회 → 3j 정렬 블록 npz와 계수·delete 값이 같아야 한다(상대 차이 < 1e-9) / main = 36회 + 검정
 판정(k = 5): (i)·(ii) 같은 부호이고 양측 P < 0.05 → "소수 블록에 의존하지 않음", 같은 부호·P ≥ 0.05 → "약화", 부호 바뀜 → "소수 블록에 의존"; (iii) Q P < 0.05 → "유지".
출력: exp3m/runs/<형질>__T3all__<집합>.npz, results/m04b_tests.txt, m04b_summary.csv
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, os, re, io, json, time, contextlib
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); B = X/"GRCh38"; Y = Path(f"{DATA}/exp3d"); F = Path(f"{DATA}/exp3f")
J = Path(f"{DATA}/exp3j"); M = Path(f"{DATA}/exp3m"); RUN = M/"runs"; SSD = M/"sumstats"; RUN.mkdir(parents=True, exist_ok=True); SSD.mkdir(exist_ok=True)
BV = Path(f"{DATA}/tools/ldsc_bv"); SS = Path(f"{DATA}/exp3b/sumstats")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
TRAITS = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "MONO": Path(f"{DATA}/exp3b/monocyte_count.sumstats.gz"), "EA": SS/"PASS_Years_of_Education2.sumstats",
          "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats", "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats",
          "IBD": SS/"PASS_IBD.sumstats", "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "LDL": SS/"PASS_LDL.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats"}   # 3j j02와 같다
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal"}
LAY = ["CAGE", "DNASE", "H3K27ac"]
def refs(t): g = OWN[t]; return [str(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.") for m in LAY]
# --- ldsc 불러오기와 패치 (3j j02와 같고, 제외 블록만 다르다) ---
sys.path.insert(0, str(BV))
import ldscore.ldsc as LD
from ldscore import sumstats as S, jackknife as jk, parse as ps
BL = np.load(J/"blocks.npz"); BOUNDS = BL["bounds"]; NB = len(BOUNDS) + 1
KEY = pd.read_parquet(J/"snp_key.parquet").drop_duplicates("SNP").set_index("SNP").key
ST = dict(sep=None, names=None, hsq=None, excl=())
_ld = ps.ldscore; CACHE = {}
def ld_cached(fh, num=None):
    if "baselineLD_v2.2" in fh or "weights.hm3" in fh:
        if (fh, num) not in CACHE: CACHE[(fh, num)] = _ld(fh, num)
        return CACHE[(fh, num)].copy()
    return _ld(fh, num)
ps.ldscore = ld_cached
def separators(snps):
    k = KEY.reindex(snps).to_numpy(); assert not np.isnan(k).any(), "위치 없는 SNP"
    assert (np.diff(k) >= 0).all(), "회귀 SNP가 위치 순서가 아니다"
    c = np.bincount(np.searchsorted(BOUNDS, k.astype(np.int64), side="right"), minlength=NB); ex = np.zeros(NB, bool); ex[list(ST["excl"])] = True
    assert len(c) == NB and (c[ex] == 0).all() and (c[~ex] > 0).all(), "제외 블록에 SNP가 남았거나 빈 블록이 있다"
    return np.concatenate([[0], np.cumsum(c[~ex])]).astype(int)
_rls = S._read_ld_sumstats
def rls(args, lg, fh, alleles=False, dropna=True):
    out = _rls(args, lg, fh, alleles=alleles, dropna=dropna); M_annot, w_ld_cname, ref_ld_cnames, ss, novar = out
    ST["names"] = list(ref_ld_cnames); cm = max(0.001 * ss.N.max(), 80)
    ST["sep"] = separators(ss.SNP.to_numpy()[np.ravel(ss.Z ** 2 < cm)]); return out
S._read_ld_sumstats = rls
_LJF = jk.LstsqJackknifeFast
class LJF(_LJF):
    def __init__(self, x, y, n_blocks=None, separators=None):
        separators = list(ST["sep"]); assert separators[-1] == x.shape[0]; _LJF.__init__(self, x, y, None, separators)
jk.LstsqJackknifeFast = LJF
_eh = S.estimate_h2
def eh(args, lg): r = _eh(args, lg); ST["hsq"] = r; return r
S.estimate_h2 = eh
LOGS = []; _li = LD.Logger.__init__
def li(self, fh): _li(self, fh); LOGS.append(self)
LD.Logger.__init__ = li
def filtered(t, setname, excl):
    f = SSD/f"{t}__{setname}.sumstats.gz"
    if not f.exists():
        d = pd.read_csv(TRAITS[t], sep="\t"); k = KEY.reindex(d.SNP).to_numpy(); b = np.full(len(d), -1); ok = ~np.isnan(k)
        b[ok] = np.searchsorted(BOUNDS, k[ok].astype(np.int64), side="right"); d = d[~np.isin(b, list(excl))]
        tmp = SSD/f"{t}__{setname}.tmp.gz"; d.to_csv(tmp, sep="\t", index=False, compression="gzip"); tmp.replace(f)
    return f
def run(t, setname, excl):
    out = RUN/f"{t}__T3all__{setname}"
    if Path(f"{out}.npz").exists(): return dict(np.load(f"{out}.npz", allow_pickle=False))
    ST.update(sep=None, names=None, hsq=None, excl=tuple(excl)); ss = TRAITS[t] if not excl else filtered(t, setname, excl)
    sys.argv = ["ldsc", "--h2", str(ss), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + refs(t)), "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--out", str(out)]
    with contextlib.redirect_stdout(io.StringIO()): LD.main()
    for l in LOGS: l.log_fh.close()
    LOGS.clear(); h = ST["hsq"]
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    r = dict(names=np.array(ST["names"]), coef=np.asarray(h.coef).ravel(), part_delete=np.asarray(h.part_delete_values), tot=float(np.asarray(h.tot).ravel()[0]),
             h2_log=h2, n_snp=int(ST["sep"][-1]), sep=np.asarray(ST["sep"]))
    assert r["part_delete"].shape == (NB - len(excl), len(r["coef"])), r["part_delete"].shape
    np.savez(f"{out}.npz", **r); return r
rel = lambda a, b: float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))) / max(np.max(np.abs(np.asarray(b, float))), 1e-300))
t0 = time.time()
if STAGE == "check":
    a = run("HEIGHT", "none", ()); z = np.load(J/"runs"/"HEIGHT__T3all.npz")
    r1, r2 = rel(a["coef"], z["coef"]), rel(a["part_delete"], z["part_delete"])
    P(f"check 키 T3all 제외 없음: 계수 {r1:.1e}, delete 값 {r2:.1e} ({(time.time() - t0) / 60:.1f}분)"); assert r1 < 1e-9 and r2 < 1e-9, "재현 실패 — 멈춘다"; P("확인 통과")
    (EXP/"results"/"m04b_check.txt").write_text("\n".join(log), encoding="utf-8"); sys.exit()
SETS = json.loads((EXP/"results"/"m04_exclusion_sets.json").read_text(encoding="utf-8"))
JOBS = [("AD vs height (c1)", "ADH", ["AD", "HEIGHT"]), ("Brain traits vs height (c1)", "BRH", ["EA", "SCZ", "BIP", "NEURO", "HEIGHT"]), ("Q (11 map traits)", "Q11", list(OWN))]
n_all = sum(len(ts) for _, _, ts in JOBS) * 2; i = 0
for lab, code, ts in JOBS:
    for kk in ("k5", "k10"):
        for t in ts:
            run(t, f"{code}_{kk}", SETS[lab][kk]); i += 1
            msg = f"[{i}/{n_all}] {t} {code}_{kk} 완료 ({(time.time() - t0) / 60:.1f}분)"; P("  " + msg); (M/"PIPELINE_STATUS.txt").write_text(msg, encoding="utf-8")
# --- 검정 (m_common의 τ*·대비·Q 계산을 그대로 쓴다) ---
sys.path.insert(0, str(Path(__file__).parent)); import m_common as mc
rows = []
P("\n=== 제외 재적합 결과 (T3all, 전체 트랙) ===")
full = {t: mc.prof(t, ["T3all"]) for t in OWN}
for lab, code, ts in JOBS:
    if code == "Q11": Q0 = mc.Qtest(full, mc.MAP)
    else: g0 = [g for l, g, _ in mc.CON if l == lab][0]; e0, s0, p0, _ = mc.comb(full, g0, 0)
    for kk in ("k5", "k10"):
        PR = {t: mc.prof(t, [f"T3all__{code}_{kk}"], runs=RUN) for t in ts}
        if code == "Q11":
            Q, df, p, pf = mc.Qtest(PR, mc.MAP); v = "유지" if p < .05 else "유지 안 됨"
            P(f"[{lab}] {kk}: Q {Q0[0]:.1f} (P {Q0[2]:.2g}) → {Q:.1f} (자유도 {df}, P {p:.2g}, F 근사 {pf:.2g}) → {v if kk == 'k5' else '(부, 기술)'}")
            rows.append(dict(contrast=lab, k=kk, full=Q0[0], full_p=Q0[2], excl=Q, excl_se=np.nan, excl_p=p, verdict=v if kk == "k5" else ""))
        else:
            e, s, p, _ = mc.comb(PR, g0, 0)
            v = ("소수 블록에 의존하지 않음" if p < .05 else "약화") if np.sign(e) == np.sign(e0) else "소수 블록에 의존"
            P(f"[{lab}] {kk}: {e0:+.2f} ± {s0:.2f} (P {p0:.2g}) → {e:+.2f} ± {s:.2f} (양측 P {p:.2g}; α 0.017 {'충족' if p < .0167 else '미달'}) → {v if kk == 'k5' else '(부, 기술)'}")
            rows.append(dict(contrast=lab, k=kk, full=e0, full_p=p0, excl=e, excl_se=s, excl_p=p, verdict=v if kk == "k5" else ""))
pd.DataFrame(rows).to_csv(EXP/"results"/"m04b_summary.csv", index=False)
(EXP/"results"/"m04b_tests.txt").write_text("\n".join(log), encoding="utf-8")
