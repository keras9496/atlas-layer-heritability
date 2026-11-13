"""실험 3t t01: 유전율 비율·농축 배수와 τ* 크기 비교 (PROTOCOL_exp3t_2026-11-10 T1).
 단계 (python t01_enrichment.py check | fit | extract | scale):
   check   = 같은 실행기로 원래 키 GWAS의 T3all을 다시 돌려 exp3f/sldsc/HEIGHT__T3all.results와 계수가 같아야 한다(상대 차이 < 1e-9).
   fit     = 더 큰 GWAS 6개 × (계통 주석 M1, T3all) = 12회. 표준 LDSC(기본 블록, --overlap-annot). 재실행 안전.
   extract = 기존 .results와 새 .results에서 Prop._SNPs, Prop._h2, Enrichment, τ*를 모은다 → results/t01_enrichment.csv
   scale   = T3all 모형에서 최상위 층위 주석의 τ*가 기준선 97개 주석과 비교해 몇 번째인지 → results/t01_scale.csv
 읽기만: exp3, exp3b, exp3c, exp3d, exp3e, exp3f, exp3n, exp3r.  쓰기: exp3t/sldsc, results/t01_*.
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, os, re, subprocess, time
import numpy as np, pandas as pd
from pathlib import Path
STAGE = sys.argv[1]
DATA = Path(DATA); EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; RES.mkdir(exist_ok=True)
X = DATA/"exp3"; B = X/"GRCh38"; C = DATA/"exp3c"; Y = DATA/"exp3d"; E = DATA/"exp3e"; F = DATA/"exp3f"; T = DATA/"exp3t"
SL = T/"sldsc"; SL.mkdir(parents=True, exist_ok=True); BV = DATA/"tools"/"ldsc_bv"; SS = DATA/"exp3b"/"sumstats"
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
MAIN = {"AD": X/"bellenguez_noAPOE.sumstats.gz", "MONO": DATA/"exp3b"/"monocyte_count.sumstats.gz", "EA": SS/"PASS_Years_of_Education2.sumstats",
        "SCZ": SS/"PASS_Schizophrenia.sumstats", "BIP": SS/"PASS_Bipolar_Disorder.sumstats", "NEURO": SS/"PASS_Neuroticism.sumstats", "BMI": SS/"PASS_BMI1.sumstats",
        "IBD": SS/"PASS_IBD.sumstats", "RA": SS/"PASS_Rheumatoid_Arthritis.sumstats", "LDL": SS/"PASS_LDL.sumstats", "HEIGHT": SS/"PASS_Height1.sumstats",
        "BMD": F/"sumstats"/"BMD.sumstats.gz", "NEUT": F/"sumstats"/"NEUT.sumstats.gz", "TC": F/"sumstats"/"TC.sumstats.gz"}
LARGER = {"HEIGHT": DATA/"exp3n"/"sumstats"/"HEIGHT.sumstats.gz", "BMI": DATA/"exp3n"/"sumstats"/"BMI.sumstats.gz", "SCZ": DATA/"exp3n"/"sumstats"/"SCZ.sumstats.gz",
          "BIP": DATA/"exp3n"/"sumstats"/"BIP.sumstats.gz", "LDL": DATA/"exp3n"/"sumstats"/"LDL.sumstats.gz", "AD": DATA/"exp3r"/"sumstats"/"N_AD.sumstats.gz"}
OWN = {"AD": "myeloid", "MONO": "myeloid", "EA": "CNS_tissue", "SCZ": "CNS_tissue", "BIP": "CNS_tissue", "NEURO": "CNS_tissue", "BMI": "CNS_tissue",
       "IBD": "lymphoid", "RA": "lymphoid", "LDL": "liver", "HEIGHT": "stromal", "BMD": "stromal", "NEUT": "myeloid", "TC": "liver"}
MAP = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT"]; VAL = ["BMD", "NEUT", "TC"]
TOP = {t: ("DNASE" if t in ("HEIGHT", "BMD") else "CAGE") for t in OWN}
LAY = ["CAGE", "DNASE", "H3K27ac"]
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
strip = lambda n: re.sub(r"L2_\d+$", "", n)
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
def refs(g, model):
    if model == "T3all": return [str(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.") for m in LAY]
    if model == "M1": return [str(C/"annot"/f"L5_{g}"/f"L5_{g}.")]
    raise ValueError(model)
def fit(ss, out, g, model):
    if not Path(f"{out}.results").exists():
        code = f"import sys; sys.path.insert(0, r'{BV}'); sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()"
        cmd = ["python", "-I", "-X", "utf8", "-W", "ignore", "-c", code, "--h2", str(ss), "--ref-ld-chr", ",".join([str(B/"baselineLD_v2.2"/"baselineLD.")] + refs(g, model)),
               "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--overlap-annot", "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."),
               "--print-coefficients", "--print-delete-vals", "--out", str(out)]
        subprocess.run(cmd, cwd=str(T), env=env, capture_output=True, text=True)
        if not Path(f"{out}.results").exists(): raise SystemExit(f"{out} 실패: {open(str(out) + '.log', encoding='utf-8').read()[-1500:]}")
    return out
def h2_of(out): return float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
t0 = time.time()
if STAGE == "check":
    out = fit(MAIN["HEIGHT"], SL/"CHECK_HEIGHT__T3all", "stromal", "T3all")
    a = pd.read_csv(f"{out}.results", sep="\t"); b = pd.read_csv(F/"sldsc"/"HEIGHT__T3all.results", sep="\t")
    r = float(np.max(np.abs(a.Coefficient - b.Coefficient)) / np.max(np.abs(b.Coefficient))); r2 = float(np.max(np.abs(a["Prop._h2"] - b["Prop._h2"])))
    P(f"확인: HEIGHT T3all 다시 적합 — 계수 상대 차이 {r:.1e}, Prop._h2 최대 차이 {r2:.1e} ({(time.time() - t0) / 60:.1f}분)"); assert r < 1e-9 and list(a.Category) == list(b.Category), "재현 실패 — 멈춘다"
    P("확인 통과"); (RES/"t01_check.txt").write_text("\n".join(log), encoding="utf-8")
elif STAGE == "fit":
    jobs = [(t, m) for t in LARGER for m in ("M1", "T3all")]
    for i, (t, m) in enumerate(jobs):
        fit(LARGER[t], SL/f"L_{t}__{m}", OWN[t], m); msg = f"[{i + 1}/{len(jobs)}] 더 큰 GWAS {t} {m} 완료 ({(time.time() - t0) / 60:.1f}분)"; P("  " + msg)
        (T/"PIPELINE_STATUS_t01.txt").write_text(msg, encoding="utf-8")
    (RES/"t01_fit.txt").write_text("\n".join(log), encoding="utf-8")
elif STAGE == "extract":
    SD = pd.concat([pd.read_csv(p, index_col=0).sd_common for p in [C/"annot_sd_common.csv", Y/"annot_sd_common.csv", E/"measured_sd_common.csv", Y/"baseline_sd_common.csv"]])
    SD = SD[~SD.index.duplicated()]
    def take(out, gwas, t, model, keep):
        R = pd.read_csv(f"{out}.results", sep="\t").rename(columns={"Prop._SNPs": "ps", "Prop._h2": "ph", "Prop._h2_std_error": "phse"}); R["name"] = R.Category.map(strip); h2 = h2_of(out); rows = []
        for r in R[R.name.isin(keep)].itertuples():
            f = SD[r.name] * M_TOT / h2
            rows.append(dict(trait=t, gwas=gwas, role="map" if t in MAP else "validation", lineage=OWN[t], model=model, annot=r.name, h2_obs=h2, prop_snps=r.ps, prop_h2=r.ph, prop_h2_se=r.phse,
                             enrichment=r.Enrichment, enrichment_se=r.Enrichment_std_error, enrichment_p=r.Enrichment_p, tau_star=r.Coefficient * f, tau_star_se=r.Coefficient_std_error * f))
        assert len(rows) == len(keep), (out, keep, [x["annot"] for x in rows]); return rows
    rows = []
    for t in MAP + VAL:
        g = OWN[t]; m1 = C/"sldsc"/f"{t}__M1__{g}" if t in MAP else F/"sldsc"/f"{t}__M1_{g}"
        rows += take(m1, "main", t, "lineage", [f"L5_{g}"])
        rows += take(F/"sldsc"/f"{t}__T3all", "main", t, "three layers, joint", [f"L5_{g}_{m}" for m in LAY] + ["GTEx_eQTL_MaxCPP"])
        if t in MAP:
            for m in LAY: rows += take(E/"sldsc"/f"{t}__M1_{m}", "main", t, "measured, single", [f"M_{g}_{m}"])
    for t in LARGER:
        g = OWN[t]; rows += take(SL/f"L_{t}__M1", "larger", t, "lineage", [f"L5_{g}"])
        rows += take(SL/f"L_{t}__T3all", "larger", t, "three layers, joint", [f"L5_{g}_{m}" for m in LAY] + ["GTEx_eQTL_MaxCPP"])
    D = pd.DataFrame(rows); D.to_csv(RES/"t01_enrichment.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
    P("=== T1 유전율 비율과 농축 (Prop._h2 ± s.e.; Enrichment ± s.e.) ===")
    for gw in ("main", "larger"):
        P(f"\n[{gw} GWAS] 계통 주석(공통 SNP의 약 5%)")
        for r in D[(D.gwas == gw) & (D.model == "lineage")].itertuples():
            P(f"  {r.trait:7s} {r.lineage:11s} SNP {100 * r.prop_snps:.1f}% → h2 {100 * r.prop_h2:.0f}% (± {100 * r.prop_h2_se:.0f}), {r.enrichment:.1f}배 (± {r.enrichment_se:.1f}), P = {r.enrichment_p:.1g}; τ* {r.tau_star:.2f}")
        P(f"[{gw} GWAS] 세 층위 결합 모형의 층위별 주석")
        for t in D[D.gwas == gw].trait.unique():
            d = D[(D.gwas == gw) & (D.trait == t) & (D.model == "three layers, joint")]
            P(f"  {t:7s} " + "; ".join(f"{r.annot.split('_')[-1]} h2 {100 * r.prop_h2:.0f}% ({r.enrichment:.1f}배 ± {r.enrichment_se:.1f})" for r in d.itertuples()))
    d = D[(D.gwas == "main") & (D.model == "lineage")]
    P(f"\n요약(주 GWAS, 계통 주석, 14형질): 유전율 비율 {100 * d.prop_h2.min():.0f}–{100 * d.prop_h2.max():.0f}% (중앙 {100 * d.prop_h2.median():.0f}%), 농축 {d.enrichment.min():.1f}–{d.enrichment.max():.1f}배 (중앙 {d.enrichment.median():.1f})")
    d = D[(D.gwas == "main") & (D.model == "lineage") & (D.role == "map")]
    P(f"요약(지도 형질 11개): 유전율 비율 {100 * d.prop_h2.min():.0f}–{100 * d.prop_h2.max():.0f}% (중앙 {100 * d.prop_h2.median():.0f}%), 농축 {d.enrichment.min():.1f}–{d.enrichment.max():.1f}배 (중앙 {d.enrichment.median():.1f})")
    P("\n[실측 주석, 단일 모형, 지도 형질] 층위별 유전율 비율(농축)")
    for t in MAP:
        d = D[(D.trait == t) & (D.model == "measured, single")]
        P(f"  {t:7s} " + "; ".join(f"{r.annot.split('_')[-1]} SNP {100 * r.prop_snps:.1f}% → h2 {100 * r.prop_h2:.0f}% ({r.enrichment:.1f}배)" for r in d.itertuples()))
    (RES/"t01_extract.txt").write_text("\n".join(log), encoding="utf-8")
elif STAGE == "scale":
    # 기준선 97개 주석의 표준편차 (MAF ≥ 5% SNP)
    sdf = T/"baseline_sd_common_all.csv"
    if not sdf.exists():
        n = 0; s1 = None; s2 = None
        for c in range(1, 23):
            A = pd.read_csv(B/"baselineLD_v2.2"/f"baselineLD.{c}.annot.gz", sep="\t"); fr = pd.read_csv(B/"plink_files"/f"1000G.EUR.hg38.{c}.frq", sep=r"\s+")
            assert len(A) == len(fr) and (A.SNP.to_numpy() == fr.SNP.to_numpy()).all(), f"chr{c}: 주석과 frq의 SNP 순서가 다르다"
            cols = [x for x in A.columns if x not in ("CHR", "BP", "SNP", "CM")]; k = (fr.MAF.to_numpy() >= 0.05); V = A[cols].to_numpy(float)[k]
            if s1 is None: s1 = np.zeros(V.shape[1]); s2 = np.zeros(V.shape[1])
            n += V.shape[0]; s1 += V.sum(0); s2 += (V ** 2).sum(0); P(f"  chr{c}: 공통 SNP {int(k.sum()):,} ({(time.time() - t0) / 60:.1f}분)")
        sd = np.sqrt(s2 / n - (s1 / n) ** 2); pd.Series(sd, index=cols, name="sd_common").to_csv(sdf); P(f"공통 SNP 합계 {n:,} (M_5_50 합 {M_TOT:,.0f})")
    SDB = pd.read_csv(sdf, index_col=0).sd_common; old = pd.read_csv(Y/"baseline_sd_common.csv", index_col=0).sd_common
    r = float(abs(SDB["GTEx_eQTL_MaxCPP"] - old["GTEx_eQTL_MaxCPP"]) / old["GTEx_eQTL_MaxCPP"]); P(f"확인: GTEx_eQTL_MaxCPP 표준편차 {SDB['GTEx_eQTL_MaxCPP']:.8f} (기존 {old['GTEx_eQTL_MaxCPP']:.8f}), 상대 차이 {r:.1e}")
    assert r < 1e-6, "기준선 표준편차 재현 실패 — 멈춘다"
    SDA = pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common; SD = pd.concat([SDB, SDA]); rows = []
    jobs = [("main", t, F/"sldsc"/f"{t}__T3all") for t in MAP + VAL] + [("larger", t, SL/f"L_{t}__T3all") for t in LARGER]
    for gw, t, out in jobs:
        R = pd.read_csv(f"{out}.results", sep="\t"); R["name"] = R.Category.map(strip); h2 = h2_of(out); g = OWN[t]
        R = R[R.name != "base"].copy(); assert R.name.isin(SD.index).all(), list(R.name[~R.name.isin(SD.index)])
        R["tau_star"] = R.Coefficient * R.name.map(SD) * M_TOT / h2; R["z"] = R["Coefficient_z-score"]
        top = f"L5_{g}_{TOP[t]}"; me = R[R.name == top].iloc[0]; base = R[~R.name.str.startswith("L5_")].sort_values("tau_star", ascending=False)
        rk = int((R.tau_star > me.tau_star).sum()) + 1; rk_abs = int((R.tau_star.abs() > abs(me.tau_star)).sum()) + 1; rkz = int((R.z > me.z).sum()) + 1
        rows.append(dict(gwas=gw, trait=t, top_annot=top, tau_star=me.tau_star, z=me.z, n_annot=len(R), rank_tau_star=rk, rank_abs_tau_star=rk_abs, rank_z=rkz,
                         best_baseline=base.name.iloc[0], best_baseline_tau_star=base.tau_star.iloc[0], second_baseline=base.name.iloc[1], second_baseline_tau_star=base.tau_star.iloc[1],
                         gtex_maxcpp_tau_star=float(R[R.name == "GTEx_eQTL_MaxCPP"].tau_star.iloc[0])))
    D = pd.DataFrame(rows); D.to_csv(RES/"t01_scale.csv", index=False)
    P("\n=== T1 크기 비교: 세 층위 결합 모형(T3all)에서 최상위 층위 주석의 τ*와 기준선 주석(기준 'base' 제외 96개 + 예측 3개) ===")
    for r in D.itertuples():
        P(f"  [{r.gwas}] {r.trait:7s} {r.top_annot}: τ* {r.tau_star:.2f} (z {r.z:.1f}) — τ* 순위 {r.rank_tau_star}/{r.n_annot}, |τ*| 순위 {r.rank_abs_tau_star}, z 순위 {r.rank_z}"
          f" | 기준선 최대 {r.best_baseline} {r.best_baseline_tau_star:.2f}, 둘째 {r.second_baseline} {r.second_baseline_tau_star:.2f} | GTEx MaxCPP {r.gtex_maxcpp_tau_star:.2f}")
    d = D[D.gwas == "main"]; P(f"주 GWAS 14형질: τ* 순위 1위 {int((d.rank_tau_star == 1).sum())}개, z 순위 1위 {int((d.rank_z == 1).sum())}개; 지도 형질 11개 중 τ* 1위 {int((d[d.trait.isin(MAP)].rank_tau_star == 1).sum())}개")
    (RES/"t01_scale.txt").write_text("\n".join(log), encoding="utf-8")
