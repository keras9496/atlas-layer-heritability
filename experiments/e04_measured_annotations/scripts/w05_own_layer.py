"""실험 3e w05 (사후 기술, 2026-10-11; 설계안 v3 보강 1): R5를 형질별 고유 층위로 다시 정리한다. 판정 없음.
 고유 층위 = 실험 3f 1단계의 시드 평균 1위 층위(results/v03_main_layers.csv 'top'). 키 = DNASE, 나머지 10개 = CAGE.
 형질 t(계통 g, 고유 층위 m)마다:
   (a) Atlas 단독 τ*: 3d D1 (baselineLD + L5_{g,m})
   (b) eQTL 대비: 같은 D1 모형 안의 GTEx_eQTL_MaxCPP τ*, 차이 SE(part_delete); eQTL 감쇠 = 3d-A B0 → D1
   (c) 같은 층위 실측 대비: 3e S_m (Atlas τ*, 실측 τ*, 차이)
   (d) 실측 세 층위 조건부: CAGE형 = 3e J1. 키 DNASE는 새로 1회 실행한다(BM_stromal + L5_stromal_DNASE, 'J1own_DNASE').
출력: results/w05_own_layer.csv, results/w05_own_layer.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os, re
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); W = Path(f"{DATA}/exp3e")
BV = Path(f"{DATA}/tools/ldsc_bv"); B = X/"GRCh38"; SS = Path(f"{DATA}/exp3b/sumstats")
R3F = EXP.parent/"e05_layer_map"/"results"; env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
OWN = pd.read_csv(R3F/"v03_main_layers.csv").set_index("trait")[["lineage", "top", "cls"]]
SD = pd.concat([pd.read_csv(Y/"baseline_sd_common.csv", index_col=0).sd_common, pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common, pd.read_csv(W/"measured_sd_common.csv", index_col=0).sd_common])
M_TOT = sum(float(open(B/"baselineLD_v2.2"/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
KO = {"AD": "AD", "MONO": "단핵구 수", "EA": "교육 연수", "SCZ": "조현병", "BIP": "양극성장애", "NEURO": "신경증", "BMI": "BMI", "IBD": "IBD", "RA": "RA", "LDL": "LDL", "HEIGHT": "키"}
def parse(out):
    R = pd.read_csv(f"{out}.results", sep="\t"); R["name"] = R.Category.str.replace(r"L2_\d+$", "", regex=True)
    h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
    D = pd.read_csv(f"{out}.part_delete", sep=r"\s+", header=None).to_numpy()
    def ts(n):
        k = int(R.index[R.name == n][0]); f = SD[n] * M_TOT / h2; return R.Coefficient[k] * f, R.Coefficient_std_error[k] * f, D[:, k] * f
    return ts
def jse(d): n = len(d); return np.sqrt((n - 1) / n * ((d - d.mean()) ** 2).sum())
# (d) 키 DNASE의 실측 세 층위 조건부 (새 1회)
out_new = W/"sldsc"/"HEIGHT__J1own_DNASE"
if not Path(f"{out_new}.results").exists():
    refs = [str(B/"baselineLD_v2.2"/"baselineLD."), str(W/"annot"/"BM_stromal"/"BM_stromal."), str(Y/"annot"/"L5_stromal_DNASE"/"L5_stromal_DNASE.")]
    cmd = ["python", "-W", "ignore", "-c", "import sys; sys.argv=['ldsc']+sys.argv[1:]; from ldscore.ldsc import main; main()", "--h2", str(SS/"PASS_Height1.sumstats"), "--ref-ld-chr", ",".join(refs),
           "--w-ld-chr", str(B/"weights"/"weights.hm3_noMHC."), "--overlap-annot", "--frqfile-chr", str(B/"plink_files"/"1000G.EUR.hg38."), "--print-coefficients", "--print-delete-vals", "--out", str(out_new)]
    subprocess.run(cmd, cwd=BV, env=env, capture_output=True, text=True)
    if not Path(f"{out_new}.results").exists(): raise SystemExit(open(f"{out_new}.log", encoding="utf-8").read()[-1500:])
rows = []
for t, r in OWN.iterrows():
    g, m = r.lineage, r.top; a = f"L5_{g}_{m}"
    d1 = parse(Y/"sldsc"/f"{t}__D1__{a}"); b0 = parse(Y/"sldsc"/f"{t}__B0__baseline")
    at, ase, ad = d1(a); et, ese, ed = d1("GTEx_eQTL_MaxCPP"); e0, _, _ = b0("GTEx_eQTL_MaxCPP")
    dd = ad - ed; dse = jse(dd)
    s = parse(W/"sldsc"/f"{t}__S_{m}"); sa, sase, sad = s(a); sm, smse, smd = s(f"M_{g}_{m}"); sdse = jse(sad - smd)
    j = parse(W/"sldsc"/(f"{t}__J1" if m == "CAGE" else f"{t}__J1own_{m}")); ja, jase, _ = j(a)
    rows.append(dict(trait=t, lineage=g, own_layer=m, cls=r.cls, atlas_D1=at, atlas_D1_se=ase, eqtl_D1=et, eqtl_D1_se=ese, atlas_minus_eqtl=at - et, atlas_minus_eqtl_se=dse,
                     atlas_minus_eqtl_p2=2 * stats.norm.sf(abs((at - et) / dse)), eqtl_B0=e0, eqtl_atten=(et - e0) / abs(e0),
                     atlas_S=sa, atlas_S_se=sase, meas_S=sm, meas_S_se=smse, atlas_minus_meas=sa - sm, atlas_minus_meas_se=sdse, atlas_minus_meas_p2=2 * stats.norm.sf(abs((sa - sm) / sdse)),
                     atlas_J1=ja, atlas_J1_se=jase, atlas_J1_p1=stats.norm.sf(ja / jase), retention_J1=ja / at, j1_source="3e J1" if m == "CAGE" else "w05 새 실행(사후)"))
T = pd.DataFrame(rows); T.to_csv(EXP/"results"/"w05_own_layer.csv", index=False)
P("=== R5 형질별 고유 층위 (고유 층위 = 3f 시드 평균 1위) ===")
P("형질        계통        층위   Atlas 단독       eQTL (같은 모형)   Atlas−eQTL (p)    eQTL 감쇠 | 같은 층위 실측 (S): Atlas / 실측 / 차이 (p) | 실측 세 층위 조건부 (p1, 유지율)")
for r in T.itertuples():
    P(f"{KO[r.trait]:9s}  {r.lineage:10s}  {r.own_layer:5s}  {r.atlas_D1:+5.2f} ± {r.atlas_D1_se:.2f}   {r.eqtl_D1:+5.2f} ± {r.eqtl_D1_se:.2f}   {r.atlas_minus_eqtl:+5.2f} ({r.atlas_minus_eqtl_p2:.2g})   {r.eqtl_atten * 100:+4.0f} % | "
      f"{r.atlas_S:+5.2f} / {r.meas_S:+5.2f} / {r.atlas_minus_meas:+5.2f} ({r.atlas_minus_meas_p2:.2g}) | {r.atlas_J1:+5.2f} ± {r.atlas_J1_se:.2f} ({r.atlas_J1_p1:.2g}, {r.retention_J1 * 100:.0f} %)")
P(f"\n요약 (형질 11개, 각자 고유 층위)")
P(f"  Atlas > eQTL: 점추정 {int((T.atlas_minus_eqtl > 0).sum())}/11, 차이 p < 0.05 {int(((T.atlas_minus_eqtl > 0) & (T.atlas_minus_eqtl_p2 < .05)).sum())}/11; eQTL 감쇠 중앙값 {T.eqtl_atten.median() * 100:+.0f} %")
P(f"  Atlas > 같은 층위 실측: 점추정 {int((T.atlas_minus_meas > 0).sum())}/11, 차이 p < 0.05 {int(((T.atlas_minus_meas > 0) & (T.atlas_minus_meas_p2 < .05)).sum())}/11 (실측이 유의하게 큼 {int(((T.atlas_minus_meas < 0) & (T.atlas_minus_meas_p2 < .05)).sum())})")
P(f"  실측 세 층위 조건부 후 Atlas τ* > 0 (단측 p < 0.05): {int((T.atlas_J1_p1 < .05).sum())}/11, 유지율 중앙값 {T.retention_J1.median() * 100:.0f} %")
P("  주의: 키 DNASE의 실측 세 층위 조건부(J1own_DNASE)는 이번에 새로 1회 실행한 사후 분석이다. 나머지는 기존 실행이다.")
(EXP/"results"/"w05_own_layer.txt").write_text("\n".join(log), encoding="utf-8")
