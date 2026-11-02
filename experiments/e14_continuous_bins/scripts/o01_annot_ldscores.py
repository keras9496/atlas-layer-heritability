"""실험 3o o01: 연속 백분위 주석과 겹치지 않는 구간 주석 (PROTOCOL_exp3o_2026-10-21, 2026-10-21 고정) → 염색체별 LD 점수.
 점수와 백분위는 3g u01과 같다: 칸(계통 g × 층위 m, 5 × 3)마다 전체 트랙 |ALT − REF| 최댓값(3d L5와 같은 트랙 묶음) → 대상 SNP 백분위 R(rank pct).
 연속 C_{g}_{m} = R (MAF < 5 % 또는 Atlas 값 없음 = 0).
 구간 (겹치지 않음, 원고 역치에 맞춤): B1 R ≥ 0.999 (상위 0–0.1 %), B2 0.995 ≤ R < 0.999, B3 0.99 ≤ R < 0.995, B4 0.95 ≤ R < 0.99.
 확인(어긋나면 멈춘다): 구간 합의 sd가 3d L5(B1–B4), 3g L0p1(B1), L0p5(B1–B2), L1p0(B1–B3)과 같다(차이 < 1e-9).
 묶음: C3_{g} = [C_{g}_CAGE, C_{g}_DNASE, C_{g}_H3K27ac], B12_{g} = [B1–B4 × 3층위]. LD 점수는 3g와 같은 l2 호출(1 cM, --thin-annot, HapMap3 출력).
 기술: 같은 계통 안 층위 간 연속 주석 상관(대상 SNP), 구간별 Jaccard.
출력: exp3o/annot/{C3,B12}_{g}/..., exp3o/annot_sd_common.csv, results/o01_annot.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); C3 = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d")
G = Path(f"{DATA}/exp3g"); O = Path(f"{DATA}/exp3o"); AN = O/"annot"; AN.mkdir(parents=True, exist_ok=True); (EXP/"results").mkdir(exist_ok=True)
PF = X/"GRCh38"/"plink_files"; LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LING = ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]; MODS = ["CAGE", "DNASE", "H3K27ac"]
BINS = {"B1": (1 - .001, 2.0), "B2": (1 - .005, 1 - .001), "B3": (1 - .01, 1 - .005), "B4": (.95, 1 - .01)}   # 3g의 R >= 1 − q, 3d의 R >= .95와 같은 식
TL = pd.read_csv(C3/"track_lineage.csv"); TL["col"] = TL.scorer + "::" + TL.key
CELLS = {(g, m): TL[(TL.lineage == g) & (TL["mod"] == m)].col.tolist() for g in LING for m in MODS}
P("칸별 트랙 수 (전체 트랙, 3d L5와 같음)\n" + pd.DataFrame({m: {g: len(CELLS[(g, m)]) for g in LING} for m in MODS}).to_string())
need = sorted(set(c for v in CELLS.values() for c in v)); parts = sorted((C3/"atlas").glob("part_*.parquet")); P(f"트랙 열 {len(need)}개, 조각 {len(parts)}개 읽기")
V = []
for f in parts:
    have = set(pq.read_schema(f).names); cols = [c for c in need if c in have]
    D = pd.read_parquet(f, columns=["atlas_variant"] + cols).set_index("atlas_variant").astype(np.float32)
    V.append(pd.DataFrame({f"{g}_{m}": (D[[c for c in v if c in have]].max(1) if any(c in have for c in v) else np.nan) for (g, m), v in CELLS.items()}, index=D.index))
V = pd.concat(V); V = V[~V.index.duplicated()]
T = pd.read_parquet(X/"targets.parquet")[["rsid", "atlas_variant"]]; T = T.join(V, on="atlas_variant")
P(f"대상 {len(T):,}, Atlas 값 있음 {int(T[f'{LING[0]}_{MODS[0]}'].notna().sum()):,}")
R = T[[f"{g}_{m}" for g in LING for m in MODS]].rank(pct=True)
ann = pd.DataFrame({"rsid": T.rsid})
for g in LING:
    for m in MODS: ann[f"C_{g}_{m}"] = R[f"{g}_{m}"].fillna(0).astype(np.float32)
for b, (lo, hi) in BINS.items():
    for g in LING:
        for m in MODS: r = R[f"{g}_{m}"]; ann[f"{b}_{g}_{m}"] = ((r >= lo) & (r < hi)).astype(np.int8)
ann = ann.drop_duplicates("rsid").set_index("rsid"); NAMES = list(ann.columns)
# --- 확인: 구간 합 = 기존 주석 ---
r3d = pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common; r3g = pd.read_csv(G/"annot_sd_common.csv", index_col=0).sd_common; worst = 0.0
for g in LING:
    for m in MODS:
        B = {b: ann[f"{b}_{g}_{m}"] for b in BINS}
        for ref, ref_sd, parts_ in [("L5", r3d[f"L5_{g}_{m}"], ["B1", "B2", "B3", "B4"]), ("L0p1", r3g[f"L0p1_{g}_{m}"], ["B1"]), ("L0p5", r3g[f"L0p5_{g}_{m}"], ["B1", "B2"]), ("L1p0", r3g[f"L1p0_{g}_{m}"], ["B1", "B2", "B3"])]:
            u = sum(B[b] for b in parts_); assert u.max() <= 1, "구간이 겹친다"; worst = max(worst, abs(float(u.std(ddof=0)) - float(ref_sd)))
P(f"\n확인: 구간 합의 sd 대 3d L5·3g L0p1·L0p5·L1p0 sd — 최대 차이 {worst:.2e}")
if worst >= 1e-9: raise SystemExit("재현 실패 — 점수 또는 백분위 방식이 3d·3g와 다르다. 멈춘다.")
P("\n구간 크기 (대상 SNP 대비 %; 층위별 같은 값이어야 한다): " + ", ".join(f"{b} {ann[f'{b}_myeloid_CAGE'].mean() * 100:.3f}" for b in BINS))
P("\n같은 계통 안 층위 간 연속 주석 상관 (CAGE–DNase / CAGE–H3K27ac / DNase–H3K27ac)")
for g in LING:
    c = ann[[f"C_{g}_{m}" for m in MODS]].corr().to_numpy(); P(f"  {g:11s} {c[0, 1]:.2f} / {c[0, 2]:.2f} / {c[1, 2]:.2f}")
def jac(a, b): return (a & b).sum() / max((a | b).sum(), 1)
P("\n구간별 층위 간 Jaccard (CAGE–DNase / CAGE–H3K27ac / DNase–H3K27ac)")
for g in LING:
    P(f"  {g:11s} " + "   ".join(f"{b} " + " / ".join(f"{jac(ann[f'{b}_{g}_{x}'], ann[f'{b}_{g}_{y}']):.2f}" for x, y in [("CAGE", "DNASE"), ("CAGE", "H3K27ac"), ("DNASE", "H3K27ac")]) for b in BINS))
ann.std(ddof=0).rename("sd_common").to_csv(O/"annot_sd_common.csv")
# --- LD 점수 ---
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
groups = {**{f"C3_{g}": [f"C_{g}_{m}" for m in MODS] for g in LING}, **{f"B12_{g}": [f"{b}_{g}_{m}" for b in BINS for m in MODS] for g in LING}}
CONT = [n for n in NAMES if n.startswith("C_")]
for c in range(1, 23):
    out_all = AN/"_all"/f"all.{c}"; (AN/"_all").mkdir(exist_ok=True)
    b = pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, names=["CHR", "SNP", "CM", "BP", "a1", "a2"])
    Ac = ann.reindex(b.SNP).fillna(0).reset_index(drop=True)
    for n in NAMES: Ac[n] = Ac[n].astype(np.float32) if n in CONT else Ac[n].astype(np.int8)
    if not Path(f"{out_all}.l2.ldscore.gz").exists():
        Ac.to_csv(f"{out_all}.annot.gz", sep="\t", index=False, compression="gzip", float_format="%.6g")
        r = subprocess.run(["python", "-W", "ignore", str(LDSC), "--l2", "--bfile", str(PF/f"1000G.EUR.hg38.{c}"), "--ld-wind-cm", "1", "--annot", f"{out_all}.annot.gz", "--thin-annot",
                            "--out", str(out_all), "--print-snps", str(X/"hm3_no_MHC.list.txt")], cwd=LDSC.parent, env=env, capture_output=True, text=True)
        if not Path(f"{out_all}.l2.ldscore.gz").exists(): raise SystemExit(f"chr{c} l2 실패\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
    L = pd.read_csv(f"{out_all}.l2.ldscore.gz", sep="\t"); M = open(f"{out_all}.l2.M").read().split(); M5 = open(f"{out_all}.l2.M_5_50").read().split()
    for gname, cols in groups.items():
        d = AN/gname; d.mkdir(exist_ok=True); idx = [NAMES.index(x) for x in cols]
        L[["CHR", "SNP", "BP"] + [f"{x}L2" for x in cols]].to_csv(d/f"{gname}.{c}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        Ac[cols].to_csv(d/f"{gname}.{c}.annot.gz", sep="\t", index=False, compression="gzip", float_format="%.6g")
        open(d/f"{gname}.{c}.l2.M", "w").write("\t".join(M[i] for i in idx)); open(d/f"{gname}.{c}.l2.M_5_50", "w").write("\t".join(M5[i] for i in idx))
    P(f"chr{c} 완료"); (O/"PIPELINE_STATUS.txt").write_text(f"o01 chr{c} 완료", encoding="utf-8")
(EXP/"results"/"o01_annot.txt").write_text("\n".join(log), encoding="utf-8")
