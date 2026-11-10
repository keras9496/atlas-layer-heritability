"""실험 3g u01: 역치 민감도 주석 (PROTOCOL_exp3g_2026-10-11) → 염색체별 LD 점수 (실험 3f v01과 같은 l2 호출).
 칸(계통 g × 층위 m, 5 × 3)마다 전체 트랙 |ALT − REF| 최댓값 (3d L5_{g}_{m}과 같은 트랙 묶음) → 대상 SNP 백분위 → 상위 q.
 q = 0.1 % (L0p1), 0.5 % (L0p5), 1 % (L1p0). 묶음 G3_{q}_{g} = [CAGE, DNASE, H3K27ac] (15개). MAF < 5 % 또는 Atlas 값 없는 SNP = 0 (3d와 같다).
 확인: 같은 방식의 상위 5 % sd가 3d annot_sd_common.csv의 L5_{g}_{m}과 같아야 한다(차이 < 1e-9, 아니면 멈춘다).
 기술: 주석 크기, 같은 계통 안 층위 간 Jaccard, 주석 SNP가 걸친 1 Mb 구간 수.
 시험: U01_TEST=1 → exp3c/test 조각, chr22 만, exp3g/test/annot 에 쓴다 (5 % 확인 생략).
출력: exp3g/annot/G3_{q}_{g}/..., exp3g/annot_sd_common.csv, results/u01_annot.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); C3 = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d"); G = Path(f"{DATA}/exp3g")
PF = X/"GRCh38"/"plink_files"; TEST = os.environ.get("U01_TEST") == "1"; LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py")
SRC = C3/"test" if TEST else C3/"atlas"; AN = G/"test"/"annot" if TEST else G/"annot"; AN.mkdir(parents=True, exist_ok=True); (EXP/"results").mkdir(exist_ok=True)
CHRS = [22] if TEST else range(1, 23)
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LING = ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]; MODS = ["CAGE", "DNASE", "H3K27ac"]
QS = {"L0p1": .001, "L0p5": .005, "L1p0": .01}
TL = pd.read_csv(C3/"track_lineage.csv"); TL["col"] = TL.scorer + "::" + TL.key
CELLS = {(g, m): TL[(TL.lineage == g) & (TL["mod"] == m)].col.tolist() for g in LING for m in MODS}
P("칸별 트랙 수 (전체 트랙, 3d L5와 같음)\n" + pd.DataFrame({m: {g: len(CELLS[(g, m)]) for g in LING} for m in MODS}).to_string())
need = sorted(set(c for v in CELLS.values() for c in v)); parts = sorted(SRC.glob("part_*.parquet")); P(f"트랙 열 {len(need)}개, 조각 {len(parts)}개 읽기")
V = []
for f in parts:
    have = set(pq.read_schema(f).names); cols = [c for c in need if c in have]
    D = pd.read_parquet(f, columns=["atlas_variant"] + cols).set_index("atlas_variant").astype(np.float32)
    V.append(pd.DataFrame({f"{g}_{m}": (D[[c for c in v if c in have]].max(1) if any(c in have for c in v) else np.nan) for (g, m), v in CELLS.items()}, index=D.index))
V = pd.concat(V); V = V[~V.index.duplicated()]
T = pd.read_parquet(X/"targets.parquet")[["rsid", "atlas_variant"]]; T = T.join(V, on="atlas_variant")
n_val = int(T[f"{LING[0]}_{MODS[0]}"].notna().sum()); P(f"대상 {len(T):,}, Atlas 값 있음 {n_val:,}")
R = T[[f"{g}_{m}" for g in LING for m in MODS]].rank(pct=True)
ann = pd.DataFrame({"rsid": T.rsid})
for q, thr in QS.items():
    for g in LING:
        for m in MODS: ann[f"{q}_{g}_{m}"] = (R[f"{g}_{m}"] >= 1 - thr).astype(np.int8)
chk = pd.DataFrame({"rsid": T.rsid})   # 5 % 재현 확인용 (LD 점수 안 만듦)
for g in LING:
    for m in MODS: chk[f"L5_{g}_{m}"] = (R[f"{g}_{m}"] >= .95).astype(np.int8)
ann = ann.drop_duplicates("rsid").set_index("rsid"); chk = chk.drop_duplicates("rsid").set_index("rsid")
NAMES = list(ann.columns)
if not TEST:
    ref = pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common; mine = chk.std(ddof=0); dif = (mine - ref[mine.index]).abs()
    P(f"\n확인: 다시 만든 상위 5 % sd 대 3d L5 sd — 최대 차이 {dif.max():.2e}")
    if dif.max() >= 1e-9: raise SystemExit("상위 5 % 재현 실패 — 트랙 묶음 또는 백분위 방식이 3d와 다르다. 멈춘다.")
P("\n주석 크기 (SNP 수 / 대상 대비 %)")
for q in QS: P(f"  {q}: " + ", ".join(f"{g}_{m} {int(ann[f'{q}_{g}_{m}'].sum()):,} ({ann[f'{q}_{g}_{m}'].mean() * 100:.3f} %)" for g in LING for m in MODS))
P("\n같은 계통 안 층위 간 Jaccard (CAGE–DNase / CAGE–H3K27ac / DNase–H3K27ac); 5 % 참조 포함")
def jac(a, b): return (a & b).sum() / max((a | b).sum(), 1)
for g in LING:
    s = []
    for q, src in [("L0p1", ann), ("L0p5", ann), ("L1p0", ann), ("L5", chk)]:
        s.append(f"{q} " + " / ".join(f"{jac(src[f'{q}_{g}_{a}'], src[f'{q}_{g}_{b}']):.2f}" for a, b in [("CAGE", "DNASE"), ("CAGE", "H3K27ac"), ("DNASE", "H3K27ac")]))
    P(f"  {g:11s} " + "   ".join(s))
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
groups = {f"G3_{q}_{g}": [f"{q}_{g}_{m}" for m in MODS] for q in QS for g in LING}
bins = {n: 0 for n in NAMES}
for c in CHRS:
    out_all = AN/"_all"/f"all.{c}"; (AN/"_all").mkdir(exist_ok=True)
    b = pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, names=["CHR", "SNP", "CM", "BP", "a1", "a2"])
    Ac = ann.reindex(b.SNP).fillna(0).astype(np.int8).reset_index(drop=True)
    for n in NAMES: bins[n] += int(np.unique(b.BP.to_numpy()[Ac[n].to_numpy() == 1] // 1_000_000).size)
    if not Path(f"{out_all}.l2.ldscore.gz").exists():
        Ac.to_csv(f"{out_all}.annot.gz", sep="\t", index=False, compression="gzip")
        r = subprocess.run(["python", "-W", "ignore", str(LDSC), "--l2", "--bfile", str(PF/f"1000G.EUR.hg38.{c}"), "--ld-wind-cm", "1", "--annot", f"{out_all}.annot.gz", "--thin-annot",
                            "--out", str(out_all), "--print-snps", str(X/"hm3_no_MHC.list.txt")], cwd=LDSC.parent, env=env, capture_output=True, text=True)
        if not Path(f"{out_all}.l2.ldscore.gz").exists(): raise SystemExit(f"chr{c} l2 실패\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
    L = pd.read_csv(f"{out_all}.l2.ldscore.gz", sep="\t"); M = open(f"{out_all}.l2.M").read().split(); M5 = open(f"{out_all}.l2.M_5_50").read().split()
    for gname, cols in groups.items():
        d = AN/gname; d.mkdir(exist_ok=True); idx = [NAMES.index(x) for x in cols]
        L[["CHR", "SNP", "BP"] + [f"{x}L2" for x in cols]].to_csv(d/f"{gname}.{c}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        Ac[cols].to_csv(d/f"{gname}.{c}.annot.gz", sep="\t", index=False, compression="gzip")
        open(d/f"{gname}.{c}.l2.M", "w").write("\t".join(M[i] for i in idx)); open(d/f"{gname}.{c}.l2.M_5_50", "w").write("\t".join(M5[i] for i in idx))
    P(f"chr{c} 완료")
P("\n주석 SNP가 걸친 1 Mb 구간 수 (참조 1000G EUR bim 기준)")
for q in QS: P(f"  {q}: " + ", ".join(f"{g}_{m} {bins[f'{q}_{g}_{m}']}" for g in LING for m in MODS))
ann.std(ddof=0).rename("sd_common").to_csv(G/("test_annot_sd_common.csv" if TEST else "annot_sd_common.csv"))
(EXP/"results"/("u01_annot_test.txt" if TEST else "u01_annot.txt")).write_text("\n".join(log), encoding="utf-8")
