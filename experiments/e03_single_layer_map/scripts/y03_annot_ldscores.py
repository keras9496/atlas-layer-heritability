"""실험 3d y03: 실험 3c 트랙별 Atlas 값 → 계통 × 모달리티 주석 (PROTOCOL_exp3d_2026-10-02) → 염색체별 LD 점수 (실험 3c z03과 같은 방식).
 L5_{g}_{m} (72): 계통 g 모달리티 m 트랙 최댓값 → 백분위 → 상위 5 %.   L5_all_{m} (6): 12계통 전체 트랙.
 L5k_{g}_{m} (72): 칸마다 트랙 최대 3개 무작위(seed 20261002), 같은 방식 (트랙 수 민감도).
 MAF < 5 % 또는 REF/ALT 미결정 SNP = 0. 조각마다 필요한 열만 읽는다.
 시험: Y03_TEST=1 → exp3c/test/part_*.parquet 만, chr22 만, exp3d/test/annot 에 쓴다.
출력: exp3d/annot/<이름>/<이름>.<chr>.{annot.gz,l2.ldscore.gz,l2.M,l2.M_5_50}, exp3d/annot_sd_common.csv, exp3d/cell_tracks.csv, results/y03_annot.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); C3 = Path(f"{DATA}/exp3c"); Y = Path(f"{DATA}/exp3d")
PF = X/"GRCh38"/"plink_files"; TEST = os.environ.get("Y03_TEST") == "1"
SRC = C3/"test" if TEST else C3/"atlas"; AN = Y/"test"/"annot" if TEST else Y/"annot"; AN.mkdir(exist_ok=True, parents=True); (EXP/"results").mkdir(exist_ok=True)
CHRS = [22] if TEST else range(1, 23); LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
MODS = ["DNASE", "H3K27ac", "H3K4me1", "H3K4me3", "CAGE", "RNA"]
TL = pd.read_csv(C3/"track_lineage.csv"); TL["col"] = TL.scorer + "::" + TL.key
rng = np.random.default_rng(20261002)
GROUPS, cells = {}, []
for g in LIN:
    for m in MODS:
        cols = TL[(TL.lineage == g) & (TL["mod"] == m)].col.tolist()
        sub = sorted(rng.choice(cols, 3, replace=False).tolist()) if len(cols) > 3 else cols
        GROUPS[f"L5_{g}_{m}"] = cols; GROUPS[f"L5k_{g}_{m}"] = sub; cells.append(dict(lineage=g, mod=m, n_tracks=len(cols), n_matched=len(sub)))
for m in MODS: GROUPS[f"L5_all_{m}"] = TL[TL.lineage.isin(LIN) & (TL["mod"] == m)].col.tolist()
CT = pd.DataFrame(cells); CT.to_csv(Y/"cell_tracks.csv", index=False)
P("계통 × 모달리티 트랙 수\n" + CT.pivot(index="lineage", columns="mod", values="n_tracks").reindex(index=LIN, columns=MODS).to_string())
P(f"트랙 0개 칸: {[k for k, v in GROUPS.items() if not v] or '없음'}")
NAMES = list(GROUPS); need = sorted(set(c for v in GROUPS.values() for c in v))
parts = sorted(SRC.glob("part_*.parquet")); P(f"주석 {len(NAMES)}개, 조각 {len(parts)}개 읽기")
V = []
for f in parts:
    have = set(pq.read_schema(f).names); cols = [c for c in need if c in have]
    D = pd.read_parquet(f, columns=["atlas_variant"] + cols).set_index("atlas_variant").astype(np.float32)
    V.append(pd.DataFrame({n: (D[[c for c in v if c in have]].max(1) if any(c in have for c in v) else np.nan) for n, v in GROUPS.items()}, index=D.index))
V = pd.concat(V); V = V[~V.index.duplicated()]
T = pd.read_parquet(X/"targets.parquet")[["rsid", "atlas_variant"]]; T = T.join(V, on="atlas_variant")
P(f"대상 {len(T):,}, Atlas 값 있음 {int(T['L5_all_DNASE'].notna().sum()):,}")
ann = pd.DataFrame({"rsid": T.rsid})
for n in NAMES: ann[n] = (T[n].rank(pct=True) >= .95).astype(np.int8)
# 같은 계통 안 모달리티 간 겹침(Jaccard), 계통별 평균
for g in LIN:
    cs = [f"L5_{g}_{m}" for m in MODS]
    J = [((ann[a] & ann[b]).sum() / max((ann[a] | ann[b]).sum(), 1)) for i, a in enumerate(cs) for b in cs[i + 1:]]
    P(f"  {g:24s} 모달리티 간 Jaccard 평균 {np.mean(J):.2f} (범위 {min(J):.2f}–{max(J):.2f})")
J = [((ann[f"L5_{g}_{m}"] & ann[f"L5k_{g}_{m}"]).sum() / max((ann[f"L5_{g}_{m}"] | ann[f"L5k_{g}_{m}"]).sum(), 1)) for g in LIN for m in MODS]
P(f"L5 vs L5k(트랙 맞춤) Jaccard: 중앙값 {np.median(J):.2f}, 최소 {min(J):.2f}")
ann = ann.drop_duplicates("rsid").set_index("rsid")
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
for c in CHRS:
    out_all = AN/"_all"/f"all.{c}"; (AN/"_all").mkdir(exist_ok=True)
    if not Path(f"{out_all}.l2.ldscore.gz").exists():
        b = pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, names=["CHR", "SNP", "CM", "BP", "a1", "a2"])
        ann.reindex(b.SNP).fillna(0).astype(np.int8).reset_index(drop=True).to_csv(f"{out_all}.annot.gz", sep="\t", index=False, compression="gzip")
        r = subprocess.run(["python", "-W", "ignore", str(LDSC), "--l2", "--bfile", str(PF/f"1000G.EUR.hg38.{c}"), "--ld-wind-cm", "1", "--annot", f"{out_all}.annot.gz", "--thin-annot",
                            "--out", str(out_all), "--print-snps", str(X/"hm3_no_MHC.list.txt")], cwd=LDSC.parent, env=env, capture_output=True, text=True)
        if not Path(f"{out_all}.l2.ldscore.gz").exists(): raise SystemExit(f"chr{c} l2 실패\n{r.stdout[-1500:]}")
    L = pd.read_csv(f"{out_all}.l2.ldscore.gz", sep="\t"); Aall = pd.read_csv(f"{out_all}.annot.gz", sep="\t")
    M = open(f"{out_all}.l2.M").read().split(); M5 = open(f"{out_all}.l2.M_5_50").read().split()
    groups = {n: [n] for n in NAMES}
    for g in LIN: groups[f"J_{g}"] = [f"L5_{g}_{m}" for m in MODS]; groups[f"Jk_{g}"] = [f"L5k_{g}_{m}" for m in MODS]   # D2 / D2k 동시 모형용
    for gname, cols in groups.items():
        d = AN/gname; d.mkdir(exist_ok=True); idx = [NAMES.index(x) for x in cols]
        L[["CHR", "SNP", "BP"] + [f"{x}L2" for x in cols]].to_csv(d/f"{gname}.{c}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        Aall[cols].to_csv(d/f"{gname}.{c}.annot.gz", sep="\t", index=False, compression="gzip")
        open(d/f"{gname}.{c}.l2.M", "w").write("\t".join(M[i] for i in idx)); open(d/f"{gname}.{c}.l2.M_5_50", "w").write("\t".join(M5[i] for i in idx))
    P(f"chr{c} 완료")
ann[NAMES].std(ddof=0).rename("sd_common").to_csv(Y/("test_annot_sd_common.csv" if TEST else "annot_sd_common.csv"))
if not TEST: (EXP/"results"/"y03_annot.txt").write_text("\n".join(log), encoding="utf-8")
