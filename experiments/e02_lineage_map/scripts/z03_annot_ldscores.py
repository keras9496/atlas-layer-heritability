"""실험 3c z03: 트랙별 Atlas 값 → 계통 주석(PROTOCOL_exp3c) → 염색체별 LD 점수 (실험 3 x03과 같은 방식).
 v_{g,m} = 계통 g 의 모달리티 m 트랙 |ALT − REF| 최댓값 (m ∈ DNASE, H3K27ac, H3K4me1, CAGE). 대상 SNP 전체에서 백분위.
 S_g = 4개 백분위 평균 → L5_g = 1[S_g 백분위 ≥ 0.95] (12개). L5_all = 12계통 합집합 트랙으로 같은 방식.
 MAF < 5 % 또는 REF/ALT 미결정 SNP = 0. 조각마다 필요한 열만 읽는다(메모리 절약).
 시험: Z03_TEST=1 → exp3c/test/part_*.parquet 만 쓰고 chr22 만, exp3c/test/annot 에 쓴다.
출력: exp3c/annot/<이름>/<이름>.<chr>.{annot.gz,l2.ldscore.gz,l2.M,l2.M_5_50}, exp3c/lineage_scores.parquet (S_g 백분위, 나중 분석용),
      exp3c/annot_sd_common.csv, results/z03_annot.txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import subprocess, os
import numpy as np, pandas as pd
import pyarrow.parquet as pq
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3c"); PF = X/"GRCh38"/"plink_files"
TEST = os.environ.get("Z03_TEST") == "1"; SRC = Y/"test" if TEST else Y/"atlas"; AN = Y/"test"/"annot" if TEST else Y/"annot"; AN.mkdir(exist_ok=True, parents=True)
CHRS = [22] if TEST else range(1, 23); LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py")
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
LIN = ["CNS_tissue", "neural_cells", "myeloid", "lymphoid", "liver", "GI", "cardiovascular", "skeletal_muscle", "adipose", "stromal", "reproductive_endocrine", "epithelial_barrier"]
CORE = ["DNASE", "H3K27ac", "H3K4me1", "CAGE"]
TL = pd.read_csv(Y/"track_lineage.csv"); TL["col"] = TL.scorer + "::" + TL.key
GROUPS = {(g, m): TL[(TL.lineage == g) & (TL["mod"] == m)].col.tolist() for g in LIN for m in CORE}
for m in CORE: GROUPS[("all", m)] = TL[TL.lineage.isin(LIN) & (TL["mod"] == m)].col.tolist()
empty = [k for k, v in GROUPS.items() if not v]; P(f"계통 × 모달리티 {len(GROUPS)}칸, 트랙 없는 칸: {empty or '없음'}")
need = sorted(set(c for v in GROUPS.values() for c in v))
parts = sorted(SRC.glob("part_*.parquet")); P(f"조각 {len(parts)}개 읽기")
V = []
for f in parts:
    have = set(pq.read_schema(f).names); cols = [c for c in need if c in have]
    D = pd.read_parquet(f, columns=["atlas_variant"] + cols).set_index("atlas_variant").astype(np.float32)
    V.append(pd.DataFrame({f"{g}|{m}": (D[[c for c in v if c in have]].max(1) if any(c in have for c in v) else np.nan) for (g, m), v in GROUPS.items()}, index=D.index))
V = pd.concat(V); V = V[~V.index.duplicated()]
T = pd.read_parquet(X/"targets.parquet")[["rsid", "atlas_variant"]]; T = T.join(V, on="atlas_variant")
P(f"대상 {len(T):,}, Atlas 값 있음 {int(T['all|DNASE'].notna().sum()):,}")
ann = pd.DataFrame({"rsid": T.rsid}); sc = pd.DataFrame({"rsid": T.rsid})
for g in LIN + ["all"]:
    S = pd.concat([T[f"{g}|{m}"].rank(pct=True) for m in CORE], axis=1).mean(1); Sp = S.rank(pct=True)
    sc[g] = Sp.astype(np.float32); ann[f"L5_{g}"] = (Sp >= .95).astype(np.int8)
sc.to_parquet(Y/("test_lineage_scores.parquet" if TEST else "lineage_scores.parquet"), index=False)
NAMES = [c for c in ann.columns if c != "rsid"]
J = pd.DataFrame([[((ann[a] & ann[b]).sum() / max((ann[a] | ann[b]).sum(), 1)) for b in NAMES] for a in NAMES], index=NAMES, columns=NAMES).round(2)
P("상위 5 % 주석 간 겹침(Jaccard)\n" + J.to_string())
if not TEST:   # 실험 3 ct5_myeloid 와의 겹침(기술)
    old = pd.read_parquet(X/"targets.parquet")[["atlas_variant"]].join(pd.concat([pd.read_parquet(f) for f in sorted((X/"atlas").glob("part_*.parquet"))]).drop_duplicates("atlas_variant").set_index("atlas_variant"), on="atlas_variant")
    o = (pd.concat([old[f"{m}_myeloid"].astype(float).rank(pct=True) for m in CORE], axis=1).mean(1).rank(pct=True) >= .95)
    P(f"L5_myeloid(3c) vs ct5_myeloid(실험 3) Jaccard {((ann.L5_myeloid == 1) & o).sum() / ((ann.L5_myeloid == 1) | o).sum():.2f}")
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
    for i, n in enumerate(NAMES):
        d = AN/n; d.mkdir(exist_ok=True)
        L[["CHR", "SNP", "BP", f"{n}L2"]].to_csv(d/f"{n}.{c}.l2.ldscore.gz", sep="\t", index=False, compression="gzip")
        Aall[[n]].to_csv(d/f"{n}.{c}.annot.gz", sep="\t", index=False, compression="gzip")
        open(d/f"{n}.{c}.l2.M", "w").write(M[i]); open(d/f"{n}.{c}.l2.M_5_50", "w").write(M5[i])
    P(f"chr{c} 완료")
ann[NAMES].std(ddof=0).rename("sd_common").to_csv(Y/("test_annot_sd_common.csv" if TEST else "annot_sd_common.csv"))
if not TEST: (EXP/"results"/"z03_annot.txt").write_text("\n".join(log), encoding="utf-8")
