"""실험 3k k01: 층위 주석의 근위·원위 분리 (PROTOCOL_exp3k_2026-10-15) → 염색체별 LD 점수 (3d·3i와 같은 l2 호출).
 TSS = GENCODE v46 종합 주석, 단백질 암호화·lncRNA 유전자의 모든 전사체 시작 위치(+ start, − end), 1–22번. 근위 = 가장 가까운 TSS ≤ 1,000 bp.
 P_{g}_{m} = 3d L5_{g}_{m} ∧ 근위, D_{g}_{m} = L5 ∧ 원위 (다시 순위를 매기지 않는다). 묶음 PD_{g} = [P × 3층위, D × 3층위].
 확인(어긋나면 멈춘다): P + D = L5 (모든 SNP), 다시 계산한 L5의 sd = 3d annot_sd_common (차이 < 1e-9).
 시험: K01_TEST=1 → chr22만, exp3k/test/annot.
출력: exp3k/annot/PD_{g}/..., exp3k/annot_sd_common.csv, exp3k/target_tss_dist.parquet, results/k01_annot[_test].txt
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import os, re, gzip, subprocess
import numpy as np, pandas as pd
from pathlib import Path
EXP = Path(__file__).resolve().parents[1]; X = Path(f"{DATA}/exp3"); Y = Path(f"{DATA}/exp3d"); K = Path(f"{DATA}/exp3k")
PF = X/"GRCh38"/"plink_files"; LDSC = Path(f"{DATA}/tools/ldsc/ldsc.py"); GTF = Path(f"{DATA}/ref/gencode.v46.annotation.gtf.gz")
TEST = os.environ.get("K01_TEST") == "1"; AN = K/"test"/"annot" if TEST else K/"annot"; AN.mkdir(parents=True, exist_ok=True)
CHRS = [22] if TEST else list(range(1, 23)); W = 1000
LING = ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]; MODS = ["CAGE", "DNASE", "H3K27ac"]
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
# --- TSS
rows = []
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.split("\t", 9)
        if c[2] != "transcript" or c[0] not in {f"chr{i}" for i in range(1, 23)}: continue
        gt = re.search(r'gene_type "([^"]+)"', c[8]).group(1)
        if gt in ("protein_coding", "lncRNA"): rows.append((int(c[0][3:]), int(c[3]) if c[6] == "+" else int(c[4]), gt))
TSS = pd.DataFrame(rows, columns=["chr", "tss", "gene_type"]).drop_duplicates(["chr", "tss"])
P(f"TSS (GENCODE v46, 단백질 암호화·lncRNA 전사체, 1–22번): 고유 위치 {len(TSS):,} (단백질 암호화 {int((TSS.gene_type == 'protein_coding').sum()):,}, lncRNA {int((TSS.gene_type == 'lncRNA').sum()):,})")
# --- 대상 SNP의 TSS 거리
T = pd.read_parquet(X/"targets.parquet", columns=["chr", "pos", "rsid"]); T["chr"] = T["chr"].astype(int); T = T.drop_duplicates("rsid")
T["dist"] = np.inf
for c, g in T.groupby("chr"):
    t = np.sort(TSS.tss[TSS.chr == c].to_numpy()); p = g.pos.to_numpy(); i = np.searchsorted(t, p)
    lo = np.abs(p - t[np.clip(i - 1, 0, len(t) - 1)]); hi = np.abs(t[np.clip(i, 0, len(t) - 1)] - p); T.loc[g.index, "dist"] = np.minimum(lo, hi)
T["prox"] = T.dist <= W
if not TEST: T[["rsid", "chr", "pos", "dist"]].to_parquet(K/"target_tss_dist.parquet", index=False)
P(f"대상 SNP {len(T):,}: 근위(≤ {W:,} bp) {int(T.prox.sum()):,} ({T.prox.mean() * 100:.1f} %)")
PROX = T.set_index("rsid").prox
# --- 염색체별 분리 + LD 점수
NAMES = [f"{k}_{g}_{m}" for g in LING for k in ["P", "D"] for m in MODS]
groups = {f"PD_{g}": [f"{k}_{g}_{m}" for k in ["P", "D"] for m in MODS] for g in LING}
env = dict(os.environ, PYTHONIOENCODING="utf-8", MSYS_NO_PATHCONV="1")
TV = {n: [] for n in NAMES + [f"L5_{g}_{m}" for g in LING for m in MODS]}   # 대상 SNP 수준 값 (sd 계산용)
for c in CHRS:
    b = pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, names=["CHR", "SNP", "CM", "BP", "a1", "a2"])
    px = PROX.reindex(b.SNP).fillna(False).to_numpy().astype(bool)
    Ac = pd.DataFrame(index=b.index)
    for g in LING:
        for m in MODS:
            L5 = pd.read_csv(Y/"annot"/f"L5_{g}_{m}"/f"L5_{g}_{m}.{c}.annot.gz", sep="\t").iloc[:, 0].to_numpy().astype(np.int8); assert len(L5) == len(b)
            Ac[f"P_{g}_{m}"] = (L5 * px).astype(np.int8); Ac[f"D_{g}_{m}"] = (L5 * ~px).astype(np.int8)
            assert ((Ac[f"P_{g}_{m}"] + Ac[f"D_{g}_{m}"]).to_numpy() == L5).all(), "P + D ≠ L5"
            Ac[f"L5_{g}_{m}"] = L5
    tr = T[T.chr == c].rsid; V = Ac.assign(SNP=b.SNP).drop_duplicates("SNP").set_index("SNP").reindex(tr)
    assert V.isna().sum().sum() == 0, "bim에 없는 대상 SNP"
    for n in TV: TV[n].append(V[n].to_numpy())
    Ac = Ac[NAMES]; out_all = AN/"_all"/f"all.{c}"; (AN/"_all").mkdir(exist_ok=True)
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
A = pd.DataFrame({n: np.concatenate(v) for n, v in TV.items()}); sd = A.std(ddof=0)
if not TEST:
    ref = pd.read_csv(Y/"annot_sd_common.csv", index_col=0).sd_common; dif = (sd[[f"L5_{g}_{m}" for g in LING for m in MODS]] - ref.reindex([f"L5_{g}_{m}" for g in LING for m in MODS])).abs()
    P(f"L5 sd 재현: 최대 차이 {dif.max():.1e}"); assert dif.max() < 1e-9, "sd 재현 실패 — 멈춘다"
    sd[NAMES].rename("sd_common").to_csv(K/"annot_sd_common.csv")
P("\n주석 크기 (대상 SNP 중 %) 와 L5 주석의 근위 비율")
for g in LING:
    P(f"  {g:11s} " + "  ".join(f"{m} P {A[f'P_{g}_{m}'].mean() * 100:.2f} / D {A[f'D_{g}_{m}'].mean() * 100:.2f} (근위 {A[f'P_{g}_{m}'].sum() / max(A[f'L5_{g}_{m}'].sum(), 1) * 100:.0f} %)" for m in MODS))
(EXP/"results"/("k01_annot_test.txt" if TEST else "k01_annot.txt")).write_text("\n".join(log), encoding="utf-8")
