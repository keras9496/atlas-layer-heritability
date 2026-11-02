"""실험 3m M3: 공선성과 계수 안정성 (PROTOCOL_exp3m_2026-10-17, 기술; 새 실행 없음).
 계통 5개(골수계, 뇌 조직, 림프계, 간, 결합조직) × 주석 두 벌(전체 트랙 L5 = 주 지도 T3all, 세포형 일치 M3 = 3i).
 (a) 대상 SNP 5,828,323개에서 세 이진 주석의 φ 상관과 Jaccard 겹침.
 (b) 회귀 SNP(baselineLD ∩ HapMap3 가중치, 1,186,879개)에서 세 층위 LD 점수 상관: 원래 값, baselineLD 97개(+ 상수)로 회귀한 잔차(편상관).
 (c) VIF: 전체(층위 j를 기준선 97 + 다른 두 층위로 회귀한 R²) / 층위 사이(잔차 상관 행렬 역행렬의 대각). 표준화 설계 행렬 조건수(기준선만, 기준선 + 3층위).
 (d) 공동 모형 τ* 사이 jackknife 상관(3j 정렬 블록, T3all / M3).
 (e) 단일 층위 모형(3d D1) → 동시 모형(T3all) τ* 변화, 지도 형질 11개.
 읽는 규칙(미리 정함): 층위 사이 VIF > 5 또는 편상관 |r| > 0.9 → 그 계통의 동시 계수는 잘 분리되지 않는다.
출력: results/m03_collinearity.csv, m03_coef_corr.csv, m03_single_vs_joint.csv, m03_tests.txt
"""
import sys, re
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent)); from m_common import *
EXP = Path(__file__).resolve().parents[1]; PF = X/"GRCh38"/"plink_files"; BL = X/"GRCh38"/"baselineLD_v2.2"; W = X/"GRCh38"/"weights"
LIN = ["myeloid", "CNS_tissue", "lymphoid", "liver", "stromal"]
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
tg = set(pd.read_parquet(X/"targets.parquet", columns=["rsid"]).rsid)
def ann_cols(g, kind):
    if kind == "L5": return [(Y/"annot"/f"L5_{g}_{m}", f"L5_{g}_{m}", [f"L5_{g}_{m}"]) for m in LAY]
    return [(I/"annot"/f"M3_{g}", f"M3_{g}", [f"M5_{g}_{m}" for m in LAY])]
# ---------- SNP 단위 주석 (a) ----------
def annot_matrix(g, kind):
    out = []
    for c in range(1, 23):
        b = pd.read_csv(PF/f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, usecols=[1], names=["SNP"]).SNP
        keep = b.isin(tg).to_numpy(); cols = []
        for d, stem, names in ann_cols(g, kind):
            A = pd.read_csv(d/f"{stem}.{c}.annot.gz", sep="\t"); cols.append(A[names].to_numpy()[keep])
        out.append(np.hstack(cols))
    return np.vstack(out).astype(np.float32)
# ---------- 회귀 SNP 단위 LD 점수 (b, c) ----------
wsnp = pd.concat([pd.read_csv(W/f"weights.hm3_noMHC.{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"]) for c in range(1, 23)]).SNP.drop_duplicates()
P("기준선 LD 점수 읽는 중…")
BLD = pd.concat([pd.read_csv(BL/f"baselineLD.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP")
BLD = BLD[BLD.SNP.isin(set(wsnp))].reset_index(drop=True); bcols = [c for c in BLD.columns if c not in ("CHR", "SNP", "BP")]; assert len(bcols) == 97, len(bcols)
P(f"회귀 SNP {len(BLD):,} (3j 기준 SNP 1,186,879와 같아야 함)"); assert len(BLD) == 1186879
Xb = BLD[bcols].to_numpy(np.float64); Xb1 = np.column_stack([np.ones(len(Xb)), Xb])
def cond(M): Z = (M - M.mean(0)) / M.std(0); ev = np.linalg.eigvalsh(np.corrcoef(Z, rowvar=False)); return float(np.sqrt(ev.max() / max(ev.min(), 1e-300)))
cond_b = cond(Xb); P(f"기준선 97개만의 조건수 {cond_b:.0f}")
def ld_matrix(g, kind):
    cols = []
    for d, stem, names in ann_cols(g, kind):
        L = pd.concat([pd.read_csv(d/f"{stem}.{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]).drop_duplicates("SNP").set_index("SNP")
        cols.append(L.loc[BLD.SNP, [f"{n}L2" for n in names]].to_numpy(np.float64))
    return np.hstack(cols)
rows = []
for kind in ["L5", "M3"]:
    for g in LIN:
        A = annot_matrix(g, kind); phi = np.corrcoef(A, rowvar=False)
        jac = {(i, j): float(((A[:, i] > 0) & (A[:, j] > 0)).sum() / ((A[:, i] > 0) | (A[:, j] > 0)).sum()) for i, j in [(0, 1), (0, 2), (1, 2)]}
        Lm = ld_matrix(g, kind); rawr = np.corrcoef(Lm, rowvar=False)
        beta = np.linalg.lstsq(Xb1, Lm, rcond=None)[0]; Rz = Lm - Xb1 @ beta; pr = np.corrcoef(Rz, rowvar=False)
        vif_lay = np.diag(np.linalg.inv(pr))
        vif_full = []
        for j in range(3):
            oth = np.column_stack([Xb1, np.delete(Lm, j, axis=1)]); e = Lm[:, j] - oth @ np.linalg.lstsq(oth, Lm[:, j], rcond=None)[0]
            vif_full.append(float(Lm[:, j].var() / e.var()))
        cnd = cond(np.column_stack([Xb, Lm])); r2b = 1 - Rz.var(0) / Lm.var(0)
        flag = bool((vif_lay > 5).any() or (np.abs(pr[np.triu_indices(3, 1)]) > .9).any())
        rec = dict(kind=kind, lineage=g, size_CAGE=A[:, 0].mean(), size_DNASE=A[:, 1].mean(), size_H3K27ac=A[:, 2].mean(),
                   phi_CD=phi[0, 1], phi_CH=phi[0, 2], phi_DH=phi[1, 2], jac_CD=jac[(0, 1)], jac_CH=jac[(0, 2)], jac_DH=jac[(1, 2)],
                   ldr_CD=rawr[0, 1], ldr_CH=rawr[0, 2], ldr_DH=rawr[1, 2], pr_CD=pr[0, 1], pr_CH=pr[0, 2], pr_DH=pr[1, 2],
                   r2base_CAGE=r2b[0], r2base_DNASE=r2b[1], r2base_H3K27ac=r2b[2],
                   viflay_CAGE=vif_lay[0], viflay_DNASE=vif_lay[1], viflay_H3K27ac=vif_lay[2], viffull_CAGE=vif_full[0], viffull_DNASE=vif_full[1], viffull_H3K27ac=vif_full[2],
                   cond_base=cond_b, cond_full=cnd, flag=flag)
        rows.append(rec)
        P(f"[{kind} {g}] 크기 C/D/H {A[:, 0].mean():.3f}/{A[:, 1].mean():.3f}/{A[:, 2].mean():.3f} | φ CD {phi[0, 1]:.2f} CH {phi[0, 2]:.2f} DH {phi[1, 2]:.2f} | Jaccard CD {jac[(0, 1)]:.2f} CH {jac[(0, 2)]:.2f} DH {jac[(1, 2)]:.2f}")
        P(f"    LD 상관 원래 CD {rawr[0, 1]:.2f} CH {rawr[0, 2]:.2f} DH {rawr[1, 2]:.2f} | 기준선 뺀 편상관 CD {pr[0, 1]:.2f} CH {pr[0, 2]:.2f} DH {pr[1, 2]:.2f} | 기준선 R² C/D/H {r2b[0]:.2f}/{r2b[1]:.2f}/{r2b[2]:.2f}")
        P(f"    VIF 층위 사이 C/D/H {vif_lay[0]:.2f}/{vif_lay[1]:.2f}/{vif_lay[2]:.2f} | 전체 {vif_full[0]:.1f}/{vif_full[1]:.1f}/{vif_full[2]:.1f} | 조건수 기준선 {cond_b:.0f} → +3층위 {cnd:.0f} | 규칙 표시 {'예' if flag else '아니오'}")
pd.DataFrame(rows).to_csv(EXP/"results"/"m03_collinearity.csv", index=False)
# ---------- (d) 공동 모형 τ* 사이 jackknife 상관 ----------
P("\n(d) 공동 모형 τ* 사이 jackknife 상관 (정렬 블록)")
cr = []
for kind, model in [("L5", lambda t: "T3all"), ("M3", lambda t: f"M3_{OWN[t]}")]:
    for t in ALL:
        L = layers(t, [model(t)]); Dm = np.column_stack([L[m][1] for m in LAY]); C = jcov(Dm); s = np.sqrt(np.diag(C)); R = C / np.outer(s, s)
        cr.append(dict(kind=kind, trait=t, lineage=OWN[t], **{f"tau_{m}": L[m][0] for m in LAY}, **{f"se_{m}": s[i] for i, m in enumerate(LAY)}, r_CD=R[0, 1], r_CH=R[0, 2], r_DH=R[1, 2]))
        P(f"  [{kind}] {KO[t]:6s} τ* C/D/H {L['CAGE'][0]:+.2f}/{L['DNASE'][0]:+.2f}/{L['H3K27ac'][0]:+.2f} | r(C,D) {R[0, 1]:+.2f} r(C,H) {R[0, 2]:+.2f} r(D,H) {R[1, 2]:+.2f}")
pd.DataFrame(cr).to_csv(EXP/"results"/"m03_coef_corr.csv", index=False)
# ---------- (e) 단일 층위 → 동시 모형 ----------
P("\n(e) 단일 층위 모형(3d D1) → 동시 모형(T3all) τ*")
sv = []
for t in MAP:
    g = OWN[t]; J3 = layers(t, ["T3all"])
    for m in LAY:
        out = Y/"sldsc"/f"{t}__D1__L5_{g}_{m}"; R = pd.read_csv(f"{out}.results", sep="\t"); R["n"] = [strip(n) for n in R.Category]
        h2 = float(re.search(r"Total Observed scale h2: ([\d.eE-]+)", open(f"{out}.log", encoding="utf-8").read()).group(1))
        s1 = float(R.set_index("n").loc[f"L5_{g}_{m}", "Coefficient"]) * SD[f"L5_{g}_{m}"] * M_TOT / h2
        sv.append(dict(trait=t, lineage=g, layer=m, single=s1, joint=J3[m][0], joint_se=jse(J3[m][1]), sign_change=bool(np.sign(s1) != np.sign(J3[m][0]) and abs(J3[m][0]) > 0)))
S = pd.DataFrame(sv); S.to_csv(EXP/"results"/"m03_single_vs_joint.csv", index=False)
for t in MAP:
    s = S[S.trait == t].set_index("layer")
    P(f"  {KO[t]:6s} " + " | ".join(f"{m} {s.loc[m, 'single']:+.2f} → {s.loc[m, 'joint']:+.2f} ({s.loc[m, 'joint_se']:.2f}){' 부호 바뀜' if s.loc[m, 'sign_change'] else ''}" for m in LAY))
P(f"  부호 바뀜 {int(S.sign_change.sum())}/{len(S)}; 동시 모형 음수 τ* {int((S.joint < 0).sum())}개, 그중 |τ*| > 1.96 SE {int(((S.joint < 0) & (S.joint.abs() > 1.96 * S.joint_se)).sum())}개")
(EXP/"results"/"m03_tests.txt").write_text("\n".join(log), encoding="utf-8")
