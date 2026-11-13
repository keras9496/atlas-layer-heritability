"""실험 3w w08: 106개 GWAS의 농축된 쌍(다섯 묶음, 117개)을 TSS 거리로 나눈 균형 (PROTOCOL_exp3w_추가 A; 사후 분석).
 모형 PD_<묶음>: 전체 트랙의 세 층위 주석을 TSS 1 kb 안(근위)과 밖(원위)으로 나눈 여섯 주석을 함께 적합(3k의 주석, 정렬 블록).
 쌍마다 원위 균형 = atan2(원위 DNase τ*, 원위 CAGE τ*), 근위도 같은 식. 집계는 w05·w07과 같다(역분산 가중 평균, 공동 jackknife).
실행: python -I -X utf8 w08_pd.py
출력: results/w08_pairs.csv, w08_summary.txt
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy import stats
from wlib_panel import *
assert PD_READY, "근위·원위 적합이 아직 다 끝나지 않았다"
log = []
def P(*x): s = " ".join(str(i) for i in x); print(s, flush=True); log.append(s)
def pz(d, se): return f"{d:+.1f}° ± {se:.1f}° (P = {2 * stats.norm.sf(abs(d / se)):.1e})"
F5 = [(i, g) for i, g in PAIRS if g in FIVE]; assert len(F5) == 117
rows = []
for i, g in F5:
    rec = dict(id=i, trait=TR.trait[i], category=cat(i), group=g)
    for lab in ("all", "dist", "prox"):
        x = D[(i, g, lab)]; rec.update({f"{lab}_cage": x["c"], f"{lab}_dnase": x["d"], f"{lab}_c1": x["c"] - x["d"], f"{lab}_c1_se": jse(x["cdv"] - x["ddv"]), f"{lab}_phiD": x["e"], f"{lab}_phiD_se": jse(x["dv"]), f"{lab}_valid": x["ok"]})
    rows.append(rec)
PR = pd.DataFrame(rows); PR.to_csv(RES/"w08_pairs.csv", index=False)
P(f"다섯 묶음의 농축된 쌍 {len(PR)}개. 균형이 정의되는 쌍: 전체 트랙 {int(PR.all_valid.sum())}, 원위 {int(PR.dist_valid.sum())}, 근위 {int(PR.prox_valid.sum())}")
for lab in ("all", "dist", "prox"): P(f"  쌍의 s.e. 중앙값 [{NAME[lab]}]: 균형 {PR[PR[lab + '_valid']][lab + '_phiD_se'].median():.1f}°, CAGE − DNase {PR[lab + '_c1_se'].median():.2f}")
# ---- 조직 묶음 주석마다 ----
P("\n== 조직 묶음 주석마다의 균형: 역분산 가중 평균 ± s.e. (중앙값; 쌍 수)")
A = {lab: gmeans(lab, FIVE) for lab in ("all", "dist", "prox")}
for g in FIVE:
    P(f"  {KO[g]:8s} " + " | ".join(f"{NAME[lab]} {A[lab][g][0]:5.1f}° ± {jse(A[lab][g][1]):.1f}° ({np.median([D[(i, gg, lab)]['e'] for i, gg in valid(lab, [g])]):.1f}°; {A[lab][g][2]})" for lab in ("all", "dist", "prox")))
P("이질성 Q (공동 jackknife):")
for lab in ("all", "dist", "prox"):
    for gs, note in ((FIVE, "다섯 주석"), (NON, "stromal을 뺀 네 주석")):
        Q, p, pc, pf = Qtest([A[lab][g][0] for g in gs], np.column_stack([A[lab][g][1] for g in gs])); P(f"  [{NAME[lab]}] {note}: Q = {Q:.1f} (자유도 {p}), P = {pc:.1e} (F 근사 {pf:.1e})")
P("stromal − 다른 주석 (공동 jackknife):")
def others(lab): return float(np.mean([A[lab][g][0] for g in NON])), np.mean([A[lab][g][1] for g in NON], axis=0)
for lab in ("all", "dist", "prox"):
    s = [f"{KO[g]} " + pz(A[lab]["stromal"][0] - A[lab][g][0], jse(A[lab]["stromal"][1] - A[lab][g][1])) for g in NON]; om, odv = others(lab)
    P(f"  [{NAME[lab]}] " + "; ".join(s)); P(f"  [{NAME[lab]}] stromal − 네 주석 평균의 평균: " + pz(A[lab]["stromal"][0] - om, jse(A[lab]["stromal"][1] - odv)))
dm = {lab: (A[lab]["stromal"][0] - others(lab)[0], A[lab]["stromal"][1] - others(lab)[1]) for lab in ("dist", "prox")}
P("  원위의 차이 − 근위의 차이 (stromal − 네 주석 평균): " + pz(dm["dist"][0] - dm["prox"][0], jse(dm["dist"][1] - dm["prox"][1])))
for g in FIVE: P(f"  {KO[g]}: 원위 − 근위 = " + pz(A["dist"][g][0] - A["prox"][g][0], jse(A["dist"][g][1] - A["prox"][g][1])))
# ---- 같은 GWAS 안 ----
P("\n== 같은 GWAS 안: stromal − 그 GWAS의 다른 묶음(다섯 가운데) 평균")
W_ = {}
for lab in ("all", "dist", "prox"):
    v = valid(lab, FIVE); items = []
    for i in sorted({i for i, g in v if g == "stromal"}):
        oth = [g for j, g in v if j == i and g != "stromal"]
        if oth: items.append((D[(i, "stromal", lab)]["e"] - np.mean([D[(i, g, lab)]["e"] for g in oth]), D[(i, "stromal", lab)]["dv"] - np.mean([D[(i, g, lab)]["dv"] for g in oth], axis=0)))
    m, dv = ivw(items); W_[lab] = (m, dv); P(f"  [{NAME[lab]}] GWAS {len(items)}개: " + pz(m, jse(dv)))
P("  원위 − 근위: " + pz(W_["dist"][0] - W_["prox"][0], jse(W_["dist"][1] - W_["prox"][1])))
# ---- CAGE − DNase 대비 (균형이 정의되지 않는 쌍도 들어간다) ----
P("\n== CAGE − DNase 대비의 역분산 가중 평균 (모든 쌍 117개)")
def c1means(lab):
    out = {}
    for g in FIVE:
        m, dv = ivw([(D[(i, g, lab)]["c"] - D[(i, g, lab)]["d"], D[(i, g, lab)]["cdv"] - D[(i, g, lab)]["ddv"]) for i, gg in F5 if gg == g]); out[g] = (m, dv)
    return out
C = {lab: c1means(lab) for lab in ("all", "dist", "prox")}
for g in FIVE: P(f"  {KO[g]:8s} " + " | ".join(f"{NAME[lab]} {C[lab][g][0]:+.2f} ± {jse(C[lab][g][1]):.2f}" for lab in ("all", "dist", "prox")))
for lab in ("all", "dist", "prox"):
    om = float(np.mean([C[lab][g][0] for g in NON])); odv = np.mean([C[lab][g][1] for g in NON], axis=0); d = C[lab]["stromal"][0] - om; se = jse(C[lab]["stromal"][1] - odv)
    Q, p, pc, pf = Qtest([C[lab][g][0] for g in FIVE], np.column_stack([C[lab][g][1] for g in FIVE])); Q4, p4, pc4, _ = Qtest([C[lab][g][0] for g in NON], np.column_stack([C[lab][g][1] for g in NON]))
    P(f"  [{NAME[lab]}] stromal − 네 주석 평균: {d:+.2f} ± {se:.2f} (P = {2 * stats.norm.sf(abs(d / se)):.1e}); 다섯 주석 Q = {Q:.1f} (P = {pc:.1e}); 네 주석 Q = {Q4:.1f} (P = {pc4:.1e})")
# ---- 쌍 수준: 부호 ----
P("\n== 쌍 수준: DNase τ* > CAGE τ*인 쌍의 수 (stromal / 다른 네 묶음)")
for lab in ("all", "dist", "prox"):
    a = PR[PR.group == "stromal"]; b = PR[PR.group != "stromal"]; P(f"  [{NAME[lab]}] stromal {int((a[lab + '_c1'] < 0).sum())}/{len(a)}; 다른 네 묶음 {int((b[lab + '_c1'] < 0).sum())}/{len(b)}")
(RES/"w08_summary.txt").write_text("\n".join(log), encoding="utf-8")
