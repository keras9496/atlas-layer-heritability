"""실험 3u u00: 형질 표와 입력 점검 (PROTOCOL_exp3u). 회귀는 하지 않는다.
 입력: exp3u/raw/sumstats_107 (Zenodo 10515792 sumstats_indep107.tgz, md5 8406e1370e5807a11f9f04174ce9f76e)
 출력: <실험 폴더>/traits_exp3u.csv (형질 묶음, 예측, 포함 여부; 실행 전 고정), results/u00_inputs.csv (파일 점검)
 실행: python -I -X utf8 u00_inputs.py
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import sys, gzip
import numpy as np, pandas as pd
from pathlib import Path
DATA = Path(DATA); EXP = Path(__file__).resolve().parents[1]; RES = EXP/"results"; RES.mkdir(exist_ok=True)
RAW = DATA/"exp3u"/"raw"/"sumstats_107"; J = DATA/"exp3j"
# 묶음: S = 골격, BI = 혈액·면역, LL = 지질, BR = 뇌·행동, U = 예측 없음.  D = 예측을 세운 14형질과 같은 형질
S = ["PASS.Knee_Osteoarthritis.Boer2021", "PASS.Total_Hip_Replacement.Boer2021", "PASS.Thumb_Osteoarthritis.Boer2021", "PASS.Height.Yengo2022", "UKB_460K.bmd_HEEL_TSCOREz"]
BI = ["GBMI.Asthma", "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP", "PASS.Inflammatory_Bowel_Disease.deLange2017", "PASS.Rheumatoid_Arthritis.Ishigaki2022",
      "PASS.Rheumatoid_Arthritis.Saevarsdottir2022", "PASS.Type_1_Diabetes.Chiou2021", "PASS.Primary_Biliary_Cirrhosis.Cordell2015", "PASS.Lupus.Bentham2015", "PASS.Celiac.Dubois2010",
      "UKB_460K.blood_WHITE_COUNT", "UKB_460K.blood_RED_COUNT", "PASS.Eosinophil_Percentage.Vuckovic2020", "UKB_460K.blood_RBC_DISTRIB_WIDTH", "PASS.Reticulocyte_Count.Vuckovic2020",
      "PASS.Platelet_Distribution_Width.Vuckovic2020", "PASS.Monocyte_Percentage.Vuckovic2020", "PASS.Basophil_Percentage.Vuckovic2020"]
LL = ["UKB_460K.biochemistry_Cholesterol", "PASS.HDL.Teslovich2010", "PASS.LDL.Teslovich2010"]
BR = ["PASS.Schizophrenia.Trubetskoy2022", "PASS.Depression.Howard2019", "PASS.Bipolar_Disorder_ALL.Mullins2021", "PASS.Attention_Deficit_Hyperactivity_Disorder_ADHD.Demontis2023",
      "PASS.Insomnia.Watanabe2022", "PASS.Anorexia_Nervosa.Watson2019", "PASS.Alzheimers_Disease.Wightman2021", "PASS.Genetic_Generalized_Epilepsy.Berkovic2022", "PASS.Parkinson.Nalls2019",
      "PASS.Autism.PGC2013", "PASS.Tourette_Syndrome.Yu2019", "PASS.BMI.Yengo2018", "PASS.Intelligence.SavageJansen2018", "PASS.Daytime_Napping.Dashti2021", "UKB_460K.other_MORNINGPERSON",
      "PASS.Reaction_Time.Davies2018", "PASS.Sleep_Duration.Dashti2019", "PASS.General_Risk_Tolerance.KarlssonLinner2019", "PASS.Drinks_Per_Week.Liu2019", "PASS.Alcohol_use.SanchezRoige2018",
      "PASS.Education_Years.Rietveld2013", "PASS.Cigarettes_Per_Day.Liu2019", "PASS.Cannabis_Use_Disorder.Johnson2020", "PASS.Ever_Smoked.TAGC2010", "PASS.Structural_Connectivity_Global.Tissink2023"]
# 예측을 세운 14형질과 같은 형질 (표본이 같거나 겹치거나, 같은 형질의 독립 표본)
D = {"PASS.Height.Yengo2022": "HEIGHT", "UKB_460K.bmd_HEEL_TSCOREz": "BMD", "PASS.BMI.Yengo2018": "BMI", "PASS.Schizophrenia.Trubetskoy2022": "SCZ", "PASS.Bipolar_Disorder_ALL.Mullins2021": "BIP",
     "PASS.Education_Years.Rietveld2013": "EA", "PASS.LDL.Teslovich2010": "LDL", "UKB_460K.biochemistry_Cholesterol": "TC", "PASS.Inflammatory_Bowel_Disease.deLange2017": "IBD",
     "PASS.Rheumatoid_Arthritis.Ishigaki2022": "RA", "PASS.Rheumatoid_Arthritis.Saevarsdottir2022": "RA", "PASS.Alzheimers_Disease.Wightman2021": "AD", "PASS.Monocyte_Percentage.Vuckovic2020": "MONO"}
EXCL = {"PASS.Multiple_Sclerosis.IMSGC2019": "파일 없음(IMSGC 승인 신청 필요)"}
# 폐활량·서혜부 탈장·하지정맥류 2개는 제외가 아니라 U 묶음(예측 없음)이다.
# 같은 형질의 독립 표본 쌍 (재현 기술용)
PAIR = {"PASS.Rheumatoid_Arthritis.Ishigaki2022": "RA", "PASS.Rheumatoid_Arthritis.Saevarsdottir2022": "RA", "PASS.Atrial_Fibrillation.Nielsen2018": "AF", "PANUKBB.Atrial_Fibrillation_And_Flutter": "AF",
        "PASS.Glaucoma.Craig2020": "GLAUCOMA", "GBMI.Primary_Apen_Angle_Glaucoma": "GLAUCOMA", "PASS.Cancer_Prostate.Schumacher2018": "PROSTATE", "UKB_460K.cancer_PROSTATE": "PROSTATE",
        "PASS.Varicose_Veins.Ahmed2022": "VARICOSE", "PANUKBB.Varicose_Vein_Surgery": "VARICOSE"}
M = pd.read_excel(RAW/"traits_indep107.xlsx", header=0); M.columns = ["trait", "file", "link", "mean_N", "h2_z", "note"]; M["file"] = M.file.str.strip()
assert len(M) == 107 and M.file.is_unique
for name, lst in [("S", S), ("BI", BI), ("LL", LL), ("BR", BR)]: assert set(lst) <= set(M.file), (name, set(lst) - set(M.file))
assert set(D) <= set(M.file) and set(EXCL) <= set(M.file) and set(PAIR) <= set(M.file)
grp = {**{f: "S" for f in S}, **{f: "BI" for f in BI}, **{f: "LL" for f in LL}, **{f: "BR" for f in BR}}; assert len(grp) == len(S) + len(BI) + len(LL) + len(BR)
M["group"] = M.file.map(grp).fillna("U")
M["prediction"] = M.group.map({"S": "accessibility", "BI": "initiation", "LL": "initiation", "BR": "initiation", "U": "none"})
M["discovery_trait"] = M.file.map(D).fillna(""); M["status"] = np.where(M.discovery_trait != "", "discovery", "new")
M["pair"] = M.file.map(PAIR).fillna(""); M["excluded"] = M.file.map(EXCL).fillna(""); M["include"] = M.excluded == ""
M["disease"] = M.index < 53   # 표의 앞 53행이 질병, 뒤 54행이 정량 형질
M["id"] = [f"T{i + 1:03d}" for i in range(len(M))]
M.loc[~M.include, ["group", "prediction"]] = ["excluded", "none"]
M[["id", "file", "trait", "disease", "mean_N", "h2_z", "group", "prediction", "status", "discovery_trait", "pair", "include", "excluded", "link"]].to_csv(EXP/"traits_exp3u.csv", index=False, encoding="utf-8-sig")
I = M[M.include]; print(f"포함 {len(I)}개, 제외 {len(M) - len(I)}개"); print(I.groupby(["group", "status"]).size().unstack(fill_value=0).to_string())
# 입력 점검: 열, 행 수, N, 평균 카이제곱, 정렬 블록 열쇠(3j)와의 겹침
KEY = set(pd.read_parquet(J/"snp_key.parquet").SNP); rows = []
for r in I.itertuples():
    p = RAW/f"{r.file}.sumstats.gz"; assert p.exists(), p
    d = pd.read_csv(p, sep=r"\s+"); assert {"SNP", "A1", "A2", "N", "Z"} <= set(d.columns), (r.file, list(d.columns))
    d = d.dropna(subset=["Z", "N"]); rows.append(dict(id=r.id, file=r.file, columns=" ".join(d.columns), n_snp=len(d), dup=int(d.SNP.duplicated().sum()), n_median=float(d.N.median()), n_max=float(d.N.max()),
                                                    mean_chi2=float((d.Z ** 2).mean()), frac_in_key=float(d.SNP.isin(KEY).mean()), has_alleles=bool(d.A1.notna().all() and d.A2.notna().all())))
Q = pd.DataFrame(rows); Q.to_csv(RES/"u00_inputs.csv", index=False)
print(Q.describe().loc[["min", "50%", "max"], ["n_snp", "n_median", "mean_chi2", "frac_in_key"]].to_string()); print("중복 SNP가 있는 파일:", int((Q.dup > 0).sum()), "| 대립유전자 빠진 파일:", int((~Q.has_alleles).sum()))
print(Q.sort_values("n_snp").head(5)[["file", "n_snp", "n_median", "mean_chi2", "frac_in_key"]].to_string(index=False))
