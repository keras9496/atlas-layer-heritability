%title Sequence-based variant-effect predictions reveal trait-specific regulatory-layer profiles of complex-trait heritability
%authors ByungHyug Choi^1,2^ and Jaeseung Jeong^1,*^
%corr ^*^ Corresponding author: Jaeseung Jeong (jsjeong@kaist.ac.kr)
%contents Supplementary Tables S1–S19; Supplementary Figures S1–S5. Reference numbers refer to the reference list of the main text.

# Supplementary Tables

::table Supplementary Table S1|GWAS summary statistics and lineages.|landscape
*N* is the median per-SNP sample size after processing (cases plus controls for case–control traits). *h*^2^~g~ (observed scale), its *Z* score, mean *χ*^2^, the LD score regression intercept and the number of regression SNPs are from the main layer-map model (baselineLD v2.2 and the three all-track layer annotations of the trait's lineage); *h*^2^~g~ varies little between models. Intercepts below 1 for height, BMI, LDL cholesterol, educational attainment and neuroticism suggest that genomic control had been applied to these summary statistics; this lowers *h*^2^~g~ estimates but has little effect on τ\*, which is scaled by *h*^2^~g~. Lineage τ\* is from the single-lineage model (Fig. 2) for the lineage used in the layer analyses, with its rank among 12 lineages. For the validation traits, the lineage was predicted before analysis.
| Trait | Role | Study | Source | *N* | *h*^2^~g~ (s.e.) | *Z* | Mean *χ*^2^ | Intercept (s.e.) | Regression SNPs | Lineage used | Lineage τ\* (rank) |
| Alzheimer's disease | Map | Bellenguez 2022 [@bellenguez2022] | GCST90027158 (stage I) | 487,511 | 0.026 (0.002) | 10.7 | 1.265 | 1.053 (0.012) | 1,162,232 | Myeloid | +2.98 (1) |
| Educational attainment | Map | Okbay 2016 [@okbay2016ea] | Price laboratory | 328,917 | 0.110 (0.004) | 26.8 | 1.645 | 0.966 (0.011) | 1,167,140 | Brain tissue | +0.87 (1) |
| Schizophrenia | Map | PGC 2014 [@pgc2014] | Price laboratory | 70,100 | 0.458 (0.026) | 17.9 | 1.697 | 1.061 (0.018) | 1,082,608 | Brain tissue | +1.55 (1) |
| Bipolar disorder | Map | PGC 2011 [@pgc2011] | Price laboratory | 16,731 | 0.392 (0.069) | 5.7 | 1.180 | 1.031 (0.014) | 750,389 | Brain tissue | +1.91 (2) |
| Neuroticism | Map | Okbay 2016 [@okbay2016neuro] | Price laboratory | 170,911 | 0.096 (0.015) | 6.6 | 1.307 | 0.991 (0.022) | 1,114,925 | Brain tissue | +1.42 (1) |
| BMI | Map | Speliotes 2010 [@speliotes2010] | Price laboratory | 123,822 | 0.146 (0.009) | 16.4 | 1.108 | 0.768 (0.011) | 1,026,547 | Brain tissue | +0.53 (1) |
| Monocyte count | Map | Vuckovic 2020 [@vuckovic2020] | GCST90002393 | 408,112 | 0.244 (0.017) | 14.8 | 2.943 | 1.277 (0.043) | 1,156,723 | Myeloid | +2.82 (1) |
| IBD | Map | Jostins 2012 [@jostins2012] | Price laboratory | 34,652 | 0.358 (0.035) | 10.2 | 1.270 | 1.058 (0.011) | 1,077,495 | Lymphoid | +3.43 (1) |
| Rheumatoid arthritis | Map | Okada 2014 [@okada2014] | Price laboratory | 38,242 | 0.183 (0.025) | 7.2 | 1.146 | 1.028 (0.011) | 1,124,653 | Lymphoid | +2.48 (1) |
| LDL cholesterol | Map | Teslovich 2010 [@teslovich2010] | Price laboratory | 95,386 | 0.102 (0.012) | 8.5 | 1.103 | 0.932 (0.012) | 1,017,468 | Liver | +3.21 (1) |
| Height | Map | Lango Allen 2010 [@langoallen2010] | Price laboratory | 133,719 | 0.221 (0.013) | 17.7 | 1.330 | 0.786 (0.015) | 1,025,325 | Stromal | +1.54 (3) |
| Bone mineral density (heel) | Validation | Morris 2019 [@morris2019] | GCST006979 | 426,824 | 0.284 (0.017) | 16.8 | 3.577 | 1.204 (0.042) | 1,122,478 | Stromal | +1.85 (3) |
| Neutrophil count | Validation | Vuckovic 2020 [@vuckovic2020] | GCST90002398 | 408,112 | 0.197 (0.010) | 20.3 | 2.581 | 1.111 (0.026) | 1,164,547 | Myeloid | +2.04 (1) |
| Total cholesterol | Validation | Zoodsma 2025 [@zoodsma2025] | GCST90500025 | 434,646 | 0.083 (0.008) | 10.5 | 1.710 | 1.073 (0.017) | 1,160,758 | Liver | +1.97 (1) |
::end

::table Supplementary Table S2|AlphaGenome Atlas tracks per lineage and layer.
Number of tracks and of biosamples (ontology terms with any layer) in each lineage after excluding cell lines (372 tracks) and unassigned biosamples (64 tracks). CAGE tracks are strand specific, two per CAGE library (hCAGE or LQhCAGE); some biosamples have both library types, so that, for example, the 36 myeloid CAGE tracks come from 18 libraries of 13 biosamples.
| Lineage | Biosamples | CAGE | DNase | H3K27ac | H3K4me1 | H3K4me3 | RNA-seq | Examples |
| Brain tissue | 38 | 52 | 12 | 8 | 10 | 11 | 24 | prefrontal and temporal cortex, caudate, hippocampus, cerebellum, spinal cord |
| Neural cells | 17 | 12 | 9 | 6 | 6 | 9 | 9 | cultured neurons, neural progenitors, astrocytes, oligodendrocyte precursors, Schwann cells |
| Myeloid | 17 | 36 | 4 | 3 | 3 | 3 | 2 | CD14^+^ monocytes, macrophages, dendritic cells, granulocytes, mast cells, myeloid progenitors |
| Lymphoid | 32 | 18 | 27 | 13 | 12 | 13 | 30 | T, B and NK cells, thymus, spleen, tonsil |
| Liver | 4 | 4 | 4 | 3 | 3 | 3 | 7 | hepatocytes, liver tissue |
| Gastrointestinal | 46 | 44 | 15 | 19 | 21 | 22 | 31 | oesophagus, stomach, colon, small intestine, pancreas |
| Cardiovascular | 59 | 58 | 29 | 13 | 11 | 14 | 49 | heart, aorta, endothelial cells, vascular smooth muscle |
| Skeletal muscle | 21 | 22 | 12 | 7 | 8 | 9 | 17 | skeletal muscle tissue, myoblasts, myotubes |
| Adipose | 17 | 24 | 3 | 1 | 2 | 2 | 8 | adipose tissue, adipocytes, preadipocytes |
| Stromal/mesenchymal | 67 | 84 | 23 | 7 | 6 | 12 | 27 | fibroblasts, mesenchymal stem cells, chondrocytes, osteoblasts |
| Reproductive/endocrine | 33 | 46 | 10 | 13 | 11 | 11 | 29 | reproductive organs, placenta, thyroid, adrenal gland |
| Epithelial/barrier | 82 | 86 | 40 | 17 | 14 | 22 | 57 | skin, lung and airway, kidney, breast, prostate |
| Total | 433 | 486 | 188 | 110 | 107 | 131 | 290 | |
::end

::table Supplementary Table S3|Comparison at each trait's own layer.
Own layer: the top layer in the layer map. Atlas alone: τ\* of the predicted own-layer annotation conditional on baselineLD. vs eQTL: difference between the predicted and the GTEx eQTL (MaxCPP) τ\* in the same model (two-sided *P*). vs measured: difference between the predicted and the measured τ\* of the same layer fitted together (two-sided *P*). Conditioned: predicted τ\* after conditioning on measured CAGE, DNase and H3K27ac annotations of the same lineage (one-sided *P*; retention, the conditioned τ\* as a percentage of the τ\* alone). The conditioned model for height was run after the layer map was known (†).
| Trait | Own layer | Atlas alone | vs eQTL (*P*) | vs measured (*P*) | Conditioned (*P*; retention) |
| Alzheimer's disease | CAGE | +3.05 | +1.56 (0.065) | +2.91 (0.004) | +2.17 (0.001; 71%) |
| Monocyte count | CAGE | +3.81 | +3.23 (5 × 10^−8^) | +2.28 (0.001) | +2.49 (1 × 10^−7^; 65%) |
| Educational attainment | CAGE | +1.21 | +0.97 (6 × 10^−4^) | +0.63 (0.035) | +0.84 (2 × 10^−4^; 69%) |
| Schizophrenia | CAGE | +1.51 | +0.95 (0.002) | +0.63 (0.093) | +1.07 (4 × 10^−5^; 71%) |
| Bipolar disorder | CAGE | +4.37 | +3.82 (8 × 10^−5^) | +3.55 (0.003) | +4.02 (8 × 10^−7^; 92%) |
| Neuroticism | CAGE | +0.97 | +0.47 (0.22) | +0.14 (0.77) | +0.52 (0.047; 54%) |
| BMI | CAGE | +0.74 | +0.38 (0.32) | +0.28 (0.53) | +0.49 (0.064; 67%) |
| IBD | CAGE | +2.94 | +2.48 (6 × 10^−4^) | +2.42 (0.010) | +2.16 (2 × 10^−5^; 73%) |
| Rheumatoid arthritis | CAGE | +2.95 | +2.80 (8 × 10^−4^) | +2.20 (0.030) | +2.04 (8 × 10^−4^; 69%) |
| LDL cholesterol | CAGE | +4.33 | +4.21 (1 × 10^−5^) | +5.24 (6 × 10^−4^) | +4.23 (3 × 10^−7^; 98%) |
| Height | DNase | +1.78 | +1.64 (5 × 10^−4^) | +0.92 (0.035) | +1.27 (9 × 10^−4^; 72%)† |
::end

::table Supplementary Table S4|Primary test (specified before the comparison): predicted CAGE conditional on measured annotations.
τ\* of the predicted CAGE annotation of each trait's lineage, alone and conditional on measured CAGE, DNase and H3K27ac annotations of the same lineage (full model). *q*: Benjamini–Hochberg across 11 traits. Retention: the conditioned τ\* as a percentage of the τ\* alone. Sensitivity analyses used three CAGE tracks matched in number to the other layers (track-matched) or measured peaks without scorer-window extension (peaks only).
| Trait | Lineage | Alone | Conditioned (s.e.) | *q* | Retention | Track-matched (*q*) | Peaks only (*q*) |
| Alzheimer's disease | Myeloid | +3.05 | +2.17 (0.70) | 0.001 | 71% | +2.69 (< 0.001) | +1.96 (0.003) |
| Monocyte count | Myeloid | +3.81 | +2.49 (0.49) | 2 × 10^−6^ | 65% | +2.37 (< 0.001) | +2.69 (< 0.001) |
| Educational attainment | Brain tissue | +1.21 | +0.84 (0.24) | 4 × 10^−4^ | 69% | +0.80 (< 0.001) | +0.93 (< 0.001) |
| Schizophrenia | Brain tissue | +1.51 | +1.07 (0.27) | 1 × 10^−4^ | 71% | +1.15 (< 0.001) | +1.15 (< 0.001) |
| Bipolar disorder | Brain tissue | +4.37 | +4.02 (0.84) | 3 × 10^−6^ | 92% | +3.70 (< 0.001) | +4.17 (< 0.001) |
| Neuroticism | Brain tissue | +0.97 | +0.52 (0.31) | 0.057 | 54% | +0.22 (0.24) | +0.62 (0.030) |
| BMI | Brain tissue | +0.74 | +0.49 (0.32) | 0.071 | 67% | +0.48 (0.052) | +0.57 (0.044) |
| IBD | Lymphoid | +2.94 | +2.16 (0.53) | 7 × 10^−5^ | 73% | +1.64 (0.002) | +2.11 (< 0.001) |
| Rheumatoid arthritis | Lymphoid | +2.95 | +2.04 (0.64) | 0.001 | 69% | +2.01 (0.003) | +1.86 (0.003) |
| LDL cholesterol | Liver | +4.33 | +4.23 (0.85) | 2 × 10^−6^ | 98% | +4.54 (< 0.001) | +4.20 (< 0.001) |
| Height | Stromal | +0.42 | +0.07 (0.33) | 0.42 | 16% | +0.83 (0.015) | +0.09 (0.40) |
::end

::table Supplementary Table S5|Leave-one-trait-out heterogeneity of layer profiles (post hoc).
Multivariate Cochran's *Q* on the (CAGE − DNase, CAGE − H3K27ac) contrasts of the layer map with one trait removed (18 degrees of freedom; 20 with all 11 traits), computed with the joint jackknife across traits; *P* from the χ^2^ distribution and from the *F* approximation (Methods). Values with traits treated as independent are in Supplementary Table S10.
| Trait removed | *Q* | *P* | *P* (*F* approximation) |
| None (all 11 traits) | 74.3 | 3.6 × 10^−8^ | 7.6 × 10^−6^ |
| Alzheimer's disease | 72.3 | 1.9 × 10^−8^ | 3.2 × 10^−6^ |
| Monocyte count | 56.3 | 8.0 × 10^−6^ | 1.9 × 10^−4^ |
| Educational attainment | 73.7 | 1.0 × 10^−8^ | 2.2 × 10^−6^ |
| Schizophrenia | 73.9 | 9.9 × 10^−9^ | 2.1 × 10^−6^ |
| Bipolar disorder | 64.5 | 3.8 × 10^−7^ | 2.4 × 10^−5^ |
| Neuroticism | 65.9 | 2.2 × 10^−7^ | 1.7 × 10^−5^ |
| BMI | 71.1 | 3.0 × 10^−8^ | 4.4 × 10^−6^ |
| IBD | 74.3 | 8.5 × 10^−9^ | 1.9 × 10^−6^ |
| Rheumatoid arthritis | 70.0 | 4.5 × 10^−8^ | 5.7 × 10^−6^ |
| LDL cholesterol | 65.7 | 2.4 × 10^−7^ | 1.8 × 10^−5^ |
| Height | 38.5 | 0.0033 | 0.014 |
::end

::table Supplementary Table S6|Analysis sequence and status.
Each stage was specified in a dated protocol before it was run. Because the stages were sequential, each protocol was written with knowledge of the results listed under Results already seen. Prospective: predictions fixed before any data for those traits were analysed. Post hoc: analyses specified after the results they follow up, without a criterion fixed in advance.
| Stage | Protocol date | Content | Results already seen | Status |
| Initial AD analysis | 2026-09-25 | Monocyte/macrophage and other cell-type annotations; thresholds 1%, 5%, 10% and continuous; AVI conditioning; *APOE* region included | Locus-level AD study [@choi2026] | Protocol before analysis |
| AD robustness and control traits | 2026-09-28 | Control traits for the myeloid and neuron annotations | Initial AD analysis | Protocol before analysis |
| | 2026-09-28 | Removal of genome-wide significant loci; leave-one-block-out | Initial AD analysis | Post hoc |
| Lineage map | 2026-09-28 | 12 lineages, 11 traits; AD myeloid, brain traits in brain tissue, AD brain tissue | Initial AD analysis; control traits | Protocol before analysis |
| Single-layer map | 2026-10-02 | 72 lineage–layer annotations; joint six-layer models | Lineage map | Protocol before analysis |
| eQTL comparison | 2026-10-04 | Predicted CAGE versus GTEx MaxCPP in the same model | Lineage and single-layer maps | Post hoc (descriptive) |
| Measured annotations | 2026-10-08 | FANTOM5 and ENCODE annotations; predicted CAGE conditional on three measured layers | Single-layer map; eQTL comparison | Protocol before analysis |
| Layer map with resampling | 2026-10-09 | Three-layer joint models, 20 resamples; heterogeneity, pairwise contrasts, classification | Single-layer map and six-layer models; eQTL and measured-annotation comparisons | Protocol before analysis (exploratory classification) |
| Validation predictions | 2026-10-09 (traits, rule); 2026-10-10 (predictions recorded) | Three new traits; predicted lineage and layer | As above and, when the predictions were recorded, the layer map; no data for the new traits | Prospective |
| Leave-one-trait-out heterogeneity | 2026-10-10 | Supplementary Table S5 | Layer map | Post hoc |
| Own-layer comparison | 2026-10-11 | Supplementary Table S3; one new conditional model for height | Layer map; measured-annotation comparison | Post hoc |
| Threshold sensitivity | 2026-10-11 | Top 1%, 0.5% and 0.1% layer maps; criterion at 0.1% fixed before analysis (Supplementary Table S7) | All of the above | Protocol before analysis |
| Cell-type-matched layer map | 2026-10-13 | One unit per layer from matched cell types or regions (Supplementary Table S8) | All of the above | Protocol before analysis |
| Comparison within shared lineages | 2026-10-13 | All-track joint models in each enriched lineage; skeletal versus lipid and blood traits (Supplementary Table S9) | All of the above | Protocol before analysis |
| Joint jackknife across traits | 2026-10-14 | All layer models refitted with genomic blocks aligned across traits; tests between traits with their covariance (Supplementary Table S10) | All of the above | Protocol before analysis |
| Proximal and distal annotations | 2026-10-15 | Main-map annotations split at 1 kb from the nearest annotated TSS; proximal and distal parts fitted jointly (Supplementary Table S11) | All of the above, including the joint jackknife | Protocol before analysis |
| Diagnostics and influential regions | 2026-10-17 | Threshold contrasts, collinearity and resampling convergence from existing outputs (Supplementary Fig. S4 and Supplementary Table S13); refitting without the most influential blocks (Supplementary Table S14) | All of the above | Protocol before analysis (block removal); descriptive (others) |
| Activity conditioning | 2026-10-18 | Predicted-activity annotations from the same tracks and windows; effect + activity models (Supplementary Table S15) | All of the above, including the diagnostics | Protocol before analysis |
| Larger GWAS | 2026-10-21 | Recent GWAS for height, BMI, schizophrenia, bipolar disorder and LDL cholesterol; layer, lineage and activity models with lineages, annotations and tests unchanged (Supplementary Table S16) | All of the above | Protocol before analysis |
| Continuous and non-overlapping annotations | 2026-10-21 | Continuous percentile annotations and four non-overlapping bins of the three layers (Supplementary Table S17, Supplementary Fig. S4b) | All of the above, including the larger-GWAS results | Protocol before analysis |
| Activity constructions | 2026-10-23 | Continuous activity covariates and activity-stratified effect annotations; criterion that the binary pattern is reproduced in both (Supplementary Table S18) | All of the above, including the activity conditioning | Protocol before analysis |
| Effects within predicted active elements | 2026-10-23 | Activity-matched larger- and smaller-effect halves within the top 10% (and 20%) of predicted activity (Supplementary Table S18) | Activity constructions | Post hoc |
| Larger AD GWAS | 2026-10-25 | Layer, lineage and activity models with the 2026 AD meta-analysis without proxy cases (Supplementary Tables S16 and S18) | All of the above | Protocol before analysis |
| Cell-type eQTL conditioning | 2026-10-25 | FunGen-xQTL MaxCPP annotations for five cell-type groups added to the effect and effect + activity models (Supplementary Table S19) | All of the above | Protocol before analysis |
::end

::table Supplementary Table S7|Layer map at stricter annotation thresholds.
Top layer and its τ\* in the three-layer joint model (all tracks of the trait's lineage, single run) at each threshold. Asterisk: the top layer exceeds the second with one-sided *P* < 0.05. The 5% column is the all-track version of the main layer map. By the criterion fixed before the analysis (at 0.1%, heterogeneity *P* < 0.05 and a top layer matching the main map in at least 9 of 11 map traits, with DNase for height), the map was not reproduced at 0.1%; at 1% the ranking part was met but heterogeneity was not significant. The summary rows give the tests of the main layer map repeated at each threshold, with the joint jackknife across traits: the pairwise contrasts are two-sided with α = 0.017, and the validation contrast is the one-sided difference in CAGE − DNase between the traits predicted to be CAGE type and DNase type. Brain traits: educational attainment, schizophrenia, bipolar disorder and neuroticism.
| Trait | Top 0.1% | Top 0.5% | Top 1% | Top 5% |
| Educational attainment | DNase +0.58 | CAGE +0.79 | DNase +0.73 | CAGE +1.22\* |
| Schizophrenia | CAGE +0.31 | CAGE +0.54 | CAGE +0.39 | CAGE +1.53\* |
| Bipolar disorder | H3K27ac +1.89 | CAGE +2.52 | CAGE +2.84 | CAGE +4.79\* |
| Neuroticism | DNase +0.62\* | DNase +0.47 | DNase +0.59 | CAGE +0.91 |
| BMI | CAGE +0.36 | DNase +0.11 | CAGE +0.57 | CAGE +0.88 |
| Alzheimer's disease | CAGE +3.49 | CAGE +2.26 | CAGE +2.39 | CAGE +2.91\* |
| Monocyte count | H3K27ac +1.61 | H3K27ac +2.36 | CAGE +3.11 | CAGE +3.47\* |
| IBD | CAGE +2.89 | CAGE +2.72 | CAGE +2.13 | CAGE +2.32 |
| Rheumatoid arthritis | DNase +1.95 | CAGE +2.90 | CAGE +2.36 | CAGE +2.98\* |
| LDL cholesterol | DNase +3.67 | CAGE +2.74 | CAGE +3.57 | CAGE +3.91\* |
| Height | DNase +0.85 | DNase +0.86 | DNase +1.47 | DNase +1.65\* |
| Bone mineral density | DNase +0.96 | DNase +0.98 | DNase +1.47 | DNase +1.51 |
| Neutrophil count | CAGE +1.95\* | CAGE +1.74 | CAGE +1.61 | CAGE +1.35 |
| Total cholesterol | DNase +1.57 | CAGE +2.66\* | CAGE +2.82 | CAGE +3.70\* |
| Summary across the 11 map traits | | | | |
| Top layer matching the main map | 5 of 11 | 8 of 11 | 9 of 11 | 11 of 11 |
| Traits with top layer exceeding the second | 1 | 0 | 0 | 8 |
| Heterogeneity *Q* (20 d.f.), *P* | 48.1, 4 × 10^−4^ | 40.0, 0.005 | 31.0, 0.055 | 70.1, 2 × 10^−7^ |
| AD vs height, CAGE − DNase, *P* | 0.028 | 0.13 | 0.057 | 0.008 |
| Brain traits vs height, CAGE − DNase, *P* | 0.71 | 0.26 | 0.025 | 8 × 10^−6^ |
| IBD and RA vs AD and monocyte count, CAGE − H3K27ac, *P* | 0.10 | 0.005 | 0.38 | 0.76 |
| Validation contrast, one-sided *P* | 0.018 | 2 × 10^−4^ | 3 × 10^−5^ | 4 × 10^−6^ |
::end

::table Supplementary Table S8|Layer map with cell-type-matched layers.
Three-layer joint models (τ\*, standard error) with annotations built from matched cell types or anatomical regions, one unit per layer from each matched group (single run). Matched groups: myeloid, CD14^+^ monocytes (CAGE, CD14^+^CD16^−^ classical monocytes); brain tissue, spinal cord, hippocampus, caudate nucleus and dorsolateral prefrontal cortex (CAGE, middle frontal gyrus); stromal, dermal fibroblasts, lung fibroblasts and mesenchymal stem cells; liver, hepatocytes and liver tissue; lymphoid, natural killer cells, naive CD4^+^ T cells, regulatory T cells, spleen and thymus. Units per layer: myeloid 1, brain tissue 4, stromal 3, liver 2, lymphoid 5. Asterisk: the top layer exceeds the second with one-sided *P* < 0.05. Main map: top layer in the main layer map (validation traits: predicted layer). By the criterion fixed before the analysis (an unchanged top layer in at least 7 of the 8 traits with a confirmed type, including DNase for height and CAGE for AD), the map was retained (8 of 8). For height, the stromal CAGE − DNase contrast with cell-type-matched layers was −1.92 ± 0.74 (one-sided *P* = 0.005; Table 1).
| Trait | Lineage | CAGE | DNase | H3K27ac | Top layer | Main map |
| Alzheimer's disease | Myeloid | +3.52 (0.84) | +0.83 (0.51) | −0.37 (0.64) | CAGE\* | CAGE |
| Monocyte count | Myeloid | +3.78 (0.63) | +1.20 (0.40) | +0.07 (0.41) | CAGE\* | CAGE |
| Educational attainment | Brain tissue | +1.28 (0.25) | +0.10 (0.25) | −0.08 (0.23) | CAGE\* | CAGE |
| Schizophrenia | Brain tissue | +1.66 (0.28) | −0.05 (0.26) | −0.18 (0.24) | CAGE\* | CAGE |
| Bipolar disorder | Brain tissue | +4.93 (0.89) | +2.25 (0.88) | −2.72 (0.88) | CAGE\* | CAGE |
| Neuroticism | Brain tissue | +1.06 (0.34) | −0.15 (0.41) | +0.12 (0.33) | CAGE\* | CAGE |
| BMI | Brain tissue | +1.09 (0.33) | −0.05 (0.34) | −0.48 (0.34) | CAGE\* | CAGE |
| IBD | Lymphoid | +2.35 (0.54) | +0.51 (0.67) | +1.29 (0.73) | CAGE | CAGE |
| Rheumatoid arthritis | Lymphoid | +2.95 (0.80) | +0.32 (0.74) | −0.26 (0.84) | CAGE\* | CAGE |
| LDL cholesterol | Liver | +3.93 (0.82) | +0.60 (0.75) | +1.16 (0.71) | CAGE\* | CAGE |
| Height | Stromal | +0.16 (0.43) | +2.08 (0.49) | +0.34 (0.40) | DNase\* | DNase |
| Bone mineral density | Stromal | +0.18 (0.32) | +2.41 (0.43) | +0.54 (0.40) | DNase\* | DNase |
| Neutrophil count | Myeloid | +2.14 (0.41) | +0.83 (0.34) | +0.27 (0.33) | CAGE\* | CAGE |
| Total cholesterol | Liver | +3.75 (0.88) | +0.41 (0.45) | −0.10 (0.73) | CAGE\* | CAGE |
::end

::table Supplementary Table S9|Comparison of skeletal and lipid or blood traits within shared lineages.
All-track three-layer joint models were fitted for each trait in every lineage in which it was enriched in the lineage map (one-sided *P* < 0.05). Shown are the lineages in which both a DNase-type trait and a CAGE-type trait were enriched; within each lineage all traits share the same annotations. Values are inverse-variance-weighted CAGE − DNase contrasts (τ\*, standard error) of the enriched traits of each type. Differences, their standard errors and *P* values use the joint jackknife across traits and lineages; the group values use the standard errors of each trait. With the conservative standard error of the protocol, which assumed perfect correlation between lineages, the pooled difference was +2.29 (0.84), one-sided *P* = 0.003, and by the criterion fixed before the analysis (pooled one-sided *P* < 0.05 with this standard error and a positive difference in at least four of six lineages), the difference was retained. Top layers depended on the lineage: H3K27ac ranked first for height in the adipose, cardiovascular, epithelial/barrier and skeletal muscle lineages, and for LDL or total cholesterol in the gastrointestinal, adipose, cardiovascular, lymphoid and reproductive/endocrine lineages.
| Lineage | CAGE-type traits | CAGE − DNase | DNase-type traits | CAGE − DNase | Difference | One-sided *P* |
| Reproductive/endocrine | LDL, total cholesterol | +2.35 (0.71) | Height, BMD | −1.58 (0.35) | +3.94 (0.87) | 3 × 10^−6^ |
| Liver | LDL, total cholesterol | +3.29 (0.72) | BMD | +0.44 (0.51) | +2.85 (0.90) | 8 × 10^−4^ |
| Gastrointestinal | LDL, total cholesterol | +2.35 (0.85) | BMD | −0.31 (0.56) | +2.66 (1.12) | 0.009 |
| Epithelial/barrier | LDL | +2.51 (1.10) | Height, BMD | −0.08 (0.35) | +2.58 (1.09) | 0.009 |
| Cardiovascular | Total cholesterol | +1.62 (0.91) | Height, BMD | −0.57 (0.37) | +2.19 (0.98) | 0.013 |
| Adipose | Monocyte count, LDL, total cholesterol | +1.10 (0.50) | Height, BMD | +0.32 (0.38) | +0.78 (0.68) | 0.13 |
| Pooled | | | | | +2.29 (0.63) | 1.5 × 10^−4^ |
::end

::table Supplementary Table S10|Tests between traits with traits treated as independent and with the joint jackknife.
Original: jackknife blocks of equal SNP number within each trait, traits treated as independent (the values of earlier versions). Aligned: genomic blocks aligned across traits, traits treated as independent. Joint: aligned blocks with the covariance between traits (the values reported in the text). Contrasts are two-sided except the validation contrast and the shared-lineage difference (one-sided). Brain traits: educational attainment, schizophrenia, bipolar disorder and neuroticism. Error correlations of the CAGE − DNase contrast between the 14 traits had a median of +0.03 (range −0.15 to +0.39; largest for monocyte and neutrophil counts, +0.39, and LDL and total cholesterol, +0.30).
| Analysis | Test | Original | Aligned | Joint |
| Layer map (Fig. 3) | Heterogeneity *Q* (20 d.f.), *P* | 66.2, *P* = 7.4 × 10^−7^ | 62.3, *P* = 3.2 × 10^−6^ | 74.3, *P* = 3.6 × 10^−8^ (*F*: 7.6 × 10^−6^) |
|  | AD vs height, CAGE − DNase | +4.51 ± 1.47, *P* = 0.0022 | ± 1.50, *P* = 0.0026 | ± 1.51, *P* = 0.0028 |
|  | Brain traits vs height, CAGE − DNase | +2.75 ± 0.55, *P* = 6.6 × 10^−7^ | ± 0.58, *P* = 1.7 × 10^−6^ | ± 0.58, *P* = 2.4 × 10^−6^ |
|  | IBD and RA vs AD and monocyte count, CAGE − H3K27ac | −1.88 ± 0.83, *P* = 0.023 | ± 0.85, *P* = 0.026 | ± 0.90, *P* = 0.037 |
| Validation (Fig. 4) | CAGE-type vs DNase-type predictions, CAGE − DNase | +3.42 ± 0.70, *P* = 5.1 × 10^−7^ | ± 0.70, *P* = 4.6 × 10^−7^ | ± 0.67, *P* = 1.6 × 10^−7^ |
| Thresholds (Supplementary Table S7) | Heterogeneity, top 0.1% | 41.4, *P* = 0.0033 | 40.9, *P* = 0.0038 | 48.1, *P* = 4.2 × 10^−4^ (*F*: 0.0039) |
|  | Heterogeneity, top 0.5% | 33.2, *P* = 0.032 | 33.6, *P* = 0.029 | 40.0, *P* = 0.005 (*F*: 0.022) |
|  | Heterogeneity, top 1% | 30.4, *P* = 0.064 | 29.7, *P* = 0.075 | 31.0, *P* = 0.055 (*F*: 0.13) |
|  | Heterogeneity, top 5% (all tracks) | 61.9, *P* = 3.7 × 10^−6^ | 58.6, *P* = 1.2 × 10^−5^ | 70.1, *P* = 1.7 × 10^−7^ (*F*: 2.1 × 10^−5^) |
| Cell-type-matched (Supplementary Table S8) | Heterogeneity | 56.8, *P* = 2.2 × 10^−5^ | 53.4, *P* = 7.0 × 10^−5^ | 56.3, *P* = 2.6 × 10^−5^ (*F*: 5.9 × 10^−4^) |
|  | AD vs height | +4.62 ± 1.29, *P* = 3.3 × 10^−4^ | ± 1.32, *P* = 4.8 × 10^−4^ | ± 1.34, *P* = 5.7 × 10^−4^ |
|  | Brain traits vs height | +3.36 ± 0.77, *P* = 1.4 × 10^−5^ | ± 0.81, *P* = 3.2 × 10^−5^ | ± 0.81, *P* = 3.4 × 10^−5^ |
|  | Validation contrast | +4.17 ± 0.74, *P* = 8.6 × 10^−9^ | ± 0.73, *P* = 5.9 × 10^−9^ | ± 0.76, *P* = 1.9 × 10^−8^ |
| Shared lineages (Supplementary Table S9) | Pooled difference (original: conservative standard error of the protocol) | +2.29 ± 0.84, *P* = 0.0033 | – | ± 0.63, *P* = 1.5 × 10^−4^ |
::end

::table Supplementary Table S11|Layer map at distal and proximal SNPs.
Each all-track 5% layer annotation of the main map was split at 1 kb from the nearest TSS of a GENCODE v46 protein-coding or lncRNA transcript, and the proximal and distal parts of the three layers were fitted jointly (τ\*, standard error; single run per trait in its own lineage). Distal annotations covered 4.0–4.2% (CAGE), 4.6–4.8% (DNase) and 4.4–4.6% (H3K27ac) of common SNPs and proximal annotations 0.2–1.0%. Asterisk: the top layer exceeds the second with one-sided *P* < 0.05. Across the 11 map traits (joint jackknife), the distal layer profiles differed between traits (*Q* = 63.7, 20 d.f., *P* = 1.9 × 10^−6^; *F* approximation 1.0 × 10^−4^), with CAGE − DNase contrasts of AD versus height of +4.41 (1.49), *P* = 0.003, and of the brain traits versus height of +2.59 (0.65), *P* = 7 × 10^−5^, and a validation contrast of +3.51 (0.70), one-sided *P* = 2 × 10^−7^. The proximal profiles did not differ (*Q* = 18.7, *P* = 0.54). By the criterion fixed before the analysis (at distal SNPs, CAGE first for at least 5 of the 7 CAGE-type traits, DNase for height and heterogeneity *P* < 0.05), the layer map was retained at distal SNPs. For height, the stromal CAGE − DNase contrast at distal SNPs was −1.68 ± 0.62 (one-sided *P* = 0.003; Table 1).
| Trait | Lineage | CAGE, distal | DNase, distal | H3K27ac, distal | Top layer, distal | Top layer, proximal |
| Alzheimer's disease | Myeloid | +2.86 (0.86) | +0.14 (0.76) | +0.23 (0.58) | CAGE\* | CAGE |
| Monocyte count | Myeloid | +2.63 (0.53) | −0.05 (0.46) | +1.40 (0.43) | CAGE\* | CAGE\* |
| Educational attainment | Brain tissue | +0.94 (0.23) | +0.29 (0.25) | −0.17 (0.23) | CAGE\* | CAGE |
| Schizophrenia | Brain tissue | +1.44 (0.26) | +0.19 (0.27) | −0.23 (0.26) | CAGE\* | H3K27ac |
| Bipolar disorder | Brain tissue | +4.45 (0.92) | +2.41 (0.87) | −2.64 (0.92) | CAGE\* | CAGE\* |
| Neuroticism | Brain tissue | +0.83 (0.35) | +0.19 (0.37) | +0.07 (0.34) | CAGE | CAGE |
| BMI | Brain tissue | +0.96 (0.37) | +0.21 (0.39) | −0.55 (0.38) | CAGE | CAGE |
| IBD | Lymphoid | +1.78 (0.56) | +1.21 (0.71) | +0.70 (0.65) | CAGE | CAGE |
| Rheumatoid arthritis | Lymphoid | +2.66 (0.80) | +0.99 (0.89) | −0.42 (0.79) | CAGE | CAGE |
| LDL cholesterol | Liver | +4.24 (0.92) | +0.67 (0.84) | +0.62 (0.65) | CAGE\* | DNase |
| Height | Stromal | −0.18 (0.37) | +1.50 (0.44) | +0.30 (0.42) | DNase | DNase |
| Bone mineral density | Stromal | +0.11 (0.30) | +1.46 (0.35) | +1.04 (0.33) | DNase | DNase |
| Neutrophil count | Myeloid | +1.18 (0.45) | +0.16 (0.42) | +1.08 (0.35) | CAGE | DNase |
| Total cholesterol | Liver | +3.61 (0.78) | −0.06 (0.51) | +0.11 (0.66) | CAGE\* | DNase |
::end

::table Supplementary Table S12|Robustness of the myeloid enrichment of Alzheimer's disease heritability.
τ\* with standard errors and one-sided *P* values for τ\* > 0, all conditional on baselineLD v2.2 and with the *APOE* region excluded unless stated. Rows marked † were specified after the main results were seen. The removal of loci and the leave-one-block-out analysis apply to the initial myeloid lineage annotation, not to the layer contrast. The last three rows give the AD results of the block-removal and activity analyses, which used the all-track myeloid layer annotations (Supplementary Tables S14 and S15).
| Analysis | Annotation | τ\* (s.e.) | *P* |
| Lineage map (Fig. 2) | Myeloid lineage, top 5% | +2.98 (0.68) | 6.4 × 10^−6^ |
| Lineage-specific model | Myeloid lineage, conditional on all-lineage annotation | +3.82 (0.92) | 1.6 × 10^−5^ |
| Single layer | Myeloid CAGE | +3.05 (0.72) | 1.0 × 10^−5^ |
| Layer map (Fig. 3) | Myeloid CAGE, joint with DNase and H3K27ac, 20 resamples | +3.22 (0.85) | first in 20/20 |
| Measured conditioning (Fig. 5) | Myeloid CAGE, conditional on measured CAGE, DNase and H3K27ac | +2.17 (0.70) | 9.5 × 10^−4^ |
| Initial analysis | Monocyte and macrophage tracks, top 5% | +2.75 (0.65) | 1.3 × 10^−5^ |
|   Threshold | top 1% | +3.13 (0.91) | 2.8 × 10^−4^ |
|   Threshold | top 10% | +2.33 (0.54) | 8.3 × 10^−6^ |
|   Threshold | continuous score | +1.39 (0.41) | 3.5 × 10^−4^ |
|   Conditional on overall variant impact | top 5%, conditional on AVI top 5% | +2.67 (0.67) | 3.6 × 10^−5^ |
|   *APOE* region included | top 5% | +2.83 (0.67) | 1.2 × 10^−5^ |
|   Loci removed† | top 5%, ±1 Mb of genome-wide significant loci removed (162 Mb) | +1.88 (0.63) | 1.5 × 10^−3^ |
|   Leave one block out† | top 5%, minimum over 200 blocks | +2.46 | all blocks > 0 |
| Lineage map (Fig. 2) | Brain-tissue lineage, top 5% | +0.88 (0.59) | 0.068 |
| Influential blocks removed | Myeloid CAGE, joint with DNase and H3K27ac, five blocks most influential for the AD–height contrast removed | +1.70 (0.65) | 0.005 |
| Activity conditioning (Fig. 5) | Myeloid CAGE allelic effect, conditional on predicted activity | +0.42 (0.94) | 0.33 |
|   Same model | Myeloid CAGE predicted activity | +1.91 (0.67) | 0.002 |
::end

::table Supplementary Table S13|Overlap and collinearity of the three layer annotations.
For each analysed lineage: overlap of the CAGE (C), DNase (D) and H3K27ac (H) annotations across target SNPs (Jaccard index); for regression SNPs, the correlation of their LD scores after regressing out the 97 baseline LD scores (partial *r*), the variance inflation factor of each layer given the other two after the baseline (VIF) and the condition number of the standardized design matrix (118.2 for the baseline alone). All tracks: annotations of the all-track models; matched: cell-type-matched annotations (Supplementary Table S8). The flags fixed before the analysis (VIF > 5 or |*r*| > 0.9) were met in no lineage. Raw correlations of the layer LD scores were 0.54–0.92, mostly because of shared baseline structure (baseline *R*^2^ 0.86–0.95). Of the 33 joint-model τ\* of the map traits (all tracks), five changed sign relative to the single-layer models, all to values near zero, and the only τ\* more than 1.96 standard errors below zero was H3K27ac for bipolar disorder (−2.73 ± 0.93; −0.34 in the single-layer model; −2.72 ± 0.88 with cell-type-matched annotations), in brain tissue, where the VIF between layers was 1.13–1.27; it was absent with a larger bipolar disorder GWAS (+0.16 ± 0.27; Supplementary Table S16).
| Lineage | Annotations | Jaccard C–D, C–H, D–H | Partial *r* C–D, C–H, D–H | VIF C, D, H | Condition number |
| Myeloid | All tracks | 0.30, 0.25, 0.31 | 0.55, 0.41, 0.55 | 1.46, 1.74, 1.47 | 121.7 |
| Brain tissue | All tracks | 0.15, 0.21, 0.22 | 0.21, 0.31, 0.39 | 1.12, 1.19, 1.27 | 122.1 |
| Lymphoid | All tracks | 0.15, 0.28, 0.21 | 0.27, 0.46, 0.41 | 1.28, 1.21, 1.42 | 122.2 |
| Liver | All tracks | 0.10, 0.20, 0.19 | 0.08, 0.37, 0.34 | 1.16, 1.13, 1.31 | 123.6 |
| Stromal | All tracks | 0.14, 0.18, 0.36 | 0.24, 0.31, 0.57 | 1.12, 1.50, 1.56 | 121.6 |
| Myeloid | Matched | 0.24, 0.19, 0.25 | 0.42, 0.33, 0.50 | 1.24, 1.49, 1.37 | 121.8 |
| Brain tissue | Matched | 0.12, 0.22, 0.18 | 0.13, 0.34, 0.34 | 1.13, 1.13, 1.26 | 122.4 |
| Lymphoid | Matched | 0.12, 0.28, 0.16 | 0.20, 0.49, 0.38 | 1.31, 1.17, 1.47 | 121.9 |
| Liver | Matched | 0.08, 0.20, 0.18 | 0.03, 0.37, 0.25 | 1.16, 1.07, 1.24 | 124.8 |
| Stromal | Matched | 0.15, 0.17, 0.38 | 0.33, 0.32, 0.59 | 1.15, 1.59, 1.57 | 122.3 |
::end

::table Supplementary Table S14|Layer contrasts without the most influential genomic blocks (all-track models).
For each test, the 200 aligned jackknife blocks were ranked by their influence on the test (the full estimate minus the estimate without the block; for *Q*, its decrease), and the regression SNPs of the 5 (primary) or 10 most influential blocks were removed before refitting. Criterion fixed before the refitting (5 blocks): same sign and two-sided *P* < 0.05, not dependent on a few blocks; same sign and *P* ≥ 0.05, weakened; sign change, dependent; for *Q*, *P* < 0.05, retained. The lower rows decompose the AD–height contrast for the blocks removed for that test, with single-trait jackknife standard errors (two-sided *P* for contrasts, one-sided for τ\*). The five blocks most influential for the AD–height contrast were chr11:47.2–69.1 Mb (containing the *MS4A* locus), chr19:46.9 Mb–chr20:0.6 Mb (*LILRB5*, *SIGLEC11*), chr14:86.5–97.8 Mb (*RIN3*), chr15:95.3 Mb–chr16:5.6 Mb (a height signal) and chr10:6.2–15.6 Mb (*USP6NL*); those for the brain-trait contrast and for *Q* lay mostly in regions with height signals. With all blocks, these all-track values differ slightly from the resample averages of the main layer map. Brain traits: educational attainment, schizophrenia, bipolar disorder and neuroticism.
| Test | All blocks (*P*) | 5 blocks removed (*P*) | 10 blocks removed (*P*) | Verdict |
| AD versus height, CAGE − DNase | +4.41 ± 1.65 (0.008) | +2.21 ± 1.20 (0.064) | +1.23 ± 1.16 (0.29) | Weakened |
| Brain traits versus height, CAGE − DNase | +2.81 ± 0.63 (8.3 × 10^−6^) | +2.06 ± 0.56 (2.2 × 10^−4^) | +1.55 ± 0.53 (0.004) | Not dependent on a few blocks |
| Heterogeneity *Q* (11 map traits, 20 degrees of freedom) | 70.1 (1.7 × 10^−7^) | 69.3 (2.4 × 10^−7^) | 62.3 (3.1 × 10^−6^) | Retained |
|   AD, CAGE − DNase | +2.70 ± 1.57 (0.087) | +0.75 ± 1.12 (0.50) | +0.10 ± 1.11 (0.93) | |
|   AD, myeloid CAGE τ\* | +2.91 ± 0.92 (7.5 × 10^−4^) | +1.70 ± 0.65 (0.005) | +1.35 ± 0.64 (0.018) | |
|   AD, myeloid DNase τ\* | +0.22 ± 0.83 (0.40) | +0.95 ± 0.69 (0.085) | +1.25 ± 0.71 (0.038) | |
|   Height, CAGE − DNase | −1.71 ± 0.60 (0.004) | −1.46 ± 0.58 (0.012) | −1.13 ± 0.53 (0.034) | |
|   Height, stromal DNase τ\* | +1.65 ± 0.44 (9.2 × 10^−5^) | +1.55 ± 0.45 (2.5 × 10^−4^) | +1.25 ± 0.40 (8.2 × 10^−4^) | |
::end

::table Supplementary Table S15|Predicted allelic effects conditional on predicted activity.|landscape
All tracks of each trait's lineage, own layer from the layer map. Effect only: τ\* of the own-layer allelic-effect annotation in the joint model of the three effect layers. Effect + activity: the same τ\* after adding the three predicted-activity annotations of the same tracks and windows (primary test; one-sided, *q* by Benjamini–Hochberg across the 11 traits). Retention: effect + activity τ\* as a percentage of effect-only τ\*. + measured: the effect τ\* with the measured CAGE, DNase and H3K27ac annotations also in the model (*q* across the 11 traits). Activity: τ\* (s.e.) of the three activity annotations in the effect + activity model. Criterion fixed before the analysis: *q* < 0.05 in at least 6 of 11 traits, retained; 3–5, partly explained by activity; 2 or fewer, largely explained by activity. Result: 5 of 11 (partly explained by activity). With activity in the model, the effect contrasts still differed between traits (*Q* = 53.1, 20 degrees of freedom, *P* = 7.9 × 10^−5^); the three pairwise contrasts of the layer map were +2.01 ± 1.63 for AD versus height (*P* = 0.22), +2.63 ± 0.64 for the four brain traits versus height (*P* = 4.5 × 10^−5^) and +0.30 ± 1.21 for IBD and rheumatoid arthritis versus AD and monocyte count (CAGE − H3K27ac, *P* = 0.81), and the validation contrast was +1.48 ± 0.80 (one-sided *P* = 0.032). The activity annotations alone also differed between traits (*Q* = 41.1, *P* = 0.0036). Conditioning with continuous activity covariates and with activity-stratified effect annotations, and the comparison within predicted active elements, are in Supplementary Table S18. Effect and activity annotations of the same layer overlapped little (Jaccard index 0.05–0.17), and the variance inflation factors of the six LD scores after the baseline were 1.11–2.07.
| Trait | Lineage | Own layer | Effect only | Effect + activity (s.e.) | *q* | Retention | + measured (*q*) | Activity CAGE | Activity DNase | Activity H3K27ac |
| Educational attainment | Brain tissue | CAGE | +1.22 | +0.55 (0.26) | 0.038 | 45% | +0.47 (0.079) | +0.38 (0.10) | +0.91 (0.18) | −0.13 (0.10) |
| Schizophrenia | Brain tissue | CAGE | +1.53 | +0.86 (0.27) | 0.004 | 56% | +0.75 (0.010) | +0.57 (0.12) | +0.94 (0.24) | −0.34 (0.14) |
| Bipolar disorder | Brain tissue | CAGE | +4.79 | +4.89 (1.01) | 7.4 × 10^−6^ | 102% | +4.69 (1.3 × 10^−5^) | +0.21 (0.33) | +0.88 (0.74) | −0.98 (0.46) |
| Neuroticism | Brain tissue | CAGE | +0.91 | +0.34 (0.36) | 0.28 | 37% | +0.21 (0.40) | +0.29 (0.12) | +1.09 (0.35) | −0.24 (0.17) |
| BMI | Brain tissue | CAGE | +0.88 | +0.35 (0.38) | 0.28 | 40% | +0.33 (0.36) | +0.31 (0.16) | +0.43 (0.30) | +0.05 (0.18) |
| Alzheimer's disease | Myeloid | CAGE | +2.91 | +0.42 (0.94) | 0.40 | 14% | +0.52 (0.40) | +1.91 (0.67) | +0.36 (0.76) | +0.36 (0.49) |
| Monocyte count | Myeloid | CAGE | +3.47 | +0.02 (0.73) | 0.49 | 1% | +0.04 (0.49) | +0.50 (0.37) | +1.60 (0.80) | +2.27 (0.42) |
| IBD | Lymphoid | CAGE | +2.32 | +0.08 (0.72) | 0.49 | 3% | +0.01 (0.49) | +0.53 (0.36) | +2.47 (0.81) | +0.15 (0.57) |
| Rheumatoid arthritis | Lymphoid | CAGE | +2.98 | +0.40 (0.90) | 0.40 | 14% | +0.35 (0.42) | +0.52 (0.42) | +1.66 (0.77) | +1.17 (0.58) |
| LDL cholesterol | Liver | CAGE | +3.91 | +2.46 (1.06) | 0.028 | 63% | +2.51 (0.024) | +1.45 (0.59) | +0.98 (0.78) | +0.10 (0.60) |
| Height | Stromal | DNase | +1.65 | +1.32 (0.43) | 0.004 | 80% | +1.22 (0.010) | +0.65 (0.25) | +0.76 (0.32) | +0.78 (0.23) |
::end

::table Supplementary Table S16|Layer map with larger GWAS for six traits.|landscape
For each trait, the first row is the GWAS of the main map and the second a recent, larger GWAS (European ancestry; 23andMe excluded; UK Biobank included for height, BMI, bipolar disorder and LDL cholesterol; for AD, the 2026 consensus meta-analysis in its release without proxy cases). The larger GWAS include the participants of the earlier ones, so this is a test of robustness, not an independent replication. *N*: median per-SNP sample size after processing (effective sample size for the larger schizophrenia, bipolar disorder and AD GWAS). *h*^2^~g~ (observed scale, with its *Z* score), mean *χ*^2^, the intercept (s.e.) and the number of regression SNPs are from the all-track layer model. Lineage τ\*: single-lineage model for the lineage used, with its rank among 12 lineages. CAGE, DNase and H3K27ac: τ\* in the three-layer joint model, averaged over 20 resamples; CAGE − DNase with its jackknife standard error; layer type by the classification rule of the main map. Effect given activity: τ\* of the own-layer effect annotation (DNase for height, CAGE otherwise) conditional on binary predicted activity (all tracks), with its retention relative to the effect-only model. With the six traits replaced, layer profiles still differed between traits (*Q* = 64.8, 20 degrees of freedom, *P* = 1.3 × 10^−6^; main map 74.3; 68.4 with the first five replaced), and the contrasts of four brain traits and of AD with height held (+1.66 ± 0.41, *P* = 4.6 × 10^−5^; +4.36 ± 1.42, *P* = 0.002; criterion two-sided *P* < 0.017; with the new AD GWAS and the original height GWAS, +5.19 ± 1.47). The criterion set for height (DNase first and CAGE − DNase < 0 with one-sided *P* < 0.05) was not met (one-sided *P* = 0.054). For AD, the lineage order (myeloid, lymphoid, brain tissue) and the effect τ\* under each activity construction were as in the original GWAS (Supplementary Table S18).
| Trait | GWAS | *N* | *h*^2^~g~ (*Z*) | Mean *χ*^2^ | Intercept (s.e.) | Regression SNPs | Lineage τ\* (rank) | CAGE | DNase | H3K27ac | CAGE − DNase (s.e.) | Layer type | Effect given activity (retention) |
| Height | Lango Allen 2010 [@langoallen2010] | 133,719 | 0.221 (17.7) | 1.330 | 0.786 (0.015) | 1,025,325 | Stromal +1.54 (3) | +0.49 | +1.92 | +0.04 | −1.43 (0.55) | DNase type | +1.32 (80%) |
|  | Yengo 2022 [@yengo2022] | 1,591,205 | 0.417 (23.0) | 14.665 | 2.193 (0.187) | 1,038,098 | Stromal +1.35 (3) | +0.88 | +1.48 | +0.12 | −0.61 (0.38) | Unresolved (DNase modal, 16/20) | +1.18 (79%) |
| BMI | Speliotes 2010 [@speliotes2010] | 123,822 | 0.146 (16.4) | 1.108 | 0.768 (0.011) | 1,026,547 | Brain tissue +0.53 (1) | +0.90 | +0.06 | −0.24 | +0.84 (0.50) | CAGE type | +0.35 (40%) |
|  | Yengo 2018 [@yengo2018] | 688,775 | 0.217 (32.4) | 3.930 | 1.067 (0.030) | 1,008,852 | Brain tissue +0.83 (1) | +0.97 | +0.05 | +0.11 | +0.92 (0.30) | CAGE type | +0.54 (55%) |
| Schizophrenia | PGC 2014 [@pgc2014] | 70,100 | 0.458 (17.9) | 1.697 | 1.061 (0.018) | 1,082,608 | Brain tissue +1.55 (1) | +1.53 | +0.07 | +0.19 | +1.46 (0.28) | CAGE type | +0.86 (56%) |
|  | Trubetskoy 2022 [@trubetskoy2022] | 58,749 | 0.821 (28.4) | 2.016 | 1.095 (0.014) | 1,165,447 | Brain tissue +1.42 (1) | +1.40 | +0.20 | +0.28 | +1.20 (0.24) | CAGE type | +0.83 (57%) |
| Bipolar disorder | PGC 2011 [@pgc2011] | 16,731 | 0.392 (5.7) | 1.180 | 1.031 (0.014) | 750,389 | Brain tissue +1.91 (2) | +4.26 | +0.47 | −1.17 | +3.79 (0.99) | CAGE type | +4.89 (102%) |
|  | O'Connell 2025 [@oconnell2025] | 163,367 | 0.172 (22.1) | 1.614 | 1.074 (0.011) | 1,146,171 | Brain tissue +1.76 (1) | +1.38 | +0.71 | +0.16 | +0.67 (0.33) | CAGE type | +0.79 (51%) |
| LDL cholesterol | Teslovich 2010 [@teslovich2010] | 95,386 | 0.102 (8.5) | 1.103 | 0.932 (0.012) | 1,017,468 | Liver +3.21 (1) | +3.90 | +0.83 | +1.13 | +3.07 (1.07) | CAGE type | +2.46 (63%) |
|  | Graham 2021 [@graham2021] | 1,230,202 | 0.093 (10.5) | 3.041 | 1.087 (0.044) | 1,165,597 | Liver +2.65 (1) | +3.54 | +0.28 | +1.09 | +3.26 (0.87) | CAGE type | +1.95 (55%) |
| Alzheimer's disease | Bellenguez 2022 [@bellenguez2022] | 487,511 | 0.026 (10.7) | 1.265 | 1.053 (0.012) | 1,162,232 | Myeloid +2.98 (1) | +3.22 | +0.14 | +0.01 | +3.08 (1.37) | CAGE type | +0.42 (14%) |
|  | 2026 meta-analysis [@ad2026meta] | 197,002 | 0.089 (13.1) | 1.350 | 1.057 (0.013) | 1,162,861 | Myeloid +2.33 (1) | +3.28 | −0.48 | +0.54 | +3.76 (1.37) | CAGE type | +0.66 (23%) |
::end

::table Supplementary Table S17|Layer contrasts with continuous annotations and non-overlapping bins.
CAGE − DNase difference in τ\* (jackknife standard error) in each trait's lineage, all tracks. Top 5%: the all-track three-layer model of the main map. Continuous: the percentiles of the three layers fitted jointly as continuous annotations. Bins: four non-overlapping bins of the percentile with the boundaries of the threshold analysis (top 0–0.1%, 0.1–0.5%, 0.5–1% and 1–5% of SNPs), all 12 bin annotations fitted jointly. Direction: stable when the contrast has the sign of the top-5% model with continuous annotations and in all four bins (criterion fixed before the analysis), with the number of these five criteria with the same sign. Heterogeneity between the 11 map traits: continuous *Q* = 43.6 (20 degrees of freedom, *P* = 0.002); bins *Q* = 58.1 (*P* = 1.4 × 10^−5^), 31.2 (*P* = 0.052), 25.8 (*P* = 0.17) and 58.0 (*P* = 1.4 × 10^−5^) from the top bin down. The validation contrast held in every criterion (one-sided *P* ≤ 0.02). Variance inflation factors of the baseline-adjusted LD scores were 1.18–3.35 for the continuous and 1.02–1.53 for the bin annotations. The lower rows use the larger GWAS of Supplementary Table S16, with their own top-5% model as the reference.
| Trait | Top 5% | Continuous | 0–0.1% | 0.1–0.5% | 0.5–1% | 1–5% | Direction |
| Educational attainment | +0.92 (0.33) | +0.92 (0.31) | −0.09 (0.43) | +0.34 (0.41) | −0.77 (0.38) | +1.11 (0.36) | Sensitive (3/5) |
| Schizophrenia | +1.39 (0.37) | +1.53 (0.33) | +0.09 (0.43) | +0.51 (0.50) | −0.35 (0.40) | +1.51 (0.38) | Sensitive (4/5) |
| Bipolar disorder | +2.50 (1.23) | +2.44 (0.95) | −0.50 (1.52) | +2.10 (1.51) | +0.38 (1.44) | +2.17 (1.27) | Sensitive (4/5) |
| Neuroticism | +0.72 (0.55) | +0.57 (0.45) | −1.17 (0.58) | −0.21 (0.56) | −0.01 (0.57) | +0.96 (0.54) | Sensitive (2/5) |
| BMI | +0.66 (0.54) | +0.87 (0.49) | +0.57 (0.53) | −0.66 (0.54) | +0.78 (0.68) | +0.60 (0.59) | Sensitive (4/5) |
| Alzheimer's disease | +2.70 (1.57) | +1.03 (1.08) | +4.15 (1.93) | +0.77 (1.38) | +0.77 (1.35) | +2.06 (1.29) | Stable (5/5) |
| Monocyte count | +3.78 (0.97) | +1.79 (0.64) | +2.99 (0.85) | +2.81 (1.14) | +2.22 (1.12) | +2.79 (0.74) | Stable (5/5) |
| IBD | +1.03 (0.90) | +1.32 (0.60) | +0.78 (1.48) | +0.92 (1.15) | −0.92 (1.09) | +0.87 (0.98) | Sensitive (4/5) |
| Rheumatoid arthritis | +2.18 (1.23) | +1.50 (0.79) | −1.51 (1.84) | +0.99 (1.76) | +0.21 (1.28) | +1.91 (1.17) | Sensitive (4/5) |
| LDL cholesterol | +2.97 (1.12) | +1.38 (0.82) | −3.30 (1.66) | +1.33 (1.40) | +2.14 (1.92) | +2.25 (0.98) | Sensitive (4/5) |
| Height | −1.71 (0.60) | −0.71 (0.52) | −0.72 (0.86) | −0.68 (0.73) | −1.37 (0.52) | −1.34 (0.59) | Stable (5/5) |
| Bone mineral density | −1.34 (0.45) | −1.68 (0.46) | −0.61 (0.60) | −1.15 (0.55) | −1.12 (0.42) | −0.84 (0.44) | Stable (5/5) |
| Neutrophil count | +1.04 (0.78) | +0.84 (0.42) | +2.56 (0.76) | +1.60 (1.02) | +0.13 (0.72) | +0.54 (0.71) | Stable (5/5) |
| Total cholesterol | +3.59 (1.00) | +1.13 (0.76) | −0.71 (0.91) | +2.22 (1.15) | +1.02 (1.19) | +2.95 (0.77) | Sensitive (4/5) |
| Height, larger GWAS | −1.47 (0.45) | −0.85 (0.36) | +0.10 (0.49) | −0.93 (0.53) | −0.85 (0.29) | −1.16 (0.44) | Sensitive (4/5) |
| BMI, larger GWAS | +0.78 (0.35) | +0.69 (0.25) | +0.04 (0.25) | −0.29 (0.26) | +0.15 (0.29) | +0.87 (0.36) | Sensitive (4/5) |
| Schizophrenia, larger GWAS | +0.93 (0.33) | +1.22 (0.26) | +0.09 (0.35) | +0.40 (0.38) | −0.33 (0.31) | +1.00 (0.34) | Sensitive (4/5) |
| Bipolar disorder, larger GWAS | +0.48 (0.44) | +0.78 (0.34) | −0.43 (0.49) | +0.07 (0.52) | −0.97 (0.50) | +0.92 (0.40) | Sensitive (3/5) |
| LDL cholesterol, larger GWAS | +3.44 (0.95) | +1.20 (0.71) | +0.32 (1.12) | +0.98 (0.82) | +1.31 (1.24) | +2.88 (0.79) | Stable (5/5) |
::end

::table Supplementary Table S18|Predicted allelic effects conditional on predicted activity: alternative constructions and effects within predicted active elements.|landscape
All tracks of each trait's lineage, own layer from the layer map (DNase for height and BMD, CAGE otherwise). Effect only: τ\* of the own-layer effect annotation in the joint model of the three effect layers. Binary: with the three top-5% activity annotations as covariates (Supplementary Table S15). Continuous: with the three activity percentiles as continuous covariates. Stratified: τ\* of the own-layer effect annotation built as the top 5% of effects within each of 100 activity strata, three layers fitted jointly. *q*: one-sided, Benjamini–Hochberg across the 11 map traits (for the validation traits and the larger GWAS, one-sided *P* without correction). Criterion fixed before the analysis: the lineage pattern of the binary model (retained for educational attainment, schizophrenia, bipolar disorder, LDL cholesterol and height; explained by activity for AD, monocyte count, IBD and rheumatoid arthritis) was reproduced if at least 4 of the 5 retained traits and at most 1 of the 4 explained traits had *q* < 0.05; it was reproduced in neither construction (continuous, 3 of 5 and 3 of 4; stratified, 0 of 5 and 0 of 4). Within top 10% (20%): Δ = τ\*(larger effect) − τ\*(smaller effect) between activity-matched halves of the SNPs in the top 10% (20%) of predicted activity of the same lineage and layer (post hoc; jackknife s.e.). Direction: number of the five annotation criteria of Supplementary Table S17 and the four activity constructions (binary, continuous, stratified, binary with measured annotations) in which the CAGE − DNase contrast had the sign of the top-5% all-track model (for the larger GWAS, eight criteria). The lower rows give, for each model, the heterogeneity of the (CAGE − DNase, CAGE − H3K27ac) contrasts across the 11 map traits and the contrasts with height (joint jackknife), and the same quantities for the larger GWAS of Supplementary Table S16. Variance inflation factors of the baseline-adjusted LD scores were at most 2.6 (continuous), 1.2 (stratified) and 4.7 (within active elements).
| Trait | Own layer | Effect only | Binary (*q*) | Continuous (*q*) | Stratified (*q*) | Within top 10% Δ (s.e.; *q*) | Within top 20% Δ (s.e.; *q*) | Direction |
| Educational attainment | CAGE | +1.22 | +0.55 (0.038) | +0.37 (0.094) | +0.06 (0.63) | +0.32 (0.32; 0.30) | +0.41 (0.32; 0.22) | 6/9 |
| Schizophrenia | CAGE | +1.53 | +0.86 (0.004) | +0.39 (0.094) | -0.19 (0.95) | +1.08 (0.41; 0.022) | +0.62 (0.35; 0.11) | 8/9 |
| Bipolar disorder | CAGE | +4.79 | +4.89 (7.4 × 10^−6^) | +4.08 (1.9 × 10^−5^) | +1.72 (0.25) | +2.73 (1.12; 0.027) | +3.03 (1.04; 0.010) | 7/9 |
| Neuroticism | CAGE | +0.91 | +0.34 (0.28) | +0.12 (0.36) | +0.15 (0.63) | +0.19 (0.45; 0.53) | -0.05 (0.43; 0.60) | 5/9 |
| BMI | CAGE | +0.88 | +0.35 (0.28) | +0.24 (0.27) | -0.59 (0.97) | +0.03 (0.47; 0.62) | +0.38 (0.45; 0.37) | 6/9 |
| Alzheimer's disease | CAGE | +2.91 | +0.42 (0.40) | +2.33 (0.012) | -0.21 (0.86) | -0.19 (1.66; 0.62) | -0.50 (1.17; 0.67) | 7/9 |
| Monocyte count | CAGE | +3.47 | +0.02 (0.49) | +2.20 (0.007) | -0.80 (0.97) | +0.97 (1.00; 0.30) | +0.62 (1.29; 0.44) | 8/9 |
| IBD | CAGE | +2.32 | +0.08 (0.49) | +1.30 (0.031) | +0.46 (0.63) | +1.43 (0.84; 0.12) | +1.33 (0.76; 0.11) | 4/9 |
| Rheumatoid arthritis | CAGE | +2.98 | +0.40 (0.40) | +1.64 (0.065) | +0.19 (0.63) | -0.57 (1.21; 0.68) | +0.17 (1.16; 0.54) | 5/9 |
| LDL cholesterol | CAGE | +3.91 | +2.46 (0.028) | +3.28 (0.004) | +0.96 (0.47) | -0.23 (1.47; 0.62) | +0.62 (1.30; 0.44) | 7/9 |
| Height | DNase | +1.65 | +1.32 (0.004) | +1.10 (0.014) | +0.49 (0.33) | +2.28 (0.49; 2 × 10^−5^) | +2.45 (0.49; 3.4 × 10^−6^) | 9/9 |
| *q* < 0.05 | | | 5 of 11 | 6 of 11 | 0 of 11 | 3 of 11 | 2 of 11 | |
| Bone mineral density | DNase | +1.51 | +1.23 | +1.02 | +0.63 | +2.37 | +2.25 | 9/9 |
| Neutrophil count | CAGE | +1.35 | -0.70 | +0.59 | -0.35 | -0.13 | -0.22 | 6/9 |
| Total cholesterol | CAGE | +3.70 | +2.72 | +3.18 | -0.04 | +2.49 | +1.72 | 8/9 |
| Heterogeneity *Q* (*P*) | | 70.1 (1.7 × 10^−7^) | 53.1 (7.9 × 10^−5^) | 60.5 (5.9 × 10^−6^) | 24.0 (0.24) | 57.8 (1.5 × 10^−5^) | 57.5 (1.7 × 10^−5^) | |
| Brain traits vs height (*P*) | | +2.81 ± 0.63 (8 × 10^−6^) | +2.63 ± 0.64 (4.5 × 10^−5^) | +1.58 ± 0.65 (0.015) | +0.83 ± 0.50 (0.094) | +3.01 ± 0.83 (2.8 × 10^−4^) | +2.82 ± 0.85 (9 × 10^−4^) | |
| AD vs height (*P*) | | +4.41 ± 1.65 (0.008) | +2.01 ± 1.63 (0.22) | +4.01 ± 1.61 (0.013) | +0.59 (0.60) | +0.68 ± 2.28 (0.77) | +1.11 ± 1.79 (0.53) | |
| Validation contrast (one-sided *P*) | | +3.35 ± 0.74 (3 × 10^−6^) | +1.48 ± 0.80 (0.032) | +2.22 (0.0017) | +0.03 (0.48) | +2.91 ± 1.05 (0.003) | +1.94 ± 1.00 (0.026) | |
| Height, larger GWAS | DNase | +1.48 | +1.18 (2.6 × 10^−4^) | +1.00 (0.002) | +0.56 (0.013) | +1.84 (0.44; 1.4 × 10^−5^) | +1.85 (0.38; 6.7 × 10^−7^) | 7/8 |
| BMI, larger GWAS | CAGE | +0.98 | +0.54 (0.002) | +0.46 (0.008) | -0.01 (0.53) | +0.25 (0.24; 0.14) | +0.35 (0.23; 0.061) | 7/8 |
| Schizophrenia, larger GWAS | CAGE | +1.46 | +0.83 (1.5 × 10^−4^) | +0.42 (0.033) | +0.01 (0.47) | +0.74 (0.33; 0.013) | +0.28 (0.29; 0.16) | 5/8 |
| Bipolar disorder, larger GWAS | CAGE | +1.56 | +0.79 (0.003) | +0.57 (0.019) | +0.04 (0.44) | +0.90 (0.43; 0.019) | +1.03 (0.41; 0.006) | 3/8 |
| LDL cholesterol, larger GWAS | CAGE | +3.53 | +1.95 (0.014) | +2.88 (0.005) | +0.04 (0.48) | +1.42 (1.14; 0.11) | +1.01 (1.18; 0.20) | 7/8 |
| Alzheimer's disease, larger GWAS | CAGE | +2.91 | +0.66 (0.23) | +2.16 (0.005) | -0.41 (0.79) | +0.35 (1.44; 0.40) | +0.00 (0.97; 0.50) | — |
::end

::table Supplementary Table S19|Predicted annotations conditional on fine-mapped cell-type eQTLs.|landscape
Five continuous MaxCPP annotations from FunGen-xQTL [@fungen2026] (microglia, neurons, other glia, bulk brain, monocytes; Methods) were added to the effect-only model and to the effect + binary-activity model of each trait's lineage, all tracks. Shown are the own-layer effect τ\* without and with the cell-type eQTL annotations, with its retention (primary test for AD and the five brain traits: one-sided *q* by Benjamini–Hochberg across the six; the other traits descriptive, one-sided *P*), the same for the effect + activity model, the myeloid CAGE activity τ\* for AD (primary test), and the τ\* (one-sided *P*, uncorrected) of each cell-type eQTL annotation in the effect + eQTL model, which also contains the GTEx and BLUEPRINT MaxCPP annotations of the baseline. Numbers of target SNPs with MaxCPP > 0.1: microglia 103,619; neurons 34,318; other glia 27,708; bulk brain 71,398; monocytes 3,444. With the cell-type eQTL annotations in the model, the layer contrasts still differed between the 11 map traits (*Q* = 70.3, 20 degrees of freedom, *P* = 1.6 × 10^−7^; 70.1 without), and the contrasts of AD and of the brain traits with height were +4.32 ± 1.60 (*P* = 0.007) and +2.74 ± 0.61 (*P* = 7 × 10^−6^); the validation contrast was +3.29 ± 0.74 (one-sided *P* = 5 × 10^−6^). Variance inflation factors of the baseline-adjusted LD scores were at most 2.1, and the partial correlations between the effect and eQTL annotations at most 0.15. The lower rows use the larger GWAS of Supplementary Table S16.
| Trait | Own layer | Effect only | + cell-type eQTLs (s.e.) | Retention (*q* or *P*) | Effect + activity | + cell-type eQTLs (s.e.) | Retention (*P*) | Microglia | Neurons | Other glia | Bulk brain | Monocytes |
| Educational attainment | CAGE | +1.22 | +1.14 (0.23) | 94% (*q* = 6.6 × 10^−7^) | +0.55 | +0.52 (0.26) | 94% (0.022) | +0.22 (0.045) | +0.09 (0.24) | +0.21 (0.052) | -0.06 (0.70) | +0.07 (0.32) |
| Schizophrenia | CAGE | +1.53 | +1.36 (0.25) | 89% (*q* = 1.1 × 10^−7^) | +0.86 | +0.76 (0.27) | 88% (0.002) | +0.18 (0.11) | +0.57 (3.1 × 10^−4^) | +0.30 (0.034) | -0.10 (0.76) | +0.17 (0.17) |
| Bipolar disorder | CAGE | +4.79 | +4.57 (0.90) | 95% (*q* = 5.9 × 10^−7^) | +4.89 | +4.61 (1.00) | 94% (1.9 × 10^−6^) | +0.15 (0.37) | +1.05 (0.011) | -0.90 (0.98) | +0.81 (0.056) | +0.53 (0.19) |
| Neuroticism | CAGE | +0.91 | +0.87 (0.35) | 97% (*q* = 0.007) | +0.34 | +0.33 (0.37) | 99% (0.18) | +0.19 (0.19) | +0.47 (0.004) | +0.14 (0.21) | -0.54 (0.98) | +0.05 (0.40) |
| BMI | CAGE | +0.88 | +0.82 (0.35) | 94% (*q* = 0.009) | +0.35 | +0.31 (0.39) | 90% (0.21) | -0.16 (0.78) | -0.02 (0.54) | +0.15 (0.25) | +0.26 (0.11) | +0.06 (0.40) |
| Alzheimer's disease | CAGE | +2.91 | +2.73 (0.88) | 94% (*q* = 0.002) | +0.42 | +0.37 (0.91) | 88% (0.34) | +0.44 (0.11) | -0.07 (0.58) | +1.12 (0.016) | +0.43 (0.15) | +0.54 (0.13) |
| Alzheimer's disease, myeloid CAGE activity | | | | | +1.91 | +1.87 (0.65) | 98% (0.002) | | | | | |
| Monocyte count | CAGE | +3.47 | +3.29 (0.62) | 95% (*P* = 5.5 × 10^−8^) | +0.02 | +0.00 (0.71) | 16% (0.50) | +0.36 (0.058) | -0.46 (0.98) | -0.61 (1.00) | +0.52 (0.065) | +1.37 (1.6 × 10^−4^) |
| IBD | CAGE | +2.32 | +2.33 (0.59) | 101% (*P* = 3.7 × 10^−5^) | +0.08 | +0.08 (0.74) | 106% (0.46) | -0.34 (0.85) | +0.18 (0.32) | -1.17 (1.00) | +0.21 (0.31) | +1.14 (0.013) |
| Rheumatoid arthritis | CAGE | +2.98 | +2.91 (0.80) | 98% (*P* = 1.4 × 10^−4^) | +0.40 | +0.38 (0.89) | 93% (0.34) | -0.01 (0.51) | -0.95 (0.98) | -0.19 (0.66) | +1.07 (0.012) | +0.45 (0.22) |
| LDL cholesterol | CAGE | +3.91 | +3.76 (0.84) | 96% (*P* = 3.8 × 10^−6^) | +2.46 | +2.36 (1.04) | 96% (0.012) | +0.01 (0.49) | -1.15 (0.99) | +0.91 (0.11) | +0.66 (0.17) | +1.22 (0.10) |
| Height | DNase | +1.65 | +1.62 (0.43) | 98% (*P* = 7.2 × 10^−5^) | +1.32 | +1.30 (0.42) | 99% (9.5 × 10^−4^) | +0.41 (0.023) | +0.33 (0.090) | +0.09 (0.34) | +0.18 (0.22) | -0.12 (0.69) |
| Bone mineral density | DNase | +1.51 | +1.51 (0.35) | 100% (*P* = 9.2 × 10^−6^) | +1.23 | +1.23 (0.35) | 100% (1.8 × 10^−4^) | -0.07 (0.63) | +0.21 (0.16) | -0.02 (0.55) | +0.25 (0.047) | -0.06 (0.64) |
| Neutrophil count | CAGE | +1.35 | +1.28 (0.45) | 95% (*P* = 0.002) | -0.70 | -0.73 (0.46) | 104% (0.94) | -0.16 (0.82) | -0.38 (0.90) | +0.12 (0.28) | +0.35 (0.14) | +1.17 (0.009) |
| Total cholesterol | CAGE | +3.70 | +3.51 (0.84) | 95% (*P* = 1.5 × 10^−5^) | +2.72 | +2.60 (1.09) | 96% (0.009) | +0.61 (0.15) | -0.16 (0.62) | +0.36 (0.26) | -0.05 (0.55) | +0.50 (0.25) |
| Height, larger GWAS | DNase | +1.48 | +1.48 (0.35) | 100% (*P* = 1.1 × 10^−5^) | | | | +0.14 (0.21) | +0.24 (0.061) | -0.14 (0.83) | +0.18 (0.14) | -0.17 (0.92) |
| BMI, larger GWAS | CAGE | +0.98 | +0.93 (0.18) | 95% (*P* = 2.1 × 10^−7^) | | | | +0.00 (0.49) | +0.09 (0.23) | +0.06 (0.32) | +0.20 (0.047) | +0.00 (0.49) |
| Schizophrenia, larger GWAS | CAGE | +1.46 | +1.34 (0.22) | 92% (*P* = 2.6 × 10^−10^) | | | | +0.14 (0.13) | +0.43 (2.6 × 10^−4^) | +0.16 (0.12) | -0.00 (0.52) | +0.14 (0.18) |
| Bipolar disorder, larger GWAS | CAGE | +1.56 | +1.33 (0.25) | 85% (*P* = 3.3 × 10^−8^) | | | | +0.31 (0.014) | +0.81 (1.4 × 10^−6^) | +0.18 (0.17) | +0.20 (0.17) | +0.14 (0.25) |
| LDL cholesterol, larger GWAS | CAGE | +3.53 | +3.39 (0.89) | 96% (*P* = 6.5 × 10^−5^) | | | | +0.36 (0.24) | -1.02 (0.99) | +0.44 (0.19) | +0.64 (0.13) | +0.49 (0.18) |
| Alzheimer's disease, larger GWAS | CAGE | +2.91 | +2.64 (0.87) | 91% (*P* = 0.001) | +1.63 | +1.54 (0.56) | 94% (0.003) | +0.95 (7.3 × 10^−4^) | -0.31 (0.86) | +0.81 (0.020) | +0.28 (0.19) | +0.61 (0.089) |
::end

# Supplementary Figures

::figure figS1.png|Supplementary Figure S1|Lineage-specific heritability enrichment.
τ\* of each lineage annotation in models that also contain an all-lineage annotation of SNPs with large predicted effects in any lineage (baselineLD + all-lineage + one lineage annotation). Bold values with an asterisk: false discovery rate < 0.05 across the 132 cells.
::end

::figure figS2.png|Supplementary Figure S2|Single-layer annotations by lineage.
τ\* of the 72 lineage–layer annotations (six layers within each of 12 lineages), each fitted alone with baselineLD v2.2. Layers are labelled above the columns (K27ac, K4me1 and K4me3: H3K27ac, H3K4me1 and H3K4me3). Dots: false discovery rate < 0.05 across the 792 cells (54 cells).
::end

::figure figS3.png|Supplementary Figure S3|Layer map at stricter annotation thresholds.
Top layer of each trait in the three-layer joint model at the top 0.1%, 0.5%, 1% and 5% of SNPs (all tracks). Dark cells with an asterisk: the top layer exceeds the second with one-sided *P* < 0.05. The heterogeneity *P* across the 11 map traits (joint jackknife) is shown under each threshold. Validation traits are below the dotted line.
::end

::figure figS4.png|Supplementary Figure S4|CAGE − DNase contrasts across annotation criteria.
**a**, Difference between the CAGE and DNase τ\* (95% confidence intervals) in the three-layer joint model with all tracks of each trait's lineage at nested thresholds (top 0.1%, 0.5%, 1% and 5% of SNPs; Supplementary Table S7). **b**, The same contrast with continuous percentile annotations of the three layers fitted jointly and in four non-overlapping bins (top 0–0.1%, 0.1–0.5%, 0.5–1% and 1–5%), with all 12 bin annotations fitted jointly (Supplementary Table S17). Standard errors are from the jackknife over aligned blocks. The 5% estimates are from the all-track model and differ slightly from the resample averages of Fig. 3. The validation traits are below the line.
::end

::figure figS5.png|Supplementary Figure S5|Size comparisons of predicted, measured and eQTL annotations.
Each trait is shown at its own layer from the layer map (in parentheses). **a**, τ\* of the predicted and the measured annotation of the same layer, fitted together; measured annotations were built from FANTOM5 (CAGE) or ENCODE (DNase, H3K27ac) data for the same ontology terms as the Atlas tracks. The predicted annotation had the larger τ\* in all 11 traits, and adding it reduced the measured τ\* by a median of 29% for CAGE (80% for AD), 5% for DNase and 2% for H3K27ac. **b**, τ\* of the predicted own-layer annotation and of the fine-mapped GTEx eQTL annotation (MaxCPP) in the same model. In **a** and **b**, asterisks mark differences with two-sided *P* < 0.05, uncorrected for multiple testing. **c**, Predicted minus measured τ\* for each layer across the 11 traits (points) with the median (bar). Measured annotations covered 0.5–1.2% (CAGE), 5–23% (DNase) and 12–17% (H3K27ac) of common SNPs, against 5% for the predictions; because τ\* scales with the standard deviation of an annotation, these differences partly reflect annotation size and are not a measure of relative performance. The difference was largest for CAGE and absent for H3K27ac; besides annotation size, differences between layers in resolution and between lineages in track number and prediction accuracy could contribute, and we did not distinguish these. Error bars are 95% confidence intervals.
::end
