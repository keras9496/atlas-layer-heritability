# Code index

Paths are relative to `experiments/`. Scripts are listed in the order in which they are run (`run_all.sh` gives the
exact order and the dependencies between steps). Every folder also holds `results/` with the result files the manuscript used. "Not included" marks the AlphaGenome Atlas query scripts,
which are left out of the repository; their outputs are inputs (see `data_layout.md`).

Shared helpers: `e12_diagnostics/scripts/m_common.py` (τ*, joint jackknife, contrasts, Cochran's *Q*; imported by
e13, e14, e17) and `e15_activity_constructions/scripts/p_common.py` (the same for the activity models; imported by
e16, e18). `e02_lineage_map/scripts/lineage_rules.py` holds the 12-lineage assignment rules verbatim.

## Experiments and scripts

| Folder (stage) | Script | What it does | Main outputs (`results/`) |
|---|---|---|---|
| `e01_targets_sldsc` (3: target SNPs and initial AD analysis; 3b: AD robustness) | `x00_prepare_sumstats.py` | Bellenguez AD summary statistics → LDSC format, APOE removed, HapMap3 SNPs, `munge_sumstats` | — |
| | `x01_targets.py` | MAF ≥ 5% SNPs of the 1000 Genomes EUR hg38 reference; REF/ALT from the high-coverage release; Atlas query list | `x01_targets.txt` |
| | `x02_atlas_query.py` | Atlas query of myeloid tracks (**not included**) | — |
| | `x03_annot_ldscores.py` | Cell-type annotations from Atlas scores → LD scores | `x03_annot.txt` |
| | `x04_sldsc.py`, `x05_summary.py` | S-LDSC conditional on baselineLD v2.2; summary | `x04_sldsc.csv` |
| | `x06_prepare_3b.py`, `x07_sldsc_3b.py` | 3b: summary statistics for controls and loci removal (locus lists in `inputs/`); robustness models (thresholds, AVI conditioning, *APOE*, loci removal, leave-one-block-out) | `x07_3b*.csv` |
| `e02_lineage_map` (3c) | `z01_lineages.py` | Track metadata and lineage assignment (**not included**; rules in `lineage_rules.py`) | `z01_lineages.txt` |
| | `z02_atlas_query.py` | Atlas query of all tracks (**not included**) | — |
| | `z03_annot_ldscores.py` | Lineage annotations → LD scores | `z03_annot.txt` |
| | `z04_sldsc.py`, `z05_summary.py` | 11 traits × 12 lineages (models M0–M2) | `z04_sldsc.csv`, `z05_map_tau_star.csv` |
| `e03_single_layer_map` (3d; 3d-A) | `y03_annot_ldscores.py` | Lineage × layer annotations (72) → LD scores (used by most later stages) | — |
| | `y04_sldsc.py A`, `B`, `y05_summary.py` | Single-layer (A) and six-layer joint (B) models | `y04_B.csv` and others |
| | `y06_eqtl_compare.py` | 3d-A: predicted CAGE and RNA versus fine-mapped GTEx eQTL annotations (post hoc, descriptive) | `y06_*` |
| `e04_measured_annotations` (3e) | `w01_measured_sources.py` | FANTOM5 and ENCODE sources matched to Atlas ontology terms (downloads the files if absent) | `w01_sources.txt` |
| | `w02_annot_ldscores.py` | Measured annotations → LD scores | `w02_annot.txt` |
| | `w03_sldsc.py`, `w04_summary.py` | Predicted CAGE conditional on measured layers (primary test) | `w03_sldsc.csv`, `w04_*.csv` |
| | `w05_own_layer.py` | Own-layer comparisons (post hoc; Fig. 5a) | `w05_own_layer.csv` |
| `e05_layer_map` (3f) | `v01_resample_annot.py` | Track-balanced resampled annotations (20 resamples) → LD scores | — |
| | `v02_sldsc.py dry|main|val`, `v03_tests.py` | Three-layer joint models; heterogeneity, pairwise contrasts, layer types; validation traits | `v03_main_layers.csv`, `v03_val_layers.csv` |
| | `v04_validation_sumstats.py` | Summary statistics of the three validation traits | — |
| | `v05_posthoc_loo.py` | Leave-one-trait-out heterogeneity (post hoc) | — |
| `e06_thresholds` (3g) | `u01_threshold_annot.py`, `u02_sldsc.py`, `u03_tests.py` | Top 1%, 0.5% and 0.1% layer annotations; layer map per threshold | `u03_main_layers.csv` |
| `e07_shared_lineages` (3h) | `h01_sldsc.py`, `h02_tests.py` | All-track models in every enriched lineage; skeletal versus lipid and blood traits | `h02_main_*.csv` |
| `e08_cell_type_matched` (3i) | `i01_matched_annot.py`, `i02_sldsc.py`, `i03_tests.py` | One unit per layer from matched cell types or regions | `i03_main_layers.csv` |
| `e09_joint_jackknife` (3j) | `j01_blocks.py` | 200 genomic blocks aligned across traits | — |
| | `j02_sldsc.py`, `j03_tests.py` | Refits with aligned blocks; tests between traits with their covariance (main *Q*, contrasts with height, validation contrast) | `j03_main_summary.csv` |
| `e10_proximal_distal` (3k) | `k01_split_annot.py`, `k02_sldsc.py`, `k03_tests.py` | Annotations split at 1 kb from annotated TSSs | `k03_main_layers.csv` |
| `e11_activity` (3l) | `l01_atlas_active.py` | Atlas activity-scorer query (**not included**) | — |
| | `l02_annot_ldscores.py` | Predicted-activity annotations → LD scores | `l02_raw_scale.csv` |
| | `l03_sldsc.py`, `l04_tests.py` | Effect + activity models | `l04_main_tests.txt`, `l04_main_summary.csv` |
| | `l05_fig_export.py` | Values for Supplementary Tables S12 and S15 (no new fits) | `l05_values.csv` |
| `e12_diagnostics` (3m) | `m_common.py` | Shared τ*, joint-jackknife and contrast functions (as in 3j) | — |
| | `m01_trait_info.py` | Trait information for Supplementary Table S1 | `m01_trait_info.csv` |
| | `m02_threshold_contrasts.py` | Contrasts per nested threshold | `m02_threshold_contrasts.csv` |
| | `m03_collinearity.py` | Overlap, baseline-adjusted LD-score correlation, VIF | `m03_collinearity.csv` |
| | `m04a_block_influence.py`, `m04b_exclusion_refit.py` | Influential blocks; refits without them (Table 1, Supplementary Table S14) | `m04a_top_blocks.csv`, `m04b_summary.csv` |
| | `m05_resample_convergence.py` | Convergence of resample averages | `m05_convergence.csv` |
| | `m06_supp_export.py` | Values for Supplementary Table S14 (no new fits) | `m06_decomp.csv` |
| `e13_larger_gwas` (3n) | `n00_download.sh` | Downloads the public GWAS files | — |
| | `n01_munge.py`, `n02_sldsc.py`, `n03_tests.py` | Processing; layer, lineage and activity models; tests | `n03_*.csv`, `n03_tests.txt` |
| `e14_continuous_bins` (3o) | `o01_annot_ldscores.py`, `o02_sldsc.py`, `o03_tests.py` | Continuous percentile annotations and four non-overlapping bins; tests | `o03_*.csv`, `o03_tests.txt` |
| `e15_activity_constructions` (3p) | `p_common.py` | Shared functions for the activity models | — |
| | `p01_annot_ldscores.py`, `p02_sldsc.py`, `p03_tests.py` | Continuous activity covariates and activity-stratified effect annotations; tests | `p03_main_summary.csv`, `p03_main_tests.txt` |
| | `p04_fig_export.py`, `p06_direction_summary.py` | Values for Fig. 5b and Supplementary Table S18; direction summary across criteria | `p04_values.csv`, `p04_contrasts.csv`, `p06_direction_summary.csv` |
| `e16_within_active` (3q) | `q01_annot_ldscores.py`, `q02_sldsc.py`, `q03_tests.py` | Activity-matched halves within the top 10% (20%) of predicted activity (post hoc; Fig. 5c) | `q03_main_summary.csv`, `q03_values.csv` |
| `e17_larger_ad_gwas` (3r) | `r00_download.sh`, `r01_munge.py`, `r02_sldsc.py`, `r03_tests.py` | 2026 AD consensus meta-analysis without proxy cases: lineage, layer, activity and within-active models | `r03_summary.csv` |
| `e18_cell_type_eqtl` (3s) | `s00_download.sh`, `s01_annot_ldscores.py`, `s02_sldsc.py`, `s03_tests.py` | FunGen-xQTL fine-mapped cell-type eQTL annotations (MaxCPP); effect models conditional on them | `s03_summary.csv`, `s03_values.csv` |

## Manuscript figures and tables

All figures are drawn by `paper/make_figures.py`, which reads only the result files listed here and writes
`paper/numbers.json` with the numbers quoted in the text. The text of the preprint is not part of this repository.

| Item | Source experiment and result file |
|---|---|
| Fig. 1 | Design schematic (`make_figures.py`) |
| Fig. 2, Supplementary Fig. S1 | e02 `z04_sldsc.csv` (models M1, M2) |
| Fig. 3 | e05 `v03_main_layers.csv`; *Q* from e09; larger GWAS from e13 `n03_layers.csv` and e17 `r03_summary.csv` |
| Fig. 4 | e05 `v03_val_layers.csv`; validation contrast from e09 `j03_main_summary.csv` |
| Fig. 5a | e04 `w05_own_layer.csv` |
| Fig. 5b | e15 `p04_values.csv`, `p03_main_summary.csv` |
| Fig. 5c | e16 `q03_main_summary.csv` |
| Table 1 | e09 `j03_main_summary.csv`; e08, e07, e10, e12 `m04b_summary.csv`, e13, e15, e18 |
| Supplementary Fig. S2 | e03 `y04_B.csv` |
| Supplementary Fig. S3 | e06 `u03_main_layers.csv` |
| Supplementary Fig. S4 | e12 `m02_threshold_contrasts.csv` (a), e14 `o03_contrasts_long.csv` (b) |
| Supplementary Fig. S5 | e04 `w05_own_layer.csv`, `w04_Sm.csv` |
| Supplementary Table S1 | e12 `m01_trait_info.csv` |
| Supplementary Table S2 | e02 `z01_lineages.txt` (rules in `lineage_rules.py`) |
| Supplementary Tables S3, S4 | e04 `w05_own_layer.csv`, `w04_*` |
| Supplementary Table S5 | e05 `v05_posthoc_loo.py`, recomputed with the joint jackknife in e09 |
| Supplementary Table S7 | e06, with the joint jackknife of e09 |
| Supplementary Table S8 | e08 `i03_main_layers.csv` |
| Supplementary Table S9 | e07 `h02_main_*.csv` |
| Supplementary Table S10 | e09 `j03_main_summary.csv` |
| Supplementary Table S11 | e10 `k03_main_layers.csv` |
| Supplementary Table S12 | e01 (`x04`, `x07`), e12 `m06_decomp.csv`, e11 `l05_values.csv` |
| Supplementary Table S13 | e12 `m03_collinearity.csv` |
| Supplementary Table S14 | e12 `m04b_summary.csv`, `m04a_top_blocks.csv`, `m06_decomp.csv` |
| Supplementary Table S15 | e11 `l04_main_summary.csv`, `l05_values.csv` |
| Supplementary Table S16 | e13 `n03_layers.csv`, `n03_trait_info.csv`, `n03_lineage.csv`, `n03_activity.csv`, `n03_contrasts.csv`; e17 `r03_summary.csv` |
| Supplementary Table S17 | e14 `o03_contrasts_long.csv`, `o03_stability.csv` |
| Supplementary Table S18 | e15 `p03_main_summary.csv`, `p04_values.csv`, `p04_contrasts.csv`; e16 `q03_main_summary.csv` |
| Supplementary Table S19 | e18 `s03_summary.csv`, `s03_values.csv` |
