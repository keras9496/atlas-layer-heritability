# Code index

Paths are relative to `experiments/`. Scripts are listed in the order in which they are run (`run_all.sh` gives the
exact order and the dependencies between steps). Every folder also has a `results/` folder into which its scripts write; the result files are not tracked in this repository.
"Not included" marks the AlphaGenome Atlas query scripts, which are left out of the repository; their outputs are inputs
(see `data_layout.md`). Stages e01–e18 analyse 11 traits and three further traits with the GWAS first used; stages e19–e22 add
the later GWAS of six traits to a benchmark of 14 traits and analyse a public panel of 106 GWAS. The scripts and file names
call the tissue groups lineages.

Shared helpers: `e12_diagnostics/scripts/m_common.py` (τ*, joint jackknife, contrasts, Cochran's *Q*; imported by
e13, e14, e17) and `e15_activity_constructions/scripts/p_common.py` (the same for the activity models; imported by
e16, e18). `e02_lineage_map/scripts/lineage_rules.py` holds the 12-lineage assignment rules verbatim.
`e22_benchmark14/scripts/wlib_panel.py` loads the enriched pairs of the panel (imported by w08–w10; w09b and w10 import the
model of `w09_weighted.py`).

## Analyses

The folders follow the order in which the analyses were run, not the order of the manuscript. For each analysis, in the
order of the Results, this table names the scripts that write the reported values and the S-LDSC fits and annotations that those scripts
read. The result file behind every figure and table is listed under *Figures and tables* at the end of this file.

| Analysis | Figures and tables | Scripts that write the reported values | Fits and annotations they read |
|---|---|---|---|
| Tissue-group map of the 14 traits and of the panel | Fig. 1; Supplementary Figs. S1, S2; Supplementary Tables S1–S4; Supplementary Data 1 | e22 `w04_bench.py` (part 1 of `w04_summary.txt`: tissue-group map of the 14 traits, without and with the all-group annotation); e19 `t01_enrichment.py` (heritability share and enrichment); e20 `u02_lineage.py` (tissue-group map of the panel); table of the GWAS: e12 `m01_trait_info.py`, e13 `n01_munge.py`, `n03_tests.py`, e17 `r01_munge.py`, `r03_tests.py`; tracks per tissue group and layer: e02 `inputs/track_lineage.csv`, `inputs/z01_lineages.txt` | Tissue-group models: e22 `w01_fit.py fit bench` (single-group models of eight traits; models with the all-group annotation of all 14), e13 `n02_sldsc.py` and e17 `r02_sldsc.py` (single-group models of the six later GWAS), e20 `u01_fit.py fit` (panel). Heritability share: e19 `t01_enrichment.py fit` and the LDSC outputs of e02, e04 and e05. Annotations: e02 `z03_annot_ldscores.py`, e03 `y03_annot_ldscores.py` |
| Map of tissue groups and layers (72 annotations, each fitted alone); types of signal | Fig. 2; Supplementary Fig. S3; Supplementary Tables S5–S7; Supplementary Note 1; Supplementary Data 3 | e22 `w04_bench.py` (part 2: the 72 annotations and the types of signal of the 14 traits); e21 `v02_types.py` (the 106 GWAS) | Each of the 72 annotations fitted alone: e22 `w02_fit72.py` (nine GWAS), e03 `y04_sldsc.py B` (five GWAS), e21 `v01_fit72.py` (panel). Annotations: e03 `y03_annot_ldscores.py` |
| Joint model of three layers (CAGE, DNase, H3K27ac) within tissue groups; CAGE–DNase balance of the 14 traits | Fig. 3; Supplementary Tables S8–S13; Supplementary Note 2 | e22 `w04_bench.py` (parts 3–5 and 9: joint model of the three layers, heterogeneity and contrasts between traits, CAGE–DNase balance, every enriched tissue group, cell-type-matched and proximal/distal layers, influential blocks); overlap of the layer annotations: e12 `m03_collinearity.py`; coordinates of the influential blocks: e12 `m04a_block_influence.py` | Three units per layer, 20 resamples: e09 `j02_sldsc.py` (eight traits), e13 `n02_sldsc.py`, e17 `r02_sldsc.py`. All tracks, every tissue group: e09 `j02_sldsc.py`, e19 `t08_matrix.py fit`, e22 `w01_fit.py fit bench`. Cell-type-matched and proximal/distal: e09 `j02_sldsc.py`, e10 `k02_sldsc.py`, e22 `w01_fit.py fit bench`. Annotations: e03 `y03_annot_ldscores.py` (all tracks), e05 `v01_resample_annot.py`, e08 `i01_matched_annot.py`, e10 `k01_split_annot.py`; blocks: e09 `j01_blocks.py` |
| CAGE–DNase balance of the enriched GWAS–tissue-group pairs of the panel; mixed model; annotation constructions | Fig. 4; Supplementary Fig. S4; Supplementary Tables S14–S16; Supplementary Note 2; Supplementary Data 2 | e22 `w05_panel.py`, `w07_central.py`, `w08_pd.py`, `w09_weighted.py`, `w09b_total.py`, `w10_constructions.py` (balance of the enriched GWAS–tissue-group pairs, mixed model, constructions, trait categories); e20 `u03_pairs.py` (two GWAS of the same trait) | e20 `u01_fit.py fit` (all tracks, 106 GWAS × 12 tissue groups), `fitR`, `fitS`; e22 `w01_fit.py fit panel`, `fit panelpd`. Annotations: e02 `z03_annot_ldscores.py`, e03 `y03_annot_ldscores.py`, e05 `v01_resample_annot.py`, e08 `i01_matched_annot.py`, e10 `k01_split_annot.py` |
| Measured annotations, fine-mapped eQTLs, predicted activity and bins of effect size | Fig. 5; Supplementary Fig. S5; Supplementary Tables S17–S22 | e22 `w04_bench.py` (parts 6–8: measured annotations, predicted activity, bins of effect size), `w03_actbin.py`, `w06_actbin.py`; e18 `s03_tests.py` (cell-type eQTL annotations) and `figures/fig5.py` (their *P* values and false discovery rate across the 14 traits, `eqtl14.csv`) | Measured annotations: e22 `w01_fit.py fit bench` with the annotations of e04 `w02_annot_ldscores.py`; the τ\* of each annotation fitted alone comes from the single-annotation fits (e22 `w02_fit72.py`, e03 `y04_sldsc.py B`). eQTL annotations: e18 `s01_annot_ldscores.py`, `s02_sldsc.py`. Predicted activity: e11 `l03_sldsc.py`, e13 `n02_sldsc.py`, e16 `q02_sldsc.py`, e17 `r02_sldsc.py`, e22 `w01_fit.py fit actbin` (annotations: e11 `l02_annot_ldscores.py`, e16 `q01_annot_ldscores.py`, e22 `w03_actbin.py`). Bins of effect size: e14 `o02_sldsc.py`, e22 `w01_fit.py fit bench` (annotations: e14 `o01_annot_ldscores.py`) |

Two tables belong to the Methods rather than to the analyses above: Supplementary Table S23 (earlier and later GWAS
of six traits) is laid out from e05 `v03_tests.py`, e12 `m01_trait_info.py`, e13
`n03_tests.py` and e17 `r03_tests.py`, and Supplementary Table S24 (robustness analyses of the earlier
AD GWAS) from e01 `x04_sldsc.py`, `x06_prepare_3b.py`, `x07_sldsc_3b.py`, e02 `z04_sldsc.py`, e03 `y04_sldsc.py`,
e04 `w05_own_layer.py`, e05 `v03_tests.py`, e11 `l05_fig_export.py` and e12 `m06_supp_export.py`.

The scripts that are not named above prepare the inputs of those that are (summary statistics, annotations, LD scores,
S-LDSC fits), or are the checks, dry runs and tests of the stage in which they were written. The test and summary scripts
of e01–e18 label their comparisons with the numbering in use when each stage was run (H1, H2 and so on, with a verdict
for each); the manuscript does not use this numbering. It reports the 14 traits through `w04_bench.py`, which reads the
fits of those stages and applies one set of definitions to all traits.

## Experiments and scripts

| Folder (stage) | Script | What it does | Main outputs (`results/`) |
|---|---|---|---|
| `e01_targets_sldsc` (3: target SNPs and initial AD analysis; 3b: AD robustness) | `x00_prepare_sumstats.py` | Bellenguez AD summary statistics → LDSC format, APOE removed, HapMap3 SNPs, `munge_sumstats` | — |
| | `x01_targets.py` | MAF ≥ 5% SNPs of the 1000 Genomes EUR hg38 reference; REF/ALT from the high-coverage release; Atlas query list | `x01_targets.txt` |
| | `x02_atlas_query.py` | Atlas query of myeloid tracks (**not included**) | — |
| | `x03_annot_ldscores.py` | Cell-type annotations from Atlas scores → LD scores | `x03_annot.txt` |
| | `x04_sldsc.py`, `x05_summary.py` | S-LDSC conditional on baselineLD v2.2; summary | `x04_sldsc.csv` |
| | `x06_prepare_3b.py`, `x07_sldsc_3b.py` | 3b: summary statistics for controls and loci removal (locus lists in `inputs/`); robustness models (thresholds, AVI conditioning, *APOE*, loci removal, leave-one-block-out) | `x07_3b*.csv` |
| `e02_lineage_map` (3c) | `z01_lineages.py` | Track metadata and lineage assignment (**not included**; rules in `lineage_rules.py`) | `inputs/z01_lineages.txt` (kept as an input) |
| | `z02_atlas_query.py` | Atlas query of all tracks (**not included**) | — |
| | `z03_annot_ldscores.py` | Lineage annotations → LD scores | `z03_annot.txt` |
| | `z04_sldsc.py`, `z05_summary.py` | 11 traits × 12 lineages (models M0–M2) | `z04_sldsc.csv`, `z05_map_tau_star.csv` |
| `e03_single_layer_map` (3d; 3d-A) | `y03_annot_ldscores.py` | Lineage × layer annotations (72) → LD scores (used by most later stages) | — |
| | `y04_sldsc.py A`, `B`, `y05_summary.py` | `A`: one layer over all tissue groups, the six layers of a tissue group fitted jointly, and track-matched joint models; `B`: each of the 72 annotations fitted alone (the single-annotation map) | `y04_B.csv` and others |
| | `y06_eqtl_compare.py` | 3d-A: predicted CAGE and RNA versus fine-mapped GTEx eQTL annotations (post hoc, descriptive) | `y06_*` |
| `e04_measured_annotations` (3e) | `w01_measured_sources.py` | FANTOM5 and ENCODE sources matched to Atlas ontology terms (downloads the files if absent) | `w01_sources.txt` |
| | `w02_annot_ldscores.py` | Measured annotations → LD scores | `w02_annot.txt` |
| | `w03_sldsc.py`, `w04_summary.py` | Predicted CAGE conditional on measured layers (primary test) | `w03_sldsc.csv`, `w04_*.csv` |
| | `w05_own_layer.py` | Own-layer comparisons (post hoc) | `w05_own_layer.csv` |
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
| | `l05_fig_export.py` | Values collected for tables (no new fits) | `l05_values.csv` |
| `e12_diagnostics` (3m) | `m_common.py` | Shared τ*, joint-jackknife and contrast functions (as in 3j) | — |
| | `m01_trait_info.py` | Heritability, mean χ², intercept and regression SNPs of each GWAS | `m01_trait_info.csv` |
| | `m02_threshold_contrasts.py` | Contrasts per nested threshold | `m02_threshold_contrasts.csv` |
| | `m03_collinearity.py` | Overlap, baseline-adjusted LD-score correlation, VIF | `m03_collinearity.csv` |
| | `m04a_block_influence.py`, `m04b_exclusion_refit.py` | Influential blocks; refits without them | `m04a_top_blocks.csv`, `m04b_summary.csv` |
| | `m05_resample_convergence.py` | Convergence of resample averages | `m05_convergence.csv` |
| | `m06_supp_export.py` | Values collected for tables (no new fits) | `m06_decomp.csv` |
| `e13_larger_gwas` (3n) | `n00_download.sh` | Downloads the public GWAS files | — |
| | `n01_munge.py`, `n02_sldsc.py`, `n03_tests.py` | Processing; layer, lineage and activity models; tests | `n03_*.csv`, `n03_tests.txt` |
| `e14_continuous_bins` (3o) | `o01_annot_ldscores.py`, `o02_sldsc.py`, `o03_tests.py` | Continuous percentile annotations and four non-overlapping bins; tests | `o03_*.csv`, `o03_tests.txt` |
| `e15_activity_constructions` (3p) | `p_common.py` | Shared functions for the activity models | — |
| | `p01_annot_ldscores.py`, `p02_sldsc.py`, `p03_tests.py` | Continuous activity covariates and activity-stratified effect annotations; tests | `p03_main_summary.csv`, `p03_main_tests.txt` |
| | `p04_fig_export.py`, `p06_direction_summary.py` | Values collected for figures and tables; direction summary across criteria | `p04_values.csv`, `p04_contrasts.csv`, `p06_direction_summary.csv` |
| `e16_within_active` (3q) | `q01_annot_ldscores.py`, `q02_sldsc.py`, `q03_tests.py` | Activity-matched halves within the top 10% (20%) of predicted activity (post hoc) | `q03_main_summary.csv`, `q03_values.csv` |
| `e17_larger_ad_gwas` (3r) | `r00_download.sh`, `r01_munge.py`, `r02_sldsc.py`, `r03_tests.py` | 2026 AD consensus meta-analysis without proxy cases: lineage, layer, activity and within-active models | `r03_summary.csv` |
| `e18_cell_type_eqtl` (3s) | `s00_download.sh`, `s01_annot_ldscores.py`, `s02_sldsc.py`, `s03_tests.py` | FunGen-xQTL fine-mapped cell-type eQTL annotations (MaxCPP); effect models conditional on them | `s03_summary.csv`, `s03_values.csv` |
| `e19_enrichment_all_groups` (3t) | `t01_enrichment.py check`, `fit`, `extract` | Share of heritability and enrichment of the tissue-group and layer annotations (LDSC with `--overlap-annot`; the six later GWAS are fitted here, the fits of the other GWAS are read from e02, e04 and e05). The script has a further stage, `scale`, which is not run: its comparison with the baseline annotations is not reported in the manuscript | `t01_enrichment.csv` |
| | `t08_matrix.py fit` | All-track joint models of the three layers for the 14 GWAS of e02–e05 in every tissue group (aligned blocks; the models already fitted in e09 are read). The fits are read by `w04_bench.py`, which writes the table of every enriched tissue group. The `test` stage of the script is not run: it tabulates the same fits for the earlier GWAS only, and the manuscript does not report it | `t08_fit.txt` |
| `e20_gwas_panel` (3u) | `u00_inputs.py` | Table of the 106 GWAS of the public panel and checks of the input files. Columns of the table: `group`, the trait category (S skeletal, BI blood and immune, LL lipid, BR brain and behaviour, U other); `prediction`, a label derived from the category and not read by any script in this repository; `status` and `discovery_trait`, whether the trait is one of the 14 traits; `pair`, independent GWAS of the same trait | `traits_exp3u.csv` (in the stage folder), `u00_inputs.csv` |
| | `u01_fit.py prep`, `check`, `fit`, `fitR`, `fitS` | Cleaned summary statistics; S-LDSC with aligned blocks: tissue-group model and all-track joint model of the three layers for 106 GWAS × 12 tissue groups (`fit`); joint models with three units per layer, 20 resamples, for two subsets of GWAS (`fitR`: GWAS whose tissue group chosen by `u02` is one of the five with resampled annotations; `fitS`: GWAS enriched in the stromal group). The resampled models of the remaining enriched pairs are fitted in e22 | `u01_*.txt`, `u01_prep.csv` |
| | `u02_lineage.py` | Tissue-group map of the panel (τ\*, one-sided *P*, FDR across the 1,272 combinations); for each GWAS the tissue group with the largest τ\* and the tissue group that defines the job list of `fitR` (the stromal group for the skeletal category, the top group otherwise) | `u02_lineage.csv`, `u02_analysis_lineage.csv` |
| | `u03_pairs.py` | Two independent GWAS of the same trait in the panel: CAGE–DNase balance of each GWAS in the tissue group with its largest τ\*, from the all-track joint model of the three layers, and the difference between the two GWAS with the block-jackknife standard error of the difference. Reported with the balance of the enriched pairs of the panel and in Supplementary Table S16, for the pairs in which the tissue group with the largest τ\* is the same for both GWAS and both are enriched in it. The script compares each pair in that tissue group only, not in the other tissue groups in which both GWAS are enriched | `u03_pairs.csv`, `u03_pairs.txt` |
| `e21_panel_single_layer_map` (3v) | `v01_fit72.py prepld`, `check`, `fit` | The 72 annotations (12 tissue groups × six layers) fitted one at a time for the 106 GWAS (LDSC default blocks, as in e03); `check` reproduces e03 for two traits | `v01_*.txt` |
| | `v02_types.py all` | Numbers of significant annotations, of tissue groups and of layers per tissue group; types of signal (focal, layer-spanning, tissue-spanning, tissue- and layer-spanning); relation to power; other FDR thresholds; the same indices for the 11 traits of e03 (`v02_main11.*`) | `v02_gwas106.csv`, `v02_matrix106.csv`, `v02_pairs106.csv`, `v02_types.txt`, `v02_main11.csv`, `v02_main11.txt` |
| `e22_benchmark14` (3w) | `w01_fit.py check`, `fit bench`, `fit panel`, `fit panelpd`, `fit actbin` | S-LDSC with aligned blocks. `bench`: models of the 14 traits not fitted in earlier stages (tissue-group models with and without the all-group annotation, all-track three-layer models in every tissue group, cell-type-matched, proximal/distal, measured and measured + predicted models, continuous and bin models of the later AD GWAS). `panel`: enriched GWAS–tissue-group pairs of the panel in five tissue groups, three units per layer × 20 resamples and cell-type-matched layers. `panelpd`: the same pairs split at 1 kb from TSSs. `actbin`: effect annotations with four bins of predicted activity | `w01_*.txt` |
| | `w02_fit72.py check`, `fit` | The 72 annotations one at a time for the nine GWAS without such fits in e03 | `w02_*.txt` |
| | `w03_actbin.py check`, `main` | Annotations of four bins of predicted activity (five tissue groups × three layers × four bins) → LD scores; overlap of the effect and activity annotations; properties of SNPs by size of the predicted effect | `w03_overlap.csv`, `w03_bins.csv` |
| | `w04_bench.py` | Tables of the 14 traits (no new fits): tissue-group map, 72-annotation map and types of signal, joint model of the three layers (20 resamples), heterogeneity and contrasts between traits, layers in every enriched tissue group, cell-type-matched and proximal/distal layers, measured annotations, predicted activity, bins of effect size, influential blocks | `w04_*.csv`, `w04_summary.txt` |
| | `w05_panel.py` | CAGE–DNase balance of the enriched GWAS–tissue-group pairs of the panel under three constructions; inverse-variance-weighted mean and median per tissue-group annotation; comparisons within GWAS and within annotation | `w05_pairs.csv`, `w05_annotations.csv`, `w05_summary.txt` |
| | `w06_actbin.py` | Effect annotations fitted with four bins of predicted activity | `w06_actbin.csv`, `w06_summary.txt` |
| | `w07_central.py` | Heterogeneity between tissue-group annotations under three constructions; differences within GWAS; two-factor models; variation between GWAS (post hoc) | `w07_summary.txt`, `w07_gwas_resid.csv` |
| | `wlib_panel.py` | Shared loader of the enriched pairs and their balance under five constructions (imported by w08–w10) | — |
| | `w08_pd.py` | Balance of the enriched pairs for SNPs within and beyond 1 kb of TSSs (post hoc) | `w08_pairs.csv`, `w08_summary.txt` |
| | `w09_weighted.py`, `w09b_total.py` | Mixed model with the known jackknife covariance of the estimates: annotation means, variation between GWAS and between pairs (REML; post hoc) | `w09_model.csv`, `w09_summary.txt`, `w09b_summary.txt` |
| | `w10_constructions.py` | Model means of five annotations under five constructions and their contrasts; model of the 12 all-track annotations | `w10_means.csv`, `w10_contrasts.csv`, `w10_summary.txt` |

## Figures and tables

The figures are drawn by the scripts in `figures/` (`python figures/make_figures.py`, or one script at a time). They read only
the result files listed here and write, next to each figure, the values drawn in it (`figures/out/numbers_<figure>.json`).
The tables and the Supplementary Data files of the manuscript are laid out from the result files in the right-hand column
(selection of rows and columns, rounding, names of traits); the text of the manuscript and the scripts that lay out its
tables are not part of this repository. Figure and table numbers are those of the current manuscript; docstrings of the
scripts of e01–e18 that name a figure or table use the numbering of the earlier version.

| Item | Script | Source stage and result file |
|---|---|---|
| Fig. 1 | `figures/fig1.py` | e22 `w04_tissue.csv`; track counts from e02 `inputs/z01_lineages.txt` |
| Fig. 2 | `figures/fig2.py` | e22 `w04_map72.csv`, `w04_types.csv`; e21 `v02_gwas106.csv` |
| Fig. 3 | `figures/fig3.py` | e22 `w04_layers.csv`, `w04_controls.csv`, `w04_summary.txt` |
| Fig. 4 | `figures/fig4.py` | e22 `w05_pairs.csv`, `w08_pairs.csv`, `w10_means.csv`, `w10_contrasts.csv` |
| Fig. 5 | `figures/fig5.py` | e22 `w04_map72.csv`, `w04_measured.csv`, `w04_activity.csv`, `w06_actbin.csv`; e18 `s03_values.csv` |
| Supplementary Fig. S1 | `figures/figS1.py` | e22 `w04_tissue.csv` (model with the all-group annotation) |
| Supplementary Fig. S2 | `figures/figS2.py` | e20 `u02_lineage.csv`, `u02_analysis_lineage.csv`; display names in `figures/panel_names.json` |
| Supplementary Fig. S3 | `figures/figS3.py` | e22 `w04_map72.csv`, `w04_types.csv` |
| Supplementary Fig. S4 | `figures/figS4.py` | e22 `w05_pairs.csv` |
| Supplementary Fig. S5 | `figures/figS5.py` | e22 `w04_bins.csv`, `w04_layers.csv` |
| Supplementary Table S1 | | e02 `inputs/track_lineage.csv`, `inputs/z01_lineages.txt` (rules in `lineage_rules.py`) |
| Supplementary Table S2 | | e12 `m01_trait_info.csv`; e13 `n01_munge.txt`, `n03_trait_info.csv`; e17 `r01_munge.txt`, `r03_tests.txt`; e22 `w04_tissue.csv`. *N* is the median of the N column of each processed summary-statistics file (logged by the munging scripts for the later GWAS) |
| Supplementary Table S3 | | e22 `w04_tissue.csv` |
| Supplementary Table S4 | | e19 `t01_enrichment.csv` |
| Supplementary Table S5 | | e22 `w04_types.csv`, `w04_map72.csv` |
| Supplementary Tables S6, S7, Supplementary Note 1 | | e21 `v02_gwas106.csv`, `v02_pairs106.csv`, `v02_types.txt` |
| Supplementary Table S8 | | e12 `m03_collinearity.csv` |
| Supplementary Tables S9, S10 | | e22 `w04_layers.csv`, `w04_contrasts.csv`, `w04_summary.txt` |
| Supplementary Table S11 | | e22 `w04_controls.csv`, `w04_summary.txt` |
| Supplementary Table S12 | | e22 `w04_blocks.csv` (block indices); coordinates of the blocks from e09 `j01_blocks.py`, as listed in e12 `m04a_top_blocks.csv` |
| Supplementary Table S13 | | e22 `w04_allgroups.csv` |
| Supplementary Table S14 | | e22 `w10_means.csv`, `w10_contrasts.csv`, `w05_annotations.csv`; the `w05`, `w07`, `w08` and `w09` summaries |
| Supplementary Table S15 | | e22 `w10_means.csv`, `w05_annotations.csv`, `w05_summary.txt`, `w09_summary.txt` |
| Supplementary Table S16, Supplementary Note 2 | | e22 `w07_summary.txt`, `w09_summary.txt`, `w09b_summary.txt`, `w05_summary.txt`; two GWAS of the same trait: e20 `u03_pairs.csv`, `u03_pairs.txt` |
| Supplementary Tables S17, S18 | | e22 `w04_measured.csv`, `w04_summary.txt` |
| Supplementary Table S19 | | e18 `s03_values.csv`, `s03_main_tests.txt`, `s01_check.txt` (matching of the fine-mapped variants); `figures/out/eqtl14.csv` (written by `fig5.py`) |
| Supplementary Table S20 | | e22 `w03_overlap.csv`, `w03_bins.csv` |
| Supplementary Table S21 | | e22 `w04_activity.csv`, `w06_actbin.csv` |
| Supplementary Table S22 | | e22 `w04_bins.csv`, `w04_layers.csv` |
| Supplementary Table S23 | | e12 `m01_trait_info.csv`; e05 `v03_main_layers.csv`; e13 `n01_munge.txt`, `n03_trait_info.csv`, `n03_layers.csv`, `n03_lineage.csv`, `n03_activity.csv`; e17 `r01_munge.txt`, `r03_summary.csv`, `r03_tests.txt` (*N* as in Supplementary Table S2) |
| Supplementary Table S24 | | e01 `x04_sldsc.csv`, `x06_prepare_3b.txt`, `x07_3b*.csv`; e02 `z04_sldsc.csv`; e03 `y04_B.csv`; e04 `w05_own_layer.csv`; e05 `v03_main_layers.csv`; e11 `l05_values.csv`; e12 `m06_decomp.csv` |
| Supplementary Data 1 | | e20 `u02_lineage.csv` |
| Supplementary Data 2 | | e22 `w05_pairs.csv`, `w08_pairs.csv` |
| Supplementary Data 3 | | e21 `v02_matrix106.csv` |
