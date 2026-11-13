# Reproduction check of stages e19–e22 and of the figures

Reference: `experiments/e19_…` to `e22_…` and `figures/out/` of this repository  
Tolerance: rtol 1e-06, atol 1e-09

One run from a copy made with `tools/make_run_copy.py` (empty `results/` folders, no `figures/out/`): every step of e19–e22 and the figure scripts, with `run_all.sh`, in the
data root that the full re-run of e01–e18 had written (inputs and the outputs of e01–e18 only). Every annotation, LD score and S-LDSC fit of these stages was therefore computed
anew. The new fits differ from the stored ones in the last digits (relative differences below 1e-9 in the tables computed from them and below 1e-6 in the mixed-model tables), so
tables are identical in content and not always byte for byte.

| Kind of file | Status | Files |
|---|---|---|
| result table | match | 19 |
| result table | same bytes | 13 |
| text output | match | 16 |
| text output | match (timing/path lines only) | 7 |
| log of a piece of a split step | match | 6 |
| log of a piece of a split step | match (timing/path lines only) | 80 |
| figure (PNG) | same bytes | 10 |
| values drawn in the figures | match | 5 |
| values drawn in the figures | same bytes | 6 |

| File | Status | Detail |
|---|---|---|
| e19_enrichment_all_groups/results/t01_check.txt | match (timing/path lines only) | 2 lines |
| e19_enrichment_all_groups/results/t01_enrichment.csv | same bytes |  |
| e19_enrichment_all_groups/results/t01_extract.txt | match | 63 lines |
| e19_enrichment_all_groups/results/t01_fit.txt | match (timing/path lines only) | 12 lines |
| e19_enrichment_all_groups/results/t08_fit.txt | match (timing/path lines only) | 9 lines |
| e20_gwas_panel/results/u00_inputs.csv | same bytes |  |
| e20_gwas_panel/results/u01_check.txt | match | 3 lines |
| e20_gwas_panel/results/u01_fitR_0of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitR_10of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitR_11of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitR_1of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitR_2of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitR_3of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitR_4of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitR_5of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitR_6of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitR_7of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitR_8of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitR_9of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_0of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_10of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_11of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_1of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_2of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_3of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_4of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_5of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_6of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_7of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_8of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fitS_9of12.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fit_0of8.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fit_1of8.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fit_2of8.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fit_3of8.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fit_4of8.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fit_5of8.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fit_6of8.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_fit_7of8.txt | match (timing/path lines only) | 1 lines |
| e20_gwas_panel/results/u01_prep.csv | same bytes |  |
| e20_gwas_panel/results/u01_prep.txt | match | 1 lines |
| e20_gwas_panel/results/u02_analysis_lineage.csv | match | 106 rows x 23 columns |
| e20_gwas_panel/results/u02_lineage.csv | match | 1272 rows x 13 columns |
| e20_gwas_panel/results/u02_lineage.txt | match | 74 lines |
| e20_gwas_panel/results/u03_pairs.csv | match | 5 rows x 24 columns |
| e20_gwas_panel/results/u03_pairs.txt | match | 7 lines |
| e20_gwas_panel/traits_exp3u.csv | same bytes |  |
| e21_panel_single_layer_map/results/v01_check.txt | match (timing/path lines only) | 3 lines |
| e21_panel_single_layer_map/results/v01_fit_0of10.txt | match (timing/path lines only) | 1 lines |
| e21_panel_single_layer_map/results/v01_fit_1of10.txt | match (timing/path lines only) | 1 lines |
| e21_panel_single_layer_map/results/v01_fit_2of10.txt | match (timing/path lines only) | 1 lines |
| e21_panel_single_layer_map/results/v01_fit_3of10.txt | match (timing/path lines only) | 1 lines |
| e21_panel_single_layer_map/results/v01_fit_4of10.txt | match (timing/path lines only) | 1 lines |
| e21_panel_single_layer_map/results/v01_fit_5of10.txt | match (timing/path lines only) | 1 lines |
| e21_panel_single_layer_map/results/v01_fit_6of10.txt | match (timing/path lines only) | 1 lines |
| e21_panel_single_layer_map/results/v01_fit_7of10.txt | match (timing/path lines only) | 1 lines |
| e21_panel_single_layer_map/results/v01_fit_8of10.txt | match (timing/path lines only) | 1 lines |
| e21_panel_single_layer_map/results/v01_fit_9of10.txt | match (timing/path lines only) | 1 lines |
| e21_panel_single_layer_map/results/v01_prepld.txt | match (timing/path lines only) | 1 lines |
| e21_panel_single_layer_map/results/v02_gwas106.csv | same bytes |  |
| e21_panel_single_layer_map/results/v02_main11.csv | same bytes |  |
| e21_panel_single_layer_map/results/v02_main11.txt | match | 15 lines |
| e21_panel_single_layer_map/results/v02_matrix106.csv | same bytes |  |
| e21_panel_single_layer_map/results/v02_pairs106.csv | same bytes |  |
| e21_panel_single_layer_map/results/v02_types.txt | match | 103 lines |
| e22_benchmark14/results/w01_check.txt | match | 4 lines |
| e22_benchmark14/results/w01_fit_actbin_0of7.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_actbin_1of7.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_actbin_2of7.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_actbin_3of7.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_actbin_4of7.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_actbin_5of7.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_actbin_6of7.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_bench_0of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_bench_1of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_bench_2of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_bench_3of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_bench_4of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_bench_5of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_bench_6of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_bench_7of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_bench_8of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_bench_9of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panel_0of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panel_1of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panel_2of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panel_3of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panel_4of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panel_5of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panel_6of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panel_7of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panel_8of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panel_9of10.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panelpd_0of6.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panelpd_1of6.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panelpd_2of6.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panelpd_3of6.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panelpd_4of6.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w01_fit_panelpd_5of6.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w02_check.txt | match (timing/path lines only) | 2 lines |
| e22_benchmark14/results/w02_fit_0of5.txt | match (timing/path lines only) | 2 lines |
| e22_benchmark14/results/w02_fit_1of5.txt | match (timing/path lines only) | 2 lines |
| e22_benchmark14/results/w02_fit_2of5.txt | match (timing/path lines only) | 2 lines |
| e22_benchmark14/results/w02_fit_3of5.txt | match (timing/path lines only) | 2 lines |
| e22_benchmark14/results/w02_fit_4of5.txt | match (timing/path lines only) | 1 lines |
| e22_benchmark14/results/w03_bins.csv | same bytes |  |
| e22_benchmark14/results/w03_check.txt | match (timing/path lines only) | 57 lines |
| e22_benchmark14/results/w03_main_0of6.txt | match | 4 lines |
| e22_benchmark14/results/w03_main_1of6.txt | match | 4 lines |
| e22_benchmark14/results/w03_main_2of6.txt | match | 4 lines |
| e22_benchmark14/results/w03_main_3of6.txt | match | 4 lines |
| e22_benchmark14/results/w03_main_4of6.txt | match | 3 lines |
| e22_benchmark14/results/w03_main_5of6.txt | match | 3 lines |
| e22_benchmark14/results/w03_overlap.csv | same bytes |  |
| e22_benchmark14/results/w04_activity.csv | match | 14 rows x 19 columns |
| e22_benchmark14/results/w04_allgroups.csv | match | 168 rows x 20 columns |
| e22_benchmark14/results/w04_bins.csv | match | 56 rows x 8 columns |
| e22_benchmark14/results/w04_blocks.csv | match | 14 rows x 8 columns |
| e22_benchmark14/results/w04_contrasts.csv | same bytes |  |
| e22_benchmark14/results/w04_controls.csv | match | 42 rows x 15 columns |
| e22_benchmark14/results/w04_layers.csv | match | 14 rows x 25 columns |
| e22_benchmark14/results/w04_map72.csv | same bytes |  |
| e22_benchmark14/results/w04_measured.csv | match | 14 rows x 23 columns |
| e22_benchmark14/results/w04_summary.txt | match | 171 lines |
| e22_benchmark14/results/w04_tissue.csv | match | 168 rows x 13 columns |
| e22_benchmark14/results/w04_types.csv | same bytes |  |
| e22_benchmark14/results/w05_annotations.csv | match | 12 rows x 13 columns |
| e22_benchmark14/results/w05_pairs.csv | match | 288 rows x 27 columns |
| e22_benchmark14/results/w05_summary.txt | match | 74 lines |
| e22_benchmark14/results/w06_actbin.csv | match | 14 rows x 20 columns |
| e22_benchmark14/results/w06_summary.txt | match | 16 lines |
| e22_benchmark14/results/w07_gwas_resid.csv | match | 294 rows x 8 columns |
| e22_benchmark14/results/w07_summary.txt | match | 52 lines |
| e22_benchmark14/results/w08_pairs.csv | match | 117 rows x 25 columns |
| e22_benchmark14/results/w08_summary.txt | match | 52 lines |
| e22_benchmark14/results/w09_model.csv | match | 117 rows x 11 columns |
| e22_benchmark14/results/w09_summary.txt | match | 62 lines |
| e22_benchmark14/results/w09b_summary.txt | match | 14 lines |
| e22_benchmark14/results/w10_contrasts.csv | match | 39 rows x 8 columns |
| e22_benchmark14/results/w10_means.csv | match | 37 rows x 8 columns |
| e22_benchmark14/results/w10_summary.txt | match | 29 lines |
| figures/out/eqtl14.csv | same bytes |  |
| figures/out/fig1.png | same bytes |  |
| figures/out/fig2.png | same bytes |  |
| figures/out/fig3.png | same bytes |  |
| figures/out/fig4.png | same bytes |  |
| figures/out/fig5.png | same bytes |  |
| figures/out/figS1.png | same bytes |  |
| figures/out/figS2.png | same bytes |  |
| figures/out/figS3.png | same bytes |  |
| figures/out/figS4.png | same bytes |  |
| figures/out/figS5.png | same bytes |  |
| figures/out/numbers_fig1.json | same bytes |  |
| figures/out/numbers_fig2.json | same bytes |  |
| figures/out/numbers_fig3.json | match | 352 numbers |
| figures/out/numbers_fig4.json | match | 837 numbers |
| figures/out/numbers_fig5.json | match | 416 numbers |
| figures/out/numbers_figS1.json | same bytes |  |
| figures/out/numbers_figS2.json | same bytes |  |
| figures/out/numbers_figS3.json | same bytes |  |
| figures/out/numbers_figS4.json | match | 222 numbers |
| figures/out/numbers_figS5.json | match | 146 numbers |
