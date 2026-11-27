# PPR_RAD_Korea v9.10 — code for the wildlife-borne PPR incursion model

Code, small public inputs and a runnable synthetic test for:

> Ahmad S, Park J, Jeong J, Jang H, Yoo DS. Wildlife-Borne Incursion of Peste des Petits Ruminants into the Republic
> of Korea, Conditional on Water-Deer Competence: A Reaction–Advection–Diffusion Scenario Analysis. *Microorganisms*
> (under review, manuscript microorganisms-4668700).

The model simulates the spread of peste des petits ruminants (PPR) across China, through the Democratic People's
Republic of Korea (DPRK) and across the fenced Military Demarcation Line into the Republic of Korea (ROK) with a
reaction–advection–diffusion equation on an 8-km grid (4 and 2 km as sensitivity runs), and ranks ROK areas with a
spatial multi-criteria risk score. Every equation, parameter value and option is in `code/PPR_model_v9_10.py`;
`python code/PPR_model_v9_10.py --help` lists the options.

## Contents

| Path | What it is |
|---|---|
| `code/PPR_model_v9_10.py` | The model (single script; md5 204bd491b2d7ec787a0b502d1c57f51d). Point locations are placed in the grid cell that contains them (default); `--point-cell round` reproduces runs made before 26 Nov 2026. |
| `code/run_model.py` | Runs the model with your input folder: `python code/run_model.py --data-root <inputs> [options]`. |
| `code/reproduce_all.sh` | Every run of the paper in a fixed order with fixed random seeds, then the external test and its audit. |
| `code/ppr_external_test.py` | External test of the hold-out run against the outbreaks of 2015–2025 (Section 2.5.6; Table S8). |
| `code/ppr_auc_audit.py` | Recomputes every count and the AUC of the external test from the raster and Table S13, plus the AUC among reached cells only. |
| `code/ppr_sa_survival.py`, `code/ppr_sa_survival_checks.py` | Censoring-aware analyses of the Latin hypercube runs (Table S10). |
| `code/ppr_partial_by_grid.py` | What the risk score adds to distance on each grid (Table S9 d). |
| `code/rok_goat_fill.py` | Builds the ROK part of the goat layer from the national register (the register itself is not shared; Methods S1.10). |
| `data/` | Public inputs of the external test: EMPRES-i outbreak records of China (cleaned), dated seed table, China boundary (Natural Earth 1:50 m), MaxEnt training points before 2015, and the 38 test records of Table S13. |
| `tests/` | Synthetic China–Korea world, its generator and an end-to-end test (`python run_synthetic_test.py`). |
| `environment.yml`, `requirements.txt` | Software environment (Python 3.9.20, NumPy 1.26.4, SciPy 1.13.1, GDAL 3.8.4, Matplotlib 3.9.2). |
| `INPUT_DATA.md` | The input layers needed for the full runs, their file names and sources. |

## Quick check (no restricted data, about 5–10 minutes)

    conda env create -f environment.yml
    conda activate ppr_rad
    cd tests
    python run_synthetic_test.py

It runs the full 8-km baseline (continental leg, peninsular leg, fence and corridor gates, risk score, internal
evaluation, maps and tables) on the synthetic world and ends with `TEST PASSED` when the timeline matches the stored
values. Outputs are in `tests/run/`.

## Full reproduction

1. Put the input layers in one folder as described in `INPUT_DATA.md`.
2. From an empty working folder: `bash /path/to/code/reproduce_all.sh /path/to/inputs` (several hours; the 2-km grid
   is the longest step). Tables and figures of each run are copied to `copies/`; logs are in `logs/`.
3. The run folders map to the paper as follows: `_g_pcfloor` (8-km run of record, Tables 1–4, Figures 7–9, Morris and
   Latin hypercube sample), `_r2_g_pcfloor` (comparator hosts, Table 5; gate-zone radii, Table S12; surveillance
   rules; sensitivity analysis, Figure 10, Table S10), `_f1ci*_g_pcfloor` (Table S9), `_g_pcfloor_grid4km` and
   `_grid2km` (Table S11), `_r2_4km_g_pcfloor` (Table S5 b), `_g_pcfloor_pubfigs` (Figures 1, 3–6),
   `_holdout2015_pcfloor` and its `external_test/` (Table S8).

## Data that are not included

The national register of goat holdings of the ROK contains personal information and is not shared; the goat-density
layer derived from it is available from the corresponding author on request. Without it the full runs cannot use
goats in the ROK (the synthetic test includes a synthetic goat layer).

## Licence

Code: MIT (see `LICENSE`). Files in `data/`: the EMPRES-i records are from FAO EMPRES-i (cite the source); the China
boundary is from Natural Earth (public domain); the other files are released under CC BY 4.0.

## Contact

Dae Sung Yoo (corresponding author), College of Veterinary Medicine, Chonnam National University, Gwangju, Republic of Korea.
