"""Round 2 (Reviewer 2, point 1): censoring-aware checks of the Latin hypercube sensitivity ranking
(Supplementary Table S10 a, Methods S1.12).

For the month of first arrival in the ROK, with the runs that never reached the ROK handled in several ways:
  1. Cox proportional-hazards model, not-reached runs censored at the month in which they stopped (as in Table S10);
  2. the same Cox model with every not-reached run censored at month 498, the end of follow-up (a run in which the
     epidemic died out or stalled cannot reach the ROK later);
  3. log-normal and 4. Weibull accelerated-failure-time (AFT) models with the censoring of 2, which do not assume
     proportional hazards;
  5. the proportional-hazards test of model 1 (scaled Schoenfeld residuals, rank time).
Parameters are rescaled to [0, 1] over their sampled ranges (gamma on a log scale). AFT z values are reported with the
sign reversed so that, as for the Cox model, a positive value means earlier arrival.

python ppr_sa_survival_checks.py PPR_v9.10_SA_lhs_runs.csv [out.csv]
Needs numpy, pandas and lifelines (pip install lifelines).
"""
import sys
import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, LogNormalAFTFitter, WeibullAFTFitter
from lifelines.statistics import proportional_hazard_test

warnings.filterwarnings("ignore")
PARAMS = ["gamma", "beta_level", "w_h", "diff_max_x", "beta_amp", "psi", "tv_amp", "a_star"]


def main(path, out=None):
    df = pd.read_csv(path, comment="#", encoding="utf-8-sig")
    df = df[df["ok"].astype(str) == "True"].copy()
    keys = [k for k in PARAMS if k in df.columns]
    arrived = df["first_rok"].notna()
    cap = float(df["end_month"].max() + 1)                      # 498: 480 months after the DPRK introduction
    Z = df[keys].astype(float).copy()
    if "gamma" in Z:
        Z["gamma"] = np.log(Z["gamma"])
    Z = (Z - Z.min()) / (Z.max() - Z.min())
    A = Z.copy()
    A["T"] = np.where(arrived, df["first_rok"], df["end_month"]) - df["intro_month"]
    A["E"] = arrived.astype(int)
    B = Z.copy()
    B["T"] = np.where(arrived, df["first_rok"], cap) - df["intro_month"]
    B["E"] = arrived.astype(int)
    c1 = CoxPHFitter().fit(A, duration_col="T", event_col="E")
    c2 = CoxPHFitter().fit(B, duration_col="T", event_col="E")
    ln = LogNormalAFTFitter().fit(B, duration_col="T", event_col="E")
    wb = WeibullAFTFitter().fit(B, duration_col="T", event_col="E")
    ph = proportional_hazard_test(c1, A, time_transform="rank").summary["p"]
    tab = pd.DataFrame({
        "cox_stop_log_hr": c1.summary["coef"].reindex(keys), "cox_stop_z": c1.summary["z"].reindex(keys),
        "cox_cap_z": c2.summary["z"].reindex(keys),
        "lognormal_aft_z_earlier": -ln.summary.loc["mu_"]["z"].reindex(keys),
        "weibull_aft_z_earlier": -wb.summary.loc["lambda_"]["z"].reindex(keys),
        "ph_test_p_rank": ph.reindex(keys)}, index=keys)
    tab.index.name = "parameter"
    pd.set_option("display.width", 200)
    print(f"runs {len(df)}; reached the ROK {int(arrived.sum())}; not reached {int((~arrived).sum())}; "
          f"end of follow-up month {cap:g}")
    print(f"concordance: Cox (censored at stop) {c1.concordance_index_:.3f}; Cox (censored at {cap:g}) "
          f"{c2.concordance_index_:.3f}; log-normal AFT {ln.concordance_index_:.3f}; Weibull AFT "
          f"{wb.concordance_index_:.3f}; AIC log-normal {ln.AIC_:.1f}, Weibull {wb.AIC_:.1f}")
    print(tab.round(3).to_string())
    for c in ("cox_stop_z", "cox_cap_z", "lognormal_aft_z_earlier", "weibull_aft_z_earlier"):
        print(f"ranking by |{c}|:", list(tab[c].abs().sort_values(ascending=False).index))
    if out:
        tab.to_csv(out, encoding="utf-8-sig")
        print("written", out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
