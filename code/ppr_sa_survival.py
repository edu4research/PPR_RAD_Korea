"""Round 2: censoring-aware re-analysis of the Latin hypercube sensitivity runs.

The PRCC of the manuscript sets runs that never reached the ROK to one common month (the latest end + 1).
This script checks whether the parameter ranking depends on that choice:
  1. PRCC on the runs that reached the ROK only (no censoring);
  2. a Cox proportional-hazards model of the time from introduction into the DPRK to first arrival in the ROK,
     with runs that never arrived treated as right-censored at the end of the run;
  3. a logistic model of whether the ROK was reached at all.
Parameters are scaled to [0, 1] over their sampled ranges (γ on a log scale).

python ppr_sa_survival.py PPR_v9.10_SA_lhs_runs.csv [out.csv]
Needs numpy, pandas, scipy and lifelines (pip install lifelines) or statsmodels.
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

PARAMS = ["gamma", "beta_level", "w_h", "diff_max_x", "beta_amp", "psi", "tv_amp", "a_star"]
LOG = {"gamma"}


def prcc(X, y):
    R = np.column_stack([stats.rankdata(X[:, j]) for j in range(X.shape[1])])
    ry = stats.rankdata(y)
    out = []
    for j in range(X.shape[1]):
        Z = np.column_stack([np.ones(len(ry)), np.delete(R, j, axis=1)])
        rx = R[:, j] - Z @ np.linalg.lstsq(Z, R[:, j], rcond=None)[0]
        res = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
        out.append(float(np.corrcoef(rx, res)[0, 1]))
    return np.array(out)


def main(path, out=None):
    df = pd.read_csv(path, comment="#", encoding="utf-8-sig")
    df = df[df["ok"].astype(str) == "True"].copy()
    keys = [k for k in PARAMS if k in df.columns]
    arrived = df["first_rok"].notna()
    end = df["end_month"].astype(float)
    df["T"] = np.where(arrived, df["first_rok"], end) - df["intro_month"]
    df["E"] = arrived.astype(int)
    Z = df[keys].astype(float).copy()
    for k in keys:
        if k in LOG:
            Z[k] = np.log(Z[k])
        Z[k] = (Z[k] - Z[k].min()) / (Z[k].max() - Z[k].min())
    tab = pd.DataFrame(index=keys)
    cens = float(end.max() + 1)
    tab["prcc_censored_at_common_month"] = prcc(df[keys].to_numpy(float),
                                                np.where(arrived, df["first_rok"], cens))
    sub = df[arrived]
    tab["prcc_arrived_only"] = prcc(sub[keys].to_numpy(float), sub["first_rok"].to_numpy(float))
    try:
        from lifelines import CoxPHFitter
        cph = CoxPHFitter().fit(pd.concat([Z, df[["T", "E"]]], axis=1), duration_col="T", event_col="E")
        tab["cox_log_hr_per_range"] = cph.summary["coef"].reindex(keys).values
        tab["cox_z"] = cph.summary["z"].reindex(keys).values
        conc = cph.concordance_index_
    except ImportError:
        import statsmodels.api as sm
        m = sm.PHReg(df["T"].to_numpy(float), Z.to_numpy(float), status=df["E"].to_numpy(int)).fit()
        tab["cox_log_hr_per_range"] = m.params
        tab["cox_z"] = m.tvalues
        conc = float("nan")
    try:
        import statsmodels.api as sm
        lg = sm.Logit(df["E"], sm.add_constant(Z)).fit(disp=0)
        tab["logit_reached_coef"] = lg.params[keys].values
        tab["logit_reached_z"] = lg.tvalues[keys].values
    except Exception:
        pass
    pd.set_option("display.width", 200)
    print(f"runs {len(df)}; reached the ROK {int(arrived.sum())}; censored {int((~arrived).sum())}; "
          f"Cox concordance {conc:.3f}")
    print(tab.round(3).to_string())
    for c in ("prcc_censored_at_common_month", "prcc_arrived_only", "cox_z"):
        print(f"ranking by |{c}|:", list(tab[c].abs().sort_values(ascending=False).index))
    if out:
        tab.round(4).to_csv(out, encoding="utf-8-sig")
    return tab


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
