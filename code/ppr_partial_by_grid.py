"""Round 2 (Reviewer 2, F1 and point 2): what the risk score adds to distance on each grid (Supplementary Table S9 d).

Partial Spearman correlation of the risk score (and of goat density) with the arrival month, with the Euclidean distance
from the entry routes held fixed: the Pearson correlation of the residuals of the ranks, each regressed on the ranks of
distance (average ranks for ties). Over the ROK cells reached, from the per-cell table written by
`PPR_model_v9_10.py --batch targeting` (PPR_v9.10_R2_targeting_cells.csv). 95% intervals from 300 bootstrap resamples
of 32 x 32 km blocks of cells (seed 20261123), as for Table S7 (b).

python ppr_partial_by_grid.py CELLS_CSV CELL_KM [CELLS_CSV CELL_KM ...]
e.g. python ppr_partial_by_grid.py v9_10_f1ci_g/PPR_v9.10_R2_targeting_cells.csv 8 \
                                   v9_10_f1ci4_g/PPR_v9.10_R2_targeting_cells.csv 4 \
                                   v9_10_f1ci2_g/PPR_v9.10_R2_targeting_cells.csv 2
"""
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def resid(y, x):
    z = np.column_stack([np.ones(y.size), x])
    return y - z @ np.linalg.lstsq(z, y, rcond=None)[0]


def partial(a, b, c):
    rc = rankdata(c)
    ea, eb = resid(rankdata(a), rc), resid(rankdata(b), rc)
    return float((ea * eb).sum() / np.sqrt((ea * ea).sum() * (eb * eb).sum()))


def run(path, cell_km, n_boot=300, block_km=32.0, seed=20261123):
    c = pd.read_csv(path, comment="#", encoding="utf-8-sig")
    c = c[c["arrival_month"].notna()].reset_index(drop=True)
    y, d = c["arrival_month"].to_numpy(float), c["distance_km"].to_numpy(float)
    sc = np.nan_to_num(c["risk_score"].to_numpy(float), nan=-1e9)
    gd = c["goat_density"].to_numpy(float)
    b = max(1, int(round(block_km / float(cell_km))))
    blk = (c["row"].to_numpy() // b).astype(np.int64) * 1_000_003 + (c["col"].to_numpy() // b)
    _u, inv = np.unique(blk, return_inverse=True)
    members = [np.flatnonzero(inv == j) for j in range(_u.size)]
    rng = np.random.default_rng(seed)
    bs_r, bs_g = [], []
    for _ in range(n_boot):
        k = np.concatenate([members[j] for j in rng.integers(0, len(members), len(members))])
        bs_r.append(partial(sc[k], y[k], d[k]))
        bs_g.append(partial(gd[k], y[k], d[k]))
    lo_r, hi_r = np.percentile(bs_r, [2.5, 97.5])
    lo_g, hi_g = np.percentile(bs_g, [2.5, 97.5])
    print(f"{cell_km:g}-km grid: {y.size} ROK cells reached, {len(members)} blocks; partial rho, risk score | distance "
          f"{partial(sc, y, d):+.3f} (95% CI {lo_r:+.3f} to {hi_r:+.3f}); goat density | distance "
          f"{partial(gd, y, d):+.3f} (95% CI {lo_g:+.3f} to {hi_g:+.3f})")


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a or len(a) % 2:
        raise SystemExit(__doc__)
    for i in range(0, len(a), 2):
        run(a[i], float(a[i + 1]))
