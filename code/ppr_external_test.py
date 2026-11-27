#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PPR v9.10 — external test of the risk maps against LATER, independent outbreaks
(item 7 of the revision plan; Reviewer 2, comment M4).

Independence first
  The run of record is seeded with every point of the outbreak file: all EMPRES-i PPR
  outbreaks in China 2007–2025 (242 locations) plus 24 locations in Mongolia. Its maps
  cannot be tested against those outbreaks. The test therefore uses a hold-out run that
  is seeded only with the outbreaks reported before 1 January 2015,

      python PPR_model_v9_10.py --preset holdout2015

  and scores that run's maps against the outbreaks of 2015 on (EMPRES-i: 38 records at
  33 locations), leaving out those within 10 km of a hold-out seed point.

In a notebook:

  from ppr_external_test import run_external_test
  run_external_test(points="EMPRESi_PPR_China_2007_2026.csv",
                    results_dir="results/ppr_v9_10_review_holdout2015")

  The cut-off date, the seed points to keep away from and the rasters are read from the
  run folder. The raw EMPRES-i export works as it is: the lines above its header are
  skipped, and each outbreak is dated by its observation date, else its report date.
  The background is the part of the grid where the outbreaks could have been recorded:
  China (china_boundary_ne50m.geojson next to this script; Natural Earth 1:50m). The
  whole grid is reported as a sensitivity row.

From a terminal:

  python ppr_external_test.py --points EMPRESi_PPR_China_2007_2026.csv \
      --results-dir results/ppr_v9_10_review_holdout2015

Reported per raster (external_test_summary.csv)
  AUC (95% CI)    outbreak pixels vs all valid pixels (Mann–Whitney; 0.5 = random);
                  bootstrap over outbreaks; one-sided p-value
  Boyce (95% CI)  continuous Boyce index (Hirzel et al. 2006, Ecol. Model. 199:142):
                  +1 = the predicted/expected ratio rises steadily with the score
  benchmark       the same AUC / Boyce for 'distance to the nearest seed point' on the
                  same pixels, and the paired AUC difference (map − distance) with its CI
  dates           Spearman rho between the simulated arrival month and the outbreak date
                  (arrival rasters only)
  classes         % of outbreaks and % of area in each risk class (P/E ratio)

Beyond distance (external_test_beyond_distance.csv, external_test_mixture.csv; round 2, Reviewer 2 F1)
  auc_within_distance_strata  each outbreak pixel compared only with background pixels in its own band of
                  distance to the nearest seed (deciles of the background); 0.5 = nothing beyond distance
  lr_chi2, lr_p   likelihood-ratio test of the map added to distance (presence/background logistic model,
                  both as percentile ranks, distance as a cubic); coef_map = log-odds across the map's whole
                  range (negative = the map points the wrong way once distance is known)
  best_lambda     weight of the map in the mixture (1-λ)·distance + λ·map with the highest AUC; mixture_gain
                  = its AUC minus that of distance (in-sample upper bound; CI with λ re-chosen per replicate)
  rho_date_*      arrival rasters: Spearman rho of outbreak date with arrival month (with its 95% CI) and with
                  distance, their difference, and the partial rho of date with arrival month given distance
Each outbreak pixel counts once (the earliest report) unless all_points=True. The same
statistics are repeated for the exclusion distances in sensitivity_km.
"""
import argparse
import csv
import datetime as _dt
import json
import math
from pathlib import Path

import numpy as np

try:
    from osgeo import gdal, ogr, osr
    gdal.UseExceptions()
    try:
        osr.UseExceptions()
        ogr.UseExceptions()
    except AttributeError:                                          # pragma: no cover
        pass
except ImportError as _e:                                           # pragma: no cover
    raise ImportError("ppr_external_test needs the GDAL Python bindings — run it in the same "
                      "environment as the model.") from _e

TAG = "PPR_v9.10"
EARTH_KM = 6371.0088


# ── statistics ────────────────────────────────────────────────────────────────────
def rank_avg(x):
    x = np.asarray(x, dtype=np.float64)
    _u, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    csum = np.cumsum(cnt)
    return ((csum - cnt + csum + 1) / 2.0)[inv]


def spearman(a, b):
    """Spearman rho with average ranks for ties; (rho, n)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if a.size < 3:
        return float("nan"), int(a.size)
    ra, rb = rank_avg(a), rank_avg(b)
    ra -= ra.mean()
    rb -= rb.mean()
    den = math.sqrt(float((ra * ra).sum()) * float((rb * rb).sum()))
    return (float((ra * rb).sum() / den) if den > 0 else float("nan")), int(a.size)


def _norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def spearman_p(rho, n):
    """Two-sided p for Spearman rho (t approximation)."""
    if not np.isfinite(rho) or n < 4:
        return float("nan")
    if abs(rho) >= 1.0:
        return 0.0
    t = rho * math.sqrt((n - 2) / (1.0 - rho * rho))
    try:
        from scipy.stats import t as _t
        return float(2.0 * _t.sf(abs(t), n - 2))
    except Exception:                                               # pragma: no cover
        return float(2.0 * _norm_sf(abs(t)))


def auc_parts(pres, back_sorted):
    """Per outbreak: share of background pixels with a lower score (ties count ½).
    Their mean is the AUC (the Mann–Whitney U / n1·n2)."""
    lo = np.searchsorted(back_sorted, pres, "left")
    hi = np.searchsorted(back_sorted, pres, "right")
    return (lo + 0.5 * (hi - lo)) / float(back_sorted.size)


def mann_whitney_p(pres, back_sorted):
    """One-sided p that outbreak pixels score higher than the background (normal
    approximation with tie correction and continuity correction)."""
    n1, n2 = int(pres.size), int(back_sorted.size)
    if n1 == 0 or n2 == 0:
        return float("nan")
    u = float(auc_parts(pres, back_sorted).sum()) * n2
    _v, t = np.unique(np.concatenate([pres, back_sorted]), return_counts=True)
    n = n1 + n2
    tie = float((t.astype(np.float64) ** 3 - t).sum()) / (n * (n - 1.0))
    var = n1 * n2 / 12.0 * ((n + 1.0) - tie)
    if var <= 0:
        return float("nan")
    return _norm_sf((u - n1 * n2 / 2.0 - 0.5) / math.sqrt(var))


def boyce_windows(back_sorted, window=0.10, step=0.01, min_back=10):
    """Moving windows of the continuous Boyce index over the background score range."""
    lo, hi = float(back_sorted[0]), float(back_sorted[-1])
    if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
        return None
    w = window * (hi - lo)
    k = int(math.floor((1.0 - window) / step + 1e-9)) + 1
    starts = lo + np.arange(k) * step * (hi - lo)
    ends = starts + w
    e_n = np.searchsorted(back_sorted, ends, "right") - np.searchsorted(back_sorted, starts, "left")
    return dict(starts=starts, ends=ends, mid=starts + w / 2.0, e_n=e_n, ok=e_n >= min_back,
                E=e_n / float(back_sorted.size))


def boyce_index(win, pres_sorted):
    """(Boyce index, presences per window, P, F)."""
    p_n = np.searchsorted(pres_sorted, win["ends"], "right") - np.searchsorted(pres_sorted, win["starts"], "left")
    P = p_n / float(max(pres_sorted.size, 1))
    F = np.where(win["E"] > 0, P / np.where(win["E"] > 0, win["E"], 1.0), np.nan)
    ok = win["ok"]
    if pres_sorted.size == 0 or ok.sum() < 3:
        return float("nan"), p_n, P, F
    return spearman(win["mid"][ok], F[ok])[0], p_n, P, F


def _ci(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size < 20:
        return float("nan"), float("nan")
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


# ── beyond distance (round 2, Reviewer 2 F1) ─────────────────────────────────────────────────────
MIX_LAMBDAS = np.round(np.linspace(0.0, 1.0, 21), 2)


def partial_spearman(a, b, c):
    """Spearman correlation of a and b with c held fixed: Pearson correlation of the residuals of
    rank(a) and rank(b), each regressed on rank(c)."""
    a, b, c = (np.asarray(v, float) for v in (a, b, c))
    ok = np.isfinite(a) & np.isfinite(b) & np.isfinite(c)
    if ok.sum() < 5:
        return float("nan")
    rc = rank_avg(c[ok])
    z = np.column_stack([np.ones(rc.size), rc])

    def res(v):
        r = rank_avg(v)
        return r - z @ np.linalg.lstsq(z, r, rcond=None)[0]
    ea, eb = res(a[ok]), res(b[ok])
    den = math.sqrt(float((ea * ea).sum()) * float((eb * eb).sum()))
    return float((ea * eb).sum() / den) if den > 0 else float("nan")


def stratified_auc_parts(pres, pres_cov, back, back_cov, n_strata=10):
    """Covariate-matched AUC parts: for every outbreak pixel, the share of background pixels in the same
    covariate band (quantile bins of the background covariate) that score lower (ties ½). Their mean is
    the AUC with the covariate held fixed; 0.5 = no information beyond the covariate."""
    edges = np.unique(np.quantile(back_cov, np.linspace(0.0, 1.0, int(n_strata) + 1)))
    inner = edges[1:-1]
    sb = np.searchsorted(inner, back_cov, "right")
    sp = np.searchsorted(inner, pres_cov, "right")
    parts = np.full(pres.size, np.nan)
    for k in np.unique(sp):
        bk = np.sort(back[sb == k])
        if bk.size:
            m = sp == k
            parts[m] = auc_parts(pres[m], bk)
    return parts, int(inner.size + 1)


def mixture_parts(pres_a, pres_b, back_a, back_b, lambdas=MIX_LAMBDAS):
    """Per-outbreak AUC parts of the mixture (1-λ)·u_a + λ·u_b for each λ, with u the percentile rank of a
    value among the background pixels (higher = higher risk for both a and b)."""
    sa, sb_ = np.sort(back_a), np.sort(back_b)
    ua_p, ub_p = auc_parts(pres_a, sa), auc_parts(pres_b, sb_)
    ua_b, ub_b = auc_parts(back_a, sa), auc_parts(back_b, sb_)
    return [auc_parts((1.0 - lam) * ua_p + lam * ub_p, np.sort((1.0 - lam) * ua_b + lam * ub_b))
            for lam in lambdas]


def _loglik(eta, y):
    return float((y * eta - np.logaddexp(0.0, eta)).sum())


def _logistic(X, y, iters=100, offset=None, beta0=None):
    """Maximum-likelihood logistic regression (Newton–Raphson with step halving):
    (beta, covariance, log-likelihood). offset: a fixed term added to the linear predictor."""
    off = np.zeros(y.size) if offset is None else np.asarray(offset, float)
    if beta0 is None:
        pbar = min(max(float(y.mean()), 1e-9), 1 - 1e-9)
        beta = np.zeros(X.shape[1])
        beta[0] = math.log(pbar / (1.0 - pbar))
    else:
        beta = np.array(beta0, float)
    ll = _loglik(X @ beta + off, y)
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-np.clip(X @ beta + off, -35.0, 35.0)))
        H = X.T @ (X * (p * (1.0 - p))[:, None])
        step = np.linalg.solve(H, X.T @ (y - p))
        for _h in range(30):                                  # step halving keeps the likelihood rising
            ll_new = _loglik(X @ (beta + step) + off, y)
            if ll_new >= ll - 1e-10:
                break
            step = step / 2.0
        beta, gain, ll = beta + step, ll_new - ll, ll_new
        if float(np.abs(step).max()) < 1e-8 or abs(gain) < 1e-10:
            break
    p = 1.0 / (1.0 + np.exp(-np.clip(X @ beta + off, -35.0, 35.0)))
    H = X.T @ (X * (p * (1.0 - p))[:, None])
    return beta, np.linalg.inv(H), ll


def _profile_ci(x0, u, y, b_hat, ll1, se, beta0, crit=3.841458820694124, tol=1e-3):
    """95% profile-likelihood interval of the coefficient of u in logit = x0·beta + b·u: the values of b
    at which 2·(ll1 - ll(b)) = 3.84. Consistent with the likelihood-ratio test, unlike a Wald interval
    when u and distance are strongly correlated."""
    def dev(b, start):
        bb, _c, ll = _logistic(x0, y, offset=b * u, beta0=start)
        return 2.0 * (ll1 - ll), bb
    out = []
    for side in (-1.0, 1.0):
        step = float(se) if np.isfinite(se) and se > 0 else 1.0
        inner, start = float(b_hat), np.array(beta0, float)
        outer, found = None, False
        for _k in range(25):                                  # bracket the crossing
            outer = b_hat + side * step
            d_, bb = dev(outer, start)
            if d_ >= crit:
                found = True
                break
            inner, start, step = outer, bb, step * 2.0
        if not found:
            out.append(float("nan"))
            continue
        for _k in range(60):                                  # bisection
            mid = 0.5 * (inner + outer)
            d_, bb = dev(mid, start)
            if d_ >= crit:
                outer = mid
            else:
                inner, start = mid, bb
            if abs(outer - inner) < tol:
                break
        out.append(0.5 * (inner + outer))
    return out[0], out[1]


def _chi2_sf1(x):
    """Upper tail of chi-square with 1 degree of freedom."""
    if not np.isfinite(x):
        return float("nan")
    return math.erfc(math.sqrt(max(x, 0.0) / 2.0))


def beyond_distance(pres, pres_dist, back, back_dist, rng, n_boot=2000, n_strata=10,
                    raw=None, dates=None):
    """What a map adds to distance to the nearest seed point. pres/back: map scores (higher = higher
    risk); pres_dist/back_dist: distance (km) on the same pixels; raw/dates: the raster values and
    outbreak dates of the outbreak pixels, for the timing test (arrival rasters)."""
    out = {}
    if pres.size == 0:
        return out, None
    # 1. distance-matched AUC
    sp, ns = stratified_auc_parts(pres, pres_dist, back, back_dist, n_strata)
    ok = np.isfinite(sp)
    out["auc_within_distance_strata"] = float(sp[ok].mean()) if ok.any() else float("nan")
    out["n_strata"] = ns
    # 2. mixtures with distance (u = percentile rank among the background; distance: closer = higher)
    parts = mixture_parts(-pres_dist, pres, -back_dist, back)
    aucs = np.array([float(p.mean()) for p in parts])
    j = int(np.argmax(aucs))
    out["auc_distance_check"] = float(aucs[0])
    out["best_lambda"] = float(MIX_LAMBDAS[j])
    out["auc_best_mixture"] = float(aucs[j])
    out["mixture_gain"] = float(aucs[j] - aucs[0])
    # 3. likelihood-ratio test, presence/background logistic model on percentile ranks; distance enters
    #    as a cubic so that a curved distance–risk relation is not credited to the map
    sd, sm_ = np.sort(-back_dist), np.sort(back)
    ud = np.concatenate([auc_parts(-pres_dist, sd), auc_parts(-back_dist, sd)])
    um = np.concatenate([auc_parts(pres, sm_), auc_parts(back, sm_)])
    y = np.concatenate([np.ones(pres.size), np.zeros(back.size)])
    x0 = np.column_stack([np.ones(y.size), ud, ud ** 2, ud ** 3])
    try:
        _b0, _c0, ll0 = _logistic(x0, y)
        b1, c1, ll1 = _logistic(np.column_stack([x0, um]), y)
        out["lr_chi2"] = max(0.0, 2.0 * (ll1 - ll0))
        out["lr_p"] = _chi2_sf1(out["lr_chi2"])
        se = math.sqrt(max(float(c1[-1, -1]), 0.0))
        out["coef_map"] = float(b1[-1])
        out["coef_map_lo"], out["coef_map_hi"] = _profile_ci(x0, um, y, float(b1[-1]), ll1, se, b1[:-1])
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        out["lr_chi2"] = out["lr_p"] = float("nan")
    # 4. timing (arrival rasters with dated outbreaks)
    if raw is not None and dates is not None:
        t = np.asarray(dates, float)
        okd = np.isfinite(t) & np.isfinite(raw) & np.isfinite(pres_dist)
        if okd.sum() >= 5:
            r_m, _n = spearman(raw[okd], t[okd])
            r_d, _n = spearman(pres_dist[okd], t[okd])
            out.update(rho_date_map=r_m, rho_date_distance=r_d, rho_date_diff=r_m - r_d,
                       partial_rho_date_map_given_distance=partial_spearman(raw[okd], t[okd], pres_dist[okd]),
                       n_dated=int(okd.sum()))
    # 5. bootstrap over outbreak pixels
    if n_boot and pres.size > 1:
        idx = rng.integers(0, pres.size, size=(int(n_boot), pres.size))
        out["auc_within_distance_strata_lo"], out["auc_within_distance_strata_hi"] = _ci(
            np.nanmean(sp[idx], axis=1))
        P = np.vstack(parts)                                     # λ x outbreaks
        g = np.array([max(0.0, P[:, i].mean(axis=1).max() - P[0, i].mean()) for i in idx])
        out["mixture_gain_lo"], out["mixture_gain_hi"] = _ci(g)
        if "n_dated" in out:
            ii = np.flatnonzero(okd)
            bd, bp, bm = [], [], []
            for _ in range(int(n_boot)):
                k = ii[rng.integers(0, ii.size, ii.size)]
                if np.unique(t[k]).size < 4:
                    continue
                bm.append(spearman(raw[k], t[k])[0])                       # 26 Nov: CI of rho itself (same resamples)
                bd.append(spearman(raw[k], t[k])[0] - spearman(pres_dist[k], t[k])[0])
                bp.append(partial_spearman(raw[k], t[k], pres_dist[k]))
            out["rho_date_map_lo"], out["rho_date_map_hi"] = _ci(bm)
            out["rho_date_diff_lo"], out["rho_date_diff_hi"] = _ci(bd)
            out["partial_rho_date_lo"], out["partial_rho_date_hi"] = _ci(bp)
    return out, aucs


# ── input ─────────────────────────────────────────────────────────────────────────
LON_NAMES = ("longitude", "lon", "long", "lng", "decimallongitude", "x_lon")
LAT_NAMES = ("latitude", "lat", "decimallatitude", "y_lat")
DATE_NAMES = ("date", "report_date", "reportingdate", "observation_date", "observationdate",
              "outbreak_start_date", "start_date", "event_date", "onset")
EMPRES_OBS = "observation.date..dd.mm.yyyy."
EMPRES_REP = "report.date..dd.mm.yyyy."
ID_NAMES = ("event_id", "event.id", "seed_id", "id", "fid")


def parse_date(s):
    """A date from a table cell: ISO (2015-03-30) or day first (30/03/2015, as in the
    EMPRES-i exports); Excel / pandas dates too. None when empty or unreadable."""
    if s is None:
        return None
    if isinstance(s, _dt.datetime):
        return s.date()
    if isinstance(s, _dt.date):
        return s
    if hasattr(s, "to_pydatetime"):                                 # pandas Timestamp
        try:
            return s.to_pydatetime().date()
        except Exception:
            return None
    if isinstance(s, float) and not np.isfinite(s):
        return None
    t = str(s).strip().split(" ")[0].split("T")[0]
    if not t or t.lower() in ("nan", "nat", "none"):
        return None
    for fmt in ("%Y-%m-%d", "%d/%m/%Y", "%Y/%m/%d", "%d.%m.%Y", "%d-%m-%Y", "%Y%m%d", "%Y-%m", "%Y"):
        try:
            return _dt.datetime.strptime(t, fmt).date()
        except ValueError:
            continue
    return None


def srs_from(wkt_or_epsg):
    s = osr.SpatialReference()
    if isinstance(wkt_or_epsg, int):
        s.ImportFromEPSG(wkt_or_epsg)
    else:
        s.ImportFromWkt(wkt_or_epsg)
    try:
        s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except AttributeError:                                          # pragma: no cover
        pass
    return s


def _table_rows(path):
    """(column names, row dicts) from a CSV or an Excel sheet. The header is the first
    line holding both a longitude and a latitude column, so the metadata lines above
    the header of an EMPRES-i export (and the '#' lines of the v9.10 CSVs) are skipped."""
    path = Path(path)

    def header_index(lines):
        for i, cells in enumerate(lines):
            if cells and str(cells[0]).lstrip().startswith("#"):
                continue
            low = {str(c).strip().lower() for c in cells}
            if low & set(LON_NAMES) and low & set(LAT_NAMES):
                return i
        return 0
    if path.suffix.lower() in (".xlsx", ".xls"):
        import pandas as pd                                         # only needed for Excel input
        raw = pd.read_excel(path, header=None, dtype=object)
        h = header_index([list(raw.iloc[i]) for i in range(min(len(raw), 60))])
        cols = [str(v).strip() for v in raw.iloc[h]]
        return cols, [dict(zip(cols, list(raw.iloc[i]))) for i in range(h + 1, len(raw))]
    text = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    h = header_index([next(csv.reader([ln]), []) for ln in text[:60]])
    rd = csv.DictReader(text[h:])
    rows = list(rd)
    return list(rd.fieldnames or []), rows


def _pick(low, wanted, names):
    """Column matching `wanted` (case-insensitive), else the first known alias."""
    if wanted and str(wanted).lower() != "auto":
        return low.get(str(wanted).strip().lower())
    return next((low[a] for a in names if a in low), None)


def _date_columns(low, date_col):
    """Columns to read the date from, in order (the first non-empty one counts)."""
    if date_col is None:
        return [], "no dates"
    if isinstance(date_col, (list, tuple)):
        want = [str(c) for c in date_col]
    elif str(date_col).lower() == "auto":
        if "date" in low:
            return [low["date"]], "'date'"
        if EMPRES_OBS in low or EMPRES_REP in low:
            cols = [low[c] for c in (EMPRES_OBS, EMPRES_REP) if c in low]
            return cols, "EMPRES-i observation date, else report date"
        c = next((low[a] for a in DATE_NAMES if a in low), None)
        return ([c], f"'{c}'") if c else ([], "no date column found")
    else:
        want = str(date_col).split("|")
    cols = [low[w.strip().lower()] for w in want if w.strip().lower() in low]
    return cols, (" else ".join(f"'{c}'" for c in cols) if cols else f"no column '{date_col}'")


def read_points(path, lon_col="auto", lat_col="auto", date_col="auto"):
    """Outbreak points as dicts (lon, lat, date, date_col, id) in WGS84, from a CSV /
    Excel table or any point layer (shapefile, GeoPackage …); plus a note on the dates."""
    path = Path(path)
    out = []
    if path.suffix.lower() in (".csv", ".txt", ".xlsx", ".xls"):
        names, rows = _table_rows(path)
        low = {str(n).strip().lower(): n for n in names}
        cx, cy = _pick(low, lon_col, LON_NAMES), _pick(low, lat_col, LAT_NAMES)
        if cx is None or cy is None:
            raise ValueError(f"{path.name}: no longitude/latitude columns among {names}; "
                             f"pass lon_col= / lat_col=")
        dcols, dnote = _date_columns(low, date_col)
        cid = next((low[a] for a in ID_NAMES if a in low), None)
        for i, row in enumerate(rows, 1):
            try:
                lon, lat = float(row[cx]), float(row[cy])
            except (KeyError, TypeError, ValueError):
                continue
            if not (np.isfinite(lon) and np.isfinite(lat)) or abs(lat) > 90 or abs(lon) > 180:
                continue
            d, used = None, ""
            for c in dcols:
                d = parse_date(row.get(c))
                if d is not None:
                    used = c
                    break
            out.append(dict(lon=lon, lat=lat, date=d, date_col=used,
                            id=str(row.get(cid, "")).strip() if cid else str(i)))
        return out, dnote
    ds = ogr.Open(str(path))
    lyr = ds.GetLayer()
    src = lyr.GetSpatialRef()
    tf = None
    if src is not None:
        try:
            src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:                                      # pragma: no cover
            pass
        tf = osr.CoordinateTransformation(src, srs_from(4326))
    dfield = None
    if date_col and str(date_col).lower() != "auto" and lyr.GetLayerDefn().GetFieldIndex(str(date_col)) >= 0:
        dfield = str(date_col)
    for f in lyr:
        g = f.GetGeometryRef()
        if g is None:
            continue
        g = g.Clone()
        if tf is not None:
            g.Transform(tf)
        pts = ([g.GetGeometryRef(i) for i in range(g.GetGeometryCount())]
               if "MULTI" in g.GetGeometryName().upper() else [g])
        d = parse_date(f.GetField(dfield)) if dfield else None
        for p in pts:
            out.append(dict(lon=p.GetX(), lat=p.GetY(), date=d, date_col=dfield or "", id=str(f.GetFID())))
    ds = None
    return out, (f"'{dfield}'" if dfield else "no dates")


def read_raster(path):
    ds = gdal.Open(str(path))
    b = ds.GetRasterBand(1)
    a = b.ReadAsArray().astype(np.float64)
    nd = b.GetNoDataValue()
    if nd is not None:
        a[np.isclose(a, nd, rtol=0, atol=1e-6)] = np.nan
    a[~np.isfinite(a)] = np.nan
    return a, ds.GetGeoTransform(), ds.GetProjection()


def _unit(lon, lat):
    lo, la = np.radians(np.asarray(lon, float)), np.radians(np.asarray(lat, float))
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def nearest_km(lon, lat, seed_lon, seed_lat, chunk=100_000):
    """Great-circle distance (km) from each point to the nearest seed point."""
    P = _unit(lon, lat)
    S = _unit(seed_lon, seed_lat)
    out = np.empty(P.shape[0])
    for i in range(0, P.shape[0], chunk):
        dots = np.clip((P[i:i + chunk] @ S.T).max(axis=1), -1.0, 1.0)
        out[i:i + chunk] = EARTH_KM * np.arccos(dots)
    return out


def pixel_centres_lonlat(gt, shape, wkt):
    """(lon, lat) of every pixel centre of a raster grid."""
    nr, nc = shape
    cc, rr = np.meshgrid(np.arange(nc) + 0.5, np.arange(nr) + 0.5)
    x = gt[0] + cc * gt[1] + rr * gt[2]
    y = gt[3] + cc * gt[4] + rr * gt[5]
    tf = osr.CoordinateTransformation(srs_from(wkt), srs_from(4326))
    xy = np.column_stack([x.ravel(), y.ravel()])
    res = np.asarray(tf.TransformPoints(xy.tolist()), float)
    return res[:, 0].reshape(shape), res[:, 1].reshape(shape)


def polygon_rings(path):
    """Polygons of a vector file in WGS84 as [(exterior, [holes])] lon/lat arrays."""
    ds = ogr.Open(str(path))
    lyr = ds.GetLayer()
    src = lyr.GetSpatialRef()
    tf = None
    if src is not None:
        try:
            src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:                                      # pragma: no cover
            pass
        if not src.IsSame(srs_from(4326)):
            tf = osr.CoordinateTransformation(src, srs_from(4326))
    out = []
    for f in lyr:
        g0 = f.GetGeometryRef()
        if g0 is None:
            continue
        g = g0.Clone()
        if tf is not None:
            g.Transform(tf)
        polys = ([g.GetGeometryRef(i) for i in range(g.GetGeometryCount())]
                 if g.GetGeometryName().upper() == "MULTIPOLYGON" else [g])
        for pg in polys:
            rings = [np.array(pg.GetGeometryRef(k).GetPoints())[:, :2] for k in range(pg.GetGeometryCount())]
            if rings:
                out.append((rings[0], rings[1:]))
    ds = None
    return out


def polygon_mask(rings, lon_g, lat_g):
    """Grid cells whose centre lies in the polygons, plus the cells touching them
    (one-cell dilation, so outbreaks on the border are kept)."""
    from matplotlib.path import Path as _MPath
    pts = np.column_stack([lon_g.ravel(), lat_g.ravel()])
    ok = np.isfinite(pts).all(axis=1)
    inside = np.zeros(pts.shape[0], bool)
    for ext, holes in rings:
        bb = ok & (pts[:, 0] >= ext[:, 0].min()) & (pts[:, 0] <= ext[:, 0].max()) \
            & (pts[:, 1] >= ext[:, 1].min()) & (pts[:, 1] <= ext[:, 1].max())
        idx = np.flatnonzero(bb)
        if not idx.size:
            continue
        m = _MPath(ext).contains_points(pts[idx])
        for h in holes:
            m &= ~_MPath(h).contains_points(pts[idx])
        inside[idx[m]] = True
    m = inside.reshape(lon_g.shape)
    d = m.copy()
    d[1:, :] |= m[:-1, :]
    d[:-1, :] |= m[1:, :]
    d[:, 1:] |= m[:, :-1]
    d[:, :-1] |= m[:, 1:]
    d[1:, 1:] |= m[:-1, :-1]
    d[:-1, :-1] |= m[1:, 1:]
    d[1:, :-1] |= m[:-1, 1:]
    d[:-1, 1:] |= m[1:, :-1]
    return d


def to_pixels(pts, gt, wkt, shape):
    """Row/col of each point on a raster (-1 when outside)."""
    tf = osr.CoordinateTransformation(srs_from(4326), srs_from(wkt))
    rc = np.full((len(pts), 2), -1, int)
    if not pts:
        return rc
    xy = np.asarray(tf.TransformPoints([(p["lon"], p["lat"]) for p in pts]), float)
    det = gt[1] * gt[5] - gt[2] * gt[4]
    dx, dy = xy[:, 0] - gt[0], xy[:, 1] - gt[3]
    col = np.floor((gt[5] * dx - gt[2] * dy) / det).astype(int)
    row = np.floor((-gt[4] * dx + gt[1] * dy) / det).astype(int)
    ok = (row >= 0) & (row < shape[0]) & (col >= 0) & (col < shape[1])
    rc[ok, 0], rc[ok, 1] = row[ok], col[ok]
    return rc


def write_csv(path, header, rows, meta=()):
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        for m in meta:
            fh.write(f"# {m}\n")
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            vals = []
            for k in header:
                v = r.get(k)
                if v is None or (isinstance(v, float) and not np.isfinite(v)):
                    vals.append("")
                elif isinstance(v, float) and k.startswith("p_"):
                    vals.append(f"{v:.3g}")                         # p-values
                elif isinstance(v, float):
                    vals.append(f"{v:.4f}")
                else:
                    vals.append(v)
            w.writerow(vals)


# ── finding the model outputs ─────────────────────────────────────────────────────
DEFAULT_POINT_FILES = ("EMPRESi_PPR_China_2007_2026.csv", "EMPRESi_PPR_China_2007_2026_clean.csv",
                       "later_outbreaks.csv", "later_outbreaks.xlsx")
RASTER_NOTES = {
    "continental_arrival_month": "simulated arrival month, continental leg (low = early = high risk)",
    "korea_leg_arrival_month": "simulated arrival month, peninsular leg (low = early = high risk)",
    "smcda_all_score": "SMCDA risk score, whole domain (ENTRY measured from this run's seed points)",
    "ppr_suitability_on_grid": "MaxEnt suitability (input map): independent only if the test outbreaks "
                               "were not among its training points",
}


BACKGROUND_FILES = ("china_boundary_ne50m.geojson",)


def find_background(setting):
    """Polygon file that limits the background: 'auto' = china_boundary_ne50m.geojson next to
    this script or in the current folder; None / 'grid' = every valid pixel."""
    if setting is None or str(setting).lower() in ("none", "grid", "all"):
        return None
    if str(setting).lower() == "auto":
        here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
        for base in (Path.cwd(), here):
            for nm in BACKGROUND_FILES:
                if (base / nm).exists():
                    return base / nm
        return None
    return Path(str(setting))


def find_v910_rasters(results_dir=None):
    """The raster folder of a v9.10 run: results_dir itself or results_dir/PPR_v9.10_rasters;
    else results/ppr_v9_10_review*/PPR_v9.10_rasters under the current folder (a hold-out
    run is preferred, because the run of record cannot be tested against its own seeds)."""
    cands = []
    if results_dir:
        rd = Path(results_dir)
        cands += [rd / f"{TAG}_rasters", rd]
    else:
        here = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
        for base in (Path.cwd(), here):
            found = sorted((base / "results").glob(f"ppr_v9_10_review*/{TAG}_rasters"))
            cands += [c for c in found if "holdout" in c.parent.name] + found
    for c in cands:
        if c.is_dir() and any(c.glob(f"{TAG}_*.tif")):
            return c
    return None


def default_rasters(folder):
    folder = Path(folder)
    spec = []
    for name, low in (("continental_arrival_month", True), ("smcda_all_score", False),
                      ("ppr_suitability_on_grid", False)):
        f = folder / f"{TAG}_{name}.tif"
        if f.exists():
            spec.append(str(f) + (":low" if low else ""))
    cls = folder / f"{TAG}_smcda_all_class.tif"
    return spec, (str(cls) if cls.exists() else None)


def read_seeding(run_dir):
    """What seeded the run: PPR_v9.10_seeding.json, else the 'seeding' entry of the run summary."""
    if run_dir is None:
        return {}
    for fn, key in ((f"{TAG}_seeding.json", None), (f"{TAG}_run_summary.json", "seeding")):
        f = Path(run_dir) / fn
        if f.exists():
            try:
                d = json.loads(f.read_text(encoding="utf-8"))
                d = d.get(key, {}) if key else d
                if d:
                    return d
            except Exception:
                pass
    return {}


HOW_TO = """
What this test needs
  1. The rasters of a v9.10 HOLD-OUT run (seeded only with the outbreaks before a cut-off):
         python PPR_model_v9_10.py --preset holdout2015
     -> results/ppr_v9_10_review_holdout2015/PPR_v9.10_rasters/
     (needs PPR_seed_points_dated.csv next to the model script or in the folder you run from)
  2. The outbreak table: the EMPRES-i export as downloaded (EMPRESi_PPR_China_2007_2026.csv),
     or any CSV / Excel file with longitude, latitude (WGS84) and a date.
  3. china_boundary_ne50m.geojson next to this script (limits the background to China).

Then, in a notebook cell:
    from ppr_external_test import run_external_test
    run_external_test(points="EMPRESi_PPR_China_2007_2026.csv",
                      results_dir="results/ppr_v9_10_review_holdout2015")
"""


# ── the test ──────────────────────────────────────────────────────────────────────
def _evaluate(pres, back_sorted, win, rng, n_boot, bench=None):
    """AUC, Boyce and their bootstrap CIs for one raster and one set of outbreak pixels;
    `bench` = (presence scores, sorted background) of the distance benchmark."""
    res = dict(auc=float("nan"), auc_lo=float("nan"), auc_hi=float("nan"), p_auc=float("nan"),
               boyce=float("nan"), boyce_lo=float("nan"), boyce_hi=float("nan"),
               auc_distance=float("nan"), boyce_distance=float("nan"),
               auc_diff=float("nan"), auc_diff_lo=float("nan"), auc_diff_hi=float("nan"))
    if pres.size == 0:
        return res, None
    u = auc_parts(pres, back_sorted)
    res["auc"] = float(u.mean())
    res["p_auc"] = mann_whitney_p(pres, back_sorted)
    curve = None
    if win is not None:
        res["boyce"], p_n, P, F = boyce_index(win, np.sort(pres))
        curve = dict(p_n=p_n, P=P, F=F)
    idx = rng.integers(0, pres.size, size=(int(n_boot), pres.size)) if n_boot and pres.size > 1 else None
    if idx is not None:
        res["auc_lo"], res["auc_hi"] = _ci(u[idx].mean(axis=1))
        if win is not None:
            res["boyce_lo"], res["boyce_hi"] = _ci([boyce_index(win, np.sort(pres[i]))[0] for i in idx])
    if bench is not None:
        bp, bb_sorted, bwin = bench
        ub = auc_parts(bp, bb_sorted)
        res["auc_distance"] = float(ub.mean())
        if bwin is not None:
            res["boyce_distance"] = boyce_index(bwin, np.sort(bp))[0]
        res["auc_diff"] = res["auc"] - res["auc_distance"]
        if idx is not None:
            res["auc_diff_lo"], res["auc_diff_hi"] = _ci((u[idx] - ub[idx]).mean(axis=1))
    return res, curve


def run_external_test(points=None, results_dir=None, rasters=None, classes=None,
                      class_labels="Very low,Low,Moderate,High,Very high",
                      lon_col="auto", lat_col="auto", date_col="auto",
                      after="auto", before=None, exclude="auto", exclude_km=10.0,
                      sensitivity_km=(0, 25), all_points=False, out=None,
                      n_boot=2000, seed=20261109, background="auto", n_strata=10):
    """Score the risk rasters of a model run against later outbreaks. Returns the summary
    rows (list of dicts), or None when an input is missing (a message says what to supply).

    points        outbreak table (CSV / Excel, EMPRES-i exports as they are) or point layer
    results_dir   folder of the run (e.g. results/ppr_v9_10_review_holdout2015)
    rasters       list of GeoTIFFs instead of the run's own; add ':low' when low = high risk
    after         keep outbreaks dated on/after this day ('auto' = the run's seeding cut-off)
    before        keep outbreaks dated before this day (optional)
    exclude       seed points to keep away from ('auto' = the run's PPR_v9.10_seed_points_used.csv)
    exclude_km    leave out outbreaks within this distance of a seed point (main result)
    sensitivity_km  further exclusion distances reported in the summary
    background    polygon file limiting the background pixels to where the outbreaks were
                  recorded ('auto' = china_boundary_ne50m.geojson; None = the whole grid,
                  which is also reported as a sensitivity row)
    n_strata      distance bands of the distance-matched AUC (quantiles of the background)"""
    # ---- the run and its rasters --------------------------------------------------------
    run_dir = Path(results_dir) if results_dir else None
    if not rasters:
        folder = find_v910_rasters(results_dir)
        if folder is None:
            print(f"No v9.10 rasters found{' in ' + str(results_dir) if results_dir else ''}." + HOW_TO)
            return None
        run_dir = folder.parent if folder.name == f"{TAG}_rasters" else folder
        rasters, cls_default = default_rasters(folder)
        classes = classes or cls_default
    seeding = read_seeding(run_dir)
    if points is None:
        points = next((p for p in DEFAULT_POINT_FILES if Path(p).exists()), None)
        if points is None:
            print("No outbreak table given (points=...)." + HOW_TO)
            return None
    if not Path(str(points)).exists():
        print(f"Outbreak table not found: {points} (relative to {Path.cwd()})" + HOW_TO)
        return None
    if str(after).lower() == "auto":
        after = seeding.get("before")
    if str(exclude).lower() == "auto":
        sl = (run_dir / f"{TAG}_seed_points_used.csv") if run_dir else None
        exclude = str(sl) if (sl is not None and sl.exists()) else None
    missing = [s[:-4] if s.endswith(":low") else s for s in map(str, rasters)
               if not Path(s[:-4] if s.endswith(":low") else s).exists()]
    if classes and not Path(str(classes)).exists():
        missing.append(str(classes))
    if exclude and not Path(str(exclude)).exists():
        missing.append(str(exclude))
    bg_path = find_background(background)
    if bg_path is not None and not bg_path.exists():
        missing.append(str(bg_path))
    if missing:
        print("File(s) not found — check the path(s), relative to " + str(Path.cwd()) + ":\n   "
              + "\n   ".join(missing))
        return None
    out = Path(out) if out else ((run_dir / "external_test") if run_dir else Path("external_test"))
    out.mkdir(parents=True, exist_ok=True)
    rng_seed = int(seed)

    # ---- run information ---------------------------------------------------------------
    print(f"run folder: {run_dir if run_dir else '(rasters given directly)'}")
    if seeding.get("source") == "table" and seeding.get("before"):
        print(f"   hold-out run: seeded with {seeding.get('seed_rows', '?')} points "
              f"({seeding.get('seed_cells', '?')} cells) dated before {seeding['before']}")
    elif seeding:
        print(f"   seeded from {Path(str(seeding.get('shapefile') or seeding.get('table'))).name} "
              f"({seeding.get('seed_rows', '?')} points, no cut-off) — the run of record")
    else:
        print("   (no seeding record in this folder — runs made before 9 Nov 2026 do not write one)")

    # ---- outbreaks: date window, distance to the seed points --------------------------------
    pts, dnote = read_points(points, lon_col, lat_col, date_col)
    d_after, d_before = parse_date(after), parse_date(before)
    print(f"outbreak table: {Path(str(points)).name} — {len(pts)} rows with coordinates; dates: {dnote}")
    for p in pts:
        p["status"] = "used"
        if (d_after or d_before) and p["date"] is None:
            p["status"] = "no date"
        elif d_after and p["date"] < d_after:
            p["status"] = f"before {d_after.isoformat()}"
        elif d_before and p["date"] >= d_before:
            p["status"] = f"on/after {d_before.isoformat()}"
    window = [p for p in pts if p["status"] == "used"]
    seeds = []
    if exclude:
        seeds, _sn = read_points(exclude, date_col=None)
    if seeds and pts:
        dk = nearest_km([p["lon"] for p in pts], [p["lat"] for p in pts],
                        [s["lon"] for s in seeds], [s["lat"] for s in seeds])
        for p, d in zip(pts, dk):
            p["seed_km"] = float(d)
    else:
        for p in pts:
            p["seed_km"] = float("nan")
    wtxt = (f"{d_after.isoformat() if d_after else '…'} to {d_before.isoformat() if d_before else 'end'}"
            if (d_after or d_before) else "all dates")
    print(f"window {wtxt}: {len(window)} outbreaks"
          + (f" ({sum(1 for p in window if p['date_col'] and p['date_col'].lower() == EMPRES_REP)} dated by "
             f"the report date)" if any(p["date_col"].lower() == EMPRES_REP for p in window) else ""))
    if not (d_after or d_before):
        print("   ⚠ no date window: every outbreak in the table is scored")
    if seeds:
        n_same = sum(1 for p in window if p["seed_km"] <= 1.0)
        print(f"   seed points: {Path(str(exclude)).name} ({len(seeds)} points); "
              f"{n_same} of the {len(window)} outbreaks {'lies' if n_same == 1 else 'lie'} within 1 km of one")
        if n_same and n_same >= 0.5 * max(len(window), 1):
            print("   ⚠⚠ most of these outbreaks seeded this run — the test is NOT independent. "
                  "Use a hold-out run (--preset holdout2015) and its folder as results_dir.")
        elif n_same and not float(exclude_km) > 0:
            print("   ⚠ some outbreaks coincide with seed points and exclude_km = 0")
    else:
        print("   ⚠ seed points of this run unknown (exclude=None): independence cannot be checked. "
              "If the outbreaks seeded the run, the test is circular.")
    kms = [float(exclude_km)] + [float(k) for k in (sensitivity_km or ()) if float(k) != float(exclude_km)]
    if not seeds:
        kms = [0.0]

    def test_set(km):
        if km > 0 and seeds:
            return [p for p in window if not (p["seed_km"] <= km)]
        return list(window)
    for km in kms:
        print(f"   exclusion {km:g} km: {len(test_set(km))} outbreaks")
    for p in window:
        if seeds and kms[0] > 0 and p["seed_km"] <= kms[0]:
            p["status"] = f"within {kms[0]:g} km of a seed point"
    # earliest first, so a pixel counted once keeps its first outbreak
    order = sorted(range(len(pts)), key=lambda i: (pts[i]["date"] is None, pts[i]["date"] or _dt.date.max))
    if not seeding and not seeds:
        print("   ⚠ this folder has no seeding record. If it is the run of record, its seeds include every "
              "EMPRES-i outbreak in China 2007–2025 and the test is circular — use the hold-out run.")

    # ---- rasters -------------------------------------------------------------------------
    rings = polygon_rings(bg_path) if bg_path is not None else None
    main_scope = f"inside {bg_path.name}" if rings else "whole grid"
    if rings:
        print(f"background: grid cells inside or touching {bg_path.name} (whole grid as a sensitivity row)")
    else:
        print("background: every valid pixel of each raster")
    grid_cache = {}

    def grid_info(gt, shape, wkt):
        """Per raster grid: distance to the nearest seed point and the background mask."""
        key = (tuple(gt), tuple(shape), wkt)
        if key not in grid_cache:
            gi = {}
            if seeds or rings:
                lon_g, lat_g = pixel_centres_lonlat(gt, shape, wkt)
                if seeds:
                    d = nearest_km(lon_g.ravel(), lat_g.ravel(),
                                   [s_["lon"] for s_ in seeds], [s_["lat"] for s_ in seeds]).reshape(shape)
                    if not np.isfinite(d).all():                    # pixels the CRS cannot map: farthest
                        d = np.where(np.isfinite(d), d, np.nanmax(d) if np.isfinite(d).any() else 0.0)
                    gi["dist"] = d
                if rings:
                    gi["mask"] = polygon_mask(rings, lon_g, lat_g)
            grid_cache[key] = gi
        return grid_cache[key]

    summary, curves, values, in_scope = [], [], {}, {}
    bd_rows, mix_rows = [], []
    for spec in map(str, rasters):
        low = spec.endswith(":low")
        path = Path(spec[:-4] if low else spec)
        name = path.name
        arr, gt, wkt = read_raster(path)
        gi = grid_info(gt, arr.shape, wkt)
        sgn = -1.0 if low else 1.0
        rc = to_pixels(pts, gt, wkt, arr.shape)
        val = np.full(len(pts), np.nan)
        inside = rc[:, 0] >= 0
        val[inside] = arr[rc[inside, 0], rc[inside, 1]]
        values[name] = val
        dist = gi.get("dist")
        note = next((v for k, v in RASTER_NOTES.items() if k in name), "")
        scopes = [(main_scope, np.isfinite(arr) & gi["mask"] if rings else np.isfinite(arr), kms)]
        if rings:
            scopes.append(("whole grid", np.isfinite(arr), kms[:1]))
        for scope, valid, kms_here in scopes:
            ok_pt = np.zeros(len(pts), bool)
            ok_pt[inside] = valid[rc[inside, 0], rc[inside, 1]]
            if scope == main_scope:
                in_scope[name] = ok_pt
            back_sorted = np.sort(sgn * arr[valid])
            if back_sorted.size < 10:
                print(f"{name}: fewer than 10 background pixels in '{scope}' — skipped")
                continue
            win = boyce_windows(back_sorted)
            bench_back = np.sort(-dist[valid]) if dist is not None else None
            bench_win = boyce_windows(bench_back) if bench_back is not None else None
            for km in kms_here:
                ids = {id(p) for p in test_set(km)}
                sel, seen = [], set()
                for i in order:
                    if id(pts[i]) not in ids or not ok_pt[i]:
                        continue
                    pix = (int(rc[i, 0]), int(rc[i, 1]))
                    if not all_points and pix in seen:
                        continue
                    seen.add(pix)
                    sel.append(i)
                sel = np.array(sel, int)
                pres = sgn * val[sel] if sel.size else np.array([])
                bench = None
                if dist is not None and sel.size:
                    bench = (-dist[rc[sel, 0], rc[sel, 1]], bench_back, bench_win)
                res, curve = _evaluate(pres, back_sorted, win, np.random.default_rng(rng_seed), n_boot, bench)
                if dist is not None and sel.size:                     # round 2: value beyond distance
                    _dts = (np.array([pts[i]["date"].toordinal() if pts[i]["date"] else np.nan for i in sel], float)
                            if low else None)
                    bd, aucs = beyond_distance(pres, dist[rc[sel, 0], rc[sel, 1]], sgn * arr[valid], dist[valid],
                                               np.random.default_rng(rng_seed + 1), n_boot, n_strata,
                                               raw=(val[sel] if low else None), dates=_dts)
                    bd_rows.append(dict(background=scope, exclude_km=km, raster=name, presence_pixels=int(sel.size),
                                        background_pixels=int(back_sorted.size), auc=res["auc"],
                                        auc_distance=res["auc_distance"], auc_diff=res["auc_diff"], **bd))
                    if aucs is not None and scope == main_scope and km == kms[0]:
                        mix_rows += [dict(raster=name, weight_of_map=float(l_), auc=float(a_))
                                     for l_, a_ in zip(MIX_LAMBDAS, aucs)]
                rho_d, n_d = float("nan"), 0
                if low and sel.size:
                    dd = np.array([pts[i]["date"].toordinal() if pts[i]["date"] else np.nan for i in sel], float)
                    rho_d, n_d = spearman(val[sel], dd)
                summary.append(dict(
                    background=scope, exclude_km=km, raster=name, note=note,
                    direction="low = high risk" if low else "high = high risk", outbreaks=len(ids),
                    outbreaks_on_raster=int(sum(1 for i in order if id(pts[i]) in ids and ok_pt[i])),
                    presence_pixels=int(sel.size), background_pixels=int(back_sorted.size),
                    rho_arrival_vs_date=rho_d, p_rho=spearman_p(rho_d, n_d), n_dated=n_d, **res))
                main = (scope == main_scope and km == kms[0])
                if main and curve is not None and win is not None:
                    lo_, hi_ = float(back_sorted[0]), float(back_sorted[-1])
                    for j in np.flatnonzero(win["ok"]):
                        a_, b_ = sgn * win["starts"][j], sgn * win["ends"][j]          # raster units
                        curves.append(dict(raster=name, window_low=min(a_, b_), window_high=max(a_, b_),
                                           position_pct=100.0 * (win["mid"][j] - lo_) / (hi_ - lo_),
                                           presences=int(curve["p_n"][j]), background=int(win["e_n"][j]),
                                           P=float(curve["P"][j]), E=float(win["E"][j]), F=float(curve["F"][j])))
                if main:
                    ci = (lambda a, b: f" ({a:.2f}–{b:.2f})" if np.isfinite(a) else "")
                    txt = (f"{name}: {sel.size} outbreak pixels — AUC {res['auc']:.3f}"
                           f"{ci(res['auc_lo'], res['auc_hi'])}, p {res['p_auc']:.2g}; Boyce {res['boyce']:+.2f}"
                           f"{ci(res['boyce_lo'], res['boyce_hi'])}")
                    if np.isfinite(res["auc_distance"]):
                        txt += (f"; distance-to-seed AUC {res['auc_distance']:.3f}, difference "
                                f"{res['auc_diff']:+.3f}{ci(res['auc_diff_lo'], res['auc_diff_hi'])}")
                    if n_d:
                        txt += f"; rho(arrival, date) {rho_d:+.2f} (n={n_d}, p {spearman_p(rho_d, n_d):.2g})"
                    print(txt)
    for i, p in enumerate(pts):
        if p["status"] != "used":
            continue
        if values and all(not np.isfinite(values[vn][i]) for vn in values):
            p["status"] = "outside the rasters / no data"
        elif rings and in_scope and not any(in_scope[vn][i] for vn in in_scope):
            p["status"] = "outside the background area"
    extra = [r for r in summary if not (r["background"] == main_scope and r["exclude_km"] == kms[0])]
    if extra:
        print("sensitivity (AUC):")
        for name in dict.fromkeys(r["raster"] for r in summary):
            print("   " + name + ": " + ", ".join(
                f"{r['exclude_km']:g} km{'' if r['background'] == main_scope else ', whole grid'} "
                f"{r['auc']:.3f} (n={r['presence_pixels']})"
                for r in summary if r["raster"] == name))

    # ---- value beyond distance (main exclusion distance and background) ----------------------
    if bd_rows:
        print("beyond distance (main test set):")
        for r_ in bd_rows:
            if not (r_["background"] == main_scope and r_["exclude_km"] == kms[0]):
                continue
            ci = (lambda a, b: f" ({a:.2f}–{b:.2f})" if np.isfinite(a) else "")
            txt = (f"   {r_['raster']}: AUC within distance bands {r_['auc_within_distance_strata']:.3f}"
                   f"{ci(r_.get('auc_within_distance_strata_lo', np.nan), r_.get('auc_within_distance_strata_hi', np.nan))}"
                   f"; map added to distance: LR chi2 {r_.get('lr_chi2', np.nan):.2f}, p {r_.get('lr_p', np.nan):.2g}, "
                   f"coefficient {r_.get('coef_map', np.nan):+.2f}"
                   f"{ci(r_.get('coef_map_lo', np.nan), r_.get('coef_map_hi', np.nan))}"
                   f"{' (wrong direction)' if r_.get('coef_map', 0) < 0 else ''}"
                   f"; best mixture weight {r_['best_lambda']:.2f}, AUC gain {r_['mixture_gain']:.3f}"
                   f"{ci(r_.get('mixture_gain_lo', np.nan), r_.get('mixture_gain_hi', np.nan))}")
            if "n_dated" in r_:
                txt += (f"; dates (n={r_['n_dated']}): rho(arrival) {r_['rho_date_map']:+.2f} vs rho(distance) "
                        f"{r_['rho_date_distance']:+.2f}, difference {r_['rho_date_diff']:+.2f}"
                        f"{ci(r_.get('rho_date_diff_lo', np.nan), r_.get('rho_date_diff_hi', np.nan))}"
                        f", partial given distance {r_['partial_rho_date_map_given_distance']:+.2f}"
                        f"{ci(r_.get('partial_rho_date_lo', np.nan), r_.get('partial_rho_date_hi', np.nan))}")
            print(txt)

    # ---- risk classes (main exclusion distance and background) ------------------------------
    crow = []
    if classes:
        cls, gtc, wktc = read_raster(classes)
        cmask = np.isfinite(cls) & (cls >= 0)
        if rings:
            cmask &= grid_info(gtc, cls.shape, wktc)["mask"]
        labels = [s_.strip() for s_ in str(class_labels).split(",")]
        rc = to_pixels(pts, gtc, wktc, cls.shape)
        ids = {id(p) for p in test_set(kms[0])}
        cv, seen = [], set()
        for i in order:
            if id(pts[i]) not in ids or rc[i, 0] < 0 or not cmask[rc[i, 0], rc[i, 1]]:
                continue
            pix = (int(rc[i, 0]), int(rc[i, 1]))
            if not all_points and pix in seen:
                continue
            seen.add(pix)
            cv.append(cls[rc[i, 0], rc[i, 1]])
        cv = np.array(cv, float)
        allc = cls[cmask]
        for k in range(int(np.nanmax(allc)) + 1 if allc.size else 0):
            p_ = float((cv == k).mean()) if cv.size else float("nan")
            e_ = float((allc == k).mean()) if allc.size else float("nan")
            crow.append(dict(cls=k, label=labels[k] if k < len(labels) else f"class {k}",
                             outbreaks=int((cv == k).sum()), pct_outbreaks=100 * p_, pct_area=100 * e_,
                             ratio_P_E=(p_ / e_ if e_ > 0 else float("nan"))))
        print(f"risk classes ({Path(str(classes)).name}; exclusion {kms[0]:g} km; % of area {main_scope}):")
        for r_ in crow:
            print(f"   {r_['label']:<10s} {r_['outbreaks']:4d} outbreaks ({r_['pct_outbreaks']:5.1f}%) on "
                  f"{r_['pct_area']:5.1f}% of the area, P/E {r_['ratio_P_E']:.2f}")

    # ---- files -----------------------------------------------------------------------------
    meta = [f"{TAG} external test — outbreaks: {Path(str(points)).name} (dates: {dnote}); window {wtxt}",
            f"run: {run_dir}; seeding: " + (f"hold-out, seeds dated before {seeding.get('before')}"
                                            if seeding.get("before") else (seeding.get("source") or "unknown")),
            f"seed points kept away from: {Path(str(exclude)).name if exclude else 'none'}; "
            f"main exclusion {kms[0]:g} km; background {main_scope}; each outbreak pixel counted once "
            f"(earliest)" + ("" if not all_points else " — all reports counted"),
            f"AUC: outbreak pixels vs all valid pixels; 95% CI = bootstrap over outbreaks ({n_boot} resamples, "
            f"seed {rng_seed}); p_auc one-sided Mann–Whitney; auc_distance = AUC of minus the distance "
            f"to the nearest seed point on the same pixels; auc_diff = map − distance (paired bootstrap CI)"]
    write_csv(out / "external_test_summary.csv",
              ["background", "exclude_km", "raster", "note", "direction", "outbreaks", "outbreaks_on_raster",
               "presence_pixels", "background_pixels", "auc", "auc_lo", "auc_hi", "p_auc", "boyce",
               "boyce_lo", "boyce_hi", "auc_distance", "boyce_distance", "auc_diff", "auc_diff_lo",
               "auc_diff_hi", "rho_arrival_vs_date", "p_rho", "n_dated"], summary, meta)
    if bd_rows:
        write_csv(out / "external_test_beyond_distance.csv",
                  ["background", "exclude_km", "raster", "presence_pixels", "background_pixels", "auc",
                   "auc_distance", "auc_diff", "auc_within_distance_strata", "auc_within_distance_strata_lo",
                   "auc_within_distance_strata_hi", "n_strata", "lr_chi2", "lr_p", "coef_map", "coef_map_lo",
                   "coef_map_hi", "best_lambda", "auc_best_mixture", "mixture_gain", "mixture_gain_lo",
                   "mixture_gain_hi", "rho_date_map", "rho_date_map_lo", "rho_date_map_hi", "rho_date_distance", "rho_date_diff", "rho_date_diff_lo",
                   "rho_date_diff_hi", "partial_rho_date_map_given_distance", "partial_rho_date_lo",
                   "partial_rho_date_hi", "n_dated"], bd_rows,
                  meta[:3] + [
                      "Value beyond distance to the nearest seed point (round 2, Reviewer 2 F1). "
                      "auc_within_distance_strata: each outbreak pixel compared only with background pixels in "
                      f"its own distance band ({n_strata} quantile bands of the background); 0.5 = nothing beyond "
                      "distance. lr_chi2/lr_p: likelihood-ratio test of the map added to distance in a "
                      "presence/background logistic model (percentile ranks; distance as a cubic); coef_map = "
                      "log-odds across the map's range (negative = wrong direction), 95% profile-likelihood "
                      "interval. best_lambda/mixture_gain: mixture (1-λ)·distance + λ·map with the highest "
                      "AUC and its gain over distance (in-sample upper bound). rho_date_*: Spearman rho of the "
                      "outbreak date with the arrival month and with distance; partial = given distance. "
                      f"Other 95% CIs: bootstrap over outbreak pixels ({n_boot} resamples, seed {rng_seed + 1})."])
        write_csv(out / "external_test_mixture.csv", ["raster", "weight_of_map", "auc"], mix_rows,
                  meta[:3] + ["AUC of the mixture (1-w)·distance + w·map (percentile ranks among the background); "
                              "w = 0 is distance alone, w = 1 the map alone; main test set"])
    write_csv(out / "external_test_boyce_curve.csv",
              ["raster", "window_low", "window_high", "position_pct", "presences", "background", "P", "E", "F"],
              curves, meta[:1] + ["window_low/high in raster units; position_pct = window midpoint in the "
                                  "risk range (0 = lowest risk, 100 = highest; arrival months reversed)"])
    if crow:
        write_csv(out / "external_test_classes.csv",
                  ["cls", "label", "outbreaks", "pct_outbreaks", "pct_area", "ratio_P_E"], crow, meta[:3])
    vnames = list(values)
    prow = []
    for i in order:
        p = pts[i]
        r_ = dict(id=p["id"], date=p["date"].isoformat() if p["date"] else "", date_column=p["date_col"],
                  longitude=f"{p['lon']:.6f}", latitude=f"{p['lat']:.6f}",
                  nearest_seed_km=p["seed_km"], status=p["status"])
        for vn in vnames:
            r_[vn] = float(values[vn][i])
        prow.append(r_)
    write_csv(out / "external_test_points.csv",
              ["id", "date", "date_column", "longitude", "latitude", "nearest_seed_km", "status"] + vnames,
              prow, meta[:3])
    (out / "external_test_settings.json").write_text(json.dumps(dict(
        points=str(points), dates=dnote, results_dir=str(run_dir) if run_dir else None, seeding=seeding,
        rasters=list(map(str, rasters)), classes=classes, after=str(d_after) if d_after else None,
        before=str(d_before) if d_before else None, exclude=exclude, exclude_km=kms[0], sensitivity_km=kms[1:],
        background=str(bg_path) if bg_path else None,
        all_points=bool(all_points), n_boot=int(n_boot), seed=rng_seed,
        created=_dt.datetime.now().strftime("%Y-%m-%d %H:%M")), indent=2, default=str), encoding="utf-8")
    try:
        # object-oriented matplotlib: a notebook's plotting backend is left alone
        from matplotlib.figure import Figure
        main_rows = [r for r in summary if r["exclude_km"] == kms[0] and r["background"] == main_scope]
        fig = Figure(figsize=(12 if crow else 6.5, 4.8))
        ax = fig.subplots(1, 2 if crow else 1, squeeze=False)
        for r_ in main_rows:
            rr = [b for b in curves if b["raster"] == r_["raster"]]
            if rr:
                lab = r_["raster"].replace(f"{TAG}_", "").replace(".tif", "")
                ax[0][0].plot([b["position_pct"] for b in rr], [b["F"] for b in rr], lw=1.8,
                              label=f"{lab} (Boyce {r_['boyce']:+.2f})")
        ax[0][0].axhline(1.0, color="#777", lw=0.8, ls="--")
        ax[0][0].set_yscale("symlog", linthresh=1.0)                # F spans 0 to ≫1
        ax[0][0].set_xlabel("Position in the risk range (%; 100 = highest risk)")
        ax[0][0].set_ylabel("Predicted / expected ratio")
        ax[0][0].legend(fontsize=8, frameon=False)
        ax[0][0].text(-0.1, 1.02, "A", transform=ax[0][0].transAxes, fontsize=13, fontweight="bold")
        if crow:
            x = np.arange(len(crow))
            ax[0][1].bar(x - 0.2, [c_["pct_area"] for c_ in crow], 0.4, color="#bdbdbd", label="% of area")
            ax[0][1].bar(x + 0.2, [c_["pct_outbreaks"] for c_ in crow], 0.4, color="#d6604d",
                         label="% of outbreaks")
            ax[0][1].set_xticks(x)
            ax[0][1].set_xticklabels([c_["label"] for c_ in crow], fontsize=9)
            ax[0][1].set_ylabel("%")
            ax[0][1].legend(fontsize=9, frameon=False)
            ax[0][1].text(-0.1, 1.02, "B", transform=ax[0][1].transAxes, fontsize=13, fontweight="bold")
        fig.tight_layout()
        fig.savefig(out / "external_test.png", dpi=200, bbox_inches="tight")
    except Exception as e:                                          # pragma: no cover
        print(f"(figure skipped: {e})")
    print(f"written to {out.resolve()}")
    return summary


# ── command line (also safe inside Jupyter: the kernel's own arguments are ignored) ──
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--points", default=None, help="outbreak table (CSV / Excel, EMPRES-i as exported) or point layer")
    ap.add_argument("--results-dir", default=None, help="run folder, e.g. results/ppr_v9_10_review_holdout2015")
    ap.add_argument("--raster", action="append", default=None,
                    help="risk raster (repeatable) instead of the run's own; append ':low' when low = high risk")
    ap.add_argument("--classes", default=None, help="class raster (0 = lowest class, -1 = no data)")
    ap.add_argument("--class-labels", default="Very low,Low,Moderate,High,Very high")
    ap.add_argument("--lon-col", default="auto")
    ap.add_argument("--lat-col", default="auto")
    ap.add_argument("--date-col", default="auto", help="'auto', a column name, or 'a|b' (first non-empty)")
    ap.add_argument("--after", default="auto", help="keep outbreaks on/after this day (auto = the run's cut-off)")
    ap.add_argument("--before", default=None, help="keep outbreaks before this day")
    ap.add_argument("--exclude", default="auto", help="seed points to keep away from (auto = the run's seed list)")
    ap.add_argument("--exclude-km", type=float, default=10.0)
    ap.add_argument("--sensitivity-km", default="0,25", help="further exclusion distances, comma-separated")
    ap.add_argument("--all-points", action="store_true", help="count every report, not every pixel once")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--background", default="auto",
                    help="polygon file limiting the background (auto = china_boundary_ne50m.geojson; 'grid' = all pixels)")
    ap.add_argument("--out", default=None, help="output folder (default: <run folder>/external_test)")
    ap.add_argument("--n-strata", type=int, default=10, help="distance bands of the distance-matched AUC")
    a, _unknown = ap.parse_known_args(argv)
    sens = tuple(float(v) for v in str(a.sensitivity_km).split(",") if v.strip())
    return run_external_test(points=a.points, results_dir=a.results_dir, rasters=a.raster, classes=a.classes,
                             class_labels=a.class_labels, lon_col=a.lon_col, lat_col=a.lat_col,
                             date_col=a.date_col, after=a.after, before=a.before, exclude=a.exclude,
                             exclude_km=a.exclude_km, sensitivity_km=sens, all_points=a.all_points,
                             out=a.out, n_boot=a.n_boot, background=a.background, n_strata=a.n_strata)


if __name__ == "__main__":
    main()
