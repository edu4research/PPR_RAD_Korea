"""Audit of the external-test AUC (Reviewer 2, point 3; 27 Nov 2026).

Recomputes the AUC of the hold-out run's simulated arrival month from the raster and the list of test records
(Supplementary Table S13), prints every count that enters it, and adds one sensitivity check: the AUC among
reached cells only (outbreak cells and background cells that the 72-month hold-out run reached), which does not
depend on how never-reached cells (coded month 73, tied) are ranked.

    conda activate rapids_env
    cd /mnt/c/Users/USER/Downloads/PPR_Results/v9_10_code
    python ppr_auc_audit.py ~/stomoxys_model/results/ppr_v9_10_review_holdout2015_pcfloor

Needs ppr_external_test.py next to it; the record list and the China boundary are found here or in ../data.
Writes <run folder>/external_test/auc_audit.csv.
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ppr_external_test as X  # noqa: E402


def auc(pres, back):
    return float(X.auc_parts(np.asarray(pres, float), np.sort(np.asarray(back, float))).mean())


def _find(name):
    for base in (HERE, HERE.parent / "data", Path.cwd()):
        if (base / name).exists():
            return base / name
    sys.exit(f"{name} not found next to this script, in ../data or in the current folder")


def main(run_dir, records=None, boundary=None, out=None):
    records = records or _find("TableS13_external_test_records.csv")
    boundary = boundary or _find("china_boundary_ne50m.geojson")
    out = Path(out) if out else Path(run_dir).expanduser() / "external_test" / "auc_audit.csv"
    run_dir = Path(run_dir).expanduser()
    hits = sorted(run_dir.rglob("*continental_arrival_month.tif")) if run_dir.is_dir() else []
    if not hits:                                            # look for the hold-out run anywhere under results/
        base = run_dir if run_dir.is_dir() else Path("~/stomoxys_model/results").expanduser()
        hits = sorted(p for p in base.parent.rglob("*continental_arrival_month.tif")
                      if "holdout2015" in str(p) and "pcfloor" in str(p))
    if not hits:
        sys.exit(f"no *continental_arrival_month.tif under {run_dir} - run: "
                 f"find ~/stomoxys_model -name '*continental_arrival_month.tif'")
    print("arrival raster:", hits[0])
    arr, gt, wkt = X.read_raster(hits[0])
    last = np.nanmax(arr)                                   # 73 = never reached in the 72-month run
    lon_g, lat_g = X.pixel_centres_lonlat(gt, arr.shape, wkt)
    china = X.polygon_mask(X.polygon_rings(boundary), lon_g, lat_g) & np.isfinite(arr)
    back = arr[china]
    with open(records, encoding="utf-8-sig") as f:
        rec = list(csv.DictReader(f))
    rows = []
    for km in (0.0, 10.0, 25.0):
        kept = [r for r in rec if km == 0 or float(r["nearest_seed_km"]) > km]
        pts = [{"lon": float(r["longitude"]), "lat": float(r["latitude"])} for r in kept]
        rc = X.to_pixels(pts, gt, wkt, arr.shape)
        cells = sorted({(int(a), int(b)) for a, b in rc if a >= 0 and china[a, b]})
        pv = np.array([arr[a, b] for a, b in cells])
        reached_p, reached_b = pv[pv < last], back[back < last]
        row = dict(exclusion_km=km, records=len(kept),
                   locations=len({(r["longitude"], r["latitude"]) for r in kept}),
                   outbreak_cells=len(cells), outbreak_cells_not_reached=int((pv >= last).sum()),
                   background_cells=int(back.size), background_not_reached=int((back >= last).sum()),
                   pct_background_not_reached=round(100.0 * (back >= last).mean(), 1),
                   auc=round(auc(-pv, -back), 4),
                   auc_reached_only=round(auc(-reached_p, -reached_b), 4) if reached_p.size else "",
                   outbreak_cells_reached=int(reached_p.size), background_reached=int(reached_b.size))
        rows.append(row)
        print(row)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("written", out)


if __name__ == "__main__":
    main(sys.argv[1])
