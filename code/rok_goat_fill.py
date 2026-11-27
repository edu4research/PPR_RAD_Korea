"""Fill the ROK part of the GLW4 goat layer with the national goat register.

The model's goat layer (GLW4GTS_1km_5179.tif, head per 1-km pixel = head/km2) was clipped
to China and the DPRK, so every ROK pixel is empty and the model reads it as zero. This
script writes a merged copy:
  - ROK territory (sigungu polygons of sig.shp): goat density from the register, head/km2;
    pixels without goats are 0, not NaN;
  - everywhere else: the GLW4 values, unchanged.

Register input, detected from its columns:
  (a) farm points: x/y in EPSG:5179 or lon/lat in WGS84, plus a head count -> heads are
      added to the 1-km pixel holding the farm;
  (b) farms or sigungu rows with a 5-digit sigungu code (SIG_CD, or the first five digits
      of a 10-digit legal-dong code) and a head count -> heads per sigungu spread evenly
      over the sigungu's land pixels.

python rok_goat_fill.py GLW.tif sig.shp REGISTER.(csv|xlsx) OUT.tif [--mode points|sigungu]
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import shapefile
from rasterio.features import rasterize
from shapely.geometry import shape as shp_shape

HEAD_KEYS = ["head", "heads", "n_head", "count", "num", "사육두수", "두수", "마릿수", "총두수", "사육수", "goats",
             "goat", "염소", "total"]
LON_KEYS = ["lon", "longitude", "long", "경도", "x_wgs84"]
LAT_KEYS = ["lat", "latitude", "위도", "y_wgs84"]
X_KEYS = ["x", "x5179", "xcoord", "x_coord", "tm_x", "utmk_x", "point_x"]
Y_KEYS = ["y", "y5179", "ycoord", "y_coord", "tm_y", "utmk_y", "point_y"]
SIG_KEYS = ["sig_cd", "sigcd", "sig_code", "시군구코드", "sgg_cd", "sgg_code", "adm_cd"]
SNAP_PX = 5
DONG_KEYS = ["bjd_cd", "법정동코드", "ldong_cd", "emd_cd", "adm_dr_cd", "법정동", "lawd_cd"]


def norm(s):
    return re.sub(r"[\s_\-()/.]", "", str(s)).lower()


def pick(cols, keys, contains=False):
    nc = {norm(c): c for c in cols}
    for k in keys:
        if norm(k) in nc:
            return nc[norm(k)]
    if contains:
        for k in keys:
            for n, c in nc.items():
                if norm(k) in n:
                    return c
    return None


def read_table(path):
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(p)
    for enc in ("utf-8-sig", "cp949", "euc-kr", "utf-8"):
        try:
            return pd.read_csv(p, encoding=enc, low_memory=False)
        except UnicodeDecodeError:
            continue
    raise SystemExit(f"cannot read {p}")


def sig_raster(sig_shp, profile):
    r = shapefile.Reader(str(sig_shp), encoding="cp949")
    names = [f[0] for f in r.fields[1:]]
    i_cd = names.index("SIG_CD")
    shapes, codes = [], []
    for sr in r.iterShapeRecords():
        code = int(sr.record[i_cd])
        shapes.append((shp_shape(sr.shape.__geo_interface__), code))
        codes.append(code)
    arr = rasterize(shapes, out_shape=(profile["height"], profile["width"]),
                    transform=profile["transform"], fill=0, dtype="int32", all_touched=False)
    return arr, sorted(set(codes))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("glw")
    ap.add_argument("sig")
    ap.add_argument("register")
    ap.add_argument("out")
    ap.add_argument("--mode", choices=["auto", "points", "sigungu"], default="auto")
    ap.add_argument("--head-col")
    ap.add_argument("--filter-col", help="column naming the species, e.g. 축종")
    ap.add_argument("--filter-val", help="value to keep, e.g. 염소")
    a = ap.parse_args()

    with rasterio.open(a.glw) as src:
        glw = src.read(1).astype(np.float64)
        prof = src.profile.copy()
        crs = src.crs
    sig, codes = sig_raster(a.sig, prof)
    rok = sig > 0
    print(f"GLW4: {np.isfinite(glw).sum():,} valid pixels; ROK pixels from sig.shp: {rok.sum():,} "
          f"({len(codes)} sigungu); GLW values already inside the ROK: {np.isfinite(glw[rok]).sum():,} "
          f"pixels, {np.nansum(glw[rok]):,.0f} head")

    df = read_table(a.register)
    print(f"register: {len(df):,} rows; columns: {list(df.columns)}")
    if a.filter_col:
        df = df[df[a.filter_col].astype(str).str.contains(a.filter_val, na=False)].copy()
        print(f"  kept {len(df):,} rows where {a.filter_col} contains {a.filter_val!r}")
    hcol = a.head_col or pick(df.columns, HEAD_KEYS) or pick(df.columns, HEAD_KEYS, contains=True)
    if hcol is None:
        raise SystemExit("no head-count column found; pass --head-col")
    heads = pd.to_numeric(df[hcol], errors="coerce").fillna(0.0).clip(lower=0.0)
    print(f"  head column: {hcol!r}; total {heads.sum():,.0f} head on {int((heads > 0).sum()):,} rows")

    dens = np.zeros(glw.shape, np.float64)
    lon, lat = pick(df.columns, LON_KEYS), pick(df.columns, LAT_KEYS)
    x, y = pick(df.columns, X_KEYS), pick(df.columns, Y_KEYS)
    sigc = pick(df.columns, SIG_KEYS)
    dong = pick(df.columns, DONG_KEYS)
    mode = a.mode
    if mode == "auto":
        mode = "points" if ((lon and lat) or (x and y)) else "sigungu"
    if mode == "points":
        if lon and lat:
            from pyproj import Transformer
            tr = Transformer.from_crs("EPSG:4326", crs.to_string(), always_xy=True)
            xs, ys = tr.transform(pd.to_numeric(df[lon], errors="coerce").values,
                                  pd.to_numeric(df[lat], errors="coerce").values)
        else:
            xs = pd.to_numeric(df[x], errors="coerce").values
            ys = pd.to_numeric(df[y], errors="coerce").values
        T = prof["transform"]
        cols = np.floor((xs - T.c) / T.a).astype("float")
        rows = np.floor((ys - T.f) / T.e).astype("float")
        ok = np.isfinite(cols) & np.isfinite(rows) & (cols >= 0) & (rows >= 0) & \
            (cols < prof["width"]) & (rows < prof["height"])
        r_, c_ = rows[ok].astype(int), cols[ok].astype(int)
        inside = rok[r_, c_]
        hv = heads.values[ok]
        np.add.at(dens, (r_[inside], c_[inside]), hv[inside])
        # farms whose 1-km pixel centre falls outside the sigungu polygons (coast, islands):
        # move them to the nearest ROK land pixel within SNAP_PX pixels
        from scipy.ndimage import distance_transform_edt
        dist, (ir, ic) = distance_transform_edt(~rok, return_indices=True)
        out_r, out_c, out_h = r_[~inside], c_[~inside], hv[~inside]
        near = dist[out_r, out_c] <= SNAP_PX
        np.add.at(dens, (ir[out_r[near], out_c[near]], ic[out_r[near], out_c[near]]), out_h[near])
        lost = heads.values[~ok].sum() + out_h[~near].sum()
        print(f"  points: {int(ok.sum()):,} on the grid, {int(inside.sum()):,} inside ROK land, "
              f"{int(near.sum()):,} moved to the nearest ROK land pixel (≤ {SNAP_PX} km; {out_h[near].sum():,.0f} head); "
              f"{int((~near).sum()) + int((~ok).sum()):,} farms ({lost:,.0f} head) dropped")
        # heads per 1-km pixel = head/km2 (pixel area 1 km2)
        dens /= (abs(T.a) * abs(T.e)) / 1e6
    else:
        if sigc:
            cd = df[sigc].astype(str).str.extract(r"(\d{5})")[0]
        elif dong:
            cd = df[dong].astype(str).str.extract(r"(\d{5})")[0]
        else:
            raise SystemExit("no coordinates and no sigungu/legal-dong code column")
        tot = heads.groupby(cd.astype(float)).sum()
        px = pd.Series(np.bincount(sig[rok].ravel()), dtype=float)
        unmatched = [int(k) for k in tot.index if k not in codes]
        if unmatched:
            print(f"  {len(unmatched)} register code(s) not in sig.shp (head "
                  f"{tot.loc[unmatched].sum():,.0f}); e.g. {unmatched[:8]}")
        area_km2 = abs(prof["transform"].a * prof["transform"].e) / 1e6
        lut = np.zeros(int(sig.max()) + 1, np.float64)
        for k, h in tot.items():
            k = int(k)
            if k < lut.size and k in codes and px.get(k, 0) > 0:
                lut[k] = h / (px[k] * area_km2)
        dens = np.where(rok, lut[sig], 0.0)
    out = glw.copy()
    out[rok] = dens[rok]
    print(f"ROK after fill: {np.nansum(out[rok]):,.0f} head over {rok.sum():,} km2 "
          f"(mean {np.nanmean(out[rok]):.2f} head/km2, max {np.nanmax(out[rok]):.1f}); "
          f"outside the ROK unchanged: {np.allclose(np.nan_to_num(out[~rok]), np.nan_to_num(glw[~rok]))}")
    prof.update(dtype="float32", nodata=np.nan, compress="lzw")
    with rasterio.open(a.out, "w", **prof) as dst:
        dst.write(out.astype(np.float32), 1)
        dst.update_tags(note="GLW4 goats (head/km2) with the ROK filled from the national goat register",
                        rok_head_total=f"{np.nansum(out[rok]):.0f}", register=Path(a.register).name)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    sys.exit(main())
