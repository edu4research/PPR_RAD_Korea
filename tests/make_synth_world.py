"""Synthetic China–Korea world (EPSG:5179) that exercises every path of the model:
MaxEnt base raster, water-deer raster, MDL line with the six real gate positions,
ROK sigungu polygons (UTF-8 DBF), an ROK island next to the DPRK coast, an
unfenced estuary contact west of the fence, outbreak points in China and three
environmental rasters. Run with python3.12 (GDAL bindings)."""
import sys
from pathlib import Path
import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("world")
NATIVE_KM = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
(OUT / "PPR_env").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(7)

P = osr.SpatialReference()
# Korea 2000 / Unified CS written without EPSG axis order (as ArcGIS/ESRI WKT1 does)
P.ImportFromProj4("+proj=tmerc +lat_0=38 +lon_0=127.5 +k=0.9996 +x_0=1000000 +y_0=2000000 "
                  "+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
G = osr.SpatialReference(); G.ImportFromEPSG(4326)
for s in (P, G):
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
T = osr.CoordinateTransformation(G, P)


def poly(coords):
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for lon, lat in coords + [coords[0]]:
        ring.AddPoint_2D(float(lon), float(lat))
    g = ogr.Geometry(ogr.wkbPolygon)
    g.AddGeometry(ring)
    g.Transform(T)
    return g


def line(coords):
    g = ogr.Geometry(ogr.wkbLineString)
    for lon, lat in coords:
        g.AddPoint_2D(float(lon), float(lat))
    g.Transform(T)
    return g


CHINA = poly([(113, 33), (120.5, 33), (120.5, 36.5), (122.2, 37.3), (121.0, 38.5), (121.5, 39.5),
              (123.2, 39.6), (124.30, 40.0), (125.45, 40.72), (126.95, 41.62), (128.45, 42.12),
              (129.95, 42.62), (130.65, 42.52), (131.3, 43.0), (131.3, 46.5), (113, 46.5)])
PEN = poly([(124.4, 39.9), (125.5, 40.6), (127.0, 41.5), (128.5, 42.0), (130.0, 42.5), (130.6, 42.35),
            (129.7, 41.0), (128.5, 39.5), (128.4, 38.6), (129.4, 36.0), (129.3, 35.2), (127.5, 34.6),
            (126.3, 34.6), (126.4, 35.5), (126.5, 36.8), (126.55, 37.3), (126.25, 37.55),
            (126.15, 37.80), (125.5, 37.85), (125.1, 38.0), (124.7, 38.05), (124.8, 38.3),
            (125.3, 38.7), (125.0, 39.3), (124.6, 39.6)])
ISLAND = poly([(124.93, 37.91), (125.05, 37.91), (125.05, 37.975), (124.93, 37.975)])
MDL_LL = [(126.55, 37.80), (126.73, 37.95), (127.08, 38.02), (127.31, 38.24), (127.68, 38.25),
          (128.00, 38.28), (128.25, 38.45), (128.45, 38.62), (128.65, 38.66)]
ROK_REGION = poly([(124.0, 37.80)] + MDL_LL + [(131, 38.66), (131, 33), (124, 33)])
ROK = PEN.Intersection(ROK_REGION).Union(ISLAND)
LAND = CHINA.Union(PEN).Union(ISLAND)

# grid of the MaxEnt raster (native resolution NATIVE_KM)
x0, x1, y0, y1 = -150_000.0, 1_500_000.0, 1_350_000.0, 2_850_000.0
px = NATIVE_KM * 1000.0
nc, nr = int(round((x1 - x0) / px)), int(round((y1 - y0) / px))
gt = (x0, px, 0.0, y1, 0.0, -px)
wkt = P.ExportToWkt()


def burn(geom, nr_, nc_, gt_):
    drv = ogr.GetDriverByName("Memory")
    ds = drv.CreateDataSource("m")
    lyr = ds.CreateLayer("l", srs=P, geom_type=ogr.wkbMultiPolygon)
    f = ogr.Feature(lyr.GetLayerDefn()); f.SetGeometry(geom); lyr.CreateFeature(f)
    r = gdal.GetDriverByName("MEM").Create("", nc_, nr_, 1, gdal.GDT_Byte)
    r.SetGeoTransform(gt_); r.SetProjection(wkt)
    gdal.RasterizeLayer(r, [1], lyr, burn_values=[1])
    return r.GetRasterBand(1).ReadAsArray().astype(bool)


def smooth_field(nr_, nc_, scale):
    from scipy.ndimage import gaussian_filter
    a = gaussian_filter(rng.random((nr_, nc_)), scale)
    return (a - a.min()) / (a.max() - a.min() + 1e-12)


def write_tif(path, arr, gt_, nodata=-9999.0):
    d = gdal.GetDriverByName("GTiff").Create(str(path), arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32)
    d.SetGeoTransform(gt_); d.SetProjection(wkt)
    b = d.GetRasterBand(1); b.WriteArray(np.where(np.isfinite(arr), arr, nodata).astype(np.float32))
    b.SetNoDataValue(nodata); d = None


land = burn(LAND, nr, nc, gt)
mx = smooth_field(nr, nc, 6)
write_tif(OUT / "PPR_MaxEnt_Final_Map_studyarea.tif", np.where(land, 0.15 + 0.8 * mx, np.nan), gt)
# water-deer raster at half the pixel size, peninsula + island only
gt2 = (x0, px / 2, 0.0, y1, 0.0, -px / 2)
pen2 = burn(PEN.Union(ISLAND), nr * 2, nc * 2, gt2)
gv = smooth_field(nr * 2, nc * 2, 10)
write_tif(OUT / "gorani_habitat_suitability_cloglog_donut_500m_10km.tif",
          np.where(pen2, 0.1 + 0.85 * gv, np.nan), gt2)
# environmental rasters (land only)
yy = np.linspace(1, 0, nr)[:, None] * np.ones((1, nc))
elev = 50 + 900 * smooth_field(nr, nc, 8) * (0.5 + yy)
write_tif(OUT / "PPR_env" / "wc2.1_30s_elev_studyarea_1km_5179.tif", np.where(land, elev, np.nan), gt)
write_tif(OUT / "PPR_env" / "GLW4GTS_1km_5179.tif", np.where(land, 5 + 200 * smooth_field(nr, nc, 3), np.nan), gt)
write_tif(OUT / "PPR_env" / "wc2.1_30s_bio_12_studyarea_1km_5179.tif",
          np.where(land, 500 + 900 * smooth_field(nr, nc, 10), np.nan), gt)


def write_shp(path, geoms, fields, srs, gtype):
    drv = ogr.GetDriverByName("ESRI Shapefile")
    if Path(path).exists():
        drv.DeleteDataSource(str(path))
    ds = drv.CreateDataSource(str(path))
    lyr = ds.CreateLayer(Path(path).stem, srs=srs, geom_type=gtype, options=["ENCODING=UTF-8"])
    for name, typ in fields:
        lyr.CreateField(ogr.FieldDefn(name, typ))
    for g, attrs in geoms:
        f = ogr.Feature(lyr.GetLayerDefn())
        for k, v in attrs.items():
            f.SetField(k, v)
        f.SetGeometry(g)
        lyr.CreateFeature(f)
    ds = None


write_shp(OUT / "MDL.shp", [(line(MDL_LL), {"NAME": "MDL"})], [("NAME", ogr.OFTString)], P, ogr.wkbLineString)
# ROK sigungu = 0.3° boxes clipped to the ROK, codes with real province prefixes
counties, k = [], 0
for lon in np.arange(124.0, 131.0, 0.3):
    for lat in np.arange(33.0, 38.7, 0.3):
        box = poly([(lon, lat), (lon + 0.3, lat), (lon + 0.3, lat + 0.3), (lon, lat + 0.3)])
        g = ROK.Intersection(box)
        if g is None or g.IsEmpty() or g.GetArea() < 1e5:
            continue
        k += 1
        if g.Intersects(ISLAND):
            code, en, kr = "28720", "Ongjin-gun", "옹진군"
        else:
            pre = ("41" if (lat >= 37.0 and lon < 127.6) else "51" if lat >= 37.0 else
                   "44" if (lat >= 36.0 and lon < 127.5) else "43" if lat >= 36.0 else
                   "46" if lon < 127.5 else "48")
            code, en, kr = f"{pre}{k:03d}", f"County-{k:03d}", f"군{k:03d}"
        counties.append((g, {"SIG_CD": code, "SIG_ENG_NM": en, "SIG_KOR_NM": kr}))
write_shp(OUT / "sig.shp", counties, [("SIG_CD", ogr.OFTString), ("SIG_ENG_NM", ogr.OFTString),
                                       ("SIG_KOR_NM", ogr.OFTString)], P, ogr.wkbPolygon)
# outbreak points in China (WGS84)
pts = []
for lon, lat in list(zip(rng.uniform(123.2, 125.5, 8), rng.uniform(41.2, 42.6, 8))) + \
        list(zip(rng.uniform(117.5, 120.0, 10), rng.uniform(35, 39.5, 10))):
    p = ogr.Geometry(ogr.wkbPoint); p.AddPoint_2D(float(lon), float(lat))
    pts.append((p, {"ID": len(pts) + 1}))
write_shp(OUT / "ppr_occ_points.shp", pts, [("ID", ogr.OFTInteger)], G, ogr.wkbPoint)
print(f"world: {nc}x{nr} px at {NATIVE_KM} km, {len(counties)} sigungu, ROK area "
      f"{ROK.GetArea() / 1e6:,.0f} km², island {ISLAND.GetArea() / 1e6:.0f} km² -> {OUT}")
