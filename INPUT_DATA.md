# Input layers for the full runs

Put these in one folder (`--data-root`); the model finds them by name. All rasters are read and resampled to the
grid of the PPR suitability map (EPSG:5179, Korea 2000 / Unified CS).

| File (name pattern) | Layer | Source (Supplementary Table S4) |
|---|---|---|
| `PPR_MaxEnt_Final_Map_*.tif` | PPR environmental suitability Hp; defines the grid | MaxEnt model fitted by the authors to 266 EMPRES-i outbreak locations |
| `gorani_habitat_suitability_cloglog_donut_500m_10km.tif` | Water-deer habitat suitability Hg | MaxEnt model of water-deer occurrence (authors) |
| `MDL.shp` (+ .dbf, .shx, .prj) | Military Demarcation Line | Digitised line, rasterised on the model grid |
| `sig.shp` (+ .dbf, .shx, .prj, .cpg) | Sigungu boundaries of the ROK (250 polygons) | Administrative boundaries of the ROK |
| point shapefile whose name contains `occ` (folder `PPR_occ/`) | PPR outbreak locations used as seeds | FAO EMPRES-i |
| `PPR_env/GLW4GTS_1km_5179.tif` | Goat density, China and the DPRK | Gridded Livestock of the World, version 4 |
| `PPR_env/GLW4GTS_1km_5179_ROKfilled.tif` | Same, completed in the ROK from the national register (`code/rok_goat_fill.py`) | Available on request (personal data in the register) |
| `PPR_env/wc2.1_30s_elev_studyarea_1km_5179.tif` | Elevation | WorldClim 2.1 |
| `PPR_env/wc2.1_30s_bio_12_studyarea_1km_5179.tif` | Annual precipitation | WorldClim 2.1 |

The corridor coordinates and weights (Table S6), the monthly dispersal radius from 44,618 road-kill records and all
parameter values are written in the code. The hold-out run (`--preset holdout2015`) reads `data/PPR_seed_points_dated.csv`.
`tests/world/` holds a synthetic set of the same files and shows the expected layout.
