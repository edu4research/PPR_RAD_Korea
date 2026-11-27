



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PPR MODEL v9.10 — DYNAMIC REACTION–ADVECTION–DIFFUSION EPIDEMIC (persistent S/I/R)
  v9.10 = v9.9 + the re-analysis for the Microorganisms review (search "v9.10"):
    1. ROK/DPRK split from the ROK sigungu sub-cells (+ land share at the border);
       ROK exclaves on the grid ('report' | 'mask'); estuary policy ('open' | 'sealed');
       establishment judged on ROK area; buffer in months; extra report months
    2. SMCDA check: distance-only benchmarks, ENTRY from gates vs gates + estuary,
       entropy / equal / rank-sum / leave-one-out weights; Table 2 CSV (shares add to 100)
    3. fence wait per gate: DPRK side vs ROK side of the fence (H1)
    4. Table 1 (timeline + state of the ROK), Table 3 (V1/V2 on ROK territory),
       first sigungu per province, Figure 9 (maps A–C + time series D, no title)
    5. global sensitivity analysis on the peninsular leg: Morris -> LHS + PRCC (seeded)
    6. Table 4 scenarios (Yalu entry, uniform fence, sealed estuary, γ×2/×4, no
       seasonality, A* 0.01/0.10, corridor variants, exclaves) + 4-km grid preset
    7. GeoTIFFs of scores / arrival months; temporal hold-out run (--preset holdout2015:
       seeds reported before 2015 only) for ppr_external_test.py
    8. round 2 (Reviewer 2, F1): what each ranking adds beyond distance from the entry routes
       (partial Spearman rho, gain in R², best distance + ranking mixture; spatial block bootstrap)
       -> Table S7b, smcda_beyond_distance()
    9. round 2 (F1): two-step surveillance rankings (distance bands first, risk score or goat density within a
       band) in the targeting comparison (Table S9); --batch targeting runs that part alone
   10. round 2 (F1): spatial block-bootstrap intervals for the rho of each surveillance ranking and for its
       difference from distance alone (R2_targeting_rho_ci.csv), and a per-cell table of the ROK cells used
       (R2_targeting_cells.csv)
   11. point locations (seeds, corridor gates before snapping, Yalu entry, centroids of very small sigungu) are placed
       in the cell that contains them, as for every raster (default since 26 Nov 2026; the paper reports these runs);
       --point-cell round (nearest cell corner, up to half a cell off) reproduces the runs made before that date
  Command line: see the end of this file (python PPR_model_v9_10.py --help).
  v9.9 = v9.8 + review fixes (search "v9.9 FIX"):
    A. seasonal drivers of the Korea leg start at the NK-introduction month (season_offset)
    B. a leg never stops before its stated minimum horizon
    C. DPRK/ROK split from the national (sigungu) boundary, not one grid row
    D. ENTRY for South-Korean domains = corridor gates + unfenced border (western terminus)
    E. Korea-leg vaccination scenarios V1 (combined host, upper bound) and V2 (ROK livestock,
       triggered by the first ROK detection, coverage capped by the livestock share)
              + CORRIDOR-GATED MDL FENCE (water-deer crossing points)
              + PLAIN-WHITE SOUTH KOREA PANELS
              + SOUTH-KOREA COUNTY (SIGUNGU) CSV EXPORT
              + SPATIAL MCDA RISK-AREA THRESHOLDS, PER DOMAIN
              + SOUTH-KOREA-ONLY RISK LEVELS + SK FRAME CROP  ← NEW IN v9.8
=============================================================================
WHAT IS NEW IN v9.8 (everything else is v9.7, unchanged)
-----------------------------------------------------------------------------
SOUTH KOREAN RISK LEVELS ARE NOW READ OFF A SOUTH KOREAN SCALE.

The problem this fixes was visible on the v9.7 "korea_reached" panel: North
Korea came out red, and essentially ALL of South Korea came out blue. That is
not a statement about South Korea. It is what a shared threshold scale does.
NK and SK were standardised and thresholded together, NK scores higher on the
criteria carrying most of the weight, so the whole South Korean distribution
sat below the upper Jenks cuts and collapsed into the two lowest classes. Every
South Korean county then reads "Very low / Low", and a map that cannot separate
one South Korean county from another cannot allocate surveillance effort
between them — which is the only thing a surveillance design needs from it.
That domain also FAILED its arrival-order validation in that run (arrival order
not monotonic, rho = -0.128), for the related reason that a single ENTRY origin
cannot describe both legs of the incursion.

TWO CHANGES, BOTH NARROW:

  1. THE DEFAULT DOMAIN LIST IS NOW SOUTH KOREA ONLY.
         SMCDA_DOMAINS = ["sk_reached", "sk"]
     Both are kept because they answer different questions:
         "sk_reached"  CONDITIONAL  — within the area the modelled incursion
                       actually covers, where does surveillance go first? This
                       is the primary domain and the one the county CSV, the
                       mobility field and the JSON summary are built from.
         "sk"          PROSPECTIVE  — across all of South Korea, whether or not
                       the model takes PPR there: if PPR enters through the
                       gates, where is risk highest? This is the domain for
                       pre-incursion planning.
     Nothing was deleted. "korea", "korea_reached", "nk" and "all" are still
     fully supported and are the correct domains for a NK-vs-SK or a
     China-vs-Korea comparison; put them back in SMCDA_DOMAINS to get them.
     What has changed is that they are no longer in the DEFAULT run, so a South
     Korean risk level cannot be read off a peninsula-scaled map by accident.

  2. THE KOREA PANEL CROPS TO SOUTH KOREA (Config.SMCDA_KOREA_PANEL_CROP).
     An SK-only classification drawn on the whole-peninsula frame spent about
     half the figure on grey North Korean land and shrank South Korea to where
     county-level differences could not be read. "auto" (the default) crops the
     frame to South Korea whenever the domain being drawn contains no North
     Korean cells, and keeps the full peninsula otherwise. The MDL and the
     corridor gates are deliberately pulled into frame along the top edge: they
     are the introduction route the ENTRY criterion measures from.

     THIS IS A CHANGE OF VIEW ONLY. No cell is reclassified, no threshold is
     recomputed, no statistic changes. The classification was already SK-only
     in v9.7 for the "sk" and "sk_reached" domains — v9.7 drew it on too wide a
     frame. Set SMCDA_KOREA_PANEL_CROP="peninsula" for the v9.7 figure.

WHAT THIS DOES NOT FIX — READ BEFORE THE CLASSES FEED THE SURVEILLANCE DESIGN
  * The risk-based surveillance design consumes a BINARY class (high risk /
    low risk) and strata built from it. The SMCDA produces FIVE ordered levels.
    No rule for collapsing five into two is implemented here, because the choice
    of cut-point is a design decision with a direct effect on how many counties
    and farms enter the high-risk stratum, and inventing one silently would put
    an unstated assumption underneath the sample size. Decide the cut, state it,
    then apply it.
  * The AHP weights are still CONSTRUCTED from an ordinal priority vector, not
    elicited from experts. Every run says so, in the log, in the CSV header and
    in the JSON.
  * The Jenks breaks are still descriptions of one map's histogram. They are
    reported per domain and are NOT comparable between domains — including
    between "sk" and "sk_reached", which are cut through different
    distributions of the same score.

WHAT IS NEW IN v9.7 (everything else is v9.3, unchanged)
-----------------------------------------------------------------------------
THE MODEL COULD SAY WHERE PPR ARRIVES AND WHEN, BUT NOT WHAT LEVEL OF RISK AN
AREA CARRIES. Those are different questions. The first is the epidemic. The
second is a DECISION: it combines several things that matter — host density,
distance from the entry route, transmission potential, habitat suitability,
wildlife density, landscape connectivity — measured in different units and not
mattering equally. A quantile cut on one layer cannot answer it, and the
habitat/accessible/barrier map in v9.3 was exactly that.

v9.6 adds a SPATIAL MULTI-CRITERIA DECISION ANALYSIS (SMCDA) layer that produces
EXPLICIT, DEFENSIBLE THRESHOLDS between named risk levels, and then tests them.

THE CHAIN (each step is separable and separately auditable)
  1. CRITERIA + STANDARDISATION. Each criterion layer is mapped onto [0,1] where
     1 = maximum contribution to risk, by min-max, percentile rank, or a fuzzy
     membership between two stated control points. Cost criteria (distance from
     the entry route) are inverted after scaling. Objective, given the layers.
  2. WEIGHTS — AHP. Pairwise judgements -> principal eigenvector -> weights, with
     Saaty's CONSISTENCY RATIO reported every run. CR > 0.10 means the judgements
     contradict each other and the weights are not defensible; the code says so
     and can be set to abort. Rank-sum, Shannon-entropy and equal weighting are
     available as alternatives, and the entropy weights are ALWAYS printed beside
     the AHP weights as an independent cross-check.
  3. AGGREGATION — WLC or OWA. WLC is fully compensatory: a very high score on
     one criterion can carry a cell. OWA takes an ANDness parameter, so you can
     require a cell to score well on most criteria before it is called high risk.
     andness = 0.5 reproduces WLC exactly. Non-compensatory VETOES are supported
     for constraints that genuinely cannot be traded away.
  4. THRESHOLDS — Fisher-Jenks natural breaks, quantile, equal-interval, standard
     deviation, or YOUR OWN absolute cut values. This is the deliverable.

ONE CLASSIFICATION PER DOMAIN — AND WHY THAT IS NOT A CONVENIENCE
  A risk CLASS is relative to the area it was computed over. "High risk" on a
  China-wide grid and "High risk" in South Korea are different statements,
  because the score is standardised and thresholded within whatever domain it is
  given. So v9.7 computes, maps, validates and exports a SEPARATE classification
  for each domain in Config.SMCDA_DOMAINS: "sk" (all of South Korea), "korea"
  (the peninsula, NK+SK on one scale), "nk", "all" — and the two domains added
  in v9.7:

  THE INVADED FOOTPRINT — "sk_reached" and "korea_reached"
  Restricted to the area PPR ACTUALLY REACHED by Config.SMCDA_REACHED_BY_MONTH
  (default: the last month of the Korea leg; set 209 for exactly 18 years). This
  is the coloured part of the progression figure, not the whole country, and it
  answers the operational question: "within the area the incursion actually
  covers, where does surveillance effort go first?"

  ‼ WHAT A _reached DOMAIN COSTS. Every other domain is defined WITHOUT
  reference to the epidemic, which is what makes validating it against that
  epidemic a real test. A _reached domain is defined BY the epidemic outcome.
  The score VALUES are still computed from static layers alone and still never
  see the simulation, but the DOMAIN is conditioned on it, and the
  standardisation ranges and Jenks breaks are computed inside it. So:
    * AUC against "reached by end of run" becomes UNDEFINED — every cell in the
      domain is reached by construction. The code says so instead of printing a
      number. Arrival-order monotonicity and Spearman rho remain valid, and are
      now the HARDER test: the score can no longer earn credit for the easy part
      of the job (separating reached from unreached).
    * Its class thresholds are NOT comparable with the unrestricted domain's.
    * It is a CONDITIONAL risk map — "given the incursion got this far" — not a
      prospective one. For pre-incursion planning, read "sk".
  Both maps are produced, so the two questions can be answered side by side.

  Each domain also has its OWN INTRODUCTION ROUTE, and the ENTRY criterion must
  measure travel cost from the route that actually applies:
      South Korea       entered through the water-deer corridor gates;
      North Korea       entered from China across the northern border;
      Korean peninsula  entered from China — the gates are an INTERNAL barrier
                        on the way south, not the entry to the peninsula;
      continental grid  entered at the outbreak seed points themselves.
  Config.SMCDA_ENTRY_SOURCE="auto" selects this per domain and logs its choice.
  Pointing ENTRY at the wrong origin makes the classification rank the epidemic
  BACKWARDS: the cells nearest the assumed entry are the ones reached last. The
  validation pass detects exactly that, and the per-criterion diagnostic names
  the offending criterion instead of leaving you to guess.

WHAT MAKES THE THRESHOLDS MEAN SOMETHING
  Jenks breaks describe one map's histogram. They will manufacture a "very high
  risk" class out of a uniformly low-risk surface, and they move whenever the
  layers or the domain change. So every domain is VALIDATED against the
  mechanistic RAD epidemic — and against the RIGHT ONE: the Korean domains are
  tested against the Korea leg (China border -> North Korea -> corridor gates ->
  South Korea), the continental domain against the continental run. Asking
  whether a South-Korea classification predicts arrival order in a run where
  South Korea is reached last, all at once, answers nothing.
    * is mean ARRIVAL MONTH monotonic across the classes (higher class reached
      earlier)?
    * AUC of the continuous score against "reached by the end of the run",
      which is independent of where the breaks fall;
    * Spearman rho between score and arrival month (must be NEGATIVE);
    * per-class reached fraction and mean final attack rate;
    * WHEN MONOTONICITY FAILS: the Spearman rho of EACH criterion on its own
      against arrival month, so the criterion pointing the wrong way over that
      domain is named, with its weight beside it.
  The score is built from STATIC layers only and never sees the epidemic, so
  this is a real test rather than a restatement of a fit. A one-at-a-time weight
  sensitivity analysis (+20% per criterion) reports what percentage of cells
  change class, i.e. how much the answer rests on the weighting judgement.

NEW OUTPUTS (all in Config.OUT_DIR, and copied to Config.DL_DIR):
  1. The RISK-AREA THRESHOLD TABLE per domain, printed to the log: class, score
     range, cell count, % of domain, fraction reached by the epidemic, mean
     arrival month and area in km², followed by the validation verdicts and, on
     failure, the per-criterion diagnostic.
  2. PPR_v9.7_smcda_thresholds.json — every domain, machine-readable: weights,
     criteria used, criteria DROPPED for want of a layer, the consistency ratio,
     the ENTRY origin, which epidemic each domain was validated against, the
     full sensitivity table and the diagnostic.
  3. PPR_v9.7_smcda_risk_areas_<domain>.png — the classified map on the model
     grid beside the score histogram with the class boundaries drawn on it, so
     the thresholds are visible as cuts in the real distribution.
  4. PPR_v9.7_smcda_korea_panel_<domain>.png — for the Korean domains, the risk
     classes on the SAME white-canvas peninsula frame as the progression figure
     (county lines, waterways, MDL, corridor gates, China introduction seeds),
     beside the epidemic's own arrival-month map on the identical frame. The
     left map never saw the right one; putting them side by side at one scale is
     the whole point.
  5. PPR_v9.7_smcda_county_risk_levels_<domain>.csv — ONE ROW PER SOUTH-KOREAN
     COUNTY: area-weighted mean score, the share of the county in each risk
     class, the assigned level (dominant area) and the highest level present
     anywhere in it. The header carries the validation verdict, so a table whose
     classes failed the arrival-order test says so on its own face.

WHAT THIS DOES NOT CHANGE
  Nothing about the epidemic. The SMCDA layer reads the model's static surfaces
  and writes its own outputs. USE_SMCDA_CLASS_IN_TRANSPORT exists for the case
  where you WANT the classes to drive movement, and it is OFF by default — with
  it on, the simulated spread changes and the validation stops being
  independent, both of which the code warns about at run time.

WHERE JUDGEMENT ENTERS, STATED PLAINLY
  The default AHP matrix is CONSTRUCTED from an ordinal priority vector
  (Config.SMCDA_PRIORITY_VECTOR), not elicited from an expert panel. Because it
  comes from a single vector it is near-consistent by construction, so its low CR
  confirms the arithmetic and says nothing about whether the priorities are
  epidemiologically correct. The run log says this every time. Replace
  SMCDA_AHP_MATRIX with your panel's pairwise judgements before these classes
  allocate surveillance effort.
-----------------------------------------------------------------------------
WHAT IS NEW IN v9.3 (everything else is v9.2, unchanged)
-----------------------------------------------------------------------------
THE MODEL PRODUCED MAPS BUT NO TABLE. A raster panel cannot be joined to a
surveillance database, a movement-control order or a vaccination budget — those
are all written per COUNTY. v9.3 adds an administrative zonal layer that turns
the modelled cumulative attack-rate fields into a per-sigungu table for SOUTH
KOREA, written the moment PPR has crossed the MDL and the Korea leg finishes.

NEW OUTPUTS (all in Config.OUT_DIR, and copied to Config.DL_DIR):
  1. PPR_v9.3_south_korea_counties_korea_leg_summary.csv
        ONE ROW PER SOUTH-KOREAN COUNTY. Columns include the county code
        (SIG_CD), Korean and English names, province, the month PPR arrived
        (both as a simulation month and as a nominal calendar month), the
        final mean / max cumulative attack rate, the fraction of the county's
        area past the arrival threshold, the peak infected density and the
        month it peaked, and the NEAREST WATER-DEER CORRIDOR GATE with its
        distance — i.e. which DMZ crossing most plausibly delivered PPR there.
  2. PPR_v9.3_south_korea_counties_korea_leg_monthly_wide.csv
        counties × months matrix of mean cumulative attack (%), for plotting.
  3. PPR_v9.3_south_korea_counties_korea_leg_long.csv
        tidy/long form (one row per county-month) for pandas / R / Excel
        pivot tables.
  4. The same three files for the CONTINENTAL run
        (…_main_run_…), when SK_COUNTY_ALSO_MAIN_RUN is on.
  5. PPR_v9.3_SK_county_arrival.png — choropleth of arrival month and final
        attack rate BY COUNTY, so the CSV is visually auditable.

HOW THE ZONAL STATISTICS ARE COMPUTED (this is the part that matters)
  * sig.shp (or any admin-2 layer: GADM KOR_2, sigungu, adm2 …) is reprojected
    onto the model CRS and rasterized at SK_COUNTY_OVERSAMPLE× the model
    resolution. Each county therefore gets a WEIGHT per model cell equal to the
    number of sub-cells it covers there, so a county that covers a third of a
    cell contributes a third of that cell — instead of the all-or-nothing
    assignment you get from rasterizing straight onto a multi-km grid.
  * Counties too small to capture even one sub-cell (the Seoul gu districts are
    ~10-25 km² against a 10-15 km cell) fall back to POINT-ON-SURFACE sampling
    and are flagged `sampling=centroid` in the CSV, so you can see exactly which
    rows are grid-limited rather than silently trusting them.
  * Korean shapefiles are almost always CP949/EUC-KR, not UTF-8. The encoding is
    AUTO-DETECTED by scoring how much Hangul each candidate decodes, and the CSV
    is written UTF-8 with a BOM so Excel opens the Korean names correctly.
  * Two arrival definitions are reported side by side, because they answer
    different questions: `arrival_month` (the county MEAN attack crosses the
    threshold — "the county is affected") and `arrival_month_first_cell` (ANY
    cell in the county crosses it — "PPR is present somewhere in the county").
  * GRID-RESOLUTION CAVEAT, logged loudly at run time: with BASE_RASTER="maxent"
    the cell is ~10-15 km, which is coarser than many South-Korean counties. For
    county-level conclusions run BASE_RASTER="gorani" (~5 km at FACTOR=10) or
    lower FACTOR; the CSV reports model_cells and sampling per row so grid-
    limited counties are identifiable.
-----------------------------------------------------------------------------
WHAT WAS NEW IN v9.2 (retained)
-----------------------------------------------------------------------------
THE SOUTH-KOREA PANELS WERE BROWN EVEN AT 0% ATTACK. In v9.1 the Korea-crop
renderer (draw_korea) painted the habitat-suitability underlay across the WHOLE
peninsula:
    base = np.where(KP_k, SUIT_k, np.nan)      # KP_k = North AND South Korea
    ax.imshow(..., cmap=_mx_cmap, alpha=0.80, ...)
SUIT_k is GOR_EFF (habitat suitability) drawn with the MaxEnt ramp
(#2c7bb6 -> #abd9e9 -> #ffffbf -> #fdae61 -> #d7191c). South Korea sits around
0.4-0.6 suitability, which lands on the pale-orange/tan part of that ramp — so
the south read as "infected" when its cumulative attack rate was exactly zero.
Two secondary contributors: sig.shp county lines added grey texture, and the
single BAND_MIN_PCT_VISIBLE = 10% floor was applied to South Korea too, hiding
genuine low-% arrivals while the tan underlay stayed visible.
v9.2 fixes it:
  1. The suitability underlay is confined to NORTH KOREA (the source region).
  2. SOUTH KOREA is a BLANK WHITE canvas carrying ONLY water bodies / waterways
     (water.shp, clipped to south of the MDL) and the PPR % attack bands. A thin
     coastline keeps the white landmass legible.
  3. NK and SK get SEPARATE visibility floors (BAND_MIN_PCT_VISIBLE for NK,
     SK_BAND_MIN_PCT_VISIBLE = 1% for SK), applied inside draw_korea instead of
     being pre-filtered by the callers.
  4. The lowest SK band colour is darkened (#edf8e9 -> #c7e9c0) so a 1-20%
     arrival is actually visible against white.
  5. County outlines in the Korea panels are OFF by default
     (SK_PANEL_SHOW_COUNTIES).
All of it is switchable from Config (SK_PANEL_*, NK_PANEL_SHOW_SUITABILITY).
Both PPR_v9.3_transboundary_china_nk_sk.png and every frame of
PPR_v9.3_NK_to_SK_progression.png route through draw_korea, so all Korea panels
are fixed by the one change.
-----------------------------------------------------------------------------
WHAT WAS NEW IN v9.1 (retained)
-----------------------------------------------------------------------------
1. THE MDL IS NO LONGER A UNIFORMLY LEAKY WALL.
   v9.0 applied Han et al.'s fence pass-through phi = 0.531 to EVERY cell of the
   digitized MDL, so gorani crossed the whole 250-km line equally. v9.1 replaces
   that scalar with a PER-CELL PERMEABILITY FIELD built from the ecologically
   justified corridor table:
       pass(cell) = sealed + (phi_max * weight * slowdown - sealed) * taper(d)
       taper(d)   = exp(-0.5 * (d / (gate_radius_km/2))^2)
   The fence is SEALED (MDL_SEALED_FACTOR) everywhere except at the named
   crossings, where it opens in proportion to the corridor's weight:
       Cheorwon         1.00   Iron Triangle — flat lowland agriculture + wetlands
       Paju_Imjin       0.85   Imjin reed marsh / Han estuary; floodgate gaps
       Yeoncheon_Hantan 0.80   Continuous Hantan River riparian path
       Hwacheon_Bukhan  0.70   Bukhan River water route through steep terrain
       Yanggu           0.70   Punchbowl basin; high rim limits lateral spread
       Goseong_Coast    0.60   Narrow lowland coastal strip
       Inje_Seorak      0.40   Steep rocky Taebaek (NO COORDINATE SUPPLIED — skipped)
   Coordinates are WGS84 lat/lon, reprojected to the model CRS and SNAPPED onto
   the nearest digitized MDL cell, so a point a few km off MDL.shp still lands
   on the fence. Every gate is logged with its snap distance.
2. SEASONAL GATING. Gate permeability is scaled by the monthly gorani dispersal
   radius (CORRIDOR_SEASONAL_GATING), so crossings concentrate in the
   late-spring/early-summer dispersal peak instead of running all winter.
3. THE NORTH -> SOUTH CROSSING NOW TAKES LONGER, by design. Four levers, all in
   Config, all documented:
       (a) the gated fence itself — ~6 narrow gates instead of a 250-km leaky
           line, at a reduced peak pass-through (CORRIDOR_MAX_PASS,
           MDL_CROSSING_SLOWDOWN);
       (b) a THICKER fence band (MDL_DILATE_ITER) so the front must traverse
           several attenuated cells rather than one;
       (c) slower transport overall (RAD_DIFF_MAX / RAD_ADV_MAX / RAD_ADV_GAIN);
       (d) a longer horizon (SIM_YEARS, plus KOREA_EXTRA_YEARS for the Korea leg)
           so the delayed crossing is still visible in the outputs.
4. NEW DIAGNOSTICS + MAPS.
       * per-gate log line (grid cell, snap distance, habitat at the gate);
       * MDL gate-permeability map (PPR_v9.3_MDL_gates.png);
       * CORRIDOR CROSSING REPORT — which gate delivered PPR into South Korea,
         in which month, and the peak attack rate on its southern probe.
5. TWO LATENT CRASHES FIXED (they only fired with BASE_RASTER="gorani"):
       * `_floor` was used in the NK/SK band branch of plot_map but only defined
         in the single-colormap branch  -> NameError;
       * `risk_nk` was plotted but never built -> NameError.
-----------------------------------------------------------------------------
v9.0 (retained): the static per-month exponential-kernel RISK engine is replaced
by a true mechanistic spatio-temporal EPIDEMIC: infection accumulates through
time, infected cells become new sources, and the front propagates month-by-month
across the landscape (no re-seeding each month).
  dI/dt = ∇·(D(x,y)∇I) − ∇·(v(x,y) I) + β_eff·S·I − γ·I + Source
where diffusion D and advection v are ecological (built from habitat suitability,
conductance, resistance, rivers and least-cost connectivity), the reaction is
local SIR transmission, and the MDL fence / roads act as transport barriers
(Han et al. 2022 φ, δ). The engine is inlined below (movement physics + disease
dynamics + numerical solver), is numerically stable (flux-form finite
differences, upwind advection, adaptive CFL sub-stepping), and is UDE-READY: the
ecological closure is a single swappable function that a learned neural closure
can later replace — no Neural ODE is implemented now.
Retained from v8.1: percent-spread bands (the cumulative attack rate 1−S of the
DYNAMIC epidemic), high-density gorani hotspots, NK habitat-risk underlay, clean
white SK background with light water/sig.shp lines, and the Han landscape
calibration below. Outputs: epidemic time-series, summary maps (max extent /
time-to-arrival / final recovered) and an animated spread sequence.
  * PERCENT-OF-SPREAD BANDS. Each cell is coloured by the % of its local
    population PPR has reached — the cumulative attack rate (I + R = 1 - S,
    since INITIAL_SUSCEPTIBLE = 1.0), binned into 0-20 / 20-40 / 40-60 /
    60-80 / 80-100 % classes. North Korea (source) uses a warm ramp; South
    Korea (reached via the cross-MDL gorani corridors) uses a blue-green ramp.
  * GORANI SPREADER DOTS. Yellow points mark high-density gorani cells (top
    GORANI_DOT_QUANTILE of the blended density surface).
  * NK -> SK cross-MDL narrative preserved and now corridor-explicit.
v8.0 added a landscape-epidemiology calibration layer anchored to the only
published LOCAL parameter-estimation study for a wildlife disease spreading
across the exact same region (Gyeonggi + Gangwon, spanning the DMZ):
  Han J-H, Yoo D-S, Pak S-I, Kim E-T (2022). "Understanding the transmission
  of African swine fever in wild boars of South Korea: A simulation study for
  parameter estimation." Transbound Emerg Dis 69:e1101-e1112.
  doi:10.1111/tbed.14403
WHAT TRANSFERS:
  * How DISCRETE LINEAR LANDSCAPE FEATURES cut spatial transmission:
        roads   delta = 0.362   (~64% blocked)
        rivers  nu    = 0.379   (~62% blocked)   <-- see CAVEAT below
        fences  phi   = 0.531   (~47% blocked; FENCES ARE LEAKY)
    plus a different-habitat factor pi = 0.624 and decay alpha = 0.650 /km.
  * The qualitative finding that the fence was a WEAK barrier. In v9.1 phi is
    the CEILING at the best corridor rather than the value everywhere.
  * The finding that spatial transmission is SLOW (see the decay diagnostic).
WHAT DOES NOT TRANSFER:
  * ASF-specific machinery (case fatality, carcass transmission, boar demography).
  * RIVERS AS A BARRIER. Han's nu applies to WILD BOAR. Gorani are strongly
    water-associated / good swimmers, which is why water is an ATTRACTOR here
    (WATER_ATTRACT, RIV) and HAN_APPLY_RIVER_BARRIER defaults to False.
Author: Saleem Ahmad (base v2.2/v4.2/v5.0/v6.0) — v7.0 hybrid neural
        — v8.0 Han et al. (2022) landscape calibration
        — v9.0 dynamic RAD epidemic
        — v9.1 corridor-gated MDL + slower N->S crossing
        — v9.2 South Korea panels: plain white canvas, water only,
          NK-confined suitability underlay, per-region % floors
        — v9.3 South-Korea county (sigungu) zonal CSV export + county choropleth
Date: August 24, 2026
"""
import sys, io, warnings, os
import csv
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap, PowerNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from scipy.ndimage import binary_dilation, distance_transform_edt, gaussian_filter, maximum_filter
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra as sp_dijkstra
from osgeo import gdal, ogr as _ogr, osr as _osr
import shutil
# PyTorch imports -- required for neural components
try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
# torchdiffeq for Neural ODE integration
# Install: pip install torchdiffeq --break-system-packages
try:
    from torchdiffeq import odeint
    TORCHDIFFEQ_AVAILABLE = True
except ImportError:
    TORCHDIFFEQ_AVAILABLE = False
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")
gdal.UseExceptions()
LOG_FILE = Path("ppr_v9_10_run.log")
logging.basicConfig(
    level=logging.DEBUG,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE, encoding="utf-8"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)
def warp_raster_to_grid(path, gt, nr, nc, wkt, resample=gdal.GRA_Average):
    """Warp any raster onto a target grid (gt/nr/nc/wkt); nodata -> NaN.
    Returns an all-NaN array if the file is missing or the warp fails."""
    if not Path(str(path)).exists():
        return np.full((nr, nc), np.nan, np.float32)
    try:
        wo = gdal.WarpOptions(
            format="MEM", dstSRS=wkt,
            outputBounds=(gt[0], gt[3] + nr * gt[5], gt[0] + nc * gt[1], gt[3]),
            width=nc, height=nr, resampleAlg=resample)
        dsw = gdal.Warp("", str(path), options=wo)
        if dsw is None:
            return np.full((nr, nc), np.nan, np.float32)
        a = dsw.GetRasterBand(1).ReadAsArray().astype(np.float32)
        nd = dsw.GetRasterBand(1).GetNoDataValue()
        if nd is not None:
            a[a == nd] = np.nan
        a[~np.isfinite(a)] = np.nan
        dsw = None
        return a
    except Exception as e:
        logger.warning(f"  ⚠ warp_raster_to_grid failed for {path}: {e}")
        return np.full((nr, nc), np.nan, np.float32)
# ══════════════════════════════════════════════════════════════════════════════
# INLINED DMZ CORRIDOR MODULE (v9.1) — water-deer crossing gates in the MDL
# ------------------------------------------------------------------------------
# Turns the ecologically justified corridor table into three concrete things the
# RAD engine consumes:
#   1. a PER-CELL MDL permeability field (the gates),
#   2. seed entries on the NK side of each gate (fallback seeding only),
#   3. diagnostics — which gate delivered PPR into South Korea, and when.
# Kept in-file so the script still runs as a single paste-able cell.
# ══════════════════════════════════════════════════════════════════════════════
# (name, lat, lon, weight, ecological basis).  weight = relative permeability of
# that gate to water deer (1.00 = best crossing on the line).
DMZ_CORRIDOR_CROSSINGS: List[tuple] = [
    ("Cheorwon",         38.24, 127.31, 1.00,
     "Iron Triangle — flat lowland agriculture + wetlands; maximum water-deer habitat"),
    ("Paju_Imjin",       37.95, 126.73, 0.85,
     "Imjin reed marsh / Han estuary; floodgate gaps; strong swimmers"),
    ("Yeoncheon_Hantan", 38.02, 127.08, 0.80,
     "Continuous Hantan River riparian path; high local water-deer density"),
    ("Hwacheon_Bukhan",  38.25, 127.68, 0.70,
     "Bukhan River as a direct water route despite steep surrounding terrain"),
    ("Yanggu",           38.28, 128.00, 0.70,
     "Punchbowl basin supports populations; high rim limits lateral spread"),
    ("Goseong_Coast",    38.60, 128.43, 0.60,
     "Lowland coastal strip — suitable but narrow"),
    # The weight table also lists Inje_Seorak (0.40, steep rocky Taebaek —
    # suboptimal for water deer, suits roe deer / goral), but NO COORDINATE was
    # supplied for it. It is kept here with lat/lon = None so the loader logs it
    # loudly and skips it rather than inventing a location. Put the real
    # crossing lat/lon in to activate it.
    ("Inje_Seorak",      None, None, 0.40,
     "Steep rocky Taebaek terrain — suboptimal for water deer (COORDINATE MISSING)"),
]
CORRIDOR_CRS_EPSG = 4326          # the table above is WGS84 lat/lon
def make_lonlat_transform(wkt: str, epsg: int = CORRIDOR_CRS_EPSG):
    """Transform lat/lon (EPSG:4326) -> model grid CRS. Axis order is forced to
    traditional GIS order so TransformPoint takes (lon, lat) on GDAL 2 and 3."""
    src = _osr.SpatialReference(); src.ImportFromEPSG(int(epsg))
    dst = _osr.SpatialReference(); dst.ImportFromWkt(wkt)
    for srs in (src, dst):
        try:
            srs.SetAxisMappingStrategy(_osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:
            pass                                   # GDAL 2: already lon/lat
    return _osr.CoordinateTransformation(src, dst)
def make_grid_to_lonlat_transform(wkt: str, epsg: int = CORRIDOR_CRS_EPSG):
    """INVERSE of make_lonlat_transform: model grid CRS -> (lon, lat) WGS84.
    Used by the v9.3 county export to report a lon/lat centroid per county."""
    src = _osr.SpatialReference(); src.ImportFromWkt(wkt)
    dst = _osr.SpatialReference(); dst.ImportFromEPSG(int(epsg))
    for srs in (src, dst):
        try:
            srs.SetAxisMappingStrategy(_osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:
            pass
    return _osr.CoordinateTransformation(src, dst)
def point_cell(x: float, y: float, gt):
    """(row, col) of a map point. Config.POINT_TO_CELL = "floor" (default): the cell that contains the point, the
    convention of every raster of the model; "round": nearest cell corner, up to half a cell off (runs before 26 Nov)."""
    fc = (float(x) - gt[0]) / gt[1]
    fr = (gt[3] - float(y)) / abs(gt[5])
    if str(getattr(Config, "POINT_TO_CELL", "floor")).lower() == "floor":
        return int(math.floor(fr)), int(math.floor(fc))
    return int(round(fr)), int(round(fc))


def lonlat_to_rc(lon: float, lat: float, gt, nr: int, nc: int, tfm):
    """(lon, lat) -> (row, col, inside, x_map, y_map) on the model grid."""
    x, y = tfm.TransformPoint(float(lon), float(lat))[:2]
    r, c = point_cell(x, y, gt)
    return r, c, (0 <= r < nr and 0 <= c < nc), x, y
def build_corridors(mdl_boundary, gt, wkt, nr, nc, cell_km,
                    kp=None, crossings=None, snap_max_km=25.0,
                    gate_radius_km=8.0, logger=None):
    """Reproject each corridor, snap it onto the nearest digitized MDL cell, and
    return a list of dicts: name, weight, lat, lon, row, col, snap_km, note."""
    def _log(msg, warn=False):
        if logger is None:
            print(msg)
        elif warn:
            logger.warning(msg)
        else:
            logger.info(msg)
    crossings = list(crossings if crossings is not None else DMZ_CORRIDOR_CROSSINGS)
    tfm = make_lonlat_transform(wkt)
    if mdl_boundary is not None and bool(np.any(mdl_boundary)):
        _, idx = distance_transform_edt(~mdl_boundary.astype(bool), return_indices=True)
        ir, ic = idx[0], idx[1]
    else:
        ir = ic = None
        _log("  ⚠ corridors: the MDL raster is empty — corridor points cannot be "
             "snapped onto the fence line.", warn=True)
    out: List[Dict] = []
    for name, lat, lon, w, note in crossings:
        if lat is None or lon is None:
            _log(f"  ⚠ corridor {name}: NO COORDINATE SUPPLIED (weight {w:.2f}) — skipped. "
                 f"Add a lat/lon to DMZ_CORRIDOR_CROSSINGS to activate this gate.", warn=True)
            continue
        r, c, inside, x, y = lonlat_to_rc(lon, lat, gt, nr, nc, tfm)
        if not inside:
            _log(f"  ⚠ corridor {name}: ({lat:.3f}N,{lon:.3f}E) -> grid ({r},{c}) falls "
                 f"OUTSIDE the {nr}x{nc} grid — skipped (check raster extent / CRS).", warn=True)
            continue
        snap_km = 0.0
        if ir is not None:
            rs, cs = int(ir[r, c]), int(ic[r, c])
            snap_km = math.hypot(rs - r, cs - c) * cell_km
            if snap_km > snap_max_km:
                _log(f"  ⚠ corridor {name}: nearest MDL cell is {snap_km:.1f} km away "
                     f"(> {snap_max_km:.0f} km). Snapping anyway, but verify this "
                     f"coordinate against MDL.shp before publishing.", warn=True)
            r, c = rs, cs
        in_domain = True if kp is None else bool(kp[r, c])
        out.append(dict(name=name, weight=float(w), lat=float(lat), lon=float(lon),
                        row=int(r), col=int(c), snap_km=float(snap_km),
                        in_domain=in_domain, note=note))
        _log(f"  · {name:<17s} w={w:.2f}  ({lat:.3f}N,{lon:.3f}E) -> cell ({r},{c})  "
             f"snap={snap_km:5.1f} km  in-domain={in_domain}   [{note}]")
    if not out:
        _log("  ⚠ corridors: NONE resolved — the MDL fence stays uniform (Han phi).", warn=True)
    else:
        _log(f"  ✓ {len(out)} DMZ corridor gates resolved on the MDL line "
             f"(cell size {cell_km:.2f} km).")
        if cell_km > gate_radius_km:
            _log(f"  ⚠ cell size ({cell_km:.1f} km) exceeds the gate radius "
                 f"({gate_radius_km:.0f} km): gates are barely resolved and neighbouring "
                 f"corridors may merge into one cell. Use the finer gorani grid "
                 f"(BASE_RASTER='gorani') for corridor-level conclusions.", warn=True)
    return out
def build_mdl_fence_field(corridors, mdl_boundary, shape, cell_km,
                          gate_radius_km=8.0, max_pass=0.35, sealed=0.02,
                          slowdown=1.0, logger=None):
    """Per-cell MDL pass-through in (0,1]: SEALED everywhere except a Gaussian
    gate at each corridor whose peak is `max_pass * weight * slowdown`.
    Consumed by the engine's build_barrier() through env['MDL_FENCE_FIELD'].
    Returns None when there are no corridors (caller keeps the uniform phi)."""
    if not corridors:
        return None
    mdl = mdl_boundary.astype(bool)
    ff = np.full(shape, float(sealed), dtype=np.float64)
    sig_cells = max(gate_radius_km / 2.0, 1e-6) / max(cell_km, 1e-6)
    for co in corridors:
        pt = np.ones(shape, dtype=bool)
        pt[co["row"], co["col"]] = False
        d_cells = distance_transform_edt(pt)
        taper = np.exp(-0.5 * (d_cells / sig_cells) ** 2)
        peak = float(np.clip(max_pass * co["weight"] * slowdown, sealed, 1.0))
        ff = np.maximum(ff, sealed + (peak - sealed) * taper)
    ff = np.where(mdl, np.clip(ff, 1e-5, 1.0), 1.0)
    if logger is not None:
        on = ff[mdl]
        open_frac = float((on > sealed * 2.0).mean()) if on.size else 0.0
        logger.info(f"  ✓ MDL fence field: {int(mdl.sum()):,} fence cells, "
                    f"{100*open_frac:.1f}% of them inside a corridor gate; "
                    f"pass-through min={on.min():.3f} max={on.max():.3f} "
                    f"(sealed={sealed}, ceiling={max_pass}×weight×{slowdown}).")
    return ff.astype(np.float64)
def corridor_seed_entries(corridors, nk_mask, offset_cells=2,
                          seed_weight_from_corridor=True):
    """(row, col, label, 0, weight) entries just NORTH of each gate, for the RAD
    seeder. Used only as the fallback when no outbreak shapefile drives the run."""
    ent = []
    nr = nk_mask.shape[0]
    for co in corridors:
        r = max(0, co["row"] - int(offset_cells))
        while r > 0 and not nk_mask[r, co["col"]]:
            r -= 1                                  # walk north into NK
        if 0 <= r < nr and nk_mask[r, co["col"]]:
            w = co["weight"] if seed_weight_from_corridor else 1.0
            ent.append((int(r), int(co["col"]), co["name"], 0, float(w)))
    return ent
def plot_corridor_markers(ax, corridors, row_offset=0, col_offset=0,
                          label=True, size_scale=90.0, zorder=9.2):
    """Draw the gates; marker area scales with the ecological weight."""
    if not corridors:
        return
    rr = np.array([co["row"] - row_offset for co in corridors], float) + 0.5
    cc = np.array([co["col"] - col_offset for co in corridors], float) + 0.5
    ww = np.array([co["weight"] for co in corridors], float)
    ax.scatter(cc, rr, s=size_scale * (0.35 + 0.65 * ww), marker="v",
               c="#00e5ff", edgecolors="#003b46", linewidths=0.8,
               alpha=0.95, zorder=zorder)
    if label and getattr(Config, "PUB_FIGS", False):
        # 11 Nov (publication figures): full names, staggered labels with leader lines
        _full = {"Paju_Imjin": "Paju–Imjin", "Yeoncheon_Hantan": "Yeoncheon–Hantan",
                 "Hwacheon_Bukhan": "Hwacheon–Bukhan", "Goseong_Coast": "Goseong",
                 "Inje_Seorak": "Inje–Seorak"}
        _ord = list(np.argsort(cc))
        _xs = np.linspace(0.10, 0.80, len(_ord)) if len(_ord) > 1 else [0.45]
        for k, i in enumerate(_ord):
            co = corridors[i]
            ax.annotate(f"{_full.get(co['name'], co['name'])} ({co['weight']:.2f})",
                        (cc[i], rr[i]), xycoords="data", xytext=(_xs[k], 0.93 if k % 2 == 0 else 0.83),
                        textcoords="axes fraction",
                        ha="center", va="center", fontsize=10, color="#003b46",
                        zorder=zorder + 0.1,
                        bbox=dict(fc="white", ec="#003b46", lw=0.6, alpha=0.9, pad=2.0),
                        arrowprops=dict(arrowstyle="-", color="#003b46", lw=0.8,
                                        shrinkA=0, shrinkB=3))
    elif label:
        for co, r_, c_ in zip(corridors, rr, cc):
            ax.annotate(f"{co['name'].split('_')[0]} ({co['weight']:.2f})",
                        (c_, r_), xytext=(4, -8), textcoords="offset points",
                        fontsize=6.5, color="#003b46", zorder=zorder + 0.1,
                        bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.6))
def corridor_crossing_report(corridors, monthly_cum, sk_mask, cell_km,
                             probe_km=15.0, threshold=0.05,
                             row_offset=0, col_offset=0, month_offset=0,
                             logger=None):
    """For each gate: the first month the South-Korean side immediately below it
    exceeds `threshold` cumulative attack — i.e. WHICH corridor delivered PPR
    into South Korea, and WHEN."""
    rows = []
    nr, nc = sk_mask.shape
    rad = max(1, int(round(probe_km / max(cell_km, 1e-6))))
    for co in corridors:
        r0 = co["row"] - row_offset
        c0 = co["col"] - col_offset
        rs, re = max(0, r0), min(nr, r0 + 2 * rad)          # probe SOUTH of the gate
        cs, ce = max(0, c0 - rad), min(nc, c0 + rad + 1)
        probe = np.zeros((nr, nc), bool)
        if re > rs and ce > cs:
            probe[rs:re, cs:ce] = True
        probe &= sk_mask
        month, peak = None, 0.0
        if probe.any():
            for m, cum in enumerate(monthly_cum):
                v = float(np.nanmean(np.nan_to_num(cum)[probe]))
                peak = max(peak, v)
                if month is None and v >= threshold:
                    month = m + 1 + int(month_offset)
        rows.append(dict(name=co["name"], weight=co["weight"],
                         arrival_month=month, peak_attack=peak,
                         probe_cells=int(probe.sum())))
    if logger is not None:
        logger.info(f"  [CORRIDOR CROSSING REPORT] South-Korean side within {probe_km:.0f} km "
                    f"of each gate (arrival = mean attack ≥ {100*threshold:.0f}%)")
        for r in sorted(rows, key=lambda d: (d["arrival_month"] is None,
                                             d["arrival_month"] or 0)):
            am = f"month {r['arrival_month']:>3d}" if r["arrival_month"] else "not reached"
            logger.info(f"    {r['name']:<17s} w={r['weight']:.2f}  {am}  "
                        f"peak attack={100*r['peak_attack']:5.1f}%  "
                        f"({r['probe_cells']} SK probe cells)")
    return rows
# ══════════════════════════════════════════════════════════════════════════════
# INLINED SOUTH-KOREA COUNTY (SIGUNGU) ZONAL-EXPORT MODULE — NEW IN v9.3
# ------------------------------------------------------------------------------
# The model produces raster fields; veterinary decisions are made per COUNTY.
# This module is the bridge. It rasterizes an administrative polygon layer
# (sig.shp / sigungu / GADM KOR_2 / any adm-2 layer) onto the model grid at a
# finer sub-resolution, builds AREA WEIGHTS for every county × model-cell pair,
# and turns the monthly cumulative-attack stack into:
#     · a per-county summary table (arrival month, final attack, area reached…),
#     · a counties × months matrix,
#     · a tidy long table,
#     · rasters for a county choropleth.
#
# THREE THINGS THIS MODULE IS CAREFUL ABOUT, because each one silently corrupts
# a county table if ignored:
#   1. ENCODING. Korean shapefiles are CP949/EUC-KR far more often than UTF-8.
#      Reading them as UTF-8 gives mojibake county names that will not join to
#      any Korean database. The encoding is chosen by scoring how much valid
#      Hangul each candidate decodes, and CSVs are written UTF-8-with-BOM so
#      Excel renders them correctly on a Korean Windows machine.
#   2. SUB-CELL AREA. Rasterizing 250 counties straight onto a 10-15 km grid
#      throws away every county smaller than a cell and biases the rest toward
#      whichever polygon happened to own the cell centre. Rasterizing at
#      OVERSAMPLE× and counting sub-cells gives each county a fractional weight
#      in every cell it touches.
#   3. HONESTY ABOUT RESOLUTION. Counties still too small to catch a single
#      sub-cell are sampled at their point-on-surface and FLAGGED in the CSV
#      (sampling="centroid"), never silently merged away.
# ══════════════════════════════════════════════════════════════════════════════
# Attribute names seen in the wild: the Korean national sigungu layer (sig.shp),
# GADM level-2, and a few other admin-2 products.
SK_NAME_KR_FIELDS = ["SIG_KOR_NM", "SIGUNGU_NM", "SIG_NM", "NL_NAME_2", "ADM_NM",
                     "SGG_NM", "KOR_NM", "NAME_KOR", "CTY_KOR_NM"]
SK_NAME_EN_FIELDS = ["SIG_ENG_NM", "NAME_2", "ENG_NM", "VARNAME_2", "ADM_ENG_NM",
                     "SIGUNGU_EN", "NAME"]
SK_CODE_FIELDS    = ["SIG_CD", "SIGUNGU_CD", "SGG_CD", "ADM_CD", "ADM_SECT_C",
                     "GID_2", "CC_2", "HC_2", "CODE"]
SK_PROV_FIELDS    = ["CTP_KOR_NM", "CTPRVN_NM", "SIDO_NM", "NAME_1", "SIDO",
                     "PROVINCE", "CTP_ENG_NM"]
def _hangul_score(s: str) -> float:
    """How Korean does this decoded string look? Hangul syllables score +2,
    replacement characters and private-use junk score -5. Used to pick the
    shapefile's DBF encoding without asking the user."""
    if not s:
        return 0.0
    han = sum(1 for ch in s if "가" <= ch <= "힣")
    bad = sum(1 for ch in s if ch == "�" or 0xE000 <= ord(ch) <= 0xF8FF)
    return 2.0 * han - 5.0 * bad
def _pick_field(field_names: Sequence[str], candidates: Sequence[str]):
    """Case-insensitive first match of `candidates` inside `field_names`."""
    up = {str(n).upper(): n for n in field_names}
    for cand in candidates:
        if cand.upper() in up:
            return up[cand.upper()]
    return None
def detect_shape_encoding(shp_path, candidates=("", "CP949", "EUC-KR", "UTF-8"),
                          probe=40, logger=None):
    """Return the SHAPE_ENCODING value that decodes the most valid Hangul.
    "" means 'let GDAL use the .cpg / DBF language driver ID'."""
    best, best_score = "", -1e18
    for enc in candidates:
        old = gdal.GetConfigOption("SHAPE_ENCODING")
        score = -1e17
        try:
            gdal.SetConfigOption("SHAPE_ENCODING", enc)
            ds = _ogr.Open(str(shp_path))
            if ds is not None:
                lyr = ds.GetLayer()
                defn = lyr.GetLayerDefn()
                names = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
                fld = _pick_field(names, SK_NAME_KR_FIELDS)
                score = 0.0
                if fld is not None:
                    lyr.ResetReading()
                    for k, ft in enumerate(lyr):
                        if k >= probe:
                            break
                        try:
                            score += _hangul_score(ft.GetFieldAsString(fld))
                        except Exception:
                            score -= 1.0
                ds = None
        except Exception:
            score = -1e17
        finally:
            gdal.SetConfigOption("SHAPE_ENCODING", old)
        if score > best_score:
            best_score, best = score, enc
    if logger is not None:
        logger.info(f"  · DBF encoding auto-detected: "
                    f"{best or '(GDAL default / .cpg)'} (Hangul score {best_score:.0f})")
    return best
class CountyZoneError(RuntimeError):
    """Raised when the admin-2 layer cannot be turned into usable zones. Carries
    a `reason` code so the caller can print advice specific to what went wrong
    instead of a catch-all 'not found or could not be rasterized'."""
    def __init__(self, reason, message):
        super().__init__(message)
        self.reason = reason
def build_county_zones(shp_path, wkt, gt, nr, nc, oversample=6,
                       max_subcells=120_000_000, encoding=None,
                       drop_north_of_row=None, logger=None):
    """Rasterize an admin-2 layer onto the model grid and return AREA-WEIGHTED
    county zones.
    Returns
    -------
    dict with:
      zones        : list of per-county dicts. Each carries `rows`, `cols`
                     (model-grid indices) and `weights` (sub-cell counts), plus
                     the attributes needed for the CSV.
      zone_raster  : (nr, nc) int32, the DOMINANT county id per model cell
                     (0 = no county). Used for the choropleth.
      fields       : which shapefile attributes were used for each column.
      encoding     : the DBF encoding actually used.
      oversample   : the sub-resolution factor actually used.
      stats        : counts of what was read, kept and skipped, and why.
    Implementation note (v9.3.1): the working copy of the layer is an OGR
    Memory layer this function CREATES and fills itself, rather than the output
    of gdal.VectorTranslate. VectorTranslate is convenient but its result is not
    reliably writable across GDAL builds, and adding the zone-id field to it is
    the one step that must not fail. Reprojection, geometry filtering and the
    id assignment therefore all happen here, per feature, where each can be
    counted and reported.
    """
    def _log(msg, warn=False):
        if logger is None:
            print(msg)
        elif warn:
            logger.warning(msg)
        else:
            logger.info(msg)
    os_ = max(1, int(oversample))
    while os_ > 1 and (nr * os_) * (nc * os_) > int(max_subcells):
        os_ -= 1
    enc = encoding if encoding is not None else detect_shape_encoding(shp_path, logger=logger)
    old_enc = gdal.GetConfigOption("SHAPE_ENCODING")
    gdal.SetConfigOption("SHAPE_ENCODING", enc or "")
    zones: List[Dict] = []
    zone_raster = np.zeros((nr, nc), np.int32)
    fields: Dict = {}
    stats = dict(features_read=0, kept=0, skipped_not_polygon=0,
                 skipped_north_of_mdl=0, skipped_outside_grid=0,
                 skipped_no_geometry=0)
    try:
        src_ds = _ogr.Open(str(shp_path))
        if src_ds is None:
            raise CountyZoneError("unreadable", f"OGR could not open {shp_path} "
                                                f"(missing .dbf/.shx, or not a vector file)")
        src_lyr = src_ds.GetLayer()
        stats["features_read"] = int(src_lyr.GetFeatureCount())
        if stats["features_read"] == 0:
            raise CountyZoneError("empty", f"{Path(str(shp_path)).name} contains 0 features")
        # ── CRS ─────────────────────────────────────────────────────────────
        src_srs = src_lyr.GetSpatialRef()
        dst_srs = _osr.SpatialReference(); dst_srs.ImportFromWkt(wkt)
        for _s in (src_srs, dst_srs):
            if _s is not None:
                try:
                    _s.SetAxisMappingStrategy(_osr.OAMS_TRADITIONAL_GIS_ORDER)
                except AttributeError:
                    pass
        tfm = None
        if src_srs is None:
            _log(f"  ⚠ {Path(str(shp_path)).name} has NO .prj — assuming it is already in "
                 f"the model CRS. If the county names come out right but every county "
                 f"reports 'no modelled cells', this is why: supply the .prj.", warn=True)
        elif not src_srs.IsSame(dst_srs):
            tfm = _osr.CoordinateTransformation(src_srs, dst_srs)
        # ── attribute mapping ───────────────────────────────────────────────
        defn = src_lyr.GetLayerDefn()
        fnames = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
        f_kr = _pick_field(fnames, SK_NAME_KR_FIELDS)
        f_en = _pick_field(fnames, SK_NAME_EN_FIELDS)
        f_cd = _pick_field(fnames, SK_CODE_FIELDS)
        f_pv = _pick_field(fnames, SK_PROV_FIELDS)
        fields = dict(name_kr=f_kr, name_en=f_en, code=f_cd, province=f_pv,
                      available=fnames)
        _log(f"  · {stats['features_read']} features; attribute mapping: "
             f"name_kr={f_kr}  name_en={f_en}  code={f_cd}  province={f_pv}")
        if f_kr is None and f_en is None:
            _log(f"  ⚠ no recognizable county-name field in this layer. Fields present: "
                 f"{fnames}. Rows will be labelled by index — add the right field name to "
                 f"SK_NAME_KR_FIELDS / SK_NAME_EN_FIELDS at the top of the script.", warn=True)
        # ── our own writable Memory layer ───────────────────────────────────
        mem_drv = _ogr.GetDriverByName("Memory")
        if mem_drv is None:
            raise CountyZoneError("no_mem_driver", "OGR 'Memory' driver unavailable")
        mds = mem_drv.CreateDataSource("ppr_county_zones")
        mlyr = mds.CreateLayer("zones", srs=dst_srs, geom_type=_ogr.wkbMultiPolygon)
        zid_field = "ZID_"
        mlyr.CreateField(_ogr.FieldDefn(zid_field, _ogr.OFTInteger))
        mdefn = mlyr.GetLayerDefn()
        tfm_ll = None
        try:
            tfm_ll = make_grid_to_lonlat_transform(wkt)
        except Exception:
            pass
        gx0, gy0 = gt[0], gt[3]
        gx1 = gt[0] + nc * gt[1]
        gy1 = gt[3] + nr * gt[5]
        gxmin, gxmax = min(gx0, gx1), max(gx0, gx1)
        gymin, gymax = min(gy0, gy1), max(gy0, gy1)
        zid = 0
        sample_names: List[str] = []
        src_lyr.ResetReading()
        for ft in src_lyr:
            geom = ft.GetGeometryRef()
            if geom is None:
                stats["skipped_no_geometry"] += 1
                continue
            g = geom.Clone()
            if tfm is not None:
                try:
                    g.Transform(tfm)
                except Exception:
                    stats["skipped_no_geometry"] += 1
                    continue
            if "POLYGON" not in (g.GetGeometryName() or "").upper():
                stats["skipped_not_polygon"] += 1
                continue
            env = g.GetEnvelope()                    # (xmin, xmax, ymin, ymax)
            if env[1] < gxmin or env[0] > gxmax or env[3] < gymin or env[2] > gymax:
                stats["skipped_outside_grid"] += 1
                continue
            try:
                pt = g.PointOnSurface()
                cx, cy = float(pt.GetX()), float(pt.GetY())
            except Exception:
                cx, cy = 0.5 * (env[0] + env[1]), 0.5 * (env[2] + env[3])
            crow = (gt[3] - cy) / abs(gt[5])
            if drop_north_of_row is not None and crow < float(drop_north_of_row):
                stats["skipped_north_of_mdl"] += 1
                continue
            zid += 1
            nf = _ogr.Feature(mdefn)
            nf.SetGeometry(g)
            nf.SetField(zid_field, zid)
            mlyr.CreateFeature(nf)
            nf = None
            try:
                area_km2 = float(g.GetArea()) / 1e6
            except Exception:
                area_km2 = float("nan")
            lon = lat = float("nan")
            if tfm_ll is not None:
                try:
                    lon, lat = tfm_ll.TransformPoint(cx, cy)[:2]
                except Exception:
                    pass
            nm_kr = ft.GetFieldAsString(f_kr) if f_kr else ""
            nm_en = ft.GetFieldAsString(f_en) if f_en else ""
            if len(sample_names) < 6:
                sample_names.append(nm_kr or nm_en or f"#{zid}")
            zones.append(dict(
                zid=zid, name_kr=nm_kr, name_en=nm_en,
                code=(ft.GetFieldAsString(f_cd) if f_cd else ""),
                province=((ft.GetFieldAsString(f_pv) if f_pv else "")
                          or province_from_sig_cd(ft.GetFieldAsString(f_cd) if f_cd else "")[0]),  # v9.10
                province_kr=province_from_sig_cd(ft.GetFieldAsString(f_cd) if f_cd else "")[1],
                poly_area_km2=area_km2, cx=cx, cy=cy, lon=lon, lat=lat))
        stats["kept"] = zid
        src_ds = None
        _log(f"  · features read={stats['features_read']}, kept={stats['kept']}; "
             f"skipped: not-polygon={stats['skipped_not_polygon']}, "
             f"outside-grid={stats['skipped_outside_grid']}, "
             f"north-of-MDL={stats['skipped_north_of_mdl']}, "
             f"no-geometry={stats['skipped_no_geometry']}")
        if sample_names:
            _log(f"  · first county names as decoded: {', '.join(sample_names)}   "
                 f"(if these look like mojibake, set Config.SK_COUNTY_SHP_ENCODING)")
        if zid == 0:
            if stats["skipped_outside_grid"] == stats["features_read"]:
                raise CountyZoneError(
                    "outside_grid",
                    "every county polygon falls OUTSIDE the model grid extent — the "
                    "layer's CRS and the base raster's CRS disagree, or this is the "
                    "wrong shapefile. Check the [DATA INVENTORY] extents in the log.")
            if stats["skipped_north_of_mdl"] == stats["features_read"]:
                raise CountyZoneError(
                    "all_north",
                    "every county polygon sits NORTH of the MDL row, so the "
                    "South-Korea-only filter removed them all. Set "
                    "Config.SK_COUNTY_DROP_NORTH_OF_MDL = False if this layer is "
                    "genuinely South Korea and the MDL row is misplaced.")
            if stats["skipped_not_polygon"] == stats["features_read"]:
                raise CountyZoneError(
                    "not_polygon",
                    "this layer contains no polygons (points or lines only) — zonal "
                    "statistics need an admin-2 POLYGON layer such as sig.shp or "
                    "GADM KOR_2.")
            raise CountyZoneError("no_usable", "no usable county polygons after filtering")
        # ── rasterize at the sub-resolution ─────────────────────────────────
        gt_os = (gt[0], gt[1] / os_, 0.0, gt[3], 0.0, gt[5] / os_)
        mem = gdal.GetDriverByName("MEM").Create("", nc * os_, nr * os_, 1, gdal.GDT_Int32)
        mem.SetGeoTransform(gt_os)
        mem.SetProjection(wkt)
        mem.GetRasterBand(1).Fill(0)
        err = gdal.RasterizeLayer(mem, [1], mlyr, options=[f"ATTRIBUTE={zid_field}"])
        if err not in (0, None, gdal.CE_None):
            _log(f"  ⚠ RasterizeLayer returned {err}; continuing with whatever burned.", warn=True)
        zarr = mem.GetRasterBand(1).ReadAsArray()
        mem = None
        mds = None
        sub_r, sub_c = np.nonzero(zarr)
        if sub_r.size:
            vals = zarr[sub_r, sub_c].astype(np.int64)
            flat = (sub_r // os_).astype(np.int64) * int(nc) + (sub_c // os_).astype(np.int64)
            order = np.lexsort((flat, vals))
            vals_s, flat_s = vals[order], flat[order]
            zids = np.arange(1, zid + 1, dtype=np.int64)
            starts = np.searchsorted(vals_s, zids, side="left")
            ends = np.searchsorted(vals_s, zids, side="right")
        else:
            starts = ends = np.zeros(zid, np.int64)
            _log("  ⚠ the county layer burned ZERO sub-cells. Every county will fall back "
                 "to centroid sampling. This usually means the grid and the layer do not "
                 "actually overlap.", warn=True)
        best_cnt = np.zeros(nr * nc, np.float64)
        zr_flat = zone_raster.ravel()
        n_centroid = n_zonal = n_lost = 0
        for z in zones:
            i0, i1 = int(starts[z["zid"] - 1]), int(ends[z["zid"] - 1])
            if i1 > i0:
                fz, wz = np.unique(flat_s[i0:i1], return_counts=True)
                z["cells_flat"] = fz.astype(np.int64)
                z["weights"] = wz.astype(np.float64)
                z["sampling"] = "zonal"
                n_zonal += 1
            else:
                rr, cc = point_cell(z["cx"], z["cy"], gt)
                if 0 <= rr < nr and 0 <= cc < nc:
                    z["cells_flat"] = np.array([rr * nc + cc], np.int64)
                    z["weights"] = np.array([1.0], np.float64)
                    z["sampling"] = "centroid"
                    n_centroid += 1
                else:
                    z["cells_flat"] = np.zeros(0, np.int64)
                    z["weights"] = np.zeros(0, np.float64)
                    z["sampling"] = "outside-grid"
                    n_lost += 1
            z["rows"] = (z["cells_flat"] // int(nc)).astype(np.int64)
            z["cols"] = (z["cells_flat"] % int(nc)).astype(np.int64)
            if z["cells_flat"].size and z["sampling"] == "zonal":
                better = z["weights"] > best_cnt[z["cells_flat"]]
                idx = z["cells_flat"][better]
                zr_flat[idx] = z["zid"]
                best_cnt[idx] = z["weights"][better]
        zone_raster = zr_flat.reshape(nr, nc)
        stats.update(zonal=n_zonal, centroid=n_centroid, outside=n_lost)
        cell_km = abs(gt[1]) / 1000.0
        _log(f"  ✓ county zones built: {len(zones)} polygons on a {nc}×{nr} grid at "
             f"{os_}× sub-resolution ({cell_km / os_:.2f} km sub-cells vs {cell_km:.2f} km "
             f"model cells); area-weighted={n_zonal}, centroid-sampled={n_centroid}, "
             f"outside grid={n_lost}.")
        if n_centroid:
            _log(f"  ⚠ {n_centroid} counties are SMALLER THAN ONE SUB-CELL and were sampled "
                 f"at a single point (flagged sampling='centroid' in the CSV). At "
                 f"{cell_km:.1f} km cells this is expected for the metropolitan gu "
                 f"districts. For county-level conclusions use BASE_RASTER='gorani' or a "
                 f"lower FACTOR / higher SK_COUNTY_OVERSAMPLE.", warn=True)
    finally:
        gdal.SetConfigOption("SHAPE_ENCODING", old_enc)
    return dict(zones=zones, zone_raster=zone_raster, fields=fields,
                encoding=enc, oversample=os_, stats=stats,
                shp_path=str(shp_path))
def _sim_month_to_calendar(global_month_1based, start_year=None, start_month=1):
    """Nominal calendar label for a 1-based simulation month. The model has no
    absolute clock — this is a labelling convenience anchored on
    Config.SIM_START_YEAR / SIM_START_MONTH, so state it as such in any figure."""
    if start_year is None or global_month_1based is None:
        return ""
    t = int(start_year) * 12 + (int(start_month) - 1) + (int(global_month_1based) - 1)
    return f"{t // 12:04d}-{t % 12 + 1:02d}"
def county_attack_table(zones, monthly_cum, domain_mask, cell_km,
                        arrival_threshold=0.05, month_offset=0,
                        row_offset=0, col_offset=0, sk_mask=None,
                        monthly_I=None, corridors=None, oversample=None,
                        start_year=None, start_month=1, logger=None):
    """Per-county epidemic summary + the month-by-month series behind it.
    Parameters
    ----------
    zones        : output of build_county_zones()["zones"] (FULL-GRID indices).
    monthly_cum  : list of (nr_c, nc_c) cumulative attack-rate arrays (0-1) for
                   the grid the epidemic was run on — which may be a CROP of the
                   full grid, hence row_offset / col_offset.
    domain_mask  : (nr_c, nc_c) bool. Cells the epidemic actually covers.
    month_offset : months elapsed BEFORE monthly_cum[0] (the Korea leg starts
                   when the China run reaches the border), so reported months are
                   on one continuous timeline.
    Returns (rows, series) where `series` maps zid -> dict of per-month arrays.
    """
    thr = float(arrival_threshold)
    NM = len(monthly_cum)
    nrc, ncc = domain_mask.shape
    # NOTE: the monthly stack is NOT pre-converted to float64. A 750x750 grid over
    # a 300-month adaptive run is ~675 MB as float32; casting the whole stack up
    # front would double that for no benefit. Each county indexes the cells it
    # needs and casts only those (a few dozen values per month).
    sub_per_cell = (float(oversample) ** 2) if oversample else None
    rows: List[Dict] = []
    series: Dict[int, Dict] = {}
    for z in zones:
        r = np.asarray(z["rows"], np.int64) - int(row_offset)
        c = np.asarray(z["cols"], np.int64) - int(col_offset)
        w = np.asarray(z["weights"], np.float64)
        rec = dict(
            sig_cd=z.get("code", ""), county_kr=z.get("name_kr", ""),
            county_en=z.get("name_en", ""), province=z.get("province", ""),
            province_kr=z.get("province_kr", ""),                    # v9.10
            poly_area_km2=z.get("poly_area_km2", float("nan")),
            centroid_lat=z.get("lat", float("nan")),
            centroid_lon=z.get("lon", float("nan")),
            sampling=z.get("sampling", ""), zid=z["zid"])
        ok = (r >= 0) & (r < nrc) & (c >= 0) & (c < ncc)
        r, c, w = r[ok], c[ok], w[ok]
        if r.size:
            keep = domain_mask[r, c]
            r, c, w = r[keep], c[keep], w[keep]
        if r.size == 0 or float(w.sum()) <= 0:
            rec.update(status="no modelled cells", model_cells=0, model_area_km2=0.0,
                       arrival_month=None, arrival_month_first_cell=None,
                       final_mean_attack_pct=float("nan"),
                       final_max_attack_pct=float("nan"),
                       final_pct_area_reached=float("nan"),
                       peak_mean_infected_pct=float("nan"), peak_infected_month=None,
                       pct_cells_south_of_mdl=float("nan"))
            rows.append(rec)
            series[z["zid"]] = dict(mean=np.zeros(NM), frac=np.zeros(NM), mx=np.zeros(NM))
            continue
        wsum = float(w.sum())
        n_cells = int(np.unique(r * ncc + c).size)
        s_mean = np.zeros(NM, np.float64)
        s_frac = np.zeros(NM, np.float64)
        s_max = np.zeros(NM, np.float64)
        s_inf = np.zeros(NM, np.float64) if monthly_I else None
        for m in range(NM):
            v = np.nan_to_num(np.asarray(monthly_cum[m][r, c], dtype=np.float64), nan=0.0)
            s_mean[m] = float(np.dot(w, v) / wsum)
            s_frac[m] = float(np.dot(w, (v >= thr).astype(np.float64)) / wsum)
            s_max[m] = float(v.max())
            if s_inf is not None:
                s_inf[m] = float(np.dot(
                    w, np.nan_to_num(np.asarray(monthly_I[m][r, c], dtype=np.float64), nan=0.0)) / wsum)
        i_mean = np.where(s_mean >= thr)[0]
        i_cell = np.where(s_max >= thr)[0]
        a_mean = int(i_mean[0]) + 1 + int(month_offset) if i_mean.size else None
        a_cell = int(i_cell[0]) + 1 + int(month_offset) if i_cell.size else None
        pk_month = pk_val = None
        if s_inf is not None and s_inf.size:
            k = int(np.argmax(s_inf))
            pk_month, pk_val = k + 1 + int(month_offset), float(s_inf[k])
        south_pct = float("nan")
        if sk_mask is not None:
            south_pct = 100.0 * float(np.dot(w, sk_mask[r, c].astype(np.float64)) / wsum)
        gate_name, gate_km, gate_w = "", float("nan"), float("nan")
        if corridors:
            cr = float(np.dot(w, r) / wsum)
            cc_ = float(np.dot(w, c) / wsum)
            best = None
            for co in corridors:
                d = math.hypot(co["row"] - int(row_offset) - cr,
                               co["col"] - int(col_offset) - cc_) * float(cell_km)
                if best is None or d < best[1]:
                    best = (co["name"], d, co["weight"])
            if best:
                gate_name, gate_km, gate_w = best[0], best[1], best[2]
        rec.update(
            status=("reached" if a_mean is not None
                    else ("present (sub-threshold mean)" if a_cell is not None
                          else "not reached")),
            model_cells=n_cells,
            # AREA-WEIGHTED modelled footprint: the sub-cell weights convert
            # straight to km² once you know how many sub-cells a model cell holds.
            # Blank for centroid-sampled counties, where the weight is a token 1
            # and any area derived from it would be meaningless — read
            # poly_area_km2 there instead.
            model_area_km2=(float("nan") if (sub_per_cell is None
                                             or rec["sampling"] == "centroid")
                            else (wsum / sub_per_cell) * float(cell_km) ** 2),
            arrival_month=a_mean,
            arrival_month_first_cell=a_cell,
            arrival_calendar=_sim_month_to_calendar(a_mean, start_year, start_month),
            arrival_year_of_sim=(None if a_mean is None else (a_mean - 1) // 12 + 1),
            final_mean_attack_pct=100.0 * float(s_mean[-1]),
            final_max_attack_pct=100.0 * float(s_max[-1]),
            final_pct_area_reached=100.0 * float(s_frac[-1]),
            peak_mean_infected_pct=(float("nan") if pk_val is None else 100.0 * pk_val),
            peak_infected_month=pk_month,
            nearest_gate=gate_name, nearest_gate_km=gate_km, nearest_gate_weight=gate_w,
            pct_cells_south_of_mdl=south_pct)
        rows.append(rec)
        series[z["zid"]] = dict(mean=s_mean, frac=s_frac, mx=s_max)
    # rank by arrival, unreached counties last
    reached = [r for r in rows if r.get("arrival_month")]
    reached.sort(key=lambda d: (d["arrival_month"], -d.get("final_mean_attack_pct", 0.0)))
    for i, r in enumerate(reached, 1):
        r["rank_arrival"] = i
    first = reached[0]["arrival_month"] if reached else None
    for r in rows:
        r.setdefault("rank_arrival", None)
        am = r.get("arrival_month")
        r["months_after_first_sk_arrival"] = (None if (am is None or first is None)
                                              else am - first)
    rows.sort(key=lambda d: (d["rank_arrival"] is None, d.get("rank_arrival") or 0,
                             str(d.get("sig_cd") or ""), str(d.get("county_kr") or "")))
    if logger is not None:
        n_reach = len(reached)
        logger.info(f"  [COUNTY TABLE] {len(rows)} counties; {n_reach} reached "
                    f"(mean attack ≥ {100*thr:.0f}%) by month "
                    f"{NM + int(month_offset)}.")
        for r in reached[:10]:
            logger.info(f"    {r['rank_arrival']:>3d}. {str(r['county_kr'] or r['county_en'])[:18]:<18s} "
                        f"({r['sig_cd']:<6s} {str(r['province'])[:10]:<10s}) "
                        f"month {r['arrival_month']:>3d} "
                        f"{('[' + r['arrival_calendar'] + ']') if r.get('arrival_calendar') else ''}  "
                        f"final mean={r['final_mean_attack_pct']:5.1f}%  "
                        f"area reached={r['final_pct_area_reached']:5.1f}%  "
                        f"via {r['nearest_gate']} ({r['nearest_gate_km']:.0f} km)")
        if not reached:
            logger.warning("  ⚠ NO county reached the arrival threshold — PPR did not "
                           "establish in South Korea within this run. Check the "
                           "corridor crossing report above, and raise MAX_SIM_YEARS or "
                           "SPREAD_SPEED_MULTIPLIER if the front is simply too slow.")
    return rows, series
SK_SUMMARY_COLUMNS = [
    "rank_arrival", "sig_cd", "county_kr", "county_en", "province", "province_kr",
    "status", "arrival_month", "arrival_calendar", "arrival_year_of_sim",
    "arrival_month_first_cell", "months_after_first_sk_arrival",
    "final_mean_attack_pct", "final_max_attack_pct", "final_pct_area_reached",
    "peak_mean_infected_pct", "peak_infected_month",
    "nearest_gate", "nearest_gate_km", "nearest_gate_weight",
    "poly_area_km2", "model_cells", "model_area_km2", "sampling",
    "pct_cells_south_of_mdl", "centroid_lat", "centroid_lon",
]
def _fmt(v, nd=3):
    if v is None:
        return ""
    if isinstance(v, float):
        if not np.isfinite(v):
            return ""
        return f"{v:.{nd}f}"
    return str(v)
def write_county_summary_csv(rows, out_path, meta_lines=(), logger=None):
    """One row per county. UTF-8 with BOM so Excel on a Korean machine renders
    the Hangul names; `#` comment lines carry the provenance of the run."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8-sig") as fh:
        for ln in meta_lines:
            fh.write(f"# {ln}\n")
        wtr = csv.writer(fh)
        wtr.writerow(SK_SUMMARY_COLUMNS)
        for r in rows:
            wtr.writerow([_fmt(r.get(k)) for k in SK_SUMMARY_COLUMNS])
    if logger is not None:
        logger.info(f"   ✓ county summary CSV: {out_path}  ({len(rows)} rows)")
    return out_path
def write_county_monthly_wide_csv(rows, series, out_path, month_offset=0,
                                  start_year=None, start_month=1,
                                  value="mean", meta_lines=(), logger=None):
    """counties × months matrix of mean cumulative attack (%) — one column per
    simulation month, headed `mNNN` (or `mNNN_YYYY-MM` when a nominal start date
    is configured)."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    NM = 0
    for s in series.values():
        NM = max(NM, len(s[value]))
    hdr = ["sig_cd", "county_kr", "county_en", "province", "arrival_month"]
    for m in range(NM):
        g = m + 1 + int(month_offset)
        cal = _sim_month_to_calendar(g, start_year, start_month)
        hdr.append(f"m{g:03d}" + (f"_{cal}" if cal else ""))
    with open(out_path, "w", newline="", encoding="utf-8-sig") as fh:
        for ln in meta_lines:
            fh.write(f"# {ln}\n")
        wtr = csv.writer(fh)
        wtr.writerow(hdr)
        for r in rows:
            s = series.get(r["zid"], {}).get(value)
            vals = ([f"{100.0 * float(v):.3f}" for v in s] if s is not None else [""] * NM)
            wtr.writerow([r.get("sig_cd", ""), r.get("county_kr", ""),
                          r.get("county_en", ""), r.get("province", ""),
                          _fmt(r.get("arrival_month"))] + vals)
    if logger is not None:
        logger.info(f"   ✓ county monthly-wide CSV: {out_path}  "
                    f"({len(rows)} counties × {NM} months)")
    return out_path
def write_county_long_csv(rows, series, out_path, month_offset=0,
                          start_year=None, start_month=1,
                          meta_lines=(), logger=None):
    """Tidy/long form: one row per county-month. This is the file to load into
    pandas / R for time-series plots and per-county epidemic curves."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(out_path, "w", newline="", encoding="utf-8-sig") as fh:
        for ln in meta_lines:
            fh.write(f"# {ln}\n")
        wtr = csv.writer(fh)
        wtr.writerow(["sig_cd", "county_kr", "county_en", "province",
                      "month", "calendar", "mean_attack_pct",
                      "pct_area_reached", "max_attack_pct"])
        for r in rows:
            s = series.get(r["zid"])
            if not s:
                continue
            for m in range(len(s["mean"])):
                g = m + 1 + int(month_offset)
                wtr.writerow([r.get("sig_cd", ""), r.get("county_kr", ""),
                              r.get("county_en", ""), r.get("province", ""),
                              g, _sim_month_to_calendar(g, start_year, start_month),
                              f"{100.0 * float(s['mean'][m]):.3f}",
                              f"{100.0 * float(s['frac'][m]):.3f}",
                              f"{100.0 * float(s['mx'][m]):.3f}"])
                n += 1
    if logger is not None:
        logger.info(f"   ✓ county long CSV: {out_path}  ({n:,} county-month rows)")
    return out_path
def zones_to_value_raster(rows, zone_raster, key, scale=1.0):
    """Paint a per-county value back onto the model grid using the dominant-zone
    raster, so the CSV can be checked against a map."""
    max_zid = int(zone_raster.max())
    lut = np.full(max_zid + 1, np.nan, np.float64)
    for r in rows:
        z = int(r.get("zid", 0))
        v = r.get(key)
        if 0 < z <= max_zid and v is not None:
            try:
                fv = float(v)
            except (TypeError, ValueError):
                continue
            if np.isfinite(fv):
                lut[z] = fv * scale
    lut[0] = np.nan
    return lut[np.clip(zone_raster, 0, max_zid)]
# ══════════════════════════════════════════════════════════════════════════════
# INLINED SPATIAL MCDA ENGINE (v9.7/9.8) — criteria -> weights -> score -> THRESHOLDS
#
# Kept in-file, like the RAD engine below, so the script stays a single runnable
# unit. numpy only; no new dependency.
#
# WHAT IS OBJECTIVE HERE AND WHAT IS NOT — read this before quoting any number
# this section produces:
#   OBJECTIVE      standardisation of each criterion onto [0,1]; the aggregation
#                  arithmetic; the break-finding algorithms; every validation
#                  statistic computed against the epidemic.
#   JUDGEMENT      which criteria enter at all; their directions; the WEIGHTS;
#                  and WHICH DOMAIN is classified, since a class is defined
#                  relative to the area it was computed over.
#                  The consistency ratio tests whether pairwise judgements agree
#                  WITH EACH OTHER. It cannot tell you whether they are right.
#   NEITHER        Jenks and quantile breaks are descriptions of one map's
#                  histogram. On their own they carry NO external validity: they
#                  will manufacture a "very high risk" class out of a uniformly
#                  low-risk surface. That is why validate_against_epidemic()
#                  exists, and why SMCDA_BREAK_METHOD="absolute" is the option to
#                  move to once the score construction has settled.
# ══════════════════════════════════════════════════════════════════════════════
SM_EPS = 1e-12

# Saaty's Random Index (RI): the mean consistency index of randomly generated
# reciprocal matrices of order n. Used to turn CI into the consistency RATIO.
# Saaty TL (1980) The Analytic Hierarchy Process. McGraw-Hill.
SAATY_RI = {1: 0.00, 2: 0.00, 3: 0.58, 4: 0.90, 5: 1.12, 6: 1.24,
            7: 1.32, 8: 1.41, 9: 1.45, 10: 1.49, 11: 1.51, 12: 1.48,
            13: 1.56, 14: 1.57, 15: 1.59}


# ══════════════════════════════════════════════════════════════════════════════
# 1. CRITERION DEFINITION
# ══════════════════════════════════════════════════════════════════════════════
@dataclass
class SMCDACriterion:
    """One SMCDA criterion.

    key        : identifier used in the AHP matrix and all reports
    label      : human-readable name for maps/tables
    direction  : "benefit" (higher raw value = higher risk) or
                 "cost"    (higher raw value = LOWER risk; inverted internally)
    standardise: "minmax" | "rank" | "fuzzy_lin" | "fuzzy_sig" | "none"
                 - minmax    : (x - min) / (max - min) over the domain
                 - rank      : percentile rank; immune to outliers and to the
                               units of the layer, so criteria measured on wildly
                               different scales become comparable. Preferred when
                               a layer is heavily skewed (farm density is).
                 - fuzzy_lin : linear membership between control points a and b
                 - fuzzy_sig : sigmoidal membership between a and b
                 - none      : already on 0-1 and meaningful as-is
    a, b       : control points for the fuzzy transforms (raw units)
    veto_below : optional NON-COMPENSATORY floor on the STANDARDISED value.
                 Cells below it are forced to the lowest risk class regardless
                 of every other criterion. Use for constraints that cannot be
                 traded away (e.g. zero susceptible hosts).
    note       : why this criterion is in the model — carried into the report
    """
    key: str
    label: str
    direction: str = "benefit"
    standardise: str = "rank"
    a: Optional[float] = None
    b: Optional[float] = None
    veto_below: Optional[float] = None
    note: str = ""


def standardise_criterion(x: np.ndarray, mask: np.ndarray, c: SMCDACriterion) -> np.ndarray:
    """Map a raw criterion layer onto [0,1] where 1 = maximum contribution to
    risk. Cost criteria are inverted AFTER scaling so the direction flip cannot
    interact with the choice of transform."""
    v = np.asarray(x, dtype=np.float64)
    m = mask & np.isfinite(v)
    out = np.zeros(v.shape, dtype=np.float64)
    if not np.any(m):
        return out

    vals = v[m]
    method = (c.standardise or "rank").lower()

    if method == "none":
        s = np.clip(vals, 0.0, 1.0)

    elif method == "minmax":
        lo, hi = float(np.min(vals)), float(np.max(vals))
        s = np.zeros_like(vals) if hi - lo < SM_EPS else (vals - lo) / (hi - lo)

    elif method == "rank":
        # Percentile rank (average ranks for ties). Robust to skew and outliers.
        order = np.argsort(vals, kind="mergesort")
        ranks = np.empty(vals.size, dtype=np.float64)
        ranks[order] = np.arange(vals.size, dtype=np.float64)
        # average tied ranks so plateaus (e.g. many zero-farm cells) get one value
        uniq, inv, cnt = np.unique(vals, return_inverse=True, return_counts=True)
        csum = np.cumsum(cnt)
        start = csum - cnt
        mean_rank = (start + csum - 1) / 2.0
        ranks = mean_rank[inv]
        s = ranks / max(vals.size - 1, 1)

    elif method in ("fuzzy_lin", "fuzzy_sig"):
        a = float(c.a if c.a is not None else np.nanpercentile(vals, 5))
        b = float(c.b if c.b is not None else np.nanpercentile(vals, 95))
        if abs(b - a) < SM_EPS:
            s = np.zeros_like(vals)
        else:
            t = np.clip((vals - a) / (b - a), 0.0, 1.0)
            s = t if method == "fuzzy_lin" else (np.sin(0.5 * math.pi * t) ** 2)

    else:
        raise ValueError(f"unknown standardisation '{c.standardise}' for '{c.key}'")

    if c.direction.lower() == "cost":
        s = 1.0 - s
    elif c.direction.lower() != "benefit":
        raise ValueError(f"direction must be 'benefit' or 'cost', got '{c.direction}'")

    out[m] = np.clip(s, 0.0, 1.0)
    return out


# ══════════════════════════════════════════════════════════════════════════════
# 2. WEIGHTING
# ══════════════════════════════════════════════════════════════════════════════
def ahp_weights(matrix: np.ndarray, keys: Sequence[str],
                cr_limit: float = 0.10) -> Dict:
    """Analytic Hierarchy Process weights by the principal-eigenvector method,
    with Saaty's consistency ratio.

    Returns dict(weights, lambda_max, CI, CR, consistent, n, cr_limit).

    CR <= 0.10 is Saaty's conventional acceptance threshold. A CR ABOVE it means
    the pairwise judgements contradict each other (A>B, B>C, but C>A) and the
    resulting weights should not be used until the matrix is revised. This
    function never silently "fixes" an inconsistent matrix.
    """
    A = np.asarray(matrix, dtype=np.float64)
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        raise ValueError("AHP matrix must be square")
    if len(keys) != n:
        raise ValueError(f"{len(keys)} keys for a {n}x{n} matrix")
    if not np.allclose(np.diag(A), 1.0, atol=1e-9):
        raise ValueError("AHP matrix diagonal must be exactly 1")
    recip_err = float(np.max(np.abs(A * A.T - 1.0)))
    if recip_err > 1e-6:
        raise ValueError(f"AHP matrix is not reciprocal (max |a_ij*a_ji - 1| = {recip_err:.2e})")

    eigval, eigvec = np.linalg.eig(A)
    k = int(np.argmax(eigval.real))
    lam = float(eigval.real[k])
    w = np.abs(eigvec[:, k].real)
    w = w / w.sum()

    ci = (lam - n) / (n - 1) if n > 2 else 0.0
    ri = SAATY_RI.get(n, 1.59)
    cr = ci / ri if ri > 0 else 0.0
    return dict(weights={k_: float(w_) for k_, w_ in zip(keys, w)},
                weight_vector=w, lambda_max=lam, CI=float(ci), CR=float(cr),
                consistent=bool(cr <= float(cr_limit)), n=n,
                cr_limit=float(cr_limit))


def ahp_matrix_from_priorities(priorities: Sequence[float],
                               scale_max: float = 9.0) -> np.ndarray:
    """Build a reciprocal Saaty matrix from an ORDINAL PRIORITY VECTOR:
    a_ij = clip(round(s_i / s_j), 1/scale_max, scale_max), with a_ji = 1/a_ij.

    This is a CONSTRUCTION, not an elicitation. It is offered so a runnable
    default exists and so the judgement is stated in one readable place (the
    priority vector) instead of 21 pairwise cells. Because the matrix is derived
    from a single vector it is near-consistent BY CONSTRUCTION, so a low CR here
    confirms the arithmetic and says nothing about whether the priorities are
    epidemiologically right. Replace it with a real pairwise matrix before the
    output drives resource allocation."""
    s = np.asarray(priorities, dtype=np.float64)
    if s.ndim != 1 or s.size == 0:
        raise ValueError("priority vector must be a non-empty 1-D sequence")
    if np.any(s <= 0):
        raise ValueError("priority vector entries must be strictly positive")
    n = s.size
    A = np.ones((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            r = s[i] / s[j]
            v = float(np.clip(round(r) if r >= 1.0 else 1.0 / max(round(1.0 / r), 1.0),
                              1.0 / scale_max, scale_max))
            if v <= 0:
                v = 1.0
            A[i, j] = v
            A[j, i] = 1.0 / v
    return A


def rank_sum_weights(order: Sequence[str]) -> Dict[str, float]:
    """Rank-sum weights: w_i = (n - r_i + 1) / sum(n - r_j + 1), r = 1 is most
    important. A transparent, judgement-light fallback when no defensible
    pairwise matrix exists."""
    n = len(order)
    raw = np.array([n - i for i in range(n)], dtype=np.float64)
    raw = raw / raw.sum()
    return {k: float(v) for k, v in zip(order, raw)}


def entropy_weights(std_stack: Dict[str, np.ndarray], mask: np.ndarray) -> Dict[str, float]:
    """Shannon-entropy (objective) weights: criteria that DISCRIMINATE more
    between cells get more weight. Purely data-driven — no expert input — so it
    is the right cross-check on an AHP result, not a replacement for it: entropy
    knows nothing about epidemiology and will happily up-weight a noisy layer."""
    out = {}
    ent = {}
    m = mask
    n = int(m.sum())
    if n <= 1:
        return {k: 1.0 / max(len(std_stack), 1) for k in std_stack}
    for k, arr in std_stack.items():
        p = np.clip(arr[m], SM_EPS, None)
        tot = float(p.sum())
        if tot <= SM_EPS:
            ent[k] = 1.0
            continue
        p = p / tot
        e = -float(np.sum(p * np.log(p))) / math.log(n)
        ent[k] = min(max(e, 0.0), 1.0)
    d = {k: 1.0 - v for k, v in ent.items()}
    s = sum(d.values())
    if s <= SM_EPS:
        return {k: 1.0 / len(d) for k in d}
    for k, v in d.items():
        out[k] = float(v / s)
    return out


# ══════════════════════════════════════════════════════════════════════════════
# 3. AGGREGATION
# ══════════════════════════════════════════════════════════════════════════════
def wlc_aggregate(std_stack: Dict[str, np.ndarray], weights: Dict[str, float],
                  mask: np.ndarray) -> np.ndarray:
    """Weighted Linear Combination — the standard SMCDA aggregator.
    FULLY COMPENSATORY: a deficit on one criterion is offset by a surplus on
    another. Where that is not acceptable, use OWA or a veto."""
    shape = next(iter(std_stack.values())).shape
    acc = np.zeros(shape, dtype=np.float64)
    wsum = 0.0
    for k, arr in std_stack.items():
        w = float(weights.get(k, 0.0))
        if w <= 0:
            continue
        acc += w * arr
        wsum += w
    if wsum > SM_EPS:
        acc /= wsum
    return np.where(mask, np.clip(acc, 0.0, 1.0), np.nan)


def owa_aggregate(std_stack: Dict[str, np.ndarray], weights: Dict[str, float],
                  mask: np.ndarray, andness: float = 0.5) -> np.ndarray:
    """Ordered Weighted Averaging (Yager 1988) with a risk-attitude parameter.

    andness = 0.5 -> equivalent to WLC (neutral, full trade-off)
    andness -> 1.0 -> AND-like / risk-AVERSE: the score is pulled toward the
                      WORST criteria, so a cell must score well on nearly all of
                      them. Use this when criteria should NOT compensate.
    andness -> 0.0 -> OR-like / risk-TAKING: pulled toward the BEST criteria.

    Implemented by reordering the criterion weights per cell and applying order
    weights derived from a power of the cumulative criterion weight.

    DERIVATION OF alpha (do not change this without redoing it).
      SMCDACriterion values are sorted BEST FIRST.  With the RIM quantifier
      Q(r) = r^alpha applied to the cumulative criterion weight, the order
      weight on position j is Q(cum_j) - Q(cum_{j-1}).  Yager's orness of that
      weight vector is the integral of Q over [0,1]:

          orness = INT_0^1 r^alpha dr = 1 / (alpha + 1)
          andness = 1 - orness       = alpha / (alpha + 1)
      =>  alpha  = andness / (1 - andness)

      alpha = 1   -> equal order weights          -> WLC          (andness 0.5)
      alpha -> inf-> all weight on the LAST (worst) position -> MIN / AND
      alpha -> 0  -> all weight on the FIRST (best) position  -> MAX / OR
    """
    keys = list(std_stack.keys())
    if not keys:
        raise ValueError("empty criterion stack")
    a = float(np.clip(andness, 0.0, 1.0))
    # alpha maps andness -> the exponent on the cumulative-weight curve
    if a <= SM_EPS:
        alpha = 1.0 / 50.0          # OR-like limit: weight on the best criterion
    elif a >= 1.0 - SM_EPS:
        alpha = 50.0                # AND-like limit: weight on the worst criterion
    else:
        alpha = a / (1.0 - a)

    shape = std_stack[keys[0]].shape
    S = np.stack([std_stack[k] for k in keys], axis=0)                 # (K, R, C)
    W = np.array([float(weights.get(k, 0.0)) for k in keys], dtype=np.float64)
    if W.sum() <= SM_EPS:
        W = np.ones_like(W)
    W = W / W.sum()

    flat = S.reshape(len(keys), -1)                                    # (K, N)
    order = np.argsort(-flat, axis=0)                                  # best first
    sorted_vals = np.take_along_axis(flat, order, axis=0)
    sorted_w = W[order]                                                # (K, N)

    cum = np.cumsum(sorted_w, axis=0)
    prev = np.vstack([np.zeros((1, cum.shape[1])), cum[:-1]])
    ow = np.power(cum, alpha) - np.power(prev, alpha)                  # order weights
    denom = ow.sum(axis=0, keepdims=True)
    ow = np.where(denom > SM_EPS, ow / np.maximum(denom, SM_EPS), 1.0 / len(keys))

    agg = (ow * sorted_vals).sum(axis=0).reshape(shape)
    return np.where(mask, np.clip(agg, 0.0, 1.0), np.nan)


# ══════════════════════════════════════════════════════════════════════════════
# 4. CLASSIFICATION — THE THRESHOLDS
# ══════════════════════════════════════════════════════════════════════════════
def _fisher_jenks_breaks(values: np.ndarray, k: int) -> List[float]:
    """Fisher-Jenks natural breaks by dynamic programming, minimising the sum of
    within-class squared deviations.

    Exact Fisher-Jenks is O(n^2 k) and impossible on a 700k-cell raster, so the
    values are reduced to `nb` ordered representatives (quantile bins) first and
    the DP is run on those. This is what production implementations do; the
    breaks are stable to ~1 representative width, which is far finer than any
    interpretation the classes carry."""
    v = np.sort(np.asarray(values, dtype=np.float64))
    v = v[np.isfinite(v)]
    if v.size == 0:
        return []
    k = max(2, int(k))
    if v.size <= k:
        return [float(v[0]), float(v[-1])]

    nb = min(v.size, 512)
    qs = np.linspace(0.0, 1.0, nb)
    reps = np.quantile(v, qs)
    reps = np.unique(reps)
    n = reps.size
    if n <= k:
        return [float(reps[0]), float(reps[-1])]

    # prefix sums for O(1) within-class variance
    p1 = np.concatenate([[0.0], np.cumsum(reps)])
    p2 = np.concatenate([[0.0], np.cumsum(reps ** 2)])

    def sse(i: int, j: int) -> float:
        """SSE of reps[i:j] (0-based, j exclusive)."""
        m = j - i
        if m <= 1:
            return 0.0
        s1 = p1[j] - p1[i]
        s2 = p2[j] - p2[i]
        return max(s2 - s1 * s1 / m, 0.0)

    INF = float("inf")
    D = np.full((k + 1, n + 1), INF)
    B = np.zeros((k + 1, n + 1), dtype=np.int32)
    D[0, 0] = 0.0
    for c in range(1, k + 1):
        for j in range(c, n + 1):
            best, arg = INF, c - 1
            for i in range(c - 1, j):
                if D[c - 1, i] == INF:
                    continue
                val = D[c - 1, i] + sse(i, j)
                if val < best:
                    best, arg = val, i
            D[c, j] = best
            B[c, j] = arg

    cuts, j = [], n
    for c in range(k, 0, -1):
        i = int(B[c, j])
        cuts.append(i)
        j = i
    cuts = sorted(set(cuts))
    edges = [float(reps[0])]
    for i in cuts[1:]:
        edges.append(float(reps[min(i, n - 1)]))
    edges.append(float(reps[-1]))
    edges = sorted(set(edges))
    return edges


def compute_breaks(score: np.ndarray, mask: np.ndarray, k: int = 5,
                   method: str = "jenks",
                   absolute: Optional[Sequence[float]] = None) -> Dict:
    """Return the class EDGES (length k+1) plus a description of how they were
    obtained.

    method:
      "jenks"    natural breaks — minimises within-class variance. Data-driven:
                 finds the map's own structure but is NOT transferable between
                 maps or years; re-derive it whenever the inputs change.
      "quantile" equal COUNT per class. Guarantees class sizes but forces a
                 fixed share of the domain into "very high" even if nowhere is
                 genuinely high-risk. Good for ranking, wrong for absolute
                 statements.
      "equal"    equal score INTERVAL. Comparable across runs, but classes can
                 be empty or hold almost everything.
      "sd"       mean +/- multiples of the standard deviation. Meaningful when
                 the score is roughly symmetric; the SMCDA score usually is not.
      "absolute" caller-supplied cut-offs. THE ONLY METHOD THAT IS COMPARABLE
                 ACROSS RUNS, SCENARIOS AND YEARS — use it once the operational
                 thresholds have been agreed.
    """
    v = np.asarray(score, dtype=np.float64)[mask & np.isfinite(score)]
    k = max(2, int(k))
    if v.size == 0:
        return dict(edges=[0.0, 1.0], method=method, note="empty domain")

    m = (method or "jenks").lower()
    if m == "absolute":
        if not absolute:
            raise ValueError("method='absolute' requires the cut-offs")
        edges = [0.0] + [float(x) for x in absolute] + [1.0]
        edges = sorted(set(np.clip(edges, 0.0, 1.0).tolist()))
        note = "caller-supplied absolute cut-offs (comparable across runs)"
    elif m == "quantile":
        qs = np.linspace(0.0, 1.0, k + 1)
        edges = [float(x) for x in np.quantile(v, qs)]
        note = f"equal-count quantiles ({100.0 / k:.0f}% of cells per class)"
    elif m == "equal":
        lo, hi = float(v.min()), float(v.max())
        edges = [float(x) for x in np.linspace(lo, hi, k + 1)]
        note = "equal score intervals"
    elif m == "sd":
        mu, sd = float(v.mean()), float(v.std())
        half = (k - 1) / 2.0
        edges = [float(v.min())] + [mu + (i - half) * sd for i in range(k - 1)] + [float(v.max())]
        edges = sorted(set(np.clip(edges, float(v.min()), float(v.max())).tolist()))
        note = f"mean {mu:.3f} +/- k*SD ({sd:.3f})"
    elif m == "jenks":
        edges = _fisher_jenks_breaks(v, k)
        note = "Fisher-Jenks natural breaks (within-class variance minimised)"
    else:
        raise ValueError(f"unknown break method '{method}'")

    edges = sorted(set(float(e) for e in edges))
    while len(edges) < k + 1:                      # degenerate distributions
        edges.append(min(1.0, edges[-1] + 1e-6))
        edges = sorted(set(edges))
    if len(edges) > k + 1:
        idx = np.linspace(0, len(edges) - 1, k + 1).round().astype(int)
        edges = [edges[i] for i in idx]
    edges[0] = min(edges[0], float(v.min()))
    edges[-1] = max(edges[-1], float(v.max()))
    return dict(edges=edges, method=m, note=note)


def classify(score: np.ndarray, edges: Sequence[float], mask: np.ndarray) -> np.ndarray:
    """Assign each cell a 0-based class index from the edges. NaN outside mask."""
    e = np.asarray(edges, dtype=np.float64)
    out = np.full(score.shape, np.nan, dtype=np.float32)
    v = np.asarray(score, dtype=np.float64)
    idx = np.clip(np.digitize(v, e[1:-1], right=True), 0, len(e) - 2)
    out[mask & np.isfinite(v)] = idx[mask & np.isfinite(v)].astype(np.float32)
    return out


# ══════════════════════════════════════════════════════════════════════════════
# 5. VALIDATION — do the classes actually rank the epidemic?
# ══════════════════════════════════════════════════════════════════════════════
def _auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """AUC via the Mann-Whitney U identity (rank-based; ties handled)."""
    s = np.asarray(scores, float)
    y = np.asarray(labels, bool)
    npos, nneg = int(y.sum()), int((~y).sum())
    if npos == 0 or nneg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(s.size, float)
    ranks[order] = np.arange(1, s.size + 1, dtype=float)
    uniq, inv, cnt = np.unique(s, return_inverse=True, return_counts=True)
    csum = np.cumsum(cnt)
    start = csum - cnt
    mean_rank = (start + csum + 1) / 2.0
    ranks = mean_rank[inv]
    return float((ranks[y].sum() - npos * (npos + 1) / 2.0) / (npos * nneg))


def validate_against_epidemic(score: np.ndarray, cls: np.ndarray, mask: np.ndarray,
                              time_to_arrival: Optional[np.ndarray] = None,
                              cum_attack: Optional[np.ndarray] = None,
                              arrival_threshold: float = 0.05,
                              class_labels: Optional[Sequence[str]] = None) -> Dict:
    """Test the SMCDA classification against the mechanistic epidemic.

    This is the step that gives the thresholds meaning. Jenks breaks describe one
    map's histogram; they say nothing about disease. The questions asked here:

      1. Is mean arrival month MONOTONIC across the classes (higher risk class
         reached earlier)? If not, the classification does not rank risk.
      2. AUC of the continuous score against "reached by the end of the run" —
         discrimination independent of where the breaks fall.
      3. Spearman rho between score and arrival month (expect NEGATIVE: higher
         score -> earlier arrival).
      4. Per-class reached fraction and mean final attack rate.

    A classification can be perfectly self-consistent and still fail every one of
    these. Report them next to the map."""
    out: Dict = {"n_domain": int(mask.sum())}
    m = mask & np.isfinite(cls)
    if not np.any(m):
        return out

    ncls = int(np.nanmax(cls[m])) + 1
    labels = list(class_labels) if class_labels else [f"class {i}" for i in range(ncls)]

    reached = None
    if cum_attack is not None:
        reached = np.nan_to_num(np.asarray(cum_attack, float)) >= float(arrival_threshold)

    rows = []
    for i in range(ncls):
        sel = m & (cls == i)
        n = int(sel.sum())
        row = dict(cls=i, label=labels[i] if i < len(labels) else f"class {i}",
                   cells=n, pct=100.0 * n / max(int(m.sum()), 1))
        if n:
            row["score_min"] = float(np.nanmin(score[sel]))
            row["score_max"] = float(np.nanmax(score[sel]))
            row["score_mean"] = float(np.nanmean(score[sel]))
        if time_to_arrival is not None and n:
            tta = np.asarray(time_to_arrival, float)[sel]
            fin = tta[np.isfinite(tta)]
            row["arr_mean"] = float(fin.mean()) + 1 if fin.size else None
            row["arr_median"] = float(np.median(fin)) + 1 if fin.size else None
            row["reached_frac"] = float(fin.size) / n
        if reached is not None and n:
            row["reached_frac"] = float(reached[sel].mean())
        if cum_attack is not None and n:
            row["attack_mean"] = float(np.nanmean(np.nan_to_num(cum_attack)[sel]))
        rows.append(row)
    out["classes"] = rows

    # (1) monotonicity of arrival month across classes
    arr = [r.get("arr_mean") for r in rows if r.get("arr_mean") is not None]
    if len(arr) >= 3:
        diffs = np.diff(np.asarray(arr, float))
        out["arrival_monotonic"] = bool(np.all(diffs <= 1e-9))
        out["arrival_by_class"] = arr
        out["arrival_inversions"] = int((diffs > 1e-9).sum())

    # (2) AUC on the continuous score
    if reached is not None:
        s = np.asarray(score, float)[m]
        y = reached[m]
        out["auc_reached"] = _auc(s, y)
        out["prevalence_reached"] = float(y.mean())

    # (3) Spearman rho, score vs arrival month
    if time_to_arrival is not None:
        tta = np.asarray(time_to_arrival, float)
        sel = m & np.isfinite(tta)
        if int(sel.sum()) > 10:
            a = score[sel]
            b = tta[sel]
            ra = np.argsort(np.argsort(a)).astype(float)
            rb = np.argsort(np.argsort(b)).astype(float)
            ra -= ra.mean(); rb -= rb.mean()
            den = math.sqrt(float((ra ** 2).sum()) * float((rb ** 2).sum()))
            out["spearman_score_vs_arrival"] = float((ra * rb).sum() / den) if den > 0 else float("nan")
            out["n_arrived"] = int(sel.sum())
    return out


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman rank correlation of two 1-D arrays (ties broken by ordering;
    adequate here because the inputs are continuous rasters)."""
    if a.size < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    den = math.sqrt(float((ra ** 2).sum()) * float((rb ** 2).sum()))
    return float((ra * rb).sum() / den) if den > 0 else float("nan")


def criterion_arrival_diagnostic(std_stack: Dict[str, np.ndarray],
                                 weights: Dict[str, float],
                                 mask: np.ndarray,
                                 time_to_arrival: np.ndarray,
                                 criteria: Optional[Sequence[SMCDACriterion]] = None) -> List[Dict]:
    """WHICH CRITERION IS BREAKING THE RANKING?

    When the aggregate classification fails its monotonicity test, the useful
    question is not "is the map wrong" but "which criterion disagrees with the
    epidemic". This correlates EACH standardised criterion, on its own, against
    arrival month.

    Expected sign is NEGATIVE: a criterion that genuinely raises risk should be
    high where the epidemic arrives EARLY. A criterion with a positive rho is
    pointing the wrong way over this domain — and if it also carries a large
    weight, it is the one dragging the aggregate score out of order.

    The commonest cause is a distance-to-entry criterion measured from the wrong
    origin: distance from a border crossing is the right predictor for the
    region the virus enters THROUGH it, and the wrong one for a region where the
    virus arrives by some other route, or arrives there last."""
    tta = np.asarray(time_to_arrival, dtype=np.float64)
    sel = mask & np.isfinite(tta)
    out: List[Dict] = []
    if int(sel.sum()) < 10:
        return out
    notes = {c.key: c.label for c in (criteria or [])}
    b = tta[sel]
    for k, arr in std_stack.items():
        a = np.asarray(arr, dtype=np.float64)[sel]
        rho = _spearman(a, b)
        w = float(weights.get(k, 0.0))
        out.append(dict(criterion=k, label=notes.get(k, k), weight=w,
                        rho_vs_arrival=rho,
                        # how much this criterion can pull the aggregate the
                        # wrong way: weight x how wrong-signed it is
                        misalignment=(w * rho if np.isfinite(rho) and rho > 0 else 0.0),
                        sign_ok=bool(np.isfinite(rho) and rho < 0)))
    out.sort(key=lambda d: -d["misalignment"])
    return out


# ══════════════════════════════════════════════════════════════════════════════
# 6. WEIGHT SENSITIVITY
# ══════════════════════════════════════════════════════════════════════════════
def weight_sensitivity(std_stack: Dict[str, np.ndarray], weights: Dict[str, float],
                       mask: np.ndarray, edges: Sequence[float],
                       delta: float = 0.20) -> List[Dict]:
    """One-at-a-time weight perturbation: raise each criterion's weight by
    `delta` (relative), renormalise the rest, and measure how much of the map
    changes class.

    A classification that reshuffles under a 20% weight nudge is not a robust
    basis for allocating surveillance. Report the largest class-change share
    alongside the map."""
    base = wlc_aggregate(std_stack, weights, mask)
    base_cls = classify(base, edges, mask)
    n = max(int(mask.sum()), 1)
    rows = []
    for k in std_stack:
        w2 = dict(weights)
        w2[k] = weights.get(k, 0.0) * (1.0 + float(delta))
        tot = sum(w2.values())
        if tot <= SM_EPS:
            continue
        w2 = {kk: vv / tot for kk, vv in w2.items()}
        s2 = wlc_aggregate(std_stack, w2, mask)
        c2 = classify(s2, edges, mask)
        changed = mask & np.isfinite(base_cls) & np.isfinite(c2) & (base_cls != c2)
        rows.append(dict(criterion=k,
                         weight=float(weights.get(k, 0.0)),
                         weight_perturbed=float(w2[k]),
                         pct_cells_changed_class=100.0 * float(changed.sum()) / n,
                         mean_abs_score_shift=float(np.nanmean(np.abs(
                             np.nan_to_num(s2) - np.nan_to_num(base))[mask]))))
    rows.sort(key=lambda r: -r["pct_cells_changed_class"])
    return rows


# ══════════════════════════════════════════════════════════════════════════════
# 7. TOP-LEVEL DRIVER
# ══════════════════════════════════════════════════════════════════════════════
def run_smcda(layers: Dict[str, np.ndarray],
              criteria: Sequence[SMCDACriterion],
              mask: np.ndarray,
              ahp_matrix: Optional[np.ndarray] = None,
              ahp_cr_limit: float = 0.10,
              weight_method: str = "ahp",
              aggregation: str = "wlc",
              owa_andness: float = 0.5,
              n_classes: int = 5,
              break_method: str = "jenks",
              absolute_breaks: Optional[Sequence[float]] = None,
              class_labels: Optional[Sequence[str]] = None,
              time_to_arrival: Optional[np.ndarray] = None,
              cum_attack: Optional[np.ndarray] = None,
              arrival_threshold: float = 0.05,
              sensitivity_delta: float = 0.20,
              logger=None) -> Dict:
    """Full SMCDA chain. Returns everything needed to report and to map."""
    def _log(msg, warn=False):
        if logger is None:
            print(msg)
        elif warn:
            logger.warning(msg)
        else:
            logger.info(msg)

    crit = [c for c in criteria if c.key in layers]
    missing = [c.key for c in criteria if c.key not in layers]
    if missing:
        _log(f"  ⚠ SMCDA: criteria with no layer supplied — DROPPED: {missing}. "
             f"The weights of the remaining criteria are renormalised, so the "
             f"result is NOT the model you specified. Supply the layers or edit "
             f"the criterion list.", warn=True)
    if not crit:
        raise ValueError("SMCDA: no criteria have a matching layer")

    keys = [c.key for c in crit]

    # ── standardise ──────────────────────────────────────────────────────────
    std = {}
    for c in crit:
        std[c.key] = standardise_criterion(layers[c.key], mask, c)
    _log(f"  SMCDA: {len(crit)} criteria standardised onto [0,1] "
         f"(1 = maximum contribution to risk).")
    for c in crit:
        v = std[c.key][mask]
        _log(f"    · {c.key:<10s} {c.label:<34s} {c.direction:<7s} "
             f"{c.standardise:<9s} mean={float(v.mean()):.3f} "
             f"p90={float(np.percentile(v, 90)):.3f}")

    # ── weights ──────────────────────────────────────────────────────────────
    ahp_info = None
    wm = (weight_method or "ahp").lower()
    if wm == "ahp":
        if ahp_matrix is None:
            raise ValueError("weight_method='ahp' requires ahp_matrix")
        A = np.asarray(ahp_matrix, float)
        if A.shape[0] != len(criteria):
            raise ValueError(f"AHP matrix is {A.shape[0]}x{A.shape[0]} but "
                             f"{len(criteria)} criteria were defined")
        if missing:                                  # drop rows/cols of dropped criteria
            keep = [i for i, c in enumerate(criteria) if c.key in layers]
            A = A[np.ix_(keep, keep)]
        ahp_info = ahp_weights(A, keys, cr_limit=ahp_cr_limit)
        weights = ahp_info["weights"]
        _log(f"  SMCDA weights: AHP eigenvector, lambda_max={ahp_info['lambda_max']:.4f}, "
             f"CI={ahp_info['CI']:.4f}, CR={ahp_info['CR']:.4f} "
             f"({'CONSISTENT' if ahp_info['consistent'] else 'INCONSISTENT'})")
        if not ahp_info["consistent"]:
            _log(f"  ⚠ AHP consistency ratio {ahp_info['CR']:.3f} EXCEEDS the "
                 f"{ahp_info['cr_limit']:.2f} acceptance limit (Saaty's convention is "
                 f"0.10). The pairwise judgements contradict each other, so these "
                 f"weights are not defensible. Revise the matrix before using this "
                 f"classification. (CR tests self-consistency only — it cannot tell "
                 f"you whether the judgements are epidemiologically correct.)", warn=True)
    elif wm == "rank_sum":
        weights = rank_sum_weights(keys)
        _log("  SMCDA weights: rank-sum (order of the criterion list).")
    elif wm == "entropy":
        weights = entropy_weights(std, mask)
        _log("  SMCDA weights: Shannon entropy (objective, data-driven). NOTE: entropy "
             "rewards layers that discriminate between cells; it has no knowledge of "
             "epidemiology and will up-weight a noisy layer.")
    elif wm == "equal":
        weights = {k: 1.0 / len(keys) for k in keys}
        _log("  SMCDA weights: equal.")
    else:
        raise ValueError(f"unknown weight_method '{weight_method}'")

    for k in sorted(weights, key=lambda x: -weights[x]):
        _log(f"    w[{k:<10s}] = {weights[k]:.4f}")

    # cross-check AHP against the objective entropy weights
    if wm == "ahp":
        ew = entropy_weights(std, mask)
        drift = {k: ew.get(k, 0) - weights.get(k, 0) for k in keys}
        worst = max(drift, key=lambda k: abs(drift[k]))
        _log(f"    [cross-check] entropy weights differ most on '{worst}' "
             f"(AHP {weights[worst]:.3f} vs entropy {ew[worst]:.3f}). Entropy is "
             f"data-driven and epidemiologically blind — a large gap is a prompt to "
             f"re-examine the judgement, not evidence that either is wrong.")

    # ── aggregate ────────────────────────────────────────────────────────────
    ag = (aggregation or "wlc").lower()
    if ag == "wlc":
        score = wlc_aggregate(std, weights, mask)
        _log("  SMCDA aggregation: WLC (fully compensatory).")
    elif ag == "owa":
        score = owa_aggregate(std, weights, mask, andness=owa_andness)
        _log(f"  SMCDA aggregation: OWA, andness={owa_andness:.2f} "
             f"({'risk-averse / AND-like' if owa_andness > 0.55 else 'risk-taking / OR-like' if owa_andness < 0.45 else 'neutral, ~WLC'}).")
    else:
        raise ValueError(f"unknown aggregation '{aggregation}'")

    # ── non-compensatory vetoes ──────────────────────────────────────────────
    veto_total = 0
    for c in crit:
        if c.veto_below is None:
            continue
        v = mask & (std[c.key] < float(c.veto_below))
        nv = int(v.sum())
        if nv:
            score = np.where(v, 0.0, score)
            veto_total += nv
            _log(f"  SMCDA veto: {nv:,} cells forced to score 0 because "
                 f"'{c.key}' < {c.veto_below:.2f} (non-compensatory constraint).")
    if veto_total:
        _log(f"  SMCDA: {veto_total:,} cells vetoed in total ({100.0 * veto_total / max(int(mask.sum()),1):.1f}% "
             f"of the domain). These cannot be compensated by any other criterion.")

    # ── classify ─────────────────────────────────────────────────────────────
    br = compute_breaks(score, mask, k=n_classes, method=break_method,
                        absolute=absolute_breaks)
    edges = br["edges"]
    cls = classify(score, edges, mask)
    labels = list(class_labels) if class_labels else \
        ["Very low", "Low", "Moderate", "High", "Very high"][:n_classes]
    _log(f"  SMCDA thresholds ({br['method']}): "
         + " | ".join(f"{labels[i]} [{edges[i]:.4f},{edges[i+1]:.4f}]"
                      for i in range(min(len(labels), len(edges) - 1))))
    # An EMPTY class is not a harmless cosmetic detail. With absolute cut-offs it
    # means the score never reaches that level anywhere in the domain — which is
    # a substantive finding ("nothing here is high risk") that is easy to miss
    # when the map simply shows fewer colours than the legend.
    _empty = [labels[i] for i in range(len(labels))
              if int((mask & (cls == i)).sum()) == 0]
    if _empty:
        _log(f"  ⚠ SMCDA: {len(_empty)} class(es) contain NO cells: {_empty}. "
             f"The score range over this domain is "
             f"[{float(np.nanmin(score[mask])):.4f}, {float(np.nanmax(score[mask])):.4f}] "
             f"while the class edges run to {edges[-1]:.4f}. With break_method="
             f"'{br['method']}' this is a real statement about the domain, not a bug — "
             f"but the map will show fewer levels than the legend lists, so say so "
             f"wherever the classification is reported.", warn=True)

    # ── validate ─────────────────────────────────────────────────────────────
    val = validate_against_epidemic(score, cls, mask,
                                    time_to_arrival=time_to_arrival,
                                    cum_attack=cum_attack,
                                    arrival_threshold=arrival_threshold,
                                    class_labels=labels)

    # ── sensitivity ──────────────────────────────────────────────────────────
    sens = weight_sensitivity(std, weights, mask, edges, delta=sensitivity_delta) \
        if ag == "wlc" else []

    diag = (criterion_arrival_diagnostic(std, weights, mask, time_to_arrival, crit)
            if time_to_arrival is not None else [])

    return dict(score=score.astype(np.float32), classes=cls, edges=edges,
                break_method=br["method"], break_note=br["note"],
                weights=weights, ahp=ahp_info, standardised=std,
                criteria=crit, labels=labels, validation=val,
                criterion_diagnostic=diag,
                sensitivity=sens, n_classes=n_classes)


def format_smcda_threshold_table(res: Dict, cell_km: float = 1.0) -> List[str]:
    """The headline deliverable: the class thresholds as ACTUAL VALUES, with the
    epidemiological evidence for each class beside it."""
    e = res["edges"]
    val = res.get("validation", {})
    rows = {r["cls"]: r for r in val.get("classes", [])}
    # Distinguish "no epidemic was supplied" (the table is printed before the
    # simulation runs) from "the epidemic ran and did not reach this class".
    # Printing 'not reached' for the first is a false statement about the model.
    has_epi = any(("arr_mean" in r) or ("reached_frac" in r) for r in rows.values())
    out = []
    out.append("")
    out.append("  ┌─ SMCDA RISK-AREA THRESHOLDS " + "─" * 78)
    out.append(f"  │ Break method : {res['break_method']} — {res['break_note']}")
    out.append(f"  │ Classes      : {res['n_classes']}")
    out.append("  │")
    out.append("  │  Class          Score range        Cells        % domain   "
               "Reached   Mean arrival")
    out.append("  │  " + "-" * 82)
    for i, lab in enumerate(res["labels"]):
        lo, hi = e[i], e[i + 1]
        r = rows.get(i, {})
        n = r.get("cells", 0)
        pct = r.get("pct", 0.0)
        rf = r.get("reached_frac")
        am = r.get("arr_mean")
        area = n * cell_km * cell_km
        if rf is None:
            rf_s = "    —" if has_epi else "  n/a"
        else:
            rf_s = f"{100 * rf:5.1f}%"
        if am is not None:
            am_s = f"month {am:6.1f}"
        else:
            am_s = "  not reached" if has_epi else "  not yet run"
        out.append(f"  │  {lab:<13s} {lo:6.4f} – {hi:6.4f}  {n:9,}   {pct:6.1f}%   "
                   f"{rf_s}   {am_s}   ({area:,.0f} km²)")
    out.append("  │")
    if not has_epi:
        out.append("  │  (no epidemic supplied yet — the 'Reached' and 'Mean arrival' columns are")
        out.append("  │   filled in by the validation pass after the simulation completes. Nothing")
        out.append("  │   here has been tested against the modelled spread so far.)")

    mono = val.get("arrival_monotonic")
    if mono is not None:
        if mono:
            out.append("  │  ✓ VALIDATION: mean arrival month decreases monotonically from the "
                       "lowest to the")
            out.append("  │    highest risk class — the classification ranks the modelled epidemic "
                       "correctly.")
            # A PASSING aggregate can still contain a criterion that disagrees
            # with the epidemic; the other criteria simply outvote it. Worth
            # saying, because it is a live question about the weighting rather
            # than a failure of the map.
            _wrong = [d for d in (res.get("criterion_diagnostic") or [])
                      if not d["sign_ok"] and d["weight"] > 0.05]
            if _wrong:
                out.append("  │    Note: the aggregate passes, but "
                           + ", ".join(f"{d['criterion']} (w={d['weight']:.2f}, "
                                       f"rho={d['rho_vs_arrival']:+.2f})"
                                       for d in _wrong[:3]))
                out.append("  │    correlate the WRONG WAY with arrival on their own — the other "
                           "criteria outvote")
                out.append("  │    them. Usually collinearity (a criterion that rises with "
                           "distance from the")
                out.append("  │    entry will trail the epidemic whatever it measures). Worth "
                           "checking before")
                out.append("  │    the weights are treated as settled.")
        else:
            inv = val.get("arrival_inversions", 0)
            arr_by = val.get("arrival_by_class") or []
            out.append(f"  │  ⚠ VALIDATION FAILED: arrival month is NOT monotonic across classes "
                       f"({inv} inversion(s)).")
            if arr_by:
                out.append("  │    mean arrival by class (low -> high risk): "
                           + " -> ".join(f"{a:.0f}" for a in arr_by))
            out.append("  │    A higher-risk class is reached LATER than a lower one, so these "
                       "thresholds do")
            out.append("  │    NOT rank risk over this domain. Do not allocate surveillance from "
                       "them yet.")
            diag = res.get("criterion_diagnostic") or []
            wrong = [d for d in diag if not d["sign_ok"]]
            if diag:
                out.append("  │")
                out.append("  │    WHICH CRITERION DISAGREES WITH THE EPIDEMIC "
                           "(rho of each criterion vs arrival month;")
                out.append("  │    NEGATIVE is correct — high risk should mean early arrival):")
                for d in diag:
                    flag = "" if d["sign_ok"] else "   <-- WRONG SIGN"
                    out.append(f"  │      {d['criterion']:<8s} w={d['weight']:.3f}  "
                               f"rho={d['rho_vs_arrival']:+.3f}{flag}")
            if wrong:
                worst = wrong[0]
                out.append("  │")
                out.append(f"  │    '{worst['criterion']}' carries weight {worst['weight']:.2f} and "
                           f"points the WRONG WAY over this")
                out.append("  │    domain. If it is a distance-to-entry criterion, it is almost "
                           "certainly measured")
                out.append("  │    from an entry point the virus does not actually use to reach "
                           "this domain — or")
                out.append("  │    from one it reaches LAST. Classify the domain the introduction "
                           "route belongs to,")
                out.append("  │    or re-point the criterion at the route that applies here. "
                           "Dropping the criterion")
                out.append("  │    or shrinking its weight would hide the mismatch rather than "
                           "fix it.")
    auc = val.get("auc_reached")
    prev = val.get("prevalence_reached")
    if auc is not None and np.isfinite(auc):
        q = ("excellent" if auc >= 0.90 else "good" if auc >= 0.80 else
             "fair" if auc >= 0.70 else "poor — barely better than chance"
             if auc >= 0.60 else "NO better than chance")
        out.append(f"  │  AUC (score vs 'reached by end of run') = {auc:.3f}  [{q}]"
                   f"   prevalence={100 * (prev or 0):.1f}%")
    elif auc is not None and prev is not None and (prev >= 0.999 or prev <= 0.001):
        # Not a failure and not a missing number: with only one outcome class
        # present, "discriminate reached from unreached" has nothing to
        # discriminate. Saying so beats printing nan or omitting the line.
        which = "every" if prev >= 0.999 else "no"
        out.append(f"  │  AUC: NOT APPLICABLE — {which} cell in this domain is 'reached' "
                   f"(prevalence {100 * prev:.1f}%),")
        out.append("  │    so there is no second class to discriminate against. This is "
                   "expected when the")
        out.append("  │    domain IS the invaded footprint. The real tests here are the "
                   "arrival-order")
        out.append("  │    monotonicity above and the Spearman rho below — and they are "
                   "the harder ones,")
        out.append("  │    because the score can no longer earn credit for telling "
                   "reached from unreached.")
    rho = val.get("spearman_score_vs_arrival")
    if rho is not None and np.isfinite(rho):
        out.append(f"  │  Spearman rho (score vs arrival month) = {rho:+.3f}  "
                   f"(negative is correct: higher score -> earlier arrival)"
                   f"{'' if rho < -0.1 else '   ⚠ WEAK OR WRONG-SIGNED'}")
    sens = res.get("sensitivity") or []
    if sens:
        top = sens[0]
        out.append(f"  │  Weight sensitivity (+20% one-at-a-time): most influential criterion is "
                   f"'{top['criterion']}'")
        out.append(f"  │    -> {top['pct_cells_changed_class']:.1f}% of cells change class. "
                   f"{'Classification is ROBUST.' if top['pct_cells_changed_class'] < 5 else 'Classification is SENSITIVE — treat class boundaries as soft.' if top['pct_cells_changed_class'] < 15 else '⚠ Classification is UNSTABLE under small weight changes.'}")
    out.append("  └" + "─" * 100)
    out.append("")
    return out


# ══════════════════════════════════════════════════════════════════════════════
# SELF-TEST
# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# INLINED RAD ENGINE (v9.0/9.1) — movement physics + disease dynamics + solver
# Reaction-Advection-Diffusion epidemic engine. Kept in-file (not a separate
# module) so the script runs as a single paste-able cell in notebooks. Uses the
# module-level `logger` and scipy already imported above. Modular by section and
# UDE-ready: `mechanistic_closure` is the single swappable ecological closure.
# ══════════════════════════════════════════════════════════════════════════════
RAD_EPS = 1e-9
_HAVE_SCIPY = True  # scipy.ndimage.distance_transform_edt/gaussian_filter imported above
@dataclass
class RADParams:
    """All tunable knobs for the RAD engine. Ecological defaults are sensible
    starting points, NOT fitted values — calibrate against occurrence data."""
    # Seasonal drivers (length-12 arrays), supplied by the host model.
    beta_monthly: np.ndarray                 # transmission scaler per month
    disp_km_monthly: np.ndarray              # gorani dispersal radius (km) per month
    tvir_monthly: np.ndarray                 # viral-survival scaler per month
    cell_km: float = 5.0                     # grid cell size (km)
    # Diffusion (random movement). Units: cell^2 / month at peak season.
    diff_max: float = 1.10
    diff_min: float = 0.02
    # Advection (directed movement up the habitat-attractiveness gradient).
    adv_max: float = 0.60                    # max speed, cells / month (CFL cap)
    adv_gain: float = 2.5                    # gradient-to-velocity gain
    # Reaction / recovery.
    gamma: float = 0.15                      # recovery+removal rate (per month)
    trans_hab_w: float = 0.6                 # habitat weight in transmissibility
    trans_lstk_w: float = 0.4                # livestock weight in transmissibility
    trans_base: float = 0.05                 # transmissibility floor
    # Border sources.
    seed_intensity: float = 0.60             # initial infected fraction at entries
    source_continuous: bool = True           # keep a low influx at entries each month
    source_rate: float = 0.02                # continuous influx (fraction / month)
    # Landscape barriers (Han et al. 2022 posterior medians as pass-through).
    mdl_fence_factor: float = 0.531          # φ: UNIFORM fallback if no gate field
    apply_mdl_barrier: bool = True
    road_factor: float = 0.362               # δ: fraction crossing major roads
    apply_road_barrier: bool = True
    # ── v9.1: corridor-gated MDL fence ──────────────────────────────────────
    # When env['MDL_FENCE_FIELD'] is supplied it REPLACES mdl_fence_factor with
    # a per-cell pass-through (sealed line + gates at the water-deer corridors).
    # Seasonal gating additionally scales the gates by the monthly dispersal
    # radius, so crossings concentrate in the dispersal peak.
    corridor_seasonal_gating: bool = True
    corridor_season_floor: float = 0.35      # winter gate opening, as a fraction of peak
    # Numerics.
    cfl_safety: float = 0.40
    max_substeps: int = 80
    i_cap: float = 6.0                       # hard cap on infected density (stability)
    # Diagnostics.
    arrival_threshold: float = 0.05          # attack-rate level defining "arrived"
    # ── Vaccination (intervention) ──────────────────────────────────────────
    # Moves susceptibles S -> V (immune). V never becomes infected, so the
    # cumulative attack rate is I+R (ever-infected), independent of V.
    vacc_enable: bool = False
    vacc_mode: str = "reactive_ring"         # "preemptive" | "continuous" | "reactive_ring"
    vacc_start_month: int = 6                # month index (0-based) vaccination begins
    vacc_coverage: float = 0.60              # max fraction of a cell that can be vaccinated
    vacc_rate: float = 0.30                  # per-month fraction of S vaccinated in targeted cells
    vacc_efficacy: float = 0.85              # vaccine efficacy (prob. of protection)
    vacc_ring_cells: int = 7                 # reactive-ring radius (grid cells) around infected
    vacc_inf_thr: float = 0.02               # infected density that triggers a reactive ring
    # ── v9.9 FIX A: calendar month of leg month 0 (0 = January). The Korea leg starts at the
    # NK-introduction month, so its seasonal drivers must start there, not in January.
    season_offset: int = 0
    # ── v9.9 E: optional campaign trigger and per-cell coverage cap (scenario V2).
    vacc_trigger_mask: Optional[np.ndarray] = None   # start once any cell here has A >= arrival_threshold
    vacc_cap_field: Optional[np.ndarray] = None      # per-cell multiplier (0..1) on vacc_coverage
    # ── v9.10: batch-run switches and the no-seasonality scenario ──────────────
    log_every: int = 1                               # 0 = silent (batch runs)
    store_states: bool = True                        # False = keep only the cumulative-attack stack
    disp_scale_ref: Optional[float] = None           # dispersal normaliser (None = max of the radii)
# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════
def _norm01(a: np.ndarray, mask: np.ndarray) -> np.ndarray:
    v = a[mask & np.isfinite(a)]
    if v.size == 0 or v.max() <= v.min():
        return np.where(mask, 0.0, 0.0).astype(np.float64)
    lo, hi = float(v.min()), float(v.max())
    return np.where(mask, np.clip((a - lo) / (hi - lo + RAD_EPS), 0.0, 1.0), 0.0).astype(np.float64)
def _get(env: Dict[str, np.ndarray], key: str, shape, default: float = 0.0) -> np.ndarray:
    a = env.get(key)
    if a is None:
        return np.full(shape, default, dtype=np.float64)
    return np.nan_to_num(a.astype(np.float64), nan=default)
# ══════════════════════════════════════════════════════════════════════════════
# MOVEMENT PHYSICS  (the UDE-swappable closure lives here)
# ══════════════════════════════════════════════════════════════════════════════
def mechanistic_closure(env: Dict[str, np.ndarray], p: RADParams):
    """Ecological closure: map static GIS layers -> the *shape* of the movement
    field (mobility surface M, attractiveness potential Φ) and the local
    transmissibility surface. This is the ONE function a UDE would replace:
    keep the signature (env, params) -> (M, Phi, transmissibility) and a learned
    net can supply the same three surfaces.
    Returns
    -------
    M    : (NR,NC) in [0,1] — mobility (scales BOTH diffusion & advection)
    Phi  : (NR,NC) in [0,1] — habitat attractiveness (advection climbs its gradient)
    beta_field : (NR,NC) in [0,1] — spatial transmissibility multiplier
    """
    mask = env["KP"].astype(bool)
    shape = mask.shape
    gor = _get(env, "GOR_EFF", shape)
    if not np.any(gor > 0):
        gor = _get(env, "GOR_nn", shape)
    cond = _get(env, "CONDUCTANCE", shape)
    resd = _get(env, "RESISTANCE", shape, default=1.0)
    watr = _get(env, "WATER_ATTRACT", shape)
    lstk = _get(env, "livestock_nn", shape)
    dens = _get(env, "GORANI_DENSITY", shape)
    # Mobility: animals move easily through suitable, well-connected, low-resistance
    # terrain; poorly through resistant/steep/urban cells.
    M = _norm01(0.45 * gor + 0.35 * cond + 0.20 * (1.0 - _norm01(resd, mask)), mask)
    # Attractiveness: animals are drawn to suitable habitat, connectivity and water
    # (gorani are water-associated — water ATTRACTS, unlike wild boar).
    Phi = _norm01(0.45 * gor + 0.20 * dens + 0.20 * watr + 0.15 * cond, mask)
    # Local transmissibility: habitat density of gorani + spillover with livestock.
    beta_field = np.clip(
        p.trans_hab_w * gor + p.trans_lstk_w * lstk + p.trans_base, 0.0, 1.0
    )
    beta_field = np.where(mask, beta_field, 0.0)
    return M, Phi, beta_field
def build_barrier(env: Dict[str, np.ndarray], p: RADParams, month: int = 0) -> np.ndarray:
    """Multiplicative transport barrier in (0,1]: reduces movement across the MDL
    fence and major roads (δ), grounded in Han et al. (2022).
    v9.1 — the MDL term is now a PER-CELL field when env['MDL_FENCE_FIELD'] is
    present: the line is sealed except at the water-deer corridor gates. Falls
    back to the uniform Han φ when no gate field is supplied, so the engine is
    backwards compatible."""
    mask = env["KP"].astype(bool)
    bar = np.ones(mask.shape, dtype=np.float64)
    if p.apply_mdl_barrier and env.get("MDL_BOUNDARY") is not None:
        mdl = env["MDL_BOUNDARY"].astype(bool)
        ff = env.get("MDL_FENCE_FIELD")
        if ff is not None:
            pass_f = np.clip(np.nan_to_num(np.asarray(ff, dtype=np.float64),
                                           nan=p.mdl_fence_factor), 1e-5, 1.0)
            if p.corridor_seasonal_gating:
                # Gates open with the dispersal season: water deer cross when they
                # are moving, not uniformly through the year.
                disp = np.asarray(p.disp_km_monthly, dtype=float)
                _dref = float(p.disp_scale_ref) if p.disp_scale_ref else float(disp.max())  # v9.10
                ds = float(disp[(month + int(p.season_offset)) % len(disp)] / (_dref + RAD_EPS))  # v9.9 FIX A
                s = p.corridor_season_floor + (1.0 - p.corridor_season_floor) * ds
                pass_f = pass_f * s
            bar = np.where(mdl, pass_f, bar)
        else:
            bar = np.where(mdl, p.mdl_fence_factor, bar)
    if p.apply_road_barrier and env.get("ROAD_MASK") is not None:
        bar = np.where(env["ROAD_MASK"].astype(bool), np.minimum(bar, p.road_factor), bar)
    return np.where(mask, bar, 0.0)
def build_fields(env, p: RADParams, month: int, closure: Callable):
    """Assemble the month's diffusion field D(x,y), advection (vx,vy) and the
    seasonal β_eff surface from the closure + seasonal drivers + barriers."""
    mask = env["KP"].astype(bool)
    M, Phi, beta_field = closure(env, p)
    bar = build_barrier(env, p, month)
    # Optional habitat-class mobility multiplier (Han-style habitat/accessible/
    # barrier). When absent (default) this is 1.0 and nothing changes; when the
    # host supplies it (USE_HABITAT_CLASS_IN_TRANSPORT), barrier cells get
    # near-zero diffusion/advection so the classification DRIVES the spread.
    cls_mob = env.get("CLASS_MOBILITY")
    if cls_mob is not None:
        bar = bar * np.nan_to_num(cls_mob.astype(np.float64), nan=0.0)
    # Seasonal dispersal scaling (0..1): more movement when the empirical gorani
    # dispersal radius is larger (peaks in early summer). Month index wraps every
    # 12 months so multi-year runs cycle the same seasonality.
    disp = np.asarray(p.disp_km_monthly, dtype=float)
    mi = (month + int(p.season_offset)) % len(disp)          # v9.9 FIX A: calendar phase
    _dref = float(p.disp_scale_ref) if p.disp_scale_ref else float(disp.max())  # v9.10
    disp_scale = float(disp[mi] / (_dref + RAD_EPS))
    # Diffusion field: mobility -> D, scaled by season & barriers.
    D = (p.diff_min + (p.diff_max - p.diff_min) * M) * disp_scale
    D = np.where(mask, D * bar, 0.0)
    # Advection: velocity climbs the attractiveness gradient, scaled by mobility,
    # season and barriers; clipped to the CFL speed cap.
    gy, gx = np.gradient(Phi)                     # gradient (row=y, col=x)
    vx = p.adv_gain * M * gx * disp_scale * bar
    vy = p.adv_gain * M * gy * disp_scale * bar
    spd = np.sqrt(vx * vx + vy * vy) + RAD_EPS
    scale = np.clip(p.adv_max / spd, 0.0, 1.0)
    vx = np.where(mask, vx * scale, 0.0)
    vy = np.where(mask, vy * scale, 0.0)
    # Seasonal β_eff surface.
    beta_eff = float(p.beta_monthly[mi]) * float(p.tvir_monthly[mi]) * beta_field
    beta_eff = np.where(mask, beta_eff, 0.0)
    return D.astype(np.float64), vx.astype(np.float64), vy.astype(np.float64), beta_eff.astype(np.float64)
# ══════════════════════════════════════════════════════════════════════════════
# NUMERICAL SOLVER  (flux-form FD, upwind advection, no-flux masked boundaries)
# ══════════════════════════════════════════════════════════════════════════════
def _transport_rate(I, D, vx, vy, mask):
    """dI/dt from advection + diffusion, flux form (mass-conserving), no-flux at
    the mask boundary. dx = 1 cell."""
    Im = I * mask
    # x-faces (between columns j and j+1)
    vxf_ok = mask[:, 1:] & mask[:, :-1]
    Dxf = 0.5 * (D[:, 1:] + D[:, :-1]) * vxf_ok
    vxf = 0.5 * (vx[:, 1:] + vx[:, :-1]) * vxf_ok
    Iup_x = np.where(vxf >= 0, Im[:, :-1], Im[:, 1:])
    Jx = vxf * Iup_x - Dxf * (Im[:, 1:] - Im[:, :-1])      # net +x flux at faces
    # y-faces (between rows i and i+1)
    vyf_ok = mask[1:, :] & mask[:-1, :]
    Dyf = 0.5 * (D[1:, :] + D[:-1, :]) * vyf_ok
    vyf = 0.5 * (vy[1:, :] + vy[:-1, :]) * vyf_ok
    Iup_y = np.where(vyf >= 0, Im[:-1, :], Im[1:, :])
    Jy = vyf * Iup_y - Dyf * (Im[1:, :] - Im[:-1, :])
    div = np.zeros_like(I)
    div[:, :-1] += Jx
    div[:, 1:] -= Jx
    div[:-1, :] += Jy
    div[1:, :] -= Jy
    return (-div) * mask
def _stable_dt(D, vx, vy, p: RADParams):
    dmax = float(D.max()) if D.size else 0.0
    vmax = float((np.abs(vx) + np.abs(vy)).max()) if vx.size else 0.0
    dt_diff = p.cfl_safety / (4.0 * dmax + RAD_EPS)     # 2-D explicit diffusion limit
    dt_adv = p.cfl_safety / (vmax + RAD_EPS)            # upwind advection limit
    dt = min(dt_diff, dt_adv, 1.0)
    nsub = int(min(max(math.ceil(1.0 / max(dt, RAD_EPS)), 1), p.max_substeps))
    return 1.0 / nsub, nsub
def rad_substep(S, I, R, D, vx, vy, mask, beta_eff, src, p: RADParams, dt: float):
    """One stable sub-step: transport, then reaction/recovery/source (Strang-lite
    operator split)."""
    # (1) transport of infected material
    I = I + dt * _transport_rate(I, D, vx, vy, mask)
    I = np.clip(I, 0.0, p.i_cap)
    # (2) reaction: local new infections deplete susceptibles (hazard form -> bounded)
    new_inf = S * (1.0 - np.exp(-np.clip(beta_eff * I * dt, 0.0, 50.0)))
    S = S - new_inf
    I = I + new_inf
    # (3) recovery / removal
    rec = I * (1.0 - math.exp(-p.gamma * dt))
    I = I - rec
    R = R + rec
    # (4) border source (imported infections drawn from local susceptibles)
    if src is not None:
        s_inf = np.minimum(src * dt, S)
        S = S - s_inf
        I = I + s_inf
    # bounds + mask
    S = np.clip(S, 0.0, 1.0) * mask
    I = np.clip(I, 0.0, p.i_cap) * mask
    R = np.clip(R, 0.0, 1.0) * mask
    return S, I, R
def simulate_month(S, I, R, env, p, month, closure, src_field):
    mask = env["KP"].astype(bool)
    D, vx, vy, beta_eff = build_fields(env, p, month, closure)
    dt, nsub = _stable_dt(D, vx, vy, p)
    src = src_field if (p.source_continuous and src_field is not None) else None
    for _ in range(nsub):
        S, I, R = rad_substep(S, I, R, D, vx, vy, mask, beta_eff, src, p, dt)
    return S, I, R, nsub
# ══════════════════════════════════════════════════════════════════════════════
# SEEDING
# ══════════════════════════════════════════════════════════════════════════════
def _nearest_valid(mask, r, c, max_radius=12):
    """Snap (r,c) to the nearest in-mask cell within max_radius (grid cells).
    Returns (r,c) or None if none found."""
    if 0 <= r < mask.shape[0] and 0 <= c < mask.shape[1] and mask[r, c]:
        return r, c
    for rad in range(1, max_radius + 1):
        r0, r1 = max(0, r - rad), min(mask.shape[0], r + rad + 1)
        c0, c1 = max(0, c - rad), min(mask.shape[1], c + rad + 1)
        sub = mask[r0:r1, c0:c1]
        if sub.any():
            rr, cc = np.where(sub)
            d = (rr + r0 - r) ** 2 + (cc + c0 - c) ** 2
            k = int(np.argmin(d))
            return int(rr[k] + r0), int(cc[k] + c0)
    return None
def seed_infection(shape, mask, entries, p: RADParams):
    """Initial infected field I0 at border-entry cells (row, col, ..., weight),
    and a persistent source field for optional continuous introduction. Entries
    that fall outside the valid mask are SNAPPED to the nearest valid cell so a
    slightly-off coordinate still starts the epidemic."""
    I0 = np.zeros(shape, dtype=np.float64)
    src = np.zeros(shape, dtype=np.float64)
    placed, snapped, lost = 0, 0, 0
    for e in entries or []:
        r, c = int(e[0]), int(e[1])
        w = float(e[4]) if len(e) >= 5 else 1.0
        rc = _nearest_valid(mask, r, c)
        if rc is None:
            lost += 1
            continue
        if rc != (r, c):
            snapped += 1
        rr, cc = rc
        I0[rr, cc] = max(I0[rr, cc], p.seed_intensity * w)
        src[rr, cc] = max(src[rr, cc], p.source_rate * w)
        placed += 1
    total = float(I0.sum())
    if placed == 0 or total <= 0:
        logger.warning(
            "  ⚠ RAD seeding placed NO infection (all %d entries fell outside the "
            "valid mask). The epidemic will be empty — check entry coordinates / "
            "rc_from_km conversion.", len(entries or []))
    elif getattr(p, "log_every", 1):                         # v9.10: quiet batch runs
        logger.info("  RAD seed: placed=%d (snapped=%d, lost=%d), total I0=%.3f",
                    placed, snapped, lost, total)
    return I0 * mask, src * mask
# ══════════════════════════════════════════════════════════════════════════════
# DRIVER
# ══════════════════════════════════════════════════════════════════════════════
def run_dynamic_rad_simulation(
    env: Dict[str, np.ndarray],
    p: RADParams,
    n_months: int = 12,
    closure: Optional[Callable] = None,
    stop_fn: Optional[Callable] = None,
    stop_extra_months: int = 0,
    max_months: Optional[int] = None,
    stop_label: str = "criterion",
) -> Dict:
    """Run the full dynamic epidemic. State persists month-to-month.
    Returns a dict with per-month S/I/R and cumulative attack fields, plus
    epidemic time-series and summary maps (see keys below)."""
    closure = closure or mechanistic_closure
    mask = env["KP"].astype(bool)
    shape = mask.shape
    cell_km = p.cell_km
    entries = list(env.get("DMZ_ENTRIES", [])) + list(env.get("NK_ENTRIES", []))
    I0, src = seed_infection(shape, mask, entries, p)
    S = np.where(mask, 1.0, 0.0).astype(np.float64)
    I = I0.copy()
    S = np.clip(S - I0, 0.0, 1.0)      # seed drawn from susceptibles
    R = np.zeros(shape, dtype=np.float64)
    V = np.zeros(shape, dtype=np.float64)      # vaccinated / immune compartment
    # Static spatial vaccination target (0..1), where 1 = fully eligible for the
    # campaign. Region targeting is baked into this by the host (env["VACC_FIELD"]).
    vacc_base = _get(env, "VACC_FIELD", shape, default=1.0) if p.vacc_enable else None
    monthly_I: List[np.ndarray] = []
    monthly_S: List[np.ndarray] = []
    monthly_R: List[np.ndarray] = []
    monthly_V: List[np.ndarray] = []
    monthly_cum: List[np.ndarray] = []     # attack rate = I + R (ever infected)
    ts = {k: [] for k in (
        "infected_area_km2", "susceptible_area_km2", "recovered_area_km2",
        "vaccinated_area_km2", "n_infected_cells", "growth_rate",
        "front_km", "mean_front_km", "substeps")}
    cell_area = cell_km * cell_km
    # distance-from-border field (km) for front metrics
    if _HAVE_SCIPY:
        border = np.zeros(shape, bool)
        for e in entries:
            r, c = int(e[0]), int(e[1])
            if 0 <= r < shape[0] and 0 <= c < shape[1]:
                border[r, c] = True
        if env.get("MDL_BOUNDARY") is not None:
            border |= env["MDL_BOUNDARY"].astype(bool)
        dist_border = distance_transform_edt(~border).astype(np.float64) * cell_km
    else:
        dist_border = np.zeros(shape)
    time_to_arrival = np.full(shape, np.nan, dtype=np.float64)
    cum_max = np.zeros(shape, dtype=np.float64)
    prev_cells = None
    # ── ADAPTIVE HORIZON (v9.2) ─────────────────────────────────────────────
    # n_months is the MINIMUM. When stop_fn is supplied the run keeps extending
    # (a year at a time, up to max_months) until the criterion fires, then runs
    # `stop_extra_months` further so the aftermath is visible. This is what lets
    # the model run "until PPR has actually crossed into South Korea" instead of
    # stopping at an arbitrary year and truncating the epidemic.
    m = 0
    target_months = max(1, int(n_months))
    hard_cap = max(target_months, int(max_months) if max_months else target_months)
    stop_hit = None
    vacc_on_month = None                                     # v9.9 E
    while m < target_months:
        S, I, R, nsub = simulate_month(S, I, R, env, p, m, closure, src)
        # ── Vaccination (monthly): move susceptibles S -> V in targeted cells ─
        # v9.9 E: either a fixed start month or a trigger (first detection in a target area)
        _vacc_active = False
        if p.vacc_enable and vacc_base is not None:
            if p.vacc_trigger_mask is not None:
                if vacc_on_month is None and monthly_cum:
                    _prev = np.nan_to_num(monthly_cum[-1])
                    if bool((_prev[np.asarray(p.vacc_trigger_mask, bool)] >= p.arrival_threshold).any()):
                        vacc_on_month = m
                _vacc_active = vacc_on_month is not None
            else:
                _vacc_active = m >= p.vacc_start_month
        if _vacc_active:
            if p.vacc_mode == "reactive_ring":
                inf = np.nan_to_num(I) > p.vacc_inf_thr
                ring = binary_dilation(inf, iterations=max(1, int(p.vacc_ring_cells)))
                tgt = vacc_base * ring.astype(np.float64)
            else:  # "preemptive" / "continuous" — vaccinate over the whole target field
                tgt = vacc_base
            _cov = (p.vacc_coverage if p.vacc_cap_field is None else
                    p.vacc_coverage * np.clip(np.nan_to_num(np.asarray(p.vacc_cap_field, float)), 0.0, 1.0))
            room = np.clip(_cov - V, 0.0, 1.0)                   # remaining capacity (v9.9 E: per-cell cap)
            dV = np.minimum(p.vacc_rate * p.vacc_efficacy * S * tgt, room)
            dV = np.where(mask, np.clip(dV, 0.0, S), 0.0)
            S = S - dV
            V = V + dV
        cum = np.where(mask, np.clip(I + R, 0.0, 1.0), np.nan)      # ever-infected (V excluded)
        if p.store_states:                                   # v9.10: batch runs keep cum only
            monthly_I.append(np.where(mask, I, np.nan).astype(np.float32))
            monthly_S.append(np.where(mask, S, np.nan).astype(np.float32))
            monthly_R.append(np.where(mask, R, np.nan).astype(np.float32))
            monthly_V.append(np.where(mask, V, np.nan).astype(np.float32))
        monthly_cum.append(cum.astype(np.float32))
        cum_max = np.maximum(cum_max, np.nan_to_num(cum))
        arrived = np.nan_to_num(cum) >= p.arrival_threshold
        newly = arrived & np.isnan(time_to_arrival)
        time_to_arrival[newly] = m
        inf_cells = int((np.nan_to_num(I) > 1e-3).sum())
        ts["n_infected_cells"].append(inf_cells)
        ts["infected_area_km2"].append(float(np.nan_to_num(I)[mask].sum()) * cell_area)
        ts["susceptible_area_km2"].append(float(S[mask].sum()) * cell_area)
        ts["recovered_area_km2"].append(float(R[mask].sum()) * cell_area)
        ts["vaccinated_area_km2"].append(float(V[mask].sum()) * cell_area)
        ts["substeps"].append(nsub)
        ts["growth_rate"].append(0.0 if prev_cells in (None, 0) else (inf_cells - prev_cells) / max(prev_cells, 1))
        prev_cells = inf_cells
        front = dist_border[arrived]
        ts["front_km"].append(float(front.max()) if front.size else 0.0)
        ts["mean_front_km"].append(float(front.mean()) if front.size else 0.0)
        if p.log_every and ((m + 1) % max(1, int(p.log_every)) == 0):     # v9.10
            logger.info(
                f"  [RAD month {m + 1:3d}/{target_months}] substeps={nsub:2d}  "
                f"infected_cells={inf_cells:6d}  front={ts['front_km'][-1]:6.1f} km  "
                f"attack(mean)={float(np.nan_to_num(cum)[mask].mean()):.3f}"
                f"{('  vacc(mean)=%.3f' % float(V[mask].mean())) if p.vacc_enable else ''}")
        # ── adaptive horizon bookkeeping ────────────────────────────────────
        if stop_fn is not None and stop_hit is None:
            try:
                _fired = bool(stop_fn(m, np.nan_to_num(cum)))
            except Exception as _e:
                logger.warning(f"  ⚠ stop_fn raised ({_e}); continuing to the minimum horizon.")
                _fired = False
            if _fired:
                stop_hit = m
                # v9.9 FIX B: never stop before the stated minimum horizon
                target_months = min(hard_cap, max(target_months, m + 1 + int(stop_extra_months)))
                if p.log_every:                                      # v9.10
                    logger.info(f"  ✓ {stop_label} MET at month {m + 1}; running "
                                f"{int(stop_extra_months)} more months "
                                f"(to month {target_months}) to show the aftermath.")
        m += 1
        if (stop_fn is not None and stop_hit is None
                and m >= target_months and target_months < hard_cap):
            target_months = min(hard_cap, target_months + 12)
            if p.log_every:                                          # v9.10
                logger.info(f"  … {stop_label} not met by month {m}; extending the run to "
                            f"{target_months} months (cap {hard_cap}).")
    if stop_fn is not None and stop_hit is None and p.log_every:
        logger.warning(f"  ⚠ {stop_label} NEVER met within the {hard_cap}-month cap "
                       f"({hard_cap / 12.0:.0f} yr). The epidemic did not get far enough — "
                       f"raise MAX_SIM_YEARS, or speed transport up "
                       f"(SPREAD_SPEED_MULTIPLIER / CONDUCTANCE_EXPONENT).")
    return {
        "monthly_I": monthly_I,
        "monthly_S": monthly_S,
        "monthly_R": monthly_R,
        "monthly_V": monthly_V,
        "monthly_cum": monthly_cum,          # attack rate I + R per month (ever infected)
        "time_series": {k: np.asarray(v) for k, v in ts.items()},
        "cum_max": cum_max.astype(np.float32),
        "time_to_arrival": time_to_arrival.astype(np.float32),  # month index (NaN = never)
        "final_R": np.where(mask, R, np.nan).astype(np.float32),
        "final_V": np.where(mask, V, np.nan).astype(np.float32),
        "dist_border_km": dist_border.astype(np.float32),
        "vacc_on_month": vacc_on_month,                     # v9.9 E (0-based leg month, or None)
        "stop_month": stop_hit,                             # 0-based month the stop criterion fired
    }
# ══════════════════════════════════════════════════════════════════════════════
# v9.10 — REVIEW RE-ANALYSIS HELPERS (search "v9.10" for every change)
#   1. ROK/DPRK split from the ROK sigungu sub-cells (+ land share at the border),
#      exclaves (islands joined to the DPRK coast on the grid), estuary policy,
#      establishment judged on ROK area
#   2. distance-only benchmarks and alternative weights for the SMCDA check
#   3. fence wait: front on the DPRK side of each gate vs first ROK cell
#   4. Tables 1-3 as CSV, footprint statistics, Figure 9
#   5. global sensitivity analysis: Morris screening -> LHS + PRCC (fixed seed)
#   6. Table 4 scenarios on the peninsular leg; 4-km grid preset (CLI)
#   7. GeoTIFF export of risk scores / arrival months for an external test
# Everything here is numpy/scipy only and is called from main() or the CLI.
# ══════════════════════════════════════════════════════════════════════════════
import json as _v910_json
import time as _v910_time
from scipy.ndimage import label as _nd_label

V910_TAG = "PPR_v9.10"
V910_FOUR = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
V910_EIGHT = np.ones((3, 3), dtype=bool)
V910_MONTHS = ["January", "February", "March", "April", "May", "June", "July",
               "August", "September", "October", "November", "December"]

# First two digits of the ROK sigungu code (SIG_CD) -> (Korean, English) province.
# sig.shp carries no province field, so the county CSVs had an empty column.
SIDO_BY_CODE = {
    "11": ("서울특별시", "Seoul"), "26": ("부산광역시", "Busan"),
    "27": ("대구광역시", "Daegu"), "28": ("인천광역시", "Incheon"),
    "29": ("광주광역시", "Gwangju"), "30": ("대전광역시", "Daejeon"),
    "31": ("울산광역시", "Ulsan"), "36": ("세종특별자치시", "Sejong"),
    "41": ("경기도", "Gyeonggi-do"), "42": ("강원도", "Gangwon-do"),
    "51": ("강원특별자치도", "Gangwon-do"), "43": ("충청북도", "Chungcheongbuk-do"),
    "44": ("충청남도", "Chungcheongnam-do"), "45": ("전라북도", "Jeollabuk-do"),
    "52": ("전북특별자치도", "Jeollabuk-do"), "46": ("전라남도", "Jeollanam-do"),
    "47": ("경상북도", "Gyeongsangbuk-do"), "48": ("경상남도", "Gyeongsangnam-do"),
    "50": ("제주특별자치도", "Jeju-do"),
}


def province_from_sig_cd(code):
    """(English, Korean) province from the first two digits of a sigungu code."""
    c = str(code or "").strip()
    if len(c) >= 2 and c[:2] in SIDO_BY_CODE:
        kr, en = SIDO_BY_CODE[c[:2]]
        return en, kr
    return "", ""


def v910_calendar(g):
    """'February, Year 14' for global month g (month 1 = January of Year 1)."""
    if g is None:
        return ""
    try:
        g = int(g)
    except (TypeError, ValueError):
        return ""
    if g < 1:
        return ""
    return f"{V910_MONTHS[(g - 1) % 12]}, Year {(g - 1) // 12 + 1}"


def v910_write_csv(path, header, rows, meta_lines=(), logger_=None):
    """UTF-8-with-BOM CSV (Excel-safe for Hangul); '#' lines carry provenance."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        for ln in meta_lines:
            fh.write(f"# {ln}\n")
        w = csv.writer(fh)
        w.writerow(list(header))
        for r in rows:
            if isinstance(r, dict):
                w.writerow([_fmt(r.get(k)) for k in header])
            else:
                w.writerow([_fmt(v) for v in r])
    if logger_ is not None:
        logger_.info(f"   ✓ {path.name} ({len(rows)} rows)")
    return path


def v910_copy_to_dl(paths, dl_dir, logger_=None):
    if not dl_dir:
        return
    try:
        Path(dl_dir).mkdir(parents=True, exist_ok=True)
        for p in paths:
            if p is not None and Path(p).exists():
                shutil.copy2(p, Path(dl_dir) / Path(p).name)
    except Exception as e:
        if logger_ is not None:
            logger_.warning(f"   ⚠ could not copy to {dl_dir}: {e}")


# ── ranks and correlations (tie-aware) ─────────────────────────────────────────
def v910_rank_avg(x):
    """1-based ranks with ties given their average rank."""
    x = np.asarray(x, dtype=np.float64).ravel()
    if x.size == 0:
        return x
    _u, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    csum = np.cumsum(cnt)
    start = csum - cnt
    return ((start + csum + 1) / 2.0)[inv]


def spearman_ties(a, b):
    """Spearman's rho with average ranks for ties. Arrival months are integers,
    so ties are the rule; the v9.7 helper broke them by array order."""
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if a.size < 3:
        return float("nan")
    ra = v910_rank_avg(a)
    rb = v910_rank_avg(b)
    ra -= ra.mean()
    rb -= rb.mean()
    den = math.sqrt(float((ra * ra).sum()) * float((rb * rb).sum()))
    return float((ra * rb).sum() / den) if den > 0 else float("nan")


def first_arrival_index(cum_list, thr, mask=None):
    """0-based month index at which each cell first reaches A >= thr (NaN = never)."""
    if not cum_list:
        return None
    tta = np.full(np.asarray(cum_list[0]).shape, np.nan, dtype=np.float64)
    for m, cum in enumerate(cum_list):
        hit = np.nan_to_num(np.asarray(cum, dtype=np.float64)) >= float(thr)
        new = hit & np.isnan(tta)
        if new.any():
            tta[new] = m
    if mask is not None:
        tta = np.where(np.asarray(mask, bool), tta, np.nan)
    return tta


def largest_remainder_shares(counts, decimals=1):
    """Percentages rounded to `decimals` that add up to exactly 100."""
    counts = np.asarray(counts, dtype=np.float64)
    tot = float(counts.sum())
    if tot <= 0:
        return [0.0] * len(counts)
    scale = 10 ** int(decimals)
    raw = counts / tot * 100.0 * scale
    base = np.floor(raw)
    short = int(round(100 * scale - base.sum()))
    order = np.argsort(-(raw - base), kind="mergesort")
    for i in order[:max(short, 0)]:
        base[i] += 1
    return [float(v) / scale for v in base]


# ══════════════════════════════════════════════════════════════════════════════
# ITEM 1 — NATIONAL SPLIT
# ══════════════════════════════════════════════════════════════════════════════
def rok_fraction_from_zones(zones, nr, nc, oversample):
    """Fraction (0-1) of each model cell covered by ROK sigungu polygons, from the
    sub-cell weights build_county_zones() already computed for the county CSVs
    (that build reprojects every polygon itself, so it works where
    gdal.VectorTranslate failed in the 6 Nov run)."""
    frac = np.zeros(int(nr) * int(nc), dtype=np.float64)
    sub = float(max(1, int(oversample or 1))) ** 2
    for z in zones or []:
        cf = np.asarray(z.get("cells_flat", []), dtype=np.int64)
        if cf.size == 0:
            continue
        if z.get("sampling") == "zonal":
            np.add.at(frac, cf, np.asarray(z.get("weights", []), dtype=np.float64) / sub)
        else:                      # smaller than one sub-cell: give it one sub-cell
            frac[cf] = np.maximum(frac[cf], 1.0 / sub)
    return np.clip(frac.reshape(int(nr), int(nc)), 0.0, 1.0)


def land_fraction_on_grid(tif_path, gt, nr, nc, wkt):
    """Fraction (0-1) of each model cell covered by valid pixels of a raster
    (the 500-m water-deer layer marks land). Used only to judge border cells."""
    ds = gdal.Open(str(tif_path))
    band = ds.GetRasterBand(1)
    a = band.ReadAsArray().astype(np.float64)
    nd = band.GetNoDataValue()
    ok = np.isfinite(a) & (a > -3.4e38 * 0.01) & (a >= 0)
    if nd is not None:
        ok &= ~np.isclose(a, nd, rtol=0, atol=1e-3)
    mem = gdal.GetDriverByName("MEM").Create("", ds.RasterXSize, ds.RasterYSize, 1,
                                             gdal.GDT_Float32)
    mem.SetGeoTransform(ds.GetGeoTransform())
    mem.SetProjection(ds.GetProjection() or wkt)
    mem.GetRasterBand(1).WriteArray(ok.astype(np.float32))
    ds = None
    wo = gdal.WarpOptions(format="MEM", dstSRS=wkt,
                          outputBounds=(gt[0], gt[3] + nr * gt[5], gt[0] + nc * gt[1], gt[3]),
                          width=int(nc), height=int(nr), resampleAlg=gdal.GRA_Average)
    out = gdal.Warp("", mem, options=wo)
    lf = out.GetRasterBand(1).ReadAsArray().astype(np.float64)
    out = None
    mem = None
    return np.clip(np.nan_to_num(lf, nan=0.0), 0.0, 1.0)


def split_rok_dprk(korea_region, rok_frac, land_frac=None, core_eps=0.02,
                   border_min_land_share=0.5):
    """Assign every cell of the Korean region to the ROK or the DPRK.

    The ROK share of a cell is its ROK area divided by its LAND area (from the
    land layer), so a coastal ROK cell that is mostly sea still has a share near 1.
    Without a land layer the share is the ROK share of the whole cell.
    1. Cells whose ROK share is (almost) zero are candidate DPRK land; the
       largest 8-connected component is the DPRK mainland. Other components
       (islands) join the DPRK if they lie nearer the DPRK mainland than ROK land.
    2. Border cells — not wholly ROK and touching the DPRK — go to the country
       holding most of their land (share >= border_min_land_share -> ROK).
    3. Every other cell holding ROK area is ROK.
    Returns (rok_mask, info)."""
    kr = np.asarray(korea_region, bool)
    f = np.clip(np.nan_to_num(np.asarray(rok_frac, dtype=np.float64)), 0.0, 1.0)
    info = dict(core_eps=float(core_eps), border_min_land_share=float(border_min_land_share),
                land_fraction_used=land_frac is not None)
    if land_frac is not None:
        lf = np.maximum(np.clip(np.nan_to_num(np.asarray(land_frac, dtype=np.float64)), 0.0, 1.0), f)
        share = np.where(lf > 1e-9, f / np.maximum(lf, 1e-9), 0.0)
    else:
        share = f
    cand = kr & (share < float(core_eps))
    lab, n = _nd_label(cand, structure=V910_EIGHT)
    if n == 0:
        rok = kr & (f > 0)
        info.update(n_rok=int(rok.sum()), n_dprk=0, n_border_to_dprk=0, n_islands_to_dprk=0)
        return rok, info
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    main_id = int(np.argmax(sizes))
    dprk = lab == main_id
    rok_core = kr & (share >= 0.5)
    n_isl = 0
    if n > 1:
        d_dprk = distance_transform_edt(~dprk)
        d_rok = (distance_transform_edt(~rok_core) if rok_core.any()
                 else np.full(f.shape, np.inf))
        for cid in range(1, n + 1):
            if cid == main_id:
                continue
            comp = lab == cid
            if float(d_dprk[comp].min()) <= float(d_rok[comp].min()):
                dprk |= comp
                n_isl += int(comp.sum())
    mixed = kr & ~dprk & (share < 1.0 - 1e-9)
    touch = binary_dilation(dprk, structure=V910_EIGHT) & mixed
    to_dprk = touch & (share < float(border_min_land_share))
    dprk |= to_dprk
    rok = kr & ~dprk
    info.update(n_rok=int(rok.sum()), n_dprk=int((kr & dprk).sum()),
                n_border_cells=int(touch.sum()), n_border_to_dprk=int(to_dprk.sum()),
                n_islands_to_dprk=int(n_isl), n_noncore_components=int(n - 1))
    return rok, info


def rok_exclave_cells(rok, dprk):
    """ROK cells in 4-connected ROK components other than the mainland that share
    a cell face with DPRK land. On the grid, infection reaches them only from the
    DPRK side, bypassing the fence (e.g. the northern Ongjin-gun islands)."""
    rok = np.asarray(rok, bool)
    lab, n = _nd_label(rok, structure=V910_FOUR)
    out = np.zeros(rok.shape, bool)
    if n <= 1:
        return out, []
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    main_id = int(np.argmax(sizes))
    touch = binary_dilation(np.asarray(dprk, bool), structure=V910_FOUR)
    comps = []
    for cid in range(1, n + 1):
        if cid == main_id:
            continue
        comp = lab == cid
        if (comp & touch).any():
            out |= comp
            comps.append(int(comp.sum()))
    return out, comps


def unfenced_contact(rok, dprk, fence_band, buffer_cells=1):
    """(ROK cells, DPRK cells) that share a face across the national border away
    from the fence band — the unfenced contact west of the fence terminus.
    buffer_cells=1 (reporting) also leaves out cells next to the band;
    buffer_cells=0 (sealing) keeps every face outside the band."""
    rok = np.asarray(rok, bool)
    dprk = np.asarray(dprk, bool)
    band = np.asarray(fence_band, bool) if fence_band is not None else np.zeros(rok.shape, bool)
    if int(buffer_cells) > 0 and band.any():
        near = binary_dilation(band, structure=V910_FOUR, iterations=int(buffer_cells))
    else:
        near = band
    rc = rok & binary_dilation(dprk, structure=V910_FOUR) & ~near
    dc = dprk & binary_dilation(rok, structure=V910_FOUR) & ~near
    return rc, dc


def seal_contact_band(mdl, fence, contact, sealed, default_phi):
    """Add contact cells to the fence band with the sealed pass-through.
    Both sides of the contact are sealed, so each crossing face is barred on
    both cells (the transport uses face-averaged coefficients)."""
    mdl = np.asarray(mdl, bool)
    contact = np.asarray(contact, bool)
    ff = (np.where(mdl, float(default_phi), 1.0) if fence is None
          else np.asarray(fence, dtype=np.float64).copy())
    ff = np.where(contact, float(sealed), ff)
    mdl2 = mdl | contact
    return mdl2, np.where(mdl2, np.clip(ff, 1e-5, 1.0), 1.0)


def clip_segments_to_mask(segs, mask):
    """Keep the runs of each polyline (x=col, y=row, crop units) whose vertices
    fall in `mask`; split at the gaps so no chord crosses the border."""
    mask = np.asarray(mask, bool)
    H, W = mask.shape
    out = []
    for s in segs or []:
        s = np.asarray(s, dtype=float)
        if s.ndim != 2 or s.shape[0] < 2:
            continue
        r = np.floor(s[:, 1]).astype(int)
        c = np.floor(s[:, 0]).astype(int)
        ok = (r >= 0) & (r < H) & (c >= 0) & (c < W)
        inside = np.zeros(s.shape[0], bool)
        inside[ok] = mask[r[ok], c[ok]]
        idx = np.flatnonzero(inside)
        if idx.size < 2:
            continue
        for run in np.split(idx, np.flatnonzero(np.diff(idx) != 1) + 1):
            if run.size >= 2:
                out.append(s[run])
    return out


# ══════════════════════════════════════════════════════════════════════════════
# ITEMS 3–4 — KOREA-LEG METRICS, FENCE WAIT, FIRST CROSSING
# ══════════════════════════════════════════════════════════════════════════════
def korea_leg_metrics(cum_list, nk_mask, rok_mask, rok_w=None, thr=0.05,
                      month_offset=0, est_mean=0.05, est_area=0.10, report_months=(),
                      first_exclude=None):
    """ROK/DPRK summary of one peninsular-leg run on the global month index.
    ROK means and areas are AREA-WEIGHTED by the ROK share of each cell (as in the
    sigungu tables); the *_cell variants use the cells assigned to the ROK.
    `first_exclude` (ROK exclave cells) is left out of the first-arrival months,
    which must describe a crossing of the border, not an island reached from the
    DPRK coast."""
    off = int(month_offset)
    nk = np.asarray(nk_mask, bool)
    rok = np.asarray(rok_mask, bool)
    fx = (np.zeros(rok.shape, bool) if first_exclude is None
          else np.asarray(first_exclude, bool))
    w = (rok.astype(np.float64) if rok_w is None
         else np.clip(np.nan_to_num(np.asarray(rok_w, dtype=np.float64)), 0.0, 1.0))
    wsum = float(w.sum())
    NM = len(cum_list)
    keys = ("rok_mean_w", "rok_area_w", "rok_mean_cell", "rok_area_cell", "nk_mean",
            "rok_any_w", "rok_any_cell")
    S = {k: np.full(NM, np.nan) for k in keys}
    wpos = (w > 0) & ~fx
    rok_first = rok & ~fx
    for m in range(NM):
        c = np.nan_to_num(np.asarray(cum_list[m], dtype=np.float64))
        hit = c >= float(thr)
        if wsum > 0:
            S["rok_mean_w"][m] = float((w * c).sum()) / wsum
            S["rok_area_w"][m] = float((w * hit).sum()) / wsum
        if rok.any():
            S["rok_mean_cell"][m] = float(c[rok].mean())
            S["rok_area_cell"][m] = float(hit[rok].mean())
        if nk.any():
            S["nk_mean"][m] = float(c[nk].mean())
        S["rok_any_w"][m] = float(bool((hit & wpos).any()))
        S["rok_any_cell"][m] = float(bool((hit & rok_first).any()))

    def first(cond):
        idx = np.flatnonzero(np.nan_to_num(cond.astype(float)) > 0)
        return int(idx[0]) + 1 + off if idx.size else None

    def ge(k, v):
        return np.nan_to_num(S[k], nan=-1.0) >= v

    out = dict(
        n_months=NM, end_month=NM + off, intro_month=off + 1, threshold=float(thr),
        first_rok_any=first(S["rok_any_w"] > 0),
        first_rok_cell=first(S["rok_any_cell"] > 0),
        area_1=first(ge("rok_area_w", 0.01)), area_5=first(ge("rok_area_w", 0.05)),
        area_10=first(ge("rok_area_w", 0.10)),
        mean_1=first(ge("rok_mean_w", 0.01)), mean_5=first(ge("rok_mean_w", 0.05)),
        est_areaw=first(ge("rok_mean_w", est_mean) | ge("rok_area_w", est_area)),
        est_cell=first(ge("rok_mean_cell", est_mean) | ge("rok_area_cell", est_area)),
        nk_5=first(ge("nk_mean", 0.05)),
    )

    def at(g):
        mi = int(g) - 1 - off
        if not (0 <= mi < NM):
            return None
        return dict(month=int(g),
                    rok_mean_w_pct=100.0 * float(S["rok_mean_w"][mi]),
                    rok_area_w_pct=100.0 * float(S["rok_area_w"][mi]),
                    rok_mean_cell_pct=100.0 * float(S["rok_mean_cell"][mi]),
                    rok_area_cell_pct=100.0 * float(S["rok_area_cell"][mi]),
                    nk_mean_pct=100.0 * float(S["nk_mean"][mi]))
    out["at_end"] = at(NM + off)
    out["at"] = {int(g): at(g) for g in (report_months or []) if at(g) is not None}
    if out["est_areaw"] is not None:
        out["at_est"] = at(out["est_areaw"])
    return out, S


def fence_wait_table(tta, corridors, rok, dprk, fence_band, cell_km, probe_km=15.0,
                     month_offset=0, contact_rok=None):
    """H1: how long the front waited at the fence. For each gate (and the unfenced
    estuary and the whole fence) the first and median month at which DPRK cells
    and ROK cells within `probe_km` reached A*, on the global month index.
    `tta` is the 0-based leg-month of arrival; corridors are in the same grid."""
    t = np.asarray(tta, dtype=np.float64)
    H, W = t.shape
    rr, cc = np.mgrid[0:H, 0:W]
    band = np.asarray(fence_band, bool) if fence_band is not None else np.zeros((H, W), bool)
    rok_all = np.asarray(rok, bool)
    rok = rok_all & ~band
    dprk = np.asarray(dprk, bool) & ~band
    off = int(month_offset)

    def glob(v):
        return None if v is None else int(round(v)) + 1 + off

    def first(m):
        v = t[m]
        v = v[np.isfinite(v)]
        return float(v.min()) if v.size else None

    def median_all(m):
        v = t[m]
        if v.size == 0:
            return None
        v = np.sort(np.where(np.isfinite(v), v, np.inf))
        mid = v[(v.size - 1) // 2]
        return float(mid) if np.isfinite(mid) else None

    def first_cells(m, v0):
        if v0 is None:
            return []
        return [tuple(int(x) for x in rc) for rc in np.argwhere(m & (t == v0))]

    rows = []

    def add(name, weight, zone, kind):
        north, south, bz = zone & dprk, zone & rok, zone & band
        nf, sf, bf = first(north), first(south), first(bz)
        rbf = first(zone & band & rok_all)
        nm, sm = median_all(north), median_all(south)
        rows.append(dict(
            zone=name, kind=kind, weight=weight,
            n_dprk_cells=int(north.sum()), n_band_cells=int(bz.sum()),
            n_rok_cells=int(south.sum()),
            dprk_first_month=glob(nf), dprk_median_month=glob(nm),
            band_first_month=glob(bf), rok_band_first_month=glob(rbf),
            rok_first_month=glob(sf), rok_median_month=glob(sm),
            wait_first_months=(None if (nf is None or sf is None) else int(round(sf - nf))),
            wait_median_months=(None if (nm is None or sm is None) else int(round(sm - nm))),
            rok_first_calendar=v910_calendar(glob(sf)),
            rok_first_cells=first_cells(south, sf)))

    rad = float(probe_km) / max(float(cell_km), 1e-9)
    for co in corridors or []:
        zone = ((rr - float(co["row"])) ** 2 + (cc - float(co["col"])) ** 2) <= rad * rad
        add(str(co["name"]), float(co.get("weight", float("nan"))), zone, "gate")
    if contact_rok is not None and np.any(contact_rok):
        d = distance_transform_edt(~np.asarray(contact_rok, bool))
        add("Han River estuary (unfenced contact)", float("nan"), d <= rad, "estuary")
    if band.any():
        ring = binary_dilation(band, structure=V910_FOUR) & ~band
        add("whole fence (cells next to the band)", float("nan"), ring, "fence")
    return rows


def first_crossing_location(tta, rok_main, corridors, contact_rok, cell_km, probe_km=15.0):
    """Where the first ROK cell(s) were reached: a gate name, 'estuary' or 'other'."""
    t = np.where(np.asarray(rok_main, bool), np.asarray(tta, dtype=np.float64), np.nan)
    if not np.isfinite(t).any():
        return dict(label="not reached", dist_km=None, row=None, col=None,
                    month_index=None, n_cells=0, labels="")
    m0 = float(np.nanmin(t))
    cells = np.argwhere(t == m0)
    dcon = (distance_transform_edt(~np.asarray(contact_rok, bool)) * float(cell_km)
            if (contact_rok is not None and np.any(contact_rok)) else None)
    best, labels = None, set()
    for r, c in cells:
        dg, gname = float("inf"), ""
        for co in corridors or []:
            d = math.hypot(r - float(co["row"]), c - float(co["col"])) * float(cell_km)
            if d < dg:
                dg, gname = d, str(co["name"])
        de = float(dcon[r, c]) if dcon is not None else float("inf")
        if dg <= probe_km and dg <= de:
            lab, d = gname, dg
        elif de <= probe_km:
            lab, d = "estuary", de
        else:
            lab, d = "other", min(dg, de)
        labels.add(lab)
        if best is None or d < best[1]:
            best = (lab, d, int(r), int(c))
    return dict(label=best[0], dist_km=float(best[1]), row=best[2], col=best[3],
                month_index=m0, n_cells=int(len(cells)), labels="|".join(sorted(labels)))


# ══════════════════════════════════════════════════════════════════════════════
# ITEM 2 — SMCDA BENCHMARKS AND ALTERNATIVE WEIGHTS
# ══════════════════════════════════════════════════════════════════════════════
def _smcda_aggregate_like(sm, std, weights):
    """Re-aggregate a standardised stack exactly as the reported classification
    was (same aggregator, same vetoes)."""
    mask = np.asarray(sm["mask"], bool)
    ag = str(Config.SMCDA_AGGREGATION).lower()
    if ag == "owa":
        score = owa_aggregate(std, weights, mask, andness=float(Config.SMCDA_OWA_ANDNESS))
    else:
        score = wlc_aggregate(std, weights, mask)
    for c in sm.get("criteria", []):
        if getattr(c, "veto_below", None) is None or c.key not in std:
            continue
        score = np.where(mask & (std[c.key] < float(c.veto_below)), 0.0, score)
    return score


def _bench_row(name, kind, score, sm, tta, reached, month_offset, base_cls=None,
               weights=None, note=""):
    mask = np.asarray(sm["mask"], bool)
    s = np.asarray(score, dtype=np.float64)
    sel = mask & np.isfinite(tta) & np.isfinite(s)
    row = dict(benchmark=name, kind=kind, note=note,
               rho_s=spearman_ties(s[sel], tta[sel]), n_arrived=int(sel.sum()))
    if reached is not None:
        m2 = mask & np.isfinite(s)
        y = reached[m2]
        row["auc_reached"] = (_auc(s[m2], y) if (y.any() and (~y).any()) else float("nan"))
    if base_cls is not None:
        br = compute_breaks(s, mask, k=int(sm.get("n_classes", 5)), method=sm.get("break_method", "jenks"),
                            absolute=Config.SMCDA_ABSOLUTE_BREAKS)
        cls = classify(s, br["edges"], mask)
        both = mask & np.isfinite(cls) & np.isfinite(base_cls)
        row["pct_same_class"] = (100.0 * float((cls[both] == base_cls[both]).mean())
                                 if both.any() else float("nan"))
        arr = []
        for i in range(int(sm.get("n_classes", 5))):
            sel_i = mask & (cls == i) & np.isfinite(tta)
            arr.append(float(tta[sel_i].mean()) + 1 + month_offset if sel_i.any() else None)
        fin = [a for a in arr if a is not None]
        row["mean_arrival_by_class"] = "|".join("" if a is None else f"{a:.1f}" for a in arr)
        row["monotonic"] = bool(len(fin) >= 3 and np.all(np.diff(fin) <= 1e-9))
    if weights:
        row["weights"] = "; ".join(f"{k}={float(v):.3f}" for k, v in weights.items())
    return row


def smcda_benchmarks(sm, tta, cum, thr, month_offset, origins, cost_dist, cell_km):
    """Item 2. Rho_s (tie-aware) between each ranking and simulated arrival over
    the domain, for: the reported composite; the composite with ENTRY measured
    from other origins; distance-only rankings (Euclidean and cost distance from
    each origin); each criterion alone; and alternative weights (entropy, equal,
    rank-sum, leave-one-out). A distance-only ranking that orders arrival as well
    as the composite means the extra criteria add little to the ranking."""
    rows = []
    mask = np.asarray(sm["mask"], bool)
    tta = np.asarray(tta, dtype=np.float64)
    reached = (np.nan_to_num(np.asarray(cum, dtype=np.float64)) >= float(thr)
               if cum is not None else None)
    std = dict(sm["standardised"])
    w0 = dict(sm["weights"])
    base_cls = np.asarray(sm["classes"], dtype=np.float64)
    rows.append(_bench_row("composite (AHP, as reported)", "composite", sm["score"], sm, tta,
                           reached, month_offset, base_cls, w0, sm.get("entry_basis", "")))
    entry_c = next((c for c in sm.get("criteria", []) if c.note == "_SMCDA_GATE_DIST"), None)
    for oname, cd in (cost_dist or {}).items():
        if cd is None:
            continue
        rows.append(_bench_row(f"cost distance from {oname} (alone)", "distance",
                               -np.asarray(cd, dtype=np.float64), sm, tta, reached, month_offset))
        if entry_c is not None and entry_c.key in std:
            std2 = dict(std)
            std2[entry_c.key] = standardise_criterion(np.asarray(cd, dtype=np.float64), mask, entry_c)
            rows.append(_bench_row(f"composite, ENTRY from {oname}", "composite",
                                   _smcda_aggregate_like(sm, std2, w0), sm, tta, reached,
                                   month_offset, base_cls, w0))
            rows.append(_bench_row(f"ENTRY criterion alone, from {oname}", "criterion",
                                   std2[entry_c.key], sm, tta, reached, month_offset))
    for oname, om in (origins or {}).items():
        if om is None or not np.any(om):
            continue
        d = distance_transform_edt(~np.asarray(om, bool)) * float(cell_km)
        rows.append(_bench_row(f"Euclidean distance from {oname} (alone)", "distance",
                               -d, sm, tta, reached, month_offset))
    for k in std:
        rows.append(_bench_row(f"criterion {k} alone", "criterion", std[k], sm, tta,
                               reached, month_offset))
    keys = list(std.keys())
    alts = {"equal weights": {k: 1.0 / len(keys) for k in keys},
            "entropy weights": entropy_weights(std, mask),
            "rank-sum weights (AHP order)": rank_sum_weights(sorted(keys, key=lambda k: -w0.get(k, 0)))}
    for k in keys:
        wl = {kk: vv for kk, vv in w0.items() if kk != k}
        s = sum(wl.values())
        if s > 0:
            alts[f"AHP without {k}"] = {kk: vv / s for kk, vv in wl.items()}
    for nm, ww in alts.items():
        rows.append(_bench_row(nm, "weights", _smcda_aggregate_like(sm, std, ww), sm, tta,
                               reached, month_offset, base_cls, ww))
    return rows


# ══════════════════════════════════════════════════════════════════════════════
# ROUND 2 (Reviewer 2, F1) — WHAT DOES A RANKING ADD BEYOND DISTANCE?
# ══════════════════════════════════════════════════════════════════════════════
BEYOND_LAMBDAS = np.round(np.linspace(0.0, 1.0, 21), 2)


def _rank_residual(y, x):
    """Residual of y after a straight-line fit on x (used on ranks)."""
    z = np.column_stack([np.ones(y.size), x])
    return y - z @ np.linalg.lstsq(z, y, rcond=None)[0]


def partial_spearman(a, b, c):
    """Spearman correlation of a and b with c held fixed: the Pearson correlation of the
    residuals of rank(a) and rank(b), each regressed on rank(c) (average ranks for ties)."""
    a, b, c = (np.asarray(v, dtype=np.float64).ravel() for v in (a, b, c))
    ok = np.isfinite(a) & np.isfinite(b) & np.isfinite(c)
    if ok.sum() < 5:
        return float("nan")
    rc = v910_rank_avg(c[ok])
    ea = _rank_residual(v910_rank_avg(a[ok]), rc)
    eb = _rank_residual(v910_rank_avg(b[ok]), rc)
    den = math.sqrt(float((ea * ea).sum()) * float((eb * eb).sum()))
    return float((ea * eb).sum() / den) if den > 0 else float("nan")


def rank_r2(y, *xs):
    """R² of rank(y) regressed (OLS) on the ranks of xs."""
    y = np.asarray(y, dtype=np.float64).ravel()
    ok = np.isfinite(y)
    for x in xs:
        ok &= np.isfinite(np.asarray(x, dtype=np.float64).ravel())
    if ok.sum() < 5:
        return float("nan")
    ry = v910_rank_avg(y[ok])
    z = np.column_stack([np.ones(ry.size)] + [v910_rank_avg(np.asarray(x, np.float64).ravel()[ok]) for x in xs])
    res = ry - z @ np.linalg.lstsq(z, ry, rcond=None)[0]
    tss = float(((ry - ry.mean()) ** 2).sum())
    return 1.0 - float((res * res).sum()) / tss if tss > 0 else float("nan")


def _pct_rank(x):
    """Percentile rank in (0, 1), average ranks for ties."""
    r = v910_rank_avg(x)
    return (r - 0.5) / float(max(r.size, 1))


def _best_mixture(u_ref, u_r, y, lambdas=BEYOND_LAMBDAS):
    """(lambda, rho) of the mixture (1-lambda)*u_ref + lambda*u_r with the most negative
    Spearman rho with arrival month; lambda = 0 is the reference alone."""
    best_l, best_r = 0.0, spearman_ties(u_ref, y)
    for lam in lambdas[1:]:
        r = spearman_ties((1.0 - lam) * u_ref + lam * u_r, y)
        if np.isfinite(r) and r < best_r - 1e-12:
            best_l, best_r = float(lam), float(r)
    return best_l, best_r


def _beyond_stats(y, d, r):
    """Statistics of one ranking r against arrival y, with reference distance d (km)."""
    rho = spearman_ties(r, y)
    rho_d = spearman_ties(-d, y)
    part = partial_spearman(r, y, d)
    r2_d = rank_r2(y, d)
    r2_dr = rank_r2(y, d, r)
    lam, rho_m = _best_mixture(_pct_rank(-d), _pct_rank(r), y)
    return dict(rho_s=rho, rho_s_reference=rho_d, partial_rho_given_distance=part,
                r2_reference=r2_d, r2_reference_plus_ranking=r2_dr, delta_r2=r2_dr - r2_d,
                best_lambda=lam, rho_s_best_mixture=rho_m, mixture_gain=rho_d - rho_m)


def smcda_beyond_distance(sm, tta, origins, cost_dist, cell_km, reference=None, n_boot=None,
                          block_km=None, seed=None, logger_=None):
    """Round 2, Reviewer 2 F1: what does each ranking add to distance from the entry routes
    in ordering the simulated arrival? Over the domain cells with an arrival month:
      rho_s                      Spearman rho with arrival (negative = higher-ranked reached earlier)
      partial_rho_given_distance the same with Euclidean distance from the reference origin held
                                 fixed: 0 = nothing beyond distance; negative = adds the right order;
                                 positive = pulls the order the wrong way
      delta_r2                   gain in R² of rank(arrival) when rank(ranking) joins rank(distance)
      best_lambda, mixture_gain  the mixture (1-λ)·distance + λ·ranking (percentile ranks) with the
                                 most negative rho; gain = rho(distance) - rho(mixture) ≥ 0, an
                                 optimistic in-sample bound (λ = 0, gain 0: distance alone is best)
    CIs: spatial block bootstrap over blocks of block_km, λ chosen again in every replicate.
    The reference distance is the benchmark of Supplementary Table S7."""
    lg = logger_ or logger
    reference = reference or Config.BEYOND_DISTANCE_REFERENCE
    n_boot = int(Config.BEYOND_DISTANCE_BOOT if n_boot is None else n_boot)
    block_km = float(Config.BEYOND_DISTANCE_BLOCK_KM if block_km is None else block_km)
    seed = int(Config.BEYOND_DISTANCE_SEED if seed is None else seed)
    if not origins:
        return []
    ref = reference if reference in origins else next(iter(origins))
    mask = np.asarray(sm["mask"], bool)
    tta = np.asarray(tta, dtype=np.float64)
    dist = distance_transform_edt(~np.asarray(origins[ref], bool)) * float(cell_km)
    sel = mask & np.isfinite(tta) & np.isfinite(dist)
    if int(sel.sum()) < 30:
        lg.info(f"    (beyond-distance analysis skipped: {int(sel.sum())} cells with an arrival month)")
        return []
    std = dict(sm["standardised"])
    w0 = dict(sm["weights"])
    rk = [("composite (AHP, as reported)", "composite", np.asarray(sm["score"], np.float64))]
    if sm.get("classes") is not None:                 # the mapped classification itself (ordinal, ties)
        rk.append((f"risk class ({int(sm.get('n_classes', 5))} classes, as mapped)", "classes",
                   np.asarray(sm["classes"], np.float64)))
    entry_c = next((c for c in sm.get("criteria", []) if c.note == "_SMCDA_GATE_DIST"), None)
    if entry_c is not None and entry_c.key in w0:
        wl = {k: v for k, v in w0.items() if k != entry_c.key}
        s_ = sum(wl.values())
        if s_ > 0:
            rk.append((f"composite without {entry_c.key} (other AHP weights renormalised)", "composite",
                       _smcda_aggregate_like(sm, std, {k: v / s_ for k, v in wl.items()})))
    for k in std:
        rk.append((f"criterion {k} alone", "criterion", np.asarray(std[k], np.float64)))
    for oname, cd in (cost_dist or {}).items():
        if cd is not None:
            rk.append((f"cost distance from {oname}", "distance", -np.asarray(cd, np.float64)))
    rr_, cc_ = np.nonzero(sel)
    y = tta[sel]
    d = dist[sel]
    vals = [(nm, kd, np.asarray(a, np.float64)[sel]) for nm, kd, a in rk]
    vals = [(nm, kd, v) for nm, kd, v in vals if np.isfinite(v).all() and np.ptp(v) > 0]
    rows = [dict(ranking=f"Euclidean distance from {ref} (reference)", kind="reference",
                 n_cells=int(y.size), reference=ref, rho_s=spearman_ties(-d, y),
                 rho_s_reference=spearman_ties(-d, y), r2_reference=rank_r2(y, d))]
    for nm, kd, v in vals:
        rows.append(dict(ranking=nm, kind=kd, n_cells=int(y.size), reference=ref, **_beyond_stats(y, d, v)))
    cap = int(getattr(Config, "BEYOND_DISTANCE_BOOT_MAX_CELLS", 20000))
    if n_boot > 0 and y.size > cap:
        lg.info(f"    (beyond distance: {y.size:,} cells > BEYOND_DISTANCE_BOOT_MAX_CELLS = {cap:,}; "
                f"point values without bootstrap intervals)")
        for row in rows:
            row["bootstrap"] = f"not run ({y.size} cells > {cap})"
        n_boot = 0
    if n_boot > 0:
        b = max(1, int(round(block_km / float(cell_km))))
        blk = (rr_ // b).astype(np.int64) * 1_000_003 + (cc_ // b)
        _u, inv = np.unique(blk, return_inverse=True)
        members = [np.flatnonzero(inv == j) for j in range(_u.size)]
        rng = np.random.default_rng(seed)
        boot = {nm: [] for nm, _k, _v in vals}
        for _ in range(n_boot):
            idx = np.concatenate([members[j] for j in rng.integers(0, len(members), len(members))])
            if np.unique(y[idx]).size < 3:
                continue
            for nm, _kd, v in vals:
                st = _beyond_stats(y[idx], d[idx], v[idx])
                boot[nm].append((st["partial_rho_given_distance"], st["delta_r2"], st["mixture_gain"]))
        for row in rows[1:]:
            a = np.asarray(boot.get(row["ranking"]) or [], dtype=np.float64)
            if a.ndim == 2 and a.shape[0] >= 20:
                for j, key in enumerate(("partial_rho_given_distance", "delta_r2", "mixture_gain")):
                    col = a[:, j][np.isfinite(a[:, j])]
                    if col.size >= 20:
                        row[key + "_lo"] = float(np.percentile(col, 2.5))
                        row[key + "_hi"] = float(np.percentile(col, 97.5))
        for row in rows:
            row["bootstrap"] = f"{n_boot} spatial blocks of {b} x {b} cells ({b * float(cell_km):g} km), seed {seed}"
    return rows


def smcda_class_table(sm, month_offset):
    """Table 2: class, score interval, cells, share (sums to 100.0), mean arrival
    month on the global index, reached share."""
    v = sm.get("validation") or {}
    cls_rows = v.get("classes") or []
    counts = [int(r.get("cells", 0)) for r in cls_rows]
    shares = largest_remainder_shares(counts, 1)
    edges = list(sm.get("edges", []))
    out = []
    for i, r in enumerate(cls_rows):
        am = r.get("arr_mean")
        out.append(dict(
            cls=r.get("cls"), label=r.get("label"),
            score_lower=(edges[i] if i < len(edges) else None),
            score_upper=(edges[i + 1] if i + 1 < len(edges) else None),
            cells=r.get("cells"), share_pct=shares[i] if i < len(shares) else None,
            mean_arrival_month=(None if am is None else float(am) + int(month_offset)),
            median_arrival_month=(None if r.get("arr_median") is None
                                  else float(r["arr_median"]) + int(month_offset)),
            reached_pct=(None if r.get("reached_frac") is None else 100.0 * float(r["reached_frac"])),
            mean_attack_pct=(None if r.get("attack_mean") is None else 100.0 * float(r["attack_mean"]))))
    return out


# ══════════════════════════════════════════════════════════════════════════════
# ITEM 7 — GEOTIFF EXPORT
# ══════════════════════════════════════════════════════════════════════════════
def write_geotiff(path, arr, gt, wkt, nodata=-9999.0, dtype="float32"):
    """Write a single-band GeoTIFF on the model grid (NaN -> nodata)."""
    a = np.asarray(arr)
    gdt = {"float32": gdal.GDT_Float32, "int16": gdal.GDT_Int16, "byte": gdal.GDT_Byte}[dtype]
    if dtype == "float32":
        a = np.where(np.isfinite(a), a, nodata).astype(np.float32)
    else:
        a = np.where(np.isfinite(np.asarray(a, dtype=np.float64)), a, nodata).astype(
            np.int16 if dtype == "int16" else np.uint8)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ds = gdal.GetDriverByName("GTiff").Create(str(path), a.shape[1], a.shape[0], 1, gdt,
                                              options=["COMPRESS=DEFLATE"])
    ds.SetGeoTransform(tuple(float(x) for x in gt))
    ds.SetProjection(wkt)
    b = ds.GetRasterBand(1)
    b.WriteArray(a)
    b.SetNoDataValue(float(nodata))
    b.FlushCache()
    ds = None
    return path


# ══════════════════════════════════════════════════════════════════════════════
# ITEMS 5–6 — PENINSULAR-LEG VARIANTS (scenarios and sensitivity analysis)
# ══════════════════════════════════════════════════════════════════════════════
# A "context" holds everything a peninsular-leg run needs (crop arrays, seeds,
# gates, parameters). It is built once after the baseline leg, saved to
# PPR_v9.10_korea_leg_context.npz, and can be re-loaded with
#   python PPR_model_v9_10.py --from-context <file> --batch scenarios
# so scenarios and the sensitivity analysis do not need the continental leg.
V910_CTX_ARRAYS = ("KP", "KP_NOMASK", "NK", "SK", "ROK_FRAC", "EXCL", "CONTACT_ROK",
                   "CONTACT_DPRK", "CONTACT_SEAL", "MDL", "FENCE", "GOR_EFF", "GOR_nn", "GORANI_DENSITY",
                   "MAXENT", "CONDUCTANCE", "RESISTANCE", "RIV_nn", "WATER_ATTRACT",
                   "livestock_nn", "SCORE")


def save_korea_ctx(ctx, path):
    arrays, meta = {}, {}
    for k, v in ctx.items():
        if k in V910_CTX_ARRAYS:
            if v is not None:
                arrays[f"a__{k}"] = np.asarray(v)
        elif k == "params":
            pm = {}
            for pk, pv in v.items():
                if isinstance(pv, np.ndarray):
                    arrays[f"p__{pk}"] = pv
                elif isinstance(pv, np.generic):
                    pm[pk] = pv.item()
                elif pv is None or isinstance(pv, (bool, int, float, str)):
                    pm[pk] = pv
            meta["params"] = pm
        else:
            meta[k] = v
    arrays["meta"] = np.array(_v910_json.dumps(meta, default=lambda o: o.tolist()
                                               if isinstance(o, (np.ndarray, np.generic)) else str(o)))
    np.savez_compressed(path, **arrays)
    return Path(path)


def load_korea_ctx(path):
    z = np.load(path, allow_pickle=False)
    meta = _v910_json.loads(str(z["meta"]))
    ctx = dict(meta)
    params = dict(meta.get("params", {}))
    for k in z.files:
        if k.startswith("a__"):
            ctx[k[3:]] = z[k]
        elif k.startswith("p__"):
            params[k[3:]] = z[k]
    ctx["params"] = params
    for k in V910_CTX_ARRAYS:
        ctx.setdefault(k, None)
    for k in ("KP", "KP_NOMASK", "NK", "SK", "EXCL", "CONTACT_ROK", "CONTACT_DPRK", "CONTACT_SEAL", "MDL"):
        if ctx.get(k) is not None:
            ctx[k] = np.asarray(ctx[k], bool)
    ctx["seeds"] = [tuple(s) for s in ctx.get("seeds", [])]
    ctx["yalu_seeds"] = [tuple(s) for s in ctx.get("yalu_seeds", [])]
    ctx["shape"] = tuple(ctx["shape"])
    return ctx


def v910_gor_eff(ctx, psi, kp):
    """GOR_EFF rebuilt on the crop with habitat-blend weight psi (same rule as main())."""
    h = ctx["habitat"]
    g = np.nan_to_num(np.asarray(ctx["GOR_nn"], dtype=np.float64))
    if h.get("use_composite"):
        g = np.maximum(g, float(psi) * np.nan_to_num(np.asarray(ctx["GORANI_DENSITY"], dtype=np.float64)))
    mx = ctx.get("MAXENT")
    if h.get("use_maxent") and mx is not None:
        mx = np.asarray(mx, dtype=np.float64)
        have = np.isfinite(mx) & np.asarray(kp, bool)
        if h.get("base_is_maxent"):
            g = np.where(have, np.maximum(g, np.nan_to_num(mx)), g)
        else:
            b = float(h.get("maxent_blend", 0.5))
            g = np.where(have, np.clip((1.0 - b) * g + b * np.nan_to_num(mx), 0.0, 1.0), g)
    return g


def build_korea_variant(ctx, spec):
    """(env, RADParams, a_border, threshold, corridors, notes) for one variant.
    spec keys: factors{...} (sensitivity factors, see Config.SA_PARAM_RANGES);
    entry='yalu'; fence='uniform_phimax'|'uniform_mean'|'sealed';
    estuary='open'|'sealed'; gamma_mult, beta_mult; seasonality='off';
    corridors='equal'|'x0.5'|'x1.5'|'drop:<name>'|'add_inje'; exclaves='mask'|'keep'."""
    spec = dict(spec or {})
    P = {}
    for k, v in ctx["params"].items():
        P[k] = np.array(v, dtype=np.float64) if isinstance(v, (np.ndarray, list)) else v
    P.update(vacc_enable=False, log_every=0, store_states=False,
             vacc_trigger_mask=None, vacc_cap_field=None)
    H, W = tuple(ctx["shape"])
    notes = []
    kp = np.asarray(ctx["KP"], bool).copy()
    exc = np.asarray(ctx["EXCL"], bool) if ctx.get("EXCL") is not None else np.zeros((H, W), bool)
    if spec.get("exclaves") == "mask":
        kp &= ~exc
    elif spec.get("exclaves") == "keep" and ctx.get("KP_NOMASK") is not None:
        kp |= exc & np.asarray(ctx["KP_NOMASK"], bool)
    env = {k: np.nan_to_num(np.asarray(ctx[k], dtype=np.float64)) for k in
           ("GOR_EFF", "GOR_nn", "CONDUCTANCE", "RESISTANCE", "RIV_nn", "WATER_ATTRACT",
            "livestock_nn", "GORANI_DENSITY")}
    seeds = [tuple(s) for s in ctx["seeds"]]
    a_b = int(ctx["a_border"])
    corridors = [dict(c) for c in ctx["corridors"]]
    fk = dict(ctx["fence_kw"])
    rebuild = False
    thr = float(ctx["thr"])
    beta = np.asarray(P["beta_monthly"], dtype=np.float64).copy()
    tv = np.asarray(P["tvir_monthly"], dtype=np.float64).copy()
    f = dict(spec.get("factors") or {})
    # ---- sensitivity factors ----------------------------------------------------
    if "gamma" in f:
        P["gamma"] = float(f["gamma"])
    if "beta_level" in f or "beta_amp" in f:
        bm = float(beta.mean())
        beta = np.clip(float(f.get("beta_level", 1.0)) * (bm + float(f.get("beta_amp", 1.0)) * (beta - bm)),
                       0.0, None)
    if "tv_amp" in f:
        tm = float(tv.mean())
        tv = np.clip(tm + float(f["tv_amp"]) * (tv - tm), 0.0, None)
    if "diff_max_x" in f:
        P["diff_max"] = float(P["diff_max"]) * float(f["diff_max_x"])
    if "adv_max_x" in f:
        P["adv_max"] = float(P["adv_max"]) * float(f["adv_max_x"])
    if "adv_gain_x" in f:
        P["adv_gain"] = float(P["adv_gain"]) * float(f["adv_gain_x"])
    if "kappa" in f:
        env["CONDUCTANCE"] = np.where(kp, np.exp(-float(f["kappa"]) * env["RESISTANCE"]), 1e-9)
    if "phi0" in f:
        fk["sealed"] = float(f["phi0"])
        rebuild = True
    if "phimax_x" in f:
        fk["max_pass"] = float(fk["max_pass"]) * float(f["phimax_x"])
        rebuild = True
    if "rg_km" in f:
        fk["gate_radius_km"] = float(f["rg_km"])
        rebuild = True
    if "f0" in f:
        P["corridor_season_floor"] = float(f["f0"])
    if "psi" in f:
        env["GOR_EFF"] = v910_gor_eff(ctx, float(f["psi"]), kp)
    if "w_h" in f:
        P["trans_hab_w"] = float(f["w_h"])
        P["trans_lstk_w"] = 1.0 - float(f["w_h"])
    if "a_star" in f:
        thr = float(f["a_star"])
        P["arrival_threshold"] = thr
    if "i0" in f:
        P["seed_intensity"] = float(f["i0"])
    # ---- scenario switches --------------------------------------------------------
    if spec.get("gamma_mult"):
        P["gamma"] = float(P["gamma"]) * float(spec["gamma_mult"])
    if spec.get("beta_mult"):
        beta = beta * float(spec["beta_mult"])
    if spec.get("seasonality") == "off":
        disp = np.asarray(P["disp_km_monthly"], dtype=np.float64)
        P["disp_scale_ref"] = float(disp.max())
        P["disp_km_monthly"] = np.full(disp.size, float(disp.mean()))
        bt = float((beta * tv).mean())
        tv = np.full(tv.size, float(tv.mean()))
        beta = np.full(beta.size, bt / max(float(tv.mean()), 1e-12))
        notes.append("β·T_v and dispersal flat at their annual means")
    if spec.get("entry") == "yalu":
        if not ctx.get("yalu_seeds"):
            raise RuntimeError("no Yalu seed cells in the context")
        seeds = [tuple(s) for s in ctx["yalu_seeds"]]
        a_b = int(ctx["a_border_yalu"])
        P["season_offset"] = a_b                 # the leg starts on the Yalu introduction calendar
        notes.append(ctx.get("yalu_note", ""))
    cm = str(spec.get("corridors") or "")
    if cm:
        if cm == "equal":
            wbar = float(np.mean([c["weight"] for c in corridors])) if corridors else 1.0
            for c in corridors:
                c["weight"] = wbar
            notes.append(f"all gate weights = {wbar:.3f} (the mean, so total permeability is kept)")
        elif cm.startswith("x"):
            fac = float(cm[1:])
            for c in corridors:
                c["weight"] = float(c["weight"]) * fac
        elif cm.startswith("drop:"):
            nm = cm.split(":", 1)[1]
            corridors = [c for c in corridors if str(c["name"]) != nm]
        elif cm == "add_inje":
            if not ctx.get("inje_gate"):
                raise RuntimeError("the Inje–Seorak gate could not be placed on the MDL")
            corridors.append(dict(ctx["inje_gate"]))
        rebuild = True
    if spec.get("dispersal") == "flat":                          # round 2: no water-deer seasonality
        disp = np.asarray(P["disp_km_monthly"], dtype=np.float64)
        P["disp_scale_ref"] = float(disp.max())
        P["disp_km_monthly"] = np.full(disp.size, float(disp.mean()))
        notes.append("dispersal radius (and gate opening) flat at the annual mean; β·T_v seasonal")
    hs = spec.get("host")
    if hs in ("neutral", "livestock"):                             # round 2: comparator hosts
        def _dm(a):
            a = np.nan_to_num(np.asarray(a, dtype=np.float64))
            return float(a[kp].mean()) if kp.any() else 0.0
        if hs == "neutral":
            for _k in ("GOR_EFF", "GOR_nn", "GORANI_DENSITY", "WATER_ATTRACT"):
                env[_k] = np.where(kp, _dm(env[_k]), 0.0)
            notes.append("habitat-neutral host: habitat, density and water attraction at their domain means")
        else:
            # same VALUES as the water-deer fields, re-arranged in the order of small-ruminant density
            # (quantile mapping): only the spatial arrangement of the host changes, not its abundance
            _idx = np.flatnonzero(kp.ravel())
            _Lv = np.nan_to_num(np.asarray(env["livestock_nn"], dtype=np.float64)).ravel()[_idx]
            _tie = np.random.default_rng(20261119).random(_idx.size)       # random order among equal densities
            _order = _idx[np.lexsort((_tie, _Lv))]                           # ascending small-ruminant density
            for _k in ("GOR_EFF", "GOR_nn", "GORANI_DENSITY"):
                _v = np.sort(np.nan_to_num(np.asarray(env[_k], dtype=np.float64)).ravel()[_idx])
                _f = np.zeros(kp.size, dtype=np.float64)
                _f[_order] = _v
                env[_k] = np.where(kp, _f.reshape(kp.shape), 0.0)
            env["WATER_ATTRACT"] = np.where(kp, _dm(env["WATER_ATTRACT"]), 0.0)
            notes.append("host distributed like small ruminants: the water-deer habitat and density values "
                         "re-arranged in the order of small-ruminant density (same values, different places)")
    P["beta_monthly"] = beta
    P["tvir_monthly"] = tv
    # ---- fence field ----------------------------------------------------------------
    mdl = np.asarray(ctx["MDL"], bool)
    fm = spec.get("fence")
    if fm == "uniform_phimax":
        ff = np.where(mdl, float(fk["max_pass"]), 1.0)
    elif fm == "uniform_mean":
        ff = np.where(mdl, float(ctx["fence_mean"]), 1.0)
    elif fm == "sealed":
        ff = np.where(mdl, float(fk["sealed"]), 1.0)
    elif rebuild:
        ff = build_mdl_fence_field(corridors, mdl, (H, W), float(ctx["cell_km"]),
                                   gate_radius_km=float(fk["gate_radius_km"]),
                                   max_pass=float(fk["max_pass"]), sealed=float(fk["sealed"]),
                                   slowdown=float(fk["slowdown"]), logger=None)
        if ff is None:
            ff = np.where(mdl, float(fk["sealed"]), 1.0)
    else:
        ff = None if ctx.get("FENCE") is None else np.asarray(ctx["FENCE"], dtype=np.float64)
    mdl_v = mdl
    est = spec.get("estuary", ctx.get("estuary_policy", "open"))
    if est == "sealed" and (ctx.get("CONTACT_SEAL") is not None or ctx.get("CONTACT_ROK") is not None):
        if ctx.get("CONTACT_SEAL") is not None:
            contact = np.asarray(ctx["CONTACT_SEAL"], bool) & kp
        else:
            contact = (np.asarray(ctx["CONTACT_ROK"], bool) | np.asarray(ctx["CONTACT_DPRK"], bool)) & kp
        mdl_v, ff = seal_contact_band(mdl, ff, contact, sealed=float(fk["sealed"]),
                                      default_phi=float(ctx["han_phi"]))
    env.update(KP=kp, NK=np.asarray(ctx["NK"], bool) & kp, SK=np.asarray(ctx["SK"], bool) & kp,
               MDL_BOUNDARY=mdl_v, MDL_FENCE_FIELD=ff, DMZ_ENTRIES=seeds, NK_ENTRIES=[])
    if spec.get("vacc") == "V2":                                   # round 2: V2 inside a variant
        _rok = np.asarray(ctx["SK"], bool) & kp
        _exc = np.asarray(ctx["EXCL"], bool) if ctx.get("EXCL") is not None else np.zeros_like(_rok)
        if spec.get("host") == "livestock" or spec.get("vacc_cap") == "all":
            _share = np.where(kp, 1.0, 0.0)                         # every host vaccinable
        else:
            _hab = np.clip(np.nan_to_num(env["GOR_EFF"]), 0.0, 1.0)
            _lst = np.clip(np.nan_to_num(env["livestock_nn"]), 0.0, 1.0)
            _wh, _wl = float(P["trans_hab_w"]), float(P["trans_lstk_w"])
            _share = np.where(kp, (_wl * _lst) / (_wh * _hab + _wl * _lst + 1e-9), 0.0)
        env["VACC_FIELD"] = _rok.astype(np.float32)
        P.update(vacc_enable=True, vacc_mode=Config.VACC_MODE, vacc_start_month=0,
                 vacc_coverage=float(Config.VACC_COVERAGE), vacc_rate=float(Config.VACC_RATE),
                 vacc_efficacy=float(Config.VACC_EFFICACY),
                 vacc_ring_cells=max(1, int(round(Config.VACC_RING_KM / max(float(ctx["cell_km"]), 1e-6)))),
                 vacc_inf_thr=float(Config.VACC_INFECTED_THRESHOLD),
                 vacc_trigger_mask=(_rok & ~_exc), vacc_cap_field=_share)
        if spec.get("vacc_start") == "fence":                      # FC1: start at the DPRK side of the fence
            _band = np.asarray(ctx["MDL"], bool)
            _nk = np.asarray(ctx["NK"], bool) & kp
            P["vacc_trigger_mask"] = _nk & binary_dilation(_band, structure=V910_FOUR) & ~_band
            notes.append("vaccination V2 triggered when the DPRK cells next to the fence band reach A*")
        else:
            notes.append("vaccination V2 (ROK small ruminants, triggered by the first ROK cell)")
        if spec.get("vacc_cap") == "all":
            notes.append("V2 coverage not capped by the small-ruminant share (every host vaccinable)")
    fields = set(RADParams.__dataclass_fields__)
    p = RADParams(**{k: v for k, v in P.items() if k in fields})
    return env, p, a_b, thr, corridors, notes


def variant_metrics(cum, ctx, env, a_b, thr, corridors):
    """Table-4 / sensitivity outputs of one peninsular-leg run."""
    kp = np.asarray(env["KP"], bool)
    nk = np.asarray(env["NK"], bool)
    rok = np.asarray(env["SK"], bool)
    exc = (np.asarray(ctx["EXCL"], bool) if ctx.get("EXCL") is not None
           else np.zeros(kp.shape, bool)) & kp
    rok_main = rok & ~exc
    w = np.where(kp, np.nan_to_num(np.asarray(ctx["ROK_FRAC"], dtype=np.float64)), 0.0)
    rep = int(ctx["report_month"])
    met, _S = korea_leg_metrics(cum, nk, rok, w, thr, a_b, float(ctx["est_mean"]),
                                float(ctx["est_area"]), [rep], first_exclude=exc)
    tta = first_arrival_index(cum, thr, mask=kp)
    loc = first_crossing_location(tta, rok_main, corridors,
                                  ctx.get("CONTACT_ROK"), float(ctx["cell_km"]),
                                  float(ctx["probe_km"]))
    band = np.asarray(ctx["MDL"], bool)          # the fence band itself, never the sealed contact
    ring_n = nk & binary_dilation(band, structure=V910_FOUR) & ~band
    tf = tta[ring_n] if ring_n.any() else np.array([])
    tf = tf[np.isfinite(tf)]
    fence_first = (int(tf.min()) + 1 + int(a_b)) if tf.size else None
    tb = tta[rok_main & ~band]
    tb = tb[np.isfinite(tb)]
    beyond = (int(tb.min()) + 1 + int(a_b)) if tb.size else None
    loc_b = first_crossing_location(tta, rok_main & ~band, corridors, ctx.get("CONTACT_ROK"),
                                    float(ctx["cell_km"]), float(ctx["probe_km"]))
    first_rok = met["first_rok_cell"]
    intro = int(a_b) + 1
    rho, n_rho = float("nan"), 0
    sc = ctx.get("SCORE")
    if sc is not None:
        sc = np.asarray(sc, dtype=np.float64)
        sel = rok_main & np.isfinite(tta) & np.isfinite(sc) & ((tta + 1 + a_b) <= rep)
        n_rho = int(sel.sum())
        if n_rho >= 10:
            rho = spearman_ties(sc[sel], tta[sel])
    at_rep = met["at"].get(rep) or {}
    return dict(
        intro_month=intro, end_month=met["end_month"], threshold=float(thr),
        first_rok=first_rok, first_rok_any=met["first_rok_any"],
        first_location=loc["label"], first_location_km=loc["dist_km"],
        first_location_all=loc["labels"],
        est=met["est_areaw"], est_cell=met["est_cell"],
        area_1=met["area_1"], area_5=met["area_5"], area_10=met["area_10"],
        mean_5=met["mean_5"], nk_5=met["nk_5"],
        fence_dprk_first=fence_first, first_rok_beyond_band=beyond,
        first_beyond_location=loc_b["label"], first_beyond_location_all=loc_b["labels"],
        wait_at_fence=(None if (beyond is None or fence_first is None) else beyond - fence_first),
        transit_months=(None if first_rok is None else first_rok - intro),
        transit_years=(None if first_rok is None else (first_rok - intro) / 12.0),
        pct_rok_reached_report=at_rep.get("rok_area_w_pct"),
        mean_A_rok_report=at_rep.get("rok_mean_w_pct"),
        mean_A_nk_report=at_rep.get("nk_mean_pct"),
        rho_s=rho, n_rho=n_rho)


def run_korea_variant(ctx, job):
    """Run one peninsular-leg variant to max(report month, ROK establishment),
    stopping early only if the epidemic is extinct. Never raises."""
    t0 = _v910_time.time()
    out = dict(id=job.get("id"), label=job.get("label"), group=job.get("group", ""), ok=False)
    try:
        env, p, a_b, thr, corridors, notes = build_korea_variant(ctx, job.get("spec") or {})
        rep = int(ctx["report_month"])
        n_min = max(12, rep - int(a_b))
        w = np.where(env["KP"], np.nan_to_num(np.asarray(ctx["ROK_FRAC"], dtype=np.float64)), 0.0)
        ws = float(w.sum())
        hist = []
        em, ea = float(ctx["est_mean"]), float(ctx["est_area"])
        thr_stop = max([thr] + [float(str(a).split("_", 1)[1]) for a in (job.get("post_astar") or [])])

        def _stop(_m, cum, _w=w, _ws=ws, _thr=thr_stop):
            c = np.nan_to_num(cum)
            hist.append(float(c.sum()))
            if _ws > 0 and ((float((_w * c).sum()) / _ws >= em)
                            or (float((_w * (c >= _thr)).sum()) / _ws >= ea)):
                return True
            # extinct or stalled: the total attack grew by < 0.01% over two years
            return len(hist) > 36 and abs(hist[-1] - hist[-25]) <= 1e-4 * max(1.0, hist[-1])

        res = run_dynamic_rad_simulation(env, p, n_months=n_min, closure=None, stop_fn=_stop,
                                         stop_extra_months=0, max_months=int(ctx["cap_months"]),
                                         stop_label="variant: ROK establishment or extinction")
        cum = res["monthly_cum"]
        out.update(variant_metrics(cum, ctx, env, a_b, thr, corridors))
        for lab in job.get("post_astar") or []:            # e.g. 'astar_0.10' (post hoc A*)
            out[str(lab)] = variant_metrics(cum, ctx, env, a_b, float(str(lab).split("_", 1)[1]),
                                            corridors)
        _von = res.get("vacc_on_month")
        out["vacc_start_month"] = None if _von is None else int(_von) + 1 + int(a_b)   # global month
        out.update(ok=True, n_months=len(cum), notes="; ".join(n for n in notes if n))
    except Exception as e:
        out.update(error=f"{type(e).__name__}: {e}")
    out["seconds"] = round(_v910_time.time() - t0, 2)
    return out


_V910_POOL_CTX = None


def _v910_pool_worker(job):
    return run_korea_variant(_V910_POOL_CTX, job)


def run_variants(ctx, jobs, workers=1, label="batch", logger_=None):
    """Run jobs (serially, or in forked worker processes on Linux/WSL)."""
    global _V910_POOL_CTX
    lg = logger_ or logger
    n = len(jobs)
    results = [None] * n
    t0 = _v910_time.time()
    step = max(1, n // 20)

    def _progress(i):
        if (i + 1) % step == 0 or i + 1 == n:
            el = _v910_time.time() - t0
            lg.info(f"   [{label}] {i + 1}/{n} runs done ({el / 60:.1f} min, "
                    f"~{el / (i + 1) * (n - i - 1) / 60:.1f} min left)")
    if int(workers) > 1 and n > 1:
        try:
            import multiprocessing as _mp
            _V910_POOL_CTX = ctx
            with _mp.get_context("fork").Pool(processes=int(workers)) as pool:
                for i, r in enumerate(pool.imap(_v910_pool_worker, jobs, chunksize=1)):
                    results[i] = r
                    _progress(i)
            return results
        except Exception as e:
            lg.warning(f"   ⚠ parallel run failed ({type(e).__name__}: {e}); running serially.")
        finally:
            _V910_POOL_CTX = None
    for i, job in enumerate(jobs):
        if results[i] is None:
            results[i] = run_korea_variant(ctx, job)
            _progress(i)
    return results


def v910_workers(setting):
    if isinstance(setting, int) and setting > 0:
        return setting
    try:
        return max(1, min(8, (os.cpu_count() or 2) - 1))
    except Exception:
        return 1


# ── sampling designs and PRCC ──────────────────────────────────────────────────
def morris_design(k, r, levels=4, seed=0):
    """Morris (1991) one-at-a-time trajectories on the unit cube.
    Returns X (r*(k+1), k), meta [(trajectory, step, factor moved)], delta."""
    rng = np.random.default_rng(int(seed))
    delta = levels / (2.0 * (levels - 1))
    grid = np.linspace(0.0, 1.0, int(levels))
    X, meta = [], []
    for t in range(int(r)):
        x = rng.choice(grid, size=int(k))
        X.append(x.copy())
        meta.append((t, 0, -1))
        for s, i in enumerate(rng.permutation(int(k)), 1):
            x[i] = x[i] + delta if x[i] + delta <= 1.0 + 1e-12 else x[i] - delta
            X.append(x.copy())
            meta.append((t, s, int(i)))
    return np.asarray(X), meta, delta


def morris_effects(X, Y, meta, k):
    """mu, mu* and sigma of the elementary effects of each factor (unit-cube scale)."""
    EE = [[] for _ in range(int(k))]
    Y = np.asarray(Y, dtype=np.float64)
    for j in range(1, len(meta)):
        t, s, i = meta[j]
        if s == 0 or meta[j - 1][0] != t:
            continue
        dx = X[j, i] - X[j - 1, i]
        dy = Y[j] - Y[j - 1]
        if dx != 0 and np.isfinite(dy):
            EE[i].append(dy / dx)
    out = []
    for i in range(int(k)):
        e = np.asarray(EE[i], dtype=np.float64)
        out.append(dict(n=int(e.size),
                        mu=float(e.mean()) if e.size else float("nan"),
                        mu_star=float(np.abs(e).mean()) if e.size else float("nan"),
                        sigma=float(e.std(ddof=1)) if e.size > 1 else float("nan")))
    return out


def lhs_design(n, k, seed=0):
    """Latin hypercube sample on the unit cube (one stratum per run and factor)."""
    rng = np.random.default_rng(int(seed))
    U = np.empty((int(n), int(k)))
    for j in range(int(k)):
        U[:, j] = (rng.permutation(int(n)) + rng.random(int(n))) / float(n)
    return U


def prcc(X, y):
    """Partial rank correlation coefficients (Marino et al. 2008) with p-values."""
    X = np.asarray(X, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    ok = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
    X, y = X[ok], y[ok]
    n, k = X.shape
    out = []
    if n < k + 3:
        return [(float("nan"), float("nan"))] * k
    R = np.column_stack([v910_rank_avg(X[:, j]) for j in range(k)])
    ry = v910_rank_avg(y)
    try:
        from scipy.stats import t as _t_dist
    except Exception:
        _t_dist = None
    for j in range(k):
        Z = np.column_stack([np.ones(n), np.delete(R, j, axis=1)])
        bx = np.linalg.lstsq(Z, R[:, j], rcond=None)[0]
        by = np.linalg.lstsq(Z, ry, rcond=None)[0]
        ex = R[:, j] - Z @ bx
        ey = ry - Z @ by
        den = math.sqrt(float(ex @ ex) * float(ey @ ey))
        r = float(ex @ ey) / den if den > 0 else float("nan")
        dfree = n - 2 - (k - 1)
        p = float("nan")
        if _t_dist is not None and dfree > 0 and np.isfinite(r) and abs(r) < 1:
            tstat = r * math.sqrt(dfree / (1.0 - r * r))
            p = float(2.0 * _t_dist.sf(abs(tstat), dfree))
        out.append((r, p))
    return out


V910_SA_OUTPUTS = [
    ("first_rok", "First arrival in the ROK (month)"),
    ("est", "Establishment in the ROK (month)"),
    ("pct_rok_reached_report", "ROK area reached by the report month (%)"),
    ("transit_years", "DPRK transit (years)"),
    ("wait_at_fence", "Wait at the fence (months)"),
    ("rho_s", "ρs, risk score vs arrival"),
]


def _sa_censor(results):
    """One censoring month for every run in which an event never happened: the
    latest month any run reached + 1, so a rank-based method puts all of them
    last together (an extinct run stops early and must not rank as 'early')."""
    ends = [int(r["end_month"]) for r in results if r and r.get("ok") and r.get("end_month")]
    return (max(ends) + 1) if ends else None


def _sa_values(results, key, censor=None):
    """Output vector across runs (NaN for failed runs); events that never
    happened are set to the common censoring month."""
    if censor is None:
        censor = _sa_censor(results)
    out = []
    for r in results:
        if not r or not r.get("ok"):
            out.append(float("nan"))
            continue
        v = r.get(key)
        if v is None and censor is not None and key in ("first_rok", "est", "transit_years"):
            v = (float(censor) if key != "transit_years"
                 else (float(censor) - float(r.get("intro_month", 0))) / 12.0)
        elif v is None and censor is not None and key == "wait_at_fence" and r.get("fence_dprk_first"):
            v = float(censor) - float(r["fence_dprk_first"])     # reached the fence, never crossed
        out.append(float("nan") if v is None else float(v))
    return np.asarray(out, dtype=np.float64)


# ── scenario list (Table 4) ────────────────────────────────────────────────────
def v910_scenario_jobs(ctx, names):
    """Expand Config.SCENARIO_LIST into jobs. 'drop_each' -> one job per gate."""
    est_base = ctx.get("estuary_policy", "open")
    exc_base = ctx.get("exclave_policy", "report")
    desc = {
        "baseline": ("main", {}, "baseline peninsular leg (re-run with the batch horizon rule)"),
        "entry_yalu": ("main", {"entry": "yalu"},
                       "entry into the DPRK across the lower Yalu River"),
        "fence_uniform_phimax": ("main", {"fence": "uniform_phimax"},
                                 "uniformly permeable MDL: φ = φmax at every fence cell"),
        "fence_uniform_mean": ("extra", {"fence": "uniform_mean"},
                               "uniform MDL at the mean pass-through of the gated fence"),
        "fence_sealed_no_gates": ("extra", {"fence": "sealed"},
                                  "no gates: φ = φ0 at every fence cell"),
        "estuary_toggle": ("main", {"estuary": "open" if est_base == "sealed" else "sealed"},
                           ("barrier extended along the Han River estuary (contact sealed at φ0)"
                            if est_base != "sealed" else "estuary left open (baseline seals it)")),
        "gamma_x2": ("main", {"gamma_mult": 2.0, "beta_mult": 2.0},
                     "γ × 2 with β(m) × 2 (β/γ constant); DPRK introduction month fixed"),
        "gamma_x4": ("main", {"gamma_mult": 4.0, "beta_mult": 4.0},
                     "γ × 4 with β(m) × 4 (β/γ constant); DPRK introduction month fixed"),
        "no_seasonality": ("main", {"seasonality": "off"},
                           "no seasonal forcing (β·T_v, dispersal and gate opening at annual means)"),
        "corridors_equal": ("corridor", {"corridors": "equal"}, "equal gate weights"),
        "host_neutral": ("comparator", {"host": "neutral"},
                         "habitat-neutral wildlife host: habitat, density and water attraction at their domain "
                         "means; same dispersal kernel, barrier and seeding"),
        "host_livestock": ("comparator", {"host": "livestock"},
                           "host distributed like small ruminants (the water-deer habitat and density values "
                           "re-arranged in the order of small-ruminant density); same dispersal kernel, "
                           "barrier and seeding"),
        "dispersal_flat": ("comparator", {"dispersal": "flat"},
                           "no seasonal dispersal: radius and gate opening at the annual mean"),
        "baseline_V2": ("vaccination", {"vacc": "V2"}, "baseline host, vaccination scenario V2"),
        "host_neutral_V2": ("vaccination", {"host": "neutral", "vacc": "V2"},
                            "habitat-neutral host, vaccination scenario V2"),
        "host_livestock_V2": ("vaccination", {"host": "livestock", "vacc": "V2"},
                              "small-ruminant host, vaccination scenario V2 (every host vaccinable, cap κv)"),
        "baseline_V2_allhosts": ("vaccination", {"vacc": "V2", "vacc_cap": "all"},
                                 "baseline host, V2 with every host vaccinable (no small-ruminant cap)"),
        "baseline_V2_fence": ("vaccination", {"vacc": "V2", "vacc_start": "fence"},
                              "baseline host, V2 started when infection reaches the DPRK side of the fence"),
        "baseline_V2_allhosts_fence": ("vaccination", {"vacc": "V2", "vacc_cap": "all", "vacc_start": "fence"},
                                       "baseline host, V2 with every host vaccinable, started at the fence"),
        "corridors_x0.5": ("corridor", {"corridors": "x0.5"}, "gate weights × 0.5"),
        "corridors_x1.5": ("corridor", {"corridors": "x1.5"}, "gate weights × 1.5"),
        "add_inje": ("corridor", {"corridors": "add_inje"},
                     "Inje–Seorak gate added at an estimated position (weight 0.40)"),
        "exclaves_toggle": ("extra", {"exclaves": "mask" if exc_base != "mask" else "keep"},
                            ("ROK islands joined to the DPRK coast removed from the grid"
                             if exc_base != "mask" else "ROK islands joined to the DPRK coast kept")),
    }
    jobs = []
    for nm in names:
        nm = str(nm)
        if nm == "drop_each":
            for c in ctx["corridors"]:
                jobs.append(dict(label=f"drop_{c['name']}", group="corridor",
                                 description=f"gate {c['name']} removed",
                                 spec={"corridors": f"drop:{c['name']}"}))
            continue
        if nm.startswith("astar_"):
            continue                      # derived post hoc from the baseline run
        if nm not in desc:
            continue
        grp, spec, d = desc[nm]
        job = dict(label=nm, group=grp, description=d, spec=spec)
        if nm == "baseline":
            job["post_astar"] = [str(a) for a in names if str(a).startswith("astar_")]
        jobs.append(job)
    for i, j in enumerate(jobs):
        j["id"] = i
    return jobs


V910_T4_COLUMNS = [
    "scenario", "group", "description", "ok", "intro_month", "first_rok", "first_rok_any",
    "first_location", "first_location_km", "est", "est_cell", "transit_months", "transit_years",
    "fence_dprk_first", "first_rok_beyond_band", "wait_at_fence", "pct_rok_reached_report", "mean_A_rok_report",
    "mean_A_nk_report", "rho_s", "delta_first_rok", "delta_est", "end_month", "threshold",
    "n_months", "seconds", "notes", "error"]


def v910_targeting_rho_ci(prio, t_i, ok_t, idx, shape, cell_km, ref_name, n_boot, block_km, seed):
    """Round 2 (F1): spatial block-bootstrap 95% intervals for the Spearman rho of each surveillance
    ranking with the arrival month, and for its difference from the reference ranking (distance from the
    entry routes). Paired: every replicate resamples the same blocks of cells for all rankings, so the
    difference is computed within each replicate. Cells: the ROK cells reached (t_i finite). Blocks of
    block_km x block_km on the grid, as for Table S7b. Returns (rows, description)."""
    sel = np.flatnonzero(np.asarray(ok_t, bool))
    y = np.asarray(t_i, dtype=np.float64)[sel]
    rr, cc = np.unravel_index(np.asarray(idx)[sel], shape)
    b = max(1, int(round(float(block_km) / float(cell_km))))
    blk = (rr // b).astype(np.int64) * 1_000_003 + (cc // b)
    _u, inv = np.unique(blk, return_inverse=True)
    members = [np.flatnonzero(inv == j) for j in range(_u.size)]
    names = list(prio.keys())
    vals = {nm: np.asarray(prio[nm], dtype=np.float64)[sel] for nm in names}
    point = {nm: spearman_ties(vals[nm], y) for nm in names}
    rng = np.random.default_rng(int(seed))
    boot = {nm: [] for nm in names}
    used = 0
    for _ in range(int(n_boot)):
        k = np.concatenate([members[j] for j in rng.integers(0, len(members), len(members))])
        if np.unique(y[k]).size < 3:
            continue
        used += 1
        for nm in names:
            boot[nm].append(spearman_ties(vals[nm][k], y[k]))
    ref = np.asarray(boot[ref_name], dtype=np.float64)
    desc = (f"{used} of {int(n_boot)} replicates; spatial blocks of {b} x {b} cells "
            f"({b * float(cell_km):g} km); seed {int(seed)}")

    def _q(a, p):
        a = a[np.isfinite(a)]
        return f"{float(np.percentile(a, p)):.5f}" if a.size >= 20 else None

    rows = []
    for nm in names:
        a = np.asarray(boot[nm], dtype=np.float64)
        row = dict(ranking=nm, n_cells=int(y.size), n_blocks=int(len(members)),
                   rho_s_arrival=f"{point[nm]:.5f}", rho_s_lo=_q(a, 2.5), rho_s_hi=_q(a, 97.5),
                   bootstrap=desc)
        if nm == ref_name:
            row.update(diff_vs_distance="0.00000")
        else:
            dl = a - ref
            ok = np.isfinite(dl)
            row.update(diff_vs_distance=f"{point[nm] - point[ref_name]:.5f}", diff_lo=_q(dl, 2.5),
                       diff_hi=_q(dl, 97.5),
                       pct_replicates_closer_than_distance=(f"{100.0 * float((dl[ok] < 0).mean()):.1f}"
                                                            if ok.sum() >= 20 else None))
        rows.append(row)
    return rows, desc


def v910_run_round2(ctx, out_dir, workers=1, logger_=None):
    """Round 2 (19 Nov 2026): comparator hosts with and without V2, gate-zone radii,
    and the surveillance-targeting trade-off (distance vs risk score vs host density)."""
    lg = logger_ or logger
    written = []
    rep = int(ctx["report_month"])
    # ---- (a) comparators ------------------------------------------------------------------
    jobs = ([] if getattr(Config, "R2_SKIP_COMPARATORS", False)
            else v910_scenario_jobs(ctx, list(Config.COMPARATOR_LIST)))
    if jobs:
        lg.info(f"[ROUND 2 — comparators] {len(jobs)} runs: " + ", ".join(j["label"] for j in jobs))
    else:
        lg.info("[ROUND 2] comparator runs skipped (--batch targeting): zone radii and targeting only")
    res = run_variants(ctx, jobs, workers=workers, label="comparators", logger_=lg) if jobs else []
    cols = ["scenario", "group", "description", "ok", "first_rok", "first_location", "first_location_km",
            "first_rok_beyond_band", "first_beyond_location", "first_beyond_location_all", "est",
            "transit_years", "wait_at_fence", "pct_rok_reached_report", "mean_A_rok_report",
            "mean_A_nk_report", "rho_s", "n_months", "notes", "error"]
    if getattr(Config, "R2_VAX_DECOMP", False):
        cols.insert(cols.index("mean_A_rok_report"), "vacc_start_month")
    rows = []
    for j, r in zip(jobs, res):
        r = r or {}
        row = {c: r.get(c) for c in cols}
        row.update(scenario=j["label"], group=j["group"], description=j["description"])
        rows.append(row)
    if rows:
      written.append(v910_write_csv(
        out_dir / f"{V910_TAG}_{getattr(Config, 'R2_COMPARATOR_STEM', 'R2_comparators')}.csv", cols, rows,
        [f"{V910_TAG} round 2 — comparator hosts (same seeding, dispersal kernel and barrier as the baseline)",
         f"months on the global index; report month {rep}; mean_A_rok_report = area-weighted mean attack rate "
         f"in the ROK at the report month (%)",
         "first_beyond_location = gate (within the zone radius), 'estuary' or 'other' for the first ROK cell "
         "outside the fence band"], lg))
    for row in rows:
        lg.info(f"   {str(row['scenario']):<20s} first ROK {str(row['first_rok']):>4s} "
                f"({str(row['first_location'])[:16]}); beyond band {str(row['first_rok_beyond_band']):>4s} "
                f"({str(row['first_beyond_location'])[:16]}); est {str(row['est']):>4s}; "
                f"mean A ROK@{rep} {row['mean_A_rok_report']}")
    if getattr(Config, "R2_VAX_DECOMP", False):                     # --batch vaxdecomp: part (a) only
        _b = {str(r["scenario"]): r for r in rows}
        try:
            a0 = float(_b["baseline"]["mean_A_rok_report"])
            for row in rows:
                if str(row["scenario"]).startswith("baseline_V2") and row.get("mean_A_rok_report") is not None:
                    lg.info(f"   {str(row['scenario']):<28s} V2 start {row.get('vacc_start_month')}; "
                            f"mean A ROK@{rep} reduction vs baseline "
                            f"{100.0 * (a0 - float(row['mean_A_rok_report'])) / max(a0, 1e-9):.1f}%")
        except Exception as e:
            lg.warning(f"   (reduction summary not computed: {e})")
        return written
    # ---- (b) one baseline run to the conditional month: zone radii and targeting -----------
    try:
        env, p, a_b, thr, corridors, _n = build_korea_variant(ctx, {})
        n_m = max(12, int(Config.R2_CONDITIONAL_MONTH) - int(a_b))
        cum = run_dynamic_rad_simulation(env, p, n_months=n_m, closure=None)["monthly_cum"]
        kp = np.asarray(env["KP"], bool)
        tta = first_arrival_index(cum, thr, mask=kp)
        exc = (np.asarray(ctx["EXCL"], bool) if ctx.get("EXCL") is not None else np.zeros(kp.shape, bool)) & kp
        rok_main = np.asarray(env["SK"], bool) & ~exc
        band = np.asarray(ctx["MDL"], bool)
        nk = np.asarray(env["NK"], bool)
        cell_km = float(ctx["cell_km"])
        zrows = []
        for rad in Config.R2_ZONE_RADII_KM:
            for z in fence_wait_table(tta, corridors, rok_main, nk, band, cell_km, float(rad), a_b,
                                      contact_rok=ctx.get("CONTACT_ROK")):
                z = dict(z)
                z.pop("rok_first_cells", None)
                z["radius_km"] = float(rad)
                zrows.append(z)
        zcols = ["radius_km", "zone", "kind", "weight", "dprk_first_month", "band_first_month",
                 "rok_band_first_month", "rok_first_month", "rok_median_month", "wait_first_months",
                 "n_dprk_cells", "n_band_cells", "n_rok_cells"]
        written.append(v910_write_csv(out_dir / f"{V910_TAG}_R2_zone_radius.csv", zcols, zrows,
                                      [f"{V910_TAG} round 2 — fence wait by gate zone for several zone radii; "
                                       f"cell size {cell_km:g} km"], lg))
        # targeting trade-off on ROK territory (prospective domain)
        lst = np.clip(np.nan_to_num(np.asarray(ctx["livestock_nn"], dtype=np.float64)), 0.0, None)
        w = np.where(kp, np.nan_to_num(np.asarray(ctx["ROK_FRAC"], dtype=np.float64)), 0.0)
        cells = rok_main & kp
        t_glob = np.where(np.isfinite(tta), tta + 1 + int(a_b), np.nan)
        entry = np.zeros(kp.shape, bool)
        for co in corridors:
            entry[int(co["row"]), int(co["col"])] = True
        if ctx.get("CONTACT_ROK") is not None:
            entry |= np.asarray(ctx["CONTACT_ROK"], bool)
        dist = distance_transform_edt(~entry) * cell_km
        sc = ctx.get("SCORE")
        rankings = {"distance from the entry routes": -dist, "small-ruminant density (HOST)": lst}
        if sc is not None:
            rankings["risk score (AHP weights)"] = np.nan_to_num(np.asarray(sc, dtype=np.float64), nan=-1e9)
        idx = np.flatnonzero(cells.ravel())
        tot_l = float((lst * w).ravel()[idx].sum())
        orders, prio = {}, {}
        for name, key in rankings.items():
            kv = key.ravel()[idx]
            orders[name] = idx[np.argsort(-kv, kind="mergesort")]
            prio[name] = kv                                       # rho_s with ties averaged
        # round 2 (F1): two-step rankings — distance band first (nearest first), a second criterion within the
        # band (highest first), then distance; and a weighted mix of distance and risk score (percentile ranks)
        d_i = dist.ravel()[idx]
        second = {}
        if sc is not None:
            second["risk score"] = rankings["risk score (AHP weights)"].ravel()[idx]
        second["goat density (HOST)"] = lst.ravel()[idx]
        for bkm in getattr(Config, "R2_TWO_STEP_BANDS_KM", ()) or ():
            band = np.floor(d_i / float(bkm))
            for sname, sv in second.items():
                o = np.lexsort((d_i, -sv, band))
                nm = f"two-step: {float(bkm):g}-km distance bands, then {sname}"
                orders[nm] = idx[o]
                pos = np.empty(idx.size)
                pos[o] = np.arange(idx.size, dtype=np.float64)
                prio[nm] = -pos
        lam = float(getattr(Config, "R2_TWO_STEP_WEIGHT", 0.0) or 0.0)
        if sc is not None and 0.0 < lam < 1.0:
            mix = ((1.0 - lam) * v910_rank_avg(-d_i) + lam * v910_rank_avg(second["risk score"])) / d_i.size
            o = np.lexsort((d_i, -mix))
            nm = f"weighted: {1.0 - lam:.0%} distance + {lam:.0%} risk score"
            orders[nm] = idx[o]
            prio[nm] = mix
        t_i = t_glob.ravel()[idx]
        ok_t = np.isfinite(t_i)
        trows = []
        for name, order in orders.items():
            rho_arr = spearman_ties(prio[name][ok_t], t_i[ok_t]) if ok_t.sum() >= 3 else float("nan")
            for share in Config.R2_TARGET_SHARES:
                sel = order[:max(1, int(round(share * idx.size)))]
                tg = t_glob.ravel()[sel]
                reached = np.isfinite(tg)
                trows.append(dict(
                    ranking=name, share_of_rok_cells=share, n_cells=int(sel.size),
                    first_detection_month=(int(np.nanmin(tg)) if reached.any() else None),
                    pct_reached_by_conditional_month=100.0 * float(reached.mean()),
                    median_arrival_month=(float(np.median(tg[reached])) if reached.any() else None),
                    pct_small_ruminant_density=(100.0 * float((lst * w).ravel()[sel].sum()) / tot_l
                                                if tot_l > 0 else None),
                    rho_s_arrival=rho_arr))
        written.append(v910_write_csv(
            out_dir / f"{V910_TAG}_R2_targeting_tradeoff.csv",
            ["ranking", "share_of_rok_cells", "n_cells", "first_detection_month",
             "pct_reached_by_conditional_month", "median_arrival_month", "pct_small_ruminant_density",
             "rho_s_arrival"],
            trows, [f"{V910_TAG} round 2 — cells placed under surveillance by each ranking (ROK territory); "
                    f"arrival from the baseline run to month {int(Config.R2_CONDITIONAL_MONTH)}",
                    "first_detection_month = earliest month in which any selected cell reached A*; "
                    "pct_small_ruminant_density = share of the ROK's small-ruminant density layer "
                    "(area-weighted) in the selected cells",
                    "two-step rankings (round 2, F1): cells ordered by band of distance from the entry routes "
                    "(nearest band first), then by the second criterion within the band (highest first), then by "
                    "distance; weighted = (1-w) x percentile rank of minus distance + w x percentile rank of the "
                    "risk score",
                    "rho_s_arrival = Spearman rho between the ranking (priority) and the arrival month over the "
                    "ROK cells reached (ties of a ranking value averaged; two-step rankings by position)"], lg))
        for r in trows:
            if r["share_of_rok_cells"] in (0.10, 0.20):
                lg.info(f"   targeting {r['share_of_rok_cells']:.0%} by {r['ranking']}: first detection "
                        f"{r['first_detection_month']}, median arrival {r['median_arrival_month']}, "
                        f"{r['pct_reached_by_conditional_month']:.1f}% reached, "
                        f"{r['pct_small_ruminant_density']:.1f}% of small-ruminant density; "
                        f"rho_s {r['rho_s_arrival']:+.3f}")
        # round 2 (F1): how sure is the gain of each ranking over distance alone? (paired block bootstrap)
        try:
            ref_nm = "distance from the entry routes"
            n_bt = int(getattr(Config, "R2_TARGET_BOOT", 0) or 0)
            if n_bt > 0 and int(ok_t.sum()) >= 30 and ref_nm in prio:
                ci_rows, ci_desc = v910_targeting_rho_ci(
                    prio, t_i, ok_t, idx, kp.shape, cell_km, ref_nm, n_bt,
                    float(Config.BEYOND_DISTANCE_BLOCK_KM), int(getattr(Config, "R2_TARGET_BOOT_SEED", 20261125)))
                written.append(v910_write_csv(
                    out_dir / f"{V910_TAG}_R2_targeting_rho_ci.csv",
                    ["ranking", "n_cells", "n_blocks", "rho_s_arrival", "rho_s_lo", "rho_s_hi", "diff_vs_distance",
                     "diff_lo", "diff_hi", "pct_replicates_closer_than_distance", "bootstrap"],
                    ci_rows, [f"{V910_TAG} round 2 (F1) — Spearman rho of each surveillance ranking with the arrival "
                              f"month over the ROK cells reached by month {int(Config.R2_CONDITIONAL_MONTH)} "
                              f"(outside the masked exclaves), with spatial block-bootstrap 95% intervals ({ci_desc})",
                              "diff_vs_distance = rho of the ranking minus rho of distance from the entry routes; "
                              "negative = the ranking orders arrival more closely than distance alone; paired "
                              "bootstrap (the same resampled cells for every ranking); "
                              "pct_replicates_closer_than_distance = share of replicates with a difference below 0"],
                    lg))
                for r in ci_rows:
                    msg = (f"   rho_s {r['ranking']}: {r['rho_s_arrival']} "
                           f"(95% CI {r.get('rho_s_lo')} to {r.get('rho_s_hi')})")
                    if r["ranking"] != ref_nm:
                        msg += (f"; minus distance {r['diff_vs_distance']} (95% CI {r.get('diff_lo')} to "
                                f"{r.get('diff_hi')}); closer than distance in "
                                f"{r.get('pct_replicates_closer_than_distance')}% of replicates")
                    lg.info(msg)
            if getattr(Config, "R2_TARGET_CELL_TABLE", False):
                rr_, cc_ = np.unravel_index(idx, kp.shape)
                sc_i = (np.asarray(sc, dtype=np.float64).ravel()[idx] if sc is not None
                        else np.full(idx.size, np.nan))
                l_i = lst.ravel()[idx]
                w_i = w.ravel()[idx]
                crow = [[int(rr_[j]), int(cc_[j]),                     # full precision (repr), so that
                         (str(int(t_i[j])) if np.isfinite(t_i[j]) else ""),   # every ranking can be rebuilt
                         repr(float(d_i[j])),
                         (repr(float(sc_i[j])) if np.isfinite(sc_i[j]) else ""),
                         repr(float(l_i[j])), repr(float(w_i[j]))] for j in range(idx.size)]
                written.append(v910_write_csv(
                    out_dir / f"{V910_TAG}_R2_targeting_cells.csv",
                    ["row", "col", "arrival_month", "distance_km", "risk_score", "goat_density", "rok_share"],
                    crow, [f"{V910_TAG} round 2 (F1) — ROK cells of the targeting comparison (Table S9); "
                           f"cell size {cell_km:g} km; arrival month on the global index, blank = not reached by "
                           f"month {int(Config.R2_CONDITIONAL_MONTH)}",
                           "distance_km = Euclidean distance from the nearest gate or ROK estuary contact cell; "
                           "goat_density = small-ruminant (goat) density layer of the model (head per km2); "
                           "rok_share = ROK share of the cell area (weight of the goat totals)"], lg))
        except Exception as e:
            lg.warning(f"   ⚠ round-2 targeting intervals failed ({type(e).__name__}: {e})", exc_info=True)
    except Exception as e:
        lg.warning(f"   ⚠ round-2 zone/targeting analysis failed ({type(e).__name__}: {e})", exc_info=True)
    return written


def v910_scenario_table(ctx, results, jobs, extra_rows=()):
    rows = []
    base = next((r for r in results if r and r.get("label") == "baseline" and r.get("ok")), None)
    jd = {j["id"]: j for j in jobs}
    for r in results:
        j = jd.get(r.get("id"), {})
        row = dict(r)
        row.update(scenario=r.get("label"), group=j.get("group", ""),
                   description=j.get("description", ""))
        rows.append(row)
        if r.get("label") == "baseline":
            for k, v in r.items():
                if k.startswith("astar_") and isinstance(v, dict):
                    rr = dict(v)
                    rr.update(scenario=k, group="main", ok=True,
                              description=f"arrival threshold A* = {k.split('_', 1)[1]} "
                                          f"(post hoc on the baseline run)")
                    rows.append(rr)
    rows.extend(extra_rows)
    if base:
        for row in rows:
            for key, dk in (("first_rok", "delta_first_rok"), ("est", "delta_est")):
                if row.get(key) is not None and base.get(key) is not None:
                    row[dk] = int(row[key]) - int(base[key])
    return rows


V910_SYM_TEX = {"w_h": "$w_h$", "D_max": "$D_{\\max}$", "T_v": "$T_v$", "A*": "$A^*$",
                "v_max": "$v_{\\max}$", "a_g": "$a_g$", "φ0": "φ$_0$", "φmax": "φ$_{\\max}$",
                "r_g": "$r_g$", "f0": "$f_0$", "I0": "$I_0$", "β amp.": "β amplitude"}


def v910_sym(s):
    """Parameter symbol as typeset in the manuscript (matplotlib mathtext)."""
    return V910_SYM_TEX.get(str(s), str(s))


def render_figure10(prcc_rows, lhs_vals, base_first, out_path, report_month, censor=None):
    """Figure 10: (A) PRCC per parameter for three outputs; (B) distribution of
    the month of first arrival in the ROK across the LHS runs. No figure title."""
    outs = [("first_rok", "First arrival (month)", "#1f77b4"),
            ("est", "Establishment (month)", "#ff7f0e"),
            ("pct_rok_reached_report", f"ROK reached by month {report_month} (%)", "#2ca02c")]
    params = []
    for r in prcc_rows:
        if r["parameter"] not in params:
            params.append(r["parameter"])
    look = {(r["parameter"], r["output"]): r for r in prcc_rows}
    fig, ax = plt.subplots(1, 2, figsize=(14, max(4.8, 0.42 * len(params) + 1.6)),
                           gridspec_kw=dict(width_ratios=[1.35, 1.0]), facecolor="white")
    y = np.arange(len(params))
    h = 0.8 / len(outs)
    for i, (key, lab, col) in enumerate(outs):
        vals = [look.get((p, key), {}).get("prcc", np.nan) for p in params]
        sig = [look.get((p, key), {}).get("p_value", np.nan) for p in params]
        bars = ax[0].barh(y + (i - 1) * h, np.nan_to_num(vals), height=h, color=col, label=lab)
        for b, s in zip(bars, sig):
            if np.isfinite(s) and s < 0.05:
                b.set_edgecolor("black")
                b.set_linewidth(0.8)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels([v910_sym(look.get((p, outs[0][0]), {}).get("symbol", p)) for p in params],
                          fontsize=11)
    ax[0].invert_yaxis()
    ax[0].axvline(0, color="#444", lw=0.8)
    ax[0].set_xlim(-1, 1)
    ax[0].set_xlabel("PRCC", fontsize=11)
    ax[0].legend(fontsize=9, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False)
    ax[0].text(-0.08, 1.01, "A", transform=ax[0].transAxes, fontsize=14, fontweight="bold")
    v_all = np.asarray(lhs_vals, dtype=np.float64)
    v_all = v_all[np.isfinite(v_all)]
    cens = (v_all >= float(censor)) if censor is not None else np.zeros(v_all.size, bool)
    v = v_all[~cens]
    if v.size:
        step = max(1.0, float(np.ceil((v.max() - v.min()) / 30.0)))
        bins = np.arange(np.floor(v.min()) - 0.5, np.ceil(v.max()) + 0.5 + step, step)
        ax[1].hist(v, bins=bins, color="#6baed6", edgecolor="white")
    if v_all.size:
        q5, q50, q95 = np.percentile(v if v.size else v_all, [5, 50, 95])   # over the runs that reached the ROK

        def _m(q):
            return f"≥{float(censor):.0f} (not reached)" if (censor is not None and q >= float(censor)) else f"{q:.0f}"
        ax[1].axvline(min(q50, v.max() if v.size else q50), color="#08519c", lw=1.6, label=f"median {_m(q50)}")
        ax[1].axvspan(q5, min(q95, (v.max() + 0.5) if v.size else q95), color="#c6dbef", alpha=0.45,
                      zorder=0, label=f"90% range {_m(q5)}–{_m(q95)}")
    if cens.any():
        ax[1].text(0.98, 0.70, f"{int(cens.sum())} of {v_all.size} runs did not reach the ROK by month {float(censor):.0f}",
                   transform=ax[1].transAxes, ha="right", fontsize=9)
    if base_first is not None:
        ax[1].axvline(float(base_first), color="#d62728", lw=1.6, ls="--",
                      label=f"baseline {int(base_first)}")
    ax[1].set_xlabel("Month of first arrival in the ROK", fontsize=11)
    ax[1].set_ylabel("Runs", fontsize=11)
    ax[1].legend(fontsize=9, frameon=False)
    ax[1].text(-0.1, 1.01, "B", transform=ax[1].transAxes, fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=600, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


def render_morris_figure(morris_rows, out_path):
    """Supplementary figure: Morris mu* vs sigma for the first-arrival month."""
    rows = [r for r in morris_rows if r["output"] == "first_rok" and np.isfinite(r.get("mu_star", np.nan))]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(7.5, 6), facecolor="white")
    _cut = 0.05 * max(float(r["mu_star"]) for r in rows)      # label only the influential parameters
    _small = []
    for r in rows:
        ax.scatter(r["mu_star"], np.nan_to_num(r["sigma"]), s=36, color="#1f77b4")
        if float(r["mu_star"]) >= _cut:
            _left = float(r["mu_star"]) > 17.0 * _cut           # labels of the right-most points go left
            ax.annotate(v910_sym(r.get("symbol", r["parameter"])),
                        (r["mu_star"], np.nan_to_num(r["sigma"])),
                        xytext=((-6, 4) if _left else (5, 3)), textcoords="offset points",
                        ha=("right" if _left else "left"), fontsize=11)
        else:
            _small.append(v910_sym(r.get("symbol", r["parameter"])))
    if _small:
        ax.text(0.97, 0.03, f"μ* < {_cut:.0f}: " + ", ".join(_small), transform=ax.transAxes,
                ha="right", va="bottom", fontsize=10, wrap=True,
                bbox=dict(facecolor="white", edgecolor="#bbbbbb", boxstyle="round,pad=0.3"))
    ax.set_xlabel("μ* (months per unit of the scaled range)", fontsize=11)
    ax.set_ylabel("σ", fontsize=11)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=600, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return Path(out_path)


def v910_run_batch(ctx, scenarios=True, sa_mode="off", out_dir=None, dl_dir=None,
                   workers=1, logger_=None):
    """Items 5–6: Table 4 scenarios and the global sensitivity analysis on the
    peninsular leg. Every file is written as soon as its runs finish."""
    lg = logger_ or logger
    out_dir = Path(out_dir or Config.OUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    rep = int(ctx["report_month"])
    nw = v910_workers(workers)
    lg.info("\n" + "=" * 96)
    lg.info(f"[v9.10 BATCH] peninsular leg only; DPRK introduction month {int(ctx['a_border']) + 1}; "
            f"report month {rep}; {nw} worker(s); seed {Config.SA_SEED}")
    lg.info("=" * 96)
    if getattr(Config, "RUN_R2", False):
        written += v910_run_round2(ctx, out_dir, workers=nw, logger_=lg)
    # ---- Table 4 -----------------------------------------------------------------
    if scenarios:
        jobs = v910_scenario_jobs(ctx, list(Config.SCENARIO_LIST))
        lg.info(f"[TABLE 4] {len(jobs)} scenario runs: " + ", ".join(j["label"] for j in jobs))
        res = run_variants(ctx, jobs, workers=nw, label="scenarios", logger_=lg)
        extra = []
        g4 = Config.GRID4KM_SUMMARY
        if g4 is None:
            g4 = out_dir.parent / (out_dir.name + "_grid4km") / f"{V910_TAG}_run_summary.json"
        try:
            if g4 and Path(g4).exists():
                js = _v910_json.loads(Path(g4).read_text(encoding="utf-8"))
                b = dict(js.get("baseline", {}))
                b.update(scenario="grid_4km", group="main", ok=True,
                         description=f"4-km grid (separate run: {Path(g4).parent.name})")
                extra.append(b)
            else:
                lg.info(f"   (no 4-km run found at {g4}; run with --preset grid4km to add that row)")
        except Exception as e:
            lg.warning(f"   ⚠ could not read the 4-km run summary: {e}")
        rows = v910_scenario_table(ctx, res, jobs, extra)
        meta = [f"{V910_TAG} Table 4 — peninsular-leg scenarios",
                f"DPRK introduction month {int(ctx['a_border']) + 1} (except entry_yalu: "
                f"{int(ctx.get('a_border_yalu', ctx['a_border'])) + 1}); A* = {ctx['thr']}; "
                f"report month {rep}",
                "months on the global index (month 1 = January of Year 1); ROK values area-weighted",
                "first_rok = first cell assigned to the ROK with A >= A* (fence band included, exclave "
                "islands excluded); est = mean A >= 5% or 10% of ROK area reached (area-weighted)",
                "wait_at_fence = first ROK cell beyond the fence band − first DPRK cell next to the band",
                "each run lasts to max(report month, establishment); empty = not reached"]
        written.append(v910_write_csv(out_dir / f"{V910_TAG}_Table4_scenarios.csv",
                                      V910_T4_COLUMNS, rows, meta, lg))
        for r in rows:
            lg.info(f"   {str(r.get('scenario')):<26s} first ROK {str(r.get('first_rok')):>4s} "
                    f"({str(r.get('first_location'))[:18]:<18s}) est {str(r.get('est')):>4s}  "
                    f"transit {('%.1f' % r['transit_years']) if r.get('transit_years') is not None else '—':>5s} yr  "
                    f"ROK reached@{rep} "
                    f"{('%.1f' % r['pct_rok_reached_report']) if r.get('pct_rok_reached_report') is not None else '—':>5s}%"
                    f"{('  ERROR ' + r['error']) if r.get('error') else ''}")
    # ---- sensitivity analysis ---------------------------------------------------------
    mode = str(sa_mode or "off").lower()
    if mode in ("off", "none", ""):
        v910_copy_to_dl(written, dl_dir, lg)
        return written
    ranges = dict(Config.SA_PARAM_RANGES)
    names = list(ranges.keys())
    base_f = dict(ctx["sa_baseline"])
    sym = {k: ranges[k][2] for k in names}
    written.append(v910_write_csv(
        out_dir / f"{V910_TAG}_TableS3_parameter_ranges.csv",
        ["parameter", "symbol", "low", "high", "baseline", "description"],
        [dict(parameter=k, symbol=ranges[k][2], low=ranges[k][0], high=ranges[k][1],
              baseline=base_f.get(k), description=ranges[k][3]) for k in names],
        [f"{V910_TAG} Supplementary Table S3 — ranges for the sensitivity analysis"], lg))

    def unit_to_val(k, u):
        lo, hi = float(ranges[k][0]), float(ranges[k][1])
        if k in tuple(getattr(Config, "SA_LOG_PARAMS", ()) or ()) and lo > 0 and hi > 0:
            return float(lo * (hi / lo) ** float(u))              # round 2: log-uniform
        return lo + float(u) * (hi - lo)

    def jobs_from(U, keys, tag):
        jobs = []
        for i, u in enumerate(U):
            fac = dict(base_f)
            for j, k in enumerate(keys):
                fac[k] = unit_to_val(k, u[j])
            jobs.append(dict(id=i, label=f"{tag}_{i:04d}", group=tag, spec={"factors": fac}))
        return jobs

    def run_table(jobs, res, keys, fname, extra_cols=()):
        hdr = (["run", "ok"] + list(keys) + [o for o, _ in V910_SA_OUTPUTS]
               + ["first_location", "first_rok_beyond_band", "first_beyond_location", "mean_A_rok_report",
                  "intro_month", "end_month", "seconds", "error"] + list(extra_cols))
        rows = []
        for j, r in zip(jobs, res):
            row = dict(run=j["label"], ok=r.get("ok"))
            for k in keys:
                row[k] = j["spec"]["factors"][k]
            for o, _ in V910_SA_OUTPUTS:
                row[o] = r.get(o)
            for c in ("first_location", "first_rok_beyond_band", "first_beyond_location", "mean_A_rok_report",
                      "intro_month", "end_month", "seconds", "error"):
                row[c] = r.get(c)
            rows.append(row)
        return v910_write_csv(out_dir / fname, hdr, rows,
                              [f"{V910_TAG} sensitivity runs ({len(rows)}), seed {Config.SA_SEED}"], lg)

    sel_keys = list(Config.SA_LHS_PARAMS) if Config.SA_LHS_PARAMS else None
    if mode == "oat":
        jobs, k_id = [dict(id=0, label="oat_base", group="oat", spec={"factors": dict(base_f)})], 1
        for k in names:
            for side, u in (("low", 0.0), ("high", 1.0)):
                fac = dict(base_f)
                fac[k] = unit_to_val(k, u)
                jobs.append(dict(id=k_id, label=f"oat_{k}_{side}", group="oat",
                                 spec={"factors": fac}, param=k, side=side))
                k_id += 1
        lg.info(f"[SENSITIVITY — one-at-a-time] {len(jobs)} runs (each parameter at both ends of its range)")
        res = run_variants(ctx, jobs, workers=nw, label="OAT", logger_=lg)
        base_r = res[0]
        rows = []
        for j, r in zip(jobs[1:], res[1:]):
            row = dict(parameter=j["param"], symbol=sym[j["param"]], side=j["side"],
                       value=j["spec"]["factors"][j["param"]], ok=r.get("ok"), error=r.get("error"))
            for o, _ in V910_SA_OUTPUTS:
                row[o] = r.get(o)
                if r.get(o) is not None and base_r.get(o) is not None:
                    row[f"delta_{o}"] = float(r[o]) - float(base_r[o])
            rows.append(row)
        hdr = (["parameter", "symbol", "side", "value", "ok"] + [o for o, _ in V910_SA_OUTPUTS]
               + [f"delta_{o}" for o, _ in V910_SA_OUTPUTS] + ["error"])
        written.append(v910_write_csv(out_dir / f"{V910_TAG}_SA_one_at_a_time.csv", hdr, rows,
                                      [f"{V910_TAG} one-at-a-time sensitivity; baseline row = "
                                       f"first ROK {base_r.get('first_rok')}, est {base_r.get('est')}"], lg))
        v910_copy_to_dl(written, dl_dir, lg)
        return written
    if "morris" in mode:
        k = len(names)
        X, meta, delta = morris_design(k, int(Config.SA_MORRIS_TRAJECTORIES),
                                       int(Config.SA_MORRIS_LEVELS), int(Config.SA_SEED))
        jobs = jobs_from(X, names, "morris")
        lg.info(f"[SENSITIVITY — Morris screening] {k} parameters × "
                f"{int(Config.SA_MORRIS_TRAJECTORIES)} trajectories = {len(jobs)} runs")
        res = run_variants(ctx, jobs, workers=nw, label="Morris", logger_=lg)
        written.append(run_table(jobs, res, names, f"{V910_TAG}_SA_morris_runs.csv"))
        mrows = []
        score = np.zeros(k)
        cens_m = _sa_censor(res)
        for o, olab in V910_SA_OUTPUTS:
            y = _sa_values(res, o, cens_m)
            eff = morris_effects(X, y, meta, k)
            ms = np.array([e["mu_star"] for e in eff], dtype=np.float64)
            mx = np.nanmax(ms) if np.isfinite(ms).any() else np.nan
            if o in ("first_rok", "est", "pct_rok_reached_report") and np.isfinite(mx) and mx > 0:
                score = np.maximum(score, np.nan_to_num(ms / mx))
            for i, e in enumerate(eff):
                mrows.append(dict(output=o, output_label=olab, parameter=names[i], symbol=sym[names[i]],
                                  mu=e["mu"], mu_star=e["mu_star"], sigma=e["sigma"], n_effects=e["n"],
                                  mu_star_scaled=(e["mu_star"] / mx if np.isfinite(mx) and mx > 0 else None)))
        written.append(v910_write_csv(out_dir / f"{V910_TAG}_SA_morris.csv",
                                      ["output", "output_label", "parameter", "symbol", "mu", "mu_star",
                                       "sigma", "n_effects", "mu_star_scaled"], mrows,
                                      [f"{V910_TAG} Morris elementary effects (unit-cube scale; "
                                       f"{int(Config.SA_MORRIS_LEVELS)} levels, Δ = {delta:.3f})"], lg))
        try:
            fp = render_morris_figure(mrows, out_dir / f"{V910_TAG}_FigS_morris.png")
            if fp:
                written.append(fp)
        except Exception as e:
            lg.warning(f"   ⚠ Morris figure skipped: {e}")
        order = list(np.argsort(-score))
        if sel_keys is None:
            sel_keys = [names[i] for i in order[:int(Config.SA_LHS_TOP_K)]]
        lg.info("   Morris ranking (max scaled μ* over first arrival, establishment, % reached): "
                + ", ".join(f"{sym[names[i]]}={score[i]:.2f}" for i in order))
    if "lhs" in mode:
        keys = sel_keys or names
        U = lhs_design(int(Config.SA_LHS_N), len(keys), int(Config.SA_SEED) + 1)
        jobs = jobs_from(U, keys, "lhs")
        lg.info(f"[SENSITIVITY — LHS + PRCC] {len(keys)} parameters ({', '.join(sym[k] for k in keys)}) "
                f"× {len(jobs)} runs; the others fixed at baseline")
        res = run_variants(ctx, jobs, workers=nw, label="LHS", logger_=lg)
        written.append(run_table(jobs, res, keys, f"{V910_TAG}_SA_lhs_runs.csv"))
        Xv = np.array([[j["spec"]["factors"][k] for k in keys] for j in jobs], dtype=np.float64)
        prows, urows = [], []
        cens = _sa_censor(res)
        for o, olab in V910_SA_OUTPUTS:
            y = _sa_values(res, o, cens)
            pc = prcc(Xv, y)
            for k, (r_, p_) in zip(keys, pc):
                prows.append(dict(output=o, output_label=olab, parameter=k, symbol=sym[k],
                                  prcc=r_, p_value=p_))
            fin = y[np.isfinite(y)]
            urows.append(dict(output=o, output_label=olab, n_runs=int(fin.size),
                               median=(float(np.median(fin)) if fin.size else None),
                               p05=(float(np.percentile(fin, 5)) if fin.size else None),
                               p95=(float(np.percentile(fin, 95)) if fin.size else None),
                               minimum=(float(fin.min()) if fin.size else None),
                               maximum=(float(fin.max()) if fin.size else None)))
        locs = [r.get("first_location") for r in res if r and r.get("ok")]
        for lab in sorted(set(locs)):
            urows.append(dict(output="first_location", output_label=f"first crossing at: {lab}",
                              n_runs=len(locs), median=100.0 * locs.count(lab) / max(len(locs), 1)))
        n_cens = sum(1 for r in res if r and r.get("ok") and r.get("first_rok") is None)
        meta = [f"{V910_TAG} PRCC (Marino et al. 2008), {len(jobs)} LHS runs, seed {int(Config.SA_SEED) + 1}",
                f"{n_cens} run(s) never reached the ROK: their months are set to {cens} (latest end + 1), "
                f"so they rank last together"]
        written.append(v910_write_csv(out_dir / f"{V910_TAG}_SA_prcc.csv",
                                      ["output", "output_label", "parameter", "symbol", "prcc", "p_value"],
                                      prows, meta, lg))
        arows = []                                                  # round 2: arrived runs only
        for o, olab in V910_SA_OUTPUTS:
            yy = np.array([np.nan if (not r or not r.get("ok") or r.get(o) is None) else float(r[o])
                           for r in res], dtype=np.float64)
            sel = np.isfinite(yy)
            if sel.sum() >= len(keys) + 3:
                for k, (r_, p_) in zip(keys, prcc(Xv[sel], yy[sel])):
                    arows.append(dict(output=o, output_label=olab, parameter=k, symbol=sym[k],
                                      prcc=r_, p_value=p_, n_runs=int(sel.sum())))
        written.append(v910_write_csv(out_dir / f"{V910_TAG}_SA_prcc_arrived.csv",
                                      ["output", "output_label", "parameter", "symbol", "prcc", "p_value",
                                       "n_runs"], arows,
                                      [f"{V910_TAG} PRCC on the runs in which the event happened (no censoring)"],
                                      lg))
        written.append(v910_write_csv(out_dir / f"{V910_TAG}_SA_uncertainty.csv",
                                      ["output", "output_label", "n_runs", "median", "p05", "p95",
                                       "minimum", "maximum"], urows,
                                      meta + ["first_location rows: 'median' holds the % of runs"], lg))
        try:
            base_first = (ctx.get("baseline_metrics") or {}).get("first_rok")
            written.append(render_figure10(prows, _sa_values(res, "first_rok", cens), base_first,
                                           out_dir / f"{V910_TAG}_Fig10_sensitivity.png", rep, cens))
        except Exception as e:
            lg.warning(f"   ⚠ Figure 10 skipped: {e}")
        for o, olab in V910_SA_OUTPUTS[:3]:
            top = sorted([r for r in prows if r["output"] == o and np.isfinite(r["prcc"])],
                         key=lambda r: -abs(r["prcc"]))[:4]
            lg.info(f"   {olab}: " + ", ".join(f"{r['symbol']} {r['prcc']:+.2f}" for r in top))
    v910_copy_to_dl(written, dl_dir, lg)
    return written


# ── item 7: temporal hold-out — seed only from the outbreaks of a dated table ─────────
# The run of record seeds from every point of the outbreak shapefile (all EMPRES-i
# outbreaks in China 2007–2025 plus Mongolia), so its maps cannot be tested against
# those outbreaks. A hold-out run seeds from the outbreaks reported before a cut-off
# date and leaves the later ones for ppr_external_test.py.
V910_SEED_INFO = {}
V910_SEED_TABLE_NAMES = ("PPR_seed_points_dated.csv",)
V910_SCRIPT_DIR = Path(globals().get("__file__", "PPR_model_v9_10.py")).resolve().parent
_V910_LON_NAMES = ("longitude", "lon", "long", "lng", "decimallongitude")
_V910_LAT_NAMES = ("latitude", "lat", "decimallatitude")
_V910_DATE_NAMES = ("date", "observation.date..dd.mm.yyyy.", "report.date..dd.mm.yyyy.")


def v910_parse_date(s):
    """A date from a table cell: ISO (2015-03-30) or day first (30/03/2015, as in
    EMPRES-i exports); Excel / pandas dates too. None when empty or unreadable."""
    import datetime as _dt
    if s is None:
        return None
    if isinstance(s, _dt.datetime):
        return s.date()
    if isinstance(s, _dt.date):
        return s
    if hasattr(s, "to_pydatetime"):
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


def v910_read_table(path):
    """(column names, row dicts) from a CSV or Excel file. The header is the first line
    holding both a longitude and a latitude column, so the metadata lines that EMPRES-i
    exports put above the header are skipped."""
    path = Path(path)

    def _hdr(lines):
        for i, cells in enumerate(lines):
            if cells and str(cells[0]).lstrip().startswith("#"):
                continue                                    # provenance lines of the v9.10 CSVs
            low = {str(c).strip().lower() for c in cells}
            if low & set(_V910_LON_NAMES) and low & set(_V910_LAT_NAMES):
                return i
        return 0
    if path.suffix.lower() in (".xlsx", ".xls"):
        import pandas as pd
        raw = pd.read_excel(path, header=None, dtype=object)
        h = _hdr([list(raw.iloc[i]) for i in range(min(len(raw), 60))])
        cols = [str(v).strip() for v in raw.iloc[h]]
        return cols, [dict(zip(cols, list(raw.iloc[i]))) for i in range(h + 1, len(raw))]
    text = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    h = _hdr([next(csv.reader([ln]), []) for ln in text[:60]])
    rd = csv.DictReader(text[h:])
    rows = list(rd)
    return list(rd.fieldnames or []), rows


def v910_seed_rows(path, date_col="date", before=None, undated="drop"):
    """Rows of a seed table kept for a hold-out run: dated before `before` (all rows
    when `before` is None). `date_col` may list several columns ('a|b'); the first
    non-empty one counts. Returns (rows, info)."""
    cols, rows = v910_read_table(path)
    low = {str(c).strip().lower(): c for c in cols}
    cx = next((low[a] for a in _V910_LON_NAMES if a in low), None)
    cy = next((low[a] for a in _V910_LAT_NAMES if a in low), None)
    if cx is None or cy is None:
        raise SystemExit(f"v9.10 seed table {path}: no longitude / latitude columns among {cols}")
    dcols = [low[c.strip().lower()] for c in str(date_col or "").split("|") if c.strip().lower() in low]
    if not dcols:
        dcols = [low[a] for a in _V910_DATE_NAMES if a in low]
    cut = v910_parse_date(before) if before else None
    if before and cut is None:
        raise SystemExit(f"v9.10: SEED_BEFORE '{before}' is not a date (use YYYY-MM-DD)")
    if cut is not None and not dcols:
        raise SystemExit(f"v9.10 seed table {path}: a cut-off date was given but the table has no "
                         f"date column (looked for '{date_col}', {_V910_DATE_NAMES})")
    cid = next((low[a] for a in ("seed_id", "event_id", "event.id", "id") if a in low), None)
    keep, n_bad, n_late, n_und_drop, n_und_keep = [], 0, 0, 0, 0
    for i, r in enumerate(rows):
        try:
            lon, lat = float(r[cx]), float(r[cy])
        except (KeyError, TypeError, ValueError):
            n_bad += 1
            continue
        if not (np.isfinite(lon) and np.isfinite(lat)) or abs(lat) > 90 or abs(lon) > 180:
            n_bad += 1
            continue
        d = None
        for c_ in dcols:
            d = v910_parse_date(r.get(c_))
            if d is not None:
                break
        if cut is not None:
            if d is None:
                if str(undated).lower() != "keep":
                    n_und_drop += 1
                    continue
                n_und_keep += 1
            elif d >= cut:
                n_late += 1
                continue
        keep.append(dict(lon=lon, lat=lat, date=d, id=str(r.get(cid, "") if cid else i + 1)))
    info = dict(source="table", table=str(path), date_columns=dcols,
                before=(cut.isoformat() if cut else None), undated=str(undated),
                table_rows=len(rows), kept_rows=len(keep), dropped_on_or_after_cutoff=n_late,
                undated_dropped=n_und_drop, undated_kept=n_und_keep, unreadable_rows=n_bad)
    ds_ = [k["date"] for k in keep if k["date"] is not None]
    info["kept_date_range"] = [min(ds_).isoformat(), max(ds_).isoformat()] if ds_ else None
    return keep, info


def v910_seed_entries(rows, wkt, xy_to_rc, nr, nc, weight, info, logger_=None, label="seed_table"):
    """RAD seed entries (row, col, label, 0, weight) for table rows (WGS84 lon/lat),
    with the same grid rule as the outbreak shapefile. Stops the run when nothing is
    left: a hold-out run must never fall back to the full seed file."""
    src = _osr.SpatialReference()
    src.ImportFromEPSG(4326)
    dst = _osr.SpatialReference()
    dst.ImportFromWkt(wkt)
    for s_ in (src, dst):
        try:
            s_.SetAxisMappingStrategy(_osr.OAMS_TRADITIONAL_GIS_ORDER)
        except AttributeError:
            pass
    tf = _osr.CoordinateTransformation(src, dst)
    entries, items, outside = [], [], 0
    for r in rows:
        x, y = tf.TransformPoint(float(r["lon"]), float(r["lat"]))[:2]
        r_, c_ = xy_to_rc(x, y)
        if 0 <= r_ < nr and 0 <= c_ < nc:
            entries.append((int(r_), int(c_), label, 0, float(weight)))
            items.append(dict(x=x, y=y, lon=r["lon"], lat=r["lat"], row=int(r_), col=int(c_),
                              date=r.get("date"), id=r.get("id", "")))
        else:
            outside += 1
    info.update(seed_rows=len(entries), seed_cells=len({(e[0], e[1]) for e in entries}),
                outside_grid=outside)
    V910_SEED_INFO.clear()
    V910_SEED_INFO.update(info)
    V910_SEED_INFO["_items"] = items
    if logger_ is not None:
        logger_.info(f"  ✓ v9.10 hold-out seeding from {Path(str(info.get('table'))).name}: "
                     f"{len(entries)} rows in {info['seed_cells']} cells"
                     + (f", reported before {info['before']}" if info.get("before") else "")
                     + f" ({info.get('dropped_on_or_after_cutoff', 0)} later and "
                       f"{info.get('undated_dropped', 0)} undated row(s) left out; {outside} outside the grid)")
    if not entries:
        raise SystemExit("v9.10: the seed table left no seed points inside the grid — check SEED_TABLE / "
                         "SEED_BEFORE. The run stops instead of seeding from the full outbreak file.")
    return entries


def v910_find_seed_table(setting, roots):
    """Config.SEED_TABLE: a path, or 'auto' = PPR_seed_points_dated.csv in one of `roots`."""
    if setting and str(setting).lower() != "auto":
        p = Path(str(setting))
        if not p.exists():
            raise SystemExit(f"v9.10: seed table not found: {p}")
        return p
    for root in roots:
        for nm in V910_SEED_TABLE_NAMES:
            p = Path(str(root)) / nm
            if p.exists():
                return p
    raise SystemExit("v9.10: hold-out seeding needs PPR_seed_points_dated.csv (in the folder you run "
                     "from, next to the model script, or in the data folder), or --seed-table PATH")


def v910_write_seed_list(out_dir, wkt, logger_=None):
    """<OUT_DIR>/PPR_v9.10_seed_points_used.csv + _seeding.json: the points that seeded
    this run. ppr_external_test.py reads them to keep the test set independent."""
    items = V910_SEED_INFO.get("_items") or []
    if not items:
        return []
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tf = None
    if any(it.get("lon") is None for it in items):
        src = _osr.SpatialReference()
        src.ImportFromWkt(wkt)
        dst = _osr.SpatialReference()
        dst.ImportFromEPSG(4326)
        for s_ in (src, dst):
            try:
                s_.SetAxisMappingStrategy(_osr.OAMS_TRADITIONAL_GIS_ORDER)
            except AttributeError:
                pass
        tf = _osr.CoordinateTransformation(src, dst)
    rows = []
    for k, it in enumerate(items, 1):
        lon, lat = it.get("lon"), it.get("lat")
        if lon is None and tf is not None:
            lon, lat = tf.TransformPoint(float(it["x"]), float(it["y"]))[:2]
        d = it.get("date")
        rows.append(dict(seed_no=k, seed_id=it.get("id", ""), longitude=f"{float(lon):.6f}",
                         latitude=f"{float(lat):.6f}", x=f"{float(it['x']):.2f}", y=f"{float(it['y']):.2f}",
                         grid_row=it["row"], grid_col=it["col"], date=(d.isoformat() if d else "")))
    fp = out_dir / f"{V910_TAG}_seed_points_used.csv"
    v910_write_csv(fp, ["seed_no", "seed_id", "longitude", "latitude", "x", "y", "grid_row", "grid_col", "date"],
                   rows, [f"{V910_TAG} — points that seeded this run ({len(rows)} rows, "
                          f"{len({(r['grid_row'], r['grid_col']) for r in rows})} grid cells)",
                          f"source: {V910_SEED_INFO.get('table') or V910_SEED_INFO.get('shapefile')}",
                          "x, y in the model CRS; longitude, latitude in WGS84"], logger_)
    fj = out_dir / f"{V910_TAG}_seeding.json"
    fj.write_text(_v910_json.dumps({k: v for k, v in V910_SEED_INFO.items() if not k.startswith("_")},
                                   ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    if logger_ is not None:
        logger_.info(f"  ✓ v9.10 seed list: {fp.name} ({len(rows)} points) + {fj.name}")
    return [fp, fj]


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
class Config:
    CWD = Path(os.getcwd())
    DATA_ROOT = Path("/mnt/c/Users/USER/Downloads") / "Claude Code"
    DL_DIR = Path("/mnt/c/Users/USER/Downloads") / "PPR_Results" / "v9_10"   # v9.10: own subfolder
    ENV_DIR = DATA_ROOT / "PPR_env"
    OUT_DIR = CWD / "results" / "ppr_v9_10_review"           # v9.10: never overwrite the 6 Nov run
    # Downsample factor applied to the gorani TIF for a workable grid size.
    FACTOR = 10
    MONTHLY_BETA = np.array([1.40, 1.35, 1.20, 1.10, 0.95, 0.85, 0.80, 0.85, 1.00, 1.25, 1.35, 1.42])
    DEER_RADIUS_KM = np.array([13.5, 16.4, 20.1, 23.6, 26.0, 26.7, 25.4, 22.5, 18.8, 15.3, 12.9, 12.2], dtype=float)
    TEMP_C = np.array([-3.5, -1.5, 5.0, 12.0, 17.5, 22.0, 25.0, 24.5, 18.5, 11.0, 3.5, -1.5])
    PRECIP_IDX = np.array([0.20, 0.25, 0.40, 0.55, 0.70, 0.85, 1.00, 0.95, 0.65, 0.45, 0.30, 0.20])
    MONTH_NAMES = ["January", "February", "March", "April", "May", "June",
                   "July", "August", "September", "October", "November", "December"]
    MONTH_SHORT = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    SEASONS = {
        "Winter\n(Dec–Feb)": [11, 0, 1],
        "Spring\n(Mar–May)": [2, 3, 4],
        "Summer\n(Jun–Aug)": [5, 6, 7],
        "Autumn\n(Sep–Nov)": [8, 9, 10]
    }
    SEASON_COLORS = {
        "Winter\n(Dec–Feb)": "#cfe2f3",
        "Spring\n(Mar–May)": "#d9ead3",
        "Summer\n(Jun–Aug)": "#fff2cc",
        "Autumn\n(Sep–Nov)": "#fce5cd"
    }
    # LEGACY DMZ crossing points, format (x_km, y_km, label, offset, weight).
    # These are the v2.2 map-coordinate entries. In v9.1 they are superseded by
    # DMZ_CORRIDOR_CROSSINGS (real lat/lon, snapped onto MDL.shp); they are kept
    # only as a last-resort seed fallback and for the AUTO_VALIDATE diagnostic.
    DMZ_ENTRIES = [
        (860, 1995, "Han W", -30, 0.65), (872, 1999, "Han C", -25, 0.75), (888, 2003, "Han E", -15, 0.70),
        (920, 2005, "Kaesong", -5, 0.95), (942, 2005, "Panmunjom", 0, 1.00), (950, 2020, "Imjin", 15, 0.70),
        (1000, 2008, "Cheorwon W", 8, 0.80), (1057, 2007, "Cheorwon E", 5, 0.78),
    ]
    NK_ENTRIES = [
        (835, 2210, "Yalu W", 5, 0.90), (1000, 2200, "Yalu C", 10, 0.65), (1155, 2215, "Tumen", 15, 0.55),
    ]
    AUTO_VALIDATE_ENTRIES = True
    # ══ v9.3: SOUTH-KOREA COUNTY (SIGUNGU) CSV EXPORT ════════════════════════
    # THE POINT OF THIS BLOCK: turn the modelled attack-rate rasters into an
    # administrative table for South Korea — one row per sigungu — written as
    # soon as the cross-MDL (Korea) leg finishes. See the module at the top of
    # the file for how the zonal weights are built.
    EXPORT_SK_COUNTY_CSV = True
    # EXPLICIT PATH — set this and the hint search below is skipped entirely.
    # Use it whenever the county layer does not live beside the other data (for
    # example when the script is run from a different project directory), or
    # when a hint keeps matching the wrong file. This is the unambiguous option:
    #   SK_COUNTY_SHP_PATH = "/mnt/c/Users/USER/Downloads/Claude Code/sig.shp"
    SK_COUNTY_SHP_PATH = None
    # Otherwise: the first hint that resolves to a .shp wins. sig.shp is the
    # Korean national sigungu layer, gadm*_KOR_2 the GADM equivalent.
    # CAUTION: these are SUBSTRING matches — "sig" also matches "de-sig-n.shp",
    # so the run logs every candidate and which one it chose.
    SK_COUNTY_SHP_HINTS = ["sig", "sigungu", "gadm41_kor_2", "gadm36_kor_2",
                           "kor_2", "admin2", "adm2", "county"]
    # Extra directories to search for the county layer, on top of CWD /
    # DATA_ROOT / ENV_DIR. Handy when the model is run from one project folder
    # but the GIS layers live in another.
    SK_COUNTY_EXTRA_SEARCH_DIRS = []
    # None = auto-detect the DBF encoding (Korean shapefiles are usually CP949).
    # Force it with e.g. "CP949" / "EUC-KR" / "UTF-8" if the names come out wrong.
    SK_COUNTY_SHP_ENCODING = None
    # Sub-resolution used to build the county × cell AREA WEIGHTS. 6 means each
    # model cell is split 6×6; a county covering 12 of those 36 sub-cells gets a
    # weight of 12 in that cell. Raise it for small counties on a coarse grid
    # (memory cost is (NR*os)×(NC*os) int32).
    SK_COUNTY_OVERSAMPLE = 6
    # None -> use RAD_ARRIVAL_THRESHOLD, so the CSV and the maps agree.
    SK_COUNTY_ARRIVAL_THRESHOLD = None
    SK_COUNTY_WRITE_MONTHLY_WIDE = True     # counties × months matrix
    SK_COUNTY_WRITE_LONG = True             # tidy county-month rows
    SK_COUNTY_ALSO_MAIN_RUN = True          # also export from the continental run
    # v9.9 FIX C: the county layer is ROK-only, and the old row test dropped ROK counties lying
    # north of the median MDL row (e.g. Goseong-gun). Keep every county.
    SK_COUNTY_DROP_NORTH_OF_MDL = False
    SPLIT_KOREA_BY_BOUNDARY = True          # v9.9 FIX C: DPRK/ROK from sig.shp, not one grid row
    KOREA_VACC_SCENARIOS = True             # v9.9 E: run V1/V2 vaccination on the Korea leg
    KOREA_VACC_V1_START_LEG_MONTH = 6       # V1 start, 0-based Korea-leg month (6 = 7th month after NK introduction)
    # ══ v9.10 — review re-analysis switches (search "v9.10") ═════════════════════
    # ---- item 1: national split, exclaves, estuary, run length ---------------------
    ROK_SPLIT_METHOD = "auto"             # "auto": sigungu sub-cells -> OGR rasterize -> row; "row" = old split
    ROK_SPLIT_USE_LAND_FRACTION = True    # a border cell goes to the country holding most of its LAND
    ROK_SPLIT_CORE_EPS = 0.02             # cells with < 2% ROK area are candidate DPRK land
    ROK_BORDER_MIN_LAND_SHARE = 0.50      # a border cell is ROK if the ROK holds >= 50% of its land
    ROK_EXCLAVE_POLICY = "report"         # "report" | "mask": ROK islands touching DPRK land on the grid
    KOREA_ESTUARY_POLICY = "open"         # "open" (baseline) | "sealed": unfenced contact gets φ0
    KOREA_POST_ARRIVAL_MONTHS = None      # buffer after establishment; None = KOREA_POST_ARRIVAL_YEARS × 12
    KOREA_REPORT_MONTHS = [206]           # extra global months at which Tables 1 and 3 are also given
    # ---- items 2–4: reporting -------------------------------------------------------
    BENCHMARK_ALT_WEIGHTS = True          # distance-only rankings + entropy/equal/rank-sum/leave-one-out
    # Round 2 (Reviewer 2, F1): what each ranking adds to distance from the entry routes
    BEYOND_DISTANCE_ANALYSIS = True       # writes PPR_v9.10_TableS7b_beyond_distance_<domain>.csv
    BEYOND_DISTANCE_REFERENCE = "gates + unfenced border (estuary)"   # Euclidean benchmark of the paper
    BEYOND_DISTANCE_BOOT = 300            # spatial block-bootstrap replicates (0 = no CIs)
    BEYOND_DISTANCE_BLOCK_KM = 32.0       # block size of that bootstrap
    BEYOND_DISTANCE_BOOT_MAX_CELLS = 20000   # larger domains (continental, fine grids): point values only
    BEYOND_DISTANCE_SEED = 20261123
    FENCE_WAIT_PROBE_KM = None            # zone radius around each gate (km); None = band half-width
                                          # + 2 cells (24 km on the 8-km grid), so the zone reaches
                                          # the first DPRK and ROK cells outside the fence band
    RENDER_FIGURE_9 = True
    FIGURE9_MONTH = None                  # None = end of the peninsular leg
    # ---- item 7: rasters for an external test -----------------------------------------
    EXPORT_RISK_RASTERS = True
    SMCDA_ADD_ALL_DOMAIN = True           # also classify the whole domain (China + Korea)
    # ---- items 5–6: batch on the peninsular leg ----------------------------------------
    # CLI: --batch scenarios | oat | morris | lhs | morris+lhs | all   (or set these two)
    SAVE_KOREA_CONTEXT = True             # PPR_v9.10_korea_leg_context.npz -> batch without a re-run
    RUN_SCENARIOS = False
    RUN_SENSITIVITY = "off"               # "off" | "oat" | "morris" | "lhs" | "morris+lhs"
    SA_WORKERS = "auto"                   # worker processes (fork, Linux/WSL); 1 = serial
    SA_SEED = 20261109
    SA_MORRIS_TRAJECTORIES = 12           # r; runs = r × (k + 1) = 204 for k = 16
    SA_MORRIS_LEVELS = 4
    SA_LHS_N = 250
    SA_LHS_TOP_K = 8                      # parameters carried from the Morris screen into the LHS
    SA_LHS_PARAMS = None                  # explicit list of keys overrides the Morris selection
    SA_REPORT_MONTH = 206                 # "% of the ROK reached by month ..." (global month)
    # Supplementary Table S3: name -> (low, high, symbol, description)
    SA_PARAM_RANGES = {
        "gamma":      (0.09, 0.36, "γ",      "recovery/removal rate (month⁻¹); log scale, centred on 0.18"),
        "beta_level": (0.50, 1.50, "β level", "multiplier on the mean of β(m)"),
        "beta_amp":   (0.00, 1.50, "β amp.",  "multiplier on the seasonal amplitude of β(m)"),
        "tv_amp":     (0.00, 1.00, "T_v",    "viral-survival seasonality (0 = off, 1 = as baseline)"),
        "diff_max_x": (0.50, 1.50, "D_max",  "multiplier on peak diffusivity (0.45–1.35 cell² month⁻¹)"),
        "adv_max_x":  (0.50, 1.50, "v_max",  "multiplier on the advection speed cap"),
        "adv_gain_x": (0.50, 1.50, "a_g",    "multiplier on the advection gain"),
        "kappa":      (1.00, 3.00, "κ",      "conductance exponent"),
        "phi0":       (0.00, 0.12, "φ0",     "pass-through of the sealed fence"),
        "phimax_x":   (0.50, 1.50, "φmax",   "multiplier on the gate ceiling (0.531)"),
        "rg_km":      (7.0, 28.0,  "r_g",    "gate half-width (km)"),
        "f0":         (0.25, 1.00, "f0",     "winter gate opening (fraction of the peak)"),
        "psi":        (0.60, 1.00, "ψ",      "weight of the water-deer density composite in habitat"),
        "w_h":        (0.40, 0.80, "w_h",    "habitat weight in transmissibility (w_l = 1 − w_h)"),
        "a_star":     (0.01, 0.10, "A*",     "arrival threshold"),
        "i0":         (0.30, 0.90, "I0",     "seed intensity"),
    }
    SA_LOG_PARAMS = ("gamma",)            # round 2: sampled log-uniformly (baseline 0.18 = geometric centre)
    SA_LHS_PARAMS_R2 = ["gamma", "beta_level", "w_h", "diff_max_x", "beta_amp", "psi", "tv_amp", "a_star"]
    # Round 2: comparator hosts (same seeding, dispersal kernel and barrier) and vaccination V2 under each
    COMPARATOR_LIST = ["baseline", "host_neutral", "host_livestock", "dispersal_flat", "fence_uniform_mean",
                       "baseline_V2", "host_neutral_V2", "host_livestock_V2"]
    R2_ZONE_RADII_KM = (12.0, 24.0, 36.0)  # gate-zone radii for the fence-wait sensitivity
    R2_TARGET_SHARES = (0.05, 0.10, 0.20, 0.33, 0.50)   # share of ROK cells placed under surveillance
    R2_CONDITIONAL_MONTH = 218            # end of the peninsular leg in the run of record
    R2_TWO_STEP_BANDS_KM = (16.0, 24.0, 32.0)   # distance bands of the two-step rankings (F1)
    R2_TWO_STEP_WEIGHT = 0.2              # weight of the risk score in the weighted ranking (Table S7b)
    R2_TARGET_BOOT = 1000                 # block-bootstrap replicates for rho of each ranking (0 = off)
    R2_TARGET_BOOT_SEED = 20261125        # blocks of BEYOND_DISTANCE_BLOCK_KM, as for Table S7b
    R2_TARGET_CELL_TABLE = True           # per-cell table of the ROK cells of the targeting comparison
    POINT_TO_CELL = "floor"               # point locations on the grid: "floor" = the cell that contains the point
                                          # (default, as reported) | "round" = nearest cell corner (runs before 26 Nov)
    R2_SKIP_COMPARATORS = False           # --batch targeting: zone radii and targeting only
    R2_VAX_DECOMP = False                 # --batch vaxdecomp: V2 decomposition runs only (Reviewer 2, FC1)
    R2_COMPARATOR_STEM = "R2_comparators"
    VAX_DECOMP_LIST = ["baseline", "baseline_V2", "baseline_V2_allhosts", "baseline_V2_fence",
                       "baseline_V2_allhosts_fence", "host_livestock_V2"]
    # 19 Nov (round 2): small-ruminant (goat) layer. The GLW4 raster was clipped at the MDL, so the
    # ROK had no goats; the "rokfilled" file adds the national goat register south of the MDL.
    LIVESTOCK_FILE = "auto"               # "auto" | "rokfilled" | "glw" | a file name in ENV_DIR
    LIVESTOCK_ROKFILLED_NAME = "GLW4GTS_1km_5179_ROKfilled.tif"
    LIVESTOCK_GLW_NAME = "GLW4GTS_1km_5179.tif"
    LIVESTOCK_USED = None                 # set at run time (file actually read)
    RUN_R2 = False                        # set by --batch r2 / comparators
    # Table 4 scenarios ('drop_each' = one run per gate; 'astar_*' come from the baseline run)
    SCENARIO_LIST = ["baseline", "entry_yalu", "fence_uniform_phimax", "fence_uniform_mean",
                     "fence_sealed_no_gates", "estuary_toggle", "gamma_x2", "gamma_x4",
                     "no_seasonality", "astar_0.01", "astar_0.10", "corridors_equal",
                     "corridors_x0.5", "corridors_x1.5", "drop_each", "add_inje", "exclaves_toggle"]
    YALU_ENTRY_LONLAT = (124.40, 40.10)   # lower Yalu (Sinuiju–Dandong)
    YALU_CHINA_PROBE_KM = 60.0            # China cells that time the Yalu introduction
    INJE_SEORAK_ESTIMATE_LATLON = (38.30, 128.20)   # ESTIMATED (Seohwa-myeon, Inje-gun), snapped to the MDL
    INJE_SEORAK_WEIGHT = 0.40
    GRID4KM_SUMMARY = None                # 4-km run summary JSON (None = <OUT_DIR>_grid4km/ auto-detect)
    # ---- 11 Nov: publication figures (--preset pubfigs) -------------------------------------
    PUB_FIGS = False                      # True: no embedded titles, plain-language labels, 600 dpi,
                                          # north arrows, Figure 4 zoomed with leader lines, Figure 8 stacked
    PUB_FIGS_DPI = 600
    HAB_PENINSULA_REPORT = False          # True (--preset habpeninsula): classify the peninsular domain
                                          # (Supplementary Figure S2), write the numbers, then stop
    # ---- item 6: 4-km grid (set by --preset grid4km) ---------------------------------------
    TRANSPORT_REFERENCE_CELL_KM = None    # "auto" or km: rescale transport to the reference grid
    REFERENCE_TARGET_MAX_PX = 750         # the settings that define the reference (8-km) grid
    REFERENCE_MAX_GRID_CELLS = 700_000
    KOREA_ONLY_INTRO_MONTH = None         # continental leg only this long; DPRK introduction at this month
    # ---- item 7: temporal hold-out for the external test (set by --preset holdout2015) ------
    # The run of record seeds from every outbreak in the seed shapefile (EMPRES-i China
    # 2007–2025 + Mongolia), so its maps cannot be tested against those outbreaks. A hold-out
    # run seeds only from the outbreaks reported before SEED_BEFORE; ppr_external_test.py then
    # scores its maps against the later ones.
    SEED_TABLE = None                     # None = outbreak shapefile (run of record); a CSV / Excel
                                          # path, or "auto" = PPR_seed_points_dated.csv
    SEED_TABLE_DATE_COL = "date"          # 'a|b' = first non-empty of several columns
    SEED_BEFORE = None                    # "2015-01-01": keep rows dated BEFORE this day
    SEED_UNDATED = "drop"                 # rows without a date: "drop" | "keep"
    WRITE_SEED_LIST = True                # <OUT_DIR>/PPR_v9.10_seed_points_used.csv + _seeding.json
    MAXENT_TIF_PATH = None                # explicit MaxEnt map (None = search by MAXENT_TIF_HINTS);
                                          # for a hold-out run with MaxEnt refitted on pre-cut-off points
    SK_COUNTY_MAP = True                    # choropleth PNG (arrival + final attack)
    SK_COUNTY_CSV_STEM = "PPR_v9.3_south_korea_counties"
    # Cache the Korea-leg attack rasters to .npz alongside the CSVs.
    # WHY THIS EXISTS: the adaptive horizon means the Korea leg can be a 20+ year
    # simulation. If the county layer is missing or misconfigured, the epidemic
    # itself was still computed correctly — losing it because a shapefile was in
    # the wrong folder is pure waste. With the cache, `sk_county_tools.py
    # --rebuild` regenerates every county CSV in seconds, and you can iterate on
    # the zonal settings (oversample, encoding, thresholds) without re-running
    # the model once.
    SAVE_EPIDEMIC_CACHE = True
    CACHE_INCLUDE_INFECTED = True    # also store monthly I (peak-infection columns)
    # The Korea crop is small (a few hundred cells across) so its cache is tens of
    # MB. The CONTINENTAL run is the whole China-scale grid × every month, which
    # runs to hundreds of MB — off by default; turn it on only if you need to
    # retune the county table against the continental run specifically.
    CACHE_MAIN_RUN = False
    # Nominal calendar anchor for month 1 of the simulation. The model has no
    # absolute clock — this only LABELS the months so the CSV is readable
    # ("2029-07" rather than "month 43"). Set SIM_START_YEAR = None to drop the
    # calendar columns entirely.
    SIM_START_YEAR = 2026
    SIM_START_MONTH = 1
    # ══ v9.1: DMZ CORRIDOR CROSSINGS — GATES IN THE MDL LINE ═════════════════
    # The MDL stops being a uniformly leaky wall (Han φ everywhere) and becomes a
    # SEALED barrier that is permeable only at the named water-deer corridors.
    # The corridor table itself (name, lat, lon, weight, basis) lives at the top
    # of this file in DMZ_CORRIDOR_CROSSINGS.
    USE_CORRIDOR_CROSSINGS = True
    # Gate geometry / strength.
    # ── DIAGNOSIS (why PPR was stuck in NK): the previous combination of
    #    CORRIDOR_MAX_PASS=0.35, MDL_CROSSING_SLOWDOWN=0.60, MDL_DILATE_ITER=2
    #    gave an effective Cheorwon throughput of 0.35×0.60² ≈ 0.126 per cell —
    #    two cells of compounding = ~1.6% of flux.  That is a closed wall, not
    #    a leaky fence. Restored to ecologically defensible values:
    #   -> Cheorwon (w=1.00): 0.531 * 1.00 * 1.00 = 0.531  (Han φ, as intended)
    #   -> Goseong  (w=0.60): 0.531 * 0.60 * 1.00 = 0.319
    #   -> between gates:     0.06  (near-sealed, not 0.02 / closed)
    CORRIDOR_GATE_RADIUS_KM = 14.0     # gate half-width along MDL (wider = more cross-section)
    CORRIDOR_MAX_PASS = 0.531          # Han et al. φ — the empirically estimated ceiling
    MDL_SEALED_FACTOR = 0.06           # non-corridor cells: very restricted but NOT zero
    MDL_CROSSING_SLOWDOWN = 1.00       # no extra multiplier; the gate weights alone throttle
    CORRIDOR_SNAP_MAX_KM = 25.0        # warn if a point is further than this from MDL.shp
    # Seasonal gating: gates open with the gorani dispersal season, so crossings
    # concentrate in late spring / early summer instead of running all year.
    CORRIDOR_SEASONAL_GATING = True
    # v9.2: 0.35 kept the gates near-shut for most of the year, which is a large
    # part of the 169-month NK->SK lag. Water deer still cross outside the
    # dispersal peak, so the winter floor is raised.
    CORRIDOR_SEASON_FLOOR = 0.50       # winter opening, as a fraction of the peak
    # Single-cell fence band (MDL_DILATE_ITER=1): the front crosses in one step
    # rather than compounding attenuation across multiple cells. Two or more cells
    # at sealed=0.06 give 0.06^N which approaches zero for N≥2 — that is why
    # PPR was stuck in North Korea in all previous panels.
    MDL_DILATE_ITER = 1
    CORRIDOR_SEED_FALLBACK = True      # seed from the gates when no outbreak shapefile
    SHOW_CORRIDOR_MARKERS = True       # ▼ markers, size ∝ ecological weight
    CORRIDOR_PROBE_KM = 15.0           # SK probe depth for the crossing report
    SHOW_MDL_GATE_MAP = True           # standalone gate-permeability map
    # ══ SOUTH-KOREA PANEL STYLING (v9.2) ═════════════════════════════════════
    # v9.1 painted the habitat-suitability underlay across the WHOLE peninsula in
    # draw_korea(), so South Korea rendered uniformly tan/brown even at 0% attack
    # (GOR_EFF ≈ 0.4-0.6 lands on the pale-orange part of the MaxEnt ramp). The
    # suitability underlay is now confined to North Korea (the source region);
    # South Korea is a blank white canvas carrying ONLY water and % attack bands.
    SK_PANEL_PLAIN = True                    # South Korea = white canvas
    SK_PANEL_LAND_COLOR = "#ffffff"
    NK_PANEL_SHOW_SUITABILITY = True         # keep the underlay NORTH of the MDL
    SK_PANEL_SHOW_WATER = True               # water.shp bodies / waterways only
    SK_PANEL_WATER_LINE_COLOR = "#5fa8d3"
    SK_PANEL_WATER_LINE_WIDTH = 0.35
    SK_PANEL_WATER_LINE_ALPHA = 0.75
    # water.shp is a full national hydrography layer (every stream, canal and
    # reservoir). v9.1 hid it under the tan underlay in a near-white colour; on
    # the new white canvas the layer's true density fills the whole country with
    # ink. These two knobs keep only water features that are RESOLVABLE at the
    # model's cell size — a 200 m stream is meaningless on a multi-km grid.
    SK_PANEL_WATER_MIN_EXTENT_CELLS = 2.0    # feature must span ≥ this many cells
    SK_PANEL_WATER_MAX_FEATURES = 150        # keep the N largest; 0 = no cap
    SK_PANEL_WATER_SIMPLIFY_STEP = 2         # keep every Nth vertex (1 = no thinning)
    # v9.2: South Korea is a WHITE canvas (no suitability underlay) but it is NOT
    # blank — it carries the sigungu county outlines and the water bodies, so the
    # % attack bands sit on a readable administrative/hydrographic base.
    SK_PANEL_SHOW_COUNTIES = True            # sig.shp county outlines ON
    SK_PANEL_COUNTY_LINE_COLOR = "#9a9a9a"
    SK_PANEL_COUNTY_LINE_WIDTH = 0.25
    SK_PANEL_COUNTY_LINE_ALPHA = 0.55
    SK_PANEL_OUTLINE = True                  # thin coastline so white SK is legible
    SK_PANEL_OUTLINE_COLOR = "#9aa3ab"
    SK_PANEL_OUTLINE_WIDTH = 0.6
    SK_BAND_MIN_PCT_VISIBLE = 1.0            # SK floor (NK keeps BAND_MIN_PCT_VISIBLE)
    # Extra years for the Korea (NK -> cross-MDL -> SK) sub-model, on top of
    # SIM_YEARS. With the fence gated the crossing is slow, so the peninsula leg
    # needs a longer horizon than the China source run to show the full picture.
    KOREA_EXTRA_YEARS = 2
    OUT_DIR_LEGACY = CWD / "results" / "ppr_v8_1_percent_spread_gorani_dots"
    # ── NEURAL COMPONENT SETTINGS (v7.0) ─────────────────────────────────────
    USE_NEURAL_COMPONENTS = True
    TRAIN_NEURAL = False
    OUTBREAK_CSV = CWD / "ppr_outbreak_locations.csv"
    OUTBREAK_CRS = "EPSG:5179"
    N_TRAIN_EPOCHS = 50
    LEARNING_RATE = 1e-3
    KERNEL_NET_HIDDEN = 32      # hidden channels in the convolutional kernel net
    KERNEL_NET_LAYERS = 3       # number of conv layers
    NODE_HIDDEN = 64            # hidden size of the ODE function MLP
    NODE_SOLVER = "dopri5"      # ODE solver: dopri5 (adaptive RK45), euler, rk4
    NODE_RTOL = 1e-3
    NODE_ATOL = 1e-4
    # Feature channels fed to SpatialKernelNet (order matters):
    # 0 GOR_nn, 1 ELEV_n, 2 SLOPE_N, 3 livestock_n, 4 RIV_nn,
    # 5 CONDUCTANCE, 6 dijkstra_dist_km, 7 I_proxy
    N_KERNEL_FEATURES = 8
    # Recovery/removal rate gamma (fraction of I removed per month).
    GAMMA_RECOVERY = 0.18
    # ══ REACTION-ADVECTION-DIFFUSION ENGINE ══════════════════════════════════
    RAD_ENGINE = True
    # Number of YEARS to simulate. Seasonality cycles each year; S/I/R persist
    # across the whole run. v9.1 raises this because the gated fence deliberately
    # delays the North -> South crossing.
    # ══ ELASTICITY-RANKED LANDSCAPE PARAMETERS (v9.2 sensitivity tuning) ═════
    # These were hard-coded literals in v9.1. The local sensitivity analysis
    # ranks them by mean elasticity (% change in spread per % change in input),
    # so they are exposed here and moved in that order — the highest-elasticity
    # knob buys the most acceleration per unit of parameter distortion.
    #
    #   CONDUCTANCE_EXPONENT  mean 0.549 / max 2.735   HIGH
    #   SLOPE_WEIGHT          mean 0.293 / max 1.453   MODERATE
    #   DEER_RADIUS_SCALE     mean 0.213 / max 0.452   MODERATE
    #   MONTHLY_BETA          mean 0.198 / max 0.684   LOW  <- deliberately NOT raised
    #
    # k in CONDUCTANCE = exp(-k * RESISTANCE). Lower k = a more permeable
    # landscape = higher mobility M = a faster front. v9.1 used 3.0.
    CONDUCTANCE_EXPONENT = 2.00
    # RESISTANCE weights (auto-normalised to sum 1). Slope carried 0.35 in v9.1;
    # trimmed per its 0.293 elasticity, with the freed weight spread over the
    # other three terms so total resistance keeps the same scale.
    RESISTANCE_WEIGHTS = {"slope": 0.28, "urban": 0.27, "elev": 0.22, "water": 0.23}
    # Multiplier on the empirical monthly gorani dispersal radius. Keep at 1.00
    # unless you can defend it — DEER_RADIUS_KM is observational data, not a
    # free parameter. Raise it only if the front is still too slow after the
    # transport caps below are corrected.
    DEER_RADIUS_SCALE = 1.00
    # One global lever on transport (diffusion + advection). Use this to hit a
    # target arrival time without touching the individual ecological terms;
    # report it in the methods as a transport-calibration factor.
    SPREAD_SPEED_MULTIPLIER = 1.00
    # ══ ADAPTIVE RUN LENGTH (v9.2) ═══════════════════════════════════════
    # The run is no longer a fixed horizon that can truncate the epidemic.
    # It continues until PPR has actually (1) reached the North-Korean border
    # from China and (2) crossed the gated MDL and ESTABLISHED in South Korea,
    # then keeps going for a buffer so the southward spread is visible.
    # SIM_YEARS / KOREA_EXTRA_YEARS are now only the MINIMUM horizon.
    AUTO_RUN_UNTIL_SK_ESTABLISHED = True
    SK_ESTABLISHED_MEAN_ATTACK = 0.05     # SK mean cumulative attack ≥ 5%, OR
    SK_ESTABLISHED_AREA_FRACTION = 0.10   # ≥10% of SK cells past the arrival threshold
    KOREA_POST_ARRIVAL_YEARS = 3          # keep simulating this long after establishment
    CHINA_POST_BORDER_YEARS = 1           # buffer after China reaches the NK border
    MAX_SIM_YEARS = 40                    # hard cap on the Korea leg (safety)
    MAIN_MAX_YEARS = 15                   # hard cap on the China/continental run
    SIM_YEARS = 6                         # MINIMUM years for the source run
    # v9.2: the v9.1 caps were the real reason the run needed 20 years.
    # RAD_ADV_MAX = 0.28 cells/month on a ~10-15 km grid is ~3-4 km/month of
    # directed movement, against an EMPIRICAL dispersal radius of 12-27 km/month
    # (DEER_RADIUS_KM). The front was ~5x slower than the animal data underneath
    # it. Raised to be consistent with that radius; the fence, the barriers and
    # the habitat surface still do the throttling. (v9.1: 0.45 / 0.28 / 1.80.)
    RAD_DIFF_MAX = 0.90               # peak random-movement diffusivity (cell^2/month)
    RAD_DIFF_MIN = 0.02
    RAD_ADV_MAX = 0.55                # max directed speed (cells/month; CFL cap)
    RAD_ADV_GAIN = 2.40               # habitat-gradient -> velocity gain
    RAD_TRANS_HAB_W = 0.60
    RAD_TRANS_LSTK_W = 0.40
    RAD_TRANS_BASE = 0.03             # unsuitable land barely transmits
    RAD_SEED_INTENSITY = 0.60
    RAD_SOURCE_CONTINUOUS = False     # seed ONCE and let it spread as a wave
    RAD_SOURCE_RATE = 0.02
    RAD_CFL_SAFETY = 0.40
    RAD_MAX_SUBSTEPS = 130
    RAD_ARRIVAL_THRESHOLD = 0.05      # attack-rate level defining "arrived"
    RAD_USE_LEARNED_CLOSURE = False   # UDE hook (closure is swappable)
    # ══ VACCINATION (intervention) ═══════════════════════════════════════════
    ENABLE_VACCINATION = True
    VACC_MODE = "reactive_ring"          # "reactive_ring" | "preemptive" | "continuous"
    VACC_START_MONTH = 6                 # 0-based month the campaign begins
    VACC_COVERAGE = 0.60
    VACC_RATE = 0.30
    VACC_EFFICACY = 0.85
    VACC_RING_KM = 50.0
    VACC_INFECTED_THRESHOLD = 0.02
    VACC_TARGET_REGION = "korea"         # "all" | "korea" | "nk" | "sk"
    SHOW_VACCINATION_MAP = True
    # ══ HABITAT / ACCESSIBLE / BARRIER CLASSIFICATION (Han et al. Fig.2) ══════
    SHOW_HABITAT_CLASSIFICATION = True
    HAB_CLASS_THRESHOLD_MODE = "quantile"   # "quantile" | "absolute"
    HAB_SUIT_QUANTILE = 0.45
    HAB_SUIT_THRESHOLD = 0.50
    ACCESS_COSTDIST_MAX_KM = 15.0
    CLASS_COLOR_HABITAT = "#77dd77"
    CLASS_COLOR_ACCESSIBLE = "#6db3f2"
    CLASS_COLOR_BARRIER = "#f4e04d"
    CLASS_MAP_OVERLAY_REACHED = False
    USE_HABITAT_CLASS_IN_TRANSPORT = False
    CLASS_MOBILITY_HABITAT = 1.00
    CLASS_MOBILITY_ACCESSIBLE = 0.50
    CLASS_MOBILITY_BARRIER = 0.05
    # ══ v9.7: SPATIAL MULTI-CRITERIA DECISION ANALYSIS — RISK-AREA THRESHOLDS ═
    #
    # WHY THIS EXISTS. The habitat / accessible / barrier split above is a
    # SINGLE-criterion cut on GOR_EFF at one quantile. It answers "is this cell
    # gorani habitat?". It does NOT answer "what level of PPR risk does this
    # area carry, and where do the boundaries between levels fall?" — that is a
    # MULTI-criteria question: host density, wild-reservoir density,
    # environmental suitability, landscape connectivity, distance from the entry
    # route and local transmission potential all bear on it, and they do not
    # bear on it equally.
    #
    # SMCDA answers it in four separable steps, so each can be audited alone:
    #   1. STANDARDISE  every criterion onto [0,1], 1 = maximum contribution to
    #                   risk. Objective, given the layers.
    #   2. WEIGHT       the criteria. THIS IS THE ONLY STEP WHERE JUDGEMENT
    #                   ENTERS. AHP makes the judgement explicit and testable
    #                   (Saaty consistency ratio).
    #   3. AGGREGATE    to one risk score. WLC = full trade-off between
    #                   criteria; OWA = a controllable degree of trade-off.
    #   4. THRESHOLD    the score into named risk levels.
    #                   >>> THE THRESHOLDS ARE THE DELIVERABLE OF THIS BLOCK. <<<
    #
    # The result is then VALIDATED against the mechanistic RAD epidemic (mean
    # arrival month per class, AUC, Spearman rho). The score is built ONLY from
    # static layers and never sees the epidemic, so that comparison is a real
    # test rather than a restatement of the fit.
    ENABLE_SMCDA = True
    SMCDA_EXPORT_COUNTY_CSV = True      # per-sigungu risk levels (needs sig.shp)
    SHOW_SMCDA_MAP = True
    #
    # ── 0. WHICH AREA IS BEING CLASSIFIED — AND WHY IT MATTERS ───────────────
    #
    # A risk CLASS is relative to the domain it was computed over. "High risk"
    # on a China-wide grid and "high risk" on the Korean peninsula are different
    # statements, because the score is standardised and thresholded within
    # whatever domain you give it. Classify the area the decision is about.
    #
    # A domain is also tied to an INTRODUCTION ROUTE. PPR reaches South Korea
    # through the water-deer corridor gates in the MDL; it reaches North Korea
    # from China across the northern border; it reaches the continental grid
    # from the outbreak points themselves. The ENTRY criterion has to measure
    # distance from the route that applies to the domain being classified — see
    # SMCDA_ENTRY_SOURCE. Get that wrong and the classification will rank the
    # epidemic backwards, which is exactly what the validation pass detects.
    #
    #   "all"    the whole model domain (China + Korea, one comparable scale).
    #            Continental context. Entry = the outbreak seed points.
    #   "korea"  the Korean peninsula, NK + SK on one scale. This is the domain
    #            of the NK -> corridor gates -> SK narrative and of the
    #            progression figure. Entry = the northern NK border with China.
    #   "sk"     ALL of SOUTH KOREA, whether or not PPR reaches it. A PROSPECTIVE
    #            map: "if PPR enters through the gates, where is risk highest?"
    #            Entry = the water-deer corridor gates.
    #   "nk"     North Korea only. Entry = the northern border with China.
    #
    #   "sk_reached"     SOUTH KOREA, RESTRICTED TO THE AREA PPR ACTUALLY REACHED
    #   "korea_reached"  by SMCDA_REACHED_BY_MONTH — the invaded footprint, i.e.
    #                    the coloured part of the progression map, not the whole
    #                    country. This answers a DIFFERENT and more operational
    #                    question: "within the area the incursion actually
    #                    covers, where does surveillance effort go first?"
    #
    # ‼ WHAT THE _reached DOMAINS COST YOU — READ BEFORE QUOTING THEIR STATS ‼
    # Every other domain is defined WITHOUT reference to the epidemic, which is
    # what makes validating it against that epidemic a real test. A _reached
    # domain is defined BY the epidemic outcome. The score values are still
    # computed from static layers alone and still never see the simulation — but
    # the DOMAIN is conditioned on it, and standardisation ranges and Jenks
    # breaks are computed within that domain. Three consequences:
    #   1. AUC against "reached by end of run" becomes UNDEFINED, because every
    #      cell in the domain is reached by construction. The code reports this
    #      rather than printing a misleading number. Monotonicity of arrival
    #      month and Spearman rho remain valid — and they are now the HARDER
    #      test, because the score can no longer earn credit for the easy part
    #      (telling reached from unreached).
    #   2. The class thresholds are NOT comparable with the "sk" thresholds.
    #      They are cut through a different, smaller distribution.
    #   3. It is a CONDITIONAL risk map ("given the incursion got this far"),
    #      not a prospective one. For pre-incursion planning use "sk".
    #
    # Every domain listed is computed, mapped, validated and exported. They are
    # cheap next to the epidemic, and comparing them is informative.
    #
    # ‼ v9.8 — WHY THE PENINSULA AND CONTINENTAL DOMAINS ARE OFF BY DEFAULT ‼
    # In the v9.7 run the "korea_reached" map came out with North Korea red and
    # essentially ALL of South Korea blue. That is not a finding about South
    # Korea; it is an artefact of the domain. NK and SK were standardised and
    # thresholded on ONE shared scale, NK scores higher on the criteria that
    # carry most of the weight (it is nearer the introduction route and holds
    # the higher small-ruminant densities in this model), so the entire South
    # Korean distribution falls below the upper Jenks cuts and collapses into
    # the two lowest classes. Every South Korean county then reads "Very low /
    # Low" and the map cannot discriminate WITHIN South Korea at all — which is
    # the only discrimination a South Korean surveillance design can use.
    #   That same domain ALSO failed the arrival-order validation in that run
    # (arrival order NOT monotonic, rho = -0.128), for the related reason that
    # one ENTRY origin cannot describe both legs of the incursion.
    #
    # So the DEFAULT run is now South Korea only. Both SK domains are kept
    # because they answer different questions (see above):
    #     "sk_reached"  conditional  — within the invaded footprint, where first?
    #     "sk"          prospective  — if PPR enters, where is risk highest?
    # The peninsula and continental domains are NOT deleted. They remain fully
    # supported and are the right domains for a NK-vs-SK or a China-vs-Korea
    # comparison. Re-enable them by putting them back in this list:
    #     SMCDA_DOMAINS = ["sk_reached", "sk", "korea_reached", "korea", "all"]
    # but do not read a South Korean risk level off those maps.
    SMCDA_DOMAINS = ["sk_reached", "sk"]
    # Which one drives the county CSV, the mobility field and the JSON summary.
    SMCDA_PRIMARY_DOMAIN = "sk_reached"
    # Retained so an older config still runs: when SMCDA_DOMAINS is None this
    # single value is used instead.
    SMCDA_APPLY_TO_DOMAIN = "sk_reached"
    #
    # ── 0a. WHAT COUNTS AS "REACHED", FOR THE _reached DOMAINS ───────────────
    # A cell is in the footprint when its modelled cumulative attack rate has
    # crossed RAD_ARRIVAL_THRESHOLD at or before this month, on the GLOBAL
    # simulation timeline (the same month numbers as the progression figure
    # titles, e.g. "May · Year 18 (month 209)").
    #   None = the last month of the Korea leg, i.e. the full extent of the
    #          incursion as simulated. This is the default and it matches the
    #          final panel of the progression figure.
    #   int  = a specific month. Set 209 to classify exactly the footprint at
    #          18 years, independently of how long the run was extended.
    SMCDA_REACHED_BY_MONTH = None
    # Refuse to threshold a footprint smaller than this — with a handful of
    # cells the breaks are noise, and an empty "Very high" class would be an
    # artefact of the sample size rather than a finding.
    SMCDA_REACHED_MIN_CELLS = 100
    #
    # ── 0b. THE ENTRY CRITERION'S ORIGIN ─────────────────────────────────────
    #   "auto"          pick per domain, as described above. RECOMMENDED.
    #   "gates"         cost-distance from the water-deer corridor gates.
    #   "mdl"           cost-distance from the whole MDL line (spreads entry
    #                   risk evenly along 250 km of fence instead of
    #                   concentrating it at the crossings — a different claim).
    #   "seeds"         cost-distance from the actual RAD seed points.
    #   "china_border"  cost-distance from the China–Korea border zone.
    #   "off"           drop the ENTRY criterion and renormalise the weights.
    SMCDA_ENTRY_SOURCE = "auto"
    #
    # ── 1. CRITERIA ──────────────────────────────────────────────────────────
    # (key, layer, label, direction, standardisation, a, b, veto_below)
    #   direction      "benefit" = more is riskier; "cost" = more is SAFER and is
    #                  inverted after scaling.
    #   standardise    "minmax"    linear on the observed range. Keeps the shape
    #                              of the layer but is hostage to its outliers.
    #                  "rank"      percentile. Immune to outliers and to the
    #                              layer's units; destroys the spacing between
    #                              values, so two cells one unit apart can end up
    #                              far apart if the layer is dense there.
    #                  "fuzzy_lin" linear between a and b, flat outside. Use when
    #                              a and b are DEFENSIBLE thresholds, not just
    #                              convenient numbers.
    #                  "fuzzy_sig" sigmoid with inflection a, spread b.
    #                  "none"      the layer is already on [0,1] and meaningful.
    #   veto_below     a NON-COMPENSATORY constraint. A cell scoring below this
    #                  on this criterion is forced to 0 overall no matter how
    #                  well it scores elsewhere. Leave None unless the constraint
    #                  is genuinely absolute. None are set by default.
    #
    # Every layer named here is one the model already builds. A criterion whose
    # layer is missing at run time is DROPPED with a warning and the weights are
    # renormalised over what remains — it is never silently filled with a
    # default, and the run log always states which criteria actually voted.
    SMCDA_CRITERIA = [
        dict(key="HOST",  layer="livestock_nn",    label="Small-ruminant host density",
             direction="benefit", standardise="rank"),
        dict(key="ENTRY", layer="_SMCDA_GATE_DIST", label="Cost-distance from the entry route",
             direction="cost",    standardise="rank"),
        dict(key="TRANS", layer="_SMCDA_BETA",      label="Local transmission potential (beta)",
             direction="benefit", standardise="minmax"),
        dict(key="SUIT",  layer="GOR_EFF",          label="Environmental / habitat suitability",
             direction="benefit", standardise="rank"),
        dict(key="WILD",  layer="GORANI_DENSITY",   label="Water-deer (gorani) density",
             direction="benefit", standardise="rank"),
        dict(key="CONN",  layer="CONDUCTANCE",      label="Landscape conductance",
             direction="benefit", standardise="rank"),
        dict(key="WATER", layer="WATER_ATTRACT",    label="Proximity to surface water",
             direction="benefit", standardise="minmax"),
    ]
    #
    # ── 2. WEIGHTS ───────────────────────────────────────────────────────────
    SMCDA_WEIGHT_METHOD = "ahp"         # "ahp" | "rank_sum" | "entropy" | "equal"
    #
    # ‼ READ THIS BEFORE USING THE OUTPUT OPERATIONALLY ‼
    # The pairwise matrix below is a DOCUMENTED DEFAULT, not the product of an
    # expert elicitation session. It was CONSTRUCTED from the ordinal priority
    # vector in SMCDA_PRIORITY_VECTOR by taking a_ij = round(s_i / s_j) clipped
    # to Saaty's 1-9 scale. That construction is transparent and near-consistent
    # by design — a low consistency ratio here therefore says the matrix is
    # arithmetically coherent, NOT that the judgement is right.
    # Replace it with your own panel's pairwise judgements before the classes
    # are used to allocate surveillance effort. Every run prints the AHP weights
    # beside rank-sum and entropy weights, and runs a one-at-a-time weight
    # sensitivity analysis, so the dependence on this choice is visible.
    #
    # The ordering encoded below, and the reasoning for it:
    #   HOST  9  PPR is primarily a small-ruminant disease. No hosts, no outbreak
    #            to detect — this is the one criterion that cannot be compensated
    #            for in practice.
    #   ENTRY 8  The design concerns an incursion from outside. In the early
    #            years, distance from the entry route is the dominant
    #            epidemiological gradient.
    #   TRANS 6  The model's own local transmissibility surface (habitat x
    #            livestock), i.e. what happens once the virus arrives.
    #   SUIT  5  Environmental suitability for host presence and virus survival.
    #   WILD  4  Gorani density — the wildlife bridge across the fence.
    #   CONN  3  Landscape conductance — how easily the front propagates onward.
    #   WATER 2  Water attraction; a secondary driver of gorani aggregation,
    #            already partly represented inside WILD.
    SMCDA_PRIORITY_VECTOR = {"HOST": 9.0, "ENTRY": 8.0, "TRANS": 6.0,
                             "SUIT": 5.0, "WILD": 4.0, "CONN": 3.0, "WATER": 2.0}
    # Set an explicit matrix here (list of lists, in SMCDA_CRITERIA key order) to
    # override the constructed one. None = build it from the priority vector.
    SMCDA_AHP_MATRIX = None
    SMCDA_AHP_CR_LIMIT = 0.10           # Saaty's conventional acceptance limit
    SMCDA_AHP_CR_FATAL = False          # True = abort the run on an inconsistent
    #                                     matrix instead of warning and continuing
    SMCDA_RANK_ORDER = ["HOST", "ENTRY", "TRANS", "SUIT", "WILD", "CONN", "WATER"]
    #
    # ── 3. AGGREGATION ───────────────────────────────────────────────────────
    SMCDA_AGGREGATION = "wlc"           # "wlc" | "owa"
    # OWA only. 0.5 reproduces WLC exactly. Raising andness makes the score
    # depend more on a cell's WORST criteria, i.e. it stops a single very high
    # criterion from carrying a cell into a high-risk class on its own.
    SMCDA_OWA_ANDNESS = 0.65
    #
    # ── 4. THRESHOLDS (the deliverable) ──────────────────────────────────────
    SMCDA_N_CLASSES = 5
    SMCDA_CLASS_LABELS = ["Very low", "Low", "Moderate", "High", "Very high"]
    SMCDA_BREAK_METHOD = "jenks"        # "jenks" | "quantile" | "equal" | "sd" | "absolute"
    #   jenks     Fisher-Jenks natural breaks. Minimises within-class variance;
    #             puts the cuts where the score distribution is genuinely sparse.
    #             DATA-DRIVEN: the cut values move if the layers or the domain
    #             change, so they are NOT transferable between runs or domains.
    #   quantile  equal COUNT per class. Guarantees a fixed share of the area in
    #             each class, which is convenient for budgeting and meaningless
    #             epidemiologically — it will manufacture a "very high" class
    #             even from a uniformly low-risk surface.
    #   equal     equal score WIDTH. Interpretable, but classes can be empty.
    #   sd        mean +/- k standard deviations. Assumes near-normal scores.
    #   absolute  YOUR OWN cut values from SMCDA_ABSOLUTE_BREAKS. The only option
    #             that is stable across runs and comparable between domains.
    #             Use it once the score construction is settled.
    SMCDA_ABSOLUTE_BREAKS = [0.20, 0.40, 0.60, 0.80]   # interior cuts, ascending
    #
    # ── 5. VALIDATION + SENSITIVITY ──────────────────────────────────────────
    SMCDA_VALIDATE_AGAINST_EPIDEMIC = True
    SMCDA_SENSITIVITY_DELTA = 0.20      # +20% one-at-a-time weight perturbation
    #
    # ── 6. FEEDING THE CLASSES BACK INTO THE DYNAMICS ────────────────────────
    # OFF by default, and it should stay off unless you mean it. The SMCDA score
    # is a PROSPECTIVE RISK RANKING built for surveillance allocation. Turning
    # this on makes it an INPUT to the epidemic, which (a) changes the simulated
    # spread and (b) destroys the independence that makes the validation above
    # meaningful — the map would then be partly validating itself.
    USE_SMCDA_CLASS_IN_TRANSPORT = False
    SMCDA_CLASS_MOBILITY = [0.05, 0.25, 0.50, 0.75, 1.00]   # one per class, low->high
    SMCDA_CLASS_COLORS = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
    # Draw the Korea-domain classifications on the same white-canvas peninsula
    # panel as the NK -> corridor gates -> SK progression figure.
    SHOW_SMCDA_KOREA_PANEL = True
    #
    # ── v9.8: HOW FAR THE KOREA PANEL IS CROPPED ─────────────────────────────
    # Drawing an SK-only classification on the whole-peninsula frame wastes
    # roughly half the figure on grey North Korean land and shrinks South Korea
    # to the point where county-level differences are not readable.
    #   "auto"       crop to South Korea whenever the domain being drawn
    #                contains no North Korean cells; keep the full peninsula
    #                otherwise. RECOMMENDED.
    #   "sk"         always crop to South Korea (the MDL and the corridor gates
    #                stay in frame along the top edge).
    #   "peninsula"  always draw the full peninsula — the v9.7 behaviour.
    SMCDA_KOREA_PANEL_CROP = "auto"
    # Cells of padding around the cropped area, so the coastline is not clipped
    # flush against the axes.
    SMCDA_KOREA_PANEL_CROP_MARGIN = 3
    # ══ STUDY AREA + OUTBREAK-POINT SEEDING ══════════════════════════════════
    STUDY_AREA_SHP_HINTS = ["studyarea", "study_area", "study area", "studysite"]
    SHOW_STUDY_AREA_OUTLINE = True
    EXPAND_GRID_TO_STUDY_AREA = True
    STUDY_AREA_MARGIN_CELLS = 4
    MAX_GRID_CELLS = 700_000
    SANITY_EXPAND_MAX_KM = 2500
    STUDY_AREA_DEFINES_DOMAIN = True
    CLIP_TO_STUDY_AREA = False
    STUDY_AREA_LINE_COLOR = "#111111"
    STUDY_AREA_LINE_WIDTH = 1.8
    STUDY_AREA_LINE_ALPHA = 0.95
    STUDY_AREA_MIN_CELLS = 50
    # ══ REGIONAL LAND BASEMAP (show the neighbouring countries) ══════════════
    SHOW_REGIONAL_BASEMAP = True
    BASEMAP_SHP_HINTS = ["ne_10m_admin", "ne_50m_admin", "ne_10m_land", "ne_50m_land",
                         "admin_0", "admin0", "countries", "country", "gadm", "world",
                         "land", "coast", "cntry", "china", "asia"]
    EXPAND_GRID_TO_BASEMAP = True
    EXPAND_GRID_TO_COVARIATES = True
    BASEMAP_LAND_COLOR = "#ece9e3"
    BASEMAP_LAND_ALPHA = 0.9
    BASEMAP_LINE_COLOR = "#b4ada1"
    BASEMAP_LINE_WIDTH = 0.6
    BASEMAP_LINE_ALPHA = 0.85
    BASEMAP_DEFINES_DOMAIN = False
    SEED_FROM_OUTBREAK_SHP = True
    OUTBREAK_SHP_HINTS = ["ppr_occ", "occ_unique", "occurrence", "occ",
                          "outbreak", "case", "seed", "index", "point"]
    OUTBREAK_SEED_WEIGHT = 1.0
    SHOW_OUTBREAK_SEED_MARKERS = True
    OUTBREAK_MARKER_COLOR = "#d81b60"
    OUTBREAK_MARKER_SIZE = 55
    # ══ PPR MaxEnt SUITABILITY MAP ═══════════════════════════════════════════
    MAXENT_TIF_HINTS = ["maxent", "ppr_maxent", "ppr_final", "final_map", "ndvi"]
    SHOW_MAXENT_MAP = True
    USE_MAXENT_IN_HABITAT = True
    MAXENT_HABITAT_BLEND = 0.50
    MAXENT_CMAP = ["#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"]
    SHOW_MAXENT_UNDERLAY = True
    MAXENT_UNDERLAY_ALPHA = 0.80
    # ══ BASE RASTER / MODEL DOMAIN ═══════════════════════════════════════════
    #   "gorani" — the Korea gorani donut (peninsula scale, ~5 km cells at
    #              FACTOR=10). USE THIS for corridor-level AND COUNTY-LEVEL
    #              conclusions: the gates are 30-60 km apart and South-Korean
    #              counties are ~500-1500 km², so both are properly resolved.
    #   "maxent" — the PPR MaxEnt raster (China-wide). Continental context, but
    #              the cells are coarse; adjacent gates can merge and the small
    #              metropolitan counties fall back to centroid sampling in the
    #              v9.3 CSV (each such row is flagged sampling="centroid").
    BASE_RASTER = "maxent"
    BASE_TARGET_MAX_PX = 750
    SINGLE_SPREAD_COLORMAP = False
    INITIAL_SUSCEPTIBLE = 1.0
    INTRODUCTION_INTENSITY = 0.55
    INFECTED_FLOOR = 0.005
    KERNEL_LOCAL_AGG_RADIUS_CELLS = 3
    BETA_WILDLIFE_TO_LIVESTOCK = 0.30
    BETA_LIVESTOCK_TO_LIVESTOCK = 0.20
    # ── LOCAL LANDSCAPE-EPIDEMIOLOGY CALIBRATION (Han et al. 2022) ───────────
    HAN_ROAD_REDUCTION = 0.362      # delta: crossing a major road -> ~36% passes
    HAN_RIVER_REDUCTION = 0.379     # nu:    crossing a river -> ~38% passes (boar!)
    HAN_FENCE_REDUCTION = 0.531     # phi:   fence pass-through (now the GATE CEILING)
    HAN_DIFF_HABITAT_REDUCTION = 0.624
    HAN_ALPHA_DECAY_PER_KM = 0.650
    HAN_APPLY_DMZ_FENCE = True
    HAN_APPLY_ROAD_BARRIER = True
    HAN_ROAD_RESISTANCE_ADD = 0.339
    HAN_APPLY_RIVER_BARRIER = False
    HAN_RIVER_RESISTANCE_ADD = 0.323
    HAN_LOG_DECAY_COMPARISON = True
    LOG_SIR_BALANCE = True
    # ── PERCENT-SPREAD BANDING + GORANI SPREADER DOTS ────────────────────────
    SPREAD_PCT_BANDS = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0]
    SHOW_PERCENT_BANDS = True
    # v9.2: the old first band ("#edf8e9") was invisible against the new white
    # South-Korea canvas; darkened so a 1-20% arrival actually reads on the map.
    SPREAD_BAND_COLORS_SK = ["#c7e9c0", "#74c476", "#31a354", "#2b8cbe", "#08519c"]
    SPREAD_BAND_COLORS_NK = ["#fff5eb", "#fdd0a2", "#fd8d3c", "#e6550d", "#a63603"]
    SHOW_GORANI_DOTS = True
    GORANI_DENSITY_WEIGHTS = {"gorani": 0.60, "water": 0.20,
                              "low_elev": 0.12, "low_slope": 0.08}
    GORANI_DOT_QUANTILE = 0.90
    GORANI_DOT_REQUIRE_SPREAD = False
    GORANI_DOT_SPREAD_FLOOR = 0.05
    GORANI_DOT_MAX = 2000
    GORANI_DOT_SIZE = 7
    GORANI_DOT_COLOR = "#ffe400"
    SHOW_VECTOR_SIG = True
    SHOW_VECTOR_WATER = True
    SIG_LINE_COLOR = "#7a7a7a"
    SIG_LINE_WIDTH = 0.3
    SIG_LINE_ALPHA = 0.6
    WATER_LINE_COLOR = "#cfe0f3"
    WATER_LINE_WIDTH = 0.20
    WATER_LINE_ALPHA = 0.5
    MAP_BACKGROUND = "white"
    SHOW_RIVER_RASTER = False
    SHOW_WATER_FILL = False
    SHOW_WETLAND_FILL = False
    LAND_TINT_ALPHA = 0.06
    USE_COMPOSITE_HABITAT_FOR_SPREAD = True
    HABITAT_BLEND = 0.85
    SHOW_RISK_UNDERLAY = True
    SHOW_RISK_UNDERLAY_NK = True
    SHOW_RISK_UNDERLAY_SK = False
    RISK_UNDERLAY_ALPHA = 0.45
    BAND_MIN_PCT_VISIBLE = 10.0
    BAND_ALPHA = 0.58
    BOUNDARIES_OVER_RISK = True
# ══════════════════════════════════════════════════════════════════════════════
# NEURAL NETWORK COMPONENTS (v7.0)
# ══════════════════════════════════════════════════════════════════════════════
def build_neural_components(device="cpu"):
    """Build SpatialKernelNet and ODEFunc; returns (None, None, None)
    if PyTorch is unavailable."""
    if not TORCH_AVAILABLE:
        return None, None, None
    class SpatialKernelNet(nn.Module):
        """Learned spatial kernel modulator.
        Input: (1, N_KERNEL_FEATURES, NR, NC) feature tensor.
        Output: (NR, NC) hybrid kernel weight = K_physics*(1 + alpha*K_learned).
        """
        def __init__(self, in_channels, hidden, n_layers):
            super().__init__()
            layers = []
            ch = in_channels
            for _ in range(n_layers - 1):
                layers += [nn.Conv2d(ch, hidden, 3, padding=1),
                           nn.BatchNorm2d(hidden), nn.ReLU(inplace=True)]
                ch = hidden
            layers.append(nn.Conv2d(ch, 1, 1))
            self.net = nn.Sequential(*layers)
            self.alpha = nn.Parameter(torch.tensor(0.1))
        def forward(self, features, K_physics_np):
            K_phys = torch.from_numpy(K_physics_np).float().to(features.device)
            K_phys = K_phys.unsqueeze(0).unsqueeze(0)
            K_learned = torch.sigmoid(self.net(features))
            K_hybrid = K_phys * (1.0 + self.alpha * K_learned)
            return K_hybrid.squeeze().detach().cpu().numpy().astype(np.float32)
    class ODEFunc(nn.Module):
        """Neural ODE: dZ/dt = f_theta(Z, context); Z = [S,I,R] flattened."""
        def __init__(self, hidden, gamma_base=0.15):
            super().__init__()
            self.gamma_base = gamma_base
            self.lambda_spatial = None  # set before each odeint call
            self.J_m_scalar = 1.0
            self.mlp = nn.Sequential(
                nn.Linear(4, hidden), nn.Tanh(),
                nn.Linear(hidden, hidden), nn.Tanh(),
                nn.Linear(hidden, 2))  # -> [delta_inf, delta_rem]
        def forward(self, t, Z):
            N = Z.shape[0] // 3
            S, I, R = Z[:N], Z[N:2*N], Z[2*N:]
            lam = self.lambda_spatial
            inp = torch.stack([S, I, R, lam], dim=1)
            delta = self.mlp(inp)
            d_inf = torch.tanh(delta[:, 0]) * 0.5
            d_rem = torch.tanh(delta[:, 1]) * 0.5
            eff_inf = torch.clamp(self.J_m_scalar + d_inf, min=0.0)
            eff_rem = torch.clamp(self.gamma_base + d_rem, min=0.0)
            new_inf = eff_inf * lam * S
            rem = eff_rem * I
            return torch.cat([-new_inf, new_inf - rem, rem])
    kernel_net = SpatialKernelNet(
        Config.N_KERNEL_FEATURES, Config.KERNEL_NET_HIDDEN, Config.KERNEL_NET_LAYERS
    ).to(device)
    ode_func = ODEFunc(Config.NODE_HIDDEN, Config.GAMMA_RECOVERY).to(device)
    return kernel_net, ode_func, device
def compute_kernel_features(GOR_nn, ELEV_n, SLOPE_N, livestock_n, RIV_nn,
                            CONDUCTANCE, dist_km, I_proxy, device):
    """Stack 8 feature maps -> (1, 8, NR, NC) tensor for SpatialKernelNet."""
    dist_norm = np.clip(dist_km / (np.nanmax(dist_km) + 1e-9), 0, 1).astype(np.float32)
    dist_norm[~np.isfinite(dist_norm)] = 1.0
    stack = np.stack([GOR_nn, np.nan_to_num(ELEV_n, nan=0.5), SLOPE_N,
                      np.nan_to_num(livestock_n, nan=0.0), RIV_nn, CONDUCTANCE,
                      dist_norm, I_proxy], axis=0).astype(np.float32)
    return torch.from_numpy(stack).unsqueeze(0).to(device)
def train_neural_components(kernel_net, ode_func, outbreak_df,
                            GT5, NR, NC,
                            S_init_combined, I_init_combined, R_init_combined,
                            monthly_lambda_cache, feature_cache,
                            device, logger):
    """Train SpatialKernelNet + ODEFunc against observed outbreak locations."""
    if not TORCHDIFFEQ_AVAILABLE:
        logger.warning("  ⚠ Training skipped: torchdiffeq not installed. "
                       "Run: pip install torchdiffeq --break-system-packages")
        return
    logger.info(f"  Training: {Config.N_TRAIN_EPOCHS} epochs, "
                f"lr={Config.LEARNING_RATE}, solver={Config.NODE_SOLVER}")
    optimizer = torch.optim.Adam(
        list(kernel_net.parameters()) + list(ode_func.parameters()),
        lr=Config.LEARNING_RATE, weight_decay=1e-5)
    outbreak_mask = torch.zeros(12, NR, NC, device=device)
    if outbreak_df is not None and len(outbreak_df) > 0:
        for _, row in outbreak_df.iterrows():
            try:
                col = int((float(row["lon"]) - GT5[0]) / GT5[1])
                r = int((GT5[3] - float(row["lat"])) / abs(GT5[5]))
                m_idx = int(row.get("month", 1)) - 1
                if 0 <= r < NR and 0 <= col < NC and 0 <= m_idx < 12:
                    outbreak_mask[m_idx, r, col] = 1.0
            except Exception:
                pass
        logger.info(f"  Outbreak cells: {int(outbreak_mask.sum().item())}")
    else:
        logger.warning("  ⚠ No outbreak data -- training loss will be unguided")
    N = NR * NC
    for epoch in range(Config.N_TRAIN_EPOCHS):
        optimizer.zero_grad()
        S = torch.from_numpy(S_init_combined.flatten()).float().to(device)
        I = torch.from_numpy(I_init_combined.flatten()).float().to(device)
        R = torch.from_numpy(R_init_combined.flatten()).float().to(device)
        total_loss = torch.zeros(1, device=device, requires_grad=False)
        for m in range(12):
            lam_combined, feats, J_m = monthly_lambda_cache[m], feature_cache[m], float(Config.MONTHLY_BETA[m])
            lam_t = torch.from_numpy(lam_combined.flatten()).float().to(device)
            ode_func.lambda_spatial = lam_t
            ode_func.J_m_scalar = J_m
            Z0 = torch.cat([S, I, R])
            try:
                Z1 = odeint(ode_func, Z0,
                            torch.tensor([0.0, 1.0], device=device),
                            method=Config.NODE_SOLVER,
                            rtol=Config.NODE_RTOL, atol=Config.NODE_ATOL)[-1]
            except Exception as e:
                logger.warning(f"  ⚠ ODE fail epoch={epoch} m={m}: {e}")
                continue
            S, I, R = Z1[:N].detach(), Z1[N:2*N], Z1[2*N:].detach()
            I_map = I.reshape(NR, NC).clamp(1e-8, 1.0)
            mask_m = outbreak_mask[m]
            if mask_m.sum() > 0:
                total_loss = total_loss + (-torch.log(I_map[mask_m.bool()]).mean())
        total_loss.backward()
        optimizer.step()
        if (epoch + 1) % 10 == 0 or epoch == 0:
            logger.info(f"  Epoch {epoch+1:3d}/{Config.N_TRAIN_EPOCHS}: loss={total_loss.item():.4f}")
    logger.info("  ✓ Training complete")
# ══════════════════════════════════════════════════════════════════════════════
# LOAD SPATIAL DATA
# ══════════════════════════════════════════════════════════════════════════════
def main():
    logger.info("\n" + "=" * 100)
    logger.info("PPR v9.3 — DYNAMIC RAD EPIDEMIC + CORRIDOR-GATED MDL FENCE "
                "+ PLAIN-WHITE SK PANELS + SOUTH-KOREA COUNTY CSV EXPORT")
    logger.info("=" * 100)
    # ── 1. Gorani habitat TIF ───────────────────────────────────────────────
    gorani_paths = [
        Config.CWD / "gorani_habitat_suitability_cloglog_donut_500m_10km.tif",
        Config.DATA_ROOT / "gorani_habitat_suitability_cloglog_donut_500m_10km.tif",
        Config.ENV_DIR / "gorani_habitat_suitability_cloglog_donut_500m_10km.tif",
    ]
    gorani_paths += list(Config.CWD.rglob("gorani*cloglog*donut*.tif")) if Config.CWD.exists() else []
    gorani_paths += list(Config.DATA_ROOT.rglob("gorani*cloglog*donut*.tif")) if Config.DATA_ROOT.exists() else []
    GORANI_TIF = None
    for gp in gorani_paths:
        if Path(str(gp)).exists():
            GORANI_TIF = Path(str(gp))
            logger.info(f"✓ FOUND GORANI TIF: {GORANI_TIF}")
            break
    # ── Choose the BASE raster: it defines the grid, CRS and domain ─────────
    def _find_maxent_tif():
        if getattr(Config, "MAXENT_TIF_PATH", None):        # v9.10 item 7: explicit MaxEnt map
            _p = Path(str(Config.MAXENT_TIF_PATH))
            if not _p.exists():
                raise SystemExit(f"v9.10: MAXENT_TIF_PATH not found: {_p}")
            logger.info(f"v9.10: MaxEnt map set explicitly: {_p}")
            return _p
        for root in [Config.CWD, Config.DATA_ROOT, Config.ENV_DIR]:
            if not root.exists():
                continue
            for t in root.rglob("*.tif"):
                if t.name.startswith("PPR_v9") or "PPR_v9.10_rasters" in t.parts:   # v9.10: own outputs
                    continue
                nl = t.name.lower()
                if any(h.lower() in nl for h in Config.MAXENT_TIF_HINTS):
                    return t
        return None
    MAXENT_TIF = _find_maxent_tif()
    base_is_maxent = (Config.BASE_RASTER == "maxent" and MAXENT_TIF is not None)
    if Config.BASE_RASTER == "maxent" and MAXENT_TIF is None:
        logger.warning("  ⚠ BASE_RASTER='maxent' but no MaxEnt TIF found — falling back to gorani base.")
    MAXENT_BASE = None
    if base_is_maxent:
        logger.info(f"✓ BASE raster = MaxEnt: {MAXENT_TIF.name} (grid + CRS taken from this raster)")
        dsb = gdal.Open(str(MAXENT_TIF))
        gtb = dsb.GetGeoTransform(); WKT = dsb.GetProjection()
        ndb = dsb.GetRasterBand(1).GetNoDataValue()
        rawb = dsb.GetRasterBand(1).ReadAsArray().astype(np.float64); dsb = None
        if ndb is not None:
            rawb[np.isclose(rawb, ndb, rtol=0, atol=1e-3)] = np.nan
        rawb[rawb <= -3.4e38 * 0.01] = np.nan
        rawb[~np.isfinite(rawb)] = np.nan
        hb, wb = rawb.shape
        FACTOR = max(1, int(math.ceil(max(hb, wb) / float(Config.BASE_TARGET_MAX_PX))))
        while (hb // FACTOR) * (wb // FACTOR) > Config.MAX_GRID_CELLS:
            FACTOR += 1
        rb, cb = hb // FACTOR, wb // FACTOR
        blkb = rawb[:rb * FACTOR, :cb * FACTOR].reshape(rb, FACTOR, cb, FACTOR)
        with np.errstate(invalid="ignore"):
            MAXENT_BASE = np.nanmean(blkb, axis=(1, 3)).astype(np.float32)
        _mv = MAXENT_BASE[np.isfinite(MAXENT_BASE)]
        if _mv.size and float(np.nanmax(MAXENT_BASE)) > float(np.nanmin(MAXENT_BASE)):
            _mn, _mx = float(np.nanmin(MAXENT_BASE)), float(np.nanmax(MAXENT_BASE))
            MAXENT_BASE = np.clip((MAXENT_BASE - _mn) / (_mx - _mn + 1e-9), 0, 1).astype(np.float32)
        GT5 = (gtb[0], gtb[1] * FACTOR, 0., gtb[3], 0., gtb[5] * FACTOR)
        NR, NC = MAXENT_BASE.shape; SH = (NR, NC)
        PW = abs(GT5[1]); PH = abs(GT5[5])
        if GORANI_TIF is not None:
            arr5 = warp_raster_to_grid(GORANI_TIF, GT5, NR, NC, WKT)
            arr5 = np.clip(arr5, 0, 1).astype(np.float32)
        else:
            arr5 = np.full(SH, np.nan, np.float32)
            logger.warning("  ⚠ No gorani TIF — gorani layer empty; habitat = MaxEnt only.")
        VALID = np.isfinite(MAXENT_BASE) | np.isfinite(arr5)
        logger.info(f"✓ Grid (MaxEnt base): {NC}×{NR} cells, FACTOR={FACTOR}; "
                    f"MaxEnt-valid={int(np.isfinite(MAXENT_BASE).sum()):,}; gorani/Korea warped in "
                    f"={int(np.isfinite(arr5).sum()):,} cells; domain(KP)={int(VALID.sum()):,}.")
    else:
        if GORANI_TIF is None:
            logger.error(f"✗ Gorani TIF NOT found in: {gorani_paths}")
            raise FileNotFoundError("gorani_habitat_suitability_cloglog_donut_500m_10km.tif")
        FACTOR = Config.FACTOR
        ds = gdal.Open(str(GORANI_TIF))
        gt_orig = ds.GetGeoTransform()
        WKT = ds.GetProjection()
        nd_ = ds.GetRasterBand(1).GetNoDataValue()
        raw = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
        ds = None
        G_X0, G_Y0 = gt_orig[0], gt_orig[3]
        G_PW, G_PH = abs(gt_orig[1]), abs(gt_orig[5])
        if nd_ is not None:
            raw[np.isclose(raw, nd_, rtol=0, atol=1e-3)] = np.nan
        raw[raw <= -3.4e38 * 0.01] = np.nan
        raw[raw < 0] = np.nan
        raw = np.clip(raw, 0, 1)
        raw[~np.isfinite(raw)] = np.nan
        h, w = raw.shape
        r5, c5 = h // FACTOR, w // FACTOR
        blk = raw[:r5 * FACTOR, :c5 * FACTOR].reshape(r5, FACTOR, c5, FACTOR)
        arr5 = np.nanmean(blk, axis=(1, 3)).astype(np.float32)
        vb = np.isfinite(raw[:r5 * FACTOR, :c5 * FACTOR]).reshape(r5, FACTOR, c5, FACTOR)
        arr5[~vb.any(axis=(1, 3))] = np.nan
        GT5 = (G_X0, G_PW * FACTOR, 0., G_Y0, 0., -(G_PH * FACTOR))
        NR, NC = arr5.shape
        SH = (NR, NC)
        PW = abs(GT5[1])
        PH = abs(GT5[5])
        VALID = np.isfinite(arr5)
        logger.info(f"✓ Grid: {NC}×{NR} cells ({NC * PW / 1e3:.0f}×{NR * PH / 1e3:.0f} km)")
    # Single (whole-domain) spread colormap — forced on at continental scale.
    SINGLE_MODE = bool(Config.SINGLE_SPREAD_COLORMAP or base_is_maxent)
    # ── 1b. Expand grid to cover the study area / basemap / covariates ───────
    def _find_first_shp(hints):
        for hint in hints:
            for root in [Config.CWD, Config.DATA_ROOT, Config.ENV_DIR]:
                if root.exists():
                    hits = list(root.rglob(f"*{hint}*.shp"))
                    if hits:
                        return hits[0]
        return None
    def _shp_bounds_wkt(shp_path):
        try:
            ds_b = _ogr.Open(str(shp_path))
            if ds_b is None:
                return None
            lyr_b = ds_b.GetLayer()
            ext = lyr_b.GetExtent()                 # (xmin,xmax,ymin,ymax) layer CRS
            src_b = lyr_b.GetSpatialRef()
            dst_b = _osr.SpatialReference(); dst_b.ImportFromWkt(WKT)
            if src_b and not src_b.IsSame(dst_b):
                tfm = _osr.CoordinateTransformation(src_b, dst_b)
                xs, ys = [], []
                for (x, y) in [(ext[0], ext[2]), (ext[1], ext[2]), (ext[0], ext[3]), (ext[1], ext[3])]:
                    p = tfm.TransformPoint(x, y); xs.append(p[0]); ys.append(p[1])
                b = (min(xs), min(ys), max(xs), max(ys))
            else:
                b = (ext[0], ext[2], ext[1], ext[3])
            ds_b = None
            return b
        except Exception as e:
            logger.warning(f"  ⚠ bounds unreadable for {shp_path}: {e}")
            return None
    def _raster_bounds_wkt(rpath):
        try:
            if not Path(str(rpath)).exists():
                return None
            dsr = gdal.Open(str(rpath))
            gtr = dsr.GetGeoTransform(); wr = dsr.RasterXSize; hr = dsr.RasterYSize
            xr0, yr0 = gtr[0], gtr[3]; xr1 = gtr[0] + wr * gtr[1]; yr1 = gtr[3] + hr * gtr[5]
            src_r = _osr.SpatialReference(); src_r.ImportFromWkt(dsr.GetProjection() or WKT)
            dst_r = _osr.SpatialReference(); dst_r.ImportFromWkt(WKT)
            dsr = None
            xmin, xmax = min(xr0, xr1), max(xr0, xr1); ymin, ymax = min(yr0, yr1), max(yr0, yr1)
            if src_r and not src_r.IsSame(dst_r):
                tfm = _osr.CoordinateTransformation(src_r, dst_r)
                xs, ys = [], []
                for (x, y) in [(xmin, ymin), (xmax, ymin), (xmin, ymax), (xmax, ymax)]:
                    p = tfm.TransformPoint(x, y); xs.append(p[0]); ys.append(p[1])
                return (min(xs), min(ys), max(xs), max(ys))
            return (xmin, ymin, xmax, ymax)
        except Exception:
            return None
    if (not base_is_maxent) and (Config.EXPAND_GRID_TO_STUDY_AREA or Config.EXPAND_GRID_TO_BASEMAP
                                 or Config.EXPAND_GRID_TO_COVARIATES):
        _sources = []
        if Config.EXPAND_GRID_TO_STUDY_AREA:
            _s = _find_first_shp(Config.STUDY_AREA_SHP_HINTS)
            if _s is not None:
                _sources.append((f"study:{_s.name}", _shp_bounds_wkt(_s)))
        if Config.EXPAND_GRID_TO_BASEMAP:
            _bm = _find_first_shp(Config.BASEMAP_SHP_HINTS)
            if _bm is not None:
                _sources.append((f"basemap:{_bm.name}", _shp_bounds_wkt(_bm)))
        if Config.EXPAND_GRID_TO_COVARIATES and Config.ENV_DIR.exists():
            for _rn in ["wc2.1_30s_elev_studyarea_1km_5179.tif", "GLW4GTS_1km_5179.tif",
                        "wc2.1_30s_bio_12_studyarea_1km_5179.tif"]:
                _rb = _raster_bounds_wkt(Config.ENV_DIR / _rn)
                if _rb is not None:
                    _sources.append((f"raster:{_rn}", _rb))
        PWm = PW
        g_xmin, g_ymax = GT5[0], GT5[3]
        g_xmax = GT5[0] + NC * PWm; g_ymin = GT5[3] - NR * PWm
        _sanity = Config.SANITY_EXPAND_MAX_KM * 1000.0
        t_xmin, t_xmax, t_ymin, t_ymax = g_xmin, g_xmax, g_ymin, g_ymax
        _used = []
        for _name, _b in _sources:
            if _b is None:
                continue
            if (_b[0] < g_xmin - _sanity or _b[2] > g_xmax + _sanity or
                    _b[1] < g_ymin - _sanity or _b[3] > g_ymax + _sanity):
                logger.warning(f"  ⚠ ignoring {_name} for grid extent — reaches absurdly far "
                               f"beyond the gorani grid (bad coordinates/CRS? span "
                               f"x[{_b[0]:.0f},{_b[2]:.0f}] y[{_b[1]:.0f},{_b[3]:.0f}]).")
                continue
            if (_b[0] < g_xmin - 1 or _b[2] > g_xmax + 1 or _b[1] < g_ymin - 1 or _b[3] > g_ymax + 1):
                _used.append(_name)
            t_xmin = min(t_xmin, _b[0]); t_xmax = max(t_xmax, _b[2])
            t_ymin = min(t_ymin, _b[1]); t_ymax = max(t_ymax, _b[3])
        m_ = Config.STUDY_AREA_MARGIN_CELLS * PWm
        if _used:
            t_xmin -= m_; t_xmax += m_; t_ymin -= m_; t_ymax += m_
        pad_left = max(0, int(math.ceil((g_xmin - t_xmin) / PWm)))
        pad_right = max(0, int(math.ceil((t_xmax - g_xmax) / PWm)))
        pad_top = max(0, int(math.ceil((t_ymax - g_ymax) / PWm)))
        pad_bottom = max(0, int(math.ceil((g_ymin - t_ymin) / PWm)))
        new_NR = NR + pad_top + pad_bottom
        new_NC = NC + pad_left + pad_right
        if (pad_left + pad_right + pad_top + pad_bottom) == 0:
            logger.info("  Grid-expand: no layer extends beyond the gorani grid — no change. "
                        "(To show other countries, add a land/countries shapefile — see "
                        "BASEMAP_SHP_HINTS — or a covariate raster that covers them.)")
        elif new_NR * new_NC > Config.MAX_GRID_CELLS:
            logger.warning(f"  ⚠ Grid-expand skipped: {new_NC}×{new_NR} ({new_NR*new_NC:,} cells) "
                           f"exceeds MAX_GRID_CELLS={Config.MAX_GRID_CELLS:,}. Raise the cap or FACTOR.")
        else:
            arr5 = np.pad(arr5, ((pad_top, pad_bottom), (pad_left, pad_right)),
                          mode="constant", constant_values=np.nan)
            NR, NC = arr5.shape; SH = (NR, NC)
            GT5 = (g_xmin - pad_left * PWm, PWm, 0., g_ymax + pad_top * PWm, 0., -PWm)
            PW = abs(GT5[1]); PH = abs(GT5[5])
            VALID = np.isfinite(arr5)
            logger.info(f"  ✓ Grid expanded to {NC}×{NR} cells "
                        f"({NC*PW/1e3:.0f}×{NR*PH/1e3:.0f} km); pad L{pad_left} R{pad_right} "
                        f"T{pad_top} B{pad_bottom}; driven by: {', '.join(_used)}")
    # ══ v9.10 ITEM 6: grid-resolution test (CLI --preset grid4km) ═══════════════
    # The transport coefficients are in CELL units. On a finer grid they are
    # rescaled to keep the physical rates of the reference grid: D ∝ Δx⁻²,
    # v_max ∝ Δx⁻¹, advection gain ∝ Δx⁻² (the habitat gradient per cell shrinks
    # with Δx), and the fence band keeps (to the nearest odd number of cells) its width in km.
    if Config.TRANSPORT_REFERENCE_CELL_KM:
        _ref_km = Config.TRANSPORT_REFERENCE_CELL_KM
        if str(_ref_km).lower() == "auto":
            _ref_km = None
            if base_is_maxent:
                _f_ref = max(1, int(math.ceil(max(hb, wb) / float(Config.REFERENCE_TARGET_MAX_PX))))
                while (hb // _f_ref) * (wb // _f_ref) > int(Config.REFERENCE_MAX_GRID_CELLS):
                    _f_ref += 1
                _ref_km = abs(gtb[1]) * _f_ref / 1000.0
            else:
                logger.warning("  ⚠ v9.10: TRANSPORT_REFERENCE_CELL_KM='auto' needs the MaxEnt base "
                               "grid; give the reference cell size in km instead. No rescaling.")
        if _ref_km:
            _s_grid = float(_ref_km) / (PW / 1000.0)
            if abs(_s_grid - 1.0) > 1e-6:
                _w_ref = 2 * int(Config.MDL_DILATE_ITER) + 1
                Config.RAD_DIFF_MAX *= _s_grid ** 2
                Config.RAD_DIFF_MIN *= _s_grid ** 2
                Config.RAD_ADV_MAX *= _s_grid
                Config.RAD_ADV_GAIN *= _s_grid ** 2
                Config.MDL_DILATE_ITER = max(1, int(round((_w_ref * _s_grid - 1) / 2.0)))
                logger.info(f"  v9.10 grid test: cells {PW / 1000.0:.2f} km vs reference {float(_ref_km):.2f} km "
                            f"(scale {_s_grid:.2f}) -> D_max={Config.RAD_DIFF_MAX:.3f}, "
                            f"D_min={Config.RAD_DIFF_MIN:.3f}, v_max={Config.RAD_ADV_MAX:.3f}, "
                            f"a_g={Config.RAD_ADV_GAIN:.3f} (cell units); fence band "
                            f"{2 * Config.MDL_DILATE_ITER + 1} cells "
                            f"≈ {(2 * Config.MDL_DILATE_ITER + 1) * PW / 1000.0:.0f} km.")
    # ── 2. MDL boundary ─────────────────────────────────────────────────────
    MDL_PATHS = [Config.CWD / "MDL.shp", Config.DATA_ROOT / "MDL.shp"]
    MDL_PATHS += list(Config.CWD.rglob("mdl*.shp")) if Config.CWD.exists() else []
    MDL_PATHS += list(Config.DATA_ROOT.rglob("mdl*.shp")) if Config.DATA_ROOT.exists() else []
    MDL_FILE = None
    for mp in MDL_PATHS:
        if Path(str(mp)).exists():
            MDL_FILE = Path(str(mp))
            logger.info(f"✓ FOUND MDL: {MDL_FILE}")
            break
    if MDL_FILE is None:
        logger.error(f"✗ MDL.shp NOT found in: {MDL_PATHS}")
        raise FileNotFoundError("MDL.shp")
    def rasterize_boundary(shp_path, tgt_wkt, tgt_gt, nr, nc, dilate=True, iters=1):
        """Burn a shapefile onto the model grid. `iters` controls how many cells
        wide the burned band becomes (v9.1: a thicker MDL band means the front
        must traverse several attenuated fence cells, which lengthens the
        crossing without touching the pass-through factor itself)."""
        mem_drv = gdal.GetDriverByName('MEM')
        out_ds = mem_drv.Create('', nc, nr, 1, gdal.GDT_Byte)
        out_ds.SetGeoTransform(tgt_gt)
        out_ds.SetProjection(tgt_wkt)
        out_ds.GetRasterBand(1).Fill(0)
        xmin = tgt_gt[0]
        xmax = tgt_gt[0] + nc * tgt_gt[1]
        ymax = tgt_gt[3]
        ymin = tgt_gt[3] + nr * tgt_gt[5]
        try:
            vt_ds = gdal.VectorTranslate('', str(shp_path), format='Memory', dstSRS=tgt_wkt)
            if vt_ds:
                vt_lyr = vt_ds.GetLayer()
                vt_lyr.SetSpatialFilterRect(xmin, ymin, xmax, ymax)
                gdal.RasterizeLayer(out_ds, [1], vt_lyr, burn_values=[1])
                vt_ds = None
        except Exception as e:
            logger.warning(f"  ⚠ rasterize via VectorTranslate failed: {e}")
        arr = out_ds.GetRasterBand(1).ReadAsArray().astype(bool)
        out_ds = None
        if dilate and iters > 0:
            return binary_dilation(arr, iterations=int(iters)).astype(bool)
        return arr
    MDL_BOUNDARY = rasterize_boundary(MDL_FILE, WKT, GT5, NR, NC,
                                      dilate=True, iters=max(1, int(Config.MDL_DILATE_ITER)))
    logger.info(f"✓ MDL rasterized: {MDL_BOUNDARY.sum():,} cells "
                f"(dilation={Config.MDL_DILATE_ITER} → ~{2*Config.MDL_DILATE_ITER+1} cells "
                f"≈ {(2*Config.MDL_DILATE_ITER+1)*PW/1000.0:.0f} km of fence band to traverse)")
    boundary_rows = np.where(MDL_BOUNDARY.any(axis=1))[0]
    DMZ_ROW = int(np.median(boundary_rows)) if len(boundary_rows) else NR // 2
    RI = np.tile(np.arange(NR).reshape(-1, 1), (1, NC))
    NK = VALID & (RI < DMZ_ROW)
    SK = VALID & (RI >= DMZ_ROW)
    KP = NK | SK
    logger.info(f"✓ NK/SK partition: NK={NK.sum():,}, SK={SK.sum():,}, DMZ row={DMZ_ROW}")
    # ══ v9.1: DMZ CORRIDOR GATES ═════════════════════════════════════════════
    # Reproject the corridor table, snap each point onto the digitized MDL, and
    # build the per-cell fence permeability field the RAD engine consumes.
    CORRIDORS: List[Dict] = []
    MDL_FENCE_FIELD = None
    if Config.USE_CORRIDOR_CROSSINGS:
        logger.info("\n[DMZ CORRIDOR CROSSINGS — water-deer gates in the MDL line]")
        CORRIDORS = build_corridors(
            MDL_BOUNDARY, GT5, WKT, NR, NC, cell_km=PW / 1000.0, kp=VALID,
            snap_max_km=Config.CORRIDOR_SNAP_MAX_KM,
            gate_radius_km=Config.CORRIDOR_GATE_RADIUS_KM, logger=logger)
        MDL_FENCE_FIELD = build_mdl_fence_field(
            CORRIDORS, MDL_BOUNDARY, (NR, NC), cell_km=PW / 1000.0,
            gate_radius_km=Config.CORRIDOR_GATE_RADIUS_KM,
            max_pass=Config.CORRIDOR_MAX_PASS,
            sealed=Config.MDL_SEALED_FACTOR,
            slowdown=Config.MDL_CROSSING_SLOWDOWN, logger=logger)
        if MDL_FENCE_FIELD is not None:
            logger.info(f"  Crossing is now gate-limited: peak pass-through "
                        f"{Config.CORRIDOR_MAX_PASS:.3f}×weight×{Config.MDL_CROSSING_SLOWDOWN:.2f} "
                        f"at the gates vs {Config.MDL_SEALED_FACTOR:.3f} on the rest of the line "
                        f"(v9.0 used {Config.HAN_FENCE_REDUCTION:.3f} everywhere), "
                        f"seasonal gating={'ON' if Config.CORRIDOR_SEASONAL_GATING else 'OFF'}.")
    GOR_P = np.where(KP, arr5, np.nan).astype(np.float32)
    GOR_nn = np.nan_to_num(GOR_P, nan=0.0)
    # ── 3. Optional geo layers (water/wetland/county) ───────────────────────
    logger.info("\n[SEARCHING FOR OPTIONAL GEO LAYERS]")
    ALL_SHPS = []
    for root in [Config.CWD, Config.DATA_ROOT, Config.ENV_DIR]:
        if root.exists():
            ALL_SHPS += list(root.rglob("*.shp"))
    # ── DATA INVENTORY (diagnostic) ──────────────────────────────────────────
    def _extent_in_wkt(shp_path):
        try:
            ds_i = _ogr.Open(str(shp_path))
            if ds_i is None:
                return None, None, 0
            lyr_i = ds_i.GetLayer()
            n_i = lyr_i.GetFeatureCount()
            gname = _ogr.GeometryTypeToName(lyr_i.GetGeomType())
            ext = lyr_i.GetExtent()  # (xmin,xmax,ymin,ymax) in layer CRS
            src_i = lyr_i.GetSpatialRef()
            dst_i = _osr.SpatialReference(); dst_i.ImportFromWkt(WKT)
            if src_i and not src_i.IsSame(dst_i):
                tfm = _osr.CoordinateTransformation(src_i, dst_i)
                xs, ys = [], []
                for (x, y) in [(ext[0], ext[2]), (ext[1], ext[2]), (ext[0], ext[3]), (ext[1], ext[3])]:
                    p = tfm.TransformPoint(x, y); xs.append(p[0]); ys.append(p[1])
                ext_w = (min(xs), min(ys), max(xs), max(ys))
            else:
                ext_w = (ext[0], ext[2], ext[1], ext[3])
            ds_i = None
            return gname, ext_w, n_i
        except Exception as e:
            return f"<error {e}>", None, 0
    _g_xmin, _g_ymax = GT5[0], GT5[3]
    _g_xmax = GT5[0] + NC * PW; _g_ymin = GT5[3] - NR * PH
    logger.info("[DATA INVENTORY]")
    logger.info(f"  GORANI/grid extent (m, EPSG of raster): "
                f"x[{_g_xmin:.0f}, {_g_xmax:.0f}]  y[{_g_ymin:.0f}, {_g_ymax:.0f}]  ({NC}×{NR} cells)")
    for _sp in ALL_SHPS:
        _gn, _ew, _n = _extent_in_wkt(_sp)
        if _ew is None:
            logger.info(f"  · {_sp.name}: {_gn}, {_n} feats, extent unreadable")
            continue
        _beyond = (_ew[0] < _g_xmin - 1 or _ew[2] > _g_xmax + 1 or
                   _ew[1] < _g_ymin - 1 or _ew[3] > _g_ymax + 1)
        logger.info(f"  · {_sp.name}: {_gn}, {_n} feats, "
                    f"x[{_ew[0]:.0f}, {_ew[2]:.0f}] y[{_ew[1]:.0f}, {_ew[3]:.0f}]"
                    f"{'  <-- EXTENDS BEYOND gorani grid' if _beyond else ''}")
    def load_polys(name_substr, label):
        polys = []
        candidates = [Config.CWD / f"{name_substr}.shp", Config.DATA_ROOT / f"{name_substr}.shp",
                      Config.ENV_DIR / f"{name_substr}.shp"]
        candidates += [p for p in ALL_SHPS if name_substr in p.name.lower()]
        found = None
        for c in candidates:
            if Path(str(c)).exists():
                found = Path(str(c))
                break
        if found is None:
            logger.info(f"  {label}: not found (optional, skipping)")
            return polys
        try:
            ds_l = _ogr.Open(str(found))
            lyr = ds_l.GetLayer()
            src_srs = lyr.GetSpatialRef()
            dst_srs = _osr.SpatialReference()
            dst_srs.ImportFromWkt(WKT)
            tfm = (_osr.CoordinateTransformation(src_srs, dst_srs)
                   if src_srs and not src_srs.IsSame(dst_srs) else None)
            for feat in lyr:
                geom = feat.GetGeometryRef()
                if geom is None:
                    continue
                g = geom.Clone()
                if tfm:
                    g.Transform(tfm)
                gname = g.GetGeometryName().upper()
                parts = ([g.GetGeometryRef(i) for i in range(g.GetGeometryCount())]
                         if 'MULTI' in gname else [g])
                for part in parts:
                    if part.GetGeometryCount() == 0:
                        continue
                    ring = part.GetGeometryRef(0)
                    pts = ring.GetPoints()
                    if pts and len(pts) > 2:
                        xs = np.array([pt[0] for pt in pts], dtype=np.float32)
                        ys = np.array([pt[1] for pt in pts], dtype=np.float32)
                        cs = (xs - GT5[0]) / GT5[1]
                        rs = (GT5[3] - ys) / PH
                        polys.append((cs, rs))
            ds_l = None
            logger.info(f"  ✓ {label}: loaded {len(polys)} polygons from {found.name}")
        except Exception as e:
            logger.warning(f"  ⚠ {label}: error loading {found}: {type(e).__name__}: {e}")
        return polys
    SIG_POLYS = load_polys("sig", "Counties (sig.shp)")
    WATER_POLYS = load_polys("water", "Water (water.shp)")
    WETLAND_POLYS = load_polys("wetland", "Wetland (wetland.shp)")
    def find_shp(name_substr):
        cands = [Config.CWD / f"{name_substr}.shp", Config.DATA_ROOT / f"{name_substr}.shp",
                 Config.ENV_DIR / f"{name_substr}.shp"]
        cands += [p for p in ALL_SHPS if name_substr in p.name.lower()]
        for c in cands:
            if Path(str(c)).exists():
                return Path(str(c))
        return None
    ROAD_SHP = find_shp("road") if Config.HAN_APPLY_ROAD_BARRIER else None
    if Config.HAN_APPLY_ROAD_BARRIER:
        logger.info(f"  Roads (road*.shp): "
                    f"{'found ' + ROAD_SHP.name if ROAD_SHP else 'not found (road barrier will be skipped)'}")
    # ══ v9.3: SOUTH-KOREA COUNTY ZONES ═══════════════════════════════════════
    # Build the county × model-cell AREA WEIGHTS once, here, so both the Korea
    # (cross-MDL) leg and the continental run can be tabulated against the same
    # administrative units later.
    SK_COUNTY_ZONES: List[Dict] = []
    SK_ZONE_RASTER = np.zeros(SH, np.int32)
    SK_COUNTY_SHP = None
    SK_ZONE_META = {}
    SK_ZONE_FAIL = None          # (reason, message) when the build fails
    if Config.EXPORT_SK_COUNTY_CSV:
        logger.info("\n[SOUTH-KOREA COUNTY (SIGUNGU) ZONES — for the CSV export]")
        # Log EVERY candidate the hints resolve to, not just the winner. When the
        # export fails, the usual cause is that a hint matched the wrong file
        # (note that "sig" is a substring of "de-sig-n"), and that is invisible
        # unless the alternatives are printed.
        _extra_shps = []
        for _d in (Config.SK_COUNTY_EXTRA_SEARCH_DIRS or []):
            _dp = Path(str(_d))
            if _dp.exists():
                _extra_shps += list(_dp.rglob("*.shp"))
                logger.info(f"  extra search dir: {_dp} ({len(_extra_shps)} shapefiles so far)")
            else:
                logger.warning(f"  ⚠ SK_COUNTY_EXTRA_SEARCH_DIRS entry does not exist: {_dp}")
        _cand = []
        if Config.SK_COUNTY_SHP_PATH:
            _explicit = Path(str(Config.SK_COUNTY_SHP_PATH))
            if _explicit.exists():
                logger.info(f"  using SK_COUNTY_SHP_PATH (hint search skipped): {_explicit}")
                _cand.append(("SK_COUNTY_SHP_PATH", _explicit))
            else:
                logger.warning(f"  ⚠ SK_COUNTY_SHP_PATH is set but does not exist: {_explicit} "
                               f"— falling back to the hint search.")
        if not _cand:
            for _h in Config.SK_COUNTY_SHP_HINTS:
                _f = find_shp(_h)
                if _f is None:
                    _f = next((p for p in _extra_shps if _h in p.name.lower()), None)
                if _f is not None and str(_f) not in [str(x[1]) for x in _cand]:
                    _cand.append((_h, _f))
        if _cand:
            logger.info("  candidate county layers (first one is used):")
            for _h, _f in _cand:
                logger.info(f"    hint '{_h}' -> {_f}")
        SK_COUNTY_SHP = _cand[0][1] if _cand else None
        if SK_COUNTY_SHP is None:
            _all_names = sorted({p.name for p in ALL_SHPS})
            SK_ZONE_FAIL = ("not_found", "no county/sigungu shapefile matched the hints")
            logger.warning("  ⚠ no county/sigungu shapefile found "
                           f"(hints: {Config.SK_COUNTY_SHP_HINTS}).")
            logger.warning(f"    shapefiles the model CAN see ({len(_all_names)}): "
                           f"{_all_names if len(_all_names) <= 40 else _all_names[:40] + ['…']}")
            logger.warning("    -> put sig.shp (with .dbf/.shx/.prj) in one of the searched "
                           "roots, or add the filename you do have to "
                           "Config.SK_COUNTY_SHP_HINTS.")
        else:
            logger.info(f"  ✓ county layer: {SK_COUNTY_SHP}")
            try:
                _zres = build_county_zones(
                    SK_COUNTY_SHP, WKT, GT5, NR, NC,
                    oversample=int(Config.SK_COUNTY_OVERSAMPLE),
                    encoding=Config.SK_COUNTY_SHP_ENCODING,
                    drop_north_of_row=(DMZ_ROW if Config.SK_COUNTY_DROP_NORTH_OF_MDL else None),
                    logger=logger)
                SK_COUNTY_ZONES = _zres["zones"]
                SK_ZONE_RASTER = _zres["zone_raster"]
                SK_ZONE_META = _zres
                _cell_km = PW / 1000.0
                _small = [z for z in SK_COUNTY_ZONES
                          if np.isfinite(z.get("poly_area_km2", np.nan))
                          and z["poly_area_km2"] < _cell_km ** 2]
                logger.info(f"  {len(SK_COUNTY_ZONES)} counties on a {_cell_km:.1f} km grid "
                            f"({_cell_km**2:.0f} km²/cell); {len(_small)} counties are smaller "
                            f"than one model cell.")
            except CountyZoneError as e:
                SK_ZONE_FAIL = (e.reason, str(e))
                logger.warning(f"  ⚠ county zone build FAILED [{e.reason}]: {e}")
                SK_COUNTY_ZONES = []
            except Exception as e:
                SK_ZONE_FAIL = (type(e).__name__, str(e))
                logger.warning(f"  ⚠ county zone build FAILED with an unexpected error "
                               f"({type(e).__name__}: {e}) — full traceback below.",
                               exc_info=True)
                SK_COUNTY_ZONES = []
    # ── 3b. Study area + outbreak-point seeding ──────────────────────────────
    logger.info("\n[STUDY AREA + OUTBREAK POINTS]")
    def _first_geom_name(shp_path):
        try:
            ds_ = _ogr.Open(str(shp_path))
            if ds_ is None:
                return None
            lyr_ = ds_.GetLayer()
            nm = None
            feat_ = lyr_.GetNextFeature()
            if feat_ is not None:
                g_ = feat_.GetGeometryRef()
                if g_ is not None:
                    nm = g_.GetGeometryName().upper()
            ds_ = None
            return nm
        except Exception as e:
            logger.warning(f"  ⚠ could not read geometry type of {shp_path}: {e}")
            return None
    def _load_points_map(shp_path):
        """Return outbreak points as (x_map, y_map) in the raster CRS."""
        pts = []
        try:
            ds_ = _ogr.Open(str(shp_path))
            if ds_ is None:
                return pts
            lyr_ = ds_.GetLayer()
            src_srs = lyr_.GetSpatialRef()
            dst_srs = _osr.SpatialReference(); dst_srs.ImportFromWkt(WKT)
            tfm = (_osr.CoordinateTransformation(src_srs, dst_srs)
                   if src_srs and not src_srs.IsSame(dst_srs) else None)
            for feat_ in lyr_:
                g = feat_.GetGeometryRef()
                if g is None:
                    continue
                g = g.Clone()
                if tfm:
                    g.Transform(tfm)
                nm = g.GetGeometryName().upper()
                if 'POINT' not in nm:
                    continue
                if 'MULTI' in nm:
                    for i in range(g.GetGeometryCount()):
                        pp = g.GetGeometryRef(i)
                        pts.append((pp.GetX(), pp.GetY()))
                else:
                    pts.append((g.GetX(), g.GetY()))
            ds_ = None
        except Exception as e:
            logger.warning(f"  ⚠ point load failed for {shp_path}: {e}")
        return pts
    def _xy_to_rc(x, y):
        return point_cell(x, y, (GT5[0], GT5[1], 0.0, GT5[3], 0.0, -PH))
    STUDY_AREA_MASK = None
    study_segments = []
    OUTBREAK_ENTRIES = []          # (row, col, label, 0, weight) — grid indices
    def _points_to_entries(shp_path, label):
        pts = _load_points_map(shp_path)
        entries, dropped = [], 0
        for (x, y) in pts:
            r_, c_ = _xy_to_rc(x, y)
            if 0 <= r_ < NR and 0 <= c_ < NC:
                entries.append((r_, c_, label, 0, Config.OUTBREAK_SEED_WEIGHT))
            else:
                dropped += 1
        logger.info(f"    {label}: {len(pts)} points, {len(entries)} inside grid, "
                    f"{dropped} dropped (outside grid — bad coords/CRS if many).")
        if entries:                          # v9.10 item 7: remember which points seed the run
            V910_SEED_INFO.clear()
            V910_SEED_INFO.update(source="shapefile", shapefile=str(shp_path), before=None,
                                  seed_rows=len(entries), seed_cells=len({(e[0], e[1]) for e in entries}))
            V910_SEED_INFO["_items"] = [dict(x=x, y=y, row=_xy_to_rc(x, y)[0], col=_xy_to_rc(x, y)[1])
                                        for (x, y) in pts
                                        if 0 <= _xy_to_rc(x, y)[0] < NR and 0 <= _xy_to_rc(x, y)[1] < NC]
        return entries
    STUDY_SHP = None
    for _h in Config.STUDY_AREA_SHP_HINTS:
        STUDY_SHP = find_shp(_h)
        if STUDY_SHP is not None:
            break
    study_is_points = False
    if STUDY_SHP is not None:
        gtype = _first_geom_name(STUDY_SHP) or ""
        logger.info(f"  ✓ Study-area shapefile: {STUDY_SHP.name} (geometry: {gtype or 'unknown'})")
        if 'POLYGON' in gtype:
            try:
                STUDY_AREA_MASK = rasterize_boundary(STUDY_SHP, WKT, GT5, NR, NC, dilate=False)
                logger.info(f"    study-area polygon rasterized: {int(STUDY_AREA_MASK.sum()):,} cells")
            except Exception as e:
                logger.warning(f"    ⚠ study-area rasterization failed: {e}")
            for (cs, rs) in load_polys(STUDY_SHP.stem.lower(), "Study area outline"):
                study_segments.append(np.column_stack([cs, rs]))
        elif 'POINT' in gtype:
            study_is_points = True
            logger.info("    study-area layer is POINTS — it will only be used to seed if no "
                        "dedicated occurrence/outbreak file is found.")
    else:
        logger.info("  Study-area shapefile not found (optional) — full peninsula kept.")
    # ══ v9.10 ITEM 7: temporal hold-out — seed only from the outbreaks of a dated table
    # reported before Config.SEED_BEFORE (--preset holdout2015). Never falls back to the
    # full outbreak file: v910_seed_entries() stops the run if nothing is left.
    if getattr(Config, "SEED_TABLE", None):
        _seed_tab = v910_find_seed_table(Config.SEED_TABLE,
                                         [Config.CWD, V910_SCRIPT_DIR, Config.DATA_ROOT,
                                          Path(str(Config.DATA_ROOT)) / "PPR_occ"])
        _seed_rows, _seed_info = v910_seed_rows(_seed_tab, Config.SEED_TABLE_DATE_COL,
                                                Config.SEED_BEFORE, Config.SEED_UNDATED)
        OUTBREAK_ENTRIES = v910_seed_entries(_seed_rows, WKT, _xy_to_rc, NR, NC,
                                             Config.OUTBREAK_SEED_WEIGHT, _seed_info, logger)
    if Config.SEED_FROM_OUTBREAK_SHP and not OUTBREAK_ENTRIES:
        _tried = set()
        for _h in Config.OUTBREAK_SHP_HINTS:
            op = find_shp(_h)
            if op is None or str(op) in _tried:
                continue
            if STUDY_SHP is not None and op == STUDY_SHP:
                continue  # the study-area file is handled as last resort below
            _tried.add(str(op))
            ent = _points_to_entries(op, op.stem)
            if ent:
                OUTBREAK_ENTRIES = ent
                logger.info(f"  ✓ Seeding from occurrence/outbreak file: {op.name} "
                            f"({len(ent)} valid points).")
                break
        if not OUTBREAK_ENTRIES and study_is_points and STUDY_SHP is not None:
            ent = _points_to_entries(STUDY_SHP, "studyarea_pt")
            if ent:
                OUTBREAK_ENTRIES = ent
                logger.info(f"  ✓ Seeding from study-area POINTS (no dedicated file found): "
                            f"{len(ent)} valid points.")
    if OUTBREAK_ENTRIES:
        logger.info(f"  Outbreak seed points ready: {len(OUTBREAK_ENTRIES)} (spread will start here)")
        if getattr(Config, "WRITE_SEED_LIST", True):                    # v9.10 item 7
            try:
                v910_copy_to_dl(v910_write_seed_list(Config.OUT_DIR, WKT, logger), Config.DL_DIR, logger)
            except Exception as _e:
                logger.warning(f"  ⚠ v9.10 seed list not written ({type(_e).__name__}: {_e})")
    elif Config.SEED_FROM_OUTBREAK_SHP:
        logger.info("  No valid in-grid outbreak points found — falling back to corridor-gate "
                    "(or legacy DMZ/NK) seeding.")
    # ══ v9.9 FIX C: DPRK/ROK split from the national boundary ══════════════
    # v9.8 split the domain at ONE grid row (the median MDL row): DPRK land south of that row
    # (Kaesong–Haeju–Ongjin) counted as South Korea and ROK land north of it (Goseong) as North
    # Korea. Inside the Korean region the split now follows the ROK sigungu polygons; outside it
    # (China) the legacy north/south halves are kept so continental plots are unchanged.
    KOREA_REGION = np.isfinite(arr5)              # same basis as KOREA_MASK below
    # ══ v9.10 ITEM 1: the split comes from the ROK sigungu SUB-CELLS that the
    # county CSVs already use (build_county_zones reprojects each polygon itself).
    # In the 6 Nov run rasterize_boundary(sig.shp) returned nothing and the split
    # silently fell back to the median MDL row, so DPRK land west of the fence
    # (Kaepung–Kaesong–Haeju–Ongjin peninsula) sat inside Ω_SK. Fallbacks: OGR
    # rasterize, then the row split (logged loudly).
    ROK_RASTER = None
    ROK_FRACTION = None
    LAND_FRACTION = None
    ROK_SPLIT_INFO = dict(method="row (median MDL row)")
    if Config.SPLIT_KOREA_BY_BOUNDARY and str(Config.ROK_SPLIT_METHOD).lower() != "row":
        _zs_split, _os_split = SK_COUNTY_ZONES, (SK_ZONE_META or {}).get("oversample")
        _shp_split = SK_COUNTY_SHP
        if _shp_split is None:                      # the county export may be switched off
            if Config.SK_COUNTY_SHP_PATH and Path(str(Config.SK_COUNTY_SHP_PATH)).exists():
                _shp_split = Path(str(Config.SK_COUNTY_SHP_PATH))
            else:
                for _h in Config.SK_COUNTY_SHP_HINTS:
                    _shp_split = find_shp(_h)
                    if _shp_split is not None:
                        break
        if (not _zs_split or Config.SK_COUNTY_DROP_NORTH_OF_MDL) and _shp_split is not None:
            try:
                _zres_split = build_county_zones(
                    _shp_split, WKT, GT5, NR, NC,
                    oversample=int(Config.SK_COUNTY_OVERSAMPLE),
                    encoding=Config.SK_COUNTY_SHP_ENCODING, drop_north_of_row=None, logger=logger)
                _zs_split, _os_split = _zres_split["zones"], _zres_split["oversample"]
            except Exception as _e:
                logger.warning(f"  ⚠ v9.10: sigungu sub-cells unavailable for the split "
                               f"({type(_e).__name__}: {_e})")
        if _zs_split:
            ROK_FRACTION = rok_fraction_from_zones(_zs_split, NR, NC, _os_split)
            if Config.ROK_SPLIT_USE_LAND_FRACTION and GORANI_TIF is not None:
                try:
                    LAND_FRACTION = land_fraction_on_grid(GORANI_TIF, GT5, NR, NC, WKT)
                except Exception as _e:
                    logger.warning(f"  ⚠ v9.10: land fraction not computed ({type(_e).__name__}: {_e}); "
                                   f"border cells are judged on the ROK share of the cell area.")
            ROK_RASTER, ROK_SPLIT_INFO = split_rok_dprk(
                KOREA_REGION, ROK_FRACTION, LAND_FRACTION,
                core_eps=float(Config.ROK_SPLIT_CORE_EPS),
                border_min_land_share=float(Config.ROK_BORDER_MIN_LAND_SHARE))
            ROK_SPLIT_INFO["method"] = (f"ROK sigungu sub-cells ({len(_zs_split)} polygons, "
                                        f"{_os_split}x oversample)"
                                        + (" + land share at the border" if LAND_FRACTION is not None else ""))
        elif _shp_split is not None:
            try:
                _rr_split = rasterize_boundary(_shp_split, WKT, GT5, NR, NC, dilate=False)
                if _rr_split.any():
                    ROK_RASTER = _rr_split & KOREA_REGION
                    ROK_FRACTION = ROK_RASTER.astype(np.float64)
                    ROK_SPLIT_INFO = dict(method=f"rasterize_boundary({Path(str(_shp_split)).name})")
            except Exception as _e:
                logger.warning(f"  ⚠ v9.10: rasterize_boundary failed as well ({_e})")
    if ROK_RASTER is None or not (ROK_RASTER & KOREA_REGION).any():
        ROK_RASTER = None
        ROK_FRACTION = None
        logger.warning("  ⚠⚠ v9.10: NO national split — NK/SK fall back to the median MDL row. "
                       "Ω_SK then contains DPRK land west of the fence: do NOT report ROK values "
                       "from this run. Check that sig.shp (with .dbf/.shx/.prj) is found.")
    else:
        _ri0 = np.tile(np.arange(NR).reshape(-1, 1), (1, NC))
        _dprk_south = int((KOREA_REGION & ~ROK_RASTER & (_ri0 >= DMZ_ROW)).sum())
        _rok_north = int((KOREA_REGION & ROK_RASTER & (_ri0 < DMZ_ROW)).sum())
        ROK_SPLIT_INFO.update(dprk_cells_south_of_row=_dprk_south, rok_cells_north_of_row=_rok_north)
        logger.info(f"  ✓ v9.10 DPRK/ROK split — {ROK_SPLIT_INFO['method']}: "
                    f"{int((KOREA_REGION & ROK_RASTER).sum()):,} ROK / "
                    f"{int((KOREA_REGION & ~ROK_RASTER).sum()):,} DPRK cells; border cells given to "
                    f"the DPRK: {ROK_SPLIT_INFO.get('n_border_to_dprk', '—')}. Against the old row "
                    f"split, {_dprk_south:,} DPRK cells south of row {DMZ_ROW} and {_rok_north:,} ROK "
                    f"cells north of it are now labelled correctly.")
    def split_korea(kp):
        """(NK, SK) masks: national boundary inside the Korean region, legacy row halves outside."""
        _ri = np.tile(np.arange(NR).reshape(-1, 1), (1, NC))
        nk_rows, sk_rows = kp & (_ri < DMZ_ROW), kp & (_ri >= DMZ_ROW)
        if ROK_RASTER is None:
            return nk_rows, sk_rows
        kr = KOREA_REGION
        return ((nk_rows & ~kr) | (kp & kr & ~ROK_RASTER)), ((sk_rows & ~kr) | (kp & kr & ROK_RASTER))
    NK, SK = split_korea(KP)
    if Config.STUDY_AREA_DEFINES_DOMAIN and STUDY_AREA_MASK is not None:
        _newKP = STUDY_AREA_MASK.copy()
        if int(_newKP.sum()) >= Config.STUDY_AREA_MIN_CELLS:
            KP = _newKP
            NK, SK = split_korea(KP)                     # v9.9 FIX C
            VALID = VALID | KP
            GOR_nn = np.where(KP, GOR_nn, 0.0).astype(np.float32)
            logger.info(f"  ✓ Study area DEFINES the domain: KP={int(KP.sum()):,} "
                        f"(NK/north={int(NK.sum()):,}, SK/south={int(SK.sum()):,}); "
                        f"epidemic runs across the full study area incl. other-country land.")
        else:
            logger.warning(f"  ⚠ Study-area domain has only {int(_newKP.sum())} cells "
                           f"(< {Config.STUDY_AREA_MIN_CELLS}); keeping gorani domain. Check the "
                           f"shapefile CRS/extent vs the raster.")
    elif Config.CLIP_TO_STUDY_AREA and STUDY_AREA_MASK is not None:
        _newKP = KP & STUDY_AREA_MASK
        if int(_newKP.sum()) >= Config.STUDY_AREA_MIN_CELLS:
            KP = _newKP
            NK = NK & STUDY_AREA_MASK
            SK = SK & STUDY_AREA_MASK
            VALID = VALID & STUDY_AREA_MASK
            GOR_nn = np.where(KP, GOR_nn, 0.0).astype(np.float32)
            logger.info(f"  ✓ Clipped model domain to study area: "
                        f"KP={int(KP.sum()):,} (NK={int(NK.sum()):,}, SK={int(SK.sum()):,})")
        else:
            logger.warning(f"  ⚠ Study-area clip would leave only {int(_newKP.sum())} cells "
                           f"(< {Config.STUDY_AREA_MIN_CELLS}); keeping full extent.")
    # ── 4. Environmental rasters ─────────────────────────────────────────────
    logger.info("\n[LOADING ENVIRONMENTAL RASTERS]")
    def warp(p):
        if not Path(str(p)).exists():
            logger.warning(f"  ⚠ {p} not found, using NaN fill")
            return np.full(SH, np.nan, np.float32)
        try:
            wo = gdal.WarpOptions(
                format="MEM", dstSRS=WKT,
                outputBounds=(GT5[0], GT5[3] + NR * GT5[5], GT5[0] + NC * GT5[1], GT5[3]),
                width=NC, height=NR, resampleAlg=gdal.GRA_Average)
            ds_w = gdal.Warp("", str(p), options=wo)
            if ds_w is None:
                logger.warning(f"  ⚠ Warp returned None for {p}")
                return np.full(SH, np.nan, np.float32)
            a = ds_w.GetRasterBand(1).ReadAsArray().astype(np.float32)
            nd2 = ds_w.GetRasterBand(1).GetNoDataValue()
            if nd2:
                a[a == nd2] = np.nan
            a[~np.isfinite(a)] = np.nan
            ds_w = None
            return a
        except Exception as e:
            logger.warning(f"  ⚠ Warp failed for {p}: {e}")
            return np.full(SH, np.nan, np.float32)
    def norm01(a):
        v = a[KP & np.isfinite(a)]
        if not len(v):
            return np.zeros(SH, np.float32)
        mn, mx = v.min(), v.max()
        return np.where(KP, np.clip((a - mn) / (mx - mn + 1e-9), 0, 1), np.nan).astype(np.float32)
    ELEV = warp(Config.ENV_DIR / "wc2.1_30s_elev_studyarea_1km_5179.tif") if Config.ENV_DIR.exists() else np.full(SH, 500, np.float32)
    _lv_path = v910_livestock_path()
    LIVESTOCK = warp(_lv_path) if Config.ENV_DIR.exists() else np.full(SH, 0.1, np.float32)
    Config.LIVESTOCK_USED = _lv_path.name
    if _lv_path.name == Config.LIVESTOCK_GLW_NAME:
        logger.warning(f"  ⚠ goat layer {_lv_path.name}: GLW4 clipped at the MDL — the ROK has NO "
                       f"small-ruminant data in this file (V2 and the HOST criterion are not meaningful "
                       f"there). Put {Config.LIVESTOCK_ROKFILLED_NAME} in {Config.ENV_DIR}.")
    else:
        logger.info(f"  Goat layer: {_lv_path.name} (GLW4 north of the MDL, national goat register south of it)")
    PRECIP = warp(Config.ENV_DIR / "wc2.1_30s_bio_12_studyarea_1km_5179.tif") if Config.ENV_DIR.exists() else np.full(SH, 1200, np.float32)
    logger.info(f"  Elevation: finite={np.isfinite(ELEV).sum()}/{ELEV.size}")
    logger.info(f"  Livestock: finite={np.isfinite(LIVESTOCK).sum()}/{LIVESTOCK.size}")
    logger.info(f"  Precip:    finite={np.isfinite(PRECIP).sum()}/{PRECIP.size}")
    # ── 4c. PPR MaxEnt suitability raster ────────────────────────────────────
    def _find_raster(hints):
        for root in [Config.CWD, Config.DATA_ROOT, Config.ENV_DIR]:
            if not root.exists():
                continue
            for t in root.rglob("*.tif"):
                if t.name.startswith("PPR_v9") or "PPR_v9.10_rasters" in t.parts:   # v9.10: own outputs
                    continue
                nl = t.name.lower()
                if any(h in nl for h in hints):
                    return t
        return None
    MAXENT = None
    if base_is_maxent and MAXENT_BASE is not None:
        MAXENT = MAXENT_BASE            # already on-grid (it IS the base raster)
        logger.info(f"  ✓ PPR MaxEnt map: base raster, finite={int(np.isfinite(MAXENT).sum()):,} cells")
    else:
        MAXENT_TIF2 = (Path(str(Config.MAXENT_TIF_PATH)) if getattr(Config, "MAXENT_TIF_PATH", None)   # v9.10
                       else _find_raster([h.lower() for h in Config.MAXENT_TIF_HINTS]))
        if MAXENT_TIF2 is not None:
            logger.info(f"  ✓ PPR MaxEnt map: {MAXENT_TIF2.name}")
            _mx = warp(MAXENT_TIF2)
            _mv = _mx[np.isfinite(_mx)]
            if _mv.size and float(np.nanmax(_mx)) > float(np.nanmin(_mx)):
                _mn, _mx2 = float(np.nanmin(_mx)), float(np.nanmax(_mx))
                MAXENT = np.clip((_mx - _mn) / (_mx2 - _mn + 1e-9), 0, 1).astype(np.float32)
            else:
                MAXENT = _mx
            MAXENT_TIF = MAXENT_TIF2
            logger.info(f"    MaxEnt suitability: finite={int(np.isfinite(MAXENT).sum()):,} cells, "
                        f"mean={float(np.nanmean(MAXENT)) if _mv.size else 0:.3f}")
        else:
            logger.info("  PPR MaxEnt map: not found (optional).")
    # ── 4b. Regional land basemap (neighbouring countries) ───────────────────
    LAND_BASE = np.zeros(SH, bool)
    basemap_segments = []
    BASEMAP_SHP = None
    if Config.SHOW_REGIONAL_BASEMAP:
        for _h in Config.BASEMAP_SHP_HINTS:
            BASEMAP_SHP = find_shp(_h)
            if BASEMAP_SHP is not None:
                break
        if BASEMAP_SHP is not None:
            try:
                LAND_BASE = rasterize_boundary(BASEMAP_SHP, WKT, GT5, NR, NC, dilate=False)
                for (cs, rs) in load_polys(BASEMAP_SHP.stem.lower(), "Regional basemap"):
                    basemap_segments.append(np.column_stack([cs, rs]))
                logger.info(f"  ✓ Regional basemap: {BASEMAP_SHP.name} "
                            f"({int(LAND_BASE.sum()):,} land cells, {len(basemap_segments)} outlines)")
            except Exception as e:
                logger.warning(f"  ⚠ basemap rasterization failed: {e}")
        if not LAND_BASE.any():
            LAND_BASE = np.isfinite(ELEV)
            _extra = int((LAND_BASE & ~VALID).sum())
            logger.info(f"  Regional basemap: no land/countries shapefile found — using elevation "
                        f"land cells ({int(LAND_BASE.sum()):,}; {_extra:,} beyond the gorani area).")
    if Config.BASEMAP_DEFINES_DOMAIN and LAND_BASE.any():
        KP = LAND_BASE.copy()
        NK, SK = split_korea(KP)                         # v9.9 FIX C
        VALID = VALID | KP
        GOR_nn = np.where(KP, GOR_nn, 0.0).astype(np.float32)
        logger.info(f"  Basemap DEFINES the domain: epidemic runs across ALL country land "
                    f"({int(KP.sum()):,} cells; NK/north={int(NK.sum()):,}, SK/south={int(SK.sum()):,}).")
    elev_n = norm01(ELEV)
    livestock_n = norm01(LIVESTOCK)
    precip_n = norm01(PRECIP)
    elev_nn = np.nan_to_num(ELEV, nan=500)
    RIV = np.where(
        KP,
        np.exp(-0.0025 * np.clip(elev_nn, 0, 1500)) * (0.4 + 0.6 * np.nan_to_num(precip_n, nan=0.5)),
        0.
    ).astype(np.float32)
    RIV = gaussian_filter(RIV, sigma=3).astype(np.float32)
    RIV_nn = np.nan_to_num(RIV, nan=0)
    if RIV_nn[KP].max() > 0:
        RIV_nn = RIV_nn / RIV_nn[KP].max()
    if Config.STUDY_AREA_DEFINES_DOMAIN and STUDY_AREA_MASK is not None:
        _land = np.isfinite(ELEV) | (GOR_nn > 0.0) | (LAND_BASE if LAND_BASE.any() else False)
        _kp_land = KP & _land
        if int(_kp_land.sum()) >= max(Config.STUDY_AREA_MIN_CELLS, int(0.15 * int(KP.sum()))):
            _dropped = int(KP.sum()) - int(_kp_land.sum())
            KP = _kp_land
            NK, SK = split_korea(KP)                     # v9.9 FIX C
            GOR_nn = np.where(KP, GOR_nn, 0.0).astype(np.float32)
            logger.info(f"  ✓ Domain trimmed to land: KP={int(KP.sum()):,} "
                        f"(dropped {_dropped:,} open-water cells; NK/north={int(NK.sum()):,}, "
                        f"SK/south={int(SK.sum()):,}).")
        else:
            logger.warning("  ⚠ Elevation land-coverage looks sparse; keeping the full "
                           "study-area domain (spread may include water cells).")
    # ══ v9.10 ITEM 1: ROK exclaves on the grid (e.g. the northern Ongjin islands) ══
    # At 8 km some ROK island cells share a face with DPRK coastal cells, so the
    # model lets infection walk across the sea from the DPRK and bypass the fence.
    # 'report' keeps them and flags them; 'mask' removes them from the peninsular leg.
    ROK_EXCLAVES = np.zeros(SH, bool)
    ROK_EXCLAVE_INFO = dict(policy=str(Config.ROK_EXCLAVE_POLICY).lower(), cells=0,
                            components=[], sigungu=[])
    if ROK_RASTER is not None:
        try:
            _kr_kp = KP & KOREA_REGION
            ROK_EXCLAVES, _exs = rok_exclave_cells(_kr_kp & ROK_RASTER, _kr_kp & ~ROK_RASTER)
            _ex_names = {}
            if ROK_EXCLAVES.any() and SK_COUNTY_ZONES:
                _zb_ex = {int(z["zid"]): z for z in SK_COUNTY_ZONES}
                for _z in np.unique(SK_ZONE_RASTER[ROK_EXCLAVES]):
                    if int(_z) in _zb_ex:
                        _zz = _zb_ex[int(_z)]
                        _ex_names[str(_zz.get("name_en") or _zz.get("name_kr"))] = int(
                            ((SK_ZONE_RASTER == _z) & ROK_EXCLAVES).sum())
            ROK_EXCLAVE_INFO.update(cells=int(ROK_EXCLAVES.sum()), components=_exs,
                                    sigungu=sorted(_ex_names))
            if ROK_EXCLAVES.any():
                logger.info(f"  v9.10 ROK exclaves touching DPRK land on this grid: "
                            f"{int(ROK_EXCLAVES.sum())} cells in {len(_exs)} island group(s) "
                            f"({', '.join(f'{k}: {v}' for k, v in _ex_names.items()) or 'n/a'}); "
                            f"policy '{ROK_EXCLAVE_INFO['policy']}'"
                            + (" — removed from the peninsular leg."
                               if ROK_EXCLAVE_INFO["policy"] == "mask"
                               else " — kept and flagged (ROK_EXCLAVE_POLICY='mask' removes them)."))
        except Exception as _e:
            logger.warning(f"  ⚠ v9.10 exclave check failed ({type(_e).__name__}: {_e})")
    ROK_EXCLAVES_DROP = (ROK_EXCLAVES.copy() if ROK_EXCLAVE_INFO["policy"] == "mask"
                         else np.zeros(SH, bool))
    KOREA_REGION_MODEL = KOREA_REGION & ~ROK_EXCLAVES_DROP
    # ── 5. Resistance / conductance / cost graph ────────────────────────────
    logger.info("\n[COMPUTING RESISTANCE & COST GRAPH]")
    _elev_f = np.nan_to_num(ELEV, nan=500).astype(np.float64)
    _gy, _gx = np.gradient(_elev_f, PH, PW)
    slope_deg = np.degrees(np.arctan(np.sqrt(_gx ** 2 + _gy ** 2)))
    SLOPE_N = np.clip(slope_deg / 45.0, 0, 1).astype(np.float32)
    _water_src = (RIV_nn > 0.40)
    DIST_WATER_M = distance_transform_edt(~_water_src) * PW
    WATER_ATTRACT = np.exp(-DIST_WATER_M / 5000.0).astype(np.float32)
    livestock_nn = np.nan_to_num(livestock_n, nan=0.0)
    URBAN_PROXY = np.clip((1.0 - GOR_nn) * livestock_nn, 0, 1).astype(np.float32)
    _En_nn = np.nan_to_num(elev_n, nan=0.5).astype(np.float32)
    _rw = dict(Config.RESISTANCE_WEIGHTS)
    _rw_sum = float(sum(_rw.values())) or 1.0
    RESISTANCE = np.where(
        KP,
        (_rw["slope"] * SLOPE_N + _rw["urban"] * URBAN_PROXY
         + _rw["elev"] * _En_nn + _rw["water"] * (1.0 - WATER_ATTRACT)) / _rw_sum,
        1.0
    ).astype(np.float32)
    logger.info("  Resistance weights (normalised): "
                + ", ".join(f"{k}={v/_rw_sum:.2f}" for k, v in _rw.items())
                + "   [SLOPE_WEIGHT elasticity 0.293 -> trimmed from 0.35 in v9.1]")
    ROAD_MASK = np.zeros(SH, bool)
    if Config.HAN_APPLY_ROAD_BARRIER and ROAD_SHP is not None:
        try:
            ROAD_MASK = rasterize_boundary(ROAD_SHP, WKT, GT5, NR, NC)
            RESISTANCE = np.where(
                KP & ROAD_MASK,
                np.clip(RESISTANCE + Config.HAN_ROAD_RESISTANCE_ADD, 0, 1.5),
                RESISTANCE
            ).astype(np.float32)
            logger.info(f"  ✓ Road barrier applied (Han delta={Config.HAN_ROAD_REDUCTION}): "
                        f"{int((ROAD_MASK & KP).sum()):,} road cells, "
                        f"+{Config.HAN_ROAD_RESISTANCE_ADD} resistance each")
        except Exception as e:
            logger.warning(f"  ⚠ Road-barrier rasterization failed: {e} -- skipping road barrier")
    if Config.HAN_APPLY_RIVER_BARRIER:
        _river_barrier = (RIV_nn > 0.40)
        RESISTANCE = np.where(
            KP & _river_barrier,
            np.clip(RESISTANCE + Config.HAN_RIVER_RESISTANCE_ADD, 0, 1.5),
            RESISTANCE
        ).astype(np.float32)
        logger.info(f"  ✓ River barrier applied (Han nu={Config.HAN_RIVER_REDUCTION}) "
                    f"[NOTE: enabled despite gorani water-association — verify this is intended]")
    CONDUCTANCE = np.where(KP, np.exp(-Config.CONDUCTANCE_EXPONENT * RESISTANCE),
                           1e-9).astype(np.float32)
    logger.info(f"  Conductance exponent k={Config.CONDUCTANCE_EXPONENT} (v9.1 used 3.0) — "
                f"HIGHEST-elasticity parameter (mean 0.549, max 2.735): lower k = more "
                f"permeable landscape = higher mobility = faster front; "
                f"mean CONDUCTANCE={float(np.nanmean(CONDUCTANCE[KP])):.3f}")
    MCP_COST = np.where(KP, 1.0 / (CONDUCTANCE + 1e-7), 1e7).astype(np.float64)
    _REF_COST = float(np.median(MCP_COST[KP & (MCP_COST < 1e6)]))
    logger.info(f"  Reference cost (median, used to calibrate distances to km): {_REF_COST:.2f}")
    def _build_cost_graph(cost_arr):
        nr, nc = cost_arr.shape
        idx = np.arange(nr * nc, dtype=np.int32).reshape(nr, nc)
        dirs = [(-1, 0, 1.0), (1, 0, 1.0), (0, -1, 1.0), (0, 1, 1.0),
                (-1, -1, 1.414), (-1, 1, 1.414), (1, -1, 1.414), (1, 1, 1.414)]
        rs, cs, ws = [], [], []
        for dr, dc, dist in dirs:
            r1 = np.arange(max(0, -dr), min(nr, nr - dr))
            c1 = np.arange(max(0, -dc), min(nc, nc - dc))
            r2, c2 = r1 + dr, c1 + dc
            rg1, cg1 = np.meshgrid(r1, c1, indexing='ij')
            rg2, cg2 = np.meshgrid(r2, c2, indexing='ij')
            i_ = idx[rg1, cg1].ravel()
            j_ = idx[rg2, cg2].ravel()
            ci = cost_arr[rg1, cg1].ravel()
            cj = cost_arr[rg2, cg2].ravel()
            rs.append(i_); cs.append(j_)
            ws.append(dist * (ci + cj) / 2.0)
        rs = np.concatenate(rs); cs = np.concatenate(cs); ws = np.concatenate(ws)
        return csr_matrix((ws, (rs, cs)), shape=(nr * nc, nr * nc))
    _mcp_capped = np.minimum(MCP_COST, 1e4).astype(np.float32)
    COST_GRAPH = _build_cost_graph(_mcp_capped)
    logger.info(f"  ✓ Cost graph built: {COST_GRAPH.nnz} directed edges over {NR * NC} nodes")
    def mcp_dist_km(src_r, src_c):
        src_idx = int(src_r) * NC + int(src_c)
        dist_flat = sp_dijkstra(COST_GRAPH, indices=src_idx, directed=False)
        dist_arr = dist_flat.reshape(NR, NC).astype(np.float32)
        return (dist_arr / (_REF_COST + 1e-9) * (PW / 1000.0)).astype(np.float32)
    def mcp_dist_km_multisource(mask):
        """Cost-distance (km) from the nearest TRUE cell in `mask`."""
        src_idx = np.where(mask.ravel())[0]
        if src_idx.size == 0:
            return np.full(SH, np.inf, np.float32)
        try:
            dist_flat = sp_dijkstra(COST_GRAPH, indices=src_idx, directed=False, min_only=True)
        except TypeError:
            dist_mat = sp_dijkstra(COST_GRAPH, indices=src_idx, directed=False)
            dist_flat = np.min(dist_mat, axis=0)
        dist_arr = dist_flat.reshape(NR, NC).astype(np.float32)
        return (dist_arr / (_REF_COST + 1e-9) * (PW / 1000.0)).astype(np.float32)
    def rc_from_km(x_km, y_km):
        c_ = max(0, min(NC - 1, int((x_km * 1000 - GT5[0]) / GT5[1])))
        r_ = max(0, min(NR - 1, int((GT5[3] - y_km * 1000) / PH)))
        return r_, c_
    # ── 5b. Gorani density composite (gorani TIF + covariate TIFs) ───────────
    _w = Config.GORANI_DENSITY_WEIGHTS
    _elev_nn0 = np.nan_to_num(elev_n, nan=0.5)
    _gor_density_raw = (
        _w["gorani"] * GOR_nn
        + _w["water"] * np.nan_to_num(WATER_ATTRACT, nan=0.0)
        + _w["low_elev"] * (1.0 - _elev_nn0)
        + _w["low_slope"] * (1.0 - np.nan_to_num(SLOPE_N, nan=0.5))
    ).astype(np.float32)
    _gd = _gor_density_raw[KP]
    if _gd.size and _gd.max() > _gd.min():
        GORANI_DENSITY = np.where(
            KP, np.clip((_gor_density_raw - _gd.min()) / (_gd.max() - _gd.min() + 1e-9), 0, 1), np.nan
        ).astype(np.float32)
    else:
        GORANI_DENSITY = np.where(KP, 0.0, np.nan).astype(np.float32)
    if Config.USE_COMPOSITE_HABITAT_FOR_SPREAD:
        GOR_EFF = np.maximum(GOR_nn, Config.HABITAT_BLEND * np.nan_to_num(GORANI_DENSITY, nan=0.0)).astype(np.float32)
    else:
        GOR_EFF = GOR_nn
    if Config.USE_MAXENT_IN_HABITAT and MAXENT is not None:
        _mx = np.nan_to_num(MAXENT, nan=0.0).astype(np.float32)
        _have = np.isfinite(MAXENT)
        if base_is_maxent:
            GOR_EFF = np.where(_have & KP, np.maximum(GOR_EFF, _mx), GOR_EFF).astype(np.float32)
            logger.info(f"  ✓ MaxEnt drives habitat (max with gorani) on "
                        f"{int((_have & KP).sum()):,} cells.")
        else:
            _b = float(Config.MAXENT_HABITAT_BLEND)
            GOR_EFF = np.where(_have & KP,
                               np.clip((1.0 - _b) * GOR_EFF + _b * _mx, 0.0, 1.0),
                               GOR_EFF).astype(np.float32)
            logger.info(f"  ✓ MaxEnt blended into spread habitat (weight={_b:.2f} on "
                        f"{int((_have & KP).sum()):,} cells).")
    def _cov(mask, arr, thr=0.05):
        v = arr[mask & np.isfinite(arr)]
        if not v.size:
            return (0, 0.0, 0.0, 0.0)
        return (100.0 * (v > thr).mean(), float(np.nanmean(v)), float(np.nanmax(v)),
                float(np.nanpercentile(v, 90)))
    for _nm, _msk in (("NK", NK), ("SK", SK)):
        gcov = _cov(_msk, GOR_nn); dcov = _cov(_msk, GORANI_DENSITY)
        logger.info(f"  [{_nm}] gorani-TIF GOR_nn: {gcov[0]:.0f}% cells >0.05, "
                    f"mean={gcov[1]:.3f} max={gcov[2]:.3f} p90={gcov[3]:.3f}")
        logger.info(f"  [{_nm}] density composite:  {dcov[0]:.0f}% cells >0.05, "
                    f"mean={dcov[1]:.3f} max={dcov[2]:.3f} p90={dcov[3]:.3f}")
        if _nm == "NK" and gcov[2] < 0.05:
            logger.warning("  ⚠ Gorani TIF is ~empty in NK — NK risk now comes from the "
                           "covariate habitat composite (USE_COMPOSITE_HABITAT_FOR_SPREAD).")
    # v9.1: habitat quality AT each gate — a gate on unsuitable ground will not
    # actually pass much infection, so this is worth reading next to the weights.
    if CORRIDORS:
        logger.info("  [CORRIDOR HABITAT CHECK] effective habitat (GOR_EFF) at each gate:")
        for _co in CORRIDORS:
            _g = float(np.nan_to_num(GOR_EFF)[_co["row"], _co["col"]])
            _flag = "" if _g > 0.05 else "   ⚠ near-zero habitat at this gate — verify the coordinate"
            logger.info(f"    {_co['name']:<17s} w={_co['weight']:.2f}  GOR_EFF={_g:.3f}{_flag}")
    # ── 5c. Habitat / Accessible / Barrier classification (Han et al. Fig.2) ──
    logger.info("\n[HABITAT / ACCESSIBLE / BARRIER CLASSIFICATION]")
    _ge = np.nan_to_num(GOR_EFF, nan=0.0)
    _gv = _ge[KP]
    if Config.HAB_CLASS_THRESHOLD_MODE == "absolute":
        _hab_thr = float(Config.HAB_SUIT_THRESHOLD)
    else:
        _hab_thr = float(np.nanquantile(_gv, Config.HAB_SUIT_QUANTILE)) if _gv.size else 0.5
    HABITAT_MASK = KP & (_ge >= _hab_thr)
    if HABITAT_MASK.any():
        _dist_from_hab = mcp_dist_km_multisource(HABITAT_MASK)
    else:
        _dist_from_hab = np.full(SH, np.inf, np.float32)
    ACCESS_MASK = KP & (~HABITAT_MASK) & np.isfinite(_dist_from_hab) & \
        (_dist_from_hab <= Config.ACCESS_COSTDIST_MAX_KM)
    BARRIER_MASK = KP & (~HABITAT_MASK) & (~ACCESS_MASK)
    HAB_CLASS = np.full(SH, np.nan, np.float32)
    HAB_CLASS[BARRIER_MASK] = 0.0
    HAB_CLASS[ACCESS_MASK] = 1.0
    HAB_CLASS[HABITAT_MASK] = 2.0
    _ntot = max(int(KP.sum()), 1)
    logger.info(f"  Habitat threshold ({Config.HAB_CLASS_THRESHOLD_MODE}): GOR_EFF>={_hab_thr:.3f}; "
                f"reach cutoff={Config.ACCESS_COSTDIST_MAX_KM:.0f} km (Han max_dist)")
    for _lab, _msk in (("HABITAT", HABITAT_MASK), ("ACCESSIBLE", ACCESS_MASK), ("BARRIER", BARRIER_MASK)):
        _n = int(_msk.sum())
        logger.info(f"  {_lab:<11s}: {_n:7,} cells ({100.0*_n/_ntot:5.1f}%)  "
                    f"[NK={int((_msk & NK).sum()):,}, SK={int((_msk & SK).sum()):,}]")
    CLASS_MOBILITY = np.ones(SH, np.float32)
    CLASS_MOBILITY[BARRIER_MASK] = Config.CLASS_MOBILITY_BARRIER
    CLASS_MOBILITY[ACCESS_MASK] = Config.CLASS_MOBILITY_ACCESSIBLE
    CLASS_MOBILITY[HABITAT_MASK] = Config.CLASS_MOBILITY_HABITAT
    CLASS_MOBILITY = np.where(KP, CLASS_MOBILITY, 0.0).astype(np.float32)
    # ── 6. Connectivity surfaces ─────────────────────────────────────────────
    SK_CONNECTIVITY = np.where(
        SK, np.clip(GOR_EFF * 0.35 + livestock_nn * 0.35 + 0.02, 0.02, 1.0) * WATER_ATTRACT * CONDUCTANCE, 0.
    ).astype(np.float32)
    _sc = SK_CONNECTIVITY[SK]
    if _sc.size and _sc.max() > 0:
        SK_CONNECTIVITY = (SK_CONNECTIVITY / _sc.max()).astype(np.float32)
    NK_CONNECTIVITY = np.where(
        NK, np.clip(GOR_EFF * 0.60 + 0.05, 0.05, 1.0) * WATER_ATTRACT * CONDUCTANCE, 0.
    ).astype(np.float32)
    _ncv = NK_CONNECTIVITY[NK]
    if _ncv.size and _ncv.max() > 0:
        NK_CONNECTIVITY = (NK_CONNECTIVITY / _ncv.max()).astype(np.float32)
    # ── 6b. SPATIAL MCDA RISK-AREA CLASSIFICATION (v9.7) ─────────────────────
    # Built from STATIC layers only. It never sees the epidemic, which is what
    # lets the post-run validation be a test rather than a tautology.
    #
    # One classification is produced PER DOMAIN, because a risk class is only
    # meaningful relative to the area it was computed over, and because each
    # domain has its own introduction route.
    SMCDA_BY_DOMAIN: Dict[str, Dict] = {}
    SMCDA: Dict[str, object] = {}
    SMCDA_CLASS = None
    SMCDA_MOBILITY = None
    # Filled by the Korea leg once it has run, so the Korea-domain
    # classifications can be validated against the epidemic that actually
    # delivers PPR into the peninsula rather than against the continental run.
    SMCDA_KOREA_EPI: Dict[str, object] = {}
    if Config.ENABLE_SMCDA:
        logger.info("\n[SPATIAL MULTI-CRITERIA DECISION ANALYSIS — RISK-AREA THRESHOLDS]")
        # ---- derived criterion layer: local transmissibility -----------------
        # The IDENTICAL expression used inside mechanistic_closure(), with the
        # IDENTICAL Config constants, so this criterion is the model's own
        # beta_field and not a new parameterisation.
        SMCDA_BETA = np.clip(
            Config.RAD_TRANS_HAB_W * np.nan_to_num(GOR_EFF, nan=0.0)
            + Config.RAD_TRANS_LSTK_W * np.nan_to_num(livestock_nn, nan=0.0)
            + Config.RAD_TRANS_BASE, 0.0, 1.0).astype(np.float32)

        # Domains whose extent is defined BY the epidemic, so they cannot be
        # built until the Korea leg has run. They are deferred, not skipped.
        SMCDA_EPIDEMIC_DOMAINS = ("sk_reached", "korea_reached")

        def _smcda_reached_mask():
            """(footprint mask, month cutoff, description) from the Korea leg.

            A cell is in the footprint when its cumulative attack crossed
            RAD_ARRIVAL_THRESHOLD at or before the cutoff month, on the GLOBAL
            timeline — the month numbering used in the progression-figure
            titles. time_to_arrival is stored as a 0-based index within the leg,
            so the global month of arrival is (tta + 1 + month_offset)."""
            _tta = SMCDA_KOREA_EPI.get("time_to_arrival")
            if _tta is None:
                raise RuntimeError(
                    "a '_reached' domain needs the Korea-leg epidemic, which has "
                    "not run yet or did not register")
            _off = int(SMCDA_KOREA_EPI.get("month_offset", 0))
            _nm = int(SMCDA_KOREA_EPI.get("n_months", 0))
            _end = _off + _nm
            _cut = Config.SMCDA_REACHED_BY_MONTH
            _cut = _end if _cut is None else int(_cut)
            _arr_global = np.asarray(_tta, np.float64) + 1.0 + _off
            _m = np.isfinite(_arr_global) & (_arr_global <= _cut)
            _desc = (f"reached by month {_cut}"
                     + (f" (the end of the Korea leg, year {_cut // 12}"
                        f"{'' if _cut % 12 == 0 else f' month {_cut % 12}'})"
                        if Config.SMCDA_REACHED_BY_MONTH is None
                        else f" (year {_cut // 12} month {_cut % 12}; the leg ran to "
                             f"month {_end})"))
            return _m, _cut, _desc

        def _smcda_domain_mask(dk):
            """(mask, human label) for a domain key."""
            dk = str(dk).lower()
            # v9.9 FIX C: NK/SK hold legacy row halves outside Korea, so intersect with the region
            if dk == "sk":
                return (KP & SK & KOREA_REGION_MODEL).astype(bool), "all of South Korea (classes are SK-relative)"
            if dk == "nk":
                return (KP & NK & KOREA_REGION_MODEL).astype(bool), "North Korea only (classes are NK-relative)"
            if dk == "korea":
                _k = (KP & (SK | NK) & KOREA_REGION_MODEL).astype(bool)
                return _k, "the Korean peninsula, NK + SK on one scale"
            if dk == "all":
                return KP.astype(bool), "the whole model domain (China + Korea on one scale)"
            if dk in SMCDA_EPIDEMIC_DOMAINS:
                _reach, _cut, _desc = _smcda_reached_mask()
                _base = SK if dk == "sk_reached" else (SK | NK)
                _where = "South Korea" if dk == "sk_reached" else "the Korean peninsula"
                return ((KP & _base & _reach).astype(bool),
                        f"the PPR-invaded footprint within {_where} — {_desc}")
            raise ValueError(f"unknown SMCDA domain '{dk}' (expected one of: "
                             f"sk, nk, korea, all, sk_reached, korea_reached)")

        def _smcda_entry_mask(dk, src_override=None):                # v9.10: origin override
            """The ENTRY criterion's ORIGIN for this domain, as (mask, description).

            This is the single most consequential choice in the whole block. The
            criterion claims 'risk falls with travel cost from where the virus
            comes in'. That claim is only true if the origin is the route by
            which the virus actually reaches THIS domain:

              South Korea      entered through the water-deer corridor gates.
              North Korea      entered from China across the northern border.
              Korean peninsula entered from China; the gates are an internal
                               barrier on the way south, not the entry to the
                               peninsula.
              Continental      entered at the outbreak seed points themselves.

            Point it at the wrong origin and the classification will rank the
            epidemic BACKWARDS — cells nearest the assumed entry are reached
            last. The validation pass detects exactly this, and the per-criterion
            diagnostic names ENTRY as the culprit when it happens."""
            src = str(src_override or Config.SMCDA_ENTRY_SOURCE).lower()
            dk = str(dk).lower()
            if src == "off":
                return None, "disabled (SMCDA_ENTRY_SOURCE='off')"
            if src == "auto":
                src = {"sk": "gates_terminus", "sk_reached": "gates_terminus",   # v9.9 FIX D
                       "nk": "china_border", "korea": "china_border",
                       "korea_reached": "china_border", "all": "seeds"}[dk]
            if src == "gates_terminus":                       # v9.9 FIX D
                # Gates PLUS the unfenced DPRK–ROK contact (the Han estuary west of the fence
                # terminus), which is where the simulated epidemic first entered.
                if not CORRIDORS:
                    return None, "gates requested but none were built"
                _m = np.zeros(SH, bool)
                for _co in CORRIDORS:
                    _m[int(_co["row"]), int(_co["col"])] = True
                if ROK_RASTER is None:
                    logger.warning("    ENTRY: no national split, so the unfenced border cannot be "
                                   "located; using the gates only.")
                    return _m, f"{len(CORRIDORS)} water-deer corridor gates"
                _near_fence = (binary_dilation(MDL_BOUNDARY.astype(bool), iterations=1)
                               if MDL_BOUNDARY is not None else np.zeros(SH, bool))
                _contact = (SK & KOREA_REGION_MODEL & ~ROK_EXCLAVES
                            & binary_dilation(NK & KOREA_REGION, iterations=1) & ~_near_fence)
                _m |= _contact
                return _m, (f"{len(CORRIDORS)} corridor gates + {int(_contact.sum())} unfenced "
                            f"border cells (western terminus / estuary)")
            if src == "gates":
                if not CORRIDORS:
                    return None, "gates requested but none were built"
                _m = np.zeros(SH, bool)
                for _co in CORRIDORS:
                    _m[int(_co["row"]), int(_co["col"])] = True
                return _m, f"{len(CORRIDORS)} water-deer corridor gates"
            if src == "mdl":
                if MDL_BOUNDARY is None or not MDL_BOUNDARY.any():
                    return None, "MDL line requested but not available"
                return MDL_BOUNDARY.astype(bool).copy(), \
                    ("the ENTIRE MDL line — entry risk spread evenly along the fence "
                     "rather than concentrated at the crossings")
            if src == "seeds":
                _m = np.zeros(SH, bool)
                _n = 0
                for _e in (OUTBREAK_ENTRIES or []):
                    _m[int(_e[0]), int(_e[1])] = True; _n += 1
                if _n:
                    return _m, f"{_n} outbreak seed points (the RAD introduction sites)"
                # fall through to the border when no seed layer was loaded
                src = "china_border"
            if src == "china_border":
                try:
                    _bcm = BORDER_CN_MASK          # v9.10: not yet defined during the first SMCDA pass
                except NameError:
                    _bcm = None
                if _bcm is not None and _bcm.any():
                    return _bcm.astype(bool).copy(), \
                        "the China–Korea border zone (introduction into the peninsula)"
                # the northernmost North-Korean cells: the same construction the
                # Korea leg uses to seed its introduction from China
                _m = np.zeros(SH, bool)
                _rows = np.where(NK.any(axis=1))[0]
                if _rows.size:
                    _top = int(_rows.min())
                    for _r in range(_top, min(_top + 3, SH[0])):
                        _m[_r, :] = NK[_r, :]
                if _m.any():
                    return _m, ("the northernmost North-Korean cells (introduction from "
                                "China; BORDER_CN_MASK was unavailable)")
                return None, "no China-border cells could be identified"
            raise ValueError(f"unknown SMCDA_ENTRY_SOURCE '{Config.SMCDA_ENTRY_SOURCE}'")

        def _smcda_gate_distance(entry_mask, dom_mask, label):
            """Cost-distance (km) from the entry origin, on the model's own
            conductance-based cost surface. Unreachable cells inside the domain
            are set to the maximum REACHABLE cost — 'as far from the entry as
            anything the model can reach' — never dropped and never given a
            fabricated distance. The count is logged."""
            if entry_mask is None or not entry_mask.any():
                return None
            _gd = mcp_dist_km_multisource(entry_mask)
            _fin = np.isfinite(_gd) & dom_mask
            if not _fin.any():
                return None
            _gmax = float(np.nanmax(_gd[_fin]))
            _n_un = int((dom_mask & ~_fin).sum())
            if _n_un:
                logger.info(f"    ENTRY: {_n_un:,} domain cells unreachable from "
                            f"{label}; each set to the maximum reachable cost-distance "
                            f"({_gmax:.0f} km) rather than dropped.")
            return np.where(_fin | ~dom_mask, np.where(np.isfinite(_gd), _gd, _gmax),
                            _gmax).astype(np.float32)

        def compute_smcda_for_domain(dk):
            """Full SMCDA chain for one domain. Returns the result dict or None."""
            _mask, _dom_lab = _smcda_domain_mask(dk)
            _n = int(_mask.sum())
            logger.info(f"\n  ── domain '{dk}': {_dom_lab} — {_n:,} cells ──")
            _floor = (int(Config.SMCDA_REACHED_MIN_CELLS)
                      if str(dk).lower() in SMCDA_EPIDEMIC_DOMAINS else 50)
            if _n < _floor:
                logger.warning(f"    ⚠ only {_n} cells (minimum {_floor}) — too few to "
                               f"threshold; the breaks would be noise. Domain '{dk}' "
                               f"skipped.")
                return None
            if str(dk).lower() in SMCDA_EPIDEMIC_DOMAINS:
                logger.warning(
                    "    ⚠ This domain is DEFINED BY THE EPIDEMIC OUTCOME. The score "
                    "values are still computed from static layers alone, but the extent "
                    "is conditional on the simulated incursion, so: AUC against 'reached' "
                    "is undefined here (every cell is reached by construction); the class "
                    "thresholds are NOT comparable with the unrestricted domain's; and "
                    "this is a CONDITIONAL risk map — 'given the incursion got this far' "
                    "— not a prospective one. Arrival-order monotonicity and Spearman rho "
                    "remain valid, and are the harder test.")
            _emask, _ebasis = _smcda_entry_mask(dk)
            logger.info(f"    ENTRY origin: {_ebasis}")
            _gdist = _smcda_gate_distance(_emask, _mask, _ebasis)
            if _gdist is None and str(Config.SMCDA_ENTRY_SOURCE).lower() != "off":
                logger.warning("    ⚠ ENTRY criterion DROPPED (no usable origin). The "
                               "remaining weights are renormalised, so the model actually "
                               "applied is NOT the one specified in Config.")
            _pool = {
                "livestock_nn": livestock_nn, "GOR_EFF": GOR_EFF,
                "GORANI_DENSITY": GORANI_DENSITY, "CONDUCTANCE": CONDUCTANCE,
                "RESISTANCE": RESISTANCE, "WATER_ATTRACT": WATER_ATTRACT,
                "SK_CONNECTIVITY": SK_CONNECTIVITY, "NK_CONNECTIVITY": NK_CONNECTIVITY,
                "MAXENT": MAXENT, "GOR_nn": GOR_nn,
                "_SMCDA_GATE_DIST": _gdist, "_SMCDA_BETA": SMCDA_BETA,
            }
            _objs, _layers, _dropped = [], {}, []
            for _cd in Config.SMCDA_CRITERIA:
                _arr = _pool.get(_cd["layer"])
                _objs.append(SMCDACriterion(
                    key=_cd["key"], label=_cd.get("label", _cd["key"]),
                    direction=_cd.get("direction", "benefit"),
                    standardise=_cd.get("standardise", "rank"),
                    a=_cd.get("a"), b=_cd.get("b"),
                    veto_below=_cd.get("veto_below"), note=_cd.get("layer", "")))
                if _arr is None:
                    _dropped.append((_cd["key"], _cd["layer"]))
                else:
                    _layers[_cd["key"]] = np.asarray(_arr, np.float64)
            for _k, _ly in _dropped:
                logger.warning(f"    ⚠ criterion '{_k}' DROPPED: layer '{_ly}' is not "
                               f"available. Remaining weights renormalised.")
            _A = None
            if str(Config.SMCDA_WEIGHT_METHOD).lower() == "ahp":
                if Config.SMCDA_AHP_MATRIX is not None:
                    _A = np.asarray(Config.SMCDA_AHP_MATRIX, float)
                    logger.info("    AHP matrix: from Config.SMCDA_AHP_MATRIX "
                                "(explicit pairwise judgements).")
                else:
                    _A = ahp_matrix_from_priorities(
                        [float(Config.SMCDA_PRIORITY_VECTOR[c.key]) for c in _objs])
                    logger.warning(
                        "    ⚠ AHP matrix CONSTRUCTED from Config.SMCDA_PRIORITY_VECTOR "
                        "(a_ij = round(s_i/s_j) on Saaty's 1-9 scale). These are NOT "
                        "elicited expert judgements. Because the matrix comes from a "
                        "single vector it is near-consistent BY CONSTRUCTION, so the "
                        "consistency ratio below confirms the arithmetic and says "
                        "NOTHING about whether the priorities are epidemiologically "
                        "correct. Replace SMCDA_AHP_MATRIX with your panel's pairwise "
                        "judgements before these classes allocate surveillance effort.")
            _res = run_smcda(
                layers=_layers, criteria=_objs, mask=_mask,
                ahp_matrix=_A, ahp_cr_limit=float(Config.SMCDA_AHP_CR_LIMIT),
                weight_method=str(Config.SMCDA_WEIGHT_METHOD),
                aggregation=str(Config.SMCDA_AGGREGATION),
                owa_andness=float(Config.SMCDA_OWA_ANDNESS),
                n_classes=int(Config.SMCDA_N_CLASSES),
                break_method=str(Config.SMCDA_BREAK_METHOD),
                absolute_breaks=Config.SMCDA_ABSOLUTE_BREAKS,
                class_labels=Config.SMCDA_CLASS_LABELS,
                sensitivity_delta=float(Config.SMCDA_SENSITIVITY_DELTA),
                logger=logger)
            if (str(Config.SMCDA_WEIGHT_METHOD).lower() == "ahp"
                    and _res.get("ahp") is not None and not _res["ahp"]["consistent"]
                    and Config.SMCDA_AHP_CR_FATAL):
                raise RuntimeError(
                    f"AHP consistency ratio {_res['ahp']['CR']:.3f} exceeds the "
                    f"{Config.SMCDA_AHP_CR_LIMIT:.2f} limit and SMCDA_AHP_CR_FATAL=True "
                    f"— revise the pairwise matrix.")
            _out = dict(_res)
            _out.update(domain=dk, mask=_mask, domain_label=_dom_lab,
                        entry_basis=_ebasis, dropped=_dropped)
            for _ln in format_smcda_threshold_table(_out, cell_km=PW / 1000.0):
                logger.info(_ln)
            return _out

        _dlist = [str(d).lower() for d in
                  (list(Config.SMCDA_DOMAINS) if getattr(Config, "SMCDA_DOMAINS", None)
                   else [Config.SMCDA_APPLY_TO_DOMAIN])]
        if Config.SMCDA_ADD_ALL_DOMAIN and "all" not in _dlist:   # v9.10 item 7 (external test)
            _dlist.append("all")
        _deferred = [d for d in _dlist if d in SMCDA_EPIDEMIC_DOMAINS]
        for _dk in [d for d in _dlist if d not in SMCDA_EPIDEMIC_DOMAINS]:
            try:
                _r = compute_smcda_for_domain(_dk)
                if _r is not None:
                    SMCDA_BY_DOMAIN[_dk] = _r
            except Exception as _e:
                logger.warning(f"  ⚠ SMCDA domain '{_dk}' skipped "
                               f"({type(_e).__name__}: {_e}). The run continues.")
        if _deferred:
            logger.info(f"\n  Deferred until the epidemic has run: {_deferred}. These "
                        f"domains are the PPR-invaded footprint, so their extent is not "
                        f"known until the Korea leg finishes.")

        def compute_deferred_smcda_domains():
            """Build the epidemic-conditioned domains. Called AFTER the Korea leg.

            It only ADDS to SMCDA_BY_DOMAIN, which every downstream consumer
            (validation, maps, county tables) reads directly — so nothing has to
            be rebound and the ordering stays easy to follow."""
            for _dk2 in _deferred:
                try:
                    _r2 = compute_smcda_for_domain(_dk2)
                    if _r2 is not None:
                        SMCDA_BY_DOMAIN[_dk2] = _r2
                except Exception as _e2:
                    logger.warning(f"  ⚠ SMCDA domain '{_dk2}' skipped "
                                   f"({type(_e2).__name__}: {_e2}). The run continues.")
        _prim = str(getattr(Config, "SMCDA_PRIMARY_DOMAIN",
                            Config.SMCDA_APPLY_TO_DOMAIN)).lower()
        if _prim in SMCDA_BY_DOMAIN:
            SMCDA = SMCDA_BY_DOMAIN[_prim]
            SMCDA_CLASS = SMCDA["classes"]
            logger.info(f"\n  Primary domain for the county table and the JSON summary: "
                        f"'{_prim}'.")
        elif _prim in _deferred:
            logger.info(f"\n  Primary domain '{_prim}' is epidemic-conditioned and will "
                        f"be built after the Korea leg.")
        elif SMCDA_BY_DOMAIN:
            logger.warning(f"  ⚠ SMCDA_PRIMARY_DOMAIN '{_prim}' is not among the domains "
                           f"requested; the county table will use the first that "
                           f"completed.")
        else:
            logger.warning("  ⚠ no SMCDA domain completed yet.")
        # ---- optional: let the classes drive the transport -------------------
        if Config.USE_SMCDA_CLASS_IN_TRANSPORT and SMCDA_CLASS is not None:
            if _prim != "all":
                raise ValueError(
                    f"USE_SMCDA_CLASS_IN_TRANSPORT=True with SMCDA_PRIMARY_DOMAIN="
                    f"'{_prim}'. A transport field has to be defined over the WHOLE "
                    f"model domain; a '{_prim}'-only classification leaves the rest of "
                    f"the grid unspecified. Set SMCDA_PRIMARY_DOMAIN='all' or turn the "
                    f"switch off.")
            _mobs = list(Config.SMCDA_CLASS_MOBILITY)
            if len(_mobs) != int(Config.SMCDA_N_CLASSES):
                raise ValueError(f"SMCDA_CLASS_MOBILITY has {len(_mobs)} entries but "
                                 f"SMCDA_N_CLASSES={Config.SMCDA_N_CLASSES}")
            SMCDA_MOBILITY = np.zeros(SH, np.float32)
            for _ci, _mv in enumerate(_mobs):
                SMCDA_MOBILITY[np.nan_to_num(SMCDA_CLASS, nan=-1.0) == _ci] = float(_mv)
            SMCDA_MOBILITY = np.where(KP, SMCDA_MOBILITY, 0.0).astype(np.float32)
            logger.warning(
                "  ⚠ USE_SMCDA_CLASS_IN_TRANSPORT=True: the SMCDA classes are now an "
                "INPUT to the epidemic. Two consequences. (1) The simulated spread "
                "differs from the run that produced the unweighted dynamics. (2) The "
                "validation of the classification against that epidemic is no longer "
                "independent — the map is partly validating itself. Report the "
                "validation statistics from a run with this switch OFF.")
    else:
        logger.info("\n[SPATIAL MCDA] disabled (Config.ENABLE_SMCDA=False).")
    # ── 7. Entry points + seasonal drivers ───────────────────────────────────
    logger.info("\n[VALIDATING LEGACY ENTRY POINTS]")
    if Config.AUTO_VALIDATE_ENTRIES:
        for ex, ey, label, _off, wt in Config.DMZ_ENTRIES:
            r_, c_ = rc_from_km(ex, ey)
            g = float(GOR_nn[r_, c_])
            flag = "" if g > 0.05 else "  ⚠ near-zero habitat — verify this coordinate"
            logger.info(f"  DMZ_ENTRY {label:<12s} -> ({r_},{c_}) GOR_nn={g:.3f} SK={bool(SK[r_,c_])}{flag}")
        for ex, ey, label, _off, wt in Config.NK_ENTRIES:
            r_, c_ = rc_from_km(ex, ey)
            g = float(GOR_nn[r_, c_])
            flag = "" if g > 0.05 else "  ⚠ near-zero habitat — verify this coordinate"
            logger.info(f"  NK_ENTRY  {label:<12s} -> ({r_},{c_}) GOR_nn={g:.3f} NK={bool(NK[r_,c_])}{flag}")
    T_VIRAL = np.clip(1.5 - 0.036 * Config.TEMP_C, 0.55, 1.55)
    JOINT = Config.MONTHLY_BETA * (Config.DEER_RADIUS_KM / 40.0) * T_VIRAL
    JOINT = JOINT / JOINT.max()
    cell_km = PW / 1000.0
    # ── kernel-decay reality check against Han et al. (2022) ────────────────
    if Config.HAN_LOG_DECAY_COMPARISON:
        alpha = Config.HAN_ALPHA_DECAY_PER_KM
        pi_h = Config.HAN_DIFF_HABITAT_REDUCTION
        logger.info("\n[KERNEL DECAY COMPARISON — ours (sigma) vs Han et al. alpha]")
        logger.info(f"  Han: pi*exp(-alpha*d), alpha={alpha}/km, pi={pi_h} "
                    f"-> ~{100*pi_h*np.exp(-alpha*np.sqrt(5)):.1f}% at 2.24 km, "
                    f"~{100*pi_h*np.exp(-alpha*10):.2g}% at 10 km")
        for mi in (0, 5):  # January (min radius) and June (peak radius)
            sig = float(Config.DEER_RADIUS_KM[mi])
            def _range_to(frac, s):
                return -s * np.log(frac)
            logger.info(
                f"  Ours ({Config.MONTH_SHORT[mi]}, sigma={sig:.1f} km): "
                f"{100*np.exp(-np.sqrt(5)/sig):.1f}% at 2.24 km, "
                f"{100*np.exp(-10/sig):.1f}% at 10 km; "
                f"reaches 10% at {_range_to(0.10, sig):.0f} km "
                f"(Han reaches 10% at ~{-np.log(0.10/pi_h)/alpha:.1f} km)")
        logger.info("  -> v9.1 narrows this gap from the transport side: slower diffusion/advection "
                    "and a gated fence, rather than a broad kernel crossing everywhere.\n")
    _elev_f2 = np.nan_to_num(ELEV, nan=500).astype(np.float64)
    _gy2, _gx2 = np.gradient(_elev_f2, PH, PW)
    SLOPE_N = np.clip(np.degrees(np.arctan(np.sqrt(_gx2**2 + _gy2**2))) / 45.0,
                      0, 1).astype(np.float32)
    # ── 7b. Neural component initialisation (v7.0) ───────────────────────────
    logger.info("\n[NEURAL COMPONENT INITIALISATION]")
    use_neural = Config.USE_NEURAL_COMPONENTS and TORCH_AVAILABLE
    device = "cuda" if (TORCH_AVAILABLE and torch.cuda.is_available()) else "cpu"
    if use_neural:
        kernel_net, ode_func, device = build_neural_components(device)
        logger.info(f"  ✓ SpatialKernelNet + ODEFunc built on device={device}")
        logger.info(f"    SpatialKernelNet params: {sum(p.numel() for p in kernel_net.parameters()):,}")
        logger.info(f"    ODEFunc params:          {sum(p.numel() for p in ode_func.parameters()):,}")
        if not TORCHDIFFEQ_AVAILABLE:
            logger.warning("  ⚠ torchdiffeq not found — Neural ODE will use Euler fallback. "
                           "Install: pip install torchdiffeq --break-system-packages")
    else:
        kernel_net, ode_func = None, None
        if not TORCH_AVAILABLE:
            logger.info("  PyTorch not available — running physics-only (v6.0 behaviour)")
        else:
            logger.info("  USE_NEURAL_COMPONENTS=False — running physics-only")
    outbreak_df = None
    if use_neural and Config.TRAIN_NEURAL and Config.OUTBREAK_CSV.exists():
        try:
            import pandas as pd
            outbreak_df = pd.read_csv(Config.OUTBREAK_CSV)
            logger.info(f"  Outbreak CSV loaded: {len(outbreak_df)} records from {Config.OUTBREAK_CSV}")
        except Exception as e:
            logger.warning(f"  ⚠ Could not load outbreak CSV: {e}")
    elif use_neural and Config.TRAIN_NEURAL:
        logger.warning(f"  ⚠ TRAIN_NEURAL=True but OUTBREAK_CSV not found: {Config.OUTBREAK_CSV}")
    # ══ KOREA CORRIDOR REGIONS (China → North Korea → South Korea) ═══════════
    KOREA_MASK = KP & np.isfinite(arr5) & ~ROK_EXCLAVES_DROP          # v9.10: exclave policy
    if not KOREA_MASK.any():
        KOREA_MASK = KP & (GOR_nn > 0)
    NK_KR = KOREA_MASK & NK                    # v9.9 FIX C: DPRK (national boundary)
    SK_KR = KOREA_MASK & SK                    # v9.9 FIX C: ROK (national boundary)
    CHINA_MASK = KP & (~KOREA_MASK) & ~ROK_EXCLAVES_DROP   # China landmass (the source region)
    logger.info(f"\n[KOREA CORRIDOR] China={int(CHINA_MASK.sum()):,} cells, "
                f"Korea={int(KOREA_MASK.sum()):,} (NK={int(NK_KR.sum()):,}, "
                f"SK={int(SK_KR.sum()):,}); spread path = China → NK → (corridor gates) → SK.")
    if getattr(Config, "HAB_PENINSULA_REPORT", False):
        # 11 Nov: the same classification as section 5c, over the peninsular domain only
        _kp = KOREA_MASK.astype(bool)
        _gek = np.nan_to_num(GOR_EFF, nan=0.0)
        _thr_k = float(np.nanquantile(_gek[_kp], Config.HAB_SUIT_QUANTILE)) if _kp.any() else float("nan")
        _hab_k = _kp & (_gek >= _thr_k)
        _dist_k = mcp_dist_km_multisource(_hab_k) if _hab_k.any() else np.full(SH, np.inf, np.float32)
        _acc_k = _kp & ~_hab_k & np.isfinite(_dist_k) & (_dist_k <= Config.ACCESS_COSTDIST_MAX_KM)
        _bar_k = _kp & ~_hab_k & ~_acc_k
        _nk = max(int(_kp.sum()), 1)
        _res = dict(domain="peninsular (Korean Peninsula, KOREA_MASK)", cells=int(_kp.sum()),
                    quantile=float(Config.HAB_SUIT_QUANTILE), threshold=round(_thr_k, 4),
                    habitat_pct=round(100.0 * _hab_k.sum() / _nk, 1),
                    accessible_pct=round(100.0 * _acc_k.sum() / _nk, 1),
                    barrier_pct=round(100.0 * _bar_k.sum() / _nk, 1),
                    reach_km=float(Config.ACCESS_COSTDIST_MAX_KM),
                    barrier_cells_in_DPRK=int((_bar_k & NK).sum()), barrier_cells_in_ROK=int((_bar_k & SK).sum()))
        logger.info(f"\n[HABITAT CLASSES — PENINSULAR DOMAIN] threshold GOR_EFF>={_thr_k:.3f} (45th percentile); "
                    f"habitat {_res['habitat_pct']}%, accessible {_res['accessible_pct']}%, "
                    f"barrier {_res['barrier_pct']}% of {_nk:,} cells")
        Config.OUT_DIR.mkdir(parents=True, exist_ok=True)
        _fj = Config.OUT_DIR / f"{V910_TAG}_habitat_peninsula.json"
        _fj.write_text(_v910_json.dumps(_res, indent=1), encoding="utf-8")
        try:
            _rr, _cc = np.where(_kp)
            _r0, _r1 = max(0, _rr.min() - 3), min(NR, _rr.max() + 4)
            _c0, _c1 = max(0, _cc.min() - 3), min(NC, _cc.max() + 4)
            _cls = np.full(SH, np.nan, np.float32)
            _cls[_bar_k] = 0; _cls[_acc_k] = 1; _cls[_hab_k] = 2
            fig, ax = plt.subplots(figsize=(8, 10), dpi=150, facecolor="white")
            ax.imshow(np.ma.masked_invalid(_cls[_r0:_r1, _c0:_c1]),
                      cmap=ListedColormap([Config.CLASS_COLOR_BARRIER, Config.CLASS_COLOR_ACCESSIBLE,
                                           Config.CLASS_COLOR_HABITAT]), vmin=-0.5, vmax=2.5,
                      extent=(0, _c1 - _c0, _r1 - _r0, 0), origin="upper", aspect="equal",
                      interpolation="nearest", zorder=1)
            _m = MDL_BOUNDARY[_r0:_r1, _c0:_c1]
            if _m.any():
                ax.imshow(np.ma.masked_where(~_m, np.ones_like(_m, float)), cmap=ListedColormap(["#d81b60"]),
                          extent=(0, _c1 - _c0, _r1 - _r0, 0), origin="upper", aspect="equal", zorder=3)
            if CORRIDORS:
                plot_corridor_markers(ax, CORRIDORS, row_offset=_r0, col_offset=_c0, label=False,
                                      size_scale=30.0)
            ax.set_xlim(0, _c1 - _c0); ax.set_ylim(_r1 - _r0, 0); ax.set_xticks([]); ax.set_yticks([])
            _bar = 100.0 / (PW / 1000.0)
            _x0, _y0 = (_c1 - _c0) * 0.60, (_r1 - _r0) * 0.95
            ax.plot([_x0, _x0 + _bar], [_y0, _y0], "k-", lw=3, zorder=9)
            ax.text(_x0 + _bar / 2, _y0 - (_r1 - _r0) * 0.012, "100 km", ha="center", va="bottom",
                    fontsize=10, fontweight="bold")
            ax.legend(handles=[Patch(facecolor=Config.CLASS_COLOR_HABITAT, label=f"Habitat ({_res['habitat_pct']:.0f}%)"),
                               Patch(facecolor=Config.CLASS_COLOR_ACCESSIBLE,
                                     label=f"Accessible ({_res['accessible_pct']:.0f}%)"),
                               Patch(facecolor=Config.CLASS_COLOR_BARRIER,
                                     label=f"Barrier ({_res['barrier_pct']:.0f}%)"),
                               Line2D([0], [0], color="#d81b60", lw=2, label="Military Demarcation Line")],
                      loc="upper left", fontsize=10, frameon=True, facecolor="white", edgecolor="#888")
            ax.annotate("N", xy=(0.92, 0.95), xytext=(0.92, 0.85), xycoords="axes fraction", ha="center",
                        va="center", fontsize=12, fontweight="bold",
                        arrowprops=dict(facecolor="black", edgecolor="black", width=3, headwidth=10))
            _fp2 = Config.OUT_DIR / "PPR_v9.10_FigS2_habitat_peninsula.png"
            fig.savefig(_fp2, dpi=int(Config.PUB_FIGS_DPI), bbox_inches="tight", facecolor="white")
            plt.close(fig)
            for _f in (_fj, _fp2):
                if Config.DL_DIR:
                    try:
                        Path(Config.DL_DIR).mkdir(parents=True, exist_ok=True)
                        shutil.copy2(_f, Path(Config.DL_DIR) / _f.name)
                    except Exception:
                        pass
            logger.info(f"   ✓ {_fj.name}, {_fp2.name}")
        except Exception as _e:
            logger.warning(f"   ⚠ Figure S2 not drawn: {_e}")
        logger.info("v9.10: habpeninsula preset — classification written; stopping here by design.")
        return
    # v9.2: the China cells closest to Korea. The source run is now driven until
    # the epidemic actually REACHES these (that month becomes the North-Korean
    # introduction date), instead of stopping at a fixed year.
    BORDER_CN_MASK = None
    if CHINA_MASK.any() and KOREA_MASK.any():
        _d2k = distance_transform_edt(~KOREA_MASK)
        _cdk = _d2k[CHINA_MASK]
        _thr_dk = np.percentile(_cdk, 5) if _cdk.size else 0.0
        BORDER_CN_MASK = CHINA_MASK & (_d2k <= max(_thr_dk, 1.0))
        logger.info(f"  China cells adjacent to Korea (introduction gateway): "
                    f"{int(BORDER_CN_MASK.sum()):,}")
    # ══ VACCINATION TARGET FIELD ═════════════════════════════════════════════
    VACC_FIELD = None
    if Config.ENABLE_VACCINATION:
        _reg = Config.VACC_TARGET_REGION
        _tgt = {"korea": KOREA_MASK, "nk": NK_KR, "sk": SK_KR}.get(_reg, KP)
        VACC_FIELD = np.where(_tgt & KP, 1.0, 0.0).astype(np.float32)
        logger.info(f"[VACCINATION] mode={Config.VACC_MODE}, region={_reg} "
                    f"({int((_tgt & KP).sum()):,} eligible cells), start month "
                    f"{Config.VACC_START_MONTH + 1}, coverage≤{Config.VACC_COVERAGE:.0%}, "
                    f"efficacy={Config.VACC_EFFICACY:.0%}.")
    # ══ DYNAMIC SPREAD ENGINE (Reaction-Advection-Diffusion) ═════════════════
    rad_results_holder = {}
    KOREA_LEG: Dict[str, object] = {}          # v9.10: the peninsular leg, for the new reports
    V910_STATE: Dict[str, object] = {}         # v9.10: baseline metrics + batch context
    def run_dynamic_spread_simulation():
        """Run the mechanistic RAD epidemic. S/I/R persist and each month evolves
        from the previous one (no re-seeding). Returns the four arrays the
        visualization consumes — monthly infected density (I) and monthly
        cumulative attack rate (I+R), each split into NK / SK views."""
        cell_km_local = abs(GT5[1]) / 1000.0
        def _entries_to_rc(entries):
            out = []
            for e in entries:
                ex, ey = e[0], e[1]
                r_, c_ = rc_from_km(ex, ey)
                lbl = e[2] if len(e) > 2 else ""
                off = e[3] if len(e) > 3 else 0
                wt = float(e[4]) if len(e) > 4 else 1.0
                out.append((int(r_), int(c_), lbl, off, wt))
            return out
        dmz_rc = _entries_to_rc(Config.DMZ_ENTRIES)
        nk_rc = _entries_to_rc(Config.NK_ENTRIES)
        if (Config.SEED_FROM_OUTBREAK_SHP or getattr(Config, "SEED_TABLE", None)) and OUTBREAK_ENTRIES:
            seed_dmz, seed_nk = list(OUTBREAK_ENTRIES), []
            logger.info(f"  RAD seeding: {len(seed_dmz)} OUTBREAK "
                        f"{'seed-table' if V910_SEED_INFO.get('source') == 'table' else 'shapefile'} points "
                        f"(spread starts here; DMZ/NK entries bypassed).")
        elif Config.CORRIDOR_SEED_FALLBACK and CORRIDORS:
            seed_dmz = corridor_seed_entries(CORRIDORS, NK)
            seed_nk = nk_rc
            logger.info(f"  RAD seeding: {len(seed_dmz)} DMZ CORRIDOR gates (weighted) "
                        f"+ {len(nk_rc)} NK entries.")
        else:
            seed_dmz, seed_nk = dmz_rc, nk_rc
            logger.info(f"  RAD seeding: {len(dmz_rc)} legacy DMZ + {len(nk_rc)} NK entries "
                        f"converted to grid indices (grid {NR}x{NC}).")
        env = {
            "KP": KP, "NK": NK, "SK": SK,
            "MDL_BOUNDARY": MDL_BOUNDARY,
            "MDL_FENCE_FIELD": MDL_FENCE_FIELD,     # v9.1: corridor gates
            "GOR_EFF": GOR_EFF, "GOR_nn": GOR_nn,
            "CONDUCTANCE": CONDUCTANCE, "RESISTANCE": RESISTANCE,
            "RIV_nn": RIV_nn, "WATER_ATTRACT": WATER_ATTRACT,
            "livestock_nn": livestock_nn, "GORANI_DENSITY": GORANI_DENSITY,
            "DMZ_ENTRIES": seed_dmz, "NK_ENTRIES": seed_nk,
        }
        if ROAD_MASK is not None and np.any(ROAD_MASK):
            env["ROAD_MASK"] = ROAD_MASK
        if Config.ENABLE_VACCINATION and VACC_FIELD is not None:
            env["VACC_FIELD"] = VACC_FIELD
        if Config.USE_HABITAT_CLASS_IN_TRANSPORT:
            env["CLASS_MOBILITY"] = CLASS_MOBILITY
            logger.info("  Habitat classification is DRIVING transport "
                        "(barrier cells near-impassable): USE_HABITAT_CLASS_IN_TRANSPORT=True")
        # v9.7: the SMCDA classes may drive transport INSTEAD. Both switches on
        # at once is a specification error, not something to silently resolve —
        # the two classifications disagree about which cells are passable.
        if Config.USE_SMCDA_CLASS_IN_TRANSPORT:
            if Config.USE_HABITAT_CLASS_IN_TRANSPORT:
                raise ValueError(
                    "USE_HABITAT_CLASS_IN_TRANSPORT and USE_SMCDA_CLASS_IN_TRANSPORT "
                    "are both True. They supply competing CLASS_MOBILITY fields; "
                    "choose one.")
            if SMCDA_MOBILITY is None:
                raise ValueError(
                    "USE_SMCDA_CLASS_IN_TRANSPORT=True but the SMCDA classification "
                    "did not produce a mobility field (see the SMCDA warnings above). "
                    "Refusing to run the epidemic with an unspecified transport field.")
            env["CLASS_MOBILITY"] = SMCDA_MOBILITY
            logger.warning("  SMCDA risk classes are DRIVING transport: "
                           "USE_SMCDA_CLASS_IN_TRANSPORT=True")
        p = RADParams(
            beta_monthly=np.asarray(Config.MONTHLY_BETA, float),
            disp_km_monthly=np.asarray(Config.DEER_RADIUS_KM, float) * Config.DEER_RADIUS_SCALE,
            tvir_monthly=np.asarray(T_VIRAL, float),
            cell_km=cell_km_local,
            diff_max=Config.RAD_DIFF_MAX * Config.SPREAD_SPEED_MULTIPLIER,
            diff_min=Config.RAD_DIFF_MIN,
            adv_max=Config.RAD_ADV_MAX * Config.SPREAD_SPEED_MULTIPLIER,
            adv_gain=Config.RAD_ADV_GAIN * Config.SPREAD_SPEED_MULTIPLIER,
            gamma=Config.GAMMA_RECOVERY,
            trans_hab_w=Config.RAD_TRANS_HAB_W, trans_lstk_w=Config.RAD_TRANS_LSTK_W,
            trans_base=Config.RAD_TRANS_BASE,
            seed_intensity=Config.RAD_SEED_INTENSITY,
            source_continuous=Config.RAD_SOURCE_CONTINUOUS,
            source_rate=Config.RAD_SOURCE_RATE,
            mdl_fence_factor=Config.HAN_FENCE_REDUCTION,
            apply_mdl_barrier=Config.HAN_APPLY_DMZ_FENCE,
            road_factor=Config.HAN_ROAD_REDUCTION,
            apply_road_barrier=Config.HAN_APPLY_ROAD_BARRIER,
            corridor_seasonal_gating=Config.CORRIDOR_SEASONAL_GATING,
            corridor_season_floor=Config.CORRIDOR_SEASON_FLOOR,
            cfl_safety=Config.RAD_CFL_SAFETY, max_substeps=Config.RAD_MAX_SUBSTEPS,
            arrival_threshold=Config.RAD_ARRIVAL_THRESHOLD,
            vacc_enable=Config.ENABLE_VACCINATION,
            vacc_mode=Config.VACC_MODE,
            vacc_start_month=int(Config.VACC_START_MONTH),
            vacc_coverage=float(Config.VACC_COVERAGE),
            vacc_rate=float(Config.VACC_RATE),
            vacc_efficacy=float(Config.VACC_EFFICACY),
            vacc_ring_cells=max(1, int(round(Config.VACC_RING_KM / max(cell_km_local, 1e-6)))),
            vacc_inf_thr=float(Config.VACC_INFECTED_THRESHOLD),
        )
        n_months = 12 * max(1, int(Config.SIM_YEARS))
        _stop_fn = _stop_extra = _cap = None
        _lab = "criterion"
        if Config.AUTO_RUN_UNTIL_SK_ESTABLISHED and BORDER_CN_MASK is not None \
                and BORDER_CN_MASK.any():
            _bmask = BORDER_CN_MASK
            def _stop_fn(_m, _cum, _bm=_bmask):
                return float(np.nanmean(_cum[_bm])) >= Config.RAD_ARRIVAL_THRESHOLD
            _stop_extra = 12 * max(0, int(Config.CHINA_POST_BORDER_YEARS))
            _cap = 12 * max(int(Config.SIM_YEARS), int(Config.MAIN_MAX_YEARS))
            _lab = "China epidemic reaching the North-Korean border"
        if Config.KOREA_ONLY_INTRO_MONTH:                   # v9.10: peninsular-leg-only run
            n_months = max(12, int(Config.KOREA_ONLY_INTRO_MONTH))
            _stop_fn, _stop_extra, _cap = None, 0, None
            logger.warning(f"  v9.10 KOREA-ONLY MODE: the continental leg runs {n_months} months only and "
                           f"the DPRK introduction is fixed at month {int(Config.KOREA_ONLY_INTRO_MONTH)}. "
                           f"Continental figures from this run are not for reporting.")
        logger.info(f"\n[DYNAMIC RAD EPIDEMIC — minimum {n_months} months "
                    f"({Config.SIM_YEARS} yr)"
                    f"{f', extending until {_lab} (cap {_cap // 12} yr)' if _stop_fn else ''}; "
                    f"S/I/R persists month-to-month]")
        res = run_dynamic_rad_simulation(env, p, n_months=n_months, closure=None,
                                         stop_fn=_stop_fn,
                                         stop_extra_months=_stop_extra or 0,
                                         max_months=_cap, stop_label=_lab)
        n_months = len(res["monthly_cum"])          # actual, after auto-extension
        # ── FRONT-SPEED DIAGNOSTIC (v9.2) ───────────────────────────────────
        # The single number to calibrate the run length against. Compare it to
        # an observed wildlife-disease front in the same landscape (e.g. the ASF
        # wild-boar wave of Han et al. 2022) and scale SPREAD_SPEED_MULTIPLIER
        # until they agree — then set SIM_YEARS from the distance to be covered.
        _fk = np.asarray(res["time_series"]["front_km"], float)
        if _fk.size > 13:
            _yrs = (_fk.size - 12) / 12.0
            _spd = (float(_fk[-1]) - float(_fk[11])) / max(_yrs, 1e-9)
            logger.info(
                f"  [FRONT SPEED] {_spd:.0f} km/yr mean after year 1 "
                f"(front {_fk[11]:.0f} -> {_fk[-1]:.0f} km over {_yrs:.1f} yr). "
                f"Transport caps: adv={p.adv_max:.2f} cells/mo "
                f"({p.adv_max * p.cell_km:.1f} km/mo), diff={p.diff_max:.2f} cell²/mo "
                f"(RMS {2.0 * math.sqrt(max(p.diff_max, 0.0)) * p.cell_km:.1f} km/mo) "
                f"vs an empirical dispersal radius of "
                f"{float(np.mean(p.disp_km_monthly)):.1f} km/mo. "
                f"SPREAD_SPEED_MULTIPLIER={Config.SPREAD_SPEED_MULTIPLIER:.2f}.")
            if _spd > 0 and _spd * Config.SIM_YEARS < 300:
                logger.warning(
                    f"  ⚠ At {_spd:.0f} km/yr the front covers only "
                    f"{_spd * Config.SIM_YEARS:.0f} km in the {Config.SIM_YEARS}-year run — "
                    f"less than the ~300 km length of South Korea. Raise "
                    f"SPREAD_SPEED_MULTIPLIER or SIM_YEARS, or lower CONDUCTANCE_EXPONENT.")
        rad_results_holder["res"] = res
        rad_results_holder["env"] = env
        rad_results_holder["params"] = p
        mnk, msk, cnk, csk = [], [], [], []
        for m in range(n_months):
            I = res["monthly_I"][m]
            cum = res["monthly_cum"][m]
            mnk.append(np.where(NK, I, np.nan).astype(np.float32))
            msk.append(np.where(SK, I, np.nan).astype(np.float32))
            cnk.append(np.where(NK, cum, np.nan).astype(np.float32))
            csk.append(np.where(SK, cum, np.nan).astype(np.float32))
        logger.info(f"✓ Dynamic RAD epidemic computed ({n_months} months)")
        return mnk, msk, cnk, csk
    (monthly_nk, monthly_sk,
     monthly_cum_nk, monthly_cum_sk) = run_dynamic_spread_simulation()
    if use_neural and Config.TRAIN_NEURAL and kernel_net is not None:
        logger.info("\n[NEURAL COMPONENT TRAINING]")
        logger.info("  Note: full training integration requires lambda caching across months.")
    results_by_mode = {"v9_3_corridor_gated": (monthly_nk, monthly_sk,
                                               monthly_cum_nk, monthly_cum_sk)}
    # ── 9. Pre-rasterize geo layers ──────────────────────────────────────────
    logger.info("\n[RASTERIZING GEO LAYERS]")
    def rasterize_polys(polys, nr, nc, label):
        raster = np.zeros((nr, nc), dtype=np.uint8)
        if polys:
            for xs, ys in polys:
                if np.any((xs > 0) & (xs < nc) & (ys > 0) & (ys < nr)):
                    xs_int = np.clip(xs.astype(int), 0, nc - 1)
                    ys_int = np.clip(ys.astype(int), 0, nr - 1)
                    raster[ys_int, xs_int] = 1
        logger.info(f"  {label}: {raster.sum():,} pixels from {len(polys)} polygons")
        return raster
    water_raster = rasterize_polys(WATER_POLYS, NR, NC, "Water")
    wetland_raster = rasterize_polys(WETLAND_POLYS, NR, NC, "Wetland")
    county_raster = rasterize_polys(SIG_POLYS, NR, NC, "Counties")
    # ── 10. Visualization ────────────────────────────────────────────────────
    logger.info("\n[RENDERING FULL VISUALIZATION]")
    CMAP_NK = LinearSegmentedColormap.from_list(
        "nk", ["#fff9e6", "#ffe082", "#ffa000", "#e65100", "#b71c1c", "#7f0000"], N=512)
    CMAP_SK = LinearSegmentedColormap.from_list(
        "sk", ["#f0f8ff", "#c5dff8", "#6baed6", "#2171b5", "#084594", "#041b5e"], N=512)
    PCT_BANDS = list(Config.SPREAD_PCT_BANDS)
    BAND_NORM = BoundaryNorm(PCT_BANDS, ncolors=len(PCT_BANDS) - 1, clip=True)
    CMAP_BAND_NK = ListedColormap(Config.SPREAD_BAND_COLORS_NK, name="nk_bands")
    CMAP_BAND_SK = ListedColormap(Config.SPREAD_BAND_COLORS_SK, name="sk_bands")
    CMAP_BAND_ALL = ListedColormap(["#ffffcc", "#fed976", "#fd8d3c", "#e31a1c", "#800026"],
                                   name="all_bands")
    _dens = GORANI_DENSITY
    if KP.any():
        _dl = _dens[KP & np.isfinite(_dens)]
        _gor_thr = float(np.nanquantile(_dl, Config.GORANI_DOT_QUANTILE)) if _dl.size else 1.0
    else:
        _gor_thr = 1.0
    GOR_HI = (KP & np.isfinite(_dens) & (np.nan_to_num(_dens, nan=-1) >= _gor_thr))
    _gor_hi_rows, _gor_hi_cols = np.where(GOR_HI)
    logger.info(f"  Gorani spreader dots: composite density >= {_gor_thr:.3f} "
                f"(q{Config.GORANI_DOT_QUANTILE:.2f}) -> {GOR_HI.sum():,} high-density cells "
                f"[NK={int((GOR_HI & NK).sum()):,}, SK={int((GOR_HI & SK).sum()):,}]")
    CMAP_RISK_NK = LinearSegmentedColormap.from_list(
        "risk_nk", ["#ffffff00", "#f3e6d8", "#e3c39c", "#cf9b63"], N=256)
    CMAP_RISK_SK = LinearSegmentedColormap.from_list(
        "risk_sk", ["#ffffff00", "#dfeadf", "#b8d4bf", "#8fb99b"], N=256)
    # FIX (v9.1): the NK risk underlay array was referenced but never built, which
    # raised NameError as soon as SINGLE_MODE was False (BASE_RASTER="gorani").
    RISK_NK_LAYER = np.ma.masked_invalid(
        np.where(NK, np.nan_to_num(GORANI_DENSITY, nan=np.nan), np.nan))
    RISK_SK_LAYER = np.ma.masked_invalid(
        np.where(SK, np.nan_to_num(GORANI_DENSITY, nan=np.nan), np.nan))
    sig_segments = [np.column_stack([cs, rs]) for (cs, rs) in SIG_POLYS] if Config.SHOW_VECTOR_SIG else []
    water_segments = [np.column_stack([cs, rs]) for (cs, rs) in WATER_POLYS] if Config.SHOW_VECTOR_WATER else []
    logger.info(f"  Vector overlays: sig.shp lines={len(sig_segments)}, water lines={len(water_segments)}")
    def mc(c):
        c = c.copy()
        try:
            c.set_bad(alpha=0)
        except Exception:
            pass
        return c
    def season_of(m):
        for s, ms in Config.SEASONS.items():
            if m in ms:
                return s
        return ""
    Config.OUT_DIR.mkdir(parents=True, exist_ok=True)
    # ══ v9.3: SOUTH-KOREA COUNTY CSV EXPORT DRIVER ═══════════════════════════
    # One function, called twice: once for the Korea (cross-MDL) leg — which is
    # the run that actually delivers PPR into South Korea — and once for the
    # continental run. Everything it needs is already on the model grid.
    sk_county_outputs: Dict[str, List[Path]] = {}
    sk_county_tables: Dict[str, tuple] = {}     # v9.10: (rows, series) per tag, reused by Table 1
    def save_epidemic_cache(cum_list, domain_mask, tag, row_offset, col_offset,
                            month_offset, sk_mask=None, monthly_I=None, extra=None):
        """Write the attack-rate stack + everything the county export needs to a
        single .npz, so the table can be rebuilt (or the zonal settings retuned)
        without re-running the simulation. Saved BEFORE the export is attempted."""
        if not Config.SAVE_EPIDEMIC_CACHE:
            return None
        try:
            import json as _json
            fp = Config.OUT_DIR / f"PPR_v9.3_epidemic_cache_{tag}.npz"
            payload = dict(
                cum=np.stack([np.nan_to_num(np.asarray(c, np.float32), nan=0.0)
                              for c in cum_list]),
                domain_mask=np.asarray(domain_mask, bool),
                row_offset=np.int64(row_offset), col_offset=np.int64(col_offset),
                month_offset=np.int64(month_offset),
                gt=np.asarray(GT5, np.float64),
                nr=np.int64(NR), nc=np.int64(NC),
                wkt=np.array(WKT), cell_km=np.float64(PW / 1000.0),
                tag=np.array(tag),
                corridors=np.array(_json.dumps(
                    [{k: v for k, v in co.items() if k in
                      ("name", "weight", "row", "col", "lat", "lon")} for co in CORRIDORS])),
                meta=np.array(_json.dumps(dict(
                    base_raster=Config.BASE_RASTER,
                    arrival_threshold=float(Config.RAD_ARRIVAL_THRESHOLD),
                    sim_start_year=Config.SIM_START_YEAR,
                    sim_start_month=Config.SIM_START_MONTH,
                    n_months=len(cum_list),
                    vaccination=(Config.VACC_MODE if Config.ENABLE_VACCINATION else None),
                    **(extra or {})))),
            )
            if sk_mask is not None:
                payload["sk_mask"] = np.asarray(sk_mask, bool)
            if monthly_I and Config.CACHE_INCLUDE_INFECTED:
                payload["inf"] = np.stack([np.nan_to_num(np.asarray(c, np.float32), nan=0.0)
                                           for c in monthly_I])
            np.savez_compressed(fp, **payload)
            logger.info(f"   ✓ epidemic cache: {fp} "
                        f"({fp.stat().st_size / 1e6:.1f} MB, {len(cum_list)} months) — "
                        f"rebuild the county CSVs from this with "
                        f"`python sk_county_tools.py --rebuild {fp.name}`")
            return fp
        except Exception as e:
            logger.warning(f"   ⚠ epidemic cache not written ({type(e).__name__}: {e}); "
                           f"a failed county export would then need a full re-run.")
            return None
    def export_sk_counties(cum_list, domain_mask, tag, title,
                           row_offset=0, col_offset=0, month_offset=0,
                           sk_mask=None, monthly_I=None, vacc_note=None):
        """Write the per-county CSVs for one epidemic run. Returns the list of
        files written (empty if the export could not run)."""
        if not (Config.EXPORT_SK_COUNTY_CSV and SK_COUNTY_ZONES and cum_list):
            return []
        try:
            thr = (Config.RAD_ARRIVAL_THRESHOLD
                   if Config.SK_COUNTY_ARRIVAL_THRESHOLD is None
                   else float(Config.SK_COUNTY_ARRIVAL_THRESHOLD))
            logger.info(f"\n[SOUTH-KOREA COUNTY CSV — {title}]")
            rows, series = county_attack_table(
                SK_COUNTY_ZONES, cum_list, domain_mask, cell_km=PW / 1000.0,
                arrival_threshold=thr, month_offset=month_offset,
                row_offset=row_offset, col_offset=col_offset, sk_mask=sk_mask,
                monthly_I=monthly_I, corridors=CORRIDORS,
                oversample=SK_ZONE_META.get("oversample"),
                start_year=Config.SIM_START_YEAR, start_month=Config.SIM_START_MONTH,
                logger=logger)
            sk_county_tables[tag] = (rows, series)                   # v9.10
            meta = [
                f"PPR v9.3 South-Korea county (sigungu) attack table — {title}",
                f"county layer: {SK_COUNTY_SHP.name if SK_COUNTY_SHP else 'n/a'} "
                f"(DBF encoding: {SK_ZONE_META.get('encoding') or 'GDAL default'})",
                f"model grid: {NC}x{NR} cells at {PW/1000.0:.2f} km "
                f"(BASE_RASTER={Config.BASE_RASTER}); zonal oversample "
                f"{SK_ZONE_META.get('oversample')}x",
                f"arrival threshold: cumulative attack rate >= {100*thr:.1f}% "
                f"(RAD_ARRIVAL_THRESHOLD)",
                f"months simulated in this leg: {len(cum_list)}; reported months are "
                f"offset by {month_offset} so all runs share one timeline",
                f"MDL: sealed={Config.MDL_SEALED_FACTOR}, gate peak="
                f"{Config.CORRIDOR_MAX_PASS}x weight x {Config.MDL_CROSSING_SLOWDOWN}, "
                f"{len(CORRIDORS)} water-deer corridor gates",
                (f"vaccination: {vacc_note}" if vacc_note else
                 f"vaccination: {'ON (' + Config.VACC_MODE + ', ' + Config.VACC_TARGET_REGION + ')' if Config.ENABLE_VACCINATION else 'OFF'}"),
                "calendar columns are a NOMINAL label anchored on "
                f"SIM_START_YEAR={Config.SIM_START_YEAR}-{Config.SIM_START_MONTH:02d}, "
                "not a dated forecast",
                "sampling=centroid marks counties smaller than one model sub-cell",
                f"ROK/DPRK split: {ROK_SPLIT_INFO.get('method')} (v9.10); ROK exclave policy "
                f"'{ROK_EXCLAVE_INFO.get('policy')}' ({', '.join(ROK_EXCLAVE_INFO.get('sigungu', [])) or 'none'})",
            ]
            stem = f"{Config.SK_COUNTY_CSV_STEM}_{tag}"
            written = [write_county_summary_csv(
                rows, Config.OUT_DIR / f"{stem}_summary.csv", meta_lines=meta, logger=logger)]
            if Config.SK_COUNTY_WRITE_MONTHLY_WIDE:
                written.append(write_county_monthly_wide_csv(
                    rows, series, Config.OUT_DIR / f"{stem}_monthly_wide.csv",
                    month_offset=month_offset, start_year=Config.SIM_START_YEAR,
                    start_month=Config.SIM_START_MONTH, meta_lines=meta, logger=logger))
            if Config.SK_COUNTY_WRITE_LONG:
                written.append(write_county_long_csv(
                    rows, series, Config.OUT_DIR / f"{stem}_long.csv",
                    month_offset=month_offset, start_year=Config.SIM_START_YEAR,
                    start_month=Config.SIM_START_MONTH, meta_lines=meta, logger=logger))
            if Config.SK_COUNTY_MAP:
                mp = render_sk_county_map(rows, tag, title, thr)
                if mp is not None:
                    written.append(mp)
            if Config.DL_DIR:
                try:
                    Config.DL_DIR.mkdir(parents=True, exist_ok=True)
                    for f in written:
                        shutil.copy2(f, Config.DL_DIR / Path(f).name)
                    logger.info(f"   Also copied {len(written)} county file(s) to: {Config.DL_DIR}")
                except Exception as e:
                    logger.warning(f"  ⚠ Could not copy county files to DL_DIR: {e}")
            sk_county_outputs[tag] = written
            return written
        except Exception as e:
            logger.warning(f"  ⚠ county CSV export skipped for {tag} "
                           f"({type(e).__name__}: {e})", exc_info=True)
            return []
    def render_sk_county_map(rows, tag, title, thr):
        """Two-panel county choropleth: arrival month and final mean attack.
        This is what makes the CSV auditable — if a row says month 84, the map
        must show that county lighting up at month 84."""
        try:
            arr_map = zones_to_value_raster(rows, SK_ZONE_RASTER, "arrival_month")
            att_map = zones_to_value_raster(rows, SK_ZONE_RASTER, "final_mean_attack_pct")
            occ = SK_ZONE_RASTER > 0
            if not occ.any():
                return None
            rr, cc = np.where(occ)
            pad = 3
            r0, r1 = max(0, rr.min() - pad), min(NR, rr.max() + 1 + pad)
            c0, c1 = max(0, cc.min() - pad), min(NC, cc.max() + 1 + pad)
            sl = (slice(r0, r1), slice(c0, c1))
            ext = (0, c1 - c0, r1 - r0, 0)
            fig, ax = plt.subplots(1, 2, figsize=(17, 10), facecolor="white")
            base = np.ma.masked_where(~occ[sl], np.ones((r1 - r0, c1 - c0), float))
            band_cmap = ListedColormap(Config.SPREAD_BAND_COLORS_SK)
            band_norm = BoundaryNorm(list(Config.SPREAD_PCT_BANDS),
                                     ncolors=len(Config.SPREAD_PCT_BANDS) - 1, clip=True)
            for a in ax:
                a.set_facecolor("white")
                a.imshow(base, cmap=ListedColormap(["#f2f2f2"]), extent=ext,
                         origin="upper", aspect="equal", zorder=1)
                a.set_xticks([]); a.set_yticks([])
            am = np.ma.masked_invalid(arr_map[sl])
            im0 = ax[0].imshow(am, cmap="viridis_r", extent=ext, origin="upper",
                               aspect="equal", zorder=2)
            ax[0].set_title(f"Month PPR arrives in each county\n"
                            f"(county mean cumulative attack ≥ {100*thr:.0f}%; grey = not reached)",
                            fontsize=12, fontweight="bold")
            cb0 = fig.colorbar(im0, ax=ax[0], shrink=0.65, pad=0.02)
            cb0.set_label("simulation month", fontsize=10, fontweight="bold")
            ax[1].imshow(np.ma.masked_invalid(att_map[sl]), cmap=band_cmap, norm=band_norm,
                         extent=ext, origin="upper", aspect="equal", zorder=2)
            ax[1].set_title("Final cumulative attack rate by county (%)",
                            fontsize=12, fontweight="bold")
            sm = plt.cm.ScalarMappable(cmap=band_cmap, norm=band_norm); sm.set_array([])
            cb1 = fig.colorbar(sm, ax=ax[1], boundaries=list(Config.SPREAD_PCT_BANDS),
                               ticks=list(Config.SPREAD_PCT_BANDS), shrink=0.65, pad=0.02)
            cb1.set_label("% of local population ever infected", fontsize=10, fontweight="bold")
            for a in ax:
                mdl_c = MDL_BOUNDARY[sl]
                if mdl_c.any():
                    a.contour(mdl_c.astype(float), levels=[0.5], colors="#d81b60",
                              linewidths=1.4, zorder=5)
                if Config.SHOW_CORRIDOR_MARKERS and CORRIDORS:
                    plot_corridor_markers(a, CORRIDORS, row_offset=r0, col_offset=c0,
                                          label=False, size_scale=20.0)
            n_reach = sum(1 for r in rows if r.get("arrival_month"))
            fig.suptitle(f"PPR by South-Korean county (sigungu) — {title}\n"
                         f"{n_reach} of {len(rows)} counties reached; "
                         f"source table: {Config.SK_COUNTY_CSV_STEM}_{tag}_summary.csv",
                         fontsize=13, fontweight="bold")
            fig.tight_layout(rect=(0, 0, 1, 0.94))
            fp = Config.OUT_DIR / f"PPR_v9.3_SK_county_arrival_{tag}.png"
            fig.savefig(fp, dpi=140, facecolor="white", bbox_inches="tight")
            plt.close(fig)
            logger.info(f"   ✓ county choropleth: {fp}")
            return fp
        except Exception as e:
            logger.warning(f"   county choropleth skipped: {e}")
            return None
    def render_mode(mode, monthly_nk, monthly_sk, monthly_cum_nk, monthly_cum_sk):
        logger.info(f"\n[RENDERING VISUALIZATION — mode={mode}]")
        all_nk = np.concatenate([np.nan_to_num(e, nan=0)[NK].ravel() for e in monthly_nk]) if NK.any() else np.array([0.])
        all_sk = np.concatenate([np.nan_to_num(e, nan=0)[SK].ravel() for e in monthly_sk]) if SK.any() else np.array([0.])
        NK_VMAX = float(np.nanpercentile(all_nk, 99)) if len(all_nk) else 1.
        SK_VMAX = float(np.nanpercentile(all_sk, 99)) if len(all_sk) else 1.
        logger.info(f"  NK_VMAX (99th pct)={NK_VMAX:.4g}  SK_VMAX (99th pct)={SK_VMAX:.4g}")
        _v910_pub5 = bool(getattr(Config, "PUB_FIGS", False))      # 11 Nov: 40-in canvas -> larger text
        _v910_pf5 = 3.0 if _v910_pub5 else 1.0
        fig = plt.figure(figsize=(40, 36), dpi=150, facecolor="white")
        gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.25, left=0.08, right=0.92, top=0.96, bottom=0.05)
        def plot_map(ax, m, nk_data, sk_data, cum_nk_data, cum_sk_data):
            ax.set_facecolor(Config.MAP_BACKGROUND)
            # FIX (v9.1): _floor was only defined inside the SINGLE_MODE branch,
            # so the NK/SK branch raised NameError. Define it once, up front.
            _floor = Config.BAND_MIN_PCT_VISIBLE
            # ── Regional land basemap (neighbouring countries) ───────────────
            if Config.SHOW_REGIONAL_BASEMAP and LAND_BASE.any():
                base = np.ma.masked_where(~LAND_BASE, np.ones(SH, np.float32))
                ax.imshow(base, cmap=ListedColormap([Config.BASEMAP_LAND_COLOR]),
                          vmin=0, vmax=1, extent=(0, NC, NR, 0), origin="upper",
                          aspect="auto", alpha=Config.BASEMAP_LAND_ALPHA,
                          interpolation="nearest", zorder=0.5)
                if basemap_segments:
                    ax.add_collection(LineCollection(
                        basemap_segments, colors=Config.BASEMAP_LINE_COLOR,
                        linewidths=Config.BASEMAP_LINE_WIDTH, alpha=Config.BASEMAP_LINE_ALPHA,
                        zorder=0.7))
            # ── PPR MaxEnt suitability as the base map (continental mode) ────
            if base_is_maxent and (MAXENT is not None) and Config.SHOW_MAXENT_UNDERLAY:
                _mx_cmap = LinearSegmentedColormap.from_list("maxent_u", Config.MAXENT_CMAP, N=256)
                ax.imshow(np.ma.masked_invalid(MAXENT), cmap=_mx_cmap, vmin=0, vmax=1,
                          extent=(0, NC, NR, 0), origin="upper", aspect="auto",
                          alpha=Config.MAXENT_UNDERLAY_ALPHA, interpolation="bilinear",
                          zorder=1.6)
            land = np.ma.masked_invalid(np.where(KP, 0.15, np.nan))
            ax.imshow(land, cmap=mc(plt.cm.Greys), vmin=0, vmax=1, extent=(0, NC, NR, 0),
                      origin="upper", aspect="auto", alpha=Config.LAND_TINT_ALPHA,
                      interpolation="nearest", zorder=1)
            # ── RISK UNDERLAY (habitat/variable-driven establishment risk) ───
            if Config.SHOW_RISK_UNDERLAY and Config.SHOW_RISK_UNDERLAY_NK and not SINGLE_MODE:
                ax.imshow(RISK_NK_LAYER, cmap=CMAP_RISK_NK, vmin=0, vmax=1, extent=(0, NC, NR, 0),
                          origin="upper", aspect="auto", alpha=Config.RISK_UNDERLAY_ALPHA,
                          interpolation="bilinear", zorder=2.5)
            if Config.SHOW_RISK_UNDERLAY and Config.SHOW_RISK_UNDERLAY_SK:
                ax.imshow(RISK_SK_LAYER, cmap=CMAP_RISK_SK, vmin=0, vmax=1, extent=(0, NC, NR, 0),
                          origin="upper", aspect="auto", alpha=Config.RISK_UNDERLAY_ALPHA,
                          interpolation="bilinear", zorder=2.55)
            if Config.SHOW_RIVER_RASTER:
                river_plot = np.ma.masked_invalid(np.where(KP, RIV_nn, np.nan))
                cmap_river = LinearSegmentedColormap.from_list("river", ["#ffffff00", "#dbeafe", "#93c5fd"], N=256)
                ax.imshow(river_plot, cmap=mc(cmap_river), vmin=0, vmax=1, extent=(0, NC, NR, 0),
                          origin="upper", aspect="auto", alpha=0.20, interpolation="bilinear", zorder=3)
            if Config.SHOW_WATER_FILL:
                water_plot = np.ma.masked_invalid(np.where(water_raster > 0, 1.0, np.nan))
                ax.imshow(water_plot, cmap=LinearSegmentedColormap.from_list("w", ["#ffffff00", "#9ecae1"]),
                          vmin=0, vmax=1, extent=(0, NC, NR, 0), origin="upper", aspect="auto",
                          alpha=0.25, interpolation="nearest", zorder=3.3)
            _water_z = 6.0 if Config.BOUNDARIES_OVER_RISK else 3.35
            if water_segments:
                ax.add_collection(LineCollection(
                    water_segments, colors=Config.WATER_LINE_COLOR,
                    linewidths=Config.WATER_LINE_WIDTH, alpha=Config.WATER_LINE_ALPHA,
                    zorder=_water_z))
            if Config.SHOW_WETLAND_FILL:
                wetland_plot = np.ma.masked_invalid(np.where(wetland_raster > 0, 1.0, np.nan))
                ax.imshow(wetland_plot, cmap=LinearSegmentedColormap.from_list("wet", ["#ffffff00", "#a7e8e8"]),
                          vmin=0, vmax=1, extent=(0, NC, NR, 0), origin="upper", aspect="auto",
                          alpha=0.25, interpolation="nearest", zorder=3.4)
            _sig_z = 6.1 if Config.BOUNDARIES_OVER_RISK else 3.5
            if sig_segments:
                ax.add_collection(LineCollection(
                    sig_segments, colors=Config.SIG_LINE_COLOR,
                    linewidths=Config.SIG_LINE_WIDTH, alpha=Config.SIG_LINE_ALPHA,
                    zorder=_sig_z))
            if ROAD_MASK.any():
                road_plot = np.ma.masked_invalid(np.where(ROAD_MASK & KP, 1.0, np.nan))
                ax.imshow(road_plot, cmap=LinearSegmentedColormap.from_list("rd", ["#ffffff00", "#555555"]),
                          vmin=0, vmax=1, extent=(0, NC, NR, 0), origin="upper", aspect="auto",
                          alpha=0.65, interpolation="nearest", zorder=3.6)
            mdl_plot = np.ma.masked_invalid(np.where(MDL_BOUNDARY, 1.0, np.nan))
            ax.imshow(mdl_plot, cmap='Reds', vmin=0, vmax=1, extent=(0, NC, NR, 0),
                      origin="upper", aspect="auto", alpha=0.88, interpolation="nearest", zorder=8)
            if Config.SHOW_STUDY_AREA_OUTLINE and study_segments:
                ax.add_collection(LineCollection(
                    study_segments, colors=Config.STUDY_AREA_LINE_COLOR,
                    linewidths=Config.STUDY_AREA_LINE_WIDTH, alpha=Config.STUDY_AREA_LINE_ALPHA,
                    zorder=8.5))
            if Config.SHOW_OUTBREAK_SEED_MARKERS and OUTBREAK_ENTRIES:
                orr = [e[0] for e in OUTBREAK_ENTRIES]
                occ = [e[1] for e in OUTBREAK_ENTRIES]
                ax.scatter(np.array(occ) + 0.5, np.array(orr) + 0.5,
                           s=Config.OUTBREAK_MARKER_SIZE, marker="*",
                           c=Config.OUTBREAK_MARKER_COLOR, edgecolors="white",
                           linewidths=0.6, zorder=9)
            # v9.1: the corridor gates themselves (labels off here — at
            # continental scale they collide).
            if Config.SHOW_CORRIDOR_MARKERS and CORRIDORS:
                plot_corridor_markers(ax, CORRIDORS, label=False, size_scale=18.0)
            if Config.SHOW_PERCENT_BANDS and SINGLE_MODE:
                cum_all = np.where(NK, np.nan_to_num(cum_nk_data, nan=0.0),
                                   np.where(SK, np.nan_to_num(cum_sk_data, nan=0.0), np.nan)) * 100.0
                cum_all = np.where(KP & (cum_all >= _floor), cum_all, np.nan)
                ax.imshow(np.ma.masked_invalid(cum_all), cmap=mc(CMAP_BAND_ALL),
                          norm=BAND_NORM, extent=(0, NC, NR, 0), origin="upper",
                          aspect="auto", alpha=Config.BAND_ALPHA, interpolation="nearest", zorder=4.5)
            elif Config.SHOW_PERCENT_BANDS:
                nk_pct = np.nan_to_num(cum_nk_data, nan=0.0) * 100.0
                nk_pct = np.where(NK & (nk_pct >= _floor), nk_pct, np.nan)
                ax.imshow(np.ma.masked_invalid(nk_pct), cmap=mc(CMAP_BAND_NK),
                          norm=BAND_NORM, extent=(0, NC, NR, 0), origin="upper",
                          aspect="auto", alpha=Config.BAND_ALPHA, interpolation="nearest", zorder=4)
                sk_pct = np.nan_to_num(cum_sk_data, nan=0.0) * 100.0
                sk_pct = np.where(SK & (sk_pct >= Config.SK_BAND_MIN_PCT_VISIBLE), sk_pct, np.nan)
                ax.imshow(np.ma.masked_invalid(sk_pct), cmap=mc(CMAP_BAND_SK),
                          norm=BAND_NORM, extent=(0, NC, NR, 0), origin="upper",
                          aspect="auto", alpha=Config.BAND_ALPHA, interpolation="nearest", zorder=5)
            else:
                _nk = np.nan_to_num(nk_data, nan=0)
                _nk_in = _nk[NK] if NK.any() else np.array([0.])
                _nk_max = max(float(np.nanpercentile(_nk_in, 98)) if _nk_in.size else NK_VMAX, 1e-7)
                ax.imshow(np.ma.masked_invalid(np.where(NK, _nk, np.nan)), cmap=mc(CMAP_NK),
                          norm=PowerNorm(0.40, 0, _nk_max), extent=(0, NC, NR, 0),
                          origin="upper", aspect="auto", alpha=0.88, interpolation="bilinear", zorder=4)
                _sk = np.nan_to_num(sk_data, nan=0)
                _sk_valid = _sk[SK & (_sk > 0)] if SK.any() else np.array([])
                _sk_max = max(float(np.nanpercentile(_sk_valid, 98)) if len(_sk_valid) else SK_VMAX, 1e-9)
                ax.imshow(np.ma.masked_invalid(np.where(SK, _sk, np.nan)), cmap=mc(CMAP_SK),
                          norm=PowerNorm(0.35, 0, _sk_max), extent=(0, NC, NR, 0),
                          origin="upper", aspect="auto", alpha=0.87, interpolation="bilinear", zorder=5)
            # ── high-density gorani spreader dots ────────────────────────────
            if Config.SHOW_GORANI_DOTS and (not base_is_maxent) and _gor_hi_rows.size:
                dr, dc = _gor_hi_rows, _gor_hi_cols
                rank = _dens[dr, dc]
                if Config.GORANI_DOT_REQUIRE_SPREAD:
                    cum_all = np.where(NK, np.nan_to_num(cum_nk_data, nan=0.0),
                                       np.nan_to_num(cum_sk_data, nan=0.0))
                    sel = cum_all[dr, dc] > Config.GORANI_DOT_SPREAD_FLOOR
                    dr, dc, rank = dr[sel], dc[sel], rank[sel]
                if dr.size > Config.GORANI_DOT_MAX:
                    order = np.argsort(rank)[::-1][:Config.GORANI_DOT_MAX]
                    dr, dc = dr[order], dc[order]
                if dr.size:
                    ax.scatter(dc + 0.5, dr + 0.5, s=Config.GORANI_DOT_SIZE,
                               c=Config.GORANI_DOT_COLOR, marker="o",
                               edgecolors="#5a4b00", linewidths=0.25,
                               alpha=0.95, zorder=7)
            ax.set_xlim(0, NC)
            ax.set_ylim(NR, 0)
            ax.set_xticks([])
            ax.set_yticks([])
            mm = m % 12               # month-of-year (0-11) for seasonal params
            yr = m // 12 + 1          # simulation year (1-based)
            season = season_of(mm)
            for sp in ax.spines.values():
                sp.set_edgecolor(Config.SEASON_COLORS[season])
                sp.set_linewidth(4)
            beta_m = Config.MONTHLY_BETA[mm]
            t_v = T_VIRAL[mm]
            j_m = JOINT[mm]
            rc = "#b71c1c" if j_m > 0.70 else "#e65100" if j_m > 0.45 else "#1565c0"
            if _v910_pub5:      # 15 Nov: black panel titles (the J colour code has no key once line 2 is dropped)
                rc = "black"
            _yr_txt = f"  ·  Year {yr}" if (12 * max(1, int(Config.SIM_YEARS)) > 12) else ""
            ax.set_title(
                f"{Config.MONTH_NAMES[mm]}{_yr_txt}\n"
                f"Deer r={Config.DEER_RADIUS_KM[mm]:.1f}km  β={0.55*beta_m:.3f}  T={t_v:.3f}  J={j_m:.3f}",
                fontsize=11 * _v910_pf5, fontweight="bold", pad=8 * _v910_pf5, color=rc)
        _NM = len(monthly_nk)
        if _NM <= 12:
            _snaps = list(range(_NM))
        else:
            _snaps = sorted(set(np.linspace(0, _NM - 1, 12).round().astype(int).tolist()))
        for idx, m in enumerate(_snaps):
            row, col = divmod(idx, 4)
            ax = fig.add_subplot(gs[row, col])
            plot_map(ax, m, monthly_nk[m], monthly_sk[m],
                     monthly_cum_nk[m], monthly_cum_sk[m])
        if Config.SHOW_PERCENT_BANDS and SINGLE_MODE:
            cax = fig.add_axes([0.94, 0.55, 0.012, 0.38])
            sm = plt.cm.ScalarMappable(cmap=CMAP_BAND_ALL, norm=BAND_NORM); sm.set_array([])
            cbar = fig.colorbar(sm, cax=cax, boundaries=PCT_BANDS, ticks=PCT_BANDS, spacing="proportional")
            cbar.set_label("% PPR spread (cumulative attack rate)", fontsize=9, fontweight="bold")
            cbar.ax.set_yticklabels([f"{int(b)}%" for b in PCT_BANDS])
            if _v910_pub5:
                cbar.ax.tick_params(labelsize=10 * _v910_pf5, width=2, length=8)
                cbar.ax.yaxis.label.set_fontsize(9 * _v910_pf5)
        elif Config.SHOW_PERCENT_BANDS:
            cax_nk = fig.add_axes([0.935, 0.55, 0.010, 0.38])
            sm_nk = plt.cm.ScalarMappable(cmap=CMAP_BAND_NK, norm=BAND_NORM)
            sm_nk.set_array([])
            cbar_nk = fig.colorbar(sm_nk, cax=cax_nk, boundaries=PCT_BANDS, ticks=PCT_BANDS, spacing="proportional")
            cbar_nk.set_label("NK: % PPR spread (cumulative attack rate)", fontsize=9, fontweight="bold")
            cbar_nk.ax.set_yticklabels([f"{int(b)}%" for b in PCT_BANDS])
            cax_sk = fig.add_axes([0.955, 0.55, 0.010, 0.38])
            sm_sk = plt.cm.ScalarMappable(cmap=CMAP_BAND_SK, norm=BAND_NORM)
            sm_sk.set_array([])
            cbar_sk = fig.colorbar(sm_sk, cax=cax_sk, boundaries=PCT_BANDS, ticks=PCT_BANDS, spacing="proportional")
            cbar_sk.set_label("SK: % PPR spread (via the MDL corridor gates)", fontsize=9, fontweight="bold")
            cbar_sk.ax.set_yticklabels([f"{int(b)}%" for b in PCT_BANDS])
        else:
            cax_nk = fig.add_axes([0.935, 0.55, 0.008, 0.38])
            sm_nk = plt.cm.ScalarMappable(cmap=mc(CMAP_NK), norm=PowerNorm(0.40, 0, NK_VMAX))
            sm_nk.set_array([])
            cbar_nk = fig.colorbar(sm_nk, cax=cax_nk)
            cbar_nk.set_label("NK Infected Intensity", fontsize=9, fontweight="bold")
            cax_sk = fig.add_axes([0.945, 0.55, 0.008, 0.38])
            sm_sk = plt.cm.ScalarMappable(cmap=mc(CMAP_SK), norm=PowerNorm(0.35, 0, SK_VMAX))
            sm_sk.set_array([])
            cbar_sk = fig.colorbar(sm_sk, cax=cax_sk)
            cbar_sk.set_label("SK Infected Intensity", fontsize=9, fontweight="bold")
        leg_seasons = [Patch(fc=v, ec="#777", lw=1.5, label=k.replace("\n", " "))
                       for k, v in Config.SEASON_COLORS.items()]
        leg_geo = [Line2D([0], [0], color=Config.WATER_LINE_COLOR, lw=1.4,
                          label=f'Water bodies / water.shp ({len(WATER_POLYS)})'),
                   Line2D([0], [0], color=Config.SIG_LINE_COLOR, lw=1.0,
                          label=f'Counties / sig.shp ({len(SIG_POLYS)})')]
        if CORRIDORS and Config.SHOW_CORRIDOR_MARKERS:
            leg_geo.append(Line2D([0], [0], marker='v', linestyle='none',
                                  markerfacecolor="#00e5ff", markeredgecolor="#003b46",
                                  markersize=9,
                                  label=f'DMZ water-deer corridor gates ({len(CORRIDORS)}; size ∝ weight)'))
        if Config.SHOW_RISK_UNDERLAY and Config.SHOW_RISK_UNDERLAY_NK and not (_v910_pub5 and SINGLE_MODE):
            leg_geo.append(Patch(facecolor="#cf9b63", edgecolor='#8a6033',
                                 label='NK habitat risk (gorani TIF + covariates)'))
        if Config.SHOW_RISK_UNDERLAY and Config.SHOW_RISK_UNDERLAY_SK:
            leg_geo.append(Patch(facecolor="#8fb99b", edgecolor='#5f8a6b',
                                 label='SK habitat risk (gorani TIF + covariates)'))
        if ROAD_MASK.any():
            leg_geo.append(Patch(facecolor='#555555', edgecolor='#222222',
                                 label=f'Roads (barrier δ={Config.HAN_ROAD_REDUCTION})'))
        if Config.SHOW_GORANI_DOTS and not (_v910_pub5 and (base_is_maxent or not _gor_hi_rows.size)):
            leg_geo.append(Line2D([0], [0], marker='o', linestyle='none',
                                  markerfacecolor=Config.GORANI_DOT_COLOR,
                                  markeredgecolor='#5a4b00', markersize=7,
                                  label=f'High-density gorani hotspots (q{Config.GORANI_DOT_QUANTILE:.2f})'))
        if base_is_maxent and (MAXENT is not None) and Config.SHOW_MAXENT_UNDERLAY:
            leg_geo.append(Patch(facecolor="#fdae61", edgecolor="#b03a2e",
                                 label='PPR MaxEnt suitability (base map)'))
        if Config.SHOW_REGIONAL_BASEMAP and LAND_BASE.any():
            leg_geo.append(Patch(facecolor=Config.BASEMAP_LAND_COLOR, edgecolor=Config.BASEMAP_LINE_COLOR,
                                 label='Neighbouring land (regional basemap)'))
        if Config.SHOW_STUDY_AREA_OUTLINE and study_segments:
            leg_geo.append(Line2D([0], [0], color=Config.STUDY_AREA_LINE_COLOR,
                                  lw=Config.STUDY_AREA_LINE_WIDTH, label='Study area (PPR_studyarea)'))
        if Config.SHOW_OUTBREAK_SEED_MARKERS and OUTBREAK_ENTRIES:
            leg_geo.append(Line2D([0], [0], marker='*', linestyle='none',
                                  markerfacecolor=Config.OUTBREAK_MARKER_COLOR,
                                  markeredgecolor='white', markersize=11,
                                  label=f'Outbreak seed points ({len(OUTBREAK_ENTRIES)})'))
        fig.legend(handles=leg_seasons + leg_geo, ncol=4 if _v910_pub5 else 6,
                   fontsize=8 * (3.0 if _v910_pub5 else 1.0), frameon=True,
                   **(dict(loc="upper center", bbox_to_anchor=(0.5, 0.044)) if _v910_pub5
                      else dict(loc="lower center", bbox_to_anchor=(0.5, 0.005))))
        if CORRIDORS:
            fence_note = (f"MDL GATED: {len(CORRIDORS)} corridors, pass≤"
                          f"{Config.CORRIDOR_MAX_PASS*Config.MDL_CROSSING_SLOWDOWN:.2f}×weight, "
                          f"line sealed at {Config.MDL_SEALED_FACTOR}")
        else:
            fence_note = (f"DMZ fence φ={Config.HAN_FENCE_REDUCTION} (uniform, leaky)"
                          if Config.HAN_APPLY_DMZ_FENCE else "DMZ fence: off")
        color_note = ("Colour = % infected (cumulative attack rate 1-S, 0-100% bands) | "
                      if Config.SHOW_PERCENT_BANDS else "Colour = infection intensity | ")
        vacc_note = (f" | Vaccination: {Config.VACC_MODE}, {Config.VACC_TARGET_REGION}, "
                     f"start m{Config.VACC_START_MONTH + 1}" if Config.ENABLE_VACCINATION else "")
        fig.suptitle(
            f"PPR DYNAMIC RAD EPIDEMIC (v9.3) — persistent S-I-R, corridor-gated cross-MDL spread\n"
            f"gamma={Config.GAMMA_RECOVERY} | {Config.SIM_YEARS}-year run | "
            f"Landscape calib: Han et al. 2022 [{fence_note}, roads δ={Config.HAN_ROAD_REDUCTION}] | "
            f"Water: {len(WATER_POLYS)} | Wetland: {len(WETLAND_POLYS)} | Counties: {len(SIG_POLYS)}\n"
            f"{color_note}Warm=NK source | Blue-green=SK spread through the water-deer corridors | "
            "▼ = corridor gate (size ∝ ecological weight) | "
            "Yellow dots = high-density gorani hotspots | S-I-R persists month-to-month" + vacc_note,
            fontsize=13, fontweight="bold", y=0.985)
        fp = Config.OUT_DIR / "PPR_v9.3_PERCENT_SPREAD_CORRIDOR_GATED.png"
        logger.info("  Saving figure...")
        plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        logger.info(f"  ✓ FIGURE SAVED: {fp}")
        if Config.DL_DIR:
            try:
                Config.DL_DIR.mkdir(parents=True, exist_ok=True)
                shutil.copy2(fp, Config.DL_DIR / fp.name)
                logger.info(f"   Also copied to: {Config.DL_DIR}")
            except Exception as e:
                logger.warning(f"  ⚠ Could not copy to DL_DIR: {e}")
        # ── STYLED animated spread GIF ───────────────────────────────────────
        try:
            import io as _io
            from PIL import Image
            FW, FH, FDPI = 7.6, 8.8, 110
            frames = []
            for m in range(len(monthly_nk)):
                figF = plt.figure(figsize=(FW, FH), dpi=FDPI, facecolor="white")
                axF = figF.add_axes([0.015, 0.02, 0.80, 0.95])
                plot_map(axF, m, monthly_nk[m], monthly_sk[m],
                         monthly_cum_nk[m], monthly_cum_sk[m])
                if Config.SHOW_PERCENT_BANDS and SINGLE_MODE:
                    caA = figF.add_axes([0.86, 0.55, 0.024, 0.34])
                    smA = plt.cm.ScalarMappable(cmap=CMAP_BAND_ALL, norm=BAND_NORM); smA.set_array([])
                    cbA = figF.colorbar(smA, cax=caA, boundaries=PCT_BANDS, ticks=PCT_BANDS, spacing="proportional")
                    cbA.ax.set_yticklabels([f"{int(b)}" for b in PCT_BANDS], fontsize=6)
                    cbA.set_label("% spread", fontsize=7)
                elif Config.SHOW_PERCENT_BANDS:
                    caN = figF.add_axes([0.83, 0.55, 0.022, 0.34])
                    smN = plt.cm.ScalarMappable(cmap=CMAP_BAND_NK, norm=BAND_NORM); smN.set_array([])
                    cbN = figF.colorbar(smN, cax=caN, boundaries=PCT_BANDS, ticks=PCT_BANDS, spacing="proportional")
                    cbN.ax.set_yticklabels([f"{int(b)}" for b in PCT_BANDS], fontsize=6)
                    cbN.set_label("NK %", fontsize=7)
                    caS = figF.add_axes([0.90, 0.55, 0.022, 0.34])
                    smS = plt.cm.ScalarMappable(cmap=CMAP_BAND_SK, norm=BAND_NORM); smS.set_array([])
                    cbS = figF.colorbar(smS, cax=caS, boundaries=PCT_BANDS, ticks=PCT_BANDS, spacing="proportional")
                    cbS.ax.set_yticklabels([f"{int(b)}" for b in PCT_BANDS], fontsize=6)
                    cbS.set_label("SK %", fontsize=7)
                figF.text(0.5, 0.005,
                          f"PPR dynamic spread — {Config.MONTH_NAMES[m % 12]} Year {m // 12 + 1} "
                          f"(month {m + 1}/{len(monthly_nk)}) — cumulative % infected, S-I-R persists",
                          ha="center", fontsize=8, color="#333")
                buf = _io.BytesIO()
                figF.savefig(buf, format="png", facecolor="white")  # NO tight bbox -> uniform size
                plt.close(figF)
                buf.seek(0)
                frames.append(Image.open(buf).convert("RGB"))
            if frames:
                fp_gif = Config.OUT_DIR / "PPR_v9.3_styled_spread.gif"
                frames[0].save(fp_gif, save_all=True, append_images=frames[1:],
                               duration=650, loop=0, optimize=True)
                logger.info(f"  ✓ Styled spread GIF: {fp_gif}")
                if Config.DL_DIR:
                    try:
                        shutil.copy2(fp_gif, Config.DL_DIR / fp_gif.name)
                    except Exception:
                        pass
        except Exception as e:
            logger.warning(f"  ⚠ Styled spread GIF skipped: {e}")
        return fp
    output_files = {}
    for mode, (m_nk, m_sk, c_nk, c_sk) in results_by_mode.items():
        output_files[mode] = render_mode(mode, m_nk, m_sk, c_nk, c_sk)
    # ══ v9.3: COUNTY CSV FROM THE CONTINENTAL RUN ════════════════════════════
    # The Korea leg below is the authoritative NK->SK narrative, but the
    # continental run also covers South Korea, so it is tabulated too when
    # SK_COUNTY_ALSO_MAIN_RUN is on. Both files share one month numbering.
    if Config.SK_COUNTY_ALSO_MAIN_RUN:
        _res_main = rad_results_holder.get("res")
        if _res_main is not None:
            if Config.CACHE_MAIN_RUN:
                save_epidemic_cache(_res_main["monthly_cum"], KP, "main_run", 0, 0, 0,
                                    sk_mask=SK_KR, monthly_I=_res_main.get("monthly_I"))
            export_sk_counties(
                _res_main["monthly_cum"], KP, tag="main_run",
                title="continental run (China source grid, whole domain)",
                row_offset=0, col_offset=0, month_offset=0,
                sk_mask=SK_KR, monthly_I=_res_main.get("monthly_I"))
    # ══ EXTENDED DYNAMIC-EPIDEMIC OUTPUTS ════════════════════════════════════
    def render_dynamic_outputs():
        res = rad_results_holder.get("res")
        if res is None:
            return
        try:
            ts = res["time_series"]
            months = np.arange(1, len(ts["n_infected_cells"]) + 1)
            # (1) Epidemic time-series (6 panels).
            fig, ax = plt.subplots(2, 3, figsize=(18, 9), facecolor="white")
            panels = [
                ("infected_area_km2", "Infected area (km²)", "#c0392b"),
                ("susceptible_area_km2", "Susceptible area (km²)", "#2980b9"),
                ("recovered_area_km2", "Recovered area (km²)", "#27ae60"),
                ("n_infected_cells", "Infected cells (count)", "#8e44ad"),
                ("growth_rate", "Spatial growth rate (Δcells/cell)", "#d35400"),
                ("front_km", "Epidemic front from border (km)", "#16a085"),
            ]
            for a, (key, title, col) in zip(ax.ravel(), panels):
                a.plot(months, ts[key], "-o", color=col, lw=2, ms=3)
                if key == "front_km" and "mean_front_km" in ts:
                    a.plot(months, ts["mean_front_km"], "--s", color="#7f8c8d", lw=1.4, ms=2, label="mean")
                    a.legend(fontsize=8)
                a.set_title(title, fontsize=11, fontweight="bold")
                a.set_xlabel("Month"); a.grid(alpha=0.3)
                if len(months) <= 24:
                    a.set_xticks(months)
            fig.suptitle("PPR dynamic epidemic — time-series (RAD engine, v9.3)",
                         fontsize=14, fontweight="bold")
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            fp_ts = Config.OUT_DIR / "PPR_v9.3_epidemic_timeseries.png"
            fig.savefig(fp_ts, dpi=130, facecolor="white"); plt.close(fig)
            logger.info(f"   time-series: {fp_ts}")
            # (2) Summary maps: cumulative max attack, time-to-arrival, final recovered.
            def _mask_show(ax_, arr, cmap, title, vmin=None, vmax=None):
                ax_.set_facecolor("white")
                ax_.imshow(np.ma.masked_invalid(np.where(KP, arr, np.nan)), cmap=cmap,
                           vmin=vmin, vmax=vmax, extent=(0, NC, NR, 0), origin="upper", aspect="auto")
                if MDL_BOUNDARY is not None:
                    ax_.imshow(np.ma.masked_where(~MDL_BOUNDARY, np.ones_like(arr)),
                               cmap=LinearSegmentedColormap.from_list("m", ["#7f0000", "#7f0000"]),
                               extent=(0, NC, NR, 0), origin="upper", aspect="auto", alpha=0.9)
                if Config.SHOW_CORRIDOR_MARKERS and CORRIDORS:
                    plot_corridor_markers(ax_, CORRIDORS, label=False, size_scale=16.0)
                ax_.set_title(title, fontsize=11, fontweight="bold")
                ax_.set_xticks([]); ax_.set_yticks([]); ax_.set_xlim(0, NC); ax_.set_ylim(NR, 0)
            fig, ax = plt.subplots(1, 3, figsize=(21, 8), facecolor="white")
            _mask_show(ax[0], res["cum_max"] * 100.0, "YlOrRd", "Maximum infection extent (% attack)", 0, 100)
            tta = res["time_to_arrival"]
            _tta_max = max(12, int(len(res.get("monthly_cum") or [12])))
            _mask_show(ax[1], np.where(np.isfinite(tta), tta + 1, np.nan), "viridis_r",
                       "Time-to-arrival (month)", 1, _tta_max)
            _mask_show(ax[2], res["final_R"] * 100.0, "Greens", "Final recovered (%)", 0, 100)
            fig.suptitle("PPR dynamic epidemic — summary maps (RAD engine, v9.3)",
                         fontsize=14, fontweight="bold")
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            fp_sm = Config.OUT_DIR / "PPR_v9.3_summary_maps.png"
            fig.savefig(fp_sm, dpi=130, facecolor="white"); plt.close(fig)
            logger.info(f"   summary maps: {fp_sm}")
            # (3) Annual end-of-year snapshots.
            cum_list = res.get("monthly_cum") or []
            n_years = max(1, int(round(len(cum_list) / 12.0)))
            if len(cum_list) >= 12 and n_years >= 1:
                year_idx = [min(len(cum_list) - 1, (y + 1) * 12 - 1) for y in range(n_years)]
                ncol = min(5, len(year_idx)); nrow = int(math.ceil(len(year_idx) / ncol))
                figY, axY = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 6.8 * nrow),
                                         facecolor="white", squeeze=False)
                band_cmap = ListedColormap(["#ffffcc", "#fed976", "#fd8d3c", "#e31a1c", "#800026"])
                band_norm = BoundaryNorm(list(Config.SPREAD_PCT_BANDS),
                                         ncolors=len(Config.SPREAD_PCT_BANDS) - 1, clip=True)
                for k, mi in enumerate(year_idx):
                    a = axY[k // ncol][k % ncol]
                    a.set_facecolor("white")
                    if base_is_maxent and MAXENT is not None:
                        a.imshow(np.ma.masked_invalid(MAXENT),
                                 cmap=LinearSegmentedColormap.from_list("mu", Config.MAXENT_CMAP, N=256),
                                 vmin=0, vmax=1, extent=(0, NC, NR, 0), origin="upper",
                                 aspect="auto", alpha=0.55, interpolation="bilinear", zorder=1)
                    pct = np.where(KP, np.nan_to_num(cum_list[mi], nan=0.0) * 100.0, np.nan)
                    pct = np.where(KP & (pct >= 1.0), pct, np.nan)
                    a.imshow(np.ma.masked_invalid(pct), cmap=band_cmap, norm=band_norm,
                             extent=(0, NC, NR, 0), origin="upper", aspect="auto",
                             alpha=0.72, interpolation="nearest", zorder=3)
                    if MDL_BOUNDARY is not None:
                        a.contour(MDL_BOUNDARY.astype(float), levels=[0.5], colors="#7f0000",
                                  linewidths=1.0, zorder=5)
                    if Config.SHOW_CORRIDOR_MARKERS and CORRIDORS:
                        plot_corridor_markers(a, CORRIDORS, label=False, size_scale=14.0)
                    if Config.SHOW_OUTBREAK_SEED_MARKERS and OUTBREAK_ENTRIES:
                        orr = [e[0] for e in OUTBREAK_ENTRIES]; occ = [e[1] for e in OUTBREAK_ENTRIES]
                        a.scatter(np.array(occ) + 0.5, np.array(orr) + 0.5, s=14, marker="*",
                                  c=Config.OUTBREAK_MARKER_COLOR, edgecolors="white",
                                  linewidths=0.3, zorder=6)
                    _mean = float(np.nanmean(np.nan_to_num(cum_list[mi])[KP])) * 100.0
                    if getattr(Config, "PUB_FIGS", False):
                        # 11 Nov: two lines, so that neighbouring titles no longer run into each other
                        a.title.set_text(f"End of Year {k + 1} (month {mi + 1})\n"
                                         f"mean attack rate {_mean:.1f}%")
                        a.title.set_fontsize(11); a.title.set_fontweight("bold")
                    else:
                        a.set_title(f"End of Year {k + 1}  (month {mi + 1})   mean attack={_mean:.1f}%",
                                    fontsize=11, fontweight="bold")
                    a.set_xticks([]); a.set_yticks([]); a.set_xlim(0, NC); a.set_ylim(NR, 0)
                for k in range(len(year_idx), nrow * ncol):
                    axY[k // ncol][k % ncol].axis("off")
                sm = plt.cm.ScalarMappable(cmap=band_cmap, norm=band_norm); sm.set_array([])
                cb = figY.colorbar(sm, ax=axY, boundaries=list(Config.SPREAD_PCT_BANDS),
                                   ticks=list(Config.SPREAD_PCT_BANDS), shrink=0.6, pad=0.02)
                cb.set_label("% PPR spread (cumulative attack rate)", fontsize=10, fontweight="bold")
                figY.suptitle(f"PPR spread — annual snapshots over {n_years} years (RAD engine, v9.3)",
                              fontsize=14, fontweight="bold")
                fp_y = Config.OUT_DIR / "PPR_v9.3_annual_snapshots.png"
                figY.savefig(fp_y, dpi=130, facecolor="white", bbox_inches="tight"); plt.close(figY)
                logger.info(f"   annual snapshots: {fp_y}")
                if Config.DL_DIR:
                    try:
                        shutil.copy2(fp_y, Config.DL_DIR / fp_y.name)
                    except Exception:
                        pass
        except Exception as e:
            logger.warning(f"   extended outputs skipped: {e}")
    render_dynamic_outputs()
    # ══ MDL GATE-PERMEABILITY MAP (v9.1) ═════════════════════════════════════
    # Shows the fence itself: which stretches are sealed and where the water-deer
    # corridors open it. This is the figure that makes the gating auditable.
    def render_mdl_gate_map():
        if not (Config.SHOW_MDL_GATE_MAP and CORRIDORS and MDL_FENCE_FIELD is not None):
            return
        try:
            focus = KOREA_MASK if KOREA_MASK.any() else KP
            rr, cc = np.where(MDL_BOUNDARY | focus)
            pad = 6
            if getattr(Config, "PUB_FIGS", False) and MDL_BOUNDARY.any():
                rr, cc = np.where(MDL_BOUNDARY)          # 11 Nov: zoom on the fence
                pad = 14
            r0, r1 = max(0, rr.min() - pad), min(NR, rr.max() + 1 + pad)
            c0, c1 = max(0, cc.min() - pad), min(NC, cc.max() + 1 + pad)
            sl = (slice(r0, r1), slice(c0, c1))
            ext = (0, c1 - c0, r1 - r0, 0)
            fig, ax = plt.subplots(figsize=((12, 7.5) if getattr(Config, "PUB_FIGS", False) else (15, 11)),
                                   dpi=150, facecolor="white")
            ax.set_facecolor("white")
            # land / habitat context
            suit = np.where(KP, np.nan_to_num(GOR_EFF), np.nan)[sl]
            ax.imshow(np.ma.masked_invalid(suit),
                      cmap=LinearSegmentedColormap.from_list("suit", Config.MAXENT_CMAP, N=256),
                      vmin=0, vmax=1, extent=ext, origin="upper", aspect="equal",
                      alpha=0.55, interpolation="bilinear", zorder=1)
            if sig_segments:
                ax.add_collection(LineCollection(
                    [seg - np.array([c0, r0]) for seg in sig_segments], colors="#6a6a6a",
                    linewidths=0.25, alpha=0.45, zorder=2))
            # the fence, coloured by pass-through
            ffc = np.where(MDL_BOUNDARY, MDL_FENCE_FIELD, np.nan)[sl]
            im = ax.imshow(np.ma.masked_invalid(ffc), cmap="plasma",
                           vmin=0, vmax=float(Config.CORRIDOR_MAX_PASS * Config.MDL_CROSSING_SLOWDOWN),
                           extent=ext, origin="upper", aspect="equal",
                           interpolation="nearest", zorder=4)
            plot_corridor_markers(ax, CORRIDORS, row_offset=r0, col_offset=c0,
                                  label=True, size_scale=40.0)
            ax.set_xticks([]); ax.set_yticks([])
            if getattr(Config, "PUB_FIGS", False):
                ax.set_xlim(0, c1 - c0); ax.set_ylim(r1 - r0, 0); ax.set_autoscale_on(False)
            cb = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
            cb.set_label("MDL pass-through per cell (fraction of transport that crosses)",
                         fontsize=10, fontweight="bold")
            _bar_km = 50.0; _bar_cells = _bar_km / (PW / 1000.0)
            _x0, _y0 = (c1 - c0) * 0.62, (r1 - r0) * 0.95
            ax.plot([_x0, _x0 + _bar_cells], [_y0, _y0], "k-", lw=3, zorder=9)
            ax.text(_x0 + _bar_cells / 2, _y0 - (r1 - r0) * 0.015, f"{int(_bar_km)} km",
                    ha="center", va="bottom", fontsize=10, fontweight="bold")
            ax.set_title(
                "MDL fence permeability — water-deer corridor gates\n"
                f"sealed line = {Config.MDL_SEALED_FACTOR:.2f}; gate peak = "
                f"{Config.CORRIDOR_MAX_PASS:.2f} × weight × {Config.MDL_CROSSING_SLOWDOWN:.2f}; "
                f"gate radius {Config.CORRIDOR_GATE_RADIUS_KM:.0f} km; "
                f"seasonal gating {'ON' if Config.CORRIDOR_SEASONAL_GATING else 'OFF'} "
                f"(Han φ={Config.HAN_FENCE_REDUCTION} is now the ceiling, not the default)",
                fontsize=12, fontweight="bold")
            fp_g = Config.OUT_DIR / "PPR_v9.3_MDL_gates.png"
            fig.savefig(fp_g, dpi=150, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            logger.info(f"   MDL gate-permeability map: {fp_g}")
            if Config.DL_DIR:
                try:
                    shutil.copy2(fp_g, Config.DL_DIR / fp_g.name)
                except Exception:
                    pass
        except Exception as e:
            logger.warning(f"   MDL gate map skipped: {e}")
    render_mdl_gate_map()
    # ══ HABITAT / ACCESSIBLE / BARRIER MAP (Han et al. 2022 Fig.2 analogue) ═══
    def render_habitat_classification():
        if not Config.SHOW_HABITAT_CLASSIFICATION:
            return
        try:
            from matplotlib.patches import Patch as _Patch
            fig, ax = plt.subplots(figsize=(12, 13), dpi=150, facecolor="white")
            ax.set_facecolor("white")
            cls_cmap = ListedColormap([Config.CLASS_COLOR_BARRIER,
                                       Config.CLASS_COLOR_ACCESSIBLE,
                                       Config.CLASS_COLOR_HABITAT])
            cls_norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], ncolors=3, clip=True)
            ax.imshow(np.ma.masked_invalid(HAB_CLASS), cmap=cls_cmap, norm=cls_norm,
                      extent=(0, NC, NR, 0), origin="upper", aspect="auto",
                      interpolation="nearest", zorder=2)
            if water_segments:
                ax.add_collection(LineCollection(
                    water_segments, colors=Config.WATER_LINE_COLOR,
                    linewidths=Config.WATER_LINE_WIDTH, alpha=Config.WATER_LINE_ALPHA, zorder=4))
            if sig_segments:
                ax.add_collection(LineCollection(
                    sig_segments, colors="#555555",
                    linewidths=Config.SIG_LINE_WIDTH, alpha=0.7, zorder=5))
            ax.imshow(np.ma.masked_where(~MDL_BOUNDARY, np.ones_like(HAB_CLASS)),
                      cmap=LinearSegmentedColormap.from_list("mdl", ["#ffffff", "#ffffff"]),
                      extent=(0, NC, NR, 0), origin="upper", aspect="auto", zorder=6)
            ax.contour(MDL_BOUNDARY.astype(float), levels=[0.5], colors="#d81b60",
                       linewidths=1.6, zorder=6.5)
            if Config.SHOW_CORRIDOR_MARKERS and CORRIDORS:
                plot_corridor_markers(ax, CORRIDORS, label=False, size_scale=22.0)
            if Config.SHOW_STUDY_AREA_OUTLINE and study_segments:
                ax.add_collection(LineCollection(
                    study_segments, colors=Config.STUDY_AREA_LINE_COLOR,
                    linewidths=Config.STUDY_AREA_LINE_WIDTH, alpha=Config.STUDY_AREA_LINE_ALPHA,
                    zorder=6.7))
            if Config.SHOW_OUTBREAK_SEED_MARKERS and OUTBREAK_ENTRIES:
                orr = [e[0] for e in OUTBREAK_ENTRIES]; occ = [e[1] for e in OUTBREAK_ENTRIES]
                ax.scatter(np.array(occ) + 0.5, np.array(orr) + 0.5,
                           s=Config.OUTBREAK_MARKER_SIZE, marker="*",
                           c=Config.OUTBREAK_MARKER_COLOR, edgecolors="white",
                           linewidths=0.6, zorder=7.5, label="Outbreak seed points")
            if getattr(Config, "CLASS_MAP_OVERLAY_REACHED", False):
                _res = rad_results_holder.get("res")
                if _res is not None and _res.get("monthly_cum"):
                    reached = np.nan_to_num(_res["monthly_cum"][-1]) >= Config.RAD_ARRIVAL_THRESHOLD
                    rr, cc = np.where(reached & KP)
                    if rr.size:
                        ax.scatter(cc + 0.5, rr + 0.5, s=3, c="#e53935", marker="s",
                                   linewidths=0, alpha=0.55, zorder=7,
                                   label="PPR-reached cells (final month)")
            ax.set_xlim(0, NC); ax.set_ylim(NR, 0)
            ax.set_xticks([]); ax.set_yticks([])
            _cell_km = PW / 1000.0
            _bar_km = 50.0
            _bar_cells = _bar_km / _cell_km
            _x0, _y0 = NC * 0.62, NR * 0.94
            ax.plot([_x0, _x0 + _bar_cells], [_y0, _y0], "k-", lw=3, zorder=8)
            ax.text(_x0 + _bar_cells / 2, _y0 - NR * 0.012, f"{int(_bar_km)} km",
                    ha="center", va="bottom", fontsize=10, fontweight="bold")
            legend = [_Patch(fc=Config.CLASS_COLOR_HABITAT, ec="#3a3a3a", label="Habitat cell"),
                      _Patch(fc=Config.CLASS_COLOR_ACCESSIBLE, ec="#3a3a3a", label="Accessible cell"),
                      _Patch(fc=Config.CLASS_COLOR_BARRIER, ec="#3a3a3a", label="Barrier cell"),
                      Line2D([0], [0], color="#d81b60", lw=1.6, label="MDL boundary")]
            if CORRIDORS:
                legend.append(Line2D([0], [0], marker="v", linestyle="none",
                                     markerfacecolor="#00e5ff", markeredgecolor="#003b46",
                                     markersize=9, label="Water-deer corridor gate"))
            ax.legend(handles=legend, loc="lower left", fontsize=11, frameon=True,
                      facecolor="white", edgecolor="#888")
            _nt = max(int(KP.sum()), 1)
            ax.set_title(
                "Gorani landscape classification — Habitat / Accessible / Barrier\n"
                f"(habitat GOR_EFF≥{_hab_thr:.2f}; accessible = cost-distance ≤ "
                f"{Config.ACCESS_COSTDIST_MAX_KM:.0f} km from habitat; Han et al. 2022 Fig.2 analogue)   "
                f"H={100*int(HABITAT_MASK.sum())/_nt:.0f}%  "
                f"A={100*int(ACCESS_MASK.sum())/_nt:.0f}%  "
                f"B={100*int(BARRIER_MASK.sum())/_nt:.0f}%",
                fontsize=12, fontweight="bold")
            fp_cls = Config.OUT_DIR / "PPR_v9.3_habitat_classification.png"
            fig.savefig(fp_cls, dpi=150, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            logger.info(f"   habitat classification map: {fp_cls}")
            if Config.DL_DIR:
                try:
                    shutil.copy2(fp_cls, Config.DL_DIR / fp_cls.name)
                except Exception:
                    pass
        except Exception as e:
            logger.warning(f"   habitat classification map skipped: {e}")
    render_habitat_classification()
    # ══ PPR MaxEnt SUITABILITY MAP ═══════════════════════════════════════════
    def render_maxent_map():
        if not Config.SHOW_MAXENT_MAP or MAXENT is None:
            return
        try:
            fig, ax = plt.subplots(figsize=(12, 13), dpi=150, facecolor="white")
            ax.set_facecolor("white")
            if Config.SHOW_REGIONAL_BASEMAP and LAND_BASE.any():
                ax.imshow(np.ma.masked_where(~LAND_BASE, np.ones(SH, np.float32)),
                          cmap=ListedColormap([Config.BASEMAP_LAND_COLOR]),
                          vmin=0, vmax=1, extent=(0, NC, NR, 0), origin="upper",
                          aspect="auto", alpha=Config.BASEMAP_LAND_ALPHA,
                          interpolation="nearest", zorder=0.5)
                if basemap_segments:
                    ax.add_collection(LineCollection(
                        basemap_segments, colors=Config.BASEMAP_LINE_COLOR,
                        linewidths=Config.BASEMAP_LINE_WIDTH, alpha=Config.BASEMAP_LINE_ALPHA, zorder=0.7))
            mx_cmap = LinearSegmentedColormap.from_list("maxent", Config.MAXENT_CMAP, N=256)
            im = ax.imshow(np.ma.masked_invalid(MAXENT), cmap=mx_cmap, vmin=0, vmax=1,
                           extent=(0, NC, NR, 0), origin="upper", aspect="auto",
                           interpolation="bilinear", zorder=2)
            if water_segments:
                ax.add_collection(LineCollection(
                    water_segments, colors=Config.WATER_LINE_COLOR,
                    linewidths=Config.WATER_LINE_WIDTH, alpha=Config.WATER_LINE_ALPHA, zorder=4))
            if sig_segments:
                ax.add_collection(LineCollection(
                    sig_segments, colors="#666666",
                    linewidths=Config.SIG_LINE_WIDTH, alpha=0.6, zorder=5))
            ax.contour(MDL_BOUNDARY.astype(float), levels=[0.5], colors="#d81b60",
                       linewidths=1.6, zorder=6)
            if Config.SHOW_STUDY_AREA_OUTLINE and study_segments:
                ax.add_collection(LineCollection(
                    study_segments, colors=Config.STUDY_AREA_LINE_COLOR,
                    linewidths=Config.STUDY_AREA_LINE_WIDTH, alpha=Config.STUDY_AREA_LINE_ALPHA, zorder=6.7))
            if Config.SHOW_CORRIDOR_MARKERS and CORRIDORS:
                plot_corridor_markers(ax, CORRIDORS, label=False, size_scale=22.0)
            if Config.SHOW_OUTBREAK_SEED_MARKERS and OUTBREAK_ENTRIES:
                orr = [e[0] for e in OUTBREAK_ENTRIES]; occ = [e[1] for e in OUTBREAK_ENTRIES]
                ax.scatter(np.array(occ) + 0.5, np.array(orr) + 0.5,
                           s=Config.OUTBREAK_MARKER_SIZE, marker="*",
                           c=Config.OUTBREAK_MARKER_COLOR, edgecolors="white",
                           linewidths=0.6, zorder=8, label="PPR outbreak points")
            ax.set_xlim(0, NC); ax.set_ylim(NR, 0); ax.set_xticks([]); ax.set_yticks([])
            cb = fig.colorbar(im, ax=ax, shrink=0.55, pad=0.02)
            cb.set_label("PPR MaxEnt suitability (0–1)", fontsize=10, fontweight="bold")
            _bar_km = 50.0; _bar_cells = _bar_km / (PW / 1000.0)
            _x0, _y0 = NC * 0.62, NR * 0.95
            ax.plot([_x0, _x0 + _bar_cells], [_y0, _y0], "k-", lw=3, zorder=9)
            ax.text(_x0 + _bar_cells / 2, _y0 - NR * 0.012, f"{int(_bar_km)} km",
                    ha="center", va="bottom", fontsize=10, fontweight="bold")
            legend = [Line2D([0], [0], color="#d81b60", lw=1.6, label="MDL boundary")]
            if Config.SHOW_OUTBREAK_SEED_MARKERS and OUTBREAK_ENTRIES:
                legend.append(Line2D([0], [0], marker="*", linestyle="none",
                                     markerfacecolor=Config.OUTBREAK_MARKER_COLOR,
                                     markeredgecolor="white", markersize=12,
                                     label=f"PPR outbreak points ({len(OUTBREAK_ENTRIES)})"))
            if CORRIDORS:
                legend.append(Line2D([0], [0], marker="v", linestyle="none",
                                     markerfacecolor="#00e5ff", markeredgecolor="#003b46",
                                     markersize=9, label="Water-deer corridor gate"))
            if Config.SHOW_REGIONAL_BASEMAP and LAND_BASE.any():
                legend.append(Patch(facecolor=Config.BASEMAP_LAND_COLOR,
                                    edgecolor=Config.BASEMAP_LINE_COLOR, label="Neighbouring land"))
            ax.legend(handles=legend, loc="lower left", fontsize=10, frameon=True,
                      facecolor="white", edgecolor="#888")
            _mx_name = MAXENT_TIF.name if MAXENT_TIF is not None else "MaxEnt raster"
            ax.set_title("PPR MaxEnt suitability with outbreak points and MDL corridor gates\n"
                         f"integrated on the model grid ({_mx_name})",
                         fontsize=12, fontweight="bold")
            fp_mx = Config.OUT_DIR / "PPR_v9.3_MaxEnt_suitability_outbreaks.png"
            fig.savefig(fp_mx, dpi=150, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            logger.info(f"   PPR MaxEnt map: {fp_mx}")
            if Config.DL_DIR:
                try:
                    shutil.copy2(fp_mx, Config.DL_DIR / fp_mx.name)
                except Exception:
                    pass
        except Exception as e:
            logger.warning(f"   PPR MaxEnt map skipped: {e}")
    render_maxent_map()
    # ══ TRANSBOUNDARY MAPS: China → North Korea → (corridor gates) → South Korea
    def render_transboundary_maps():
        """China → North Korea → South Korea, modelled explicitly.
        On the working grid the China (MaxEnt) raster and the Korea (gorani)
        raster are separated by open sea, so infection cannot diffuse across the
        gap. We therefore stage the Korea leg: (1) find when the China epidemic
        reaches the cells nearest the North-Korean border; (2) introduce PPR into
        North Korea at that time; (3) run the RAD spread down the CONTIGUOUS
        peninsula, where it must now find one of the WATER-DEER CORRIDOR GATES to
        get into South Korea. This is the leg the v9.1 gating slows down, and it
        is given KOREA_EXTRA_YEARS of extra horizon for that reason.
        v9.3: this is also the run that feeds the SOUTH-KOREA COUNTY CSV — the
        table is written from `cum_k` as soon as the leg finishes, before any
        plotting, so a figure failure can never cost you the data."""
        res = rad_results_holder.get("res")
        if res is None or not KOREA_MASK.any():
            logger.info("   transboundary maps skipped (no Korea cells in domain).")
            return
        try:
            cum_list = res["monthly_cum"]; NM = len(cum_list)
            # (a) China → NK-border arrival month.
            if CHINA_MASK.any():
                d_to_korea = distance_transform_edt(~KOREA_MASK)
                _cd = d_to_korea[CHINA_MASK]
                _thr_d = np.percentile(_cd, 5) if _cd.size else 0.0
                border_cn = CHINA_MASK & (d_to_korea <= max(_thr_d, 1.0))
                chn_curve = np.array([float(np.nanmean(np.nan_to_num(cum_list[m])[CHINA_MASK]))
                                      for m in range(NM)]) * 100.0
                bord_curve = np.array([float(np.nanmean(np.nan_to_num(cum_list[m])[border_cn]))
                                       if border_cn.any() else 0.0 for m in range(NM)]) * 100.0
                _bi = np.where(bord_curve >= 5.0)[0]
                a_border = int(_bi[0]) if _bi.size else 0
            else:
                chn_curve = np.zeros(NM); a_border = 0
            if Config.KOREA_ONLY_INTRO_MONTH:                   # v9.10 Korea-only mode
                a_border = max(0, int(Config.KOREA_ONLY_INTRO_MONTH) - 1)
                logger.info(f"   v9.10 Korea-only mode: DPRK introduction fixed at month {a_border + 1}")
            # (b) Korea sub-model on a cropped, contiguous peninsula grid.
            rr, cc = np.where(KOREA_MASK); pad = 4
            r0, r1 = max(0, rr.min() - pad), min(NR, rr.max() + 1 + pad)
            c0, c1 = max(0, cc.min() - pad), min(NC, cc.max() + 1 + pad)
            sl = (slice(r0, r1), slice(c0, c1))
            KP_k = KOREA_MASK[sl].copy(); NK_k = NK_KR[sl].copy(); SK_k = SK_KR[sl].copy()
            if int(KP_k.sum()) < 20:
                logger.info("   transboundary maps skipped (Korea sub-grid too small).")
                return
            env_k = {
                "KP": KP_k, "NK": NK_k, "SK": SK_k,
                "MDL_BOUNDARY": MDL_BOUNDARY[sl],
                # v9.1: the corridor gates travel with the crop, so the peninsula
                # leg crosses ONLY at the water-deer corridors.
                "MDL_FENCE_FIELD": (MDL_FENCE_FIELD[sl] if MDL_FENCE_FIELD is not None else None),
                "GOR_EFF": np.nan_to_num(GOR_EFF)[sl], "GOR_nn": GOR_nn[sl],
                "CONDUCTANCE": CONDUCTANCE[sl], "RESISTANCE": RESISTANCE[sl],
                "RIV_nn": RIV_nn[sl], "WATER_ATTRACT": WATER_ATTRACT[sl],
                "livestock_nn": livestock_nn[sl],
                "GORANI_DENSITY": np.nan_to_num(GORANI_DENSITY)[sl],
            }
            # Seed the northernmost North-Korean cells = introduction from China.
            seeds = []
            nk_rows = np.where(NK_k.any(axis=1))[0]
            if nk_rows.size:
                top = nk_rows.min()
                for _r in range(top, min(top + 3, NK_k.shape[0])):
                    cols = np.where(NK_k[_r])[0]
                    if cols.size:
                        step = max(1, len(cols) // 6)
                        for _c in cols[::step]:
                            seeds.append((int(_r), int(_c), "china_border", 0, 1.0))
            if not seeds:
                logger.info("   transboundary maps skipped (no North-Korea seed cells).")
                return
            env_k["DMZ_ENTRIES"] = seeds; env_k["NK_ENTRIES"] = []
            # ══ v9.10: ROK area weights, exclaves and the unfenced contact on the crop ══
            _band_k = MDL_BOUNDARY[sl].astype(bool)
            _W_k = (np.where(KP_k, np.clip(np.nan_to_num(ROK_FRACTION[sl]), 0.0, 1.0), 0.0)
                    if ROK_FRACTION is not None else SK_k.astype(np.float64))
            _EXCL_k = ROK_EXCLAVES[sl].copy()
            _rc_main, _dc_main = unfenced_contact(SK_k & ~_EXCL_k, NK_k, _band_k)
            # sealing needs EVERY ROK–DPRK face outside the band, including those next to it
            _rs0, _ds0 = unfenced_contact(SK_k & ~_EXCL_k, NK_k, _band_k, buffer_cells=0)
            KOREA_LEG.update(
                ROK_W=_W_k, EXCL=_EXCL_k, CONTACT_ROK=_rc_main, CONTACT_DPRK=_dc_main,
                CONTACT_SEAL=(_rs0 | _ds0),
                BAND=_band_k.copy(), MDL_ORIG=_band_k.copy(),
                FENCE_ORIG=(None if env_k["MDL_FENCE_FIELD"] is None
                            else np.asarray(env_k["MDL_FENCE_FIELD"], np.float64).copy()))
            if str(Config.KOREA_ESTUARY_POLICY).lower() == "sealed":
                _contact_k = (_rs0 | _ds0) & KP_k
                env_k["MDL_BOUNDARY"], env_k["MDL_FENCE_FIELD"] = seal_contact_band(
                    _band_k, env_k["MDL_FENCE_FIELD"], _contact_k,
                    sealed=float(Config.MDL_SEALED_FACTOR),
                    default_phi=float(Config.HAN_FENCE_REDUCTION))
                logger.info(f"   v9.10 estuary policy 'sealed': {int(_contact_k.sum())} contact cells "
                            f"added to the fence at φ0={Config.MDL_SEALED_FACTOR}.")
            else:
                logger.info(f"   v9.10 estuary policy 'open': {int(_rc_main.sum())} ROK and "
                            f"{int(_dc_main.sum())} DPRK cells meet across the unfenced border "
                            f"(Han River estuary)"
                            + (f"; {int((_EXCL_k & KP_k).sum())} ROK exclave cells kept"
                               if (_EXCL_k & KP_k).any() else "") + ".")
            # v9.1: the peninsula leg gets extra years, because the gated fence
            # deliberately delays the crossing into South Korea.
            nmk = 12 * max(1, int(Config.SIM_YEARS) + int(Config.KOREA_EXTRA_YEARS))
            if Config.KOREA_REPORT_MONTHS:                       # v9.10: reach every report month
                nmk = max(nmk, max(int(g) for g in Config.KOREA_REPORT_MONTHS) - int(a_border))
            p_k = RADParams(
                beta_monthly=np.asarray(Config.MONTHLY_BETA, float),
                disp_km_monthly=np.asarray(Config.DEER_RADIUS_KM, float) * Config.DEER_RADIUS_SCALE,
                tvir_monthly=np.asarray(T_VIRAL, float),
                cell_km=PW / 1000.0,
                diff_max=Config.RAD_DIFF_MAX * Config.SPREAD_SPEED_MULTIPLIER,
                diff_min=Config.RAD_DIFF_MIN,
                adv_max=Config.RAD_ADV_MAX * Config.SPREAD_SPEED_MULTIPLIER,
                adv_gain=Config.RAD_ADV_GAIN * Config.SPREAD_SPEED_MULTIPLIER,
                gamma=Config.GAMMA_RECOVERY,
                trans_hab_w=Config.RAD_TRANS_HAB_W, trans_lstk_w=Config.RAD_TRANS_LSTK_W,
                trans_base=Config.RAD_TRANS_BASE,
                seed_intensity=Config.RAD_SEED_INTENSITY,
                source_continuous=Config.RAD_SOURCE_CONTINUOUS, source_rate=Config.RAD_SOURCE_RATE,
                mdl_fence_factor=Config.HAN_FENCE_REDUCTION, apply_mdl_barrier=Config.HAN_APPLY_DMZ_FENCE,
                road_factor=Config.HAN_ROAD_REDUCTION, apply_road_barrier=False,
                corridor_seasonal_gating=Config.CORRIDOR_SEASONAL_GATING,
                corridor_season_floor=Config.CORRIDOR_SEASON_FLOOR,
                cfl_safety=Config.RAD_CFL_SAFETY, max_substeps=Config.RAD_MAX_SUBSTEPS,
                arrival_threshold=Config.RAD_ARRIVAL_THRESHOLD,
                vacc_enable=False,   # show the natural crossing (vaccination has its own maps)
                season_offset=int(a_border),   # v9.9 FIX A: leg month 0 = calendar month of NK introduction
            )
            # v9.2: run until PPR has crossed the MDL and ESTABLISHED in South
            # Korea (mean attack OR reached-area criterion), then keep going for
            # KOREA_POST_ARRIVAL_YEARS so the southward spread is visible. nmk is
            # only the minimum.
            _k_stop = _k_extra = _k_cap = None
            _k_lab = "criterion"
            if Config.AUTO_RUN_UNTIL_SK_ESTABLISHED and SK_k.any():
                # v9.10: establishment judged on ROK AREA (area-weighted, as in the sigungu tables)
                _wk_stop = np.asarray(KOREA_LEG["ROK_W"], np.float64)
                _ws_stop = float(_wk_stop.sum())
                def _k_stop(_m, _cum, _w=_wk_stop, _ws=_ws_stop):
                    if _ws <= 0:
                        return False
                    _c = np.nan_to_num(_cum)
                    return (float((_w * _c).sum()) / _ws >= Config.SK_ESTABLISHED_MEAN_ATTACK
                            or float((_w * (_c >= Config.RAD_ARRIVAL_THRESHOLD)).sum()) / _ws
                            >= Config.SK_ESTABLISHED_AREA_FRACTION)
                _k_extra = (int(Config.KOREA_POST_ARRIVAL_MONTHS)
                            if Config.KOREA_POST_ARRIVAL_MONTHS is not None
                            else 12 * max(0, int(Config.KOREA_POST_ARRIVAL_YEARS)))
                _k_cap = 12 * max(int(Config.MAX_SIM_YEARS),
                                  int(Config.SIM_YEARS) + int(Config.KOREA_EXTRA_YEARS))
                _k_lab = (f"PPR established in South Korea "
                          f"(mean attack ≥{100*Config.SK_ESTABLISHED_MEAN_ATTACK:.0f}% or "
                          f"≥{100*Config.SK_ESTABLISHED_AREA_FRACTION:.0f}% of SK reached)")
            _k_note = (f", extending until {_k_lab}, cap {_k_cap // 12} yr"
                       if _k_stop else "")
            logger.info(f"   [Korea sub-model] {int(KP_k.sum()):,} cells, "
                        f"{len(seeds)} NK border seeds, China reaches the NK border "
                        f"~month {a_border + 1}; minimum {nmk} months{_k_note}, "
                        f"through {len(CORRIDORS)} corridor gates...")
            res_k = run_dynamic_rad_simulation(env_k, p_k, n_months=nmk, closure=None,
                                               stop_fn=_k_stop,
                                               stop_extra_months=_k_extra or 0,
                                               max_months=_k_cap, stop_label=_k_lab)
            cum_k = res_k["monthly_cum"]
            nmk = len(cum_k)      # ACTUAL length after auto-extension
            KOREA_LEG.update(cum=cum_k, res=res_k, env=env_k, p=p_k, seeds=list(seeds), sl=sl,   # v9.10
                             r0=int(r0), c0=int(c0), a_border=int(a_border), nmk=int(nmk),
                             KP=KP_k, NK=NK_k, SK=SK_k, vacc={},
                             stop_month=(None if res_k.get("stop_month") is None
                                         else int(res_k["stop_month"]) + 1 + int(a_border)))
            # v9.1: WHICH gate delivered PPR into South Korea, and when.
            if CORRIDORS:
                corridor_crossing_report(
                    CORRIDORS, cum_k, SK_k, cell_km=PW / 1000.0,
                    probe_km=Config.CORRIDOR_PROBE_KM,
                    threshold=Config.RAD_ARRIVAL_THRESHOLD,
                    row_offset=r0, col_offset=c0, month_offset=a_border,
                    logger=logger)
            # ══ v9.3: THE SOUTH-KOREA COUNTY CSV ═════════════════════════════
            # Written HERE, immediately after the cross-MDL leg and BEFORE any
            # figure code, so the table survives even if plotting fails. Months
            # are offset by a_border so the CSV shares one timeline with the
            # China source run.
            # The cache goes FIRST: if the county layer turns out to be missing
            # or misconfigured, this 20-year simulation is still on disk and the
            # table can be rebuilt in seconds instead of re-run overnight.
            save_epidemic_cache(cum_k, KP_k, "korea_leg", r0, c0, a_border,
                                sk_mask=SK_k, monthly_I=res_k.get("monthly_I"),
                                extra=dict(nk_seed_month=int(a_border) + 1,
                                           n_gates=len(CORRIDORS)))
            export_sk_counties(
                cum_k, KP_k, tag="korea_leg",
                title=(f"Korea leg — PPR crosses the MDL through "
                       f"{len(CORRIDORS)} water-deer corridor gates "
                       f"(North Korea seeded at month {a_border + 1})"),
                row_offset=r0, col_offset=c0, month_offset=a_border,
                sk_mask=SK_k, monthly_I=res_k.get("monthly_I"),
                vacc_note="OFF in this leg (V1/V2 are separate runs: see the v9.10 Table 3)")
            # ══ v9.9 E: VACCINATION SCENARIOS ON THE KOREA LEG ═══════════════════
            # V1 = combined host on the whole peninsula from a fixed leg month (upper bound).
            # V2 = ROK only, triggered by the first ROK detection, coverage capped by the
            #      livestock share of local transmissibility (Eq. 4) — the realistic case.
            # Same length as the baseline leg, no stop rule, so the runs are directly comparable.
            if Config.KOREA_VACC_SCENARIOS:
                import dataclasses as _dc
                _rok_k = (SK_k & KP_k).astype(bool)
                _hab = np.clip(np.nan_to_num(np.asarray(env_k["GOR_EFF"], float)), 0.0, 1.0)
                _lst = np.clip(np.nan_to_num(np.asarray(env_k["livestock_nn"], float)), 0.0, 1.0)
                _share = np.where(KP_k, (Config.RAD_TRANS_LSTK_W * _lst) /
                                  (Config.RAD_TRANS_HAB_W * _hab + Config.RAD_TRANS_LSTK_W * _lst + 1e-9), 0.0)
                _scen = {
                    "V1_combined_host_peninsula": dict(field=KP_k, trigger=None, cap=None,
                                                       start=int(Config.KOREA_VACC_V1_START_LEG_MONTH)),
                    # v9.10: the campaign is triggered by a crossing, not by an exclave island
                    "V2_ROK_livestock_on_detection": dict(
                        field=_rok_k, trigger=(_rok_k & ~np.asarray(KOREA_LEG.get("EXCL", np.zeros_like(_rok_k)), bool)),
                        cap=_share, start=0),
                }
                _rows_csv = ["scenario,global_month,mean_attack_NK_pct,mean_attack_ROK_pct"]
                _base_nk = [float(np.nanmean(np.nan_to_num(c)[NK_k])) * 100 for c in cum_k]
                _base_sk = [float(np.nanmean(np.nan_to_num(c)[_rok_k])) * 100 for c in cum_k]
                for _i in range(len(cum_k)):
                    _rows_csv.append(f"baseline,{_i + 1 + a_border},{_base_nk[_i]:.3f},{_base_sk[_i]:.3f}")
                for _nm, _sc in _scen.items():
                    env_v = dict(env_k)
                    env_v["VACC_FIELD"] = np.asarray(_sc["field"], np.float32)
                    p_v = _dc.replace(
                        p_k, vacc_enable=True, vacc_mode=Config.VACC_MODE,
                        vacc_start_month=int(_sc["start"]),
                        vacc_coverage=float(Config.VACC_COVERAGE), vacc_rate=float(Config.VACC_RATE),
                        vacc_efficacy=float(Config.VACC_EFFICACY),
                        vacc_ring_cells=max(1, int(round(Config.VACC_RING_KM / max(PW / 1000.0, 1e-6)))),
                        vacc_inf_thr=float(Config.VACC_INFECTED_THRESHOLD),
                        vacc_trigger_mask=_sc["trigger"], vacc_cap_field=_sc["cap"])
                    logger.info(f"   [Korea leg — {_nm}] {len(cum_k)} months, same seeds and calendar as the baseline")
                    res_v = run_dynamic_rad_simulation(env_v, p_v, n_months=len(cum_k), closure=None)
                    cum_v = res_v["monthly_cum"]
                    _von = res_v.get("vacc_on_month")
                    KOREA_LEG.setdefault("vacc", {})[_nm] = dict(            # v9.10: for Table 3 / Fig. 9
                        cum=cum_v, von=_von, start=int(_sc["start"]), triggered=_sc["trigger"] is not None)
                    save_epidemic_cache(cum_v, KP_k, f"korea_leg_{_nm}", r0, c0, a_border,
                                        sk_mask=SK_k, monthly_I=res_v.get("monthly_I"),
                                        extra=dict(scenario=_nm,
                                                   vacc_on_global_month=(None if _von is None and _sc["trigger"] is not None
                                                                         else int((_von if _von is not None else _sc["start"])) + 1 + int(a_border))))
                    _nk = [float(np.nanmean(np.nan_to_num(c)[NK_k])) * 100 for c in cum_v]
                    _sk = [float(np.nanmean(np.nan_to_num(c)[_rok_k])) * 100 for c in cum_v]
                    for _i in range(len(cum_v)):
                        _rows_csv.append(f"{_nm},{_i + 1 + a_border},{_nk[_i]:.3f},{_sk[_i]:.3f}")
                    logger.info(f"     end of leg: mean attack NK {_base_nk[-1]:.1f}% -> {_nk[-1]:.1f}%, "
                                f"ROK {_base_sk[-1]:.1f}% -> {_sk[-1]:.1f}%; campaign started at global month "
                                f"{(int(_von) + 1 + int(a_border)) if _von is not None else ('never' if _sc['trigger'] is not None else int(_sc['start']) + 1 + int(a_border))}")
                try:
                    _fpv = Config.OUT_DIR / f"{V910_TAG}_korea_leg_vaccination_scenarios.csv"
                    _fpv.write_text("\n".join(_rows_csv) + "\n", encoding="utf-8")
                    logger.info(f"   ✓ vaccination scenario curves: {_fpv}")
                except Exception as _e:
                    logger.warning(f"   ⚠ could not write the vaccination CSV: {_e}")
            def rmean_k(msk):
                return np.array([float(np.nanmean(np.nan_to_num(cum_k[m])[msk])) if msk.any() else 0.0
                                 for m in range(nmk)]) * 100.0
            nk_curve, sk_curve = rmean_k(NK_k), rmean_k(SK_k)
            x_k = np.arange(nmk) + a_border + 1
            def arrival_k(series, off, thr=5.0):
                idx = np.where(series >= thr)[0]
                return int(idx[0]) + off + 1 if idx.size else None
            a_nk = arrival_k(nk_curve, a_border); a_sk = arrival_k(sk_curve, a_border)
            _bi2 = np.where(chn_curve >= 5.0)[0]; a_chn = int(_bi2[0]) + 1 if _bi2.size else None
            band_cmap = ListedColormap(["#ffffcc", "#fed976", "#fd8d3c", "#e31a1c", "#800026"])
            band_norm = BoundaryNorm(list(Config.SPREAD_PCT_BANDS),
                                     ncolors=len(Config.SPREAD_PCT_BANDS) - 1, clip=True)
            months_full = np.arange(1, NM + 1)
            _mx_cmap = LinearSegmentedColormap.from_list("maxent_u", Config.MAXENT_CMAP, N=256)
            SUIT_k = np.nan_to_num(GOR_EFF)[sl]
            sig_seg_k = [seg - np.array([c0, r0]) for seg in sig_segments] if sig_segments else []
            # v9.2: water bodies / waterways for the Korea crop, clipped to the
            # SOUTH side of the MDL so the South-Korea panel carries water and
            # nothing else. County lines get the same treatment (off by default).
            _mdl_row_crop = float(DMZ_ROW - r0)
            def _clip_segments_south(segs, row_min):
                """Keep only the runs of each polyline lying south of `row_min`
                (crop-row units); split at the crossings so no chord is drawn
                across the MDL when a polygon straddles it."""
                out = []
                for s in segs:
                    s = np.asarray(s, dtype=float)
                    if s.ndim != 2 or s.shape[0] < 2:
                        continue
                    idx = np.where(s[:, 1] >= row_min)[0]
                    if idx.size < 2:
                        continue
                    for run in np.split(idx, np.where(np.diff(idx) != 1)[0] + 1):
                        if run.size >= 2:
                            out.append(s[run])
                return out
            def _filter_water_segments(segs, min_extent_cells, max_features, step):
                """Keep only water features big enough to be RESOLVABLE on this
                grid, largest first, and thin their vertices.
                water.shp carries the full national hydrography (thousands of
                small streams). At multi-km cells those cannot be read as
                individual features and simply flood the panel with blue, so
                anything whose bounding box spans less than `min_extent_cells`
                is dropped."""
                scored = []
                for s in segs:
                    s = np.asarray(s, dtype=float)
                    if s.shape[0] < 2:
                        continue
                    ext = max(float(s[:, 0].max() - s[:, 0].min()),
                              float(s[:, 1].max() - s[:, 1].min()))
                    if ext < float(min_extent_cells):
                        continue
                    if step > 1 and s.shape[0] > 2 * step:
                        s = np.vstack([s[::int(step)], s[-1:]])
                    scored.append((ext, s))
                scored.sort(key=lambda t: -t[0])
                if max_features and len(scored) > int(max_features):
                    scored = scored[:int(max_features)]
                return [s for _, s in scored]
            # v9.10: clip the ROK panel's lines to ROK territory, not to the MDL row
            _rok_draw = ((SK_k | (binary_dilation(SK_k, structure=V910_FOUR) & ~NK_k))
                         if ROK_RASTER is not None else None)
            if water_segments and Config.SK_PANEL_SHOW_WATER:
                _wsk = [seg - np.array([c0, r0]) for seg in water_segments]
                _water_south = (clip_segments_to_mask(_wsk, _rok_draw) if _rok_draw is not None
                                else _clip_segments_south(_wsk, _mdl_row_crop))
            else:
                _water_south = []
            water_seg_k = _filter_water_segments(
                _water_south,
                Config.SK_PANEL_WATER_MIN_EXTENT_CELLS,
                Config.SK_PANEL_WATER_MAX_FEATURES,
                Config.SK_PANEL_WATER_SIMPLIFY_STEP)
            sig_seg_k = ((clip_segments_to_mask(sig_seg_k, _rok_draw) if _rok_draw is not None
                          else _clip_segments_south(sig_seg_k, _mdl_row_crop)) if sig_seg_k else [])
            logger.info(f"   [SK panel] water polylines south of the MDL: "
                        f"{len(_water_south)} found -> {len(water_seg_k)} drawn "
                        f"(kept features spanning ≥{Config.SK_PANEL_WATER_MIN_EXTENT_CELLS:.1f} cells "
                        f"≈ {Config.SK_PANEL_WATER_MIN_EXTENT_CELLS * PW / 1000.0:.0f} km, "
                        f"cap {Config.SK_PANEL_WATER_MAX_FEATURES}); "
                        f"counties {'ON' if Config.SK_PANEL_SHOW_COUNTIES else 'OFF'}; "
                        f"suitability underlay "
                        f"{'NK only' if Config.SK_PANEL_PLAIN else 'whole peninsula'}")
            # ══ v9.7: hand the Korea-leg epidemic to the SMCDA validator ═════
            # The Korea-domain risk classifications must be tested against THIS
            # run — the one that carries PPR from the Chinese border down through
            # North Korea and across the gated MDL — not against the continental
            # run, where the peninsula is reached last and arrival order there
            # says nothing useful about risk WITHIN Korea.
            #
            # The Korea leg is computed on a CROP of the model grid, so its
            # arrays are re-embedded on the full grid here. Cells outside the
            # crop stay NaN and are simply absent from the validation.
            try:
                _tta_full = np.full(SH, np.nan, np.float32)
                _cum_full = np.full(SH, np.nan, np.float32)
                _tta_k = res_k.get("time_to_arrival")
                if _tta_k is not None:
                    _tta_full[sl] = np.asarray(_tta_k, np.float32)
                _cum_full[sl] = np.nan_to_num(np.asarray(cum_k[-1], np.float32))
                SMCDA_KOREA_EPI.update(
                    time_to_arrival=_tta_full, cum_attack=_cum_full,
                    month_offset=int(a_border), n_months=int(nmk),
                    sl=sl, r0=int(r0), c0=int(c0),
                    KP_k=KP_k, NK_k=NK_k, SK_k=SK_k,
                    sig_seg_k=list(sig_seg_k), water_seg_k=list(water_seg_k),
                    SUIT_k=SUIT_k, mx_cmap=_mx_cmap, seeds=list(seeds))
                logger.info("   [SMCDA] Korea-leg epidemic registered for validation "
                            f"({int(np.isfinite(_tta_full).sum()):,} cells with an "
                            f"arrival month; month offset {int(a_border)}).")
            except Exception as _e:
                logger.warning(f"   [SMCDA] could not register the Korea-leg epidemic "
                               f"({type(_e).__name__}: {_e}); the Korea risk classes "
                               f"will be reported WITHOUT validation.")
            def draw_korea(ax, pct_crop, title, gate_labels=False):
                """Korea-crop frame (v9.2).
                NORTH KOREA keeps the habitat-suitability underlay (it is the
                source region). SOUTH KOREA is a BLANK WHITE canvas: the only
                colour south of the MDL is the PPR cumulative-attack band, over
                water bodies / waterways drawn as thin blue lines.
                `pct_crop` is the RAW % attack grid — the visibility floors are
                applied HERE, separately for NK and SK, so a 2% arrival in the
                south is still drawn while the northern panel keeps its 10% floor.
                """
                ax.set_facecolor("white")
                ext = (0, KP_k.shape[1], KP_k.shape[0], 0)
                ones_k = np.ones(KP_k.shape, float)
                _bn_k = BoundaryNorm(list(Config.SPREAD_PCT_BANDS),
                                     ncolors=len(Config.SPREAD_PCT_BANDS) - 1, clip=True)
                # ── (1) LAND CANVAS ─────────────────────────────────────────
                # THIS is what made South Korea brown in v9.1: the underlay was
                # drawn over KP_k (both Koreas). It is now NK-only.
                if Config.SK_PANEL_PLAIN:
                    ax.imshow(np.ma.masked_where(~SK_k, ones_k),
                              cmap=ListedColormap([Config.SK_PANEL_LAND_COLOR]),
                              vmin=0, vmax=1, extent=ext, origin="upper", aspect="equal",
                              interpolation="nearest", zorder=1.0)
                    hab_mask = NK_k if Config.NK_PANEL_SHOW_SUITABILITY else np.zeros_like(KP_k)
                else:
                    hab_mask = KP_k
                if np.any(hab_mask):
                    ax.imshow(np.ma.masked_invalid(np.where(hab_mask, SUIT_k, np.nan)),
                              cmap=_mx_cmap, vmin=0, vmax=1, extent=ext, origin="upper",
                              aspect="equal", alpha=0.80, interpolation="bilinear", zorder=1.2)
                # ── (2) NORTH KOREA % attack — warm ramp ────────────────────
                nk_pct = np.where(NK_k & (pct_crop >= Config.BAND_MIN_PCT_VISIBLE),
                                  pct_crop, np.nan)
                ax.imshow(np.ma.masked_invalid(nk_pct),
                          cmap=ListedColormap(Config.SPREAD_BAND_COLORS_NK, name="nk_k"),
                          norm=_bn_k, extent=ext, origin="upper", aspect="equal",
                          alpha=0.85, interpolation="nearest", zorder=3.0)
                # ── (3) SOUTH KOREA % attack — blue-green ramp, own floor ───
                sk_pct = np.where(SK_k & (pct_crop >= Config.SK_BAND_MIN_PCT_VISIBLE),
                                  pct_crop, np.nan)
                ax.imshow(np.ma.masked_invalid(sk_pct),
                          cmap=ListedColormap(Config.SPREAD_BAND_COLORS_SK, name="sk_k"),
                          norm=_bn_k, extent=ext, origin="upper", aspect="equal",
                          alpha=0.90, interpolation="nearest", zorder=3.1)
                # ── (4) SK coastline + WATER ONLY (counties off by default) ─
                if Config.SK_PANEL_OUTLINE and SK_k.any():
                    ax.contour(SK_k.astype(float), levels=[0.5],
                               colors=Config.SK_PANEL_OUTLINE_COLOR,
                               linewidths=Config.SK_PANEL_OUTLINE_WIDTH, zorder=4.0)
                if Config.SK_PANEL_SHOW_COUNTIES and sig_seg_k:
                    ax.add_collection(LineCollection(
                        sig_seg_k, colors=Config.SK_PANEL_COUNTY_LINE_COLOR,
                        linewidths=Config.SK_PANEL_COUNTY_LINE_WIDTH,
                        alpha=Config.SK_PANEL_COUNTY_LINE_ALPHA, zorder=4.1))
                if Config.SK_PANEL_SHOW_WATER and water_seg_k:
                    ax.add_collection(LineCollection(
                        water_seg_k, colors=Config.SK_PANEL_WATER_LINE_COLOR,
                        linewidths=Config.SK_PANEL_WATER_LINE_WIDTH,
                        alpha=Config.SK_PANEL_WATER_LINE_ALPHA, zorder=4.3))
                # ── (5) MDL ────────────────────────────────────────────────
                mdl_c = MDL_BOUNDARY[sl]
                if mdl_c.any():
                    ax.imshow(np.ma.masked_where(~mdl_c, np.ones_like(mdl_c, float)),
                              cmap=ListedColormap(["#d81b60"]), extent=ext, origin="upper",
                              aspect="equal", zorder=5)
                # ── (6) corridor gates + (7) China introduction seeds ──────
                if Config.SHOW_CORRIDOR_MARKERS and CORRIDORS:
                    plot_corridor_markers(ax, CORRIDORS, row_offset=r0, col_offset=c0,
                                          label=gate_labels, size_scale=18.0)
                sr = [s[0] for s in seeds]; sc = [s[1] for s in seeds]
                ax.scatter(np.array(sc) + 0.5, np.array(sr) + 0.5, s=18, marker="v",
                           c="#7a3b00", edgecolors="white", linewidths=0.4, zorder=6)
                # ── (8) hemisphere labels ──────────────────────────────────
                if NK_k.any():
                    ax.text(KP_k.shape[1] * 0.5, np.where(NK_k.any(axis=1))[0].min() + 1.5,
                            "NORTH KOREA", ha="center", fontsize=8,
                            fontweight="bold", color="#7a3b00")
                if SK_k.any():
                    ax.text(KP_k.shape[1] * 0.5, KP_k.shape[0] - 2, "SOUTH KOREA",
                            ha="center", fontsize=8, fontweight="bold", color="#0b3d6b")
                ax.set_xticks([]); ax.set_yticks([])
                ax.set_title(title, fontsize=11, fontweight="bold")
            KOREA_LEG["draw_korea"] = draw_korea                 # v9.10: reused for Figure 9
            # ── Figure 1: China→NK→SK time-series + final Korea zoom ─────────
            fig = plt.figure(figsize=(21, 9), facecolor="white")
            gsT = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0], wspace=0.16)
            axL = fig.add_subplot(gsT[0, 0])
            axL.plot(months_full, chn_curve, "-o", color="#d7191c", lw=2.2, ms=2.5, label="China (source)")
            axL.plot(x_k, nk_curve, "-o", color="#fd8d3c", lw=2.2, ms=2.5, label="North Korea")
            axL.plot(x_k, sk_curve, "-o", color="#2c7bb6", lw=2.2, ms=2.5, label="South Korea")
            for a_, c_ in [(a_chn, "#d7191c"), (a_nk, "#fd8d3c"), (a_sk, "#2c7bb6")]:
                if a_:
                    axL.axvline(a_, color=c_, ls="--", alpha=0.55)
            axL.set_xlabel("Month"); axL.set_ylabel("Mean cumulative attack rate (%)")
            axL.set_title("Transboundary spread: China → North Korea → (corridor gates) → South Korea",
                          fontsize=13, fontweight="bold")
            axL.grid(alpha=0.3); axL.legend(fontsize=11, loc="upper left")
            axL.set_xlim(1, max(NM, a_border + nmk))
            _fmt = lambda a: f"month {a}" if a else "not reached"
            _lag = (a_sk - a_nk) if (a_sk and a_nk) else None
            axL.text(0.02, 0.80,
                     f"Arrival (≥5% mean attack):\n  China: {_fmt(a_chn)}\n  North Korea: {_fmt(a_nk)}"
                     f"\n  South Korea: {_fmt(a_sk)}\n"
                     f"  NK→SK crossing lag: {(str(_lag) + ' months') if _lag else 'n/a'}\n"
                     f"(NK seeded when China reaches the border, month {a_border+1};\n"
                     f" MDL gated at {len(CORRIDORS)} water-deer corridors)",
                     transform=axL.transAxes, va="top", fontsize=10,
                     bbox=dict(fc="white", ec="#999", alpha=0.9))
            axR = fig.add_subplot(gsT[0, 1])
            # v9.2: pass the RAW % grid — draw_korea applies the NK/SK floors
            # separately (a single 10% pre-filter here also blanked genuine
            # low-% arrivals in South Korea).
            final_pct = np.where(KP_k, np.nan_to_num(cum_k[-1]) * 100.0, np.nan)
            draw_korea(axR, final_pct, f"Korea — cumulative attack, end of run (month {a_border + nmk})",
                       gate_labels=True)
            # Split colorbars: NK warm on the left, SK blue-green on the right.
            _bn = BoundaryNorm(list(Config.SPREAD_PCT_BANDS),
                               ncolors=len(Config.SPREAD_PCT_BANDS) - 1, clip=True)
            _fig_bbox = axR.get_position()
            cax_nk = fig.add_axes([_fig_bbox.x1 + 0.01, _fig_bbox.y0 + _fig_bbox.height * 0.55,
                                   0.013, _fig_bbox.height * 0.38])
            sm_nk = plt.cm.ScalarMappable(
                cmap=ListedColormap(Config.SPREAD_BAND_COLORS_NK, name="nk_cb"), norm=_bn)
            sm_nk.set_array([])
            cb_nk = fig.colorbar(sm_nk, cax=cax_nk, boundaries=list(Config.SPREAD_PCT_BANDS),
                                 ticks=list(Config.SPREAD_PCT_BANDS), spacing="proportional")
            cb_nk.ax.set_yticklabels([f"{int(b)}%" for b in Config.SPREAD_PCT_BANDS], fontsize=7)
            cb_nk.set_label("NK %", fontsize=8, fontweight="bold")
            cax_sk = fig.add_axes([_fig_bbox.x1 + 0.028, _fig_bbox.y0 + _fig_bbox.height * 0.55,
                                   0.013, _fig_bbox.height * 0.38])
            sm_sk = plt.cm.ScalarMappable(
                cmap=ListedColormap(Config.SPREAD_BAND_COLORS_SK, name="sk_cb"), norm=_bn)
            sm_sk.set_array([])
            cb_sk = fig.colorbar(sm_sk, cax=cax_sk, boundaries=list(Config.SPREAD_PCT_BANDS),
                                 ticks=list(Config.SPREAD_PCT_BANDS), spacing="proportional")
            cb_sk.ax.set_yticklabels([f"{int(b)}%" for b in Config.SPREAD_PCT_BANDS], fontsize=7)
            cb_sk.set_label("SK %", fontsize=8, fontweight="bold")
            fig.suptitle("PPR transboundary spread — China → North Korea (warm) → corridor gates "
                         "→ South Korea (blue-green on white, with counties + water)\n"
                         f"MDL sealed={Config.MDL_SEALED_FACTOR}, gate peak=Han φ={Config.CORRIDOR_MAX_PASS}×weight, "
                         f"radius={Config.CORRIDOR_GATE_RADIUS_KM:.0f} km",
                         fontsize=13, fontweight="bold")
            fp_tb = Config.OUT_DIR / "PPR_v9.3_transboundary_china_nk_sk.png"
            fig.savefig(fp_tb, dpi=130, facecolor="white", bbox_inches="tight"); plt.close(fig)
            logger.info(f"   transboundary China→NK→SK map: {fp_tb} "
                        f"(China={a_chn}, NK={a_nk}, SK={a_sk}, NK→SK lag={_lag})")
            if Config.DL_DIR:
                try:
                    shutil.copy2(fp_tb, Config.DL_DIR / fp_tb.name)
                except Exception:
                    pass
            # ── Figure 2: NK → SK spread progression (Korea crop) ─────────────
            snaps = sorted(set(np.linspace(0, nmk - 1, 12).round().astype(int).tolist()))
            figP, axP = plt.subplots(3, 4, figsize=(20, 16), facecolor="white")
            for k, mi in enumerate(snaps):
                a = axP[k // 4][k % 4]
                # v9.2: no pre-filter — draw_korea applies the NK/SK floors.
                pct = np.where(KP_k, np.nan_to_num(cum_k[mi]) * 100.0, np.nan)
                _gm = (mi + a_border) % 12
                _gy = (mi + a_border) // 12 + 1
                draw_korea(a, pct, f"{Config.MONTH_NAMES[_gm]} · Year {_gy}  (month {mi + a_border + 1})")
            for k in range(len(snaps), 12):
                axP[k // 4][k % 4].axis("off")
            # Two narrow colorbars on the right — NK warm, SK blue-green.
            _bn2 = BoundaryNorm(list(Config.SPREAD_PCT_BANDS),
                                ncolors=len(Config.SPREAD_PCT_BANDS) - 1, clip=True)
            cax2_nk = figP.add_axes([0.925, 0.35, 0.012, 0.55])
            sm2_nk = plt.cm.ScalarMappable(
                cmap=ListedColormap(Config.SPREAD_BAND_COLORS_NK), norm=_bn2)
            sm2_nk.set_array([])
            cb2_nk = figP.colorbar(sm2_nk, cax=cax2_nk, boundaries=list(Config.SPREAD_PCT_BANDS),
                                   ticks=list(Config.SPREAD_PCT_BANDS), spacing="proportional")
            cb2_nk.ax.set_yticklabels([f"{int(b)}%" for b in Config.SPREAD_PCT_BANDS], fontsize=7)
            cb2_nk.set_label("NK % attack (warm)", fontsize=8, fontweight="bold")
            cax2_sk = figP.add_axes([0.943, 0.35, 0.012, 0.55])
            sm2_sk = plt.cm.ScalarMappable(
                cmap=ListedColormap(Config.SPREAD_BAND_COLORS_SK), norm=_bn2)
            sm2_sk.set_array([])
            cb2_sk = figP.colorbar(sm2_sk, cax=cax2_sk, boundaries=list(Config.SPREAD_PCT_BANDS),
                                   ticks=list(Config.SPREAD_PCT_BANDS), spacing="proportional")
            cb2_sk.ax.set_yticklabels([f"{int(b)}%" for b in Config.SPREAD_PCT_BANDS], fontsize=7)
            cb2_sk.set_label("SK % attack (blue-green)", fontsize=8, fontweight="bold")
            if getattr(Config, "PUB_FIGS", False):
                # 11 Nov: DPRK ticks and label on the left of their bar, ROK on the right of theirs
                for _cb2, _side in ((cb2_nk, "left"), (cb2_sk, "right")):
                    _cb2.ax.yaxis.set_ticks_position(_side)
                    _cb2.ax.yaxis.set_label_position(_side)
                    _cb2.ax.tick_params(labelsize=10)
                    _cb2.ax.yaxis.label.set_fontsize(11)
            figP.suptitle("PPR spread down the Korean peninsula — North Korea (warm) → corridor gates "
                          "→ South Korea (blue-green on a white canvas with counties + water)\n"
                          "(▼ brown = China introduction; ▼ cyan = water-deer corridor gate, size ∝ weight; "
                          "magenta = MDL, sealed at %.2f, gates ≤ Han φ=%.2f×weight)"
                          % (Config.MDL_SEALED_FACTOR,
                             Config.CORRIDOR_MAX_PASS * Config.MDL_CROSSING_SLOWDOWN),
                          fontsize=14, fontweight="bold")
            fp_pr = Config.OUT_DIR / "PPR_v9.3_NK_to_SK_progression.png"
            figP.savefig(fp_pr, dpi=130, facecolor="white", bbox_inches="tight"); plt.close(figP)
            logger.info(f"   NK→SK progression map (Korea crop): {fp_pr}")
            if Config.DL_DIR:
                try:
                    shutil.copy2(fp_pr, Config.DL_DIR / fp_pr.name)
                except Exception:
                    pass
        except Exception as e:
            logger.warning(f"   transboundary maps skipped: {e}", exc_info=True)
    render_transboundary_maps()
    # ══ v9.7: SMCDA VALIDATION, MAPS AND COUNTY TABLE (all domains) ══════════
    # Placed AFTER the transboundary maps so the Korea leg has run: the Korean
    # domains are validated against the epidemic that actually enters the
    # peninsula, and the continental domain against the continental run.
    # The PPR-invaded footprint is only knowable now that the Korea leg has run,
    # so the epidemic-conditioned domains are built here rather than up front.
    if Config.ENABLE_SMCDA:
        try:
            compute_deferred_smcda_domains()
        except Exception as _e:
            logger.warning(f"  ⚠ the epidemic-conditioned SMCDA domains could not be "
                           f"built ({type(_e).__name__}: {_e}). The unrestricted domains "
                           f"are unaffected.")
    def _shift_segments(segs, dx, dy):
        """Translate LineCollection segments stored in crop coordinates.

        Used only by the v9.8 South-Korea frame crop. The county and waterway
        outlines are cached once, in the coordinates of the Korea-leg crop; when
        the panel is re-cropped to South Korea they have to move with it or they
        will be drawn a few hundred kilometres off the coastline. Segments are
        (x, y) = (column, row) pairs, so dx shifts columns and dy shifts rows.
        No segment is dropped — matplotlib clips whatever falls outside the
        axes, and clipping a coastline is correct where discarding a segment
        that merely starts outside the frame would leave a visible gap."""
        out = []
        for seg in segs or []:
            try:
                a = np.asarray(seg, float)
                if a.ndim != 2 or a.shape[1] != 2 or a.shape[0] < 2:
                    continue
                b = a.copy()
                b[:, 0] += float(dx)
                b[:, 1] += float(dy)
                out.append(b)
            except Exception:
                continue
        return out
    def _smcda_epidemic_for(dk):
        """(time_to_arrival, cum_attack, source label) for a domain — or None.

        WHICH EPIDEMIC A DOMAIN IS TESTED AGAINST IS PART OF THE TEST. A South-
        Korea classification asks 'within South Korea, where does PPR arrive
        first?'. Only the Korea leg answers that; in the continental run South
        Korea is reached last and almost simultaneously, so arrival order there
        carries no information about risk within the country."""
        dk = str(dk).lower()
        if (dk in ("sk", "korea", "nk", "sk_reached", "korea_reached")
                and SMCDA_KOREA_EPI.get("time_to_arrival") is not None):
            return (SMCDA_KOREA_EPI["time_to_arrival"], SMCDA_KOREA_EPI["cum_attack"],
                    "the Korea leg (China border -> North Korea -> corridor gates "
                    "-> South Korea)")
        _res = rad_results_holder.get("res")
        if _res is None or not _res.get("monthly_cum"):
            return None
        _tta = _res.get("time_to_arrival")
        if _tta is None:
            return None
        if dk in ("sk", "korea", "nk", "sk_reached", "korea_reached"):
            return (np.asarray(_tta, np.float64),
                    np.asarray(_res["monthly_cum"][-1], np.float64),
                    "the CONTINENTAL run (the Korea leg was unavailable — arrival "
                    "order within Korea is weakly informative here)")
        return (np.asarray(_tta, np.float64),
                np.asarray(_res["monthly_cum"][-1], np.float64),
                "the continental run")

    # ══ v9.10 ITEM 2: entry origins for the benchmark rankings ═══════════════════
    def v910_entry_origins(dk, dom_mask):
        """(origin masks, cost distances) used by smcda_benchmarks() for one domain."""
        dk = str(dk).lower()
        if dk in ("sk", "sk_reached"):
            srcs = (("the corridor gates", "gates"),
                    ("gates + unfenced border (estuary)", "gates_terminus"),
                    ("the MDL line", "mdl"))
        elif dk in ("korea", "korea_reached", "nk"):
            srcs = (("the China–Korea border", "china_border"), ("the corridor gates", "gates"),
                    ("gates + unfenced border (estuary)", "gates_terminus"))
        else:
            srcs = (("the outbreak seed points", "seeds"),)
        origins, costd = {}, {}
        for lab, src in srcs:
            try:
                m_, _d = _smcda_entry_mask(dk, src_override=src)
            except Exception as _e:
                logger.info(f"    (benchmark origin '{lab}' unavailable: {_e})")
                m_ = None
            if m_ is None or not np.any(m_):
                continue
            origins[lab] = m_
            if src != "mdl":
                costd[lab] = _smcda_gate_distance(m_, dom_mask, lab)
        return origins, costd

    def v910_write_smcda_tables(dk, sm):
        """Table 2 and Supplementary Table S7 for one SMCDA domain."""
        files = []
        if sm.get("table2"):
            files.append(v910_write_csv(
                Config.OUT_DIR / f"{V910_TAG}_Table2_smcda_classes_{dk}.csv",
                ["cls", "label", "score_lower", "score_upper", "cells", "share_pct",
                 "mean_arrival_month", "median_arrival_month", "reached_pct", "mean_attack_pct"],
                sm["table2"],
                [f"{V910_TAG} Table 2 — SMCDA classes, domain '{dk}' ({sm.get('domain_label')})",
                 f"ENTRY measured from: {sm.get('entry_basis')}",
                 f"ROK/DPRK split: {ROK_SPLIT_INFO.get('method')}",
                 "share_pct rounded by largest remainder so the shares add up to 100.0; "
                 "months on the global index (month 1 = January of Year 1)"], logger))
        if sm.get("benchmarks"):
            files.append(v910_write_csv(
                Config.OUT_DIR / f"{V910_TAG}_TableS7_benchmarks_{dk}.csv",
                ["benchmark", "kind", "rho_s", "auc_reached", "n_arrived", "pct_same_class",
                 "monotonic", "mean_arrival_by_class", "weights", "note"],
                sm["benchmarks"],
                [f"{V910_TAG} item 2 — rankings vs simulated arrival, domain '{dk}'",
                 "rho_s: Spearman correlation (average ranks for ties) between the ranking and the "
                 "arrival month; negative = higher-ranked cells are reached earlier",
                 "distance rankings use minus the distance; AUC only where some cells were never reached",
                 "pct_same_class: share of cells in the same class as the reported AHP map "
                 "(each alternative re-classified with the same break method)",
                 "the thresholds JSON (spearman_score_vs_arrival) breaks ties by array order; quote "
                 "the tie-averaged rho_s from this table"], logger))
            logger.info(f"    [v9.10 benchmarks, domain '{dk}'] Spearman rho with arrival month:")
            for _r in sm["benchmarks"]:
                if _r["kind"] in ("composite", "distance") or _r["benchmark"].startswith(
                        ("equal", "entropy", "rank-sum")):
                    logger.info(f"      {_r['benchmark'][:64]:<64s} rho_s = {_r['rho_s']:+.3f}")
        if sm.get("beyond_distance"):
            _bd = sm["beyond_distance"]
            files.append(v910_write_csv(
                Config.OUT_DIR / f"{V910_TAG}_TableS7b_beyond_distance_{dk}.csv",
                ["ranking", "kind", "n_cells", "reference", "rho_s", "rho_s_reference",
                 "partial_rho_given_distance", "partial_rho_given_distance_lo", "partial_rho_given_distance_hi",
                 "r2_reference", "r2_reference_plus_ranking", "delta_r2", "delta_r2_lo", "delta_r2_hi",
                 "best_lambda", "rho_s_best_mixture", "mixture_gain", "mixture_gain_lo", "mixture_gain_hi",
                 "bootstrap"],
                _bd,
                [f"{V910_TAG} round 2 (Reviewer 2, F1) — what each ranking adds to distance from the entry "
                 f"routes, domain '{dk}' ({sm.get('domain_label')}); cells with an arrival month",
                 "rho_s: Spearman rho with arrival month (average ranks; negative = higher-ranked cells reached "
                 "earlier); rho_s_reference: the same for minus the Euclidean distance from the reference origin",
                 "partial_rho_given_distance: rho_s with that distance held fixed (partial Spearman); 0 = the "
                 "ranking adds nothing to distance; negative = adds the right order; positive = the wrong order",
                 "r2_*: R² of rank(arrival) on rank(distance), and on rank(distance) + rank(ranking); delta_r2 = "
                 "their difference",
                 "best_lambda: weight of the ranking in the mixture (1-λ)·distance + λ·ranking (percentile "
                 "ranks, λ = 0, 0.05, ..., 1) with the most negative rho_s; mixture_gain = rho_s_reference - "
                 "rho_s_best_mixture ≥ 0 (in-sample, so an upper bound); λ = 0 means distance alone is best",
                 "_lo/_hi: 95% spatial block-bootstrap intervals, λ chosen again in every replicate"],
                logger))
            logger.info(f"    [round 2, beyond distance, domain '{dk}'] reference: {_bd[0]['reference']} "
                        f"(rho_s {_bd[0]['rho_s']:+.3f}, R² {_bd[0]['r2_reference']:.3f})")
            for _r in _bd[1:]:
                _ci = (f" ({_r['partial_rho_given_distance_lo']:+.2f} to {_r['partial_rho_given_distance_hi']:+.2f})"
                       if "partial_rho_given_distance_lo" in _r else "")
                logger.info(f"      {_r['ranking'][:58]:<58s} rho_s {_r['rho_s']:+.3f}; given distance "
                            f"{_r['partial_rho_given_distance']:+.3f}{_ci}; ΔR² {_r['delta_r2']:+.3f}; best "
                            f"λ {_r['best_lambda']:.2f}, gain {_r['mixture_gain']:.3f}")
        v910_copy_to_dl(files, Config.DL_DIR, logger)
        return files

    def validate_smcda_domains():
        """Re-run the validation now that the epidemics exist, and re-print each
        threshold table with the epidemiological evidence filled in.

        Every score was built from static layers BEFORE any simulation, so these
        are out-of-sample checks on the classifications, not measures of fit."""
        if not SMCDA_BY_DOMAIN:
            return
        if not Config.SMCDA_VALIDATE_AGAINST_EPIDEMIC:
            logger.info("\n[SMCDA VALIDATION] disabled "
                        "(SMCDA_VALIDATE_AGAINST_EPIDEMIC=False).")
            return
        logger.info("\n" + "=" * 96)
        logger.info("[SMCDA VALIDATION AGAINST THE MECHANISTIC EPIDEMIC]")
        logger.info("  Each risk score uses static layers only and never saw an epidemic, so")
        logger.info("  the statistics below are genuine tests of the classifications rather")
        logger.info("  than restatements of their fit.")
        if Config.USE_SMCDA_CLASS_IN_TRANSPORT:
            logger.warning(
                "  ⚠ USE_SMCDA_CLASS_IN_TRANSPORT=True, so the classes DID enter the "
                "epidemic. The independence claim above does NOT hold for this run and "
                "these statistics are inflated. Re-run with the switch off before "
                "reporting them.")
        logger.info("=" * 96)
        for _dk, _sm in SMCDA_BY_DOMAIN.items():
            try:
                _epi = _smcda_epidemic_for(_dk)
                logger.info(f"\n  ── domain '{_dk}' — {_sm['domain_label']} ──")
                if _epi is None:
                    logger.warning("    ⚠ no epidemic available; the thresholds stand but "
                                   "nothing has tested whether they rank the modelled "
                                   "spread.")
                    continue
                _tta, _cum, _srclab = _epi
                if np.asarray(_tta).shape != SH or np.asarray(_cum).shape != SH:
                    logger.warning(f"    ⚠ validation skipped: epidemic arrays are "
                                   f"{np.asarray(_tta).shape}/{np.asarray(_cum).shape} but "
                                   f"the classification grid is {SH}.")
                    continue
                _ov = int((_sm["mask"] & np.isfinite(np.asarray(_tta, float))).sum())
                logger.info(f"    Validated against: {_srclab}")
                logger.info(f"    {_ov:,} of {int(_sm['mask'].sum()):,} domain cells have "
                            f"an arrival month in that run.")
                if _ov < 50:
                    logger.warning("    ⚠ too few cells were ever reached to validate "
                                   "meaningfully. Treat the verdicts below as absent, not "
                                   "as passes.")
                _sm["validation"] = validate_against_epidemic(
                    np.asarray(_sm["score"], np.float64), _sm["classes"], _sm["mask"],
                    time_to_arrival=np.asarray(_tta, np.float64),
                    cum_attack=np.asarray(_cum, np.float64),
                    arrival_threshold=float(Config.RAD_ARRIVAL_THRESHOLD),
                    class_labels=_sm["labels"])
                _sm["criterion_diagnostic"] = criterion_arrival_diagnostic(
                    _sm["standardised"], _sm["weights"], _sm["mask"],
                    np.asarray(_tta, np.float64), _sm["criteria"])
                _sm["validated_against"] = _srclab
                # ══ v9.10 ITEM 2: Table 2 + benchmarks / alternative weights ══════
                try:
                    _moff = (int(SMCDA_KOREA_EPI.get("month_offset", 0))
                             if (str(_dk).lower() in ("sk", "korea", "nk", "sk_reached", "korea_reached")
                                 and SMCDA_KOREA_EPI.get("time_to_arrival") is not None) else 0)
                    _sm["table2"] = smcda_class_table(_sm, _moff)
                    if Config.BENCHMARK_ALT_WEIGHTS:
                        _org, _cst = v910_entry_origins(_dk, _sm["mask"])
                        _sm["benchmarks"] = smcda_benchmarks(
                            _sm, np.asarray(_tta, np.float64), np.asarray(_cum, np.float64),
                            float(Config.RAD_ARRIVAL_THRESHOLD), _moff, _org, _cst, PW / 1000.0)
                        if getattr(Config, "BEYOND_DISTANCE_ANALYSIS", False):
                            try:                          # round 2: what each ranking adds to distance
                                _sm["beyond_distance"] = smcda_beyond_distance(
                                    _sm, np.asarray(_tta, np.float64), _org, _cst, PW / 1000.0)
                            except Exception as _e2:
                                logger.warning(f"    ⚠ beyond-distance analysis for '{_dk}' failed "
                                               f"({type(_e2).__name__}: {_e2})", exc_info=True)
                    v910_write_smcda_tables(_dk, _sm)
                except Exception as _e:
                    logger.warning(f"    ⚠ v9.10 Table 2 / benchmarks for '{_dk}' failed "
                                   f"({type(_e).__name__}: {_e})", exc_info=True)
                for _ln in format_smcda_threshold_table(_sm, cell_km=PW / 1000.0):
                    logger.info(_ln)
            except Exception as _e:
                logger.warning(f"    ⚠ validation of domain '{_dk}' failed "
                               f"({type(_e).__name__}: {_e}); its thresholds are "
                               f"unaffected.")
        # ---- one machine-readable summary covering every domain --------------
        try:
            import json as _json
            _fp = Config.OUT_DIR / f"{V910_TAG}_smcda_thresholds.json"
            _doms = {}
            for _dk, _sm in SMCDA_BY_DOMAIN.items():
                _v = _sm.get("validation", {}) or {}
                _doms[_dk] = dict(
                    domain_label=_sm["domain_label"],
                    entry_basis=_sm["entry_basis"],
                    validated_against=_sm.get("validated_against"),
                    break_method=_sm["break_method"], break_note=_sm["break_note"],
                    consistency_ratio=(None if _sm.get("ahp") is None
                                       else float(_sm["ahp"]["CR"])),
                    ahp_consistent=(None if _sm.get("ahp") is None
                                    else bool(_sm["ahp"]["consistent"])),
                    weights={k: float(v) for k, v in _sm["weights"].items()},
                    criteria=[dict(key=c.key, label=c.label, direction=c.direction,
                                   standardise=c.standardise, layer=c.note)
                              for c in _sm["criteria"]],
                    dropped_criteria=[{"key": k, "layer": l}
                                      for k, l in _sm.get("dropped", [])],
                    class_labels=list(_sm["labels"]),
                    edges=[float(x) for x in _sm["edges"]],
                    thresholds=[dict(cls=i, label=lab,
                                     lower=float(_sm["edges"][i]),
                                     upper=float(_sm["edges"][i + 1]))
                                for i, lab in enumerate(_sm["labels"])],
                    validation={k: v for k, v in _v.items() if k != "classes"},
                    validation_by_class=_v.get("classes", []),
                    criterion_diagnostic=_sm.get("criterion_diagnostic", []),
                    sensitivity=_sm.get("sensitivity", []),
                    table2=_sm.get("table2", []),                        # v9.10
                    benchmarks=_sm.get("benchmarks", []),
                )
            _payload = dict(
                version="9.10", cell_km=float(PW / 1000.0),
                rok_split=ROK_SPLIT_INFO, exclaves=ROK_EXCLAVE_INFO,
                base_raster=Config.BASE_RASTER,
                primary_domain=str(getattr(Config, "SMCDA_PRIMARY_DOMAIN",
                                           Config.SMCDA_APPLY_TO_DOMAIN)).lower(),
                weight_method=str(Config.SMCDA_WEIGHT_METHOD),
                aggregation=str(Config.SMCDA_AGGREGATION),
                owa_andness=(float(Config.SMCDA_OWA_ANDNESS)
                             if str(Config.SMCDA_AGGREGATION).lower() == "owa" else None),
                entry_source=str(Config.SMCDA_ENTRY_SOURCE),
                ahp_matrix_source=(
                    "Config.SMCDA_AHP_MATRIX (explicit pairwise judgements)"
                    if Config.SMCDA_AHP_MATRIX is not None
                    else "CONSTRUCTED from SMCDA_PRIORITY_VECTOR — NOT an expert "
                         "elicitation; replace before operational use"),
                note=("Class boundaries are relative to the domain they were computed "
                      "over. Do NOT compare an edge value from one domain with an edge "
                      "value from another unless break_method is 'absolute'."),
                domains=_doms)
            _fp.parent.mkdir(parents=True, exist_ok=True)
            with open(_fp, "w", encoding="utf-8") as _fh:
                _json.dump(_payload, _fh, ensure_ascii=False, indent=2, default=float)
            logger.info(f"\n   ✓ SMCDA thresholds (all domains, machine-readable): {_fp}")
            if Config.DL_DIR:
                try:
                    shutil.copy2(_fp, Config.DL_DIR / _fp.name)
                except Exception:
                    pass
        except Exception as _je:
            logger.warning(f"   thresholds JSON not written ({_je}).")
    validate_smcda_domains()
    def _smcda_colors():
        _c = list(Config.SMCDA_CLASS_COLORS)[:int(Config.SMCDA_N_CLASSES)]
        while len(_c) < int(Config.SMCDA_N_CLASSES):
            _c.append("#999999")
        return _c
    def render_smcda_risk_map(dk, sm):
        """Two panels: the classified risk map on the full model grid, and the
        score histogram with the class boundaries drawn on it so the thresholds
        are visible as cuts in the actual distribution, not just numbers."""
        try:
            from matplotlib.patches import Patch as _Patch
            import textwrap as _tw
            _cols = _smcda_colors()
            _labels = list(sm["labels"]); _edges = list(sm["edges"])
            _cls = sm["classes"]
            _val = sm.get("validation", {}) or {}
            _rows = {r["cls"]: r for r in _val.get("classes", [])}
            fig = plt.figure(figsize=(17, 13), dpi=150, facecolor="white")
            gs = fig.add_gridspec(1, 2, width_ratios=[2.15, 1.0], wspace=0.14)
            ax = fig.add_subplot(gs[0, 0]); ax.set_facecolor("white")
            _cmap = ListedColormap(_cols)
            _norm = BoundaryNorm(np.arange(-0.5, len(_cols) + 0.5, 1.0),
                                 ncolors=len(_cols), clip=True)
            # the rest of the model domain in pale grey, so a restricted domain
            # reads as "not classified" rather than as "no data anywhere"
            _out = KP & ~sm["mask"]
            if _out.any():
                ax.imshow(np.ma.masked_where(~_out, np.ones(SH, float)),
                          cmap=ListedColormap(["#ededed"]), extent=(0, NC, NR, 0),
                          origin="upper", aspect="auto", interpolation="nearest",
                          zorder=1.5)
            ax.imshow(np.ma.masked_invalid(_cls), cmap=_cmap, norm=_norm,
                      extent=(0, NC, NR, 0), origin="upper", aspect="auto",
                      interpolation="nearest", zorder=2)
            if water_segments:
                ax.add_collection(LineCollection(
                    water_segments, colors=Config.WATER_LINE_COLOR,
                    linewidths=Config.WATER_LINE_WIDTH,
                    alpha=Config.WATER_LINE_ALPHA, zorder=4))
            if sig_segments:
                ax.add_collection(LineCollection(
                    sig_segments, colors="#555555",
                    linewidths=Config.SIG_LINE_WIDTH, alpha=0.7, zorder=5))
            ax.contour(MDL_BOUNDARY.astype(float), levels=[0.5], colors="#111111",
                       linewidths=1.6, zorder=6.5)
            if Config.SHOW_CORRIDOR_MARKERS and CORRIDORS:
                plot_corridor_markers(ax, CORRIDORS, label=False, size_scale=22.0)
            if Config.SHOW_STUDY_AREA_OUTLINE and study_segments:
                ax.add_collection(LineCollection(
                    study_segments, colors=Config.STUDY_AREA_LINE_COLOR,
                    linewidths=Config.STUDY_AREA_LINE_WIDTH,
                    alpha=Config.STUDY_AREA_LINE_ALPHA, zorder=6.7))
            ax.set_xlim(0, NC); ax.set_ylim(NR, 0)
            ax.set_xticks([]); ax.set_yticks([])
            _cell_km = PW / 1000.0
            _bar = 50.0 / _cell_km
            _x0, _y0 = NC * 0.62, NR * 0.94
            ax.plot([_x0, _x0 + _bar], [_y0, _y0], "k-", lw=3, zorder=8)
            ax.text(_x0 + _bar / 2, _y0 - NR * 0.012, "50 km", ha="center",
                    va="bottom", fontsize=10, fontweight="bold")
            _leg = []
            for _i, _lab in enumerate(_labels):
                _pc = _rows.get(_i, {}).get("pct")
                _leg.append(_Patch(fc=_cols[_i], ec="#3a3a3a",
                                   label=f"{_lab}  [{_edges[_i]:.3f}–{_edges[_i+1]:.3f}]"
                                         + (f"  {_pc:.0f}% of domain" if _pc is not None else "")))
            if _out.any():
                _leg.append(_Patch(fc="#ededed", ec="#3a3a3a", label="outside this domain"))
            _leg.append(Line2D([0], [0], color="#111111", lw=1.6, label="MDL boundary"))
            if CORRIDORS:
                _leg.append(Line2D([0], [0], marker="v", linestyle="none",
                                   markerfacecolor="#00e5ff", markeredgecolor="#003b46",
                                   markersize=9, label="Water-deer corridor gate"))
            ax.legend(handles=_leg, loc="lower left", fontsize=10, frameon=True,
                      facecolor="white", edgecolor="#888")
            _wtxt = "  ".join(f"{k}={sm['weights'][k]:.2f}"
                              for k in sorted(sm["weights"], key=lambda x: -sm["weights"][x]))
            _crtxt = ("" if sm.get("ahp") is None
                      else f" · CR={sm['ahp']['CR']:.3f}"
                           f"{'' if sm['ahp']['consistent'] else ' ⚠ INCONSISTENT'}")
            _sub = (f"{str(Config.SMCDA_WEIGHT_METHOD).upper()} weights · "
                    f"{str(Config.SMCDA_AGGREGATION).upper()}"
                    + (f" (andness {Config.SMCDA_OWA_ANDNESS:.2f})"
                       if str(Config.SMCDA_AGGREGATION).lower() == "owa" else "")
                    + f" · {sm['break_method']} breaks{_crtxt}")
            _ent = f"ENTRY measured from {sm.get('entry_basis')}"
            _ttl = (f"PPR risk areas — SMCDA over {sm['domain_label']}\n"
                    + "\n".join(_tw.wrap(_sub, 74)) + "\n"
                    + "\n".join(_tw.wrap(_ent, 74)) + "\n"
                    + "\n".join(_tw.wrap(_wtxt, 74)))
            ax.set_title(_ttl, fontsize=10.5, fontweight="bold")
            # ── right panel: the distribution the thresholds actually cut ────
            ax2 = fig.add_subplot(gs[0, 1]); ax2.set_facecolor("white")
            _sv = np.asarray(sm["score"], np.float64)[sm["mask"]]
            _sv = _sv[np.isfinite(_sv)]
            _n, _bins, _patches = ax2.hist(_sv, bins=60, color="#bdbdbd",
                                           edgecolor="white", linewidth=0.4)
            _ctr = 0.5 * (_bins[:-1] + _bins[1:])
            for _p, _c in zip(_patches, _ctr):
                _ci = int(np.clip(np.searchsorted(_edges[1:-1], _c, side="right"),
                                  0, len(_cols) - 1))
                _p.set_facecolor(_cols[_ci])
            for _e_ in _edges[1:-1]:
                ax2.axvline(_e_, color="#111111", lw=1.4, ls="--", zorder=5)
                ax2.text(_e_, ax2.get_ylim()[1] * 0.985, f" {_e_:.3f}", rotation=90,
                         va="top", ha="left", fontsize=9, fontweight="bold",
                         color="#111111")
            ax2.set_xlabel("SMCDA risk score (0–1)", fontsize=11)
            ax2.set_ylabel("model cells", fontsize=11)
            _mono = _val.get("arrival_monotonic")
            _auc = _val.get("auc_reached")
            _rho = _val.get("spearman_score_vs_arrival")
            _bits = []
            if _mono is not None:
                _bits.append("arrival order " + ("CORRECT" if _mono else "NOT monotonic"))
            if _auc is not None:
                _bits.append(f"AUC={_auc:.3f}")
            if _rho is not None:
                _bits.append(f"rho={_rho:+.3f}")
            _t2 = "Score distribution and class thresholds"
            if _bits:
                _t2 += "\n" + " · ".join(_bits)
            if _mono is False:
                _bad = [d for d in (sm.get("criterion_diagnostic") or [])
                        if not d["sign_ok"]]
                if _bad:
                    _t2 += ("\n⚠ wrong-signed: "
                            + ", ".join(f"{d['criterion']} ({d['rho_vs_arrival']:+.2f})"
                                        for d in _bad[:3]))
            ax2.set_title(_t2, fontsize=11, fontweight="bold")
            ax2.grid(axis="y", alpha=0.25, lw=0.6)
            for _s in ("top", "right"):
                ax2.spines[_s].set_visible(False)
            _fp = Config.OUT_DIR / f"PPR_v9.8_smcda_risk_areas_{dk}.png"
            fig.savefig(_fp, dpi=150, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            logger.info(f"   SMCDA risk-area map [{dk}]: {_fp}")
            if Config.DL_DIR:
                try:
                    shutil.copy2(_fp, Config.DL_DIR / _fp.name)
                except Exception:
                    pass
        except Exception as _e:
            logger.warning(f"   SMCDA risk-area map [{dk}] skipped: "
                           f"{type(_e).__name__}: {_e}")
    def render_smcda_korea_panel(dk, sm):
        """The risk classes drawn on the SAME white-canvas peninsula panel as the
        NK -> corridor gates -> SK progression figure: county lines, waterways,
        the MDL, the gates and the China introduction seeds, cropped to Korea.

        Beside it, the epidemic's own arrival-month map on the identical frame —
        so the classification and the thing it is meant to predict can be read
        against each other without switching between two figures at two scales."""
        if not Config.SHOW_SMCDA_KOREA_PANEL:
            return
        if not SMCDA_KOREA_EPI or SMCDA_KOREA_EPI.get("sl") is None:
            logger.info(f"   SMCDA Korea panel [{dk}] skipped: the Korea leg did not run.")
            return
        try:
            from matplotlib.patches import Patch as _Patch
            import textwrap as _tw
            _sl = SMCDA_KOREA_EPI["sl"]
            _r0 = SMCDA_KOREA_EPI["r0"]; _c0 = SMCDA_KOREA_EPI["c0"]
            KP_k = SMCDA_KOREA_EPI["KP_k"]; NK_k = SMCDA_KOREA_EPI["NK_k"]
            SK_k = SMCDA_KOREA_EPI["SK_k"]
            _sig_seg = list(SMCDA_KOREA_EPI.get("sig_seg_k") or [])
            _water_seg = list(SMCDA_KOREA_EPI.get("water_seg_k") or [])
            _seeds_k = list(SMCDA_KOREA_EPI.get("seeds") or [])
            _cls_k = np.asarray(sm["classes"], np.float32)[_sl]
            _msk_k = np.asarray(sm["mask"], bool)[_sl]
            _cls_k = np.where(_msk_k, _cls_k, np.nan)
            if not np.any(np.isfinite(_cls_k)):
                logger.info(f"   SMCDA Korea panel [{dk}] skipped: this domain has no "
                            f"cells inside the Korea crop.")
                return

            # ── v9.8: crop the frame to South Korea for an SK-only domain ────
            # An SK-only classification drawn on the whole-peninsula frame
            # spends half the figure on grey North Korean land and shrinks the
            # part that carries the decision. The crop below is a pure change
            # of VIEW: nothing is reclassified, no threshold is recomputed, no
            # cell changes class. It composes a second slice onto the Korea-leg
            # crop and shifts everything that is stored in crop coordinates
            # (the county and waterway segments, the seed markers) by the same
            # offset, so the two panels stay on one identical frame.
            _crop_mode = str(getattr(Config, "SMCDA_KOREA_PANEL_CROP", "auto")).lower()
            _nk_in_domain = bool((NK_k & _msk_k).any())
            _do_sk_crop = (_crop_mode == "sk" or
                           (_crop_mode == "auto" and not _nk_in_domain))
            _crop_note = ""
            if _do_sk_crop:
                # Keep the whole of South Korea in view even when only part of
                # it is classified (the _reached footprint), so the reader can
                # see how far up the country the footprint actually goes.
                _keep = SK_k | (_msk_k & np.isfinite(_cls_k))
                if not _keep.any():
                    _do_sk_crop = False
            if _do_sk_crop:
                _mrg = max(0, int(getattr(Config, "SMCDA_KOREA_PANEL_CROP_MARGIN", 3)))
                _rr = np.where(_keep.any(axis=1))[0]
                _cc = np.where(_keep.any(axis=0))[0]
                _ra = max(0, int(_rr.min()) - _mrg)
                _rb = min(KP_k.shape[0], int(_rr.max()) + 1 + _mrg)
                _ca = max(0, int(_cc.min()) - _mrg)
                _cb = min(KP_k.shape[1], int(_cc.max()) + 1 + _mrg)
                # Pull the MDL and the corridor gates into frame along the top
                # edge: they are the introduction route this map is about, and
                # a map of gate-distance risk that does not show the gates is
                # harder to read than one that does.
                _mdl_pen = MDL_BOUNDARY[_sl]
                if _mdl_pen.any():
                    _mr = np.where(_mdl_pen.any(axis=1))[0]
                    _ra = max(0, min(_ra, int(_mr.min()) - 1))
                if (_rb - _ra) >= 8 and (_cb - _ca) >= 8:
                    _sub_r = slice(_ra, _rb); _sub_c = slice(_ca, _cb)
                    KP_k = KP_k[_sub_r, _sub_c]; NK_k = NK_k[_sub_r, _sub_c]
                    SK_k = SK_k[_sub_r, _sub_c]
                    _cls_k = _cls_k[_sub_r, _sub_c]; _msk_k = _msk_k[_sub_r, _sub_c]
                    _sl = (slice(_sl[0].start + _ra, _sl[0].start + _rb),
                           slice(_sl[1].start + _ca, _sl[1].start + _cb))
                    _r0 = int(_r0) + _ra; _c0 = int(_c0) + _ca
                    _sig_seg = _shift_segments(_sig_seg, -_ca, -_ra)
                    _water_seg = _shift_segments(_water_seg, -_ca, -_ra)
                    _seeds_k = [(s[0] - _ra, s[1] - _ca) for s in _seeds_k
                                if _ra <= s[0] < _rb and _ca <= s[1] < _cb]
                    _crop_note = ("frame cropped to South Korea (North Korea is "
                                  "outside this domain and outside this view)")
                    logger.info(f"   SMCDA Korea panel [{dk}]: frame cropped to South "
                                f"Korea, rows {_ra}:{_rb}, cols {_ca}:{_cb} of the "
                                f"peninsula crop.")
                else:
                    _do_sk_crop = False
                    logger.info(f"   SMCDA Korea panel [{dk}]: South Korea crop refused "
                                f"— the resulting frame would be {_rb-_ra}x{_cb-_ca} "
                                f"cells, too small to draw. Full peninsula used.")
            if not np.any(np.isfinite(_cls_k)):
                logger.info(f"   SMCDA Korea panel [{dk}] skipped: no classified cells "
                            f"survived the frame crop.")
                return
            _cols = _smcda_colors()
            _labels = list(sm["labels"]); _edges = list(sm["edges"])
            _val = sm.get("validation", {}) or {}
            _rows = {r["cls"]: r for r in _val.get("classes", [])}
            _ext = (0, KP_k.shape[1], KP_k.shape[0], 0)
            _ones = np.ones(KP_k.shape, float)
            _pub = bool(getattr(Config, "PUB_FIGS", False))
            if _pub:
                fig, axes = plt.subplots(2, 1, figsize=(10, 13), dpi=150, facecolor="white")
            else:
                fig, axes = plt.subplots(1, 2, figsize=(15, 12), dpi=150, facecolor="white")

            def _frame(ax):
                ax.set_facecolor("white")
                # the peninsula outside this domain: pale grey land, so an
                # SK-only classification does not look like missing data in NK
                _rest = KP_k & ~_msk_k
                if _rest.any():
                    ax.imshow(np.ma.masked_where(~_rest, _ones),
                              cmap=ListedColormap(["#ededed"]), vmin=0, vmax=1,
                              extent=_ext, origin="upper", aspect="equal",
                              interpolation="nearest", zorder=1.0)
                if Config.SK_PANEL_SHOW_COUNTIES and _sig_seg:
                    ax.add_collection(LineCollection(
                        _sig_seg,
                        colors=Config.SK_PANEL_COUNTY_LINE_COLOR,
                        linewidths=Config.SK_PANEL_COUNTY_LINE_WIDTH,
                        alpha=Config.SK_PANEL_COUNTY_LINE_ALPHA, zorder=4.1))
                if Config.SK_PANEL_SHOW_WATER and _water_seg:
                    ax.add_collection(LineCollection(
                        _water_seg,
                        colors=Config.SK_PANEL_WATER_LINE_COLOR,
                        linewidths=Config.SK_PANEL_WATER_LINE_WIDTH,
                        alpha=Config.SK_PANEL_WATER_LINE_ALPHA, zorder=4.3))
                if Config.SK_PANEL_OUTLINE and SK_k.any():
                    ax.contour(SK_k.astype(float), levels=[0.5],
                               colors=Config.SK_PANEL_OUTLINE_COLOR,
                               linewidths=Config.SK_PANEL_OUTLINE_WIDTH, zorder=4.0)
                _mdl_c = MDL_BOUNDARY[_sl]
                if _mdl_c.any():
                    ax.imshow(np.ma.masked_where(~_mdl_c, np.ones_like(_mdl_c, float)),
                              cmap=ListedColormap(["#d81b60"]), extent=_ext,
                              origin="upper", aspect="equal", zorder=5)
                if Config.SHOW_CORRIDOR_MARKERS and CORRIDORS:
                    plot_corridor_markers(ax, CORRIDORS, row_offset=_r0, col_offset=_c0,
                                          label=False, size_scale=18.0)
                _sd = _seeds_k
                if _sd:
                    ax.scatter([s[1] + 0.5 for s in _sd], [s[0] + 0.5 for s in _sd],
                               s=18, marker="v", c="#7a3b00", edgecolors="white",
                               linewidths=0.4, zorder=6)
                # On the South Korea crop only a sliver of NK land is left in
                # frame above the MDL; labelling it there collides with the
                # frame edge and misleadingly implies NK is part of the map.
                if NK_k.any() and int(NK_k.sum()) >= max(40, int(0.02 * NK_k.size)):
                    ax.text(KP_k.shape[1] * 0.5,
                            np.where(NK_k.any(axis=1))[0].min() + 1.5, "NORTH KOREA",
                            ha="center", fontsize=8, fontweight="bold", color="#7a3b00")
                if SK_k.any():
                    ax.text(KP_k.shape[1] * 0.5, KP_k.shape[0] - 2, "SOUTH KOREA",
                            ha="center", fontsize=8, fontweight="bold", color="#0b3d6b")
                ax.set_xticks([]); ax.set_yticks([])
                # ── v9.8: LOCK THE AXES TO THE FRAME ─────────────────────────
                # add_collection() autoscales. The county and waterway outlines
                # extend beyond the cropped frame (they are clipped, which is
                # what we want), so without this the axes silently expand back
                # out to cover them and the crop has no visible effect at all —
                # the figure renders at the old peninsula size with the map
                # squeezed into a corner. imshow's extent is the frame; nothing
                # drawn afterwards may change it.
                ax.set_xlim(_ext[0], _ext[1])
                ax.set_ylim(_ext[2], _ext[3])
                if _pub and _msk_k.any():
                    # 11 Nov: crop to the classified domain and the fence
                    _rr, _cc = np.where(_msk_k | MDL_BOUNDARY[_sl])
                    _pd = 4
                    ax.set_xlim(max(0, _cc.min() - _pd), min(KP_k.shape[1], _cc.max() + 1 + _pd))
                    ax.set_ylim(min(KP_k.shape[0], _rr.max() + 1 + _pd), max(0, _rr.min() - _pd))
                    _xl, _yl = ax.get_xlim(), ax.get_ylim()
                    _bar = 25.0 / (PW / 1000.0)
                    _x0 = _xl[0] + 0.05 * (_xl[1] - _xl[0]); _y0 = _yl[0] - 0.05 * (_yl[0] - _yl[1])
                    ax.plot([_x0, _x0 + _bar], [_y0, _y0], "k-", lw=3, zorder=9)
                    ax.text(_x0 + _bar / 2, _y0 - 0.015 * (_yl[0] - _yl[1]), "25 km", ha="center",
                            va="bottom", fontsize=10, fontweight="bold", zorder=9)
                ax.set_autoscale_on(False)

            # ── left: the risk classes ───────────────────────────────────────
            axA = axes[0]; _frame(axA)
            axA.imshow(np.ma.masked_invalid(_cls_k),
                       cmap=ListedColormap(_cols),
                       norm=BoundaryNorm(np.arange(-0.5, len(_cols) + 0.5, 1.0),
                                         ncolors=len(_cols), clip=True),
                       extent=_ext, origin="upper", aspect="equal",
                       interpolation="nearest", alpha=0.92, zorder=3.0)
            _leg = []
            for _i, _lab in enumerate(_labels):
                _pc = _rows.get(_i, {}).get("pct")
                _leg.append(_Patch(fc=_cols[_i], ec="#3a3a3a",
                                   label=f"{_lab}  [{_edges[_i]:.3f}–{_edges[_i+1]:.3f}]"
                                         + (f"  {_pc:.0f}%" if _pc is not None else "")))
            _restricted = str(dk).lower() in ("sk_reached", "korea_reached")
            if (KP_k & ~_msk_k).any():
                _leg.append(_Patch(fc="#ededed", ec="#3a3a3a",
                                   label=("not reached by the cutoff month"
                                          if _restricted else "outside this domain")))
            if _restricted:
                _leg.append(Line2D([0], [0], color="#111111", lw=1.3,
                                   label="edge of the invaded footprint"))
            _leg += [Line2D([0], [0], color="#d81b60", lw=2.0, label="MDL"),
                     Line2D([0], [0], marker="v", linestyle="none",
                            markerfacecolor="#00e5ff", markeredgecolor="#003b46",
                            markersize=8, label="Water-deer corridor gate"),
                     Line2D([0], [0], marker="v", linestyle="none",
                            markerfacecolor="#7a3b00", markeredgecolor="white",
                            markersize=8, label="China introduction seed")]
            # v9.8: the South Korea crop makes the frame tall and narrow, so an
            # in-axes legend covers the southern third of the country — exactly
            # the part a surveillance design has to read. Put it under the map.
            if _pub:
                axA.legend(handles=_leg, loc="center left", bbox_to_anchor=(1.02, 0.5),
                           fontsize=9.5, frameon=True, facecolor="white", edgecolor="#888")
            else:
              axA.legend(handles=_leg, loc="upper center", bbox_to_anchor=(0.5, -0.01),
                       ncol=2, fontsize=8.5, frameon=True,
                       facecolor="white", edgecolor="#888")
            _wtxt = "  ".join(f"{k}={sm['weights'][k]:.2f}"
                              for k in sorted(sm["weights"], key=lambda x: -sm["weights"][x]))
            axA.set_title(
                "\n".join(_tw.wrap(f"PPR risk areas — SMCDA, {sm['domain_label']}", 62))
                + "\n"
                + "\n".join(_tw.wrap(f"ENTRY from {sm.get('entry_basis')} · "
                                     f"{sm['break_method']} breaks", 62)) + "\n"
                + "\n".join(_tw.wrap(_wtxt, 62))
                + (("\n" + "\n".join(_tw.wrap(_crop_note, 62))) if _crop_note else ""),
                fontsize=10, fontweight="bold")
            # ── right: what the epidemic actually did, same frame ────────────
            axB = axes[1]; _frame(axB)
            _tta_k = np.asarray(SMCDA_KOREA_EPI["time_to_arrival"], np.float64)[_sl]
            _off = int(SMCDA_KOREA_EPI.get("month_offset", 0))
            _arr = np.where(np.isfinite(_tta_k) & KP_k, _tta_k + 1.0 + _off, np.nan)
            if np.any(np.isfinite(_arr)):
                _im = axB.imshow(np.ma.masked_invalid(_arr), cmap="viridis_r",
                                 extent=_ext, origin="upper", aspect="equal",
                                 interpolation="nearest", alpha=0.92, zorder=3.0)
                _cb = fig.colorbar(_im, ax=axB, fraction=0.040, pad=0.02)
                _cb.set_label("Arrival month (simulation timeline)", fontsize=10)
            _never = KP_k & ~np.isfinite(_tta_k)
            if _never.any():
                axB.imshow(np.ma.masked_where(~_never, _ones),
                           cmap=ListedColormap(["#f2f2f2"]), vmin=0, vmax=1,
                           extent=_ext, origin="upper", aspect="equal",
                           interpolation="nearest", zorder=2.9)
            # On a footprint domain, draw its edge on BOTH panels so the reader
            # can see that the left map covers exactly the invaded area and no
            # more — the single most likely thing to be misread about this figure.
            if str(dk).lower() in ("sk_reached", "korea_reached"):
                for _a in (axA, axB):
                    _a.contour(_msk_k.astype(float), levels=[0.5], colors="#111111",
                               linewidths=1.3, zorder=7.0)
            _bits = []
            if _val.get("arrival_monotonic") is not None:
                _bits.append("arrival order "
                             + ("CORRECT" if _val["arrival_monotonic"] else "NOT monotonic"))
            if _val.get("auc_reached") is not None:
                _bits.append(f"AUC={_val['auc_reached']:.3f}")
            if _val.get("spearman_score_vs_arrival") is not None:
                _bits.append(f"rho={_val['spearman_score_vs_arrival']:+.3f}")
            axB.set_title("What the epidemic did — modelled arrival month\n"
                          + ("Korea leg: China border → North Korea → corridor gates "
                             "→ South Korea\n")
                          + (" · ".join(_bits) if _bits else "not validated"),
                          fontsize=10, fontweight="bold")
            # v9.8: the left title runs to six wrapped lines on a footprint
            # domain; at y=0.995 the suptitle sat on top of it. Give the panel
            # titles their own band.
            if _pub:
                for _a, _l in ((axA, "A"), (axB, "B")):
                    _a.text(-0.03, 1.02, _l, transform=_a.transAxes, fontsize=18, fontweight="bold",
                            va="bottom", ha="left")
            fig.suptitle("Risk classification (left) vs the mechanistic epidemic (right) — "
                         "the left map never saw the right one",
                         fontsize=12.5, fontweight="bold", y=1.02)
            fig.subplots_adjust(top=0.90)
            _fp = Config.OUT_DIR / f"PPR_v9.8_smcda_korea_panel_{dk}.png"
            fig.savefig(_fp, dpi=150, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            logger.info(f"   SMCDA Korea panel [{dk}]: {_fp}")
            if Config.DL_DIR:
                try:
                    shutil.copy2(_fp, Config.DL_DIR / _fp.name)
                except Exception:
                    pass
        except Exception as _e:
            logger.warning(f"   SMCDA Korea panel [{dk}] skipped: "
                           f"{type(_e).__name__}: {_e}")
    if Config.SHOW_SMCDA_MAP:
        for _dk, _sm in SMCDA_BY_DOMAIN.items():
            render_smcda_risk_map(_dk, _sm)
            if str(_dk).lower() in ("sk", "korea", "nk", "sk_reached", "korea_reached"):
                render_smcda_korea_panel(_dk, _sm)
    def export_smcda_county_csv(dk, sm):
        """Per-sigungu risk levels: the area-weighted mean score, the share of
        the county in each risk class, and the class the county is ASSIGNED —
        by dominant area, which is the defensible choice for a county-level
        control measure, with the maximum class reported beside it because a
        small very-high-risk pocket is what a surveillance team actually visits."""
        if not SK_COUNTY_ZONES:
            logger.info("\n[SMCDA COUNTY TABLE] skipped: no county zones were built "
                        "(sig.shp not found). The raster thresholds are unaffected.")
            return
        try:
            logger.info(f"\n[SMCDA COUNTY RISK TABLE — domain '{dk}']")
            _labels = list(sm["labels"]); _K = len(_labels)
            _cls = np.nan_to_num(sm["classes"], nan=-1.0)
            _score = np.asarray(sm["score"], np.float64)
            _m = np.asarray(sm["mask"], bool)
            _rows = []
            for _z in SK_COUNTY_ZONES:
                _r = np.asarray(_z["rows"], np.int64)
                _c = np.asarray(_z["cols"], np.int64)
                _w = np.asarray(_z["weights"], np.float64)
                _ok = (_r >= 0) & (_r < NR) & (_c >= 0) & (_c < NC)
                _r, _c, _w = _r[_ok], _c[_ok], _w[_ok]
                # How much of the county the classification actually covers. On a
                # restricted domain (the invaded footprint) a county can be 5%
                # inside it, and a mean score computed over that 5% must not be
                # read as a statement about the county.
                _w_all = float(_w.sum()) if _w.size else 0.0
                if _r.size:
                    _keep = _m[_r, _c]
                    _r, _c, _w = _r[_keep], _c[_keep], _w[_keep]
                _rec = dict(sig_cd=_z.get("code", ""), county_kr=_z.get("name_kr", ""),
                            county_en=_z.get("name_en", ""), province=_z.get("province", ""),
                            sampling=_z.get("sampling", ""),
                            poly_area_km2=_z.get("poly_area_km2", float("nan")),
                            centroid_lat=_z.get("lat", float("nan")),
                            centroid_lon=_z.get("lon", float("nan")))
                _wsum = float(_w.sum()) if _w.size else 0.0
                _rec["county_coverage_pct"] = (
                    100.0 * _wsum / _w_all if _w_all > 0 else float("nan"))
                if _r.size == 0 or _wsum <= 0:
                    _rec.update(status="outside the classified domain", model_cells=0,
                                mean_score=float("nan"), max_score=float("nan"),
                                risk_class=None, risk_level="", max_risk_level="")
                    for _i, _lab in enumerate(_labels):
                        _rec[f"pct_{_i}_{_lab.lower().replace(' ', '_')}"] = float("nan")
                    _rows.append(_rec)
                    continue
                _sv = _score[_r, _c]
                _cv = _cls[_r, _c].astype(int)
                _sh = np.zeros(_K, np.float64)
                for _i in range(_K):
                    _sh[_i] = float(_w[_cv == _i].sum())
                _sh = 100.0 * _sh / _wsum
                _dom = int(np.argmax(_sh)); _mx = int(_cv.max())
                _rec.update(status="classified", model_cells=int(_r.size),
                            mean_score=float(np.average(_sv, weights=_w)),
                            max_score=float(np.max(_sv)),
                            risk_class=_dom, risk_level=_labels[_dom],
                            max_risk_level=_labels[_mx])
                for _i, _lab in enumerate(_labels):
                    _rec[f"pct_{_i}_{_lab.lower().replace(' ', '_')}"] = float(_sh[_i])
                _rows.append(_rec)
            _ranked = [r for r in _rows if r.get("status") == "classified"]
            _ranked.sort(key=lambda d: -d["mean_score"])
            for _i, _r in enumerate(_ranked, 1):
                _r["rank_risk"] = _i
            for _r in _rows:
                _r.setdefault("rank_risk", None)
            _rows.sort(key=lambda d: (d["rank_risk"] is None, d.get("rank_risk") or 0))
            _cols_out = (["rank_risk", "sig_cd", "county_kr", "county_en", "province",
                          "status", "risk_level", "risk_class", "max_risk_level",
                          "mean_score", "max_score"]
                         + [f"pct_{_i}_{_l.lower().replace(' ', '_')}"
                            for _i, _l in enumerate(_labels)]
                         + ["county_coverage_pct", "model_cells", "poly_area_km2",
                            "sampling", "centroid_lat", "centroid_lon"])
            _v = sm.get("validation", {}) or {}
            _mono = _v.get("arrival_monotonic")
            _meta = [
                f"PPR v9.8 — SMCDA risk-area classification by South-Korean county (sigungu).",
                f"Classification domain: {sm['domain_label']}",
                "Class boundaries are RELATIVE TO THIS DOMAIN. Do not compare them with "
                "the boundaries of another domain unless break_method is 'absolute'.",
                f"ENTRY criterion measured from: {sm.get('entry_basis')}",
                f"Weights: {Config.SMCDA_WEIGHT_METHOD}"
                + ("" if sm.get("ahp") is None else f" (AHP CR={sm['ahp']['CR']:.4f}, "
                   f"{'consistent' if sm['ahp']['consistent'] else 'INCONSISTENT'})"),
                "AHP matrix source: " + ("Config.SMCDA_AHP_MATRIX (explicit pairwise "
                                         "judgements)" if Config.SMCDA_AHP_MATRIX is not None
                                         else "CONSTRUCTED from SMCDA_PRIORITY_VECTOR — NOT an "
                                              "expert elicitation; replace before operational use"),
                "Weight values: " + ", ".join(f"{k}={v:.4f}" for k, v in sm["weights"].items()),
                "Criteria used: " + ", ".join(f"{c.key}({c.direction},{c.standardise})"
                                              for c in sm["criteria"]),
                ("Criteria DROPPED (layer unavailable): "
                 + ", ".join(f"{k}<-{l}" for k, l in sm.get("dropped", []))
                 ) if sm.get("dropped") else "Criteria dropped: none",
                f"Aggregation: {Config.SMCDA_AGGREGATION}"
                + (f" (andness={Config.SMCDA_OWA_ANDNESS})"
                   if str(Config.SMCDA_AGGREGATION).lower() == "owa" else ""),
                f"Break method: {sm['break_method']} — {sm['break_note']}",
                "Class thresholds: " + "; ".join(
                    f"{_l} [{sm['edges'][_i]:.4f}, {sm['edges'][_i+1]:.4f}]"
                    for _i, _l in enumerate(_labels)),
                f"Validated against: {sm.get('validated_against', 'not validated')}",
                ("VALIDATION: arrival month is monotonic across classes — the "
                 "classification ranks the modelled epidemic correctly."
                 if _mono is True else
                 ("⚠ VALIDATION FAILED: arrival month is NOT monotonic across classes. "
                  "These class boundaries do not rank the modelled epidemic over this "
                  "domain — do not allocate surveillance from them until the cause is "
                  "resolved."
                  if _mono is False else
                  "VALIDATION: not available for this domain.")),
                # An undefined AUC must not be written as "nan": on a footprint
                # domain every cell is reached by construction, so there is no
                # second class to discriminate against. Say that.
                ((f"AUC (score vs reached) = {_v['auc_reached']:.3f}"
                  if np.isfinite(_v.get("auc_reached", float("nan")))
                  else ("AUC: NOT APPLICABLE — every cell in this domain is 'reached' "
                        "by construction, so there is nothing to discriminate against. "
                        "The arrival-order verdict and the Spearman rho are the tests "
                        "that apply here."))
                 if _v.get("auc_reached") is not None else "AUC: not available"),
                (f"Spearman rho (score vs arrival month) = "
                 f"{_v['spearman_score_vs_arrival']:+.3f} (negative is correct)"
                 if _v.get("spearman_score_vs_arrival") is not None
                 else "Spearman rho: not available"),
                "risk_level = the class covering the LARGEST share of the county's area. "
                "max_risk_level = the highest class present anywhere in it. They differ "
                "wherever a county holds a small high-risk pocket; the pct_* columns show "
                "how much.",
                "county_coverage_pct = how much of the county lies inside the "
                "classified domain. On a footprint domain this is often well under 100%: "
                "a county at 8% has had its mean_score computed from 8% of its area, and "
                "the row is a statement about that fragment, NOT about the county. Filter "
                "on this column before ranking. Counties at 0% appear with status "
                "'outside the classified domain'.",
                "pct_0..pct_N are shares of the CLASSIFIED part of the county and sum to "
                "100% of it, not of the county. Multiply by county_coverage_pct "
                "to get shares of the whole county.",
                f"Model cell size: {PW/1000.0:.2f} km. Counties smaller than a cell are "
                f"centroid-sampled and flagged in `sampling`.",
                f"Base raster: {Config.BASE_RASTER}.",
            ]
            _fp = Config.OUT_DIR / f"PPR_v9.8_smcda_county_risk_levels_{dk}.csv"
            _fp.parent.mkdir(parents=True, exist_ok=True)
            with open(_fp, "w", newline="", encoding="utf-8-sig") as _fh:
                for _ln in _meta:
                    _fh.write(f"# {_ln}\n")
                _wtr = csv.writer(_fh)
                _wtr.writerow(_cols_out)
                # Scores and class shares need more than the 3 decimals _fmt gives:
                # the counties are RANKED on mean_score, and at 3 dp neighbouring
                # counties tie and the ranking becomes arbitrary.
                _hi = ({"mean_score", "max_score", "county_coverage_pct"}
                       | {_c for _c in _cols_out if _c.startswith("pct_")})
                for _r in _rows:
                    _wtr.writerow([_fmt(_r.get(_k), 6) if _k in _hi else _fmt(_r.get(_k))
                                   for _k in _cols_out])
            logger.info(f"   ✓ SMCDA county risk CSV [{dk}]: {_fp}  ({len(_rows)} counties)")
            if Config.DL_DIR:
                try:
                    shutil.copy2(_fp, Config.DL_DIR / _fp.name)
                except Exception:
                    pass
            _cnt = {}
            for _r in _rows:
                if _r.get("status") == "classified":
                    _cnt[_r["risk_level"]] = _cnt.get(_r["risk_level"], 0) + 1
            _n_out = sum(1 for _r in _rows if _r.get("status") != "classified")
            _n_part = sum(1 for _r in _rows if _r.get("status") == "classified"
                          and _r.get("county_coverage_pct", 100.0) < 50.0)
            logger.info(f"  {len(_ranked)} counties intersect this domain; {_n_out} lie "
                        f"entirely outside it.")
            if _n_part:
                logger.info(f"  {_n_part} of them are less than half covered — their "
                            f"mean_score describes only the covered fragment. Filter on "
                            f"county_coverage_pct before ranking.")
            logger.info("  Counties by assigned (dominant-area) risk level:")
            for _lab in _labels:
                logger.info(f"    {_lab:<12s} {_cnt.get(_lab, 0):>4d} counties")
            logger.info("  Highest-risk counties by mean score:")
            for _r in _ranked[:20]:
                _cov = _r.get("county_coverage_pct")
                logger.info(f"    {_r['rank_risk']:>3d}. "
                            f"{str(_r['county_kr'] or _r['county_en'])[:18]:<18s} "
                            f"({_r['sig_cd']:<6s} {str(_r['province'])[:10]:<10s}) "
                            f"score={_r['mean_score']:.4f}  {_r['risk_level']:<10s} "
                            f"(max here: {_r['max_risk_level']}"
                            + (f"; {_cov:.0f}% of county covered"
                               if _cov is not None and np.isfinite(_cov) and _cov < 99.5
                               else "") + ")")
            if _mono is False:
                logger.warning("  ⚠ This table's class boundaries FAILED the arrival-order "
                               "validation over this domain. The county ranking by "
                               "mean_score is still a ranking, but the CLASS LABELS should "
                               "not be treated as risk tiers until the cause is resolved "
                               "(see the per-criterion diagnostic above).")
        except Exception as _e:
            logger.warning(f"   SMCDA county CSV [{dk}] skipped: "
                           f"{type(_e).__name__}: {_e}")
    if Config.SMCDA_EXPORT_COUNTY_CSV and SMCDA_BY_DOMAIN:
        _prim_dk = str(getattr(Config, "SMCDA_PRIMARY_DOMAIN",
                               Config.SMCDA_APPLY_TO_DOMAIN)).lower()
        for _dk, _sm in SMCDA_BY_DOMAIN.items():
            # Counties are South-Korean, so a North-Korea-only classification has
            # nothing to say about them. The continental domain does, but its
            # classes are East-Asia-relative — exported only when it is primary.
            if str(_dk).lower() == "nk":
                continue
            if str(_dk).lower() == "korea_reached" and _prim_dk != "korea_reached":
                # The peninsula footprint spans both Koreas; its classes are not
                # SK-relative, so a South-Korean county table built from them
                # would carry a scale set partly by North Korea.
                logger.info("\n[SMCDA COUNTY TABLE] domain 'korea_reached' not exported: "
                            "its class scale is set across both Koreas. Use 'sk_reached' "
                            "for a South-Korean county table.")
                continue
            if str(_dk).lower() == "all" and _prim_dk != "all":
                logger.info("\n[SMCDA COUNTY TABLE] domain 'all' not exported: its classes "
                            "are East-Asia-relative, so 'High' there is not a statement "
                            "about South Korea. Set SMCDA_PRIMARY_DOMAIN='all' if you "
                            "want it anyway.")
                continue
            export_smcda_county_csv(_dk, _sm)
    # ══ VACCINATION MAPS (coverage + resulting attack in the target region) ══
    # ══ VACCINATION MAPS (coverage + resulting attack in the target region) ══
    # ══ VACCINATION MAPS (coverage + resulting attack in the target region) ══
    # ══ VACCINATION MAPS (coverage + resulting attack in the target region) ══
    def render_vaccination_maps():
        if not (Config.ENABLE_VACCINATION and Config.SHOW_VACCINATION_MAP):
            return
        res = rad_results_holder.get("res")
        if res is None or res.get("final_V") is None:
            return
        try:
            fV = np.nan_to_num(res["final_V"]) * 100.0            # coverage %
            fcum = np.nan_to_num(res["monthly_cum"][-1]) * 100.0  # attack %
            focus = KOREA_MASK if (Config.VACC_TARGET_REGION in ("korea", "nk", "sk")
                                   and KOREA_MASK.any()) else KP
            rr, cc = np.where(focus); pad = 4
            r0, r1 = max(0, rr.min() - pad), min(NR, rr.max() + 1 + pad)
            c0, c1 = max(0, cc.min() - pad), min(NC, cc.max() + 1 + pad)
            ext = (c0, c1, r1, r0)
            fmask = focus[r0:r1, c0:c1]
            band_cmap = ListedColormap(["#ffffcc", "#fed976", "#fd8d3c", "#e31a1c", "#800026"])
            band_norm = BoundaryNorm(list(Config.SPREAD_PCT_BANDS),
                                     ncolors=len(Config.SPREAD_PCT_BANDS) - 1, clip=True)
            fig, ax = plt.subplots(1, 2, figsize=(19, 9), facecolor="white")
            for a in ax:
                a.set_facecolor("white")
                a.imshow(np.ma.masked_where(~fmask, np.ones_like(fmask, float)),
                         cmap=ListedColormap(["#ece9e3"]), extent=ext, origin="upper",
                         aspect="equal", zorder=1)
                mdl_crop = MDL_BOUNDARY[r0:r1, c0:c1]
                if mdl_crop.any():
                    a.imshow(np.ma.masked_where(~mdl_crop, np.ones_like(mdl_crop, float)),
                             cmap=ListedColormap(["#d81b60"]), extent=ext, origin="upper",
                             aspect="equal", zorder=5)
                a.set_xticks([]); a.set_yticks([]); a.set_xlim(c0, c1); a.set_ylim(r1, r0)
            vac = np.where(focus, fV, np.nan)[r0:r1, c0:c1]
            im0 = ax[0].imshow(np.ma.masked_invalid(vac), cmap="Greens", vmin=0, vmax=100,
                               extent=ext, origin="upper", aspect="equal", alpha=0.9, zorder=3)
            ax[0].set_title("Vaccination coverage at end of run (%)", fontsize=12, fontweight="bold")
            fig.colorbar(im0, ax=ax[0], shrink=0.7, pad=0.02).set_label("% vaccinated", fontsize=10)
            att = np.where(focus, fcum, np.nan)[r0:r1, c0:c1]
            ax[1].imshow(np.ma.masked_invalid(att), cmap=band_cmap, norm=band_norm, extent=ext,
                         origin="upper", aspect="equal", alpha=0.9, zorder=3)
            ax[1].set_title("Resulting cumulative attack rate (%)", fontsize=12, fontweight="bold")
            sm = plt.cm.ScalarMappable(cmap=band_cmap, norm=band_norm); sm.set_array([])
            fig.colorbar(sm, ax=ax[1], boundaries=list(Config.SPREAD_PCT_BANDS),
                         ticks=list(Config.SPREAD_PCT_BANDS), shrink=0.7, pad=0.02
                         ).set_label("% attack", fontsize=10)
            _mean_cov = float(np.nanmean(fV[focus])) if focus.any() else 0.0
            fig.suptitle(
                f"PPR vaccination scenario — mode={Config.VACC_MODE}, region={Config.VACC_TARGET_REGION}, "
                f"start month {Config.VACC_START_MONTH + 1}, efficacy {Config.VACC_EFFICACY:.0%}  "
                f"(mean coverage in region: {_mean_cov:.1f}%)",
                fontsize=13, fontweight="bold")
            fp_v = Config.OUT_DIR / "PPR_v9.3_vaccination.png"
            fig.savefig(fp_v, dpi=130, facecolor="white", bbox_inches="tight"); plt.close(fig)
            logger.info(f"   vaccination map: {fp_v}")
            if Config.DL_DIR:
                try:
                    shutil.copy2(fp_v, Config.DL_DIR / fp_v.name)
                except Exception:
                    pass
        except Exception as e:
            logger.warning(f"   vaccination map skipped: {e}")
    render_vaccination_maps()
    # ══ v9.10: REVIEW RE-ANALYSIS OUTPUTS ════════════════════════════════════════
    # Tables 1 and 3, the fence-wait table, Figure 9, GeoTIFFs, the peninsular-leg
    # context file, then (optionally) the scenario / sensitivity batch. Each step is
    # isolated, so one failure never costs the others.
    V910_VACC_LABELS = {
        "baseline": "No vaccination",
        "V1_combined_host_peninsula": "V1: combined host, whole peninsula",
        "V2_ROK_livestock_on_detection": "V2: small ruminants in the ROK, after the first detection",
    }

    def v910_probe_km():
        """Zone radius for the fence-wait table and the first-crossing label."""
        if Config.FENCE_WAIT_PROBE_KM:
            return float(Config.FENCE_WAIT_PROBE_KM)
        return (int(Config.MDL_DILATE_ITER) + 2) * PW / 1000.0

    def v910_footprint(ct, mi, thr):
        """Sigungu statistics at leg-month index mi from the county series."""
        if not ct:
            return None
        rows_c, series_c = ct
        vals, present, n_tot = [], 0, 0
        for r in rows_c:
            n_tot += 1
            s = series_c.get(r.get("zid"))
            if s is None or r.get("status") == "no modelled cells" or mi < 0 or mi >= len(s["mean"]):
                continue
            mv, mx = float(s["mean"][mi]), float(s["mx"][mi])
            if mv >= thr:
                vals.append(100.0 * mv)
            elif mx >= thr:
                present += 1
        v = np.asarray(vals, dtype=np.float64)
        return dict(n_total=n_tot, n_reached=int(v.size), n_present=int(present),
                    n_ge80=int((v >= 80.0).sum()),
                    median=float(np.median(v)) if v.size else float("nan"),
                    q25=float(np.percentile(v, 25)) if v.size else float("nan"),
                    q75=float(np.percentile(v, 75)) if v.size else float("nan"))

    def v910_korea_reports():
        """Items 3–4: fence wait, Table 1, first sigungu per province, Table 3,
        Figure 9 and the run summary, all on ROK territory."""
        if not KOREA_LEG.get("cum"):
            logger.warning("\n[v9.10] the peninsular leg did not run — Tables 1/3, the fence-wait "
                           "table and Figure 9 were skipped.")
            return
        out = Config.OUT_DIR
        files = []
        cum_k = KOREA_LEG["cum"]
        sl = KOREA_LEG["sl"]
        r0_, c0_, a_b = int(KOREA_LEG["r0"]), int(KOREA_LEG["c0"]), int(KOREA_LEG["a_border"])
        KP_k, NK_k, SK_k = KOREA_LEG["KP"], KOREA_LEG["NK"], KOREA_LEG["SK"]
        W_k = KOREA_LEG.get("ROK_W")
        if W_k is None:
            W_k = SK_k.astype(np.float64)
        EX_k = np.asarray(KOREA_LEG.get("EXCL", np.zeros(KP_k.shape, bool)), bool) & KP_k
        cell_km = PW / 1000.0
        thr = float(Config.RAD_ARRIVAL_THRESHOLD)
        em, ea = float(Config.SK_ESTABLISHED_MEAN_ATTACK), float(Config.SK_ESTABLISHED_AREA_FRACTION)
        reps = sorted({int(g) for g in (Config.KOREA_REPORT_MONTHS or [])})
        reps_all = sorted(set(reps) | {int(Config.SA_REPORT_MONTH)})
        probe = v910_probe_km()
        logger.info("\n" + "=" * 96)
        logger.info(f"[v9.10 PENINSULAR-LEG REPORTS] ROK = {ROK_SPLIT_INFO.get('method')}; "
                    f"A* = {thr}; estuary '{Config.KOREA_ESTUARY_POLICY}'; exclaves "
                    f"'{ROK_EXCLAVE_INFO.get('policy')}' ({ROK_EXCLAVE_INFO.get('cells', 0)} cells)")
        logger.info("=" * 96)
        if ROK_RASTER is None:
            logger.warning("  ⚠ the national split failed in this run: every 'ROK' value below "
                           "includes DPRK land west of the fence.")
        met, S0 = korea_leg_metrics(cum_k, NK_k, SK_k, W_k, thr, a_b, em, ea, reps_all,
                                    first_exclude=EX_k)
        V910_STATE["baseline_metrics"] = met
        tta = first_arrival_index(cum_k, thr, mask=KP_k)
        _zb = {int(z["zid"]): z for z in (SK_COUNTY_ZONES or [])}

        def county_of(r, c):
            R, C = int(r) + r0_, int(c) + c0_
            if not (0 <= R < NR and 0 <= C < NC):
                return "", "", ""
            z = _zb.get(int(SK_ZONE_RASTER[R, C]))
            if not z:
                return "", "", ""
            return (str(z.get("name_en") or ""), str(z.get("name_kr") or ""),
                    str(z.get("province") or ""))

        def county_list(cells):
            names = []
            for (r, c) in cells:
                en, kr, pv = county_of(r, c)
                if en or kr:
                    s_ = f"{en or kr}" + (f" ({pv})" if pv else "")
                    if s_ not in names:
                        names.append(s_)
            return "; ".join(names)

        corr_k = [dict(name=co["name"], weight=float(co["weight"]),
                       row=int(co["row"]) - r0_, col=int(co["col"]) - c0_) for co in CORRIDORS]
        rok_main = SK_k & ~EX_k
        # ---- item 3: fence wait ------------------------------------------------------------
        fw = fence_wait_table(tta, corr_k, rok_main, NK_k,
                              KOREA_LEG.get("BAND", MDL_BOUNDARY[sl]), cell_km,
                              probe, a_b,
                              contact_rok=KOREA_LEG.get("CONTACT_ROK"))
        for row in fw:
            row["rok_first_sigungu"] = county_list(row.pop("rok_first_cells", []) or [])
        files.append(v910_write_csv(
            out / f"{V910_TAG}_fence_wait_by_gate.csv",
            ["zone", "kind", "weight", "dprk_first_month", "dprk_median_month", "band_first_month",
             "rok_band_first_month", "rok_first_month", "rok_median_month", "wait_first_months",
             "wait_median_months",
             "rok_first_calendar", "rok_first_sigungu", "n_dprk_cells", "n_band_cells", "n_rok_cells"],
            fw,
            [f"{V910_TAG} fence wait (H1) — zone = cells within {probe:.0f} km "
             f"of each gate, or of the unfenced contact; months on the global index",
             "dprk_* / rok_* = DPRK / ROK cells OUTSIDE the fence band; band_first = any band cell; "
             "rok_band_first = band cells on ROK territory; wait = ROK first − DPRK first (H1)",
             "median = month by which half of the zone's cells had reached A*"], logger))
        logger.info("  [FENCE WAIT] DPRK side first → ROK side first (months):")
        for row in fw:
            logger.info(f"    {row['zone'][:40]:<40s} DPRK {str(row['dprk_first_month']):>4s} → ROK "
                        f"{str(row['rok_first_month']):>4s}  wait {str(row['wait_first_months']):>4s}  "
                        f"{row['rok_first_sigungu'][:60]}")
        # ---- item 4: Table 1 (timeline + state of the ROK) ---------------------------------
        loc = first_crossing_location(tta, rok_main, corr_k, KOREA_LEG.get("CONTACT_ROK"),
                                      cell_km, probe)

        def cells_at(mask, g):
            if g is None:
                return []
            return [tuple(int(x) for x in rc) for rc in np.argwhere(mask & (tta == g - 1 - a_b))]

        t1 = []

        def ev(event, month, value="", detail=""):
            t1.append(dict(event=event, month=month, calendar=v910_calendar(month),
                           value=value, detail=detail))
        ev("Introduction into the DPRK (continental stopping criterion met)", a_b + 1)
        ev("Mean attack rate in Ω_NK reaches 5%", met["nk_5"])
        whole = next((r for r in fw if r["kind"] == "fence"), None)
        if whole:
            ev("Front reaches the DPRK side of the fence (first cell next to the band)",
               whole["dprk_first_month"])
        ev("First infection in the ROK (first ROK cell, fence band included)", met["first_rok_cell"],
           detail=(county_list(cells_at(rok_main, met["first_rok_cell"]))
                   + (f"; first crossing: {loc['label']} ({loc['dist_km']:.1f} km)"
                      if loc.get("dist_km") is not None else "")))
        _band_now = np.asarray(KOREA_LEG.get("BAND", MDL_BOUNDARY[sl]), bool)
        _tb = tta[rok_main & ~_band_now]
        _tb = _tb[np.isfinite(_tb)]
        _g_beyond = (int(_tb.min()) + 1 + a_b) if _tb.size else None
        ev("First ROK cell beyond the fence band", _g_beyond,
           detail=county_list(cells_at(rok_main & ~_band_now, _g_beyond)))
        ev("First infection in any cell holding ROK land (incl. border cells)", met["first_rok_any"])
        for r in sorted([r for r in fw if r["kind"] in ("gate", "estuary")],
                        key=lambda r: (r["rok_first_month"] is None, r["rok_first_month"] or 0)):
            lab = (f"First infection south of the {r['zone']} gate (beyond the fence band)"
                   if r["kind"] == "gate" else "First infection across the unfenced Han River estuary")
            ev(lab, r["rok_first_month"],
               detail=f"{r['rok_first_sigungu']}; DPRK side first reached in month "
                      f"{r['dprk_first_month']}; wait {r['wait_first_months']} months; ROK band cells "
                      f"first reached in month {r['rok_band_first_month']}")
        ev("1% of the area of the ROK reached", met["area_1"])
        ev("Mean attack rate over the ROK reaches 1%", met["mean_1"])
        ev("5% of the area of the ROK reached", met["area_5"])
        ev("10% of the area of the ROK reached", met["area_10"])
        ev("Mean attack rate over the ROK reaches 5%", met["mean_5"])
        _ae = met.get("at_est") or {}
        ev("Establishment criterion for the ROK met (area-weighted)", met["est_areaw"],
           detail=(f"mean A {_ae.get('rok_mean_w_pct', float('nan')):.1f}%; "
                   f"{_ae.get('rok_area_w_pct', float('nan')):.1f}% of the area reached; "
                   f"criterion on ROK cells: month {met['est_cell']}"))
        ev("End of the peninsular leg", met["end_month"],
           detail=(f"the leg ends at the later of establishment (month {KOREA_LEG.get('stop_month')}) + "
                   f"{Config.KOREA_POST_ARRIVAL_MONTHS if Config.KOREA_POST_ARRIVAL_MONTHS is not None else 12 * int(Config.KOREA_POST_ARRIVAL_YEARS)}"
                   f" months and the last report month ({max(reps) if reps else '—'})"))
        ct = sk_county_tables.get("korea_leg")
        for g in [met["end_month"]] + [g for g in reps if g != met["end_month"]]:
            at_ = met["at_end"] if g == met["end_month"] else met["at"].get(g)
            if not at_:
                ev(f"(report month {g} lies beyond the end of the leg)", g)
                continue
            ev(f"State at month {g}: area of the ROK reached (A ≥ A*)", g, f"{at_['rok_area_w_pct']:.1f}%")
            ev(f"State at month {g}: mean cumulative attack rate over the ROK", g,
               f"{at_['rok_mean_w_pct']:.1f}%")
            ev(f"State at month {g}: mean cumulative attack rate in Ω_NK", g, f"{at_['nk_mean_pct']:.1f}%")
            fp_ = v910_footprint(ct, g - 1 - a_b, thr)
            if fp_:
                ev(f"State at month {g}: sigungu reached (mean A ≥ A*), of {fp_['n_total']}", g,
                   f"{fp_['n_reached']}", f"+{fp_['n_present']} with A ≥ A* in some cells only; "
                                         f"{fp_['n_ge80']} sigungu ≥ 80%")
                ev(f"State at month {g}: mean A in the reached sigungu, median (IQR)", g,
                   f"{fp_['median']:.1f}% ({fp_['q25']:.1f}–{fp_['q75']:.1f}%)")
        if ROK_EXCLAVE_INFO.get("cells"):
            ev("Note: ROK exclave cells on this grid (islands joined to the DPRK coast)", None,
               ROK_EXCLAVE_INFO.get("policy"), ", ".join(ROK_EXCLAVE_INFO.get("sigungu", [])))
        files.append(v910_write_csv(
            out / f"{V910_TAG}_Table1_timeline.csv", ["event", "month", "calendar", "value", "detail"], t1,
            [f"{V910_TAG} Table 1 — peninsular leg on ROK territory ({ROK_SPLIT_INFO.get('method')})",
             "months on one continuous index shared by both legs (month 1 = January of Year 1)",
             f"'reached' = cumulative attack rate A ≥ A* = {thr}; ROK values area-weighted by the ROK "
             f"share of each cell; sigungu values from the county table"], logger))
        for r in t1:
            logger.info(f"    {str(r['month']) if r['month'] is not None else '':>4s}  {r['event'][:78]:<78s} "
                        f"{r['value']}")
        # ---- first sigungu per province --------------------------------------------------------
        if ct:
            rows_c = ct[0]
            byp, cnt = {}, {}
            for r in rows_c:
                am = r.get("arrival_month")
                if am is None:
                    continue
                pv = r.get("province") or province_from_sig_cd(r.get("sig_cd"))[0] or "?"
                cnt[pv] = cnt.get(pv, 0) + 1
                if pv not in byp or am < byp[pv]["arrival_month"]:
                    byp[pv] = r
            prov = [dict(province=pv, first_sigungu=r.get("county_en") or r.get("county_kr"),
                         sig_cd=r.get("sig_cd"), arrival_month=r["arrival_month"],
                         calendar=v910_calendar(r["arrival_month"]), sigungu_reached=cnt.get(pv, 0))
                    for pv, r in sorted(byp.items(), key=lambda kv: kv[1]["arrival_month"])]
            files.append(v910_write_csv(
                out / f"{V910_TAG}_first_sigungu_by_province.csv",
                ["province", "first_sigungu", "sig_cd", "arrival_month", "calendar", "sigungu_reached"],
                prov, [f"{V910_TAG} first sigungu reached (mean A ≥ A*) in each province, "
                       f"by month {met['end_month']}"], logger))
        # ---- item 4: Table 3 (vaccination on ROK territory) ----------------------------------
        vac = KOREA_LEG.get("vacc") or {}
        scen = [("baseline", cum_k, None)]
        for nm, d in vac.items():
            if d.get("von") is not None:
                st = int(d["von"]) + 1 + a_b
            elif d.get("triggered"):
                st = None
            else:
                st = int(d["start"]) + 1 + a_b
            scen.append((nm, d["cum"], st))
        t3, series9, base_at = [], {}, {}
        for key, cum_s, st in scen:
            m_s, S_s = korea_leg_metrics(cum_s, NK_k, SK_k, W_k, thr, a_b, em, ea, reps_all,
                                         first_exclude=EX_k)
            series9[key] = (S_s, st)
            for g in [m_s["end_month"]] + [g for g in reps if g != m_s["end_month"]]:
                at_ = m_s["at_end"] if g == m_s["end_month"] else m_s["at"].get(g)
                if not at_:
                    continue
                if key == "baseline":
                    base_at[g] = at_
                b_ = base_at.get(g)
                t3.append(dict(
                    scenario=V910_VACC_LABELS.get(key, key), reference_month=g,
                    start_month=(st if st is not None else ("—" if key == "baseline" else "never")),
                    mean_A_nk_pct=at_["nk_mean_pct"], mean_A_rok_pct=at_["rok_mean_w_pct"],
                    reduction_nk_pct=(100.0 * (1.0 - at_["nk_mean_pct"] / b_["nk_mean_pct"])
                                      if (b_ and key != "baseline" and b_["nk_mean_pct"] > 0) else None),
                    reduction_rok_pct=(100.0 * (1.0 - at_["rok_mean_w_pct"] / b_["rok_mean_w_pct"])
                                       if (b_ and key != "baseline" and b_["rok_mean_w_pct"] > 0) else None),
                    rok_area_reached_pct=at_["rok_area_w_pct"],
                    first_rok=m_s["first_rok_cell"], first_rok_any=m_s["first_rok_any"],
                    est_rok=(m_s["est_areaw"] if m_s["est_areaw"] is not None else "not met"),
                    nk_5pct_month=m_s["nk_5"]))
        files.append(v910_write_csv(
            out / f"{V910_TAG}_Table3_vaccination.csv",
            ["scenario", "reference_month", "start_month", "mean_A_nk_pct", "mean_A_rok_pct",
             "reduction_nk_pct", "reduction_rok_pct", "rok_area_reached_pct", "first_rok",
             "first_rok_any", "est_rok", "nk_5pct_month"], t3,
            [f"{V910_TAG} Table 3 — vaccination scenarios on the peninsular leg, ROK territory",
             "reference_month = end of the leg, then each KOREA_REPORT_MONTHS entry; months global",
             "A = cumulative attack rate; ROK means/areas area-weighted; reductions vs no vaccination",
             f"V1 starts at leg month {Config.KOREA_VACC_V1_START_LEG_MONTH}; V2 starts when the first "
             f"ROK cell reaches A* (start_month = month the campaign began)"], logger))
        for r in t3:
            logger.info(f"    {r['scenario'][:46]:<46s} @{r['reference_month']}: NK {r['mean_A_nk_pct']:5.1f}%  "
                        f"ROK {r['mean_A_rok_pct']:5.1f}% (−{(r['reduction_rok_pct'] or 0):.1f}%)  "
                        f"ROK reached {r['rok_area_reached_pct']:5.1f}%  first {r['first_rok']}  est {r['est_rok']}")
        # ---- item 4: Figure 9 (maps A–C + time series D, no figure title) -----------------
        dk = KOREA_LEG.get("draw_korea")
        _rows9 = []
        for key, (S_s, st) in series9.items():
            for i in range(len(S_s["rok_mean_w"])):
                _rows9.append(dict(scenario=key, month=i + 1 + a_b,
                                   rok_mean_w_pct=100.0 * S_s["rok_mean_w"][i],
                                   rok_area_w_pct=100.0 * S_s["rok_area_w"][i],
                                   nk_mean_pct=100.0 * S_s["nk_mean"][i]))
        files.append(v910_write_csv(out / f"{V910_TAG}_Figure9_series.csv",
                                    ["scenario", "month", "rok_mean_w_pct", "rok_area_w_pct", "nk_mean_pct"],
                                    _rows9, [f"{V910_TAG} Figure 9 panel D data (ROK area-weighted)"], logger))
        if Config.RENDER_FIGURE_9 and dk is not None and vac:
            try:
                g9 = int(Config.FIGURE9_MONTH) if Config.FIGURE9_MONTH else int(met["end_month"])
                panels = [("baseline", cum_k)] + [(nm, d["cum"]) for nm, d in vac.items()]
                fig9 = plt.figure(figsize=(21, 14), facecolor="white")
                gs9 = fig9.add_gridspec(2, len(panels[:3]), height_ratios=[1.45, 1.0], hspace=0.10, wspace=0.04)
                for i, (key, cum_s) in enumerate(panels[:3]):
                    ax9 = fig9.add_subplot(gs9[0, i])
                    mi9 = int(np.clip(g9 - 1 - a_b, 0, len(cum_s) - 1))
                    pct9 = np.where(KP_k, np.nan_to_num(cum_s[mi9]) * 100.0, np.nan)
                    dk(ax9, pct9, f"{V910_VACC_LABELS.get(key, key).split(':')[0]} — month {mi9 + 1 + a_b}")
                    ax9.text(0.01, 0.99, "ABC"[i], transform=ax9.transAxes, fontsize=18,
                             fontweight="bold", va="top")
                _bb = ax9.get_position()
                _bn9 = BoundaryNorm(list(Config.SPREAD_PCT_BANDS), ncolors=len(Config.SPREAD_PCT_BANDS) - 1,
                                    clip=True)
                for _j, (_cols, _lab) in enumerate(((Config.SPREAD_BAND_COLORS_NK, "DPRK %"),
                                                    (Config.SPREAD_BAND_COLORS_SK, "ROK %"))):
                    _cax = fig9.add_axes([_bb.x1 + 0.012 + 0.050 * _j, _bb.y0 + _bb.height * 0.30,
                                          0.010, _bb.height * 0.60])
                    _sm9 = plt.cm.ScalarMappable(cmap=ListedColormap(_cols), norm=_bn9)
                    _sm9.set_array([])
                    _cb9 = fig9.colorbar(_sm9, cax=_cax, boundaries=list(Config.SPREAD_PCT_BANDS),
                                         ticks=list(Config.SPREAD_PCT_BANDS), spacing="proportional")
                    _cb9.ax.set_yticklabels([f"{int(b)}" for b in Config.SPREAD_PCT_BANDS], fontsize=9)
                    _cb9.set_label(_lab, fontsize=10, fontweight="bold")
                axD = fig9.add_subplot(gs9[1, :])
                _col9 = {"baseline": "#d62728", "V1_combined_host_peninsula": "#2ca02c",
                         "V2_ROK_livestock_on_detection": "#1f77b4"}
                for key, (S_s, st) in series9.items():
                    xs = np.arange(len(S_s["rok_mean_w"])) + a_b + 1
                    c9 = _col9.get(key, "#555555")
                    lab9 = V910_VACC_LABELS.get(key, key)
                    axD.plot(xs, 100.0 * S_s["rok_mean_w"], color=c9, lw=2.2, label=f"{lab9} — ROK")
                    axD.plot(xs, 100.0 * S_s["nk_mean"], color=c9, lw=1.4, ls="--", label=f"{lab9} — Ω_NK")
                    if st is not None and key != "baseline":
                        axD.axvline(st, color=c9, lw=1.0, ls=":", alpha=0.9)
                axD.set_xlim(a_b + 1, a_b + len(cum_k))
                axD.set_ylim(bottom=0)
                axD.set_xlabel("Month", fontsize=12)
                axD.set_ylabel("Mean cumulative attack rate (%)", fontsize=12)
                axD.grid(alpha=0.3)
                axD.legend(fontsize=9 * (1.45 if getattr(Config, "PUB_FIGS", False) else 1.0), ncol=3,
                           loc="upper center", bbox_to_anchor=(0.5, -0.13), frameon=False)
                axD.text(-0.035, 1.02, "D", transform=axD.transAxes, fontsize=18, fontweight="bold")
                fp9 = out / f"{V910_TAG}_Figure9_vaccination.png"
                fig9.savefig(fp9, dpi=200, facecolor="white", bbox_inches="tight")
                plt.close(fig9)
                files.append(fp9)
                logger.info(f"   ✓ Figure 9 (no figure title; panels A–C at month {g9}): {fp9}")
            except Exception as _e:
                logger.warning(f"   ⚠ Figure 9 skipped ({type(_e).__name__}: {_e})", exc_info=True)
        # ---- run summary (read by the scenario table of another run, e.g. the 4-km grid) ------
        whole_first = whole["dprk_first_month"] if whole else None
        at_rep = met["at"].get(int(Config.SA_REPORT_MONTH)) or {}
        # same wait definition as Table 4 / the sensitivity runs: first ROK cell beyond
        # the fence band − first DPRK cell next to the band
        summ = dict(
            version="9.10", created=_v910_time.strftime("%Y-%m-%d %H:%M"),
            out_dir=str(out), cell_km=cell_km, grid=[int(NR), int(NC)],
            rok_split=ROK_SPLIT_INFO, exclaves=ROK_EXCLAVE_INFO,
            estuary_policy=str(Config.KOREA_ESTUARY_POLICY), livestock_file=Config.LIVESTOCK_USED,
            transport_reference_cell_km=Config.TRANSPORT_REFERENCE_CELL_KM,
            baseline=dict(intro_month=a_b + 1, first_rok=met["first_rok_cell"],
                          first_rok_any=met["first_rok_any"], first_location=loc["label"],
                          first_location_km=loc.get("dist_km"), est=met["est_areaw"],
                          est_cell=met["est_cell"], end_month=met["end_month"],
                          transit_months=(None if met["first_rok_cell"] is None
                                          else met["first_rok_cell"] - (a_b + 1)),
                          transit_years=(None if met["first_rok_cell"] is None
                                         else (met["first_rok_cell"] - (a_b + 1)) / 12.0),
                          fence_dprk_first=whole_first, first_rok_beyond_band=_g_beyond,
                          wait_at_fence=(None if (whole_first is None or _g_beyond is None)
                                         else _g_beyond - whole_first),
                          pct_rok_reached_report=at_rep.get("rok_area_w_pct"),
                          mean_A_rok_report=at_rep.get("rok_mean_w_pct"),
                          mean_A_nk_report=at_rep.get("nk_mean_pct"), threshold=thr),
            seeding={k: v for k, v in V910_SEED_INFO.items() if not str(k).startswith("_")},
            metrics=met, report_month=int(Config.SA_REPORT_MONTH))
        _fps = out / f"{V910_TAG}_run_summary.json"
        _fps.write_text(_v910_json.dumps(summ, ensure_ascii=False, indent=2,
                                         default=lambda o: o.tolist() if isinstance(o, (np.ndarray, np.generic)) else str(o)),
                        encoding="utf-8")
        files.append(_fps)
        logger.info(f"   ✓ run summary: {_fps}")
        v910_copy_to_dl(files, Config.DL_DIR, logger)

    def v910_export_rasters():
        """Item 7: GeoTIFFs on the model grid — SMCDA scores and classes, arrival
        months of both legs, the ROK/DPRK split, MaxEnt — for an external test."""
        if not Config.EXPORT_RISK_RASTERS:
            return
        outd = Config.OUT_DIR / f"{V910_TAG}_rasters"
        n = 0
        for _dk, _sm in SMCDA_BY_DOMAIN.items():
            write_geotiff(outd / f"{V910_TAG}_smcda_{_dk}_score.tif",
                          np.where(_sm["mask"], _sm["score"], np.nan), GT5, WKT)
            write_geotiff(outd / f"{V910_TAG}_smcda_{_dk}_class.tif",
                          np.where(_sm["mask"], _sm["classes"], np.nan), GT5, WKT,
                          nodata=-1, dtype="int16")
            n += 2
        # arrival months: cells of the domain never reached get (last month + 1), so an
        # external test counts outbreaks there as 'latest' instead of dropping them
        if KOREA_LEG.get("cum"):
            _kpk = KOREA_LEG["KP"]
            _t = first_arrival_index(KOREA_LEG["cum"], float(Config.RAD_ARRIVAL_THRESHOLD), _kpk)
            _end = len(KOREA_LEG["cum"]) + int(KOREA_LEG["a_border"])
            _full = np.full(SH, np.nan)
            _full[KOREA_LEG["sl"]] = np.where(_kpk, np.where(np.isfinite(_t), _t + 1 + int(KOREA_LEG["a_border"]),
                                                             _end + 1), np.nan)
            write_geotiff(outd / f"{V910_TAG}_korea_leg_arrival_month.tif", _full, GT5, WKT)
            n += 1
        _res = rad_results_holder.get("res")
        if _res is not None and _res.get("time_to_arrival") is not None:
            _tc = np.asarray(_res["time_to_arrival"], np.float64)
            _nmc = len(_res.get("monthly_cum") or [])
            write_geotiff(outd / f"{V910_TAG}_continental_arrival_month.tif",
                          np.where(KP, np.where(np.isfinite(_tc), _tc + 1.0, _nmc + 1.0), np.nan), GT5, WKT)
            n += 1
        if ROK_RASTER is not None:
            _cls = np.zeros(SH, np.float64)
            _cls[KOREA_REGION & ROK_RASTER] = 1
            _cls[KOREA_REGION & ~ROK_RASTER] = 2
            _cls[ROK_EXCLAVES] = 3
            write_geotiff(outd / f"{V910_TAG}_rok_dprk_split.tif", _cls, GT5, WKT, nodata=255, dtype="byte")
            write_geotiff(outd / f"{V910_TAG}_rok_area_fraction.tif",
                          np.where(KOREA_REGION, ROK_FRACTION, np.nan), GT5, WKT)
            n += 2
        if MAXENT is not None:
            write_geotiff(outd / f"{V910_TAG}_ppr_suitability_on_grid.tif", np.where(KP, MAXENT, np.nan), GT5, WKT)
            n += 1
        logger.info(f"   ✓ {n} GeoTIFFs for the external test: {outd} "
                    f"(split raster: 1 = ROK, 2 = DPRK, 3 = ROK exclave; arrival rasters: "
                    f"last month + 1 = never reached in the run)")

    def v910_build_context():
        """Items 5–6: everything a peninsular-leg variant needs, with self-checks,
        saved to PPR_v9.10_korea_leg_context.npz."""
        if not KOREA_LEG.get("cum"):
            return
        import dataclasses as _dcl
        sl = KOREA_LEG["sl"]
        r0_, c0_, a_b = int(KOREA_LEG["r0"]), int(KOREA_LEG["c0"]), int(KOREA_LEG["a_border"])
        env_k, p_k = KOREA_LEG["env"], KOREA_LEG["p"]
        KP_k = KOREA_LEG["KP"]
        H, W = KP_k.shape
        cell_km = PW / 1000.0
        params = {f_.name: getattr(p_k, f_.name) for f_ in _dcl.fields(p_k)}
        for _k in ("vacc_trigger_mask", "vacc_cap_field"):
            params.pop(_k, None)
        corr_k = [dict(name=str(co["name"]), weight=float(co["weight"]),
                       row=int(co["row"]) - r0_, col=int(co["col"]) - c0_) for co in CORRIDORS]
        inje = None
        try:
            _lat, _lon = Config.INJE_SEORAK_ESTIMATE_LATLON
            _ij = build_corridors(MDL_BOUNDARY, GT5, WKT, NR, NC, cell_km=cell_km, kp=VALID,
                                  crossings=[("Inje_Seorak", float(_lat), float(_lon),
                                              float(Config.INJE_SEORAK_WEIGHT),
                                              "ESTIMATED position (v9.10 scenario only)")],
                                  snap_max_km=Config.CORRIDOR_SNAP_MAX_KM,
                                  gate_radius_km=Config.CORRIDOR_GATE_RADIUS_KM, logger=logger)
            if _ij:
                inje = dict(name="Inje_Seorak", weight=float(Config.INJE_SEORAK_WEIGHT),
                            row=int(_ij[0]["row"]) - r0_, col=int(_ij[0]["col"]) - c0_)
        except Exception as _e:
            logger.warning(f"   ⚠ Inje–Seorak gate not placed ({_e})")
        yalu_seeds, a_yalu, ynote = [], a_b, "lower-Yalu entry: introduction month kept at the baseline's"
        try:
            _lon, _lat = Config.YALU_ENTRY_LONLAT
            _yr, _yc, _inside, _x, _y = lonlat_to_rc(_lon, _lat, GT5, NR, NC, make_lonlat_transform(WKT))
            _nk = KOREA_LEG["NK"] & KP_k
            _rr, _cc = np.where(_nk)
            if _rr.size:
                _d = np.hypot(_rr + r0_ - _yr, _cc + c0_ - _yc)
                _ord = np.argsort(_d, kind="mergesort")[:max(1, len(KOREA_LEG["seeds"]))]
                yalu_seeds = [(int(_rr[i]), int(_cc[i]), "yalu_entry", 0, 1.0) for i in _ord]
            _res = rad_results_holder.get("res")
            if (_res is not None and CHINA_MASK.any() and not Config.KOREA_ONLY_INTRO_MONTH):
                _RRf, _CCf = np.mgrid[0:NR, 0:NC]
                _probe = CHINA_MASK & (np.hypot(_RRf - _yr, _CCf - _yc) * cell_km
                                       <= float(Config.YALU_CHINA_PROBE_KM))
                if _probe.any():
                    _curve = np.array([float(np.nanmean(np.nan_to_num(c)[_probe]))
                                       for c in _res["monthly_cum"]])
                    _hit = np.flatnonzero(_curve >= 0.05)
                    if _hit.size:
                        a_yalu = int(_hit[0])
                        ynote = (f"China cells within {Config.YALU_CHINA_PROBE_KM:.0f} km of the lower Yalu "
                                 f"reached a mean A of 5% in month {a_yalu + 1}")
                    else:
                        ynote = ("China cells near the lower Yalu never reached 5% in the continental run; "
                                 "introduction month kept at the baseline's")
            logger.info(f"   Yalu scenario: {len(yalu_seeds)} DPRK seed cells nearest "
                        f"({_lat:.2f}N, {_lon:.2f}E); {ynote}.")
        except Exception as _e:
            logger.warning(f"   ⚠ Yalu entry not prepared ({type(_e).__name__}: {_e})")
        score, sdom = None, None
        for _dk in ("sk", "sk_reached", "korea"):
            if _dk in SMCDA_BY_DOMAIN:
                _s = SMCDA_BY_DOMAIN[_dk]
                score = np.where(_s["mask"], _s["score"], np.nan)[sl].astype(np.float32)
                sdom = _dk
                break
        fence0 = KOREA_LEG.get("FENCE_ORIG")
        mdl0 = KOREA_LEG.get("MDL_ORIG", MDL_BOUNDARY[sl].astype(bool))
        fence_mean = (float(np.mean(fence0[mdl0])) if (fence0 is not None and mdl0.any())
                      else float(Config.HAN_FENCE_REDUCTION))
        _korea_nomask = KP & np.isfinite(arr5)
        bm = V910_STATE.get("baseline_metrics") or {}
        ctx = dict(
            version="9.10", shape=(int(H), int(W)), cell_km=float(cell_km), r0=r0_, c0=c0_,
            a_border=a_b, a_border_yalu=int(a_yalu), yalu_note=ynote,
            KP=KP_k, KP_NOMASK=_korea_nomask[sl], NK=(NK & _korea_nomask)[sl],
            SK=(SK & _korea_nomask)[sl],
            ROK_FRAC=(np.clip(np.nan_to_num(ROK_FRACTION[sl]), 0.0, 1.0) if ROK_FRACTION is not None
                      else (SK & _korea_nomask)[sl].astype(np.float64)),
            EXCL=ROK_EXCLAVES[sl].copy(), CONTACT_ROK=KOREA_LEG.get("CONTACT_ROK"),
            CONTACT_DPRK=KOREA_LEG.get("CONTACT_DPRK"), CONTACT_SEAL=KOREA_LEG.get("CONTACT_SEAL"),
            MDL=mdl0, FENCE=fence0, fence_mean=fence_mean,
            GOR_EFF=env_k["GOR_EFF"], GOR_nn=env_k["GOR_nn"], GORANI_DENSITY=env_k["GORANI_DENSITY"],
            MAXENT=(None if MAXENT is None else np.asarray(MAXENT, np.float32)[sl]),
            CONDUCTANCE=env_k["CONDUCTANCE"], RESISTANCE=env_k["RESISTANCE"], RIV_nn=env_k["RIV_nn"],
            WATER_ATTRACT=env_k["WATER_ATTRACT"], livestock_nn=env_k["livestock_nn"],
            SCORE=score, score_domain=sdom,
            seeds=[tuple(s) for s in KOREA_LEG["seeds"]], yalu_seeds=yalu_seeds,
            corridors=corr_k, inje_gate=inje, params=params,
            fence_kw=dict(gate_radius_km=float(Config.CORRIDOR_GATE_RADIUS_KM),
                          max_pass=float(Config.CORRIDOR_MAX_PASS),
                          sealed=float(Config.MDL_SEALED_FACTOR),
                          slowdown=float(Config.MDL_CROSSING_SLOWDOWN)),
            han_phi=float(Config.HAN_FENCE_REDUCTION),
            habitat=dict(use_composite=bool(Config.USE_COMPOSITE_HABITAT_FOR_SPREAD),
                         blend=float(Config.HABITAT_BLEND),
                         use_maxent=bool(Config.USE_MAXENT_IN_HABITAT),
                         base_is_maxent=bool(base_is_maxent),
                         maxent_blend=float(Config.MAXENT_HABITAT_BLEND)),
            thr=float(Config.RAD_ARRIVAL_THRESHOLD), est_mean=float(Config.SK_ESTABLISHED_MEAN_ATTACK),
            est_area=float(Config.SK_ESTABLISHED_AREA_FRACTION),
            report_month=int(Config.SA_REPORT_MONTH), probe_km=v910_probe_km(),
            cap_months=int(12 * int(Config.MAX_SIM_YEARS)),
            estuary_policy=str(Config.KOREA_ESTUARY_POLICY).lower(), livestock_file=Config.LIVESTOCK_USED,
            exclave_policy=str(ROK_EXCLAVE_INFO.get("policy", "report")),
            sa_baseline=dict(gamma=float(Config.GAMMA_RECOVERY), beta_level=1.0, beta_amp=1.0,
                             tv_amp=1.0, diff_max_x=1.0, adv_max_x=1.0, adv_gain_x=1.0,
                             kappa=float(Config.CONDUCTANCE_EXPONENT),
                             phi0=float(Config.MDL_SEALED_FACTOR), phimax_x=1.0,
                             rg_km=float(Config.CORRIDOR_GATE_RADIUS_KM),
                             f0=float(Config.CORRIDOR_SEASON_FLOOR), psi=float(Config.HABITAT_BLEND),
                             w_h=float(Config.RAD_TRANS_HAB_W), a_star=float(Config.RAD_ARRIVAL_THRESHOLD),
                             i0=float(Config.RAD_SEED_INTENSITY)),
            baseline_metrics=dict(first_rok=bm.get("first_rok_cell"), est=bm.get("est_areaw"),
                                  end_month=bm.get("end_month")),
            rok_split=str(ROK_SPLIT_INFO.get("method")))
        # self-checks: the variant builder must reproduce the baseline inputs exactly
        _chk = []
        _g = v910_gor_eff(ctx, float(Config.HABITAT_BLEND), KP_k)
        _chk.append(("GOR_EFF(ψ)", float(np.nanmax(np.abs(_g - np.nan_to_num(env_k["GOR_EFF"]))[KP_k]))
                     if KP_k.any() else 0.0))
        _cd = np.where(KP_k, np.exp(-float(Config.CONDUCTANCE_EXPONENT) * env_k["RESISTANCE"]), 1e-9)
        _chk.append(("CONDUCTANCE(κ)", float(np.nanmax(np.abs(_cd - env_k["CONDUCTANCE"])[KP_k]))
                     if KP_k.any() else 0.0))
        if fence0 is not None and corr_k:
            _ff = build_mdl_fence_field(corr_k, mdl0, (H, W), cell_km,
                                        gate_radius_km=Config.CORRIDOR_GATE_RADIUS_KM,
                                        max_pass=Config.CORRIDOR_MAX_PASS, sealed=Config.MDL_SEALED_FACTOR,
                                        slowdown=Config.MDL_CROSSING_SLOWDOWN, logger=None)
            _chk.append(("fence field", float(np.max(np.abs(_ff - fence0)[mdl0])) if mdl0.any() else 0.0))
        logger.info("   context self-checks (max abs difference, should be ~0): "
                    + ", ".join(f"{a} {b:.2e}" for a, b in _chk))
        if any(b > 1e-4 for _a, b in _chk):
            logger.warning("   ⚠ a variant would not reproduce the baseline inputs exactly — read the "
                           "self-check above before trusting the batch.")
        V910_STATE["ctx"] = ctx
        if Config.SAVE_KOREA_CONTEXT:
            _fpc = Config.OUT_DIR / f"{V910_TAG}_korea_leg_context.npz"
            save_korea_ctx(ctx, _fpc)
            logger.info(f"   ✓ peninsular-leg context: {_fpc} ({_fpc.stat().st_size / 1e6:.1f} MB) — "
                        f"re-run the batch without the continental leg: python {Path(globals().get('__file__', 'PPR_model_v9_10.py')).name} "
                        f"--from-context \"{_fpc}\" --batch scenarios")

    for _v910_fn in (v910_korea_reports, v910_export_rasters, v910_build_context):
        try:
            _v910_fn()
        except Exception as _e:
            logger.warning(f"  ⚠ {_v910_fn.__name__} failed ({type(_e).__name__}: {_e})", exc_info=True)
    if Config.RUN_SCENARIOS or str(Config.RUN_SENSITIVITY).lower() not in ("off", "none", ""):
        if V910_STATE.get("ctx") is None:
            logger.warning("  ⚠ v9.10 batch skipped: no peninsular-leg context (see the warnings above).")
        else:
            try:
                v910_run_batch(V910_STATE["ctx"], scenarios=bool(Config.RUN_SCENARIOS),
                               sa_mode=str(Config.RUN_SENSITIVITY), out_dir=Config.OUT_DIR,
                               dl_dir=Config.DL_DIR, workers=Config.SA_WORKERS, logger_=logger)
            except Exception as _e:
                logger.warning(f"  ⚠ v9.10 batch failed ({type(_e).__name__}: {_e})", exc_info=True)
    logger.info("\n" + "=" * 100)
    logger.info("✓ COMPLETE")
    if CORRIDORS:
        logger.info(f"   MDL gated at {len(CORRIDORS)} water-deer corridors "
                    f"(sealed {Config.MDL_SEALED_FACTOR}, gate peak "
                    f"{Config.CORRIDOR_MAX_PASS:.2f}×weight×{Config.MDL_CROSSING_SLOWDOWN:.2f}, "
                    f"radius {Config.CORRIDOR_GATE_RADIUS_KM:.0f} km)")
    logger.info(f"   South-Korea panels: "
                f"{'white canvas + counties + water' if Config.SK_PANEL_PLAIN else 'suitability underlay (v9.1 style)'}; "
                f"SK band floor {Config.SK_BAND_MIN_PCT_VISIBLE:.0f}% vs NK "
                f"{Config.BAND_MIN_PCT_VISIBLE:.0f}%")
    if Config.AUTO_RUN_UNTIL_SK_ESTABLISHED:
        logger.info(f"   Run length: ADAPTIVE — source run extended until the epidemic "
                    f"reached the North-Korean border (+{Config.CHINA_POST_BORDER_YEARS} yr), "
                    f"Korea leg extended until PPR established in South Korea "
                    f"(+{Config.KOREA_POST_ARRIVAL_YEARS} yr). Search the log for "
                    f"'MET at month' for the two dates.")
    # ── v9.3: the county tables ─────────────────────────────────────────────
    if Config.EXPORT_SK_COUNTY_CSV:
        if sk_county_outputs:
            logger.info("   SOUTH-KOREA COUNTY TABLES (the v9.3 deliverable):")
            for tag, files in sk_county_outputs.items():
                logger.info(f"     [{tag}]")
                for f in files:
                    logger.info(f"       {f}")
            logger.info("     -> open the *_summary.csv for one row per sigungu "
                        "(arrival month, final attack rate, % of county reached, "
                        "nearest DMZ corridor gate).")
        elif not SK_COUNTY_ZONES:
            _reason, _msg = SK_ZONE_FAIL or ("unknown", "no diagnosis recorded")
            logger.warning(f"   ⚠ NO county CSV written — cause: [{_reason}] {_msg}")
            _advice = {
                "not_found": ("Put sig.shp (with its .dbf/.shx/.prj siblings) under "
                              f"{Config.DATA_ROOT}, {Config.CWD} or {Config.ENV_DIR}, or add "
                              "the filename you do have to Config.SK_COUNTY_SHP_HINTS. The "
                              "list of shapefiles the model could see is logged above."),
                "outside_grid": ("The county layer and the base raster are in different "
                                 "coordinate systems, or a hint matched the wrong file. "
                                 "Compare the extents in the [DATA INVENTORY] block."),
                "all_north": ("Set Config.SK_COUNTY_DROP_NORTH_OF_MDL = False, or check the "
                              "DMZ row against MDL.shp."),
                "not_polygon": ("A hint matched a point/line layer. Name the polygon layer "
                                "explicitly at the front of Config.SK_COUNTY_SHP_HINTS."),
                "unreadable": ("A .shp without its .dbf/.shx/.prj cannot be read. Copy all "
                               "sidecar files, not just the .shp."),
            }.get(_reason)
            if _advice:
                logger.warning(f"     -> {_advice}")
            logger.warning("     -> The epidemic rasters were cached, so you do NOT need to "
                           "re-run the simulation: fix the layer, then run "
                           "`python sk_county_tools.py --rebuild` against the cache below.")
        else:
            logger.warning("   ⚠ NO county CSV written: the export ran but produced nothing — "
                           "see the warnings above.")
    for mode, fp in output_files.items():
        logger.info(f"   {mode}: {fp}")
    logger.info("=" * 100 + "\n")
# ══════════════════════════════════════════════════════════════════════════════
# v9.10 COMMAND LINE
#   python PPR_model_v9_10.py                              baseline run of record
#   python PPR_model_v9_10.py --batch scenarios            + Table 4 scenarios
#   python PPR_model_v9_10.py --batch all                  + scenarios, Morris, LHS/PRCC
#   python PPR_model_v9_10.py --from-context <npz> --batch morris+lhs
#                                                          batch only, no re-run
#   python PPR_model_v9_10.py --preset grid4km             4-km peninsular leg
#   python PPR_model_v9_10.py --preset gamma_x2            full re-run, γ and β × 2
#   python PPR_model_v9_10.py --preset holdout2015         full re-run seeded only from the
#                                                          outbreaks reported before 2015
#                                                          (external test, item 7)
#   python PPR_model_v9_10.py --seed-before 2019-01-01     the same with another cut-off
# The same switches can be given as environment variables PPR_PRESET / PPR_BATCH.
# ══════════════════════════════════════════════════════════════════════════════
V910_PRESETS = {
    "baseline": {},
    # 'baseline' = the DPRK introduction month of the 8-km run in the same results folder
    # (read from its PPR_v9.10_run_summary.json; 18, the 6 Nov value, if it is not there yet)
    "grid4km": dict(BASE_TARGET_MAX_PX=1500, MAX_GRID_CELLS=2_800_000,
                    TRANSPORT_REFERENCE_CELL_KM="auto", KOREA_ONLY_INTRO_MONTH="baseline",
                    RAD_MAX_SUBSTEPS=240),
    "korea_only": dict(KOREA_ONLY_INTRO_MONTH="baseline"),
    # round 2: grid-convergence check (8, 4 and 2 km); the MaxEnt base raster is 1 km
    "grid2km": dict(BASE_TARGET_MAX_PX=3000, MAX_GRID_CELLS=6_000_000,
                    TRANSPORT_REFERENCE_CELL_KM="auto", KOREA_ONLY_INTRO_MONTH="baseline",
                    RAD_MAX_SUBSTEPS=960, KOREA_VACC_SCENARIOS=False),
    "gamma_x2": dict(_gamma_mult=2.0),
    "gamma_x4": dict(_gamma_mult=4.0),
    "estuary_sealed": dict(KOREA_ESTUARY_POLICY="sealed"),
    "exclaves_masked": dict(ROK_EXCLAVE_POLICY="mask"),
    "buffer24": dict(KOREA_POST_ARRIVAL_MONTHS=24),
    # item 7: seeds = PPR_seed_points_dated.csv rows dated before 1 January 2015
    "holdout2015": dict(SEED_TABLE="auto", SEED_BEFORE="2015-01-01", SEED_FROM_OUTBREAK_SHP=True),
    # 11 Nov: the run of record again, drawing the manuscript figures for publication
    "pubfigs": dict(PUB_FIGS=True),
    # 11 Nov: habitat / accessible / barrier classes of the peninsular domain (Figure S2), then stop
    "habpeninsula": dict(HAB_PENINSULA_REPORT=True, PUB_FIGS=True),
}
V910_BATCH_MODES = {
    "none": (False, "off"), "scenarios": (True, "off"), "oat": (False, "oat"),
    "morris": (False, "morris"), "lhs": (False, "lhs"), "morris+lhs": (False, "morris+lhs"),
    "scenarios+oat": (True, "oat"), "all": (True, "morris+lhs"),
    "r2": (False, "morris+lhs"), "comparators": (False, "off"), "targeting": (False, "off"),
    "vaxdecomp": (False, "off"), "lhsr2": (False, "lhs"),
}


def v910_livestock_path():
    """19 Nov: which goat raster to read (see Config.LIVESTOCK_FILE)."""
    mode = str(Config.LIVESTOCK_FILE or "auto")
    filled = Config.ENV_DIR / Config.LIVESTOCK_ROKFILLED_NAME
    glw = Config.ENV_DIR / Config.LIVESTOCK_GLW_NAME
    if mode.lower() == "rokfilled":
        if not filled.exists():
            raise SystemExit(f"--goat-layer rokfilled: {filled} not found (copy it from "
                             f"PPR_Results\\v9_10_code into {Config.ENV_DIR})")
        return filled
    if mode.lower() == "glw":
        return glw
    if mode.lower() == "auto":
        return filled if filled.exists() else glw
    p = Path(mode)
    p = p if p.is_absolute() else Config.ENV_DIR / p
    if not p.exists():
        raise SystemExit(f"--goat-layer {mode}: {p} not found")
    return p


def v910_apply_preset(name, intro_month=None, out_suffix=None):
    name = str(name or "baseline").lower()
    if name not in V910_PRESETS:
        raise SystemExit(f"unknown preset '{name}' (choose from {sorted(V910_PRESETS)})")
    for k, v in V910_PRESETS[name].items():
        if k == "_gamma_mult":
            Config.GAMMA_RECOVERY = float(Config.GAMMA_RECOVERY) * float(v)
            Config.MONTHLY_BETA = np.asarray(Config.MONTHLY_BETA, float) * float(v)
        else:
            setattr(Config, k, v)
    if intro_month is not None:
        Config.KOREA_ONLY_INTRO_MONTH = int(intro_month)
    elif str(Config.KOREA_ONLY_INTRO_MONTH).lower() == "baseline":
        _fs = Config.OUT_DIR / f"{V910_TAG}_run_summary.json"
        try:
            Config.KOREA_ONLY_INTRO_MONTH = int(_v910_json.loads(_fs.read_text(encoding="utf-8"))
                                                ["baseline"]["intro_month"])
            logger.info(f"v9.10: DPRK introduction month {Config.KOREA_ONLY_INTRO_MONTH} taken from {_fs}")
        except Exception:
            Config.KOREA_ONLY_INTRO_MONTH = 18
            logger.info(f"v9.10: no 8-km run summary at {_fs}; DPRK introduction fixed at month 18 "
                        f"(the 6 Nov run). Run the baseline first, or pass --intro-month.")
    suffix = out_suffix if out_suffix is not None else ("" if name == "baseline" else f"_{name}")
    if suffix:
        Config.OUT_DIR = Config.OUT_DIR.with_name(Config.OUT_DIR.name + suffix)
        if Config.DL_DIR:
            Config.DL_DIR = Path(Config.DL_DIR).with_name(Path(Config.DL_DIR).name + suffix)
    return name


V910_PUB_FIG_FILES = {
    "PPR_v9.3_MaxEnt_suitability_outbreaks.png", "PPR_v9.3_habitat_classification.png",
    "PPR_v9.3_MDL_gates.png", "PPR_v9.3_PERCENT_SPREAD_CORRIDOR_GATED.png",
    "PPR_v9.3_annual_snapshots.png", "PPR_v9.3_NK_to_SK_progression.png",
    "PPR_v9.8_smcda_korea_panel_sk_reached.png", "PPR_v9.10_Figure9_vaccination.png",
}
V910_PUB_REPLACE = [
    ("% PPR spread (cumulative attack rate)", "Cumulative attack rate (%)"),     # 11 Nov: was doubled
    ("PPR MaxEnt suitability (base map)", "PPR suitability, $H_p$ (base map)"),
    ("Neighbouring land (regional basemap)", "Neighbouring land"),
    ("Study area (PPR_studyarea)", "Study area"),
    ("Ω_NK", r"$\Omega_{\mathrm{NK}}$"), ("DPRK %", "DPRK (%)"), ("ROK %", "ROK (%)"),
    ("NK % attack (warm)", "Cumulative attack rate, DPRK (%)"),
    ("SK % attack (blue-green)", "Cumulative attack rate, ROK (%)"),
    ("PPR MaxEnt suitability (0–1)", "PPR suitability, $H_p$ (0–1)"),
    ("MDL pass-through per cell (fraction of transport that crosses)", "Pass-through per cell, φ"),
    ("% PPR spread", "Cumulative attack rate (%)"),
    ("NORTH KOREA", "DPRK"), ("SOUTH KOREA", "ROK"), ("North Korea", "DPRK"), ("South Korea", "ROK"),
    ("South-Korean", "ROK"), ("Gorani", "Water deer"), ("gorani", "water deer"), ("GOR_EFF", "$H_e$"),
    ("Korea leg", "peninsular leg"), ("Stage 1", "Continental leg"), ("Stage 2", "Peninsular leg"),
    ("NK %", "DPRK (%)"), ("SK %", "ROK (%)"), ("MDL boundary", "Military Demarcation Line"),
    ("China introduction seed", "Introduction cells"),
]
V910_PUB_REPLACE_RE = [                       # 11 Nov: legend labels that carry file names or counts
    (r"Water bodies / water\.shp \(\d+\)", "Water bodies"),
    (r"Counties / sig\.shp \(\d+\)", "Sigungu boundaries (ROK)"),
    (r"DMZ water-deer corridor gates \(\d+; size ∝ weight\)", "Water-deer corridor gates (size ∝ weight)"),
    (r"Outbreak seed points \((\d+)\)", r"PPR outbreak locations (\1)"),
]
V910_PUB_ARROW_SCALE = {"PPR_v9.3_PERCENT_SPREAD_CORRIDOR_GATED.png": 3.0}


def v910_pubfig_mode():
    """--preset pubfigs: draw the manuscript figures without embedded titles, with the terms of the
    manuscript, a north arrow and 600 dpi. Only figure files are touched; GIF frames are not."""
    import matplotlib.figure as _mfig
    import matplotlib.axes as _maxes
    import matplotlib.text as _mtext
    _mfig.Figure.suptitle = lambda self, *a, **k: None
    _orig_title = _maxes.Axes.set_title

    import re as _re

    def _title(self, label=None, *a, **k):
        s = "" if label is None else str(label)
        if "\n" in s:
            # 11 Nov: a panel title whose first line is a date ("March  ·  Year 3") keeps that line
            _first = s.split("\n", 1)[0].strip()
            s = _first if _re.search(r"\bYear\s+\d+", _first) else ""
        if not s or len(s) > 48:
            return None
        return _orig_title(self, s, *a, **k)
    _maxes.Axes.set_title = _title
    _orig_save = _mfig.Figure.savefig

    def _save(self, fname, *a, **k):
        name = Path(str(fname)).name if isinstance(fname, (str, Path)) else ""
        if name in V910_PUB_FIG_FILES:
            for t in self.findobj(_mtext.Text):
                s = t.get_text()
                s2 = s
                for o, n in V910_PUB_REPLACE:
                    s2 = s2.replace(o, n)
                for o, n in V910_PUB_REPLACE_RE:
                    s2 = _re.sub(o, n, s2)
                if s2 != s:
                    t.set_text(s2)
            # 11 Nov: a label placed in data coordinates outside a cropped map (e.g. "ROK" below the
            # footprint of Figure 8) is not drawn, so it no longer stretches the saved canvas
            for ax in self.axes:
                _x0, _x1 = sorted(ax.get_xlim()); _y0, _y1 = sorted(ax.get_ylim())
                for t in list(ax.texts):
                    if type(t) is _mtext.Text and t.get_transform() is ax.transData:
                        _px, _py = t.get_position()
                        if not (_x0 <= _px <= _x1 and _y0 <= _py <= _y1):
                            t.set_visible(False)
            maps = [ax for ax in self.axes if ax.images and not ax.get_xticks().size]
            if maps:
                ax = maps[0]
                _ny = (0.06, 0.30, 0.18) if name == "PPR_v9.3_MDL_gates.png" else (0.94, 0.95, 0.83)
                _as = float(V910_PUB_ARROW_SCALE.get(name, 1.0))
                ax.annotate("N", xy=(_ny[0], _ny[1]), xytext=(_ny[0], _ny[2]), xycoords="axes fraction",
                            ha="center", va="center", fontsize=12 * _as, fontweight="bold", zorder=20,
                            arrowprops=dict(facecolor="black", edgecolor="black", width=3 * _as,
                                            headwidth=10 * _as, headlength=12 * _as))
            w, h = self.get_size_inches()
            k["dpi"] = int(min(float(Config.PUB_FIGS_DPI), 9000.0 / max(w, h, 1e-6)))
        return _orig_save(self, fname, *a, **k)
    _mfig.Figure.savefig = _save
    logger.info("v9.10: publication figures — no embedded titles, manuscript terms, north arrow, "
                f"{Config.PUB_FIGS_DPI} dpi (capped at 9000 px)")


def v910_cli(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="PPR model v9.10 (review re-analysis)")
    ap.add_argument("--preset", default=os.environ.get("PPR_PRESET", "baseline"),
                    help=f"one of {sorted(V910_PRESETS)}")
    ap.add_argument("--batch", default=os.environ.get("PPR_BATCH", None),
                    help=f"one of {sorted(V910_BATCH_MODES)}")
    ap.add_argument("--workers", type=int, default=None, help="parallel runs (default: CPUs − 1, max 8)")
    ap.add_argument("--zone-km", type=float, default=None,
                    help="gate-zone radius (km) for a --from-context batch (default: value in the context)")
    ap.add_argument("--from-context", default=None,
                    help="PPR_v9.10_korea_leg_context.npz: run the batch only")
    ap.add_argument("--out-suffix", default=None, help="suffix for the output folder")
    ap.add_argument("--goat-layer", default=None,
                    help="goat raster: auto (default) | rokfilled | glw | file name in the PPR_env folder")
    ap.add_argument("--dl-suffix", default=None,
                    help="extra suffix for the PPR_Results copy folder only (keeps earlier results)")
    ap.add_argument("--intro-month", type=int, default=None,
                    help="Korea-only mode: DPRK introduction month (skips most of the continental leg)")
    ap.add_argument("--seed-before", default=None,
                    help="hold-out run: seed only from outbreaks dated before this day (YYYY-MM-DD)")
    ap.add_argument("--seed-table", default=None,
                    help="dated seed table for the hold-out (default: PPR_seed_points_dated.csv)")
    ap.add_argument("--maxent", default=None,
                    help="MaxEnt map to use instead of the one found in the data folder "
                         "(e.g. refitted on the pre-cut-off outbreaks for the hold-out)")
    ap.add_argument("--point-cell", default=None, choices=["round", "floor"],
                    help="point locations on the grid: floor (default; the cell containing the point, as reported) "
                         "or round (nearest cell corner; the runs made before 26 Nov 2026)")
    args, _unknown = ap.parse_known_args(argv)
    if args.point_cell:
        Config.POINT_TO_CELL = str(args.point_cell)
    logger.info(f"v9.10: point locations (seeds, gates, Yalu entry) placed by '{Config.POINT_TO_CELL}'"
                + (" — the cell that contains each point" if Config.POINT_TO_CELL == "floor"
                   else " — the nearest cell corner (runs before 26 Nov 2026)"))
    _suffix = args.out_suffix
    if args.goat_layer:
        Config.LIVESTOCK_FILE = str(args.goat_layer)
    if args.seed_before and str(args.preset).lower() == "baseline" and _suffix is None:
        _suffix = "_holdout" + str(args.seed_before).replace("-", "")
    preset = v910_apply_preset(args.preset, args.intro_month, _suffix)
    if args.dl_suffix and Config.DL_DIR:
        Config.DL_DIR = Path(Config.DL_DIR).with_name(Path(Config.DL_DIR).name + str(args.dl_suffix))
    if getattr(Config, "PUB_FIGS", False):
        v910_pubfig_mode()
    if args.seed_before or args.seed_table:                          # item 7: hold-out seeding
        Config.SEED_TABLE = args.seed_table or Config.SEED_TABLE or "auto"
        Config.SEED_FROM_OUTBREAK_SHP = True
        if args.seed_before:
            Config.SEED_BEFORE = str(args.seed_before)
    if args.maxent:
        Config.MAXENT_TIF_PATH = str(args.maxent)
    if Config.SEED_TABLE:
        logger.info(f"v9.10: hold-out seeding — table '{Config.SEED_TABLE}', rows dated before "
                    f"{Config.SEED_BEFORE or '(no cut-off)'}; undated rows: {Config.SEED_UNDATED}")
    if args.batch:
        b = str(args.batch).lower()
        if b not in V910_BATCH_MODES:
            raise SystemExit(f"unknown --batch '{args.batch}' (choose from {sorted(V910_BATCH_MODES)})")
        Config.RUN_SCENARIOS, Config.RUN_SENSITIVITY = V910_BATCH_MODES[b]
        if b in ("r2", "comparators"):                              # round 2
            Config.RUN_R2 = True
            Config.SA_LHS_PARAMS = list(Config.SA_LHS_PARAMS_R2)
        if b == "targeting":                                         # round 2, F1: Table S9 only
            Config.RUN_R2 = True
            Config.R2_SKIP_COMPARATORS = True
        if b == "lhsr2":                                             # round 2, FC2: LHS only
            Config.SA_LHS_PARAMS = list(Config.SA_LHS_PARAMS_R2)
        if b == "vaxdecomp":                                         # round 2, FC1: V2 decomposition
            Config.RUN_R2 = True
            Config.R2_VAX_DECOMP = True
            Config.COMPARATOR_LIST = list(Config.VAX_DECOMP_LIST)
            Config.R2_COMPARATOR_STEM = "R2_vaccination_decomposition"
    if args.workers:
        Config.SA_WORKERS = int(args.workers)
    logger.info(f"v9.10: preset '{preset}', batch scenarios={Config.RUN_SCENARIOS}, "
                f"sensitivity='{Config.RUN_SENSITIVITY}', output {Config.OUT_DIR}")
    if args.from_context:
        if preset != "baseline":
            raise SystemExit("--from-context reuses the parameters stored in the context file; a preset "
                             "cannot change them. Run the preset as a full run instead.")
        ctx = load_korea_ctx(args.from_context)
        if getattr(args, "zone_km", None):                           # FC2: gate-zone radius
            logger.info(f"v9.10: gate-zone radius {float(ctx['probe_km']):g} km -> {float(args.zone_km):g} km")
            ctx["probe_km"] = float(args.zone_km)
            Config.FENCE_WAIT_PROBE_KM = float(args.zone_km)
        logger.info(f"v9.10: context goat layer: {ctx.get('livestock_file') or 'GLW4GTS_1km_5179.tif (pre-19 Nov context)'}")
        _cdir = Path(args.from_context).parent
        out = _cdir if args.out_suffix is None else _cdir.with_name(_cdir.name + args.out_suffix)
        logger.info(f"v9.10: batch from context {args.from_context} (ROK split: {ctx.get('rok_split')})")
        if not (Config.RUN_SCENARIOS or Config.RUN_R2
                or str(Config.RUN_SENSITIVITY).lower() not in ("off", "none")):
            Config.RUN_SCENARIOS = True
        v910_run_batch(ctx, scenarios=bool(Config.RUN_SCENARIOS), sa_mode=str(Config.RUN_SENSITIVITY),
                       out_dir=out, dl_dir=Config.DL_DIR, workers=Config.SA_WORKERS, logger_=logger)
        return
    if Config.DL_DIR:
        try:
            Path(Config.DL_DIR).mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logger.warning(f"v9.10: copy folder {Config.DL_DIR} not available ({e})")
    main()


if __name__ == "__main__":
    try:
        v910_cli()
    except SystemExit:
        raise
    except Exception as e:
        logger.error(f"\n✗ ERROR: {e}", exc_info=True)
        raise

