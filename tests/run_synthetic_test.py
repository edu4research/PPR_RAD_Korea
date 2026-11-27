"""End-to-end test on a small synthetic China-Korea world (no restricted data needed).

    cd tests && python run_synthetic_test.py          (about 5-10 minutes)

Runs the 8-km baseline of the model on tests/world (made by make_synth_world.py; same grid conventions, the six real
gate coordinates, an unfenced estuary, ROK sigungu polygons and outbreak points in China) and checks the timeline
against the values below, produced with this release on Python 3.12 / NumPy 1.26 / GDAL 3.8.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code"))
import PPR_model_v9_10 as M  # noqa: E402

EXPECTED = {"intro_month": 14, "first_rok": 95, "first_location": "Goseong_Coast", "est": 110,
            "transit_months": 81, "fence_dprk_first": 88, "first_rok_beyond_band": 101, "wait_at_fence": 13}

C = M.Config
RUN = HERE / "run"
RUN.mkdir(exist_ok=True)
C.CWD = RUN
C.DATA_ROOT = HERE / "world"
C.ENV_DIR = C.DATA_ROOT / "PPR_env"
C.OUT_DIR = RUN / "results" / "ppr_v9_10_review"
C.DL_DIR = RUN / "copies" / "v9_10"
C.BASE_TARGET_MAX_PX = 206              # the synthetic world is small
C.REFERENCE_TARGET_MAX_PX = 206
M.v910_cli(["--out-suffix", "_test", "--goat-layer", "rokfilled"])

summ = json.loads((RUN / "copies" / "v9_10_test" / "PPR_v9.10_run_summary.json").read_text())["baseline"]
bad = {k: (v, summ.get(k)) for k, v in EXPECTED.items() if summ.get(k) != v}
print("\nexpected vs obtained:", {k: (v, summ.get(k)) for k, v in EXPECTED.items()})
print("TEST PASSED" if not bad else f"TEST FAILED: {bad}")
sys.exit(1 if bad else 0)
