#!/usr/bin/env bash
# Reproduces every model run of the paper (8-, 4- and 2-km grids, scenarios, sensitivity analysis, comparator hosts,
# surveillance rules, publication figures, hold-out run and external test). Several hours in total; a step whose
# outputs exist is skipped, so the script can be restarted after an interruption.
#
#   conda activate ppr_rad
#   bash code/reproduce_all.sh /path/to/input_folder        (run from the folder that should hold results/)
set -u
DATA=${1:?give the input folder (see INPUT_DATA.md)}
HERE=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-python}
RUN="$PY $HERE/run_model.py --data-root $DATA --copy-dir $PWD/copies/v9_10"
R=results/ppr_v9_10_review
C8=${R}_g_pcfloor/PPR_v9.10_korea_leg_context.npz
C4=${R}_g_pcfloor_grid4km/PPR_v9.10_korea_leg_context.npz
C2=${R}_g_pcfloor_grid2km/PPR_v9.10_korea_leg_context.npz
D=copies/v9_10
mkdir -p logs

step() {   # step NAME FILE-THAT-MARKS-DONE model-arguments...
    local name=$1 done=$2; shift 2
    if [ -s "$done" ]; then echo "SKIP  $name"; return 0; fi
    echo "START $name $(date '+%F %T')"
    $RUN "$@" 2>&1 | tee "logs/$name.log" > /dev/null
    [ -s "$done" ] && echo "DONE  $name" || { echo "FAILED $name (see logs/$name.log)"; exit 1; }
}

step base8   "${D}_g_pcfloor/PPR_v9.10_run_summary.json"        --out-suffix _g_pcfloor --goat-layer rokfilled
step batch8  "${D}_g_pcfloor/PPR_v9.10_SA_lhs_runs.csv"          --from-context "$C8" --batch all --dl-suffix _g_pcfloor
step r2      "${D}_r2_g_pcfloor/PPR_v9.10_SA_lhs_runs.csv"       --from-context "$C8" --batch r2 --out-suffix _r2 --dl-suffix _g_pcfloor
step f1ci8   "${D}_f1ci_g_pcfloor/PPR_v9.10_R2_targeting_rho_ci.csv"  --from-context "$C8" --batch targeting --out-suffix _f1ci --dl-suffix _g_pcfloor
step vd      "${D}_vd_g_pcfloor/PPR_v9.10_R2_vaccination_decomposition.csv" --from-context "$C8" --batch vaxdecomp --out-suffix _vd --dl-suffix _g_pcfloor
for z in 12 36; do
  step z$z "${D}_z${z}_g_pcfloor/PPR_v9.10_SA_lhs_runs.csv" --from-context "$C8" --batch lhsr2 --zone-km $z --out-suffix _z$z --dl-suffix _g_pcfloor
done
step grid4km "${D}_g_pcfloor_grid4km/PPR_v9.10_run_summary.json" --preset grid4km --out-suffix _g_pcfloor_grid4km --goat-layer rokfilled
step r2_4km  "${D}_r2_4km_g_pcfloor/PPR_v9.10_Table4_scenarios.csv"   --from-context "$C4" --batch scenarios --out-suffix _r2_4km --dl-suffix _g_pcfloor
step f1ci4   "${D}_f1ci4_g_pcfloor/PPR_v9.10_R2_targeting_rho_ci.csv" --from-context "$C4" --batch targeting --out-suffix _f1ci4 --dl-suffix _g_pcfloor
step pubfigs "${D}_g_pcfloor_pubfigs/PPR_v9.10_run_summary.json" --preset pubfigs --out-suffix _g_pcfloor_pubfigs --goat-layer rokfilled
step grid2km "${D}_g_pcfloor_grid2km/PPR_v9.10_run_summary.json" --preset grid2km --out-suffix _g_pcfloor_grid2km --goat-layer rokfilled
step f1ci2   "${D}_f1ci2_g_pcfloor/PPR_v9.10_R2_targeting_rho_ci.csv" --from-context "$C2" --batch targeting --out-suffix _f1ci2 --dl-suffix _g_pcfloor
step holdout "${D}_holdout2015_pcfloor/PPR_v9.10_run_summary.json" --preset holdout2015 --out-suffix _holdout2015_pcfloor \
             --seed-table "$HERE/../data/PPR_seed_points_dated.csv"

# external test against the outbreaks of 2015-2025 (Section 2.5.6, Table S8) and the audit of its counts and AUC
OUT=$PWD
( cd "$HERE/../data" && $PY "$HERE/ppr_external_test.py" --points EMPRESi_PPR_China_2007_2026_clean.csv \
      --results-dir "$OUT/${R}_holdout2015_pcfloor" ) 2>&1 | tee logs/external_test.log
$PY "$HERE/ppr_auc_audit.py" "$OUT/${R}_holdout2015_pcfloor" 2>&1 | tee logs/auc_audit.log
echo "all steps finished $(date '+%F %T')"
