#!/usr/bin/env bash
# =============================================================================
# 02_run_production.sh
# Runs (or resumes) the production simulations prepared by 01_prepare.sh.
# Rerunning this script continues from md.cpt, so an interrupted 100 ns run
# picks up where it stopped.
#
# Runs are sequential (one mdrun at a time, all threads): on a single
# workstation this is faster than running replicas side by side.
#
# Usage:
#   bash scripts/02_run_production.sh                 # all force fields, all replicas
#   bash scripts/02_run_production.sh charmm36m       # one force field
#   REPS="2 3" bash scripts/02_run_production.sh      # selected replicas
#   MAXH=23.5 bash scripts/02_run_production.sh       # stop cleanly before a queue limit
#   nohup bash scripts/02_run_production.sh > production.log 2>&1 &   # detached
# =============================================================================

source "$(dirname "$0")/lib.sh"
check_gmx
parse_ff_args "$@"

if "$GMX" -quiet --version 2>/dev/null | grep -Eq '^MPI library: *MPI'; then
    THREAD_FLAGS="-ntomp $NT"
else
    THREAD_FLAGS="-ntmpi 1 -ntomp $NT"
fi

REPS="${REPS:-$(seq -s ' ' 1 "$N_REPLICAS")}"
MAXH="${MAXH:-}"

for ff in "${SELECTED_FF[@]}"; do
    for r in $REPS; do
        R="$RUN_DIR/$ff/rep$r"
        need_file "$R/md.tpr"
        if [[ -s "$R/md.gro" ]]; then
            log "[$ff rep$r] already finished (md.gro present), skipping"
            continue
        fi

        cpi=""
        [[ -s "$R/md.cpt" ]] && { cpi="-cpi md.cpt"; log "[$ff rep$r] resuming from md.cpt"; }
        maxh=""
        [[ -n "$MAXH" ]] && maxh="-maxh $MAXH"

        log "[$ff rep$r] production started ($(ff_display_name "$ff"), ${PROD_NS} ns)"
        ( cd "$R" && gmx_run "$REPO_ROOT/$R/md_stdout.log" mdrun -deffnm md $cpi $maxh \
            $THREAD_FLAGS $MDRUN_FLAGS )

        if [[ -s "$R/md.gro" ]]; then
            perf=$(grep -E '^Performance:' "$R/md.log" | tail -n 1 | awk '{print $2}')
            log "[$ff rep$r] finished (${perf:-?} ns/day)"
        else
            log "[$ff rep$r] stopped before the end (time limit); rerun to continue"
        fi
    done
done
