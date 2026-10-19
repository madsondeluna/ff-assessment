#!/usr/bin/env bash
# =============================================================================
# 02_run_production.sh
# Runs (or resumes) the production simulations prepared by 01_prepare.sh,
# strictly one at a time: the next mdrun starts only when the previous one
# has finished (or stopped at MAXH). Rerunning this script continues every
# unfinished run from its md.cpt, so an interrupted campaign picks up where
# it stopped.
#
# Order (RUN_ORDER in config.sh):
#   replica     rep1 of every force field, then rep2, ... (default: a
#               one-replica comparison is available early)
#   forcefield  all replicas of one force field, then the next
#
# Progress log: logs/production.log (one line per start/finish, with ns/day)
# Per-run mdrun output: runs/<ff>/rep<N>/md_stdout.log and md.log
# Live status of every run: bash scripts/status.sh
#
# Usage:
#   bash scripts/02_run_production.sh                 # all force fields, all replicas
#   bash scripts/02_run_production.sh charmm36m       # one force field
#   REPS="2 3" bash scripts/02_run_production.sh      # selected replicas
#   MAXH=23.5 bash scripts/02_run_production.sh       # stop cleanly before a queue limit
#   nohup bash scripts/02_run_production.sh > logs/production_stdout.log 2>&1 &
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
mkdir -p "$LOG_DIR"
PLOG="$LOG_DIR/production.log"

# plog <message>: timestamped line to the terminal and to the progress log
plog() { log "$*"; printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >>"$PLOG"; }

# Build the sequence of (ff, replica) pairs
QUEUE=()
if [[ "$RUN_ORDER" == forcefield ]]; then
    for ff in "${SELECTED_FF[@]}"; do for r in $REPS; do QUEUE+=("$ff:$r"); done; done
else
    for r in $REPS; do for ff in "${SELECTED_FF[@]}"; do QUEUE+=("$ff:$r"); done; done
fi

plog "campaign: ${#QUEUE[@]} run(s) of ${PROD_NS} ns, order=$RUN_ORDER, threads=$NT, GROMACS $GMX_VERSION${MAXH:+, maxh=$MAXH}"
for pair in "${QUEUE[@]}"; do
    ff=${pair%%:*}; r=${pair##*:}
    R="$RUN_DIR/$ff/rep$r"
    need_file "$R/md.tpr"
    if [[ -s "$R/md.gro" ]]; then
        log "[$ff rep$r] already finished (md.gro present), skipping"
        continue
    fi

    cpi=""
    [[ -s "$R/md.cpt" ]] && cpi="-cpi md.cpt"
    maxh=""
    [[ -n "$MAXH" ]] && maxh="-maxh $MAXH"

    plog "[$ff rep$r] ${cpi:+resuming}${cpi:-started} ($(ff_display_name "$ff"), ${PROD_NS} ns)"
    t0=$(date +%s)
    ( cd "$R" && gmx_run "$REPO_ROOT/$R/md_stdout.log" mdrun -deffnm md $cpi $maxh \
        $THREAD_FLAGS $MDRUN_FLAGS )
    elapsed=$(( $(date +%s) - t0 ))

    if [[ -s "$R/md.gro" ]]; then
        perf=$( { grep -E '^Performance:' "$R/md.log" || true; } | tail -n 1 | awk '{print $2}')
        plog "[$ff rep$r] finished in $(( elapsed / 3600 )) h $(( elapsed % 3600 / 60 )) min (${perf:-?} ns/day)"
    else
        plog "[$ff rep$r] stopped before the end after $(( elapsed / 60 )) min (time limit); rerun to continue"
    fi
done
plog "campaign done: all requested runs finished or stopped. Next: bash scripts/03_analysis.sh"
