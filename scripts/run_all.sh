#!/usr/bin/env bash
# =============================================================================
# run_all.sh
# The whole campaign, strictly in sequence, in one detached process:
#   get_structure -> 01_prepare -> 02_run_production -> 03_analysis -> 04_plot
# Each stage starts only when the previous one has finished, and every line
# of output is copied to logs/pipeline_<timestamp>.log. Rerunning the script
# skips what is already done (each stage checks its own outputs).
#
# Usage:
#   nohup bash scripts/run_all.sh > /dev/null 2>&1 &
#   tail -f logs/pipeline_*.log          # follow
#   bash scripts/status.sh               # table of every run
#   PROD_NS=10 N_REPLICAS=1 bash scripts/run_all.sh    # smoke test
# =============================================================================

source "$(dirname "$0")/lib.sh"
mkdir -p "$LOG_DIR"
PIPE_LOG="$LOG_DIR/pipeline_$(date '+%Y%m%d_%H%M%S').log"
ln -sfn "$(basename "$PIPE_LOG")" "$LOG_DIR/pipeline_latest.log"

# Everything below goes to the terminal and to the log
exec > >(tee -a "$PIPE_LOG") 2>&1

stage() {
    local name=$1; shift
    log "===== $name: start ====="
    local t0; t0=$(date +%s)
    if "$@"; then
        log "===== $name: done in $(( ($(date +%s) - t0) / 60 )) min ====="
    else
        die "===== $name: FAILED (see $PIPE_LOG) ====="
    fi
}

log "pipeline log: $PIPE_LOG"
log "force fields: ${FF_LABELS[*]} | replicas: $N_REPLICAS | production: ${PROD_NS} ns | threads: $NT"
stage "prepare"    bash "$REPO_ROOT/scripts/01_prepare.sh" "$@"
stage "production" bash "$REPO_ROOT/scripts/02_run_production.sh" "$@"
stage "analysis"   bash "$REPO_ROOT/scripts/03_analysis.sh" "$@"
stage "figures"    python3 "$REPO_ROOT/scripts/04_plot_comparison.py"
log "pipeline finished. Figures: $FIG_DIR/"
