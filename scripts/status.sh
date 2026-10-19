#!/usr/bin/env bash
# =============================================================================
# status.sh
# One line per force field and replica: build / equilibration / production
# state, production progress in ns and percent, performance, and whether the
# analysis is up to date. Reads md.log (last "Step ... Time" block), so it
# works while a run is in progress and does not touch the simulation.
#
# Usage:
#   bash scripts/status.sh             # table
#   watch -n 300 bash scripts/status.sh
#   tail -f logs/production.log        # event log written by 02_run_production.sh
# =============================================================================

source "$(dirname "$0")/lib.sh"

total_ps=$(fcalc "$PROD_NS * 1000")

# Last simulated time (ps) in an mdrun log: the "Step Time" header is
# followed by a line with the step number and the time.
md_time_ps() {
    [[ -s "$1" ]] || { echo 0; return; }
    awk '/^ +Step +Time/ { getline; t = $2 } END { printf "%s", (t == "" ? 0 : t) }' "$1"
}

printf '%-11s %-4s %-6s %-6s %-11s %8s %5s %9s  %s\n' \
    "forcefield" "rep" "build" "equil" "production" "ns" "%" "ns/day" "analysis"
for ff in "${FF_LABELS[@]}"; do
    B="$RUN_DIR/$ff/build"
    build="-"; [[ -s "$B/em.gro" ]] && build="ok"
    for (( r = 1; r <= N_REPLICAS; r++ )); do
        R="$RUN_DIR/$ff/rep$r"; A="$ANA_DIR/$ff/rep$r"
        equil="-"
        [[ -s "$R/nvt.gro" ]] && equil="nvt"
        [[ -s "$R/npt.gro" ]] && equil="ok"
        prod="pending"; ns="-"; pct="-"; perf="-"
        if [[ -s "$R/md.gro" ]]; then
            prod="finished"; t=$total_ps
        elif [[ -s "$R/md.log" ]]; then
            t=$(md_time_ps "$R/md.log")
            # running if md.log was modified in the last 10 minutes
            if [[ -n $(find "$R/md.log" -mmin -10 2>/dev/null) ]]; then prod="running"; else prod="stopped"; fi
        elif [[ -s "$R/md.tpr" ]]; then
            prod="ready"; t=0
        else
            t=0
        fi
        if [[ "$prod" != pending && "$prod" != ready ]]; then
            ns=$(awk "BEGIN{printf \"%.1f\", $t / 1000}")
            pct=$(awk "BEGIN{printf \"%.0f\", 100 * $t / $total_ps}")
            perf=$( { grep -E '^Performance:' "$R/md.log" 2>/dev/null || true; } | tail -n 1 | awk '{print $2}')
            [[ -z "$perf" ]] && perf=$( { grep -oE 'will finish [^,]*' "$R/md_stdout.log" 2>/dev/null || true; } | tail -n 1 | sed 's/will finish //')
            perf=${perf:--}
        fi
        ana="-"
        if [[ -s "$A/meta.txt" ]]; then
            if [[ -s "$R/md.xtc" && "$R/md.xtc" -nt "$A/md_center.xtc" ]]; then ana="stale"; else ana="ok"; fi
        fi
        printf '%-11s %-4s %-6s %-6s %-11s %8s %5s %9s  %s\n' "$ff" "$r" "$build" "$equil" "$prod" "$ns" "$pct" "$perf" "$ana"
    done
done
[[ -s "$LOG_DIR/production.log" ]] && { echo; echo "last events ($LOG_DIR/production.log):"; tail -n 5 "$LOG_DIR/production.log"; }
exit 0
