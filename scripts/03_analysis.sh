#!/usr/bin/env bash
# =============================================================================
# 03_analysis.sh
# Trajectory analyses for every force field and replica:
#   RMSD (backbone), RMSF (C-alpha), radius of gyration, intramolecular
#   H-bonds, protein-water H-bonds, SASA, residue contact map (gmx mdmat),
#   secondary structure per residue and frame (gmx dssp), plus equilibration
#   quality checks (EM, NVT, NPT, production energies).
# Results: analysis/<ff>/rep<N>/*.xvg|*.xpm  (plots: 04_plot_comparison.py)
#
# Comparability between the all-atom models (CHARMM36m, AMBER) and the
# united-atom ones (GROMOS 53A6, 54A7) is handled explicitly:
#   - RMSD/RMSF use backbone / C-alpha atoms, present in every model.
#   - SASA is reported twice: "native" (all atoms each model has) and
#     "heavy" (hydrogens removed in all). Only "heavy" is directly
#     comparable, because GROMOS lacks aliphatic hydrogens and gmx sasa uses
#     element radii, which underestimates united-atom CHn surfaces.
#   - The contact map uses heavy atoms only, for the same reason (explicit
#     hydrogens shorten residue-residue minimum distances in CHARMM/AMBER).
#   - H-bonds involve polar hydrogens, which are explicit in every model.
# Every replica is analysed on its own; 04_plot_comparison.py compares the
# replicas within a force field and the force fields with each other.
#
# Usage:
#   bash scripts/03_analysis.sh               # all force fields and replicas
#   bash scripts/03_analysis.sh gromos54a7    # one force field
#   FORCE=1 bash scripts/03_analysis.sh       # recompute existing outputs
# =============================================================================

source "$(dirname "$0")/lib.sh"
check_gmx
parse_ff_args "$@"
FORCE="${FORCE:-0}"
STALE=0   # set per replica when the trajectory is newer than the analysis
must_run() { [[ "$FORCE" == 1 || "$STALE" == 1 || ! -s "$1" ]]; }

# GROMACS 2024 rewrote gyrate and hbond; the original implementations are
# kept as *-legacy. Using the legacy tools when available keeps the algorithm
# and output identical across GROMACS 2023, 2024 and 2025.
# (help text is captured first: piping into grep -q under pipefail would make
# the test fail spuriously when gmx gets SIGPIPE)
_help_gyrate=$("$GMX" -quiet gyrate -h 2>&1 || true)
_help_hbond=$("$GMX" -quiet hbond -h 2>&1 || true)
if gmx_has_tool gyrate-legacy; then GYRATE=gyrate-legacy; GYRATE_NEW=0
elif grep -qE -- '^ *-sel' <<<"$_help_gyrate"; then GYRATE=gyrate; GYRATE_NEW=1
else GYRATE=gyrate; GYRATE_NEW=0; fi

if gmx_has_tool hbond-legacy; then HBOND=hbond-legacy; HBOND_NEW=0
elif grep -qE -- '-(\[no\])?contact' <<<"$_help_hbond"; then HBOND=hbond; HBOND_NEW=0
else HBOND=hbond; HBOND_NEW=1; fi

log "GROMACS $GMX_VERSION | gyrate: $GYRATE | hbond: $HBOND"

# Last time (ps) present in a trajectory.
traj_last_time() {
    # gmx check prints progress with carriage returns, so search the text
    "$GMX" -quiet check -f "$1" 2>&1 | tr '\r' '\n' |
        grep -oE 'Last frame +-?[0-9]+ +time +[0-9.]+' | awk '{print $NF}' | tail -n 1
}

# -----------------------------------------------------------------------------
# Equilibration quality checks (one set per force field and replica)
# -----------------------------------------------------------------------------
energy_terms() {
    # energy_terms <logfile> <edr> <out.xvg> <term...>
    local L=$1 edr=$2 out=$3; shift 3
    [[ -s "$edr" ]] || { warn "missing $edr, skipping"; return 0; }
    local input; input="$(printf '%s\\n' "$@")\n"
    gmx_run_input "$L" "$input" energy -f "$edr" -o "$out"
}

# -----------------------------------------------------------------------------
# Main analysis of one replica
# -----------------------------------------------------------------------------
analyse_replica() {
    local ff=$1 r=$2
    local B="$RUN_DIR/$ff/build" R="$RUN_DIR/$ff/rep$r" A="$ANA_DIR/$ff/rep$r"
    mkdir -p "$A"
    local L="$REPO_ROOT/$A/analysis.log"
    : >"$L"

    need_file "$R/md.tpr"
    need_file "$R/md.xtc"
    [[ -s "$R/md.gro" ]] || warn "[$ff rep$r] md.gro missing: production not finished, analysing the partial trajectory"

    # Time window for averages; fall back to the whole run if it is too short
    local tlast skip_ps
    tlast=$(traj_last_time "$R/md.xtc"); tlast=${tlast:-0}
    skip_ps=$(fcalc "$SKIP_NS * 1000")
    if awk "BEGIN{exit !($tlast <= $skip_ps)}"; then
        warn "[$ff rep$r] trajectory ends at ${tlast} ps, shorter than SKIP_NS=${SKIP_NS} ns; averages use the whole run"
        skip_ps=0
    fi
    local dssp_dt="$DSSP_DT_PS"; [[ "$dssp_dt" == 0 ]] && dssp_dt="$XTC_PS"
    printf 'force_field=%s\ndisplay_name=%s\nff_dir=%s\nwater=%s\nfamily=%s\nreplica=%s\ngromacs=%s\nlast_time_ps=%s\nskip_ps=%s\nmdmat_cutoff_nm=%s\ndssp_dt_ps=%s\n' \
        "$ff" "$(ff_display_name "$ff")" "$(ff_dir_name "$ff")" "$(ff_water "$ff")" "$(ff_family "$ff")" \
        "$r" "$GMX_VERSION" "$tlast" "$skip_ps" "$MDMAT_CUTOFF" "$dssp_dt" >"$A/meta.txt"
    log "[$ff rep$r] trajectory length ${tlast} ps; averages from ${skip_ps} ps"

    local tpr="$R/md.tpr" ndx="$A/index.ndx" trj="$A/md_center.xtc"

    # A previous analysis of a partial (still running) trajectory is outdated
    # once md.xtc grows: recompute everything for this replica.
    STALE=0
    if [[ -s "$trj" && "$R/md.xtc" -nt "$trj" ]]; then
        log "[$ff rep$r] trajectory changed since the last analysis, recomputing"
        STALE=1
    fi

    # ---- Index: default groups + Protein_heavy (protein without hydrogens) --
    if must_run "$ndx"; then
        gmx_run_input "$L" "q\n" make_ndx -f "$tpr" -o "$ndx"
        gmx_run "$L" select -s "$tpr" -on "$A/heavy.ndx" \
            -select '"Protein_heavy" group "Protein" and not name "H*"' -n "$ndx"
        cat "$A/heavy.ndx" >>"$ndx"
        rm -f "$A/heavy.ndx"
    fi

    # ---- PBC treatment: whole molecules, protein centred, compact cell ------
    if must_run "$trj"; then
        log "[$ff rep$r] trjconv (center, whole molecules)"
        gmx_run_input "$L" "Protein\nSystem\n" trjconv -s "$tpr" -f "$R/md.xtc" -n "$ndx" \
            -o "$trj" -pbc mol -center -ur compact
        # Protein-only, fitted trajectory and first frame for visualization
        gmx_run_input "$L" "Backbone\nProtein\n" trjconv -s "$tpr" -f "$trj" -n "$ndx" \
            -o "$A/protein_fit.xtc" -fit rot+trans
        gmx_run_input "$L" "Protein\n" trjconv -s "$tpr" -f "$trj" -n "$ndx" \
            -o "$A/protein_start.pdb" -dump 0
    fi

    # ---- RMSD: backbone vs energy-minimized structure -----------------------
    if must_run "$A/rmsd_backbone.xvg"; then
        log "[$ff rep$r] RMSD"
        gmx_run_input "$L" "Backbone\nBackbone\n" rms -s "$B/em.tpr" -f "$trj" -n "$ndx" \
            -o "$A/rmsd_backbone.xvg" -tu ns
    fi

    # ---- RMSF: C-alpha, per residue, equilibrated part ----------------------
    if must_run "$A/rmsf_ca.xvg"; then
        log "[$ff rep$r] RMSF"
        gmx_run_input "$L" "C-alpha\n" rmsf -s "$tpr" -f "$trj" -n "$ndx" -b "$skip_ps" \
            -res -o "$A/rmsf_ca.xvg"
    fi

    # ---- Radius of gyration (mass weighted, all protein atoms) --------------
    # Mass-weighted Rg is comparable between AA and UA: united atoms carry
    # the mass of their hydrogens.
    if must_run "$A/gyrate.xvg"; then
        log "[$ff rep$r] radius of gyration"
        if [[ $GYRATE_NEW == 1 ]]; then
            gmx_run "$L" "$GYRATE" -s "$tpr" -f "$trj" -n "$ndx" -sel 'group "Protein"' \
                -o "$A/gyrate.xvg"
        else
            gmx_run_input "$L" "Protein\n" "$GYRATE" -s "$tpr" -f "$trj" -n "$ndx" -o "$A/gyrate.xvg"
        fi
    fi

    # ---- Hydrogen bonds (r_DA <= 0.35 nm, angle HDA <= 30 deg) --------------
    if must_run "$A/hbnum_intra.xvg"; then
        log "[$ff rep$r] H-bonds protein-protein"
        if [[ $HBOND_NEW == 1 ]]; then
            gmx_run "$L" "$HBOND" -s "$tpr" -f "$trj" -n "$ndx" -r 'group "Protein"' \
                -num "$A/hbnum_intra.xvg" -tu ns
        else
            gmx_run_input "$L" "Protein\nProtein\n" "$HBOND" -s "$tpr" -f "$trj" -n "$ndx" \
                -num "$A/hbnum_intra.xvg" -tu ns
        fi
    fi
    if must_run "$A/hbnum_water.xvg"; then
        log "[$ff rep$r] H-bonds protein-water"
        local dtw=""
        [[ "$HB_WATER_DT_PS" != 0 ]] && dtw="-dt $HB_WATER_DT_PS"
        if [[ $HBOND_NEW == 1 ]]; then
            gmx_run "$L" "$HBOND" -s "$tpr" -f "$trj" -n "$ndx" -r 'group "Protein"' \
                -t 'group "SOL"' -num "$A/hbnum_water.xvg" -tu ns $dtw
        else
            gmx_run_input "$L" "Protein\nSOL\n" "$HBOND" -s "$tpr" -f "$trj" -n "$ndx" \
                -num "$A/hbnum_water.xvg" -tu ns $dtw
        fi
    fi

    # ---- SASA ----------------------------------------------------------------
    if must_run "$A/sasa_heavy.xvg"; then
        log "[$ff rep$r] SASA"
        gmx_run "$L" sasa -s "$tpr" -f "$trj" -n "$ndx" -surface 'group "Protein"' \
            -o "$A/sasa_native.xvg" -tu ns
        gmx_run "$L" sasa -s "$tpr" -f "$trj" -n "$ndx" -surface 'group "Protein_heavy"' \
            -o "$A/sasa_heavy.xvg" -tu ns
        # per-residue average over the equilibrated part
        gmx_run "$L" sasa -s "$tpr" -f "$trj" -n "$ndx" -surface 'group "Protein_heavy"' \
            -b "$skip_ps" -o "$A/sasa_heavy_window.xvg" -or "$A/sasa_heavy_res.xvg"
    fi

    # ---- Residue contact map (minimum heavy-atom distance) ------------------
    if must_run "$A/mdmat_mean.xpm"; then
        log "[$ff rep$r] contact map (mdmat)"
        local dtm=""
        [[ "$MDMAT_DT_PS" != 0 ]] && dtm="-dt $MDMAT_DT_PS"
        gmx_run_input "$L" "Protein_heavy\n" mdmat -s "$tpr" -f "$trj" -n "$ndx" \
            -b "$skip_ps" $dtm -t "$MDMAT_CUTOFF" -nlevels 100 -mean "$A/mdmat_mean.xpm" -no "$A/mdmat_ncontacts.xvg"
    fi

    # ---- Secondary structure (gmx dssp, native since GROMACS 2023) ----------
    # -hmode dssp rebuilds the amide hydrogen from the backbone geometry, as
    # the original DSSP does, so the assignment does not depend on how each
    # force field places its explicit hydrogens (all-atom vs united-atom).
    # One line of letters per frame (H/G/I helices, E/B strands, T/S/P/~
    # others); 04_plot_comparison.py reduces them to helix / strand / coil.
    if must_run "$A/dssp.dat"; then
        log "[$ff rep$r] DSSP"
        local dtd=""
        [[ "$DSSP_DT_PS" != 0 ]] && dtd="-dt $DSSP_DT_PS"
        gmx_run "$L" dssp -s "$tpr" -f "$trj" -n "$ndx" -sel 'group "Protein"' \
            -hmode dssp -o "$A/dssp.dat" $dtd
    fi
    # Reference assignment of the input (crystal) structure, once per analysis
    if [[ ! -s "$ANA_DIR/dssp_input.dat" && -s "$B/protein_clean.pdb" ]]; then
        gmx_run "$L" dssp -s "$B/protein_clean.pdb" -f "$B/protein_clean.pdb" \
            -sel 'group "Protein"' -hmode dssp -o "$ANA_DIR/dssp_input.dat"
    fi

    # ---- Equilibration and production quality checks ------------------------
    if must_run "$A/qc_md_thermo.xvg"; then
        log "[$ff rep$r] energy checks"
        [[ -s "$ANA_DIR/$ff/qc_em_potential.xvg" ]] || \
            energy_terms "$L" "$B/em.edr" "$ANA_DIR/$ff/qc_em_potential.xvg" Potential
        energy_terms "$L" "$R/nvt.edr" "$A/qc_nvt_temperature.xvg" Temperature
        energy_terms "$L" "$R/npt.edr" "$A/qc_npt_pressure.xvg" Pressure
        energy_terms "$L" "$R/npt.edr" "$A/qc_npt_density.xvg" Density
        energy_terms "$L" "$R/md.edr" "$A/qc_md_thermo.xvg" Temperature Pressure Density Potential
    fi

    log "[$ff rep$r] done -> $A"
}

for ff in "${SELECTED_FF[@]}"; do
    for (( r = 1; r <= N_REPLICAS; r++ )); do
        if [[ -s "$RUN_DIR/$ff/rep$r/md.xtc" ]]; then
            analyse_replica "$ff" "$r"
        else
            warn "[$ff rep$r] no md.xtc, skipping"
        fi
    done
done
log "done. Next: python3 scripts/04_plot_comparison.py"
