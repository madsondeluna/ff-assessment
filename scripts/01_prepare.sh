#!/usr/bin/env bash
# =============================================================================
# 01_prepare.sh
# Builds and equilibrates the test protein (ubiquitin, 1UBQ, by default; see
# PDB_ID / INPUT_PDB in config.sh) in water with every force field in FF_LABELS
# (CHARMM36m, AMBER, GROMOS 53A6, GROMOS 54A7 by default):
#   clean PDB -> pdb2gmx -> box -> solvate -> ions -> EM      (runs/<ff>/build)
#   NVT (posres) -> NPT (posres) -> production .tpr           (runs/<ff>/repN)
# Each replica gets its own NVT velocity seed (SEED_BASE + replica) and its
# own equilibration, so the N_REPLICAS production runs are independent.
# Production itself is run by 02_run_production.sh.
#
# Usage:
#   bash scripts/01_prepare.sh                    # all force fields
#   bash scripts/01_prepare.sh gromos54a7 amber   # a subset
#   FORCE=1 bash scripts/01_prepare.sh            # redo steps already completed
# =============================================================================

source "$(dirname "$0")/lib.sh"
check_gmx
parse_ff_args "$@"
FORCE="${FORCE:-0}"

# Download and clean the structure if it is missing (PDB_ID in config.sh)
if [[ ! -s "$INPUT_PDB" && -n "${PDB_ID:-}" ]]; then
    bash "$REPO_ROOT/scripts/get_structure.sh"
fi
need_file "$INPUT_PDB"
# Regenerate the .mdp files from config.sh every time (config.sh is the
# single source of truth; edit it, not the .mdp files).
bash "$REPO_ROOT/scripts/00_make_mdps.sh"

# mdrun thread flags: thread-MPI builds accept -ntmpi, real MPI builds do not.
if "$GMX" -quiet --version 2>/dev/null | grep -Eq '^MPI library: *MPI'; then
    THREAD_FLAGS="-ntomp $NT"
else
    THREAD_FLAGS="-ntmpi 1 -ntomp $NT"
fi

# done_or_force <file>: true when the step must run
must_run() { [[ "$FORCE" == 1 || ! -s "$1" ]]; }

# -----------------------------------------------------------------------------
# 1. Input structure: downloaded and cleaned by get_structure.sh when missing;
#    cleaned again here (idempotent) so any user-provided PDB is handled.
#    clean_pdb and check_gaps live in lib.sh.
# -----------------------------------------------------------------------------

# pdb2gmx histidine menu order is the same for every force field in
# GROMACS: 0 = H on ND1, 1 = H on NE2, 2 = both.
his_input() {
    local nhis=$1 s out=""
    read -r -a states <<<"$HIS_STATES"
    (( ${#states[@]} == nhis )) || die "HIS_STATES has ${#states[@]} entries but the structure has $nhis histidines"
    for s in "${states[@]}"; do
        case "${s^^}" in
            D) out+="0\n" ;; E) out+="1\n" ;; P) out+="2\n" ;;
            *) die "invalid HIS_STATES entry '$s' (use D, E or P)" ;;
        esac
    done
    printf '%s' "$out"
}

# Histidine states actually assigned by pdb2gmx, as D/E/P letters. pdb2gmx
# keeps the input name (HIS) in the coordinates and in [ atoms ], but writes
# the building block it used in a comment: "; residue   4 HIS rtp HISB ...".
# Multi-chain systems are written to topol_Protein_chain_*.itp.
his_states_from_top() {
    local files=(topol.top) f
    for f in topol_Protein*.itp; do [[ -e $f ]] && files+=("$f"); done
    awk '$1 == ";" && $2 == "residue" && $5 == "rtp" { print $6 }' "${files[@]}" |
    awk 'BEGIN { m["HSD"]="D"; m["HISA"]="D"; m["HISD"]="D"; m["HID"]="D"
                 m["HSE"]="E"; m["HISB"]="E"; m["HISE"]="E"; m["HIE"]="E"
                 m["HSP"]="P"; m["HISH"]="P"; m["HIP"]="P" }
         ($1 in m) { printf "%s ", m[$1] }'
}

# -----------------------------------------------------------------------------
# 2. System build + energy minimization
# -----------------------------------------------------------------------------
build_system() {
    local ff=$1 ffname water ffpath pos neg solvent_box
    ffname=$(ff_dir_name "$ff"); water=$(ff_water "$ff")
    ffpath=$(require_ff "$ff")
    check_ff_combo "$ff"
    solvent_box=$(ff_solvent_box "$ff")

    local B="$RUN_DIR/$ff/build" L
    mkdir -p "$B"; L="$REPO_ROOT/$B/prepare.log"
    log "[$ff] $(ff_display_name "$ff") | force field: $ffpath | water: $water"

    if must_run "$B/em.gro"; then
        : >"$L"
        clean_pdb "$INPUT_PDB" "$B/protein_clean.pdb"
        check_gaps "$B/protein_clean.pdb"

        local nhis hisflag="" hisin=""
        nhis=$(awk 'substr($0,13,4)==" CA " && substr($0,18,3)=="HIS"' "$B/protein_clean.pdb" | wc -l)
        if [[ "$HIS_STATES" != auto && $nhis -gt 0 ]]; then
            hisflag="-his"; hisin=$(his_input "$nhis")
        fi

        pushd "$B" >/dev/null
        # pdb2gmx: -ignh discards input hydrogens so each force field builds
        # its own (all-atom for CHARMM and AMBER, polar + aromatic H only for
        # GROMOS). Residue renaming for CMAP assignment (ff19SB, GROMACS >=
        # 2026, -rtpres auto) is left at its default.
        log "[$ff] pdb2gmx"
        gmx_run_input "$L" "$hisin" pdb2gmx -f protein_clean.pdb -o processed.gro \
            -p topol.top -i posre.itp -ff "$ffname" -water "$water" -ignh $hisflag
        { grep -E "^Total charge" "$L" || true; } | tail -n 1 | sed "s/^/[$ff] /"
        his_states_from_top >his_states.txt

        # Ion molecule names come from the ion .itp the topology includes
        # (water-model-specific in the built-in AMBER ff14SB/ff19SB ports).
        read -r pos neg <<<"$(ff_ion_names "$ffpath" topol.top)"
        log "[$ff] ions: $pos/$neg | solvent box: $solvent_box"

        log "[$ff] box: $BOX_TYPE, ${BOX_DIST} nm"
        gmx_run "$L" editconf -f processed.gro -o boxed.gro -c -d "$BOX_DIST" -bt "$BOX_TYPE"

        log "[$ff] solvate"
        gmx_run "$L" solvate -cp boxed.gro -cs "$solvent_box" -o solvated.gro -p topol.top

        log "[$ff] ions: neutralize + ${SALT_CONC} M NaCl"
        gmx_run "$L" grompp -f "$REPO_ROOT/$MDP_DIR/$ff/ions.mdp" -c solvated.gro \
            -p topol.top -o ions.tpr -po mdout_ions.mdp -maxwarn "$(ff_maxwarn "$ff")"
        gmx_run_input "$L" "SOL\n" genion -s ions.tpr -o ionized.gro -p topol.top \
            -pname "$pos" -nname "$neg" -neutral -conc "$SALT_CONC"

        log "[$ff] energy minimization"
        gmx_run "$L" grompp -f "$REPO_ROOT/$MDP_DIR/$ff/em.mdp" -c ionized.gro \
            -p topol.top -o em.tpr -po mdout_em.mdp -maxwarn "$(ff_maxwarn "$ff")"
        gmx_run "$L" mdrun -deffnm em $THREAD_FLAGS
        popd >/dev/null
    else
        log "[$ff] build/em.gro exists, skipping build (FORCE=1 to redo)"
    fi

    # Report EM convergence (steepest descent stops at emtol or nsteps)
    if grep -q "converged to Fmax <" "$B/em.log"; then
        log "[$ff] EM converged: $(grep -E '^Potential Energy|^Maximum force' "$B/em.log" | tr -s ' ' | paste -sd ';')"
    else
        warn "[$ff] EM did not reach emtol; inspect $B/em.log before trusting the run"
    fi
    log "[$ff] atoms: $(sed -n 2p "$B/em.gro" | tr -d ' ')"
}

# -----------------------------------------------------------------------------
# 3. Equilibration per replica and production .tpr
# -----------------------------------------------------------------------------
equilibrate_replica() {
    local ff=$1 r=$2
    local B="$RUN_DIR/$ff/build" R="$RUN_DIR/$ff/rep$r" seed=$(( SEED_BASE + r ))
    mkdir -p "$R"
    local L="$REPO_ROOT/$R/equil.log"
    local Rabs="$REPO_ROOT/$R"

    # Per-replica copies of the mdp files (the NVT seed is replica-specific)
    sed "s/@SEED@/$seed/" "$MDP_DIR/$ff/nvt.mdp" >"$R/nvt.mdp"
    cp "$MDP_DIR/$ff/npt.mdp" "$MDP_DIR/$ff/md.mdp" "$R/"

    # grompp runs inside build/ so that topol.top finds posre.itp and the
    # force field includes; outputs go to the replica directory.
    if must_run "$R/nvt.gro"; then
        log "[$ff rep$r] NVT ${NVT_PS} ps (seed $seed)"
        ( cd "$B" && gmx_run "$L" grompp -f "$Rabs/nvt.mdp" -c em.gro -r em.gro \
            -p topol.top -o "$Rabs/nvt.tpr" -po "$Rabs/mdout_nvt.mdp" -maxwarn "$(ff_maxwarn "$ff")" )
        ( cd "$R" && gmx_run "$L" mdrun -deffnm nvt $THREAD_FLAGS $MDRUN_FLAGS )
    else
        log "[$ff rep$r] nvt.gro exists, skipping"
    fi

    if must_run "$R/npt.gro"; then
        log "[$ff rep$r] NPT ${NPT_PS} ps"
        ( cd "$B" && gmx_run "$L" grompp -f "$Rabs/npt.mdp" -c "$Rabs/nvt.gro" \
            -r "$Rabs/nvt.gro" -t "$Rabs/nvt.cpt" -p topol.top -o "$Rabs/npt.tpr" \
            -po "$Rabs/mdout_npt.mdp" -maxwarn "$(ff_maxwarn "$ff")" )
        ( cd "$R" && gmx_run "$L" mdrun -deffnm npt $THREAD_FLAGS $MDRUN_FLAGS )
    else
        log "[$ff rep$r] npt.gro exists, skipping"
    fi

    if must_run "$R/md.tpr"; then
        log "[$ff rep$r] production tpr (${PROD_NS} ns)"
        ( cd "$B" && gmx_run "$L" grompp -f "$Rabs/md.mdp" -c "$Rabs/npt.gro" \
            -t "$Rabs/npt.cpt" -p topol.top -o "$Rabs/md.tpr" \
            -po "$Rabs/mdout_md.mdp" -maxwarn "$(ff_maxwarn "$ff")" )
    fi
}

# Cross-check that every built force field received the same protonation
# states and the same net charge; differences would confound the comparison.
# Every system is compared with the first one that has been built.
crosscheck() {
    local ff b h q ref_h="" ref_q="" ref_ff=""
    for ff in "${FF_LABELS[@]}"; do
        b="$RUN_DIR/$ff/build"
        [[ -s "$b/his_states.txt" && -s "$b/prepare.log" ]] || continue
        h=$(cat "$b/his_states.txt")
        q=$(grep -E "^Total charge" "$b/prepare.log" | tail -n 1 | grep -oE -- '-?[0-9]+\.[0-9]+' || true)
        log "$(printf '%-11s' "$ff") HIS: [${h% }]  net charge: ${q:-?}"
        if [[ -z $ref_ff ]]; then
            ref_ff=$ff; ref_h=$h; ref_q=$q; continue
        fi
        [[ "$h" == "$ref_h" ]] || warn "histidine states differ between $ff and $ref_ff; set HIS_STATES explicitly and rerun with FORCE=1"
        [[ "$q" == "$ref_q" ]] || warn "net protein charge differs between $ff and $ref_ff; check termini and protonation"
    done
}

# -----------------------------------------------------------------------------
log "GROMACS $GMX_VERSION | threads: $NT | replicas: $N_REPLICAS | force fields: ${SELECTED_FF[*]}"
for ff in "${SELECTED_FF[@]}"; do
    build_system "$ff"
done
crosscheck
for ff in "${SELECTED_FF[@]}"; do
    for (( r = 1; r <= N_REPLICAS; r++ )); do
        equilibrate_replica "$ff" "$r"
    done
done
log "done. Next: bash scripts/02_run_production.sh"
