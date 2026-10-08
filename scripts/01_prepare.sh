#!/usr/bin/env bash
# =============================================================================
# 01_prepare.sh
# Builds and equilibrates eIF4E in water with CHARMM36(m) and GROMOS 53A6:
#   clean PDB -> pdb2gmx -> box -> solvate -> ions -> EM      (runs/<ff>/build)
#   NVT (posres) -> NPT (posres) -> production .tpr           (runs/<ff>/repN)
# Production itself is run by 02_run_production.sh.
#
# Usage:
#   bash scripts/01_prepare.sh              # both force fields
#   bash scripts/01_prepare.sh gromos53a6   # only one
#   FORCE=1 bash scripts/01_prepare.sh      # redo steps already completed
# =============================================================================

source "$(dirname "$0")/lib.sh"
check_gmx
parse_ff_args "$@"
FORCE="${FORCE:-0}"

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
# 1. Clean the input structure (identical file for both force fields)
# -----------------------------------------------------------------------------
# Keeps only standard amino acids of the selected chain(s), first alternate
# location, renames protonation-specific residue names back to the generic
# ones (HIS, CYS) so each force field assigns its own, converts MSE to MET,
# and drops hydrogens (pdb2gmx -ignh rebuilds them for each force field).
clean_pdb() {
    local in=$1 out=$2
    awk -v chains="$CHAIN" '
    BEGIN {
        n = split("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL", aa, " ")
        for (i = 1; i <= n; i++) std[aa[i]] = 1
        split("HID HIE HIP HSD HSE HSP HISA HISB HISH HISD HISE", h, " "); for (i in h) his[h[i]] = 1
        split("CYX CYM CYS2", c, " "); for (i in c) cys[c[i]] = 1
    }
    /^(ATOM  |HETATM)/ {
        res = substr($0, 18, 4); gsub(/ /, "", res)
        if (res in his) res = "HIS"
        if (res in cys) res = "CYS"
        name = substr($0, 13, 4); gsub(/ /, "", name)
        if (res == "MSE") { res = "MET"; if (name == "SE") name = "SD" }
        if (!(res in std)) next
        alt = substr($0, 17, 1); if (alt != " " && alt != "A") next
        ch = substr($0, 22, 1); if (chains != "" && index(chains, ch) == 0) next
        elem = substr($0, 77, 2); gsub(/ /, "", elem)
        if (elem == "H" || elem == "D" || name ~ /^[0-9]*H/ || name ~ /^D/) next
        line = $0
        # rewrite record type, residue name (cols 18-20) and atom name if changed
        line = "ATOM  " substr(line, 7)
        line = substr(line, 1, 16) " " sprintf("%-3s", res) " " substr(line, 22)
        if (name == "SD" && substr($0, 13, 4) ~ /SE/) line = substr(line, 1, 12) " SD " substr(line, 17)
        print line
        next
    }
    /^TER/ { print "TER"; next }
    END { print "END" }
    ' "$in" >"$out"
    [[ $(grep -c '^ATOM' "$out") -gt 0 ]] || die "no protein atoms left after cleaning $in (check CHAIN='$CHAIN')"
}

# Detect chain breaks: consecutive CA atoms farther than 4.2 A apart.
check_gaps() {
    local pdb=$1 gaps
    gaps=$(awk '
        substr($0,13,4) == " CA " {
            ch = substr($0,22,1); r = substr($0,23,4) + 0
            x = substr($0,31,8); y = substr($0,39,8); z = substr($0,47,8)
            if (ch == pch) {
                d = sqrt((x-px)^2 + (y-py)^2 + (z-pz)^2)
                if (d > 4.2) printf "  chain %s: residue %d -> %d (CA-CA %.1f A)\n", ch, pr, r, d
            }
            pch = ch; pr = r; px = x; py = y; pz = z
        }' "$pdb")
    if [[ -n "$gaps" ]]; then
        if [[ "$ALLOW_GAPS" == 1 ]]; then
            warn "chain breaks found (continuing because ALLOW_GAPS=1):"$'\n'"$gaps"
        else
            die "chain breaks (missing residues) found:
$gaps
  pdb2gmx would bond across these gaps. Model the missing residues
  (MODELLER, PDBFixer, AlphaFold model) or trim the termini, then rerun.
  Set ALLOW_GAPS=1 only if you know what you are doing."
        fi
    fi
}

# pdb2gmx histidine menu order is identical for CHARMM36 and GROMOS in
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
    local ff=$1 ffname water ffpath pos neg
    ffname=$(ff_dir_name "$ff"); water=$(ff_water "$ff")
    ffpath=$(require_ff "$ff")
    read -r pos neg <<<"$(ff_ion_names "$ffpath")"

    local B="$RUN_DIR/$ff/build" L
    mkdir -p "$B"; L="$REPO_ROOT/$B/prepare.log"
    log "[$ff] force field: $ffpath | water: $water | ions: $pos/$neg"

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
        # its own (all-atom for CHARMM, polar + aromatic H only for GROMOS).
        log "[$ff] pdb2gmx"
        gmx_run_input "$L" "$hisin" pdb2gmx -f protein_clean.pdb -o processed.gro \
            -p topol.top -i posre.itp -ff "$ffname" -water "$water" -ignh $hisflag
        { grep -E "^Total charge" "$L" || true; } | tail -n 1 | sed "s/^/[$ff] /"
        his_states_from_top >his_states.txt

        log "[$ff] box: $BOX_TYPE, ${BOX_DIST} nm"
        gmx_run "$L" editconf -f processed.gro -o boxed.gro -c -d "$BOX_DIST" -bt "$BOX_TYPE"

        log "[$ff] solvate"
        gmx_run "$L" solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

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

# Cross-check that both force fields received the same protonation states
# and the same net charge; differences would confound the comparison.
crosscheck() {
    local a="$RUN_DIR/charmm36/build" b="$RUN_DIR/gromos53a6/build"
    [[ -s "$a/his_states.txt" && -s "$b/his_states.txt" ]] || return 0
    local ha hb qa qb
    ha=$(cat "$a/his_states.txt"); hb=$(cat "$b/his_states.txt")
    qa=$(grep -E "^Total charge" "$a/prepare.log" | tail -n 1 | grep -oE -- '-?[0-9]+\.[0-9]+' || true)
    qb=$(grep -E "^Total charge" "$b/prepare.log" | tail -n 1 | grep -oE -- '-?[0-9]+\.[0-9]+' || true)
    log "HIS states  charmm36: [${ha% }]  gromos53a6: [${hb% }]"
    log "net charge  charmm36: ${qa:-?}  gromos53a6: ${qb:-?}"
    [[ "$ha" == "$hb" ]] || warn "histidine states differ between force fields; set HIS_STATES explicitly and rerun with FORCE=1"
    [[ "$qa" == "$qb" ]] || warn "net protein charge differs between force fields; check termini and protonation"
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
