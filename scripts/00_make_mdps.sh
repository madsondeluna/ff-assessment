#!/usr/bin/env bash
# =============================================================================
# 00_make_mdps.sh
# Writes the .mdp files for both force fields into mdp/<label>/.
#
# Design of the comparison
#   The protocol follows the standard lysozyme-in-water workflow
#   (ions -> EM -> NVT 100 ps -> NPT 100 ps -> production), and every
#   protocol parameter is IDENTICAL for both force fields: integrator, time
#   step, thermostat, barostat, coupling groups, restraints, output, seeds.
#
#   Only two blocks follow each force field's own parametrization, because
#   using one force field's settings on the other is a known source of
#   artefacts and would confound the comparison:
#     1. Nonbonded interactions (cutoffs, vdW modifier, dispersion correction)
#     2. Bond constraints
#
#   CHARMM36/36m : PME, rc = 1.2 nm, LJ force-switch 1.0-1.2 nm, no DispCorr,
#                  H-bond constraints (as in the CHARMM-GUI/MacKerell inputs).
#   GROMOS 53A6  : PME, rc = 1.4 nm (the GROMOS parametrization cutoff),
#                  plain LJ cutoff, DispCorr EnerPres, all bonds constrained
#                  (GROMOS was parametrized with rigid bonds).
#
# Usage: bash scripts/00_make_mdps.sh
# =============================================================================

source "$(dirname "$0")/lib.sh"

# Steps derived from config.sh (ps -> steps with the chosen dt)
NVT_STEPS=$(fcalc "$NVT_PS / $DT")
NPT_STEPS=$(fcalc "$NPT_PS / $DT")
PROD_STEPS=$(fcalc "$PROD_NS * 1000 / $DT")
XTC_STEPS=$(fcalc "$XTC_PS / $DT")
EDR_STEPS=$(fcalc "$EDR_PS / $DT")
EQ_OUT_STEPS=$(fcalc "1.0 / $DT")          # 1 ps output during equilibration

# -----------------------------------------------------------------------------
# Force-field-specific blocks
# -----------------------------------------------------------------------------
nonbonded_block() {
    case "$1" in
        charmm36) cat <<'EOF'
; --- Nonbonded: CHARMM36 / CHARMM36m recommended settings ---------------------
cutoff-scheme            = Verlet
nstlist                  = 10        ; adjusted automatically by mdrun
coulombtype              = PME
rcoulomb                 = 1.2       ; nm
vdwtype                  = Cut-off
vdw-modifier             = Force-switch
rvdw-switch              = 1.0       ; nm, LJ forces smoothly switched to zero
rvdw                     = 1.2       ; nm
DispCorr                 = no        ; CHARMM36 was parametrized without it
pme-order                = 4
fourierspacing           = 0.16      ; nm
EOF
        ;;
        gromos53a6) cat <<'EOF'
; --- Nonbonded: GROMOS 53A6 recommended settings (Verlet scheme) -------------
; GROMOS was parametrized with a 1.4 nm cutoff (twin-range + reaction field).
; Here PME replaces the reaction field so that long-range electrostatics are
; treated identically in both force fields; the LJ cutoff stays at 1.4 nm.
cutoff-scheme            = Verlet
nstlist                  = 10        ; adjusted automatically by mdrun
coulombtype              = PME
rcoulomb                 = 1.4       ; nm
vdwtype                  = Cut-off
vdw-modifier             = Potential-shift
rvdw                     = 1.4       ; nm
DispCorr                 = EnerPres  ; long-range dispersion correction
pme-order                = 4
fourierspacing           = 0.16      ; nm
EOF
        ;;
    esac
}

constraints_block() {
    case "$1" in
        charmm36) cat <<'EOF'
; --- Constraints: CHARMM36 (bonds to hydrogen) ------------------------------
constraints              = h-bonds
constraint-algorithm     = lincs
lincs-iter               = 1
lincs-order              = 4
EOF
        ;;
        gromos53a6) cat <<'EOF'
; --- Constraints: GROMOS 53A6 (all bonds, as in its parametrization) ---------
constraints              = all-bonds
constraint-algorithm     = lincs
lincs-iter               = 1
lincs-order              = 4
EOF
        ;;
    esac
}

# -----------------------------------------------------------------------------
# Protocol blocks (identical for both force fields)
# -----------------------------------------------------------------------------
write_mdps() {
    local label=$1 out="$MDP_DIR/$1"
    mkdir -p "$out"

    # ---- ions.mdp: only used to build a .tpr for gmx genion ----------------
    # Plain cutoff electrostatics here: the system is not yet neutral, and
    # grompp (>= 2023) refuses PME with a net charge unless -maxwarn is used.
    # This .tpr is never simulated, so the electrostatics choice is irrelevant.
    {
        cat <<EOF
; ions.mdp ($label) - input for gmx genion only, never run
integrator               = steep
emtol                    = 1000.0
emstep                   = 0.01
nsteps                   = $EM_NSTEPS
pbc                      = xyz
EOF
        nonbonded_block "$label" |
            sed -E 's/^coulombtype( *)= PME.*/coulombtype\1= Cut-off   ; genion only/' |
            grep -vE '^(pme-order|fourierspacing)'
    } >"$out/ions.mdp"

    # ---- em.mdp ------------------------------------------------------------
    {
        cat <<EOF
; em.mdp ($label) - steepest descent energy minimization
integrator               = steep
emtol                    = 1000.0    ; kJ/mol/nm, stop when Fmax < emtol
emstep                   = 0.01      ; nm
nsteps                   = $EM_NSTEPS
pbc                      = xyz
EOF
        nonbonded_block "$label"
    } >"$out/em.mdp"

    # ---- nvt.mdp -----------------------------------------------------------
    {
        cat <<EOF
; nvt.mdp ($label) - ${NVT_PS} ps NVT, heavy-atom position restraints
define                   = -DPOSRES
integrator               = md        ; leap-frog
dt                       = $DT
nsteps                   = $NVT_STEPS
; Output
nstxout                  = $EQ_OUT_STEPS
nstvout                  = $EQ_OUT_STEPS
nstenergy                = $EQ_OUT_STEPS
nstlog                   = $EQ_OUT_STEPS
; Temperature coupling
tcoupl                   = V-rescale
tc-grps                  = Protein Non-Protein
tau-t                    = 0.1 0.1   ; ps
ref-t                    = $TEMP $TEMP
; Pressure coupling
pcoupl                   = no
; Velocities
continuation             = no
gen-vel                  = yes
gen-temp                 = $TEMP
gen-seed                 = @SEED@    ; replaced per replica by 01_prepare.sh
pbc                      = xyz
EOF
        constraints_block "$label"
        nonbonded_block "$label"
    } >"$out/nvt.mdp"

    # ---- npt.mdp -----------------------------------------------------------
    # C-rescale (GROMACS >= 2021) is used for equilibration: it is correct for
    # fluctuations and, unlike Parrinello-Rahman, robust far from equilibrium.
    {
        cat <<EOF
; npt.mdp ($label) - ${NPT_PS} ps NPT, heavy-atom position restraints
define                   = -DPOSRES
integrator               = md
dt                       = $DT
nsteps                   = $NPT_STEPS
; Output
nstxout                  = $EQ_OUT_STEPS
nstvout                  = $EQ_OUT_STEPS
nstenergy                = $EQ_OUT_STEPS
nstlog                   = $EQ_OUT_STEPS
; Temperature coupling
tcoupl                   = V-rescale
tc-grps                  = Protein Non-Protein
tau-t                    = 0.1 0.1
ref-t                    = $TEMP $TEMP
; Pressure coupling
pcoupl                   = C-rescale
pcoupltype               = isotropic
tau-p                    = 2.0       ; ps
ref-p                    = $PRESSURE
compressibility          = 4.5e-5    ; bar^-1, water
refcoord-scaling         = com
; Continuation from NVT
continuation             = yes
gen-vel                  = no
pbc                      = xyz
EOF
        constraints_block "$label"
        nonbonded_block "$label"
    } >"$out/npt.mdp"

    # ---- md.mdp ------------------------------------------------------------
    {
        cat <<EOF
; md.mdp ($label) - ${PROD_NS} ns production, unrestrained
integrator               = md
dt                       = $DT
nsteps                   = $PROD_STEPS
; Output: compressed coordinates of the whole system (water is needed for
; the protein-water hydrogen bond analysis)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstxout-compressed       = $XTC_STEPS
compressed-x-grps        = System
nstenergy                = $EDR_STEPS
nstlog                   = $EDR_STEPS
nstcalcenergy            = 100
; Temperature coupling
tcoupl                   = V-rescale
tc-grps                  = Protein Non-Protein
tau-t                    = 0.1 0.1
ref-t                    = $TEMP $TEMP
; Pressure coupling
pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau-p                    = 2.0
ref-p                    = $PRESSURE
compressibility          = 4.5e-5
; Continuation from NPT
continuation             = yes
gen-vel                  = no
pbc                      = xyz
EOF
        constraints_block "$label"
        nonbonded_block "$label"
    } >"$out/md.mdp"

    log "wrote $out/{ions,em,nvt,npt,md}.mdp"
}

for ff in "${FF_LABELS[@]}"; do
    write_mdps "$ff"
done

log "production: ${PROD_NS} ns = ${PROD_STEPS} steps (dt = ${DT} ps), frames every ${XTC_PS} ps"
