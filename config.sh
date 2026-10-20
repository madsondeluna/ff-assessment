# shellcheck shell=bash
# =============================================================================
# config.sh
# Central configuration for the four-way force field comparison of a protein
# (ubiquitin, PDB 1UBQ, by default):
#   CHARMM36m | AMBER (ff14SB by default) | GROMOS 53A6 | GROMOS 54A7
# Every script in scripts/ sources this file. Any variable can be overridden
# from the environment, e.g.:  PROD_NS=10 NT=8 bash scripts/01_prepare.sh
# =============================================================================

# -----------------------------------------------------------------------------
# GROMACS executable and parallelism
# -----------------------------------------------------------------------------
# GROMACS executable. Empty = autodetect: the first of gmx, gmx_mpi, gmx_d,
# gmx_mpi_d found in PATH (so a machine with a single GROMACS installation,
# from apt, conda or a module, needs no configuration). Set it explicitly
# (name or full path) only when several versions are installed.
# Requires GROMACS >= 2023 (>= 2026 for the built-in AMBER ff14SB / ff19SB).
GMX="${GMX:-}"

# OpenMP threads for mdrun. Default: all cores minus one (keeps one core free).
_ncores=$(nproc 2>/dev/null || echo 2)
NT="${NT:-$(( _ncores > 1 ? _ncores - 1 : 1 ))}"

# Extra mdrun flags, e.g. GPU offload:
#   MDRUN_FLAGS="-nb gpu -pme gpu -bonded gpu"
# Note: '-update gpu' may be refused for GROMOS because all bonds are
# constrained (large coupled LINCS clusters); leave it on auto.
MDRUN_FLAGS="${MDRUN_FLAGS:-}"

# Number of grompp warnings tolerated (keep 0 unless you know why).
MAXWARN="${MAXWARN:-0}"

# -----------------------------------------------------------------------------
# Input structure
# -----------------------------------------------------------------------------
# PDB entry to download and clean automatically when INPUT_PDB is missing
# (scripts/get_structure.sh). 1UBQ is human ubiquitin at 1.8 A: 76 residues,
# no gaps, no disulfides, one alpha helix, one 3-10 helix and a five-stranded
# mixed beta sheet, which makes it a good small alpha/beta test case for
# secondary structure (lysozyme, 1AKI, is almost all helix). Alternatives:
# 1PGA (protein G B1, 56 residues, one helix + four strands). Set PDB_ID=""
# to disable the download and provide INPUT_PDB yourself.
PDB_ID="${PDB_ID:-1UBQ}"

# Input structure (PDB format): the cleaned download above, or any
# experimental or modelled structure. Ligands, waters, ions and non-protein
# residues are removed automatically; hydrogens are rebuilt by pdb2gmx.
INPUT_PDB="${INPUT_PDB:-input/1ubq.pdb}"

# Chain(s) to keep, e.g. "A" or "AB". Empty keeps all protein chains.
CHAIN="${CHAIN-A}"

# Abort if chain breaks (missing residues) are detected. Set to 1 only if
# the gaps are intentional (e.g. you will model them separately).
ALLOW_GAPS="${ALLOW_GAPS:-0}"

# Histidine protonation. "auto" lets pdb2gmx decide from the H-bond network
# (same algorithm for every force field; the result is cross-checked).
# Otherwise give one letter per HIS, in sequence order:
#   D = H on ND1 (HSD / HID / HISA), E = H on NE2 (HSE / HIE / HISB),
#   P = both (HSP / HIP / HISH)
# Example for three histidines: HIS_STATES="E D E"
HIS_STATES="${HIS_STATES:-auto}"

# -----------------------------------------------------------------------------
# Force fields
# -----------------------------------------------------------------------------
# Labels used for directory names (mdp/, runs/, analysis/) and plots. Each
# label maps to a force field directory and a water model below; the
# force-field family (charmm, amber, gromos) decides the nonbonded and
# constraint settings written by scripts/00_make_mdps.sh.
#
# Remove a label from this list to skip that force field, override it from
# the environment (FF_LABELS="amber gromos54a7" bash scripts/run_all.sh), or
# pass labels on the command line:  bash scripts/01_prepare.sh charmm36m amber
# (bash cannot export arrays: the list is kept as the exported string
# FF_LABELS_STR for child scripts and converted to an array here)
FF_LABELS_STR="${FF_LABELS:-${FF_LABELS_STR:-charmm36m amber gromos53a6 gromos54a7}}"
export FF_LABELS_STR
read -r -a FF_LABELS <<<"$FF_LABELS_STR"

# Reference force field for the difference panels in 04_plot_comparison.py
# (every other force field is plotted as "ff minus reference").
FF_REFERENCE="${FF_REFERENCE:-charmm36m}"

# --- CHARMM36m -------------------------------------------------------------
# GROMACS port from the MacKerell lab. Place the extracted directory in
# forcefields/ (e.g. forcefields/charmm36-feb2026_cgenff-5.0.ff) and give the
# name WITHOUT the .ff suffix. Current ports contain the CHARMM36m protein
# parameters (the 'm' revision of the backbone CMAP) by default.
CHARMM_FF="${CHARMM_FF:-charmm36-feb2026_cgenff-5.0}"
CHARMM_WATER="${CHARMM_WATER:-tip3p}"      # CHARMM-modified TIP3P (LJ on H)

# --- AMBER -------------------------------------------------------------------
# Which AMBER protein force field the label "amber" refers to:
#   auto           amber14sb when the installed GROMACS (or forcefields/)
#                  has it, otherwise amber99sb-ildn, with a warning. Default.
#   amber14sb      ff14SB, TIP3P. Built into GROMACS >= 2026 (validation
#                  pending per the release notes): the usual reference AMBER
#                  force field and a 3-site water like the other systems.
#   amber19sb      ff19SB with OPC (4-site water, per-residue CMAP). Built
#                  into GROMACS >= 2026; requires AMBER_WATER=opc (or opc3).
#   amber99sb-ildn ships with every GROMACS version (AMBER_WATER=tip3p).
# A user-contributed port can be extracted into forcefields/ and named here
# the same way as the CHARMM port (directory name without .ff).
AMBER_FF="${AMBER_FF:-auto}"
AMBER_WATER="${AMBER_WATER:-tip3p}"        # standard TIP3P (no LJ on H)

# --- GROMOS 53A6 and 54A7 (ship with GROMACS) ----------------------------------
GROMOS53_FF="${GROMOS53_FF:-gromos53a6}"
GROMOS54_FF="${GROMOS54_FF:-gromos54a7}"
GROMOS_WATER="${GROMOS_WATER:-spc}"        # both were parametrized with SPC

# -----------------------------------------------------------------------------
# System setup (identical for every force field)
# -----------------------------------------------------------------------------
BOX_TYPE="${BOX_TYPE:-dodecahedron}"
BOX_DIST="${BOX_DIST:-1.2}"         # nm, solute-box distance
SALT_CONC="${SALT_CONC:-0.15}"      # mol/L NaCl on top of neutralization
TEMP="${TEMP:-300}"                 # K
PRESSURE="${PRESSURE:-1.0}"         # bar

# -----------------------------------------------------------------------------
# Protocol (identical for every force field; only nonbonded settings and
# constraints follow each family's parametrization, see 00_make_mdps.sh)
# -----------------------------------------------------------------------------
DT="${DT:-0.002}"                   # ps
EM_NSTEPS="${EM_NSTEPS:-50000}"
NVT_PS="${NVT_PS:-100}"             # ps, position-restrained
NPT_PS="${NPT_PS:-100}"             # ps, position-restrained
PROD_NS="${PROD_NS:-100}"           # ns, production
XTC_PS="${XTC_PS:-10}"              # ps between saved frames (production)
EDR_PS="${EDR_PS:-10}"              # ps between energy frames (production)

# Independent replicas per force field (different velocity seeds).
# One replica answers "what happened in this run"; three are the minimum to
# separate force field effects from run-to-run variability, so the comparison
# runs 3 x 4 = 12 production simulations by default.
N_REPLICAS="${N_REPLICAS:-3}"
SEED_BASE="${SEED_BASE:-2026}"      # replica r uses SEED_BASE + r

# -----------------------------------------------------------------------------
# Analysis
# -----------------------------------------------------------------------------
# Initial part of the production run excluded from averages, RMSF, SASA per
# residue and contact maps (time series are always plotted in full).
SKIP_NS="${SKIP_NS:-10}"

# Stride for the expensive analyses (ps). 0 = every saved frame.
HB_WATER_DT_PS="${HB_WATER_DT_PS:-0}"
MDMAT_DT_PS="${MDMAT_DT_PS:-100}"
DSSP_DT_PS="${DSSP_DT_PS:-100}"     # secondary structure (gmx dssp) per frame

# Contact-map cutoff used by gmx mdmat (nm).
MDMAT_CUTOFF="${MDMAT_CUTOFF:-1.5}"

# -----------------------------------------------------------------------------
# Execution
# -----------------------------------------------------------------------------
# Production runs are strictly sequential (one mdrun at a time). RUN_ORDER
# decides the sequence: "replica" runs replica 1 of every force field first,
# then replica 2, ... so that a partial campaign already allows a one-replica
# comparison; "forcefield" finishes all replicas of one force field before
# starting the next.
RUN_ORDER="${RUN_ORDER:-replica}"

# -----------------------------------------------------------------------------
# Paths (relative to the repository root)
# -----------------------------------------------------------------------------
MDP_DIR="${MDP_DIR:-mdp}"
RUN_DIR="${RUN_DIR:-runs}"
ANA_DIR="${ANA_DIR:-analysis}"
FIG_DIR="${FIG_DIR:-figures}"
FF_DIR="${FF_DIR:-forcefields}"
LOG_DIR="${LOG_DIR:-logs}"          # pipeline and production logs
