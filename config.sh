# shellcheck shell=bash
# =============================================================================
# config.sh
# Central configuration for the CHARMM36 vs GROMOS 53A6 comparison.
# Every script in scripts/ sources this file. Any variable can be overridden
# from the environment, e.g.:  PROD_NS=10 NT=8 bash scripts/01_prepare.sh
# =============================================================================

# -----------------------------------------------------------------------------
# GROMACS executable and parallelism
# -----------------------------------------------------------------------------
# gmx, gmx_mpi, or a full path. Requires GROMACS >= 2023.
GMX="${GMX:-gmx}"

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
# Experimental or modelled eIF4E structure (PDB format). Ligands, waters,
# ions and non-protein residues are removed automatically.
INPUT_PDB="${INPUT_PDB:-input/eif4e.pdb}"

# Chain(s) to keep, e.g. "A" or "AB". Empty keeps all protein chains.
CHAIN="${CHAIN-A}"

# Abort if chain breaks (missing residues) are detected. Set to 1 only if
# the gaps are intentional (e.g. you will model them separately).
ALLOW_GAPS="${ALLOW_GAPS:-0}"

# Histidine protonation. "auto" lets pdb2gmx decide from the H-bond network
# (same algorithm for both force fields; the result is cross-checked).
# Otherwise give one letter per HIS, in sequence order:
#   D = H on ND1 (HSD / HISA), E = H on NE2 (HSE / HISB), P = both (HSP / HISH)
# Example for three histidines: HIS_STATES="E D E"
HIS_STATES="${HIS_STATES:-auto}"

# -----------------------------------------------------------------------------
# Force fields
# -----------------------------------------------------------------------------
# Labels used for directory names and plots. Do not change the labels; change
# the force field directories below instead.
FF_LABELS=(charmm36 gromos53a6)

# CHARMM36 GROMACS port (MacKerell lab). Place the extracted directory in
# forcefields/ (e.g. forcefields/charmm36-feb2026_cgenff-5.0.ff) and give the
# name WITHOUT the .ff suffix. The port contains the CHARMM36m protein
# parameters by default.
CHARMM_FF="${CHARMM_FF:-charmm36-feb2026_cgenff-5.0}"
CHARMM_WATER="${CHARMM_WATER:-tip3p}"      # CHARMM-modified TIP3P (LJ on H)

# GROMOS 53A6 ships with GROMACS.
GROMOS_FF="${GROMOS_FF:-gromos53a6}"
GROMOS_WATER="${GROMOS_WATER:-spc}"        # GROMOS was parametrized with SPC

# -----------------------------------------------------------------------------
# System setup (identical for both force fields)
# -----------------------------------------------------------------------------
BOX_TYPE="${BOX_TYPE:-dodecahedron}"
BOX_DIST="${BOX_DIST:-1.2}"         # nm, solute-box distance
SALT_CONC="${SALT_CONC:-0.15}"      # mol/L NaCl on top of neutralization
TEMP="${TEMP:-300}"                 # K
PRESSURE="${PRESSURE:-1.0}"         # bar

# -----------------------------------------------------------------------------
# Protocol (identical for both force fields; only nonbonded settings and
# constraints follow each force field's parametrization, see 00_make_mdps.sh)
# -----------------------------------------------------------------------------
DT="${DT:-0.002}"                   # ps
EM_NSTEPS="${EM_NSTEPS:-50000}"
NVT_PS="${NVT_PS:-100}"             # ps, position-restrained
NPT_PS="${NPT_PS:-100}"             # ps, position-restrained
PROD_NS="${PROD_NS:-100}"           # ns, production
XTC_PS="${XTC_PS:-10}"              # ps between saved frames (production)
EDR_PS="${EDR_PS:-10}"              # ps between energy frames (production)

# Independent replicas per force field (different velocity seeds).
# One replica answers "what happened in this run"; at least three are needed
# to separate force field effects from run-to-run variability.
N_REPLICAS="${N_REPLICAS:-1}"
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

# Contact-map cutoff used by gmx mdmat (nm).
MDMAT_CUTOFF="${MDMAT_CUTOFF:-1.5}"

# -----------------------------------------------------------------------------
# Paths (relative to the repository root)
# -----------------------------------------------------------------------------
MDP_DIR="${MDP_DIR:-mdp}"
RUN_DIR="${RUN_DIR:-runs}"
ANA_DIR="${ANA_DIR:-analysis}"
FIG_DIR="${FIG_DIR:-figures}"
FF_DIR="${FF_DIR:-forcefields}"
