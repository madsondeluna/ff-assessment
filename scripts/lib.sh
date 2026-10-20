# shellcheck shell=bash
# =============================================================================
# scripts/lib.sh
# Shared helpers sourced by every bash script in this repository.
# =============================================================================

set -euo pipefail

# Repository root = parent directory of scripts/
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# shellcheck source=../config.sh
source "$REPO_ROOT/config.sh"

# Make force fields placed in forcefields/ visible to pdb2gmx and grompp
# without touching the GROMACS installation. GMXLIB is searched in addition
# to the default share/top directory, so built-in force fields still work.
export GMXLIB="$REPO_ROOT/$FF_DIR${GMXLIB:+:$GMXLIB}"

# Prevent GROMACS from writing #backup# files on reruns; reruns overwrite.
export GMX_MAXBACKUP=-1

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------
log()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
warn() { printf '[%s] WARNING: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; exit 1; }

need_file() { [[ -s "$1" ]] || die "missing or empty file: $1"; }

# -----------------------------------------------------------------------------
# GROMACS wrappers
# -----------------------------------------------------------------------------

# Run gmx quietly, appending stdout/stderr to a log file. On failure, show the
# tail of the log so the user sees the actual GROMACS error.
#   gmx_run <logfile> <gmx subcommand and arguments...>
gmx_run() {
    local logfile=$1; shift
    if ! "$GMX" -quiet "$@" >>"$logfile" 2>&1; then
        echo "------ last lines of $logfile ------" >&2
        tail -n 30 "$logfile" >&2
        die "gmx $1 failed (full log: $logfile)"
    fi
}

# Same as gmx_run, but feeds interactive group selections from a string.
#   gmx_run_input <logfile> "<stdin text>" <gmx subcommand and arguments...>
gmx_run_input() {
    local logfile=$1 input=$2; shift 2
    if ! printf '%b' "$input" | "$GMX" -quiet "$@" >>"$logfile" 2>&1; then
        echo "------ last lines of $logfile ------" >&2
        tail -n 30 "$logfile" >&2
        die "gmx $1 failed (full log: $logfile)"
    fi
}

# Autodetect the GROMACS executable when GMX is empty.
detect_gmx() {
    [[ -n "$GMX" ]] && return 0
    local c
    for c in gmx gmx_mpi gmx_d gmx_mpi_d; do
        if command -v "$c" >/dev/null 2>&1; then GMX=$c; return 0; fi
    done
    return 1
}

check_gmx() {
    detect_gmx || die "GROMACS not found in PATH (looked for gmx, gmx_mpi, gmx_d, gmx_mpi_d).
  Ubuntu:        sudo apt install gromacs          (2023.3 on Ubuntu 24.04)
  conda:         conda install -c conda-forge gromacs
  module system: module load gromacs
  or set GMX to the executable, e.g. GMX=/opt/gromacs/bin/gmx"
    command -v "$GMX" >/dev/null 2>&1 || die "'$GMX' not found (GMX is set explicitly; check the name or path)."
    local ver major
    ver=$("$GMX" -quiet --version 2>/dev/null | awk -F': *' '/^GROMACS version/{print $2; exit}')
    major=$(sed -E 's/^([0-9]{4}).*/\1/' <<<"$ver")
    [[ "$major" =~ ^[0-9]{4}$ ]] || die "could not parse GROMACS version from '$ver'"
    (( major >= 2023 )) || die "GROMACS >= 2023 required, found $ver"
    GMX_VERSION="$ver"
    GMX_MAJOR="$major"
    resolve_amber_ff
}

# True if a gmx subcommand exists (used to pick legacy tools in 2024+).
gmx_has_tool() {
    local cmds
    cmds=$("$GMX" -quiet help commands 2>/dev/null || true)
    awk '{print $1}' <<<"$cmds" | grep -qx -- "$1"
}

# Directory holding the built-in force fields of this GROMACS installation.
gmx_share_top() {
    local prefix
    prefix=$("$GMX" -quiet --version 2>/dev/null | awk -F': *' '/^Data prefix/{print $2; exit}')
    echo "$prefix/share/gromacs/top"
}

# -----------------------------------------------------------------------------
# Force field metadata
# -----------------------------------------------------------------------------
# Every label known to the pipeline maps to: a family (decides the nonbonded
# and constraint settings), a force field directory and a water model.
ff_family() {
    case "$1" in
        charmm36m)  echo charmm ;;
        amber)      echo amber ;;
        gromos53a6) echo gromos ;;
        gromos54a7) echo gromos ;;
        *) die "unknown force field label '$1' (use: ${FF_LABELS[*]})" ;;
    esac
}

# AMBER_FF=auto: amber14sb when available (GROMACS >= 2026 or a port in
# forcefields/), otherwise the built-in amber99sb-ildn. Resolved once.
resolve_amber_ff() {
    [[ "$AMBER_FF" == auto ]] || return 0
    [[ -n "$GMX" ]] || detect_gmx || true
    if [[ -n "$GMX" ]] && AMBER_FF=amber14sb ff_path amber >/dev/null 2>&1; then
        AMBER_FF=amber14sb
    else
        AMBER_FF=amber99sb-ildn
        [[ "${_AMBER_WARNED:-0}" == 1 ]] || warn "AMBER_FF=auto: amber14sb not available in this GROMACS; using the built-in amber99sb-ildn (set AMBER_FF to override)"
        _AMBER_WARNED=1
    fi
}

ff_dir_name() {
    case "$1" in
        charmm36m)  echo "$CHARMM_FF" ;;
        amber)      resolve_amber_ff; echo "$AMBER_FF" ;;
        gromos53a6) echo "$GROMOS53_FF" ;;
        gromos54a7) echo "$GROMOS54_FF" ;;
        *) die "unknown force field label '$1' (use: ${FF_LABELS[*]})" ;;
    esac
}

ff_water() {
    case "$1" in
        charmm36m)  echo "$CHARMM_WATER" ;;
        amber)      echo "$AMBER_WATER" ;;
        gromos53a6|gromos54a7) echo "$GROMOS_WATER" ;;
        *) die "unknown force field label '$1'" ;;
    esac
}

# Human-readable name used in logs and figure legends (04_plot_comparison.py
# reads it from analysis/<label>/rep<N>/meta.txt).
ff_display_name() {
    case "$1" in
        charmm36m)  echo "CHARMM36m" ;;
        gromos53a6) echo "GROMOS 53A6" ;;
        gromos54a7) echo "GROMOS 54A7" ;;
        amber)
            resolve_amber_ff
            case "$AMBER_FF" in
                amber14sb*)      echo "AMBER ff14SB" ;;
                amber19sb*)      echo "AMBER ff19SB" ;;
                amber99sb-ildn*) echo "AMBER ff99SB-ILDN" ;;
                amber99sb*)      echo "AMBER ff99SB" ;;
                amber03*)        echo "AMBER ff03" ;;
                *)               echo "AMBER ($AMBER_FF)" ;;
            esac ;;
        *) echo "$1" ;;
    esac
}

# Pre-equilibrated solvent box for gmx solvate. 3-site models use spc216.gro;
# 4-site models (TIP4P, TIP4P-Ew, OPC) need tip4p.gro so that the virtual
# site is present. gmx solvate only needs the geometry; the parameters come
# from the water model chosen in pdb2gmx.
ff_solvent_box() {
    case "$(ff_water "$1")" in
        tip4p*|opc)  echo tip4p.gro ;;
        tip5p)       echo tip5p.gro ;;
        *)           echo spc216.gro ;;
    esac
}

# Absolute path of a force field directory (repo forcefields/ first, then
# any GMXLIB entry, then the GROMACS installation).
ff_path() {
    local name; name="$(ff_dir_name "$1").ff"
    local d
    IFS=':' read -r -a _dirs <<<"$GMXLIB"
    for d in "${_dirs[@]}" "$(gmx_share_top)"; do
        [[ -d "$d/$name" ]] && { echo "$d/$name"; return 0; }
    done
    return 1
}

require_ff() {
    local label=$1 path
    if ! path=$(ff_path "$label"); then
        case "$label" in
            charmm36m)
                die "CHARMM36 port '$(ff_dir_name charmm36m).ff' not found.
  Download it from https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs
  and extract it into $FF_DIR/, e.g.:
    tar -xzf charmm36-feb2026_cgenff-5.0.ff.tgz -C $FF_DIR/
  or set CHARMM_FF in config.sh to the directory name you have." ;;
            amber)
                die "AMBER force field '$(ff_dir_name amber).ff' not found (GROMACS $GMX_VERSION).
  amber14sb and amber19sb are built into GROMACS >= 2026. With an older
  GROMACS either extract a port into $FF_DIR/ and set AMBER_FF to its
  directory name, or use the built-in ff99SB-ILDN:
    AMBER_FF=amber99sb-ildn AMBER_WATER=tip3p bash scripts/01_prepare.sh amber" ;;
        esac
        die "force field '$(ff_dir_name "$label").ff' not found"
    fi
    echo "$path"
}

# Consistency checks on the force field / water combination before building.
check_ff_combo() {
    local label=$1 water; water=$(ff_water "$label")
    case "$(ff_dir_name "$label")" in
        amber19sb*)
            [[ $water == opc || $water == opc3 ]] ||
                die "ff19SB was parametrized with OPC: set AMBER_WATER=opc (or opc3), not '$water'" ;;
        amber14sb*|amber99*|amber03*)
            [[ $water == opc* ]] &&
                warn "[$label] $water with $(ff_dir_name "$label"): AMBER recommends TIP3P/TIP4P-Ew for this force field" ;;
    esac
    return 0
}

# grompp warnings tolerated per force field. GROMACS >= 2023 always emits one
# WARNING for any GROMOS force field: GROMOS was parametrized with a twin-range
# multiple-time-stepping scheme that modern integrators do not reproduce, so
# properties such as density may deviate slightly. This is a documented
# limitation of GROMOS in GROMACS (gitlab issue 2884), not a setup error, and
# it is part of what this comparison measures. Exactly one extra warning is
# allowed for GROMOS; anything else still stops the run.
ff_maxwarn() {
    case "$1" in
        gromos*) echo $(( MAXWARN + 1 )) ;;
        *)       echo "$MAXWARN" ;;
    esac
}

# Positive and negative ion molecule names for gmx genion.
#   ff_ion_names <force field dir> [topol.top]
# The ion .itp actually included by the topology is inspected first: the
# built-in AMBER ff14SB/ff19SB ports (GROMACS >= 2026) ship ion parameters
# per water model (ions_tip3p.itp, ions_opc.itp, ...), which pdb2gmx
# includes according to the chosen water. Otherwise <ff>/ions.itp is used.
# Recent ports use NA/CL; older CHARMM ports use SOD/CLA; some AMBER
# conversions use Na+/Cl-.
ff_ion_names() {
    local ffpath=$1 top=${2:-} ions="" mols f
    if [[ -n $top && -f $top ]]; then
        f=$(grep -oE '#include +"[^"]*ions[^"]*\.itp"' "$top" | head -n 1 | grep -oE '"[^"]*"' | tr -d '"' || true)
        if [[ -n $f ]]; then
            if [[ -f "$ffpath/${f#*.ff/}" ]]; then ions="$ffpath/${f#*.ff/}"
            elif [[ -f "$(dirname "$top")/$f" ]]; then ions="$(dirname "$top")/$f"
            elif [[ -f "$(gmx_share_top)/$f" ]]; then ions="$(gmx_share_top)/$f"
            fi
        fi
    fi
    [[ -n $ions ]] || ions="$ffpath/ions.itp"
    [[ -f $ions ]] || { echo "NA CL"; return; }
    mols=$(awk '/^\[ *moleculetype *\]/{getline; while($0 ~ /^ *;/ || $0 ~ /^ *$/) getline; print $1}' "$ions")
    if grep -qx NA <<<"$mols" && grep -qx CL <<<"$mols"; then
        echo "NA CL"
    elif grep -qx SOD <<<"$mols" && grep -qx CLA <<<"$mols"; then
        echo "SOD CLA"
    elif grep -qx 'Na+' <<<"$mols" && grep -qx 'Cl-' <<<"$mols"; then
        echo "Na+ Cl-"
    else
        die "could not identify Na+/Cl- molecule names in $ions"
    fi
}

# -----------------------------------------------------------------------------
# Input structure cleaning (shared by get_structure.sh and 01_prepare.sh)
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

# -----------------------------------------------------------------------------
# Misc
# -----------------------------------------------------------------------------

# Resolve command line force field arguments ("all" or labels).
parse_ff_args() {
    if [[ $# -eq 0 || "$1" == "all" ]]; then
        SELECTED_FF=("${FF_LABELS[@]}")
    else
        SELECTED_FF=("$@")
        local f
        for f in "${SELECTED_FF[@]}"; do ff_dir_name "$f" >/dev/null; done
    fi
}

# Floating point helpers (bash has integer arithmetic only).
fcalc() { awk "BEGIN{ printf \"%.0f\", $* }"; }
