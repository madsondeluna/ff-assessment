#!/usr/bin/env bash
# =============================================================================
# check.sh
# Pre-flight check of everything the pipeline needs, with instructions for
# whatever is missing. Run it first on a new machine; run_all.sh runs it
# automatically and stops if something required is absent.
#
#   - GROMACS: executable (autodetected), version, GPU support
#   - force fields: one line per label in FF_LABELS, found / missing, and
#     for a missing one where to download it and the exact commands
#   - Python and its packages, curl or wget, nproc, free disk
#   - input structure (present, or will be downloaded)
#
# Usage:
#   bash scripts/check.sh            # exit 0 when everything required is there
#   bash scripts/check.sh --fix-hint # same, printed hints only (never fails)
# =============================================================================

source "$(dirname "$0")/lib.sh"

MISSING=0
ok()   { printf '  [ok]      %s\n' "$*"; }
miss() { printf '  [MISSING] %s\n' "$*"; MISSING=1; }
note() { printf '  [note]    %s\n' "$*"; }
hint() { printf '            %s\n' "$*"; }

echo "ff-assessment pre-flight check ($(date '+%Y-%m-%d %H:%M'))"
echo "repository: $REPO_ROOT"
echo

# -----------------------------------------------------------------------------
echo "GROMACS"
if detect_gmx && command -v "$GMX" >/dev/null 2>&1; then
    ver=$("$GMX" -quiet --version 2>/dev/null | awk -F': *' '/^GROMACS version/{print $2; exit}')
    major=$(sed -E 's/^([0-9]{4}).*/\1/' <<<"$ver")
    if [[ "$major" =~ ^[0-9]{4}$ ]] && (( major >= 2023 )); then
        ok "$GMX ($(command -v "$GMX")), version $ver"
    else
        miss "$GMX is version '$ver'; GROMACS >= 2023 is required (gmx dssp, C-rescale, Verlet-only)"
        hint "Ubuntu 24.04 ships 2023.3 (sudo apt install gromacs); conda-forge has newer builds."
    fi
    gpu=$("$GMX" -quiet --version 2>/dev/null | awk -F': *' '/^GPU support/{print $2; exit}')
    mpi=$("$GMX" -quiet --version 2>/dev/null | awk -F': *' '/^MPI library/{print $2; exit}')
    prec=$("$GMX" -quiet --version 2>/dev/null | awk -F': *' '/^Precision/{print $2; exit}')
    note "GPU support: ${gpu:-?} | MPI: ${mpi:-?} | precision: ${prec:-?}"
    if [[ "${gpu:-disabled}" == disabled ]]; then
        hint "No GPU support in this build: the runs use the CPU only (slower, but fine)."
        hint "Leave MDRUN_FLAGS empty. A CUDA build comes from conda-forge (gromacs=*=cuda*) or from source."
    else
        hint "To use the GPU: export MDRUN_FLAGS=\"-nb gpu -pme gpu -bonded gpu\" (check nvidia-smi first)."
    fi
    GMX_OK=1
    resolve_amber_ff      # AMBER_FF=auto -> amber14sb or amber99sb-ildn, once
else
    miss "no GROMACS executable found (gmx, gmx_mpi, gmx_d, gmx_mpi_d)"
    hint "Ubuntu/Debian: sudo apt install gromacs"
    hint "conda:         conda install -c conda-forge gromacs"
    hint "module system: module load gromacs   (then rerun this check)"
    hint "or point GMX to the executable: GMX=/path/to/gmx bash scripts/check.sh"
    GMX_OK=0
fi
echo

# -----------------------------------------------------------------------------
echo "Force fields (FF_LABELS: ${FF_LABELS[*]})"
for label in "${FF_LABELS[@]}"; do
    name=$(ff_dir_name "$label")
    water=$(ff_water "$label")
    if [[ $GMX_OK == 1 ]] && path=$(ff_path "$label" 2>/dev/null); then
        ok "$label: $(ff_display_name "$label") -> $path (water: $water)"
        case "$label" in
            charmm36m)
                [[ "$name" == *lj-pme* ]] && \
                    note "this is the LJ-PME variant of the port; the pipeline needs the standard (force-switch) one." ;;
            amber)
                [[ "$name" == amber99sb-ildn ]] && \
                    note "ff14SB is not in this GROMACS (built in from 2026); using ff99SB-ILDN, which ships with every version." ;;
        esac
        continue
    fi
    miss "$label: directory '$name.ff' not found in $FF_DIR/ nor in the GROMACS installation"
    case "$label" in
        charmm36m)
            hint "The CHARMM36 port is not distributed with GROMACS. Download it by hand:"
            hint "  1. open https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs in a browser"
            hint "     (registration with an e-mail address is required)"
            hint "  2. download the most recent 'charmm36-<month><year>.ff.tgz' (NOT the lj-pme variant)"
            hint "  3. copy it to this machine and extract it into $FF_DIR/:"
            hint "       tar -xzf ~/Downloads/charmm36-*.ff.tgz -C $REPO_ROOT/$FF_DIR/"
            hint "  4. check the directory name:  ls $FF_DIR/"
            hint "     if it is not '$name.ff', set CHARMM_FF to the name without '.ff', e.g."
            hint "       export CHARMM_FF=charmm36-jul2022    (or edit config.sh)"
            hint "  To run without CHARMM for now: FF_LABELS=\"amber gromos53a6 gromos54a7\" bash scripts/run_all.sh"
            ;;
        amber)
            hint "amber14sb / amber19sb are built into GROMACS >= 2026 only. Options:"
            hint "  a) use the built-in ff99SB-ILDN (every version):  export AMBER_FF=amber99sb-ildn"
            hint "     (AMBER_FF=auto, the default, does this by itself; you set AMBER_FF explicitly)"
            hint "  b) extract a GROMACS port of ff14SB into $FF_DIR/ and set AMBER_FF to its directory"
            hint "     name without '.ff' (ports are distributed by their authors, not by GROMACS)"
            hint "  c) install GROMACS >= 2026 (conda-forge)"
            ;;
        gromos53a6|gromos54a7)
            hint "$name ships with every GROMACS; if it is missing the installation is incomplete"
            hint "(on Ubuntu the data files are in the package 'gromacs-data')."
            ;;
    esac
done
echo

# -----------------------------------------------------------------------------
echo "Python and tools"
if command -v python3 >/dev/null 2>&1; then
    pyv=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
    ok "python3 $pyv"
    missing_pkgs=$(python3 - <<'EOF'
import importlib.util
print(" ".join(p for p in ("numpy", "matplotlib") if importlib.util.find_spec(p) is None))
EOF
)
    if [[ -z "$missing_pkgs" ]]; then
        ok "numpy, matplotlib"
    else
        miss "python packages: $missing_pkgs"
        hint "python3 -m venv ~/.venvs/ffa && source ~/.venvs/ffa/bin/activate && pip install -r requirements.txt"
        hint "(or: pip install --user -r requirements.txt)"
    fi
    python3 -c 'import scipy' 2>/dev/null && ok "scipy (p-values in summary_table.csv)" || \
        note "scipy not installed: summary_table.csv will have no p-values (optional)"
else
    miss "python3"
    hint "sudo apt install python3 python3-venv"
fi
if command -v curl >/dev/null 2>&1 || command -v wget >/dev/null 2>&1; then
    ok "curl or wget (structure download)"
else
    [[ -s "$INPUT_PDB" ]] && note "neither curl nor wget, but $INPUT_PDB is already there" || \
        { miss "curl or wget (needed to download $PDB_ID)"; hint "sudo apt install curl, or download https://files.rcsb.org/download/${PDB_ID}.pdb by hand to input/${PDB_ID,,}_raw.pdb"; }
fi
ok "$(nproc) CPU cores (mdrun will use NT=$NT threads)"
free_gb=$(df -BG --output=avail "$REPO_ROOT" 2>/dev/null | tail -n 1 | tr -dc '0-9')
need_gb=$(( ${#FF_LABELS[@]} * N_REPLICAS ))     # about 1 GB per 100 ns replica of ubiquitin
if [[ -n "$free_gb" ]]; then
    if (( free_gb >= need_gb + 2 )); then ok "${free_gb} GB free (about ${need_gb} GB needed)"; else miss "${free_gb} GB free; about ${need_gb} GB needed for ${#FF_LABELS[@]} force fields x $N_REPLICAS replicas"; fi
fi
command -v nvidia-smi >/dev/null 2>&1 && note "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -n 1)" || note "nvidia-smi not found (no NVIDIA GPU, or driver not installed)"
echo

# -----------------------------------------------------------------------------
echo "Input structure"
if [[ -s "$INPUT_PDB" ]]; then
    ok "$INPUT_PDB ($(awk 'substr($0,13,4)==" CA "' "$INPUT_PDB" | wc -l) residues)"
elif [[ -n "${PDB_ID:-}" ]]; then
    ok "$INPUT_PDB will be downloaded from the RCSB (PDB_ID=$PDB_ID) by 01_prepare.sh"
else
    miss "$INPUT_PDB not found and PDB_ID is empty"
    hint "set PDB_ID (e.g. PDB_ID=1UBQ) or copy a cleaned PDB file to $INPUT_PDB"
fi
echo

# -----------------------------------------------------------------------------
echo "Plan: ${#FF_LABELS[@]} force field(s) x $N_REPLICAS replica(s) x ${PROD_NS} ns production, one run at a time"
if [[ $MISSING == 0 ]]; then
    echo "Everything required is in place. Start with:  nohup bash scripts/run_all.sh > /dev/null 2>&1 &"
    exit 0
fi
echo "Something required is missing (see [MISSING] above). Fix it and rerun: bash scripts/check.sh"
[[ "${1:-}" == --fix-hint ]] && exit 0
exit 1
