#!/usr/bin/env bash
# =============================================================================
# get_structure.sh
# Downloads a PDB entry from the RCSB and cleans it into INPUT_PDB:
#   - keeps the protein chain(s) in CHAIN, standard amino acids only
#   - first alternate location, MSE -> MET, generic HIS/CYS names
#   - removes waters, ions, ligands and hydrogens (pdb2gmx rebuilds them per
#     force field)
#   - stops on chain breaks (missing residues) unless ALLOW_GAPS=1
# The raw download is kept next to the cleaned file for reference.
# 01_prepare.sh calls this script automatically when INPUT_PDB is missing.
#
# Usage:
#   bash scripts/get_structure.sh              # PDB_ID from config.sh (1UBQ)
#   PDB_ID=1PGA bash scripts/get_structure.sh  # another entry
#   FORCE=1 bash scripts/get_structure.sh      # download again
# =============================================================================

source "$(dirname "$0")/lib.sh"
FORCE="${FORCE:-0}"

[[ -n "${PDB_ID:-}" ]] || die "PDB_ID is empty; set it in config.sh or provide INPUT_PDB yourself"
id=$(tr '[:lower:]' '[:upper:]' <<<"$PDB_ID")
[[ "$id" =~ ^[0-9][A-Z0-9]{3}$ ]] || die "'$PDB_ID' is not a 4-character PDB identifier"

raw="$(dirname "$INPUT_PDB")/${id,,}_raw.pdb"
mkdir -p "$(dirname "$INPUT_PDB")"

if [[ "$FORCE" == 1 || ! -s "$raw" ]]; then
    url="https://files.rcsb.org/download/${id}.pdb"
    log "downloading $id from $url"
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL "$url" -o "$raw" || die "download failed: $url"
    elif command -v wget >/dev/null 2>&1; then
        wget -q "$url" -O "$raw" || die "download failed: $url"
    else
        die "neither curl nor wget found; download $url manually to $raw"
    fi
    grep -q '^ATOM' "$raw" || die "$raw contains no ATOM records (wrong identifier, or an mmCIF-only entry)"
else
    log "$raw exists, skipping download (FORCE=1 to redo)"
fi

title=$(grep -E '^TITLE ' "$raw" | cut -c11- | tr -s ' ' | paste -sd' ' | sed 's/^ //')
res=$(grep -E '^REMARK   2 RESOLUTION' "$raw" | grep -oE '[0-9]+\.[0-9]+ *ANGSTROMS' | head -n 1)
log "$id: ${title:-?} (resolution ${res:-?})"

clean_pdb "$raw" "$INPUT_PDB"
check_gaps "$INPUT_PDB"

nres=$(awk 'substr($0,13,4)==" CA "' "$INPUT_PDB" | wc -l)
chains=$(awk 'substr($0,1,4)=="ATOM" {print substr($0,22,1)}' "$INPUT_PDB" | sort -u | paste -sd,)
first=$(awk 'substr($0,13,4)==" CA " {print substr($0,18,3) substr($0,23,4)+0; exit}' "$INPUT_PDB")
last=$(awk 'substr($0,13,4)==" CA " {r=substr($0,18,3) substr($0,23,4)+0} END{print r}' "$INPUT_PDB")
removed=$(awk '/^HETATM/ {print substr($0,18,3)}' "$raw" | sort | uniq -c | sort -rn | awk '{printf "%s(%s) ", $2, $1}')
log "cleaned structure: $INPUT_PDB | chain(s) $chains | $nres residues ($first to $last)"
[[ -n "$removed" ]] && log "removed HETATM groups: $removed"
log "done. Next: bash scripts/01_prepare.sh"
