#!/bin/bash
# ──────────────────────────────────────────────────────────────────────────────
# extract_fasta.sh — extrait les jonctions de fusion Arriba au format FASTA
#
# Header  : >gene1_breakpoint1_gene2_breakpoint2
# Séquence: 26 nt en amont + 25 nt en aval du séparateur "|" = 51 nt
#
# Usage :
#   ./extract_fasta.sh fichier.tsv [sortie.fasta]   # un seul fichier
#   ./extract_fasta.sh JB_*.tsv                      # plusieurs fichiers
#   ./extract_fasta.sh -d dossier/                   # tous les JB_*.tsv d'un dossier
#
# Pour chaque entrée fichier.tsv, la sortie est fichier.fasta (sauf si un nom
# de sortie est donné explicitement avec un seul fichier en argument).
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

extract_one() {
    local input="$1"
    local output="$2"
    awk -F'\t' 'NR>1 {
        # En-tête FASTA : gene1_breakpoint1_gene2_breakpoint2
        header = $1"_"$5"_"$2"_"$6
        # Couper la colonne fusion_transcript ($28) au niveau du "|"
        split($28, parts, "|")
        left_part  = parts[1]
        right_part = parts[2]
        # 26 nt en amont du "|"
        left  = substr(left_part, length(left_part)-25, 26)
        # 25 nt en aval du "|"
        right = substr(right_part, 1, 25)
        if (length(left) == 26 && length(right) == 25) {
            seq = left right
            # Supprimer les doublons exacts (meme en-tete ET meme sequence)
            key = header SUBSEP seq
            if (key in seen_pair) next
            seen_pair[key] = 1
            # Rendre l en-tete unique en cas de collision (kmerator l exige)
            cnt[header]++
            h = (cnt[header] > 1) ? header "_" cnt[header] : header
            print ">"h
            print seq
        }
    }' "$input" > "$output"
    echo "✓ $input → $output"
}

# ── Mode dossier : -d <dossier> ───────────────────────────────────────────────
if [[ "${1:-}" == "-d" ]]; then
    dir="${2:?Usage: $0 -d <dossier>}"
    shopt -s nullglob
    files=("$dir"/JB_*.tsv)
    if [[ ${#files[@]} -eq 0 ]]; then
        echo "Aucun fichier JB_*.tsv trouvé dans : $dir" >&2
        exit 1
    fi
    for f in "${files[@]}"; do
        extract_one "$f" "${f%.tsv}.fasta"
    done
    exit 0
fi

# ── Un seul fichier avec nom de sortie explicite ──────────────────────────────
if [[ $# -eq 2 && "$1" == *.tsv && "$2" == *.fasta ]]; then
    extract_one "$1" "$2"
    exit 0
fi

# ── Un ou plusieurs fichiers (sortie auto : .tsv → .fasta) ────────────────────
if [[ $# -ge 1 ]]; then
    for f in "$@"; do
        extract_one "$f" "${f%.tsv}.fasta"
    done
    exit 0
fi

echo "Usage: $0 fichier.tsv [sortie.fasta] | $0 JB_*.tsv | $0 -d dossier/" >&2
exit 1
