#!/bin/bash
# =============================================================================
# run_pipeline.sh — pipeline fusions Chromoanagenesis-AML, étapes 1 → 7
# -----------------------------------------------------------------------------
# TOUT se règle dans le bloc CONFIG ci-dessous : si un chemin change, tu ne
# modifies QUE la variable correspondante. Le reste du script en découle.
#
# Usage :
#   bash run_pipeline.sh            # tout
#   bash run_pipeline.sh 3 4        # seulement les étapes 3 et 4
# =============================================================================
set -euo pipefail

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║ CONFIG — à adapter ici, et NULLE PART AILLEURS                             ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# ── Dossier de travail : toutes les sorties vont dessous ─────────────────────
WORKDIR="/scratch/ambre"

# ── Outils / scripts (chemins des exécutables) ───────────────────────────────
EXTRACT="/home/ambre/extract_fasta.sh"                                   # étape 1
KMERATOR="kmerator"                                                      # étape 2 (dans le PATH)
RDEER="rdeer"                                                            # étape 3 (dans le PATH)
GET_FULLNAMES="/home/chloe/2026_Projects/2026-05_kmers_AML/get_fullnames.awk"  # étape 4
MERGE="/data/nas/projects/kmer-collections/github-repository/bin/merge-kmer.py"  # étape 5
ANALYSE_R="/scratch/ambre/score-biologique/analyse_fusions_chromoAML.R"  # étape 7

# ── Entrées externes ─────────────────────────────────────────────────────────
ARRIBA_DIR="/data/nas/projects/2025/JB_GAILLARD/analysis/trimmed/starriba"  # JB_*.tsv Arriba

# ── Paramètres Reindeer / kmerator ───────────────────────────────────────────
RDEER_SRV="janis"
RDEER_PORT="12800"
INDEX_CHROMO="CHU_Montpellier-Projet_Chromoanagenesis-AML"   # cohorte patho
INDEX_WT="REMPLACER_PAR_NOM_INDEX_WT"                        # cohorte normale
KMERATOR_G="0"          # option -G de kmerator
MERGE_M="1"             # option -m de merge-kmer.py

# ── Sorties (dérivées de WORKDIR — ne toucher que si tu veux réorganiser) ────
FASTA_DIR="$WORKDIR/fasta_JB"                # étape 1 : contigs
KMER_DIR="$WORKDIR/JB_kmers"                 # étape 2 : JB_*_kmers/
RDEER_CHROMO="$WORKDIR/rdeer_results"        # étape 3a
RDEER_WT="$WORKDIR/rdeer_results_wt"         # étape 3b
FULLNAMES="$WORKDIR/rdeer_fullnames"         # étape 4
MERGE_CHROMO="$WORKDIR/JB_chromoAML_merge1"  # étape 5 (patho)
MERGE_WT="$WORKDIR/JB_wt_merge1"             # étape 5 (WT)
COMBINED="$WORKDIR/JB_chromo_wt_merge1"      # étape 6
OUT_DIR="$WORKDIR/analyse_fusions"           # étape 7

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║ PIPELINE — normalement rien à modifier en dessous                         ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Étapes à lancer (par défaut 1..7, ou celles passées en argument)
STEPS="${*:-1 2 3 4 5 6 7}"
run() { [[ " $STEPS " == *" $1 "* ]]; }

# ── 1. Extraire les contigs FASTA des fusions Arriba ─────────────────────────
if run 1; then
  echo "### Étape 1 — extract_fasta"
  mkdir -p "$FASTA_DIR"
  for f in "$ARRIBA_DIR"/JB_*.tsv; do
    case "$f" in *_discarded.tsv) continue;; esac
    nom=$(basename "${f%.tsv}")
    bash "$EXTRACT" "$f" "$FASTA_DIR/${nom}.fasta"
  done
fi

# ── 2. kmerator (k-mers spécifiques par échantillon) ─────────────────────────
if run 2; then
  echo "### Étape 2 — kmerator"
  mkdir -p "$KMER_DIR"
  for f in "$FASTA_DIR"/JB_*.fasta; do
    base=$(basename "$f" .fasta)
    "$KMERATOR" -f "$f" -G "$KMERATOR_G" -y -o "$KMER_DIR/${base}_kmers"
  done
fi

# ── 3. Requêtes Reindeer (patho + WT) ────────────────────────────────────────
if run 3; then
  echo "### Étape 3a — rdeer index patho"
  mkdir -p "$RDEER_CHROMO"
  for d in "$KMER_DIR"/JB_*_kmers; do
    base=$(basename "$d" _kmers)
    "$RDEER" query -s "$RDEER_SRV" -p "$RDEER_PORT" -q "$d/kmers.fa" \
      -o "$RDEER_CHROMO/query_result_${base}_on_chromoAML" "$INDEX_CHROMO"
  done
  echo "### Étape 3b — rdeer index WT"
  mkdir -p "$RDEER_WT"
  for d in "$KMER_DIR"/JB_*_kmers; do
    base=$(basename "$d" _kmers)
    "$RDEER" query -s "$RDEER_SRV" -p "$RDEER_PORT" -q "$d/kmers.fa" \
      -o "$RDEER_WT/query_result_${base}_wt" "$INDEX_WT"
  done
fi

# ── 4. Restaurer les noms complets (get_fullnames.awk) ───────────────────────
if run 4; then
  echo "### Étape 4 — get_fullnames"
  mkdir -p "$FULLNAMES"
  for d in "$KMER_DIR"/JB_*_kmers; do
    base=$(basename "$d" _kmers)
    c=$(ls "$RDEER_CHROMO/query_result_${base}_on_chromoAML"* 2>/dev/null | head -1 || true)
    [[ -f "$d/kmers.fa" && -n "$c" ]] && awk -f "$GET_FULLNAMES" "$d/kmers.fa" "$c" \
      > "$FULLNAMES/query_result_${base}_on_chromoAML_fullNames.tsv"
    w=$(ls "$RDEER_WT/query_result_${base}_wt"* 2>/dev/null | head -1 || true)
    [[ -f "$d/kmers.fa" && -n "$w" ]] && awk -f "$GET_FULLNAMES" "$d/kmers.fa" "$w" \
      > "$FULLNAMES/query_result_${base}_wt_fullNames.tsv"
  done
fi

# ── 5. Merge (merge-kmer.py) ─────────────────────────────────────────────────
if run 5; then
  echo "### Étape 5 — merge-kmer"
  mkdir -p "$MERGE_CHROMO" "$MERGE_WT"
  for fn in "$FULLNAMES"/*_on_chromoAML_fullNames.tsv; do
    b=$(basename "$fn"); b=${b#query_result_}; base=${b%%_on_*}
    "$MERGE" -m "$MERGE_M" -o "$MERGE_CHROMO/query_${base}_chromoAML_merge1.tsv" "$fn"
  done
  for fn in "$FULLNAMES"/*_wt_fullNames.tsv; do
    b=$(basename "$fn"); b=${b#query_result_}; base=${b%%_wt_*}
    "$MERGE" -m "$MERGE_M" -o "$MERGE_WT/query_${base}_wt_merge1.tsv" "$fn"
  done
fi

# ── 6. Combiner patho + WT par échantillon ───────────────────────────────────
if run 6; then
  echo "### Étape 6 — combine patho+WT"
  mkdir -p "$COMBINED"
  MERGE_CHROMO="$MERGE_CHROMO" MERGE_WT="$MERGE_WT" COMBINED="$COMBINED" python3 - << 'PY'
import pandas as pd, glob, os, re
cd, wd, od = os.environ["MERGE_CHROMO"], os.environ["MERGE_WT"], os.environ["COMBINED"]
wt = {re.search(r"JB_\d+", os.path.basename(f)).group(0): f for f in glob.glob(wd + "/*JB_*")}
for fc in sorted(glob.glob(cd + "/*JB_*")):
    base = re.search(r"JB_\d+", os.path.basename(fc)).group(0)
    if base not in wt:
        print("wt manquant", base); continue
    m = pd.merge(pd.read_csv(fc, sep="\t"), pd.read_csv(wt[base], sep="\t"),
                 on="seq_name", how="outer")
    cols = [c for c in m.columns if c != "seq_name"]
    m[cols] = m[cols].fillna(0).astype(int)
    m.to_csv(f"{od}/query_{base}_chromo_wt_merge1.tsv", sep="\t", index=False)
    print("ok", base, len(m))
PY
fi

# ── 7. Analyse R (spécificité + Arriba + score + figures) ────────────────────
if run 7; then
  echo "### Étape 7 — analyse R"
  Rscript "$ANALYSE_R" \
    --dir-merge  "$COMBINED" \
    --dir-arriba "$ARRIBA_DIR" \
    --dir-fasta  "$FASTA_DIR" \
    --dir-kmers  "$KMER_DIR" \
    --dir-out    "$OUT_DIR"
fi

echo "### Terminé (étapes : $STEPS)"
