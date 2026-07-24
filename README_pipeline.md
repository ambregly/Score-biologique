# Pipeline complet — `run_pipeline.sh`

Enchaîne les **7 étapes** de A à Z : des sorties Arriba jusqu'aux tableaux +
figures. Tout se règle dans un **seul bloc CONFIG** en tête de `run_pipeline.sh`
— si un chemin change, une seule ligne à modifier.

> Pour l'étape 7 (l'analyse R) utilisée seule, voir **`README_analyse.md`**.

---

## Prérequis (outils externes)

`kmerator`, `rdeer` (Reindeer), `merge-kmer.py`, le script `get_fullnames.awk`,
`python3` (+ `pandas`), et **R** avec les packages de l'analyse (cf.
`README_analyse.md`).

---

## Configuration (le seul endroit à modifier)

En tête de `run_pipeline.sh` :

```bash
WORKDIR="/scratch/ambre"          # base : toutes les sorties vont dessous

# Outils / scripts
EXTRACT="/home/ambre/extract_fasta.sh"
KMERATOR="kmerator"
RDEER="rdeer"
GET_FULLNAMES="/home/chloe/2026_Projects/2026-05_kmers_AML/get_fullnames.awk"
MERGE="/data/nas/projects/kmer-collections/github-repository/bin/merge-kmer.py"
ANALYSE_R="/scratch/ambre/score-biologique/analyse_fusions_chromoAML.R"

# Entrées externes
ARRIBA_DIR="/data/nas/projects/2025/JB_GAILLARD/analysis/trimmed/starriba"

# Reindeer / kmerator
RDEER_SRV="janis"; RDEER_PORT="12800"
INDEX_CHROMO="CHU_Montpellier-Projet_Chromoanagenesis-AML"
INDEX_WT="REMPLACER_PAR_NOM_INDEX_WT"     # ← à compléter
KMERATOR_G="0"; MERGE_M="1"
```

Les **sorties** (`FASTA_DIR`, `KMER_DIR`, `RDEER_CHROMO`, …, `OUT_DIR`) sont
dérivées de `$WORKDIR`. Changer `WORKDIR` déplace tout ; changer une seule ligne
déplace un seul dossier.

---

## Lancement

```bash
bash run_pipeline.sh          # les 7 étapes
bash run_pipeline.sh 3 4      # seulement les étapes 3 et 4
bash run_pipeline.sh 7        # relancer juste l'analyse R
```

Chaque étape crée son dossier (`mkdir -p`) et n'utilise que les variables du bloc CONFIG.

---

## Les 7 étapes

| # | Étape | Outil | Entrée → Sortie |
|---|---|---|---|
| 1 | Extraire les contigs FASTA (51 nt) | `extract_fasta.sh` | `ARRIBA_DIR/JB_*.tsv` → `FASTA_DIR` |
| 2 | K-mers spécifiques par échantillon | `kmerator` | `FASTA_DIR` → `KMER_DIR/JB_*_kmers/` |
| 3a | Requête Reindeer cohorte patho | `rdeer` | `KMER_DIR` → `RDEER_CHROMO` |
| 3b | Requête Reindeer cohorte WT | `rdeer` | `KMER_DIR` → `RDEER_WT` |
| 4 | Restaurer les noms complets | `get_fullnames.awk` | `RDEER_*` → `FULLNAMES` |
| 5 | Merge k-mers | `merge-kmer.py` | `FULLNAMES` → `MERGE_CHROMO`, `MERGE_WT` |
| 6 | Combiner patho + WT (par échantillon) | `python3/pandas` | `MERGE_*` → `COMBINED` |
| 7 | Spécificité + Arriba + score + figures | `analyse_fusions_chromoAML.R` | `COMBINED` (+ Arriba, fasta, kmers) → `OUT_DIR` |

---

## Points d'attention

- **`INDEX_WT`** est un placeholder — mets le vrai nom de ton index WT avant
  l'étape 3 (vérifie avec `rdeer list -s janis`).
- Les fichiers Arriba `_discarded` sont **ignorés** (seules les fusions
  retenues sont traitées).
- Les étapes sont indépendantes : tu peux relancer uniquement celles qui ont
  changé, tant que les sorties des étapes précédentes existent.
