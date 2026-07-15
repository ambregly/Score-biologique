# Fusions chimériques — Chromoanagenesis-AML

Outils pour caractériser les fusions géniques **spécifiques** d'une cohorte AML
(chromoanagenèse) par rapport à des témoins normaux (WT), en croisant les
comptages k-mers (Transipedia / Reindeer) avec les sorties **Arriba**, puis en
calculant un **score biologique** et en produisant des figures.

Inspiré du pipeline LAM d'origine — sans les analyses de survie ni le Cox.

---

## Contenu du dépôt

| Fichier | Rôle |
|---|---|
| `extract_fasta.sh` | Extrait les jonctions de fusion Arriba au format FASTA (51 nt) |
| `analyse_fusions_chromoAML.R` | Outil principal en ligne de commande (spécificité + matching Arriba + score + figures) |
| `analyse_fusions_chromoAML.Rmd` | Même analyse en rapport HTML/PDF (paramétrable) |
| `README.md` | Ce fichier |

---

## Prérequis

- **R ≥ 4.1** avec : `tidyverse`, `scales` (obligatoires) ;
  `pheatmap`, `ggrepel`, `karyoploteR` + `GenomicRanges` (optionnels — figures
  correspondantes simplement sautées si absents) ; `rmarkdown`, `DT` (pour le rapport).
- Outils externes du pipeline amont : `kmerator`, `rdeer` (Reindeer),
  `merge-kmer.py`, et le script `get_fullnames.awk`.

```r
install.packages(c("tidyverse","scales","pheatmap","ggrepel","rmarkdown","DT"))
# karyoploteR (Bioconductor) :
# BiocManager::install(c("karyoploteR","GenomicRanges"))
```

---

## Pipeline complet (de A à Z)

Les étapes 1→6 préparent les données ; l'étape 7 est l'analyse R de ce dépôt.

### 1. Extraire les FASTA des jonctions de fusion
```bash
chmod +x extract_fasta.sh
mkdir -p fasta_JB
for f in /data/nas/projects/2025/JB_GAILLARD/analysis/trimmed/starriba/JB_*.tsv; do
    case "$f" in *_discarded.tsv) continue;; esac
    nom=$(basename "${f%.tsv}")
    bash extract_fasta.sh "$f" "fasta_JB/${nom}.fasta"
done
```

### 2. kmerator (k-mers spécifiques par échantillon)
```bash
mkdir -p JB_kmers && cd JB_kmers
for f in ella/home/ambre/fasta_JB/JB_*.fasta; do
    base=$(basename "$f" .fasta)
    kmerator -f "$f" -G 0 -y -o "${base}_kmers"
done
cd ..
```

### 3. Requêtes Reindeer sur les deux index
```bash
# 3a. Index patho (chromoAML)
mkdir -p rdeer_results
for d in scratch/ambre/JB_kmers/JB_*_kmers; do
    base=$(basename "$d" _kmers)
    rdeer query -s janis -p 12800 -q "$d/kmers.fa" \
        -o "rdeer_results/query_result_${base}_on_chromoAML" \
        CHU_Montpellier-Projet_Chromoanagenesis-AML
done

# 3b. Index WT (normaux) — remplacer <NOM_INDEX_WT>
mkdir -p rdeer_results_wt
for d in JB_kmers/JB_*_kmers; do
    base=$(basename "$d" _kmers)
    rdeer query -s janis -p 12800 -q "$d/kmers.fa" \
        -o "rdeer_results_wt/query_result_${base}_wt" \
        <NOM_INDEX_WT>
done
```

### 4. Restaurer les noms complets (`get_fullnames.awk`)
```bash
AWK=/home/chloe/2026_Projects/2026-05_kmers_AML/get_fullnames.awk
mkdir -p rdeer_fullnames
for d in /scratch/ambre/JB_kmers/JB_*_kmers; do
    base=$(basename "$d" _kmers)
    c=$(ls rdeer_results/query_result_${base}_on_chromoAML* 2>/dev/null | head -1)
    [ -n "$c" ] && awk -f "$AWK" "$d/kmers.fa" "$c" \
        > "rdeer_fullnames/query_result_${base}_on_chromoAML_fullNames.tsv"
    c=$(ls rdeer_results_wt/query_result_${base}_wt* 2>/dev/null | head -1)
    [ -n "$c" ] && awk -f "$AWK" "$d/kmers.fa" "$c" \
        > "rdeer_fullnames/query_result_${base}_wt_fullNames.tsv"
done
```

### 5. Merge (`merge-kmer.py`)
```bash
MERGE=/data/nas/projects/kmer-collections/github-repository/bin/merge-kmer.py
mkdir -p JB_chromoAML_merge1 JB_wt_merge1
for fn in rdeer_fullnames/*_on_chromoAML_fullNames.tsv; do
    b=$(basename "$fn"); b=${b#query_result_}; base=${b%%_on_*}
    "$MERGE" -m 1 -o "JB_chromoAML_merge1/query_${base}_chromoAML_merge1.tsv" "$fn"
done
for fn in rdeer_fullnames/*_wt_fullNames.tsv; do
    b=$(basename "$fn"); b=${b#query_result_}; base=${b%%_wt_*}
    "$MERGE" -m 1 -o "JB_wt_merge1/query_${base}_wt_merge1.tsv" "$fn"
done
```

### 6. Combiner patho + WT par échantillon
```bash
python3 - << 'PY'
import pandas as pd, glob, os, re
os.makedirs("JB_chromo_wt_merge1", exist_ok=True)
wt={re.search(r"JB_\d+",os.path.basename(f)).group(0):f for f in glob.glob("JB_wt_merge1/*JB_*")}
for fc in sorted(glob.glob("JB_chromoAML_merge1/*JB_*")):
    base=re.search(r"JB_\d+",os.path.basename(fc)).group(0)
    if base not in wt: print("wt manquant",base); continue
    m=pd.merge(pd.read_csv(fc,sep="\t"),pd.read_csv(wt[base],sep="\t"),on="seq_name",how="outer")
    cols=[c for c in m.columns if c!="seq_name"]; m[cols]=m[cols].fillna(0).astype(int)
    m.to_csv(f"JB_chromo_wt_merge1/query_{base}_chromo_wt_merge1.tsv",sep="\t",index=False)
    print("ok",base,len(m))
PY
```

### 7. Analyse R (ce dépôt) — voir ci-dessous.

---

## Utilisation de l'outil R

### Deux façons de fournir les chemins

Le script lit **3 dossiers** et en écrit **1** :

| Rôle | Défaut (relatif) | Option CLI |
|---|---|---|
| Comptages mergés patho+WT | `JB_chromo_wt_merge1` | `--dir-merge` |
| Sorties Arriba `JB_*.tsv` | `starriba` | `--dir-arriba` |
| Contigs de départ `JB_*.fasta` | `fasta_JB` | `--dir-fasta` |
| Dossiers k-mers `JB_*_kmers/kmers.fa` | `JB_kmers` | `--dir-kmers` |
| Dossier de sortie | `analyse_fusions` | `--dir-out` |

Les tables de sortie incluent, pour chaque fusion : la séquence contig
(`contig_seq`, 51 nt), le nombre de k-mers (`n_kmers`) et **un k-mer par colonne**
(`kmer1`, `kmer2`, …), rattachés par `seq_name`.

**Mode A — dossier local + liens symboliques** (les défauts relatifs suffisent) :
```bash
mkdir -p ~/mon_analyse && cd ~/mon_analyse
ln -s /scratch/ambre/JB_chromo_wt_merge1 JB_chromo_wt_merge1
ln -s /data/nas/projects/2025/JB_GAILLARD/analysis/trimmed/starriba starriba
Rscript /chemin/vers/analyse_fusions_chromoAML.R
# -> lit ./JB_chromo_wt_merge1 et ./starriba, écrit ./analyse_fusions
```

**Mode B — tout en ligne de commande** (aucun symlink) :
```bash
Rscript analyse_fusions_chromoAML.R \
    --dir-merge  /scratch/ambre/JB_chromo_wt_merge1 \
    --dir-arriba /data/nas/projects/2025/JB_GAILLARD/analysis/trimmed/starriba \
    --dir-out    /scratch/ambre/analyse_fusions
```

> Le script affiche au démarrage les **chemins résolus** (absolus) qu'il utilise
> réellement, pour vérifier qu'il lit/écrit au bon endroit.

### Options de filtrage (toujours actives ; seuil 0 = valide)

| Option | Défaut | Effet |
|---|---|---|
| `--wt-min N` | 0 | fusion retenue si `max(WT) ≤ N` (défaut : absente des WT) |
| `--patho-min N` | 5 | fusion retenue si `max(patho) ≥ N` |

### Poids des composantes du score (0 = composante retirée)

| Option | Défaut | Composante |
|---|---|---|
| `--type N`  | 3 | type chimérique |
| `--conf N`  | 2 | confidence Arriba |
| `--spec N`  | 2 | spécificité WT |
| `--who N`   | 2 | fusion WHO d'intérêt |
| `--frame N` | 2 | reading frame |
| `--reads N` | 5 | couverture reads (split reads + paires discordantes) |

Score normalisé = Σ(fraction × poids) / Σ(poids). Priorités : **P1 ≥ 65 %**,
**P2 ≥ 40 %**, **P3 ≥ 20 %**. Retirer une composante (`--frame 0`) l'exclut du
calcul **et** du dénominateur.

> **Pondération** : la **couverture en reads** (`--reads`, défaut 5) prime
> désormais sur le **type chimérique** (`--type`, défaut 3). Un signal de fusion
> bien couvert est jugé plus fiable qu'un type « fort » (translocation/inversion)
> faiblement supporté. La composante reads a une granularité fine
> (≥ 5 / 10 / 25 / 50 / 100 reads → 0,2 / 0,4 / 0,6 / 0,8 / 1,0).

### Détermination du type chimérique

Le type est repris de la colonne `type` d'Arriba quand elle est renseignée. Si
elle est **absente**, il est reconstruit par géométrie — directions (`upstream`/
`downstream`) + ordre des breakpoints, reproduisant la logique `get_fusion_type`
d'Arriba (concordance 100 % sur les fusions testées). En dernier recours (ni type
ni direction), seuls les chromosomes distinguent une translocation d'un cas
« même chr » indéterminé. La colonne `type_source` (`arriba` / `géométrie` /
`chr` / `absent`) trace l'origine de chaque classification.

Autres : `--n-top N` (fusions dans les figures, défaut 30), `--help`.

### Exemples
```bash
Rscript analyse_fusions_chromoAML.R                       # réglages par défaut
Rscript analyse_fusions_chromoAML.R --who 3 --patho-min 10
Rscript analyse_fusions_chromoAML.R --frame 0 --conf 0   # score sans cadre ni confidence
```

### Sorties (dans `--dir-out`)
```
analyse_fusions/
├── fusions_all_specificite_annotees.tsv       # toutes les fusions
├── fusions_chromo_specifiques_annotees.tsv     # cible : chromo-spécifiques annotées
└── figures/
    ├── score_classement.png        # top fusions par score (noms sur l'axe)
    ├── barplot_expression.png      # top fusions par expression max
    ├── charge_par_echantillon.png  # nb de fusions par patient (chromoanagenèse)
    ├── repartition_types.png       # types chimériques / classes Rufflé
    ├── score_decomposition.png     # décomposition A..H du score
    ├── heatmap_fusions.png         (si pheatmap)
    ├── karyotype_overview.png      (si karyoploteR)
    └── karyotypes/                 # un karyotype par patient (si karyoploteR)
        └── karyotype_JB_*.png
```

---

## Rapport R Markdown (HTML / PDF)

Même analyse, sous forme de rapport avec tableaux interactifs et figures inline :
```bash
# HTML
Rscript -e 'rmarkdown::render("analyse_fusions_chromoAML.Rmd")'

# En modifiant des paramètres
Rscript -e 'rmarkdown::render("analyse_fusions_chromoAML.Rmd",
             params = list(who = 3, patho_min = 10))'
```
Pour un **PDF**, remplacer dans le YAML le bloc `html_document` par
`pdf_document` (indiqué en commentaire dans le fichier).

---

## Notes méthodologiques

- **Agrégation** des comptages : **MAX** par fusion (capture le k-mer de
  jonction le plus exprimé).
- **Spécificité** basée sur le comptage **max par cohorte**.
- **Matching Arriba** sur la **paire de gènes** (orientation indifférente,
  alias HGNC résolus) — les noms du `seq_name` et d'Arriba viennent de la même
  source, donc la combinaison de gènes suffit.
- **Type chimérique / classe Rufflé** inférés depuis chr, brins et distance
  (read-through < 300 kb).
- Références : Rufflé 2017 & 2024.
