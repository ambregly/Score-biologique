# Outil d'analyse R — `analyse_fusions_chromoAML.R`

Analyse seule (sans le pipeline amont) : à partir des comptages **déjà mergés**
patho + WT, elle calcule la spécificité tumorale, rapproche les fusions des
sorties **Arriba**, calcule un **score biologique** pondéré et produit les
figures + tableaux annotés (contig + k-mers inclus).

> Pour générer les comptages en amont (extract_fasta → kmerator → rdeer →
> fullnames → merge → combine), voir **`README_pipeline.md`**.

---

## Prérequis

- **R ≥ 4.1** : `tidyverse`, `scales` (obligatoires) ; `pheatmap`, `ggrepel`,
  `karyoploteR` + `GenomicRanges` (optionnels) ; `rmarkdown`, `DT` (rapport).

```r
install.packages(c("tidyverse","scales","pheatmap","ggrepel","rmarkdown","DT"))
# BiocManager::install(c("karyoploteR","GenomicRanges"))   # karyotypes
```

---

## Dossiers lus / écrits

| Rôle | Défaut (relatif) | Option CLI |
|---|---|---|
| Comptages mergés patho+WT | `JB_chromo_wt_merge1` | `--dir-merge` |
| Sorties Arriba `JB_*.tsv` | `starriba` | `--dir-arriba` |
| Contigs de départ `JB_*.fasta` | `fasta_JB` | `--dir-fasta` |
| Dossiers k-mers `JB_*_kmers/kmers.fa` | `JB_kmers` | `--dir-kmers` |
| Dossier de sortie | `analyse_fusions` | `--dir-out` |

Les défauts sont **relatifs** au dossier courant. Deux façons de faire :

**Mode A — symlinks** (défauts relatifs) :
```bash
cd ~/mon_analyse
ln -s /scratch/ambre/JB_chromo_wt_merge1 JB_chromo_wt_merge1
ln -s /data/nas/.../starriba starriba
ln -s /home/ambre/fasta_JB fasta_JB
ln -s /scratch/ambre/JB_kmers JB_kmers
Rscript /chemin/vers/analyse_fusions_chromoAML.R
```

**Mode B — chemins en CLI** :
```bash
Rscript analyse_fusions_chromoAML.R \
    --dir-merge  /scratch/ambre/JB_chromo_wt_merge1 \
    --dir-arriba /data/nas/.../starriba \
    --dir-fasta  /home/ambre/fasta_JB \
    --dir-kmers  /scratch/ambre/JB_kmers \
    --dir-out    /home/ambre/analyse_fusions
```
Le script affiche au démarrage les **chemins résolus** (absolus) réellement utilisés.

---

## Filtres (sur le comptage MAX par cohorte ; toujours actifs)

| Option | Défaut | Effet |
|---|---|---|
| `--wt-min N`    | 0 | fusion retenue si `max(WT) ≤ N` (défaut : absente des WT) |
| `--patho-min N` | 5 | fusion retenue si `max(patho) ≥ N` |

## Poids des composantes du score (0 = composante retirée)

| Option | Défaut | Composante |
|---|---|---|
| `--type N`  | 4 | type chimérique (Arriba) |
| `--conf N`  | 2 | confidence Arriba |
| `--spec N`  | 2 | spécificité WT |
| `--who N`   | 2 | fusion WHO d'intérêt |
| `--frame N` | 2 | reading frame |
| `--reads N` | 2 | reads Arriba |

Score normalisé = Σ(fraction × poids) / Σ(poids). Priorités **P1 ≥ 65 %**,
**P2 ≥ 40 %**, **P3 ≥ 20 %**. Autres : `--n-top N` (figures, défaut 30), `--help`.

```bash
Rscript analyse_fusions_chromoAML.R --who 3 --patho-min 10
Rscript analyse_fusions_chromoAML.R --frame 0 --conf 0   # score sans cadre ni confidence
```

---

## Sorties

```
analyse_fusions/
├── fusions_all_specificite_annotees.tsv       # toutes les fusions
├── fusions_chromo_specifiques_annotees.tsv     # cible : chromo-spécifiques
└── figures/
    ├── score_classement.png        # top fusions par score (noms sur l'axe)
    ├── barplot_expression.png      # top fusions par expression max
    ├── charge_par_echantillon.png  # nb de fusions par patient
    ├── repartition_types.png       # types chimériques / classes Rufflé
    ├── score_decomposition.png
    ├── heatmap_fusions.png         (si pheatmap)
    ├── karyotype_overview.png      (si karyoploteR)
    └── karyotypes/karyotype_JB_*.png   # un karyotype par patient
```

Colonnes ajoutées à chaque fusion : `contig_seq` (51 nt), `n_kmers`, et un
**k-mer par colonne** (`kmer1`, `kmer2`, …). Un même `seq_name` avec des contigs
différents selon l'échantillon est dupliqué en `seq_name.2`, `.3` (colonne
`fusion_id`) ; les figures restent sur l'identité dédupliquée.

---

## Rapport R Markdown (HTML / PDF)

Même analyse en rapport, tableaux interactifs `DT` + figures inline :
```bash
Rscript -e 'rmarkdown::render("analyse_fusions_chromoAML.Rmd")'
Rscript -e 'rmarkdown::render("analyse_fusions_chromoAML.Rmd", params=list(who=3, patho_min=10))'
```
PDF : remplacer `html_document` par `pdf_document` dans le YAML.

---

## Notes méthodologiques

- **Agrégation** des comptages : MAX par fusion.
- **Spécificité** sur le comptage max par cohorte.
- **Matching Arriba** sur la paire de gènes (alias HGNC résolus).
- **Type / classe Rufflé** : dérivés du champ `type` d'Arriba
  (`translocation`→Class 1, `inversion`→Class 4, `duplication`→Class 3,
  `deletion`/`deletion/read-through`→Class 2).
- Références : Rufflé 2017 & 2024.
