#!/usr/bin/env Rscript
# =============================================================================
# Analyse fusions — projet Chromoanagenesis-AML (JB_GAILLARD)
# -----------------------------------------------------------------------------
# Étape 1 (coeur) :
#   - charge les comptages mergés chromo + WT (query_JB_*_chromo_wt_merge1.tsv)
#   - sépare les fusions CHROMO-SPÉCIFIQUES (absentes des WT) des fusions
#     présentes aussi dans les WT  (Règle 1 / Règle 2 du pipeline LAM)
#   - rapproche les fusions des sorties Arriba (clé = paire de breakpoints)
#   - récupère toutes les caractéristiques Arriba (type, confidence,
#     reading_frame, sites, reads, domaines...) + déduit le type chimérique
#     et la classe Rufflé (1–4)
#
# Inspiré du pipeline LAM (score biologique) — sans survie ni Cox.
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))

# ── PARAMÈTRES ───────────────────────────────────────────────────────────────
DIR_MERGE  <- "/scratch/ambre/JB_chromo_wt_merge1"
DIR_ARRIBA <- "/data/nas/projects/2025/JB_GAILLARD/analysis/trimmed/starriba"
DIR_OUT    <- "/scratch/ambre/analyse_fusions"

COUNT_MIN        <- 0        # un échantillon est "positif" si comptage > COUNT_MIN
MAX_WT_POS       <- 2        # Règle 2 : WT positifs maximum
MIN_CHROMO_POS   <- 5        # Règle 2 : chromo positifs minimum
MAX_FREQ_WT      <- 0.05     # Règle 2 : fréquence WT maximum
READTHROUGH_DIST <- 300000   # seuil read-through (Rufflé 2024) en pb

dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

# Classe Rufflé d'après le type chimérique (cf. TYPE_CONFIDENCE du pipeline LAM)
TYPE_CLASS <- c(
  Translocation  = "Class 1",  # autre chromosome
  Inversion      = "Class 4",  # même chr, brins opposés
  Duplication    = "Class 3",
  "Délétion"     = "Class 2",
  "Read-through" = "Class 2",  # même chr, < 300 kb
  Inconnu        = "NA"
)

# ── 1. CHARGEMENT DES COMPTAGES MERGÉS ───────────────────────────────────────
merge_files <- list.files(DIR_MERGE, pattern = "chromo_wt_merge1.*\\.tsv$",
                          full.names = TRUE)
if (length(merge_files) == 0)
  stop("Aucun fichier *chromo_wt_merge1*.tsv dans : ", DIR_MERGE)
cat(length(merge_files), "fichier(s) de comptage trouvé(s)\n")

raw_stacked <- map_dfr(merge_files, ~ read_tsv(.x, show_col_types = FALSE) %>%
                         mutate(source_file = basename(.x)))

# Colonnes chromo = commencent par "JB_" ; WT = tout le reste
all_cols    <- setdiff(names(raw_stacked), c("seq_name", "source_file"))
chromo_cols <- all_cols[str_starts(all_cols, "JB_")]
wt_cols     <- setdiff(all_cols, chromo_cols)
cat(length(chromo_cols), "colonnes chromo |", length(wt_cols), "colonnes WT\n")

# Forcer les comptages en numérique
raw_stacked <- raw_stacked %>%
  mutate(across(all_of(all_cols), ~ suppressWarnings(as.numeric(.x))))

# ── 2. AGRÉGATION MAX PAR FUSION ─────────────────────────────────────────────
# Une fusion peut apparaître dans plusieurs fichiers JB (k-mers légèrement
# différents) → on garde le comptage maximal par échantillon (Rufflé 2024).
agg <- raw_stacked %>%
  group_by(seq_name) %>%
  summarise(
    across(all_of(all_cols), ~ suppressWarnings(max(.x, na.rm = TRUE))),
    n_fichiers = n_distinct(source_file),
    .groups = "drop"
  ) %>%
  mutate(across(all_of(all_cols), ~ ifelse(is.finite(.x), .x, 0)))

cat(nrow(agg), "fusions uniques après agrégation MAX\n")

# ── 3. SPÉCIFICITÉ TUMORALE (chromo vs WT) ───────────────────────────────────
mat_chromo <- as.matrix(agg[, chromo_cols])
mat_wt     <- as.matrix(agg[, wt_cols])
n_chromo   <- length(chromo_cols)
n_wt       <- length(wt_cols)

spec <- agg %>%
  transmute(seq_name, n_fichiers) %>%
  mutate(
    n_chromo_pos   = rowSums(mat_chromo > COUNT_MIN),
    n_wt_pos       = rowSums(mat_wt     > COUNT_MIN),
    freq_chromo    = n_chromo_pos / n_chromo,
    freq_wt        = n_wt_pos / n_wt,
    max_chromo     = apply(mat_chromo, 1, max),
    max_wt         = apply(mat_wt,     1, max),
    samples_chromo = apply(mat_chromo > COUNT_MIN, 1,
                           function(r) paste(chromo_cols[r], collapse = ", ")),
    specificite = case_when(
      n_wt_pos == 0 & n_chromo_pos >= 1 ~ "Chromo-spécifique",
      n_wt_pos > 0 & n_wt_pos <= MAX_WT_POS &
        n_chromo_pos >= MIN_CHROMO_POS & freq_wt <= MAX_FREQ_WT ~ "Quasi-spécifique",
      n_wt_pos > 0 ~ "Partagée WT",
      TRUE ~ "Autre"
    )
  )

# ── 4. PARSE seq_name → gènes + breakpoints ──────────────────────────────────
# Format : gene1_chr1:pos1_gene2_chr2:pos2
# On capture les 2 jetons chr:pos (le séparateur "_" n'est jamais dans un chr).
mm <- str_match(spec$seq_name,
                "^(.*?)_([A-Za-z0-9.]+:[0-9]+)_(.*)_([A-Za-z0-9.]+:[0-9]+)$")
parsed <- spec %>%
  mutate(
    gene1 = mm[, 2],
    bp1   = str_remove(mm[, 3], "^chr"),
    gene2 = mm[, 4],
    bp2   = str_remove(mm[, 5], "^chr")
  ) %>%
  mutate(bp_key = map2_chr(bp1, bp2, ~ paste(sort(c(.x, .y)), collapse = "|")))

n_unparsed <- sum(is.na(parsed$bp1))
if (n_unparsed > 0)
  warning(n_unparsed, " fusion(s) au seq_name non reconnu (breakpoints non extraits)")

# ── 5. CHARGEMENT & ANNOTATION ARRIBA ────────────────────────────────────────
arriba_files <- list.files(DIR_ARRIBA, pattern = "^JB_.*\\.tsv$", full.names = TRUE)
arriba_files <- arriba_files[!str_detect(arriba_files, "_discarded")]
if (length(arriba_files) == 0)
  stop("Aucun fichier Arriba JB_*.tsv dans : ", DIR_ARRIBA)
cat(length(arriba_files), "fichier(s) Arriba trouvé(s)\n")

arriba <- map_dfr(arriba_files, function(f) {
  read_tsv(f, show_col_types = FALSE,
           col_types = cols(.default = col_character())) %>%
    mutate(sample = str_extract(basename(f), "JB_[0-9]+"))
}) %>%
  rename(gene1 = `#gene1`)

arriba <- arriba %>%
  mutate(
    chr1 = str_remove(str_split_fixed(breakpoint1, ":", 2)[, 1], "^chr"),
    pos1 = suppressWarnings(as.numeric(str_split_fixed(breakpoint1, ":", 2)[, 2])),
    chr2 = str_remove(str_split_fixed(breakpoint2, ":", 2)[, 1], "^chr"),
    pos2 = suppressWarnings(as.numeric(str_split_fixed(breakpoint2, ":", 2)[, 2])),
    strand1_gene = str_split_fixed(`strand1(gene/fusion)`, "/", 2)[, 1],
    strand2_gene = str_split_fixed(`strand2(gene/fusion)`, "/", 2)[, 1],
    bp1_norm = paste0(chr1, ":", pos1),
    bp2_norm = paste0(chr2, ":", pos2),
    bp_key   = map2_chr(bp1_norm, bp2_norm, ~ paste(sort(c(.x, .y)), collapse = "|")),
    total_reads = suppressWarnings(as.numeric(split_reads1) + as.numeric(split_reads2)),
    conf_num = case_when(
      str_detect(confidence, regex("high",   ignore_case = TRUE)) ~ 3L,
      str_detect(confidence, regex("medium", ignore_case = TRUE)) ~ 2L,
      str_detect(confidence, regex("low",    ignore_case = TRUE)) ~ 1L,
      TRUE ~ 0L),
    # Type chimérique inféré (Rufflé 2017)
    type_chimerique = case_when(
      chr1 != chr2 ~ "Translocation",
      chr1 == chr2 & strand1_gene != strand2_gene ~ "Inversion",
      chr1 == chr2 & strand1_gene == strand2_gene & pos1 > pos2 ~ "Duplication",
      chr1 == chr2 & strand1_gene == strand2_gene &
        (pos2 - pos1) < READTHROUGH_DIST ~ "Read-through",
      chr1 == chr2 & strand1_gene == strand2_gene &
        (pos2 - pos1) >= READTHROUGH_DIST ~ "Délétion",
      TRUE ~ "Inconnu"
    )
  )

# Résumé Arriba : une ligne par paire de breakpoints, on garde le meilleur record
arriba_sum <- arriba %>%
  arrange(desc(conf_num)) %>%
  group_by(bp_key) %>%
  summarise(
    a_gene1 = first(gene1), a_gene2 = first(gene2),
    a_breakpoint1 = first(breakpoint1), a_breakpoint2 = first(breakpoint2),
    site1 = first(site1), site2 = first(site2),
    arriba_type     = first(type),
    confidence      = first(confidence),
    type_chimerique = first(type_chimerique),
    reading_frame = case_when(
      any(reading_frame == "in-frame",     na.rm = TRUE) ~ "in-frame",
      any(reading_frame == "out-of-frame", na.rm = TRUE) ~ "out-of-frame",
      TRUE ~ "."),
    split_reads1 = first(split_reads1), split_reads2 = first(split_reads2),
    total_reads  = suppressWarnings(max(total_reads, na.rm = TRUE)),
    coverage1 = first(coverage1), coverage2 = first(coverage2),
    retained_protein_domains = first(retained_protein_domains),
    transcript_id1 = first(transcript_id1), transcript_id2 = first(transcript_id2),
    direction1 = first(direction1), direction2 = first(direction2),
    tags = first(tags),
    n_arriba       = n_distinct(sample),
    samples_arriba = paste(sort(unique(sample)), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    total_reads  = ifelse(is.finite(total_reads), total_reads, NA_real_),
    class_ruffle = unname(TYPE_CLASS[type_chimerique])
  )

# ── 6. JOINTURE & SORTIES ────────────────────────────────────────────────────
annot <- parsed %>%
  left_join(arriba_sum, by = "bp_key") %>%
  mutate(arriba_matched = !is.na(confidence))

out_cols <- c(
  "seq_name", "gene1", "bp1", "gene2", "bp2", "specificite",
  "n_chromo_pos", "freq_chromo", "n_wt_pos", "freq_wt", "max_chromo", "max_wt",
  "samples_chromo", "n_fichiers",
  "arriba_matched", "type_chimerique", "class_ruffle", "reading_frame",
  "confidence", "arriba_type", "site1", "site2",
  "split_reads1", "split_reads2", "total_reads", "coverage1", "coverage2",
  "retained_protein_domains", "transcript_id1", "transcript_id2",
  "direction1", "direction2", "tags",
  "a_gene1", "a_gene2", "a_breakpoint1", "a_breakpoint2",
  "n_arriba", "samples_arriba"
)
annot_out <- annot %>% select(any_of(out_cols))

# Table complète (toutes fusions, toutes catégories)
write_tsv(annot_out, file.path(DIR_OUT, "fusions_all_specificite_annotees.tsv"))

# Table ciblée : fusions chromo-spécifiques + quasi-spécifiques, annotées Arriba
chromo_spec <- annot_out %>%
  filter(specificite %in% c("Chromo-spécifique", "Quasi-spécifique")) %>%
  arrange(specificite, desc(n_chromo_pos))
write_tsv(chromo_spec, file.path(DIR_OUT, "fusions_chromo_specifiques_annotees.tsv"))

# ── 7. RÉSUMÉ CONSOLE ────────────────────────────────────────────────────────
cat("\n=== RÉSUMÉ ===\n")
cat("Fusions totales (pool) :", nrow(annot_out), "\n")
print(count(annot_out, specificite))
cat("\nChromo/quasi-spécifiques :", nrow(chromo_spec),
    "| dont matchées Arriba :", sum(chromo_spec$arriba_matched), "\n")
cat("\nRépartition type chimérique (chromo/quasi-spécifiques matchées) :\n")
print(chromo_spec %>% filter(arriba_matched) %>% count(type_chimerique, class_ruffle))
cat("\nSorties écrites dans :", DIR_OUT, "\n")
