#!/usr/bin/env Rscript
# =============================================================================
# Analyse fusions — projet Chromoanagenesis-AML (JB_GAILLARD)
# -----------------------------------------------------------------------------
#   1. comptages mergés chromo + WT  -> spécificité tumorale (chromo vs WT)
#   2. matching Arriba (clé = paire de breakpoints) -> caractéristiques
#   3. score biologique adapté (A,B,C,[E],G,H)  -> priorités P1/P2/P3
#   4. figures : barplot, volcano, heatmap, décomposition du score, karyotype
#
# Inspiré du pipeline LAM — sans survie ni Cox.
# Lancement :  Rscript analyse_fusions_chromoAML.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ── PARAMÈTRES ───────────────────────────────────────────────────────────────
DIR_MERGE  <- "/scratch/ambre/JB_chromo_wt_merge1"
DIR_ARRIBA <- "/data/nas/projects/2025/JB_GAILLARD/analysis/trimmed/starriba"
DIR_OUT    <- "/scratch/ambre/analyse_fusions"
DIR_FIG    <- file.path(DIR_OUT, "figures")

COUNT_MIN        <- 0        # échantillon "positif" si comptage > COUNT_MIN
MAX_WT_POS       <- 2        # Règle 2 : WT positifs maximum
MIN_CHROMO_POS   <- 5        # Règle 2 : chromo positifs minimum
MAX_FREQ_WT      <- 0.05     # Règle 2 : fréquence WT maximum
READTHROUGH_DIST <- 300000   # seuil read-through (Rufflé 2024) en pb
N_TOP            <- 30        # nb de fusions affichées dans les figures

# Composante E (optionnelle) : gènes/fusions drivers d'intérêt.
# Laisser vide, ou remplir ex : c("KMT2A","RUNX1","MECOM","NUP98")
GENES_INTEREST <- character(0)

dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)
theme_set(theme_bw(base_size = 12))

# ── Constantes visuelles ─────────────────────────────────────────────────────
TYPE_COLORS <- c(
  Translocation = "#d62728", Inversion = "#9467bd", "Délétion" = "#ff7f0e",
  Duplication = "#1f77b4", "Read-through" = "#2ca02c",
  "Non confirmé Arriba" = "grey65", Inconnu = "grey50"
)
TYPE_SHAPES <- c(
  Translocation = 16, Inversion = 17, Duplication = 15, "Délétion" = 18,
  "Read-through" = 25, "Non confirmé Arriba" = 4, Inconnu = 3
)
# Classe Rufflé d'après le type chimérique
TYPE_CLASS <- c(
  Translocation = "Class 1", Inversion = "Class 4", Duplication = "Class 3",
  "Délétion" = "Class 2", "Read-through" = "Class 2",
  "Non confirmé Arriba" = "NA", Inconnu = "NA"
)

# ── 1. CHARGEMENT DES COMPTAGES MERGÉS ───────────────────────────────────────
merge_files <- list.files(DIR_MERGE, pattern = "chromo_wt_merge1.*\\.tsv$",
                          full.names = TRUE)
if (length(merge_files) == 0)
  stop("Aucun fichier *chromo_wt_merge1*.tsv dans : ", DIR_MERGE)
cat(length(merge_files), "fichier(s) de comptage trouvé(s)\n")

raw_stacked <- map_dfr(merge_files, ~ read_tsv(.x, show_col_types = FALSE) %>%
                         mutate(source_file = basename(.x)))

all_cols    <- setdiff(names(raw_stacked), c("seq_name", "source_file"))
chromo_cols <- all_cols[str_starts(all_cols, "JB_")]
wt_cols     <- setdiff(all_cols, chromo_cols)
cat(length(chromo_cols), "colonnes chromo |", length(wt_cols), "colonnes WT\n")

raw_stacked <- raw_stacked %>%
  mutate(across(all_of(all_cols), ~ suppressWarnings(as.numeric(.x))))

# ── 2. AGRÉGATION MAX PAR FUSION ─────────────────────────────────────────────
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
mat_chromo <- as.matrix(agg[, chromo_cols]); rownames(mat_chromo) <- agg$seq_name
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
mm <- str_match(spec$seq_name,
                "^(.*?)_([A-Za-z0-9.]+:[0-9]+)_(.*)_([A-Za-z0-9.]+:[0-9]+)$")
parsed <- spec %>%
  mutate(gene1 = mm[, 2], bp1 = str_remove(mm[, 3], "^chr"),
         gene2 = mm[, 4], bp2 = str_remove(mm[, 5], "^chr")) %>%
  mutate(bp_key = map2_chr(bp1, bp2, ~ paste(sort(c(.x, .y)), collapse = "|")))
n_unparsed <- sum(is.na(parsed$bp1))
if (n_unparsed > 0)
  warning(n_unparsed, " fusion(s) au seq_name non reconnu")

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
  rename(gene1 = `#gene1`) %>%
  mutate(
    chr1 = str_remove(str_split_fixed(breakpoint1, ":", 2)[, 1], "^chr"),
    pos1 = suppressWarnings(as.numeric(str_split_fixed(breakpoint1, ":", 2)[, 2])),
    chr2 = str_remove(str_split_fixed(breakpoint2, ":", 2)[, 1], "^chr"),
    pos2 = suppressWarnings(as.numeric(str_split_fixed(breakpoint2, ":", 2)[, 2])),
    strand1_gene = str_split_fixed(`strand1(gene/fusion)`, "/", 2)[, 1],
    strand2_gene = str_split_fixed(`strand2(gene/fusion)`, "/", 2)[, 1],
    bp_key   = map2_chr(paste0(chr1, ":", pos1), paste0(chr2, ":", pos2),
                        ~ paste(sort(c(.x, .y)), collapse = "|")),
    total_reads = suppressWarnings(as.numeric(split_reads1) + as.numeric(split_reads2)),
    conf_num = case_when(
      str_detect(confidence, regex("high",   ignore_case = TRUE)) ~ 3L,
      str_detect(confidence, regex("medium", ignore_case = TRUE)) ~ 2L,
      str_detect(confidence, regex("low",    ignore_case = TRUE)) ~ 1L,
      TRUE ~ 0L),
    type_chimerique = case_when(
      chr1 != chr2 ~ "Translocation",
      chr1 == chr2 & strand1_gene != strand2_gene ~ "Inversion",
      chr1 == chr2 & strand1_gene == strand2_gene & pos1 > pos2 ~ "Duplication",
      chr1 == chr2 & strand1_gene == strand2_gene &
        (pos2 - pos1) < READTHROUGH_DIST ~ "Read-through",
      chr1 == chr2 & strand1_gene == strand2_gene &
        (pos2 - pos1) >= READTHROUGH_DIST ~ "Délétion",
      TRUE ~ "Inconnu"))

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
    n_arriba = n_distinct(sample),
    samples_arriba = paste(sort(unique(sample)), collapse = ", "),
    .groups = "drop") %>%
  mutate(total_reads = ifelse(is.finite(total_reads), total_reads, NA_real_),
         class_ruffle = unname(TYPE_CLASS[type_chimerique]))

# ── 6. JOINTURE ──────────────────────────────────────────────────────────────
annot <- parsed %>%
  left_join(arriba_sum, by = "bp_key") %>%
  mutate(arriba_matched = !is.na(confidence),
         type_chimerique = replace_na(type_chimerique, "Non confirmé Arriba"))

# ── 7. SCORE BIOLOGIQUE ADAPTÉ ───────────────────────────────────────────────
MAX_SCORE <- 12L + if (length(GENES_INTEREST) > 0) 2L else 0L
annot <- annot %>%
  mutate(
    score_A = case_when(
      type_chimerique %in% c("Translocation", "Inversion") ~ 4L,
      type_chimerique %in% c("Délétion", "Duplication")    ~ 2L,
      type_chimerique == "Read-through"                    ~ 1L,
      TRUE ~ 0L),
    score_B = case_when(
      str_detect(coalesce(confidence, ""), regex("high",   ignore_case = TRUE)) ~ 2L,
      str_detect(coalesce(confidence, ""), regex("medium", ignore_case = TRUE)) ~ 1L,
      TRUE ~ 0L),
    score_C = case_when(n_wt_pos == 0 ~ 2L, n_wt_pos == 1 ~ 1L, TRUE ~ 0L),
    score_E = if_else(gene1 %in% GENES_INTEREST | gene2 %in% GENES_INTEREST, 2L, 0L),
    score_G = case_when(reading_frame == "in-frame" ~ 2L,
                        reading_frame == "out-of-frame" ~ 0L, TRUE ~ 1L),
    score_H = case_when(coalesce(total_reads, 0) >= 50 ~ 2L,
                        coalesce(total_reads, 0) >= 10 ~ 1L, TRUE ~ 0L),
    score_total = score_A + score_B + score_C + score_E + score_G + score_H,
    score_norm  = score_total / MAX_SCORE,
    priorite = factor(case_when(
      score_norm >= 0.65 ~ "P1", score_norm >= 0.40 ~ "P2",
      score_norm >= 0.20 ~ "P3", TRUE ~ "NP"),
      levels = c("P1", "P2", "P3", "NP"))
  )

# ── 8. SORTIES TABLES ────────────────────────────────────────────────────────
out_cols <- c(
  "seq_name", "gene1", "bp1", "gene2", "bp2", "specificite",
  "score_norm", "priorite", "score_A", "score_B", "score_C",
  "score_E", "score_G", "score_H",
  "n_chromo_pos", "freq_chromo", "n_wt_pos", "freq_wt", "max_chromo", "max_wt",
  "samples_chromo", "n_fichiers",
  "arriba_matched", "type_chimerique", "class_ruffle", "reading_frame",
  "confidence", "arriba_type", "site1", "site2",
  "split_reads1", "split_reads2", "total_reads", "coverage1", "coverage2",
  "retained_protein_domains", "transcript_id1", "transcript_id2",
  "direction1", "direction2", "tags",
  "a_gene1", "a_gene2", "a_breakpoint1", "a_breakpoint2",
  "n_arriba", "samples_arriba")
annot_out <- annot %>% select(any_of(out_cols)) %>% arrange(desc(score_norm))

write_tsv(annot_out, file.path(DIR_OUT, "fusions_all_specificite_annotees.tsv"))
chromo_spec <- annot_out %>%
  filter(specificite %in% c("Chromo-spécifique", "Quasi-spécifique"))
write_tsv(chromo_spec, file.path(DIR_OUT, "fusions_chromo_specifiques_annotees.tsv"))

# ── 9. FIGURES ───────────────────────────────────────────────────────────────
has_pheatmap <- requireNamespace("pheatmap",    quietly = TRUE)
has_repel    <- requireNamespace("ggrepel",     quietly = TRUE)
has_karyo    <- requireNamespace("karyoploteR", quietly = TRUE) &&
                requireNamespace("GenomicRanges", quietly = TRUE)

fig_df <- annot %>%
  filter(specificite %in% c("Chromo-spécifique", "Quasi-spécifique")) %>%
  mutate(fusion_label = paste(gene1, gene2, sep = "--"),
         type_chimerique = factor(type_chimerique, levels = names(TYPE_COLORS)))

# 9a. Barplot — top N par expression max, coloré par type chimérique
bar_df <- fig_df %>% slice_max(max_chromo, n = N_TOP, with_ties = FALSE) %>%
  mutate(fusion_ord = reorder(fusion_label, max_chromo))
p_bar <- ggplot(bar_df, aes(fusion_ord, max_chromo, fill = type_chimerique)) +
  geom_col() +
  geom_text(aes(label = paste0("n=", n_chromo_pos)), hjust = -0.1, size = 2.5) +
  coord_flip() +
  scale_fill_manual(values = TYPE_COLORS, drop = FALSE, name = "Type chimérique") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = paste0("Top ", N_TOP, " fusions chromo-spécifiques — expression"),
       subtitle = "n = nb d'échantillons JB positifs",
       x = NULL, y = "Comptage k-mers max (chromo)") +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(DIR_FIG, "barplot_fusions.png"), p_bar, width = 11, height = 9, dpi = 150)

# 9b. Volcano — score vs fréquence chromo
vol_df <- fig_df %>% mutate(categorie = case_when(
  specificite == "Chromo-spécifique" & priorite == "P1" ~ "Chromo-spéc. P1",
  specificite == "Chromo-spécifique"                    ~ "Chromo-spéc.",
  TRUE                                                  ~ "Quasi-spéc."))
p_vol <- ggplot(vol_df, aes(score_norm, freq_chromo,
                            shape = type_chimerique, color = categorie)) +
  geom_vline(xintercept = 0.65, linetype = "dashed",  color = "grey30") +
  geom_vline(xintercept = 0.40, linetype = "dotted",  color = "grey45") +
  geom_vline(xintercept = 0.20, linetype = "dotdash", color = "grey60") +
  geom_point(size = 2.6, alpha = 0.85, stroke = 0.4) +
  scale_shape_manual(values = TYPE_SHAPES, drop = FALSE, name = "Type chimérique") +
  scale_color_manual(values = c("Chromo-spéc. P1" = "#d62728",
                                "Chromo-spéc." = "#2ca02c",
                                "Quasi-spéc." = "#1f77b4"), name = "Catégorie") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Score biologique vs fréquence — fusions chromo-spécifiques",
       subtitle = paste0("Score / ", MAX_SCORE, " pts  |  seuils P1=65% P2=40% P3=20%"),
       x = "Score biologique normalisé", y = "Fréquence cohorte chromo")
if (has_repel)
  p_vol <- p_vol + ggrepel::geom_text_repel(
    data = vol_df %>% slice_max(score_norm, n = 15, with_ties = FALSE),
    aes(label = fusion_label), size = 2.5, color = "grey20",
    max.overlaps = 20, show.legend = FALSE)
ggsave(file.path(DIR_FIG, "volcano_fusions.png"), p_vol, width = 11, height = 7, dpi = 150)

# 9c. Décomposition du score — top N
COMP_LABELS <- c(score_A = "A — Type", score_B = "B — Conf.", score_C = "C — Spéc. WT",
                 score_E = "E — Gènes", score_G = "G — Cadre", score_H = "H — Reads")
dec_df <- fig_df %>% slice_max(score_norm, n = N_TOP, with_ties = FALSE) %>%
  mutate(fusion_ord = reorder(fusion_label, score_norm)) %>%
  pivot_longer(c(score_A, score_B, score_C, score_E, score_G, score_H),
               names_to = "comp", values_to = "val") %>%
  mutate(comp = factor(COMP_LABELS[comp], levels = COMP_LABELS))
p_dec <- ggplot(dec_df, aes(val, fusion_ord, fill = comp)) +
  geom_col(width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(title = "Décomposition du score biologique",
       x = paste0("Points cumulés (max = ", MAX_SCORE, ")"), y = NULL) +
  theme(axis.text.y = element_text(size = 7), legend.position = "bottom")
ggsave(file.path(DIR_FIG, "score_decomposition.png"), p_dec, width = 11, height = 9, dpi = 150)

# 9d. Heatmap présence/absence (top score) à travers les échantillons JB
if (has_pheatmap) {
  bin <- (mat_chromo > COUNT_MIN) * 1L
  top_seq <- fig_df %>% slice_max(score_norm, n = 40, with_ties = FALSE) %>% pull(seq_name)
  bin_top <- bin[rownames(bin) %in% top_seq, , drop = FALSE]
  if (nrow(bin_top) > 1) {
    lab_map <- setNames(fig_df$fusion_label, fig_df$seq_name)
    rownames(bin_top) <- make.unique(unname(lab_map[rownames(bin_top)]))
    bin_top <- bin_top[rowSums(bin_top) > 0, , drop = FALSE]
    ar <- fig_df %>% filter(seq_name %in% top_seq) %>%
      mutate(rn = make.unique(fusion_label)) %>%
      transmute(rn, Type = as.character(type_chimerique),
                Spécificité = specificite, Priorité = as.character(priorite)) %>%
      distinct(rn, .keep_all = TRUE) %>% column_to_rownames("rn")
    ar <- ar[rownames(bin_top), , drop = FALSE]
    pheatmap::pheatmap(
      bin_top, color = colorRampPalette(c("#f7f7f7", "#111827"))(2),
      annotation_row = ar, cluster_cols = TRUE, legend = FALSE,
      fontsize_row = 7, fontsize_col = 7,
      main = "Présence/absence — fusions chromo-spécifiques (top score)",
      filename = file.path(DIR_FIG, "heatmap_fusions.png"),
      width = 10, height = 11)
  }
}

# 9e. Karyotype d'ensemble (breakpoints des fusions confirmées Arriba)
if (has_karyo) {
  ks <- fig_df %>% filter(arriba_matched) %>%
    transmute(
      chr1 = paste0("chr", str_remove(str_split_fixed(bp1, ":", 2)[, 1], "^chr")),
      pos1 = suppressWarnings(as.numeric(str_split_fixed(bp1, ":", 2)[, 2])),
      chr2 = paste0("chr", str_remove(str_split_fixed(bp2, ":", 2)[, 1], "^chr")),
      pos2 = suppressWarnings(as.numeric(str_split_fixed(bp2, ":", 2)[, 2])),
      type_chimerique = as.character(type_chimerique)) %>%
    filter(!is.na(pos1), !is.na(pos2))
  if (nrow(ks) > 0) {
    png(file.path(DIR_FIG, "karyotype_overview.png"),
        width = 1300, height = 1000, res = 120)
    kp <- karyoploteR::plotKaryotype(genome = "hg38",
            main = "Breakpoints des fusions chromo-spécifiques (Arriba)")
    for (i in seq_len(nrow(ks))) {
      r <- ks[i, ]
      col <- TYPE_COLORS[r$type_chimerique]; if (is.na(col)) col <- "grey50"
      try(karyoploteR::kpPlotLinks(
        kp,
        data  = GenomicRanges::GRanges(r$chr1, IRanges::IRanges(r$pos1, width = 1)),
        data2 = GenomicRanges::GRanges(r$chr2, IRanges::IRanges(r$pos2, width = 1)),
        col   = scales::alpha(col, 0.55)), silent = TRUE)
    }
    dev.off()
  }
}

# ── 10. RÉSUMÉ CONSOLE ───────────────────────────────────────────────────────
cat("\n=== RÉSUMÉ ===\n")
cat("Fusions totales (pool) :", nrow(annot_out), "\n")
print(count(annot_out, specificite))
cat("\nChromo/quasi-spécifiques :", nrow(chromo_spec),
    "| matchées Arriba :", sum(chromo_spec$arriba_matched), "\n")
cat("\nPriorités (chromo/quasi-spécifiques) :\n")
print(count(chromo_spec, priorite))
cat("\nFigures écrites dans :", DIR_FIG, "\n")
cat("Tables  écrites dans :", DIR_OUT, "\n")
if (!has_pheatmap) cat("(pheatmap absent → heatmap non générée)\n")
if (!has_karyo)    cat("(karyoploteR/GenomicRanges absents → karyotype non généré)\n")
