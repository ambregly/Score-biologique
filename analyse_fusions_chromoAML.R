#!/usr/bin/env Rscript
# =============================================================================
# Analyse fusions — projet Chromoanagenesis-AML (JB_GAILLARD)
# -----------------------------------------------------------------------------
#   1. comptages mergés patho (JB_*) + WT  -> spécificité tumorale
#   2. matching Arriba (clé = paire de gènes) -> caractéristiques
#   3. score biologique pondéré & paramétrable -> priorités P1/P2/P3
#   4. figures : barplot, volcano, heatmap, décomposition du score, karyotype
#
# Inspiré du pipeline LAM — sans survie ni Cox.
#
# Usage : Rscript analyse_fusions_chromoAML.R [options]
#
# Chemins : par défaut RELATIFS au dossier courant (où tu lances Rscript).
# Astuce : place-toi dans un dossier de travail et crée des liens symboliques
#   ln -s /chemin/reel/JB_chromo_wt_merge1 JB_chromo_wt_merge1
#   ln -s /data/nas/.../starriba           starriba
# ...ou passe les chemins réels en CLI (--dir-merge, --dir-arriba, --dir-out).
#
# Filtres (basés sur le comptage MAX par cohorte ; toujours actifs) :
#   --wt-min N     fusion retenue si max(WT)    <= N   (défaut 0 : absente des WT)
#   --patho-min N  fusion retenue si max(patho) >= N   (défaut 5)
#
# Poids des composantes du score (0 = composante retirée du calcul) :
#   --type N   type chimérique   (défaut 4)
#   --conf N   confidence Arriba (défaut 2)
#   --spec N   spécificité WT    (défaut 2)
#   --who N    fusion WHO        (défaut 2)
#   --frame N  reading frame     (défaut 2)
#   --reads N  reads Arriba      (défaut 2)
#
# Autres : --n-top N | --dir-merge | --dir-arriba | --dir-out | --help
#
# Exemple : Rscript analyse_fusions_chromoAML.R --who 3 --patho-min 10
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ── CONFIG PAR DÉFAUT ────────────────────────────────────────────────────────
opt <- list(
  dir_merge  = "JB_chromo_wt_merge1",   # relatif au dossier courant (getwd())
  dir_arriba = "starriba",              # relatif au dossier courant (getwd())
  dir_out    = "analyse_fusions",       # relatif au dossier courant (getwd())
  wt_min = 0, patho_min = 5,          # filtres (sur le max par cohorte)
  n_top = 30,                          # figures
  w_type = 4, w_conf = 2, w_spec = 2, w_who = 2, w_frame = 2, w_reads = 2
)
READTHROUGH_DIST <- 300000   # seuil read-through (Rufflé 2024) en pb

# ── PARSER CLI ───────────────────────────────────────────────────────────────
alias <- c(type = "w_type", conf = "w_conf", confidence = "w_conf",
           spec = "w_spec", specificity = "w_spec", who = "w_who",
           frame = "w_frame", reads = "w_reads")
string_opts <- c("dir_merge", "dir_arriba", "dir_out")

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || "-h" %in% args) {
  cat("Options : --wt-min --patho-min --type --conf --spec --who --frame",
      "--reads --n-top --dir-merge --dir-arriba --dir-out\n")
  cat("Voir l'entête du script pour le détail.\n")
  quit(status = 0)
}

i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (grepl("^--", a)) {
    a2 <- sub("^--", "", a)
    if (grepl("=", a2)) {
      kv <- strsplit(a2, "=", fixed = TRUE)[[1]]; key <- kv[1]; val <- kv[2]
    } else {
      key <- a2
      if (i < length(args) && !grepl("^--", args[i + 1])) {
        val <- args[i + 1]; i <- i + 1
      } else val <- "TRUE"
    }
    key <- gsub("-", "_", key)
    if (key %in% names(alias)) key <- alias[[key]]
    if (key %in% names(opt)) {
      opt[[key]] <- if (key %in% string_opts) val else suppressWarnings(as.numeric(val))
    } else warning("Option inconnue ignorée : --", key)
  }
  i <- i + 1
}

DIR_OUT <- opt$dir_out
DIR_FIG <- file.path(DIR_OUT, "figures")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)
theme_set(theme_bw(base_size = 12))

cat("=== CONFIGURATION ===\n")
cat("Dossier courant :", getwd(), "\n")
cat("  merge  :", normalizePath(opt$dir_merge,  mustWork = FALSE), "\n")
cat("  arriba :", normalizePath(opt$dir_arriba, mustWork = FALSE), "\n")
cat("  sortie :", normalizePath(opt$dir_out,    mustWork = FALSE), "\n")
cat(sprintf("Filtres : max(WT) <= %s  |  max(patho) >= %s\n", opt$wt_min, opt$patho_min))
cat(sprintf("Poids   : type=%s conf=%s spec=%s who=%s frame=%s reads=%s\n\n",
            opt$w_type, opt$w_conf, opt$w_spec, opt$w_who, opt$w_frame, opt$w_reads))

# ── Fusions WHO d'intérêt (composante E) + alias HGNC ─────────────────────────
WHO_FUSIONS_INTEREST <- c(
  "RUNX1--RUNX1T1", "CBFB--MYH11", "PML--RARA", "KMT2A--MLLT3",
  "DEK--NUP214", "NUP98--NSD1", "BCR--ABL1", "NRIP1--MIR99AHG"
)
GENE_ALIASES <- c(
  MLL = "KMT2A", HRX = "KMT2A", ALL1 = "KMT2A", MLL1 = "KMT2A", TRX1 = "KMT2A",
  AF9 = "MLLT3", LTG9 = "MLLT3", AF6 = "MLLT4", AFDN = "MLLT4", ENL = "MLLT1",
  ETO = "RUNX1T1", MTG8 = "RUNX1T1", CBFA2T1 = "RUNX1T1", BCR1 = "BCR"
)
resolve_alias <- function(g) ifelse(g %in% names(GENE_ALIASES), GENE_ALIASES[g], g)

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
TYPE_CLASS <- c(
  Translocation = "Class 1", Inversion = "Class 4", Duplication = "Class 3",
  "Délétion" = "Class 2", "Read-through" = "Class 2",
  "Non confirmé Arriba" = "NA", Inconnu = "NA"
)

# ── 1. CHARGEMENT DES COMPTAGES MERGÉS ───────────────────────────────────────
merge_files <- list.files(opt$dir_merge, pattern = "chromo_wt_merge1.*\\.tsv$",
                          full.names = TRUE)
if (length(merge_files) == 0)
  stop("Aucun fichier *chromo_wt_merge1*.tsv dans : ", opt$dir_merge)
cat(length(merge_files), "fichier(s) de comptage trouvé(s)\n")

raw_stacked <- map_dfr(merge_files, ~ read_tsv(.x, show_col_types = FALSE) %>%
                         mutate(source_file = basename(.x)))
all_cols   <- setdiff(names(raw_stacked), c("seq_name", "source_file"))
patho_cols <- all_cols[str_starts(all_cols, "JB_")]
wt_cols    <- setdiff(all_cols, patho_cols)
cat(length(patho_cols), "colonnes patho (JB_) |", length(wt_cols), "colonnes WT\n")

raw_stacked <- raw_stacked %>%
  mutate(across(all_of(all_cols), ~ suppressWarnings(as.numeric(.x))))

# ── 2. AGRÉGATION MAX PAR FUSION ─────────────────────────────────────────────
agg <- raw_stacked %>%
  group_by(seq_name) %>%
  summarise(across(all_of(all_cols), ~ suppressWarnings(max(.x, na.rm = TRUE))),
            n_fichiers = n_distinct(source_file), .groups = "drop") %>%
  mutate(across(all_of(all_cols), ~ ifelse(is.finite(.x), .x, 0)))
cat(nrow(agg), "fusions uniques après agrégation MAX\n")

# ── 3. SPÉCIFICITÉ TUMORALE (filtres sur le max par cohorte) ─────────────────
mat_patho <- as.matrix(agg[, patho_cols]); rownames(mat_patho) <- agg$seq_name
mat_wt    <- as.matrix(agg[, wt_cols])
n_patho   <- length(patho_cols); n_wt <- length(wt_cols)

spec <- agg %>%
  transmute(seq_name, n_fichiers) %>%
  mutate(
    max_patho   = apply(mat_patho, 1, max),
    max_wt      = apply(mat_wt,    1, max),
    # positivité par individu (pour la fréquence / les figures)
    n_patho_pos = rowSums(mat_patho >= opt$patho_min),
    n_wt_pos    = rowSums(mat_wt    >  opt$wt_min),
    freq_patho  = n_patho_pos / n_patho,
    freq_wt     = n_wt_pos / n_wt,
    samples_patho = apply(mat_patho >= opt$patho_min, 1,
                          function(r) paste(patho_cols[r], collapse = ", ")),
    # sélection : strictement absente des WT ET assez exprimée en patho
    specificite = case_when(
      max_wt <= opt$wt_min & max_patho >= opt$patho_min ~ "Chromo-spécifique",
      max_wt >  opt$wt_min                              ~ "Présente WT",
      TRUE                                              ~ "Sous seuil patho"
    )
  )

# ── 4. PARSE seq_name → gènes + breakpoints + clé paire de gènes ─────────────
mm <- str_match(spec$seq_name,
                "^(.*?)_([A-Za-z0-9.]+:[0-9]+)_(.*)_([A-Za-z0-9.]+:[0-9]+)$")
parsed <- spec %>%
  mutate(gene1 = mm[, 2], bp1 = str_remove(mm[, 3], "^chr"),
         gene2 = mm[, 4], bp2 = str_remove(mm[, 5], "^chr")) %>%
  mutate(g1n = resolve_alias(gene1), g2n = resolve_alias(gene2),
         gkey = map2_chr(g1n, g2n, ~ paste(sort(c(.x, .y)), collapse = "|")))
if (sum(is.na(parsed$bp1)) > 0)
  warning(sum(is.na(parsed$bp1)), " fusion(s) au seq_name non reconnu")

# ── 5. CHARGEMENT & ANNOTATION ARRIBA ────────────────────────────────────────
arriba_files <- list.files(opt$dir_arriba, pattern = "^JB_.*\\.tsv$", full.names = TRUE)
arriba_files <- arriba_files[!str_detect(arriba_files, "_discarded")]
if (length(arriba_files) == 0)
  stop("Aucun fichier Arriba JB_*.tsv dans : ", opt$dir_arriba)
cat(length(arriba_files), "fichier(s) Arriba trouvé(s)\n")

arriba <- map_dfr(arriba_files, function(f) {
  read_tsv(f, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
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
    gkey = map2_chr(resolve_alias(gene1), resolve_alias(gene2),
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
  group_by(gkey) %>%
  summarise(
    a_gene1 = first(gene1), a_gene2 = first(gene2),
    a_breakpoint1 = first(breakpoint1), a_breakpoint2 = first(breakpoint2),
    site1 = first(site1), site2 = first(site2),
    arriba_type = first(type), confidence = first(confidence),
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

# ── 6. JOINTURE (sur la paire de gènes) ──────────────────────────────────────
annot <- parsed %>%
  left_join(arriba_sum, by = "gkey") %>%
  mutate(arriba_matched = !is.na(confidence),
         type_chimerique = replace_na(type_chimerique, "Non confirmé Arriba"),
         is_who_interest = paste(g1n, g2n, sep = "--") %in% WHO_FUSIONS_INTEREST |
                           paste(g2n, g1n, sep = "--") %in% WHO_FUSIONS_INTEREST,
         is_who_interest = replace_na(is_who_interest, FALSE))

# ── 7. SCORE BIOLOGIQUE PONDÉRÉ & PARAMÉTRABLE ───────────────────────────────
# Chaque composante = fraction dans [0,1] × poids. Poids 0 => composante retirée.
MAX_SCORE <- with(opt, w_type + w_conf + w_spec + w_who + w_frame + w_reads)
if (MAX_SCORE <= 0) stop("Tous les poids sont nuls : score impossible.")

annot <- annot %>%
  mutate(
    frac_type = case_when(
      type_chimerique %in% c("Translocation", "Inversion") ~ 1.0,
      type_chimerique %in% c("Délétion", "Duplication")    ~ 0.5,
      type_chimerique == "Read-through"                    ~ 0.25,
      TRUE ~ 0),
    frac_conf = case_when(
      str_detect(coalesce(confidence, ""), regex("high",   ignore_case = TRUE)) ~ 1.0,
      str_detect(coalesce(confidence, ""), regex("medium", ignore_case = TRUE)) ~ 0.5,
      TRUE ~ 0),
    frac_spec = case_when(n_wt_pos == 0 ~ 1.0, n_wt_pos == 1 ~ 0.5, TRUE ~ 0),
    frac_who  = if_else(is_who_interest, 1.0, 0),
    frac_frame = case_when(reading_frame == "in-frame" ~ 1.0,
                           reading_frame == "out-of-frame" ~ 0, TRUE ~ 0.5),
    frac_reads = case_when(coalesce(total_reads, 0) >= 50 ~ 1.0,
                           coalesce(total_reads, 0) >= 10 ~ 0.5, TRUE ~ 0),
    score_type  = frac_type  * opt$w_type,
    score_conf  = frac_conf  * opt$w_conf,
    score_spec  = frac_spec  * opt$w_spec,
    score_who   = frac_who   * opt$w_who,
    score_frame = frac_frame * opt$w_frame,
    score_reads = frac_reads * opt$w_reads,
    score_total = score_type + score_conf + score_spec + score_who +
                  score_frame + score_reads,
    score_norm = score_total / MAX_SCORE,
    priorite = factor(case_when(
      score_norm >= 0.65 ~ "P1", score_norm >= 0.40 ~ "P2",
      score_norm >= 0.20 ~ "P3", TRUE ~ "NP"), levels = c("P1", "P2", "P3", "NP"))
  )

# ── 8. SORTIES TABLES ────────────────────────────────────────────────────────
out_cols <- c(
  "seq_name", "gene1", "bp1", "gene2", "bp2", "specificite", "is_who_interest",
  "score_norm", "priorite",
  "score_type", "score_conf", "score_spec", "score_who", "score_frame", "score_reads",
  "n_patho_pos", "freq_patho", "n_wt_pos", "freq_wt", "max_patho", "max_wt",
  "samples_patho", "n_fichiers",
  "arriba_matched", "type_chimerique", "class_ruffle", "reading_frame",
  "confidence", "arriba_type", "site1", "site2",
  "split_reads1", "split_reads2", "total_reads", "coverage1", "coverage2",
  "retained_protein_domains", "transcript_id1", "transcript_id2",
  "direction1", "direction2", "tags",
  "a_gene1", "a_gene2", "a_breakpoint1", "a_breakpoint2", "n_arriba", "samples_arriba")
annot_out <- annot %>% select(any_of(out_cols)) %>% arrange(desc(score_norm))

write_tsv(annot_out, file.path(DIR_OUT, "fusions_all_specificite_annotees.tsv"))
chromo_spec <- annot_out %>% filter(specificite == "Chromo-spécifique")
write_tsv(chromo_spec, file.path(DIR_OUT, "fusions_chromo_specifiques_annotees.tsv"))

# ── 9. FIGURES ───────────────────────────────────────────────────────────────
has_pheatmap <- requireNamespace("pheatmap",    quietly = TRUE)
has_repel    <- requireNamespace("ggrepel",     quietly = TRUE)
has_karyo    <- requireNamespace("karyoploteR", quietly = TRUE) &&
                requireNamespace("GenomicRanges", quietly = TRUE)
N_TOP <- opt$n_top

fig_df <- annot %>%
  filter(specificite == "Chromo-spécifique") %>%
  mutate(fusion_label = paste(gene1, gene2, sep = "--"),
         type_chimerique = factor(type_chimerique, levels = names(TYPE_COLORS)))

# Palette priorités (partagée)
PRIO_COLORS <- c(P1 = "#d62728", P2 = "#ff7f0e", P3 = "#9467bd", NP = "#bdbdbd")

# 9a. Barplot — top N par expression max (noms sur l'axe Y, toujours visibles)
p_bar <- fig_df %>% slice_max(max_patho, n = N_TOP, with_ties = FALSE) %>%
  mutate(fusion_ord = reorder(fusion_label, max_patho)) %>%
  ggplot(aes(max_patho, fusion_ord, fill = type_chimerique)) +
  geom_col() +
  scale_fill_manual(values = TYPE_COLORS, drop = FALSE, name = "Type chimérique") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = paste0("Top ", N_TOP, " fusions chromo-spécifiques — expression max"),
       x = "Comptage k-mer max chez un individu (patho)", y = NULL) +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(DIR_FIG, "barplot_expression.png"), p_bar, width = 11, height = 9, dpi = 150)

# 9b. Classement par score biologique (noms sur l'axe Y = toujours visibles)
p_rank <- fig_df %>% slice_max(score_norm, n = N_TOP, with_ties = FALSE) %>%
  mutate(fusion_ord = reorder(fusion_label, score_norm),
         etiq = paste0(class_ruffle, " · ", reading_frame)) %>%
  ggplot(aes(score_norm, fusion_ord, fill = priorite)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0.65, linetype = "dashed",  color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0.40, linetype = "dotted",  color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0.20, linetype = "dotdash", color = "grey60", linewidth = 0.3) +
  geom_text(aes(label = etiq), hjust = -0.05, size = 2.2, color = "grey25") +
  scale_fill_manual(values = PRIO_COLORS, drop = FALSE, name = "Priorité") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.18), breaks = seq(0, 1, 0.2)) +
  labs(title = paste0("Top ", N_TOP, " fusions par score biologique"),
       subtitle = paste0("Étiquette = classe Rufflé · cadre  |  score / ", MAX_SCORE, " pts"),
       x = "Score normalisé", y = NULL) +
  theme(axis.text.y = element_text(size = 7),
        panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_FIG, "score_classement.png"), p_rank, width = 12, height = 9, dpi = 200)

# 9c. Charge de fusions par échantillon (signature chromoanagenèse)
present <- mat_patho[fig_df$seq_name, , drop = FALSE] >= opt$patho_min
ftype   <- setNames(as.character(fig_df$type_chimerique), fig_df$seq_name)
burden  <- as.data.frame(present) %>%
  rownames_to_column("seq_name") %>%
  pivot_longer(-seq_name, names_to = "sample", values_to = "present") %>%
  filter(present) %>%
  mutate(type = factor(ftype[seq_name], levels = names(TYPE_COLORS))) %>%
  count(sample, type, name = "n")
if (nrow(burden) > 0) {
  tot <- burden %>% group_by(sample) %>% summarise(t = sum(n), .groups = "drop")
  p_burden <- burden %>%
    mutate(sample = factor(sample, levels = tot$sample[order(tot$t)])) %>%
    ggplot(aes(n, sample, fill = type)) +
    geom_col() +
    scale_fill_manual(values = TYPE_COLORS, drop = FALSE, name = "Type chimérique") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Charge de fusions chromo-spécifiques par échantillon",
         subtitle = "Nb de fusions portées par chaque patient, par type chimérique",
         x = "Nombre de fusions", y = NULL) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(DIR_FIG, "charge_par_echantillon.png"), p_burden,
         width = 10, height = 7, dpi = 150)
}

# 9d. Répartition des types chimériques / classes Rufflé
p_type <- fig_df %>% count(type_chimerique, class_ruffle, name = "n") %>%
  ggplot(aes(n, fct_reorder(type_chimerique, n, sum), fill = class_ruffle)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "Classe Rufflé", na.value = "grey70") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Répartition des types chimériques",
       x = "Nombre de fusions chromo-spécifiques", y = NULL) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(DIR_FIG, "repartition_types.png"), p_type, width = 9, height = 5, dpi = 150)

# 9e. Décomposition du score — seulement les composantes actives (poids > 0)
comp_def <- tibble::tribble(
  ~col,          ~label,       ~weight,
  "score_type",  "Type",       opt$w_type,
  "score_conf",  "Confiance",  opt$w_conf,
  "score_spec",  "Spéc. WT",   opt$w_spec,
  "score_who",   "WHO",        opt$w_who,
  "score_frame", "Cadre",      opt$w_frame,
  "score_reads", "Reads",      opt$w_reads
) %>% filter(weight > 0)
dec_df <- fig_df %>% slice_max(score_norm, n = N_TOP, with_ties = FALSE) %>%
  mutate(fusion_ord = reorder(fusion_label, score_norm)) %>%
  pivot_longer(all_of(comp_def$col), names_to = "comp", values_to = "val") %>%
  mutate(comp = factor(setNames(comp_def$label, comp_def$col)[comp],
                       levels = comp_def$label))
p_dec <- ggplot(dec_df, aes(val, fusion_ord, fill = comp)) +
  geom_col(width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(title = "Décomposition du score biologique",
       x = paste0("Points cumulés (max = ", MAX_SCORE, ")"), y = NULL) +
  theme(axis.text.y = element_text(size = 7), legend.position = "bottom")
ggsave(file.path(DIR_FIG, "score_decomposition.png"), p_dec, width = 11, height = 9, dpi = 150)

# 9d. Heatmap présence/absence (top score)
if (has_pheatmap) {
  bin <- (mat_patho >= opt$patho_min) * 1L
  top_seq <- fig_df %>% slice_max(score_norm, n = 40, with_ties = FALSE) %>% pull(seq_name)
  bin_top <- bin[rownames(bin) %in% top_seq, , drop = FALSE]
  if (nrow(bin_top) > 1) {
    lab_map <- setNames(fig_df$fusion_label, fig_df$seq_name)
    rownames(bin_top) <- make.unique(unname(lab_map[rownames(bin_top)]))
    bin_top <- bin_top[rowSums(bin_top) > 0, , drop = FALSE]
    ar <- fig_df %>% filter(seq_name %in% top_seq) %>%
      mutate(rn = make.unique(fusion_label)) %>%
      transmute(rn, Type = as.character(type_chimerique),
                Priorité = as.character(priorite)) %>%
      distinct(rn, .keep_all = TRUE) %>% column_to_rownames("rn")
    ar <- ar[rownames(bin_top), , drop = FALSE]
    pheatmap::pheatmap(
      bin_top, color = colorRampPalette(c("#f7f7f7", "#111827"))(2),
      annotation_row = ar, cluster_cols = TRUE, legend = FALSE,
      fontsize_row = 7, fontsize_col = 7,
      main = "Présence/absence — fusions chromo-spécifiques (top score)",
      filename = file.path(DIR_FIG, "heatmap_fusions.png"), width = 10, height = 11)
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
    png(file.path(DIR_FIG, "karyotype_overview.png"), width = 1300, height = 1000, res = 120)
    kp <- karyoploteR::plotKaryotype(genome = "hg38",
            main = "Breakpoints des fusions chromo-spécifiques (Arriba)")
    for (j in seq_len(nrow(ks))) {
      r <- ks[j, ]; col <- TYPE_COLORS[r$type_chimerique]
      if (is.na(col)) col <- "grey50"
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
cat("\nChromo-spécifiques :", nrow(chromo_spec),
    "| matchées Arriba :", sum(chromo_spec$arriba_matched), "\n")
cat("\nPriorités (chromo-spécifiques) :\n")
print(count(chromo_spec, priorite))
cat("\nScore /", MAX_SCORE, "pts | Figures :", DIR_FIG, "| Tables :", DIR_OUT, "\n")
if (!has_pheatmap) cat("(pheatmap absent → heatmap non générée)\n")
if (!has_karyo)    cat("(karyoploteR/GenomicRanges absents → karyotype non généré)\n")
