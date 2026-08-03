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
#   --type N   type chimérique   (défaut 3)
#   --conf N   confidence Arriba (défaut 2)
#   --spec N   spécificité WT    (défaut 1)
#   --who N    fusion WHO        (défaut 2)
#   --frame N  reading frame     (défaut 2)
#   --reads N  couverture reads  (défaut 5 ; > type : la couverture prime)
#
# Autres : --n-top N | --dir-merge | --dir-arriba | --dir-out | --help
#   --fig-format pdf|png   format des figures (défaut pdf)
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
  dir_fasta  = "fasta_JB",              # contigs de départ JB_*.fasta
  dir_kmers  = "JB_kmers",              # dossiers JB_*_kmers/kmers.fa
  dir_out    = "analyse_fusions",       # relatif au dossier courant (getwd())
  fig_format = "pdf",                    # format des figures : "pdf" ou "png"
  wt_min = 0, patho_min = 5,          # filtres (sur le max par cohorte)
  n_top = 30,                          # figures
  # Poids : la couverture en reads (w_reads) prime désormais sur le type
  # chimérique (w_type) — un signal de fusion bien couvert est plus fiable
  # qu'un type « fort » faiblement supporté.
  w_type = 3, w_conf = 2, w_spec = 1, w_who = 2, w_frame = 2, w_reads = 5
)

# Seuil de distance read-through (bp) — délétion courte entre gènes voisins
# colinéaires (Rufflé 2024). Sert au repli géométrique du type chimérique.
READTHROUGH_DIST <- 300000L

# ── PARSER CLI ───────────────────────────────────────────────────────────────
alias <- c(type = "w_type", conf = "w_conf", confidence = "w_conf",
           spec = "w_spec", specificity = "w_spec", who = "w_who",
           frame = "w_frame", reads = "w_reads",
           fig = "fig_format", format = "fig_format")
string_opts <- c("dir_merge", "dir_arriba", "dir_fasta", "dir_kmers", "dir_out",
                 "fig_format")

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || "-h" %in% args) {
  cat("Options : --wt-min --patho-min --type --conf --spec --who --frame",
      "--reads --n-top --dir-merge --dir-arriba --dir-out --fig-format\n")
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

# ── Format des figures (pdf par défaut) ──────────────────────────────────────
FIG_EXT <- tolower(opt$fig_format)
if (!FIG_EXT %in% c("pdf", "png")) {
  warning("Format figure inconnu (", opt$fig_format, ") — pdf utilisé.")
  FIG_EXT <- "pdf"
}
# Chemin d'une figure avec l'extension choisie
fig_path <- function(name) file.path(DIR_FIG, paste0(name, ".", FIG_EXT))

cat("=== CONFIGURATION ===\n")
cat("Dossier courant :", getwd(), "\n")
cat("  merge  :", normalizePath(opt$dir_merge,  mustWork = FALSE), "\n")
cat("  arriba :", normalizePath(opt$dir_arriba, mustWork = FALSE), "\n")
cat("  fasta  :", normalizePath(opt$dir_fasta,  mustWork = FALSE), "\n")
cat("  kmers  :", normalizePath(opt$dir_kmers,  mustWork = FALSE), "\n")
cat("  sortie :", normalizePath(opt$dir_out,    mustWork = FALSE), "\n")
cat("  figures:", DIR_FIG, "(format", FIG_EXT, ")\n")
cat(sprintf("Filtres : max(WT) <= %s  |  max(patho) >= %s\n", opt$wt_min, opt$patho_min))
cat(sprintf("Poids   : type=%s conf=%s spec=%s who=%s frame=%s reads=%s\n\n",
            opt$w_type, opt$w_conf, opt$w_spec, opt$w_who, opt$w_frame, opt$w_reads))

# ── Fusions WHO d'intérêt (composante E) + alias HGNC ─────────────────────────
WHO_FUSIONS_INTEREST <- c(
  "RUNX1--RUNX1T1", "CBFB--MYH11", "PML--RARA", "KMT2A--MLLT3",
  "DEK--NUP214", "NUP98--NSD1", "BCR--ABL1", "NRIP1--MIR99AHG"
)
# Gènes dont une duplication interne en tandem (ITD) est un driver LAM connu.
# Une ITD (ex. FLT3--FLT3) est traitée comme un driver d'intérêt : elle reçoit
# le bonus WHO et une fraction de type pleine, pour la faire remonter en P1.
DRIVER_ITD_GENES <- c("FLT3", "KMT2A", "UBTF")
GENE_ALIASES <- c(
  MLL = "KMT2A", HRX = "KMT2A", ALL1 = "KMT2A", MLL1 = "KMT2A", TRX1 = "KMT2A",
  AF9 = "MLLT3", LTG9 = "MLLT3", AF6 = "MLLT4", AFDN = "MLLT4", ENL = "MLLT1",
  ETO = "RUNX1T1", MTG8 = "RUNX1T1", CBFA2T1 = "RUNX1T1", BCR1 = "BCR"
)
resolve_alias <- function(g) ifelse(g %in% names(GENE_ALIASES), GENE_ALIASES[g], g)

# ── Lecture FASTA simple : renvoie tibble(id, seq) ───────────────────────────
read_fasta <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[!grepl("^\\s*$", ln)]
  if (length(ln) == 0) return(tibble(id = character(), seq = character()))
  is_h <- startsWith(ln, ">")
  grp  <- cumsum(is_h)
  ids  <- sub("^>", "", ln[is_h])
  seqs_by_grp <- tapply(ln[!is_h], grp[!is_h], paste0, collapse = "")
  seq_vec <- rep("", length(ids))
  seq_vec[as.integer(names(seqs_by_grp))] <- as.character(seqs_by_grp)
  tibble(id = ids, seq = seq_vec)
}

# ── Type & classe à partir du type Arriba détaillé ───────────────────────────
# arriba_type ex : "translocation", "deletion/read-through", "duplication/ITD"...
# On prend le token de base (avant "/") ; read-through détecté à part.
arriba_base <- function(t) {
  tl <- str_to_lower(coalesce(t, ""))
  b  <- str_split_fixed(tl, "/", 2)[, 1]
  dplyr::case_when(
    str_detect(tl, "read-through") ~ "Read-through",
    b == "translocation"           ~ "Translocation",
    b == "inversion"               ~ "Inversion",
    b == "duplication"             ~ "Duplication",
    b == "deletion"                ~ "Délétion",
    is.na(t) | t == ""             ~ "Non confirmé Arriba",
    TRUE                           ~ "Autre")
}
ruffle_class <- function(t) {
  b <- str_split_fixed(str_to_lower(coalesce(t, "")), "/", 2)[, 1]
  dplyr::case_when(
    b == "translocation" ~ "Class 1",
    b == "inversion"     ~ "Class 4",
    b == "duplication"   ~ "Class 3",
    b == "deletion"      ~ "Class 2",
    TRUE                 ~ NA_character_)
}

# ── Constantes visuelles ─────────────────────────────────────────────────────
TYPE_COLORS <- c(
  Translocation = "#d62728", Inversion = "#9467bd", "Délétion" = "#ff7f0e",
  Duplication = "#1f77b4", "Read-through" = "#2ca02c",
  "Indéterminé (même chr)" = "grey75",
  "Non confirmé Arriba" = "grey65", Autre = "grey50"
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
    sum_patho   = rowSums(mat_patho),        # comptage cumulé sur la cohorte patho
    sumsq_patho = rowSums(mat_patho^2),      # pour l'indice de focalité
    # n_eff : nombre EFFECTIF de patients porteurs (~1 = mono-patient, ~k = diffus)
    n_eff_patho = if_else(sumsq_patho > 0, sum_patho^2 / sumsq_patho, 0),
    # foc_index : max² / somme — élevé si forte expression concentrée sur peu de
    # patients, faible si diffuse et faible. Chiffre de focalité (spécificité patient).
    foc_index   = if_else(sum_patho > 0, max_patho^2 / sum_patho, 0),
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
    gkey = map2_chr(resolve_alias(gene1), resolve_alias(gene2),
                    ~ paste(sort(c(.x, .y)), collapse = "|")),
    total_reads = suppressWarnings(as.numeric(split_reads1) + as.numeric(split_reads2)),
    # couverture de jonction = split reads + paires discordantes (support total)
    support_reads = suppressWarnings(
      as.numeric(split_reads1) + as.numeric(split_reads2) +
      coalesce(as.numeric(discordant_mates), 0)),
    conf_num = case_when(
      str_detect(confidence, regex("high",   ignore_case = TRUE)) ~ 3L,
      str_detect(confidence, regex("medium", ignore_case = TRUE)) ~ 2L,
      str_detect(confidence, regex("low",    ignore_case = TRUE)) ~ 1L,
      TRUE ~ 0L))

# Résumé Arriba : une ligne par paire de gènes, on garde le meilleur record.
# On n'affiche PAS les gènes/breakpoints d'Arriba : on garde ceux du seq_name.
arriba_sum <- arriba %>%
  arrange(desc(conf_num)) %>%
  group_by(gkey) %>%
  summarise(
    site1 = first(site1), site2 = first(site2),
    arriba_type = first(type), confidence = first(confidence),
    reading_frame = case_when(
      any(reading_frame == "in-frame",     na.rm = TRUE) ~ "in-frame",
      any(reading_frame == "out-of-frame", na.rm = TRUE) ~ "out-of-frame",
      TRUE ~ "."),
    split_reads1 = first(split_reads1), split_reads2 = first(split_reads2),
    total_reads  = suppressWarnings(max(total_reads, na.rm = TRUE)),
    support_reads = suppressWarnings(max(support_reads, na.rm = TRUE)),
    coverage1 = first(coverage1), coverage2 = first(coverage2),
    transcript_id1 = first(transcript_id1), transcript_id2 = first(transcript_id2),
    direction1 = first(direction1), direction2 = first(direction2),
    n_arriba = n_distinct(sample),
    samples_arriba = paste(sort(unique(sample)), collapse = ", "),
    .groups = "drop") %>%
  mutate(total_reads   = ifelse(is.finite(total_reads),   total_reads,   NA_real_),
         support_reads = ifelse(is.finite(support_reads), support_reads, NA_real_))

# ── 6. JOINTURE (sur la paire de gènes) ──────────────────────────────────────
# type_base (couleurs) et class_ruffle dérivés du type Arriba détaillé.
# Cascade de détermination du type chimérique :
#   (1) colonne `type` d'Arriba si renseignée (source de vérité) ;
#   (2) repli géométrique : directions Arriba + ordre des breakpoints
#       (reproduit get_fusion_type d'Arriba, ~99,6 % de concordance) ;
#   (3) repli minimal : chromosomes seuls (translocation vs même chr) ;
#   (4) sinon "Non confirmé Arriba".
annot <- parsed %>%
  left_join(arriba_sum, by = "gkey") %>%
  mutate(
    arriba_matched = !is.na(confidence),
    # coordonnées issues du seq_name (toujours présentes) pour le repli
    tmp_chr1 = str_split_fixed(bp1, ":", 2)[, 1],
    tmp_pos1 = suppressWarnings(as.numeric(str_split_fixed(bp1, ":", 2)[, 2])),
    tmp_chr2 = str_split_fixed(bp2, ":", 2)[, 1],
    tmp_pos2 = suppressWarnings(as.numeric(str_split_fixed(bp2, ":", 2)[, 2])),
    tmp_d1   = str_trim(coalesce(direction1, "")),
    tmp_d2   = str_trim(coalesce(direction2, "")),
    tmp_has_type = !is.na(arriba_type) & !arriba_type %in% c("", "."),
    tmp_has_dir  = tmp_d1 %in% c("upstream", "downstream") &
                   tmp_d2 %in% c("upstream", "downstream"),
    tmp_has_bp   = !is.na(tmp_pos1) & !is.na(tmp_pos2) &
                   tmp_chr1 != "" & tmp_chr2 != "",
    # orientation délétion vs duplication (logique get_fusion_type d'Arriba)
    tmp_del_orient = (tmp_d1 == "downstream" & tmp_pos1 < tmp_pos2) |
                     (tmp_d1 == "upstream"   & tmp_pos1 > tmp_pos2),
    # (2) type reconstruit par géométrie (directions + breakpoints)
    tmp_type_geom = case_when(
      !tmp_has_dir | !tmp_has_bp                                       ~ NA_character_,
      tmp_chr1 != tmp_chr2                                             ~ "Translocation",
      tmp_d1 == tmp_d2                                                 ~ "Inversion",
      tmp_del_orient & abs(tmp_pos2 - tmp_pos1) < READTHROUGH_DIST     ~ "Read-through",
      tmp_del_orient                                                   ~ "Délétion",
      TRUE                                                            ~ "Duplication"),
    # (3) repli minimal sur les chromosomes seuls
    tmp_type_chr = case_when(
      !tmp_has_bp          ~ NA_character_,
      tmp_chr1 != tmp_chr2 ~ "Translocation",
      TRUE                 ~ "Indéterminé (même chr)"),
    # type final + traçabilité de la source
    type_base = case_when(
      tmp_has_type              ~ arriba_base(arriba_type),
      !is.na(tmp_type_geom)     ~ tmp_type_geom,
      !is.na(tmp_type_chr)      ~ tmp_type_chr,
      TRUE                      ~ "Non confirmé Arriba"),
    type_source = case_when(
      tmp_has_type          ~ "arriba",
      !is.na(tmp_type_geom) ~ "géométrie",
      !is.na(tmp_type_chr)  ~ "chr",
      TRUE                  ~ "absent"),
    # classe Rufflé dérivée du type final (couvre aussi les replis géométriques)
    class_ruffle = case_when(
      type_base == "Translocation" ~ "Class 1",
      type_base == "Délétion"      ~ "Class 2",
      type_base == "Read-through"  ~ "Class 2",
      type_base == "Duplication"   ~ "Class 3",
      type_base == "Inversion"     ~ "Class 4",
      TRUE                         ~ NA_character_),
    is_who_interest = paste(g1n, g2n, sep = "--") %in% WHO_FUSIONS_INTEREST |
                      paste(g2n, g1n, sep = "--") %in% WHO_FUSIONS_INTEREST,
    is_who_interest = replace_na(is_who_interest, FALSE),
    # ITD driver (ex. FLT3-ITD) : duplication interne d'un même gène driver
    is_itd = str_detect(coalesce(arriba_type, ""), regex("ITD", ignore_case = TRUE)),
    is_driver_itd = replace_na(
      is_itd & g1n == g2n & g1n %in% DRIVER_ITD_GENES, FALSE),
    # driver d'intérêt = fusion WHO connue OU ITD driver
    is_driver_interest = is_who_interest | is_driver_itd) %>%
  select(-starts_with("tmp_"))

# ── 7. SCORE BIOLOGIQUE PONDÉRÉ & PARAMÉTRABLE ───────────────────────────────
# Chaque composante = fraction dans [0,1] × poids. Poids 0 => composante retirée.
MAX_SCORE <- with(opt, w_type + w_conf + w_spec + w_who + w_frame + w_reads)
if (MAX_SCORE <= 0) stop("Tous les poids sont nuls : score impossible.")

annot <- annot %>%
  mutate(
    frac_type = case_when(
      is_driver_itd                                  ~ 1.0,   # FLT3-ITD & co : signal fort
      type_base %in% c("Translocation", "Inversion") ~ 1.0,
      type_base %in% c("Délétion", "Duplication")    ~ 0.5,
      type_base == "Read-through"                    ~ 0.25,
      TRUE ~ 0),
    frac_conf = case_when(
      str_detect(coalesce(confidence, ""), regex("high",   ignore_case = TRUE)) ~ 1.0,
      str_detect(coalesce(confidence, ""), regex("medium", ignore_case = TRUE)) ~ 0.5,
      TRUE ~ 0),
    frac_spec = case_when(n_wt_pos == 0 ~ 1.0, n_wt_pos == 1 ~ 0.5, TRUE ~ 0),
    # bonus driver connu : fusions WHO d'intérêt ET ITD drivers (FLT3-ITD…)
    frac_who  = if_else(is_driver_interest, 1.0, 0),
    frac_frame = case_when(reading_frame == "in-frame" ~ 1.0,
                           reading_frame == "out-of-frame" ~ 0, TRUE ~ 0.5),
    # couverture reads : granularité fine sur le support de jonction
    # (split reads + paires discordantes). Fraction dans [0,1] × w_reads.
    frac_reads = case_when(
      coalesce(support_reads, 0) >= 100 ~ 1.00,
      coalesce(support_reads, 0) >=  50 ~ 0.80,
      coalesce(support_reads, 0) >=  25 ~ 0.60,
      coalesce(support_reads, 0) >=  10 ~ 0.40,
      coalesce(support_reads, 0) >=   5 ~ 0.20,
      TRUE                              ~ 0),
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

# ── 7b. CONTIGS (JB_*.fasta) & K-MERS (JB_*_kmers/kmers.fa) ───────────────────
# Header du contig = seq_name ; headers k-mers = "<seq_name>.kmerN".
# Un même seq_name peut avoir des contigs DIFFÉRENTS selon l'échantillon
# (fusion_transcript reconstruit à partir des reads). Chaque contig distinct
# devient une entrée : seq_name, puis seq_name.2, seq_name.3 pour les variantes.
fasta_files <- list.files(opt$dir_fasta, pattern = "^JB_.*\\.fasta$", full.names = TRUE)
kmer_files  <- list.files(opt$dir_kmers, pattern = "^kmers\\.fa$",
                          recursive = TRUE, full.names = TRUE)

# Contigs par (échantillon, seq_name)
if (length(fasta_files) > 0) {
  contigs_s <- map_dfr(fasta_files, ~ read_fasta(.x) %>%
                         mutate(sample = str_extract(basename(.x), "JB_[0-9]+"))) %>%
    filter(seq != "") %>% transmute(sample, seq_name = id, contig_seq = seq)
} else {
  contigs_s <- tibble(sample = character(), seq_name = character(), contig_seq = character())
}

# K-mers par (échantillon, seq_name), une colonne par position
if (length(kmer_files) > 0) {
  kmers_s <- map_dfr(kmer_files, ~ read_fasta(.x) %>%
                       mutate(sample = str_extract(.x, "JB_[0-9]+"))) %>%
    filter(seq != "") %>%
    mutate(seq_name = sub("\\.kmer.*$", "", id),
           kidx = as.integer(str_match(id, "\\.kmer([0-9]+)")[, 2])) %>%
    filter(!is.na(kidx))
} else {
  kmers_s <- tibble(sample = character(), seq_name = character(),
                    kidx = integer(), seq = character())
}

KMER_COLS <- if (nrow(kmers_s) > 0) paste0("kmer", seq_len(max(kmers_s$kidx))) else character(0)
if (nrow(kmers_s) > 0) {
  kmers_sw <- kmers_s %>% distinct(sample, seq_name, kidx, .keep_all = TRUE) %>%
    mutate(kcol = paste0("kmer", kidx)) %>%
    pivot_wider(id_cols = c(sample, seq_name), names_from = kcol, values_from = seq)
  nk <- kmers_s %>% distinct(sample, seq_name, kidx) %>% count(sample, seq_name, name = "n_kmers")
} else {
  kmers_sw <- tibble(sample = character(), seq_name = character())
  nk <- tibble(sample = character(), seq_name = character(), n_kmers = integer())
}

# Une entrée par contig distinct ; suffixe .2, .3 pour les variantes d'un seq_name
variants <- contigs_s %>%
  left_join(kmers_sw, by = c("sample", "seq_name")) %>%
  left_join(nk,       by = c("sample", "seq_name")) %>%
  arrange(seq_name, sample) %>%
  distinct(seq_name, contig_seq, .keep_all = TRUE) %>%   # contigs identiques fusionnés
  group_by(seq_name) %>% mutate(v = row_number()) %>% ungroup() %>%
  mutate(fusion_id = if_else(v == 1, seq_name, paste0(seq_name, ".", v)),
         n_kmers   = coalesce(n_kmers, 0L)) %>%
  select(fusion_id, seq_name, contig_seq, n_kmers, any_of(KMER_COLS))
cat(nrow(variants), "entrées contig (dont",
    sum(str_detect(variants$fusion_id, "\\.[0-9]+$")), "variantes .N)\n")

# ── 8. SORTIES TABLES ────────────────────────────────────────────────────────
base_cols <- c(
  "seq_name", "gene1", "bp1", "gene2", "bp2", "specificite",
  "is_who_interest", "is_driver_itd",
  "score_norm", "priorite",
  "score_type", "score_conf", "score_spec", "score_who", "score_frame", "score_reads",
  "n_patho_pos", "freq_patho", "n_wt_pos", "freq_wt",
  "max_patho", "sum_patho", "n_eff_patho", "foc_index", "max_wt",
  "samples_patho", "n_fichiers",
  "arriba_matched", "type_base", "type_source", "arriba_type", "class_ruffle",
  "reading_frame", "confidence", "site1", "site2",
  "split_reads1", "split_reads2", "total_reads", "support_reads",
  "coverage1", "coverage2",
  "transcript_id1", "transcript_id2", "direction1", "direction2",
  "n_arriba", "samples_arriba")
# Jointure 1:N (duplique les fusions à plusieurs contigs) — sorties uniquement
annot_out <- annot %>% select(any_of(base_cols)) %>%
  left_join(variants, by = "seq_name") %>%
  mutate(fusion_id = coalesce(fusion_id, seq_name)) %>%
  relocate(fusion_id) %>%
  select(any_of(c("fusion_id", base_cols, "contig_seq", "n_kmers", KMER_COLS))) %>%
  arrange(desc(score_norm), fusion_id)

write_tsv(annot_out, file.path(DIR_OUT, "fusions_all_specificite_annotees.tsv"))
chromo_spec <- annot_out %>% filter(specificite == "Chromo-spécifique")
write_tsv(chromo_spec, file.path(DIR_OUT, "fusions_chromo_specifiques_annotees.tsv"))

# ── 9. FIGURES ───────────────────────────────────────────────────────────────
has_pheatmap <- requireNamespace("pheatmap",    quietly = TRUE)
has_repel    <- requireNamespace("ggrepel",     quietly = TRUE)
has_karyo    <- requireNamespace("karyoploteR", quietly = TRUE) &&
                requireNamespace("GenomicRanges", quietly = TRUE)
N_TOP <- opt$n_top

# Représentations : uniquement les fusions chromo-spécifiques ANNOTÉES par Arriba
# (arriba_matched). Pour inclure aussi les chromo-spé. non retrouvées par Arriba,
# retirer le second filtre ci-dessous.
fig_df <- annot %>%
  filter(specificite == "Chromo-spécifique", arriba_matched) %>%
  mutate(fusion_label = paste(gene1, gene2, sep = "--"),
         type_base = factor(type_base, levels = names(TYPE_COLORS)))
cat(nrow(fig_df), "fusions chromo-spécifiques annotées Arriba (base des figures)\n")

# Une même paire de gènes peut avoir plusieurs breakpoints (plusieurs seq_name).
# Les figures "par fusion" (barplot, classement, décomposition, répartition)
# doivent afficher UNE barre par paire — sinon geom_col empile les variantes et
# gonfle les totaux (score cumulé > MAX_SCORE, classement > 100 %).
# fig_uni : un représentant par paire = meilleur breakpoint (score le plus élevé),
# avec l'expression max prise sur l'ensemble des variantes.
fig_uni <- fig_df %>%
  group_by(fusion_label) %>%
  mutate(max_patho = max(max_patho, na.rm = TRUE),
         sum_patho = max(sum_patho, na.rm = TRUE)) %>%
  slice_max(score_norm, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(foc_cat = factor(case_when(
           n_patho_pos <= 1 ~ "1 patient (mono)",
           n_patho_pos <= 3 ~ "2–3 patients (focale)",
           TRUE             ~ "≥ 4 patients (diffuse)"),
         levels = c("1 patient (mono)", "2–3 patients (focale)", "≥ 4 patients (diffuse)")))
cat(nrow(fig_uni), "paires de gènes uniques (base des figures par fusion)\n")

# Palette priorités (partagée)
PRIO_COLORS <- c(P1 = "#d62728", P2 = "#ff7f0e", P3 = "#9467bd", NP = "#bdbdbd")

# 9a. Barplot — top N par expression : total cohorte patho + max chez un individu
# Deux barres par fusion (dodge) : comptage k-mer cumulé sur les patients patho
# (rouge) et maximum observé chez un seul individu (orange). L'écart entre les
# deux reflète l'hétérogénéité intra-cohorte (fusion diffuse vs portée par 1-2 cas).
p_bar <- fig_uni %>% slice_max(sum_patho, n = N_TOP, with_ties = FALSE) %>%
  mutate(fusion_ord = reorder(fusion_label, sum_patho)) %>%
  pivot_longer(c(sum_patho, max_patho), names_to = "metric", values_to = "valeur") %>%
  mutate(metric = c(sum_patho = "Total cohorte patho",
                    max_patho = "Max chez un individu")[metric],
         metric = factor(metric, levels = c("Total cohorte patho", "Max chez un individu"))) %>%
  ggplot(aes(valeur, fusion_ord, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("Total cohorte patho"  = "#d62728",
                               "Max chez un individu" = "#ff7f0e"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = scales::comma_format(big.mark = " ")) +
  labs(title = paste0("Top ", N_TOP, " fusions chromo-spécifiques — expression k-mers"),
       subtitle = "Total cumulé sur la cohorte patho vs maximum chez un seul individu",
       x = "Comptage k-mer", y = NULL) +
  theme(axis.text.y = element_text(size = 7), legend.position = "top")
ggsave(fig_path("barplot_expression"), p_bar, width = 11, height = 9, dpi = 150)

# 9b. Classement par score biologique (noms sur l'axe Y = toujours visibles)
p_rank <- fig_uni %>% slice_max(score_norm, n = N_TOP, with_ties = FALSE) %>%
  mutate(fusion_ord = reorder(fusion_label, score_norm),
         etiq = paste0(coalesce(class_ruffle, "—"), " · ", coalesce(reading_frame, "—"))) %>%
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
ggsave(fig_path("score_classement"), p_rank, width = 12, height = 9, dpi = 200)

# 9c. Charge de fusions par échantillon (signature chromoanagenèse)
# On ne compte que les fusions au type confirmé par Arriba ;
# les non confirmées / inconnues sont signalées en sous-titre.
fig_conf   <- fig_df %>% filter(!type_base %in% c("Non confirmé Arriba", "Autre"))
n_non_conf <- nrow(fig_df) - nrow(fig_conf)
sub_burden <- paste0("Nb de fusions portées par chaque patient, par type (Arriba)",
  if (n_non_conf > 0) paste0("  |  ", n_non_conf,
      " fusion(s) chromo-spé. non retrouvée(s) dans Arriba (exclues)") else "")
present <- mat_patho[fig_conf$seq_name, , drop = FALSE] >= opt$patho_min
ftype   <- setNames(as.character(fig_conf$type_base), fig_conf$seq_name)
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
    scale_fill_manual(values = TYPE_COLORS, drop = TRUE, name = "Type (Arriba)") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Charge de fusions chromo-spécifiques par échantillon",
         subtitle = sub_burden, x = "Nombre de fusions", y = NULL) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(fig_path("charge_par_echantillon"), p_burden,
         width = 10, height = 7, dpi = 150)
}

# 9d. Répartition des types chimériques / classes Rufflé
p_type <- fig_uni %>% count(type_base, class_ruffle, name = "n") %>%
  ggplot(aes(n, fct_reorder(type_base, n, sum), fill = class_ruffle)) +
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "Classe Rufflé", na.value = "grey70") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Répartition des types chimériques",
       x = "Nombre de fusions chromo-spécifiques", y = NULL) +
  theme(plot.title = element_text(face = "bold"))
ggsave(fig_path("repartition_types"), p_type, width = 9, height = 5, dpi = 150)

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
dec_df <- fig_uni %>% slice_max(score_norm, n = N_TOP, with_ties = FALSE) %>%
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
ggsave(fig_path("score_decomposition"), p_dec, width = 11, height = 9, dpi = 150)

# 9f. Carte de priorisation : score biologique × expression focale ────────────
# Met en évidence les fusions potentiellement importantes :
#   x = score biologique  |  y = max de comptage chez UN patient (log)
#   couleur = focalité (exprimée chez 1 / 2-3 / >=4 patients)
#   taille  = foc_index (max^2 / somme : force × concentration).
# Cadran cible = haut-droite (score élevé) + points chauds (mono / focale).
foc_df <- fig_uni %>%
  mutate(a_labeliser = priorite == "P1")   # noms affichés seulement pour les P1
FOC_COLORS <- c("1 patient (mono)"       = "#d62728",
                "2–3 patients (focale)"  = "#ff7f0e",
                "≥ 4 patients (diffuse)" = "#1f77b4")
p_focal <- ggplot(foc_df, aes(score_norm, max_patho)) +
  annotate("rect", xmin = 0.65, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#d62728", alpha = 0.04) +
  geom_vline(xintercept = 0.65, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0.40, linetype = "dotted", color = "grey55", linewidth = 0.3) +
  geom_point(aes(color = foc_cat, size = foc_index), alpha = 0.8) +
  scale_y_log10(labels = scales::comma_format(big.mark = " ")) +
  scale_color_manual(values = FOC_COLORS, name = "Focalité") +
  scale_size_continuous(name = "Indice focalité\n(max² / somme)", range = c(1.5, 9)) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Carte de priorisation des fusions chromo-spécifiques",
       subtitle = "Score biologique × expression max par patient · focalité = spécificité à un sous-groupe",
       x = "Score biologique", y = "Comptage k-mer max chez un patient (log)") +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")
# Étiquettes des candidates : ggrepel si dispo (rien n'est masqué -> Inf),
# sinon repli sur geom_text pour que les noms s'affichent toujours.
lab_df <- dplyr::filter(foc_df, a_labeliser)
cat(nrow(lab_df), "fusion(s) étiquetée(s) sur la carte de priorisation",
    if (!has_repel) "(geom_text — installer ggrepel pour un placement propre)" else "", "\n")
if (nrow(lab_df) > 0) {
  if (has_repel) {
    p_focal <- p_focal +
      ggrepel::geom_text_repel(
        data = lab_df, aes(label = fusion_label), size = 2.6,
        max.overlaps = Inf, min.segment.length = 0, box.padding = 0.4,
        segment.color = "grey60", color = "grey20")
  } else {
    p_focal <- p_focal +
      geom_text(data = lab_df, aes(label = fusion_label),
                size = 2.4, vjust = -0.8, color = "grey20", check_overlap = TRUE)
  }
}
ggsave(fig_path("carte_priorisation"), p_focal, width = 11, height = 8, dpi = 200)

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
      transmute(rn, Type = as.character(type_base),
                Priorité = as.character(priorite)) %>%
      distinct(rn, .keep_all = TRUE) %>% column_to_rownames("rn")
    ar <- ar[rownames(bin_top), , drop = FALSE]
    pheatmap::pheatmap(
      bin_top, color = colorRampPalette(c("#f7f7f7", "#111827"))(2),
      annotation_row = ar, cluster_cols = TRUE, legend = FALSE,
      fontsize_row = 7, fontsize_col = 7,
      main = "Présence/absence — fusions chromo-spécifiques (top score)",
      filename = fig_path("heatmap_fusions"), width = 10, height = 11)
  }
}

# 9e. Karyotypes (breakpoints des fusions confirmées Arriba) — ensemble + par patient
draw_karyo <- function(ks_df, titre, fichier) {
  if (nrow(ks_df) == 0) return(invisible())
  if (FIG_EXT == "pdf") pdf(fichier, width = 11, height = 8.5)
  else                  png(fichier, width = 1300, height = 1000, res = 120)
  kp <- karyoploteR::plotKaryotype(genome = "hg38", main = titre)
  for (j in seq_len(nrow(ks_df))) {
    r <- ks_df[j, ]; col <- TYPE_COLORS[r$type_base]
    if (is.na(col)) col <- "grey50"
    try(karyoploteR::kpPlotLinks(
      kp,
      data  = GenomicRanges::GRanges(r$chr1, IRanges::IRanges(r$pos1, width = 1)),
      data2 = GenomicRanges::GRanges(r$chr2, IRanges::IRanges(r$pos2, width = 1)),
      col   = scales::alpha(col, 0.6)), silent = TRUE)
  }
  dev.off()
}

if (has_karyo) {
  ks_all <- fig_df %>% filter(arriba_matched) %>%
    transmute(seq_name,
      chr1 = paste0("chr", str_remove(str_split_fixed(bp1, ":", 2)[, 1], "^chr")),
      pos1 = suppressWarnings(as.numeric(str_split_fixed(bp1, ":", 2)[, 2])),
      chr2 = paste0("chr", str_remove(str_split_fixed(bp2, ":", 2)[, 1], "^chr")),
      pos2 = suppressWarnings(as.numeric(str_split_fixed(bp2, ":", 2)[, 2])),
      type_base = as.character(type_base)) %>%
    filter(!is.na(pos1), !is.na(pos2))

  # Vue d'ensemble (toutes cohortes confondues)
  draw_karyo(ks_all, "Breakpoints des fusions chromo-spécifiques (toutes)",
             fig_path("karyotype_overview"))

  # Un karyotype par patient : fusions présentes chez cet échantillon
  dir_kar <- file.path(DIR_FIG, "karyotypes")
  dir.create(dir_kar, showWarnings = FALSE, recursive = TRUE)
  if (nrow(ks_all) > 0) {
    present_k <- mat_patho[ks_all$seq_name, , drop = FALSE] >= opt$patho_min
    for (s in colnames(present_k)) {
      seqs  <- rownames(present_k)[present_k[, s]]
      ks_s  <- ks_all %>% filter(seq_name %in% seqs)
      if (nrow(ks_s) == 0) next
      draw_karyo(ks_s, paste0(s, " — ", nrow(ks_s), " fusion(s) chromo-spé."),
                 file.path(dir_kar, paste0("karyotype_", s, ".", FIG_EXT)))
    }
    cat("Karyotypes par patient écrits dans :", dir_kar, "\n")
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
