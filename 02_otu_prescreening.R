# =============================================================================
# Comparative Performance of Traditional Sampling and Water Environmental DNA
# for Detecting Fish Diversity in Amazonian Streams Across Hydrological Periods
#
# Authors: Fabricio dos Anjos Santa Rosa, Derlan José Ferreira Silva,
#          Cintia Oliveira Carvalho, Silvia Britto Barreto,
#          Silvana Melo Sviggum, Anna Karolina Oliveira de Queiroz,
#          Enki Gautier, Marcia da Silva Anjos,
#          Jessica Conceição Nascimento Dergan,
#          Orlando Gabriel da Silva Solano Reis,
#          Tibério Cesar Tortola Burlamaqui, Birgitte Lisbeth Graae Thorbek,
#          João Bráullio de Luna Sales, Quentin Mauvisseau, Hugo de Boer,
#          Jonathan Stuart Ready
# Contact: Fabricio dos Anjos Santa Rosa (fabriciodosanjos14@gmail.com)
#
# Script:  02_otu_prescreening.R
# Purpose: Automated pre-screening of fish OTUs (step 2): stop-codon screening,
#          abundance and sequence-divergence metrics, BOLD/GenBank validation
#          and a proposed curation decision per OTU. Decisions were reviewed
#          manually; final decisions are in sheet 2_reads_curation_decisions,
#          where 24 of 137 OTUs differ from this pre-screening.
# Input:   SuppInfo_S1.xlsx, sheets 2_OTU_sequences and 2_reads_after_controls
#          02_otu_prescreening_functions.R (BOLD/GenBank helpers)
# Output:  global_curation_output_final_water_per_subsample_upd.txt
#          presence_absence_matrix_curated_water_per_subsample_upd.txt
# Notes:   BOLD/GenBank validation requires internet access and a Chrome
#          browser (chromote).
# =============================================================================

library(Biostrings)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(readr)

# Part 1. Stop-codon screening ================================================
# 1.1 Inputs ------------------------------------------------------------------

# Working directory: folder of this script
args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
if (length(script_dir) == 1 && nzchar(script_dir)) setwd(script_dir)

workbook <- "SuppInfo_S1.xlsx"

# 1.2 Read OTU sequences ------------------------------------------------------

seq_tbl  <- readxl::read_excel(workbook, sheet = "2_OTU_sequences")
seqs_raw <- Biostrings::DNAStringSet(setNames(seq_tbl$sequence, seq_tbl$seq_id))

raw_seq_df <- tibble::tibble(
  seq_id = names(seqs_raw),
  raw_sequence = as.character(seqs_raw),
  raw_seq_length = Biostrings::width(seqs_raw)
)

# 1.3 Read OTU table ----------------------------------------------------------

metadata_raw <- readxl::read_excel(workbook, sheet = "2_reads_after_controls")

names(metadata_raw) <- trimws(names(metadata_raw))

pident_col <- which(tolower(names(metadata_raw)) == "pident")

if (length(pident_col) == 0) {
  stop("Column 'pident' not found.")
}

sample_cols <- names(metadata_raw)[(pident_col + 1):ncol(metadata_raw)]

metadata <- metadata_raw %>%
  dplyr::rename(seq_id = qseqid)

# 1.4 Keep OTUs present in both files -----------------------------------------

common_otus <- intersect(
  raw_seq_df$seq_id,
  metadata$seq_id
)

raw_seq_df <- raw_seq_df %>%
  dplyr::filter(seq_id %in% common_otus)

metadata <- metadata %>%
  dplyr::filter(seq_id %in% common_otus)

# 1.5 Abundance per OTU -------------------------------------------------------

metadata <- metadata %>%
  dplyr::mutate(
    SUM = rowSums(
      dplyr::across(
        dplyr::all_of(sample_cols),
        ~ as.numeric(.x)
      ),
      na.rm = TRUE
    )
  )

# 1.6 Join sequences and OTU table --------------------------------------------

data <- raw_seq_df %>%
  dplyr::inner_join(metadata, by = "seq_id") %>%
  dplyr::mutate(
    raw_N_count = stringr::str_count(raw_sequence, "N")
  )

# 1.7 Stop-codon functions ----------------------------------------------------

# Vertebrate mitochondrial code: TGA = W; stops = TAA, TAG, AGA, AGG.

stop_codons_mito_vertebrate <- c("TAA", "TAG", "AGA", "AGG")

clean_raw_sequence <- function(seq) {
  seq <- toupper(seq)
  seq <- gsub("-", "", seq)
  seq <- gsub("\\s+", "", seq)
  seq
}

get_codons_raw_by_frame <- function(seq, frame = 1) {
  seq <- clean_raw_sequence(seq)

  if (nchar(seq) < frame + 2) {
    return(character(0))
  }

  seq_frame <- substr(seq, frame, nchar(seq))

  remainder <- nchar(seq_frame) %% 3

  if (remainder != 0) {
    seq_frame <- substr(seq_frame, 1, nchar(seq_frame) - remainder)
  }

  if (nchar(seq_frame) < 3) {
    return(character(0))
  }

  substring(
    seq_frame,
    seq(1, nchar(seq_frame), by = 3),
    seq(3, nchar(seq_frame), by = 3)
  )
}

get_stop_codons_raw_by_frame <- function(seq, frame = 1) {
  codons <- get_codons_raw_by_frame(seq, frame = frame)

  valid_codons <- !stringr::str_detect(codons, "N")

  codons_valid <- codons[valid_codons]

  stops <- codons_valid[
    codons_valid %in% stop_codons_mito_vertebrate
  ]

  if (length(stops) == 0) {
    return("")
  }

  paste(stops, collapse = "|")
}

get_stop_positions_raw_by_frame <- function(seq, frame = 1) {
  codons <- get_codons_raw_by_frame(seq, frame = frame)

  stop_index <- which(
    !stringr::str_detect(codons, "N") &
      codons %in% stop_codons_mito_vertebrate
  )

  if (length(stop_index) == 0) {
    return("")
  }

  positions <- frame + ((stop_index - 1) * 3)

  paste(positions, collapse = "|")
}

count_stops <- function(stop_string) {
  if (is.na(stop_string) || stop_string == "") {
    return(0)
  }

  stringr::str_count(stop_string, "\\|") + 1
}

is_terminal_stop <- function(seq, stop_positions_string) {
  seq_clean <- clean_raw_sequence(seq)
  seq_len <- nchar(seq_clean)

  if (is.na(stop_positions_string) || stop_positions_string == "") {
    return(FALSE)
  }

  positions <- as.numeric(strsplit(stop_positions_string, "\\|")[[1]])

  any(
    positions == 1 |
      positions >= (seq_len - 2)
  )
}

remove_terminal_stop_positions <- function(seq, stop_positions_string) {
  seq_clean <- clean_raw_sequence(seq)
  seq_len <- nchar(seq_clean)

  if (is.na(stop_positions_string) || stop_positions_string == "") {
    return("")
  }

  positions <- as.numeric(strsplit(stop_positions_string, "\\|")[[1]])

  internal_positions <- positions[
    !(positions == 1 | positions >= (seq_len - 2))
  ]

  if (length(internal_positions) == 0) {
    return("")
  }

  paste(internal_positions, collapse = "|")
}

remove_terminal_stop_codons <- function(seq, stop_positions_string, stop_codons_string) {
  seq_clean <- clean_raw_sequence(seq)
  seq_len <- nchar(seq_clean)

  if (
    is.na(stop_positions_string) || stop_positions_string == "" ||
      is.na(stop_codons_string) || stop_codons_string == ""
  ) {
    return("")
  }

  positions <- as.numeric(strsplit(stop_positions_string, "\\|")[[1]])
  codons <- strsplit(stop_codons_string, "\\|")[[1]]

  keep <- !(positions == 1 | positions >= (seq_len - 2))

  codons_internal <- codons[keep]

  if (length(codons_internal) == 0) {
    return("")
  }

  paste(codons_internal, collapse = "|")
}

# 1.8 Trimming function -------------------------------------------------------

# Cut at the first stop in the best frame and keep the longest fragment
# (tie: fragment before the stop); the stop codon is removed.

select_longest_fragment_around_first_stop <- function(seq, frame = 1) {
  seq_clean <- clean_raw_sequence(seq)

  all_positions <- get_stop_positions_raw_by_frame(
    seq_clean,
    frame = frame
  )

  if (all_positions == "") {
    return(seq_clean)
  }

  positions <- as.numeric(strsplit(all_positions, "\\|")[[1]])

  first_stop <- positions[1]

  fragment_before <- substr(
    seq_clean,
    1,
    first_stop - 1
  )

  fragment_after <- substr(
    seq_clean,
    first_stop + 3,
    nchar(seq_clean)
  )

  length_before <- nchar(fragment_before)
  length_after <- nchar(fragment_after)

  if (length_after > length_before) {
    return(fragment_after)
  }

  fragment_before
}

get_fragment_selected_after_stop_cut <- function(seq, frame = 1) {
  seq_clean <- clean_raw_sequence(seq)

  all_positions <- get_stop_positions_raw_by_frame(
    seq_clean,
    frame = frame
  )

  if (all_positions == "") {
    return("no_stop")
  }

  positions <- as.numeric(strsplit(all_positions, "\\|")[[1]])
  first_stop <- positions[1]

  fragment_before <- substr(seq_clean, 1, first_stop - 1)
  fragment_after <- substr(seq_clean, first_stop + 3, nchar(seq_clean))

  length_before <- nchar(fragment_before)
  length_after <- nchar(fragment_after)

  if (length_after > length_before) {
    return("after_stop")
  }

  if (length_before == length_after) {
    return("before_stop_tie")
  }

  "before_stop"
}

get_first_stop_position <- function(seq, frame = 1) {
  seq_clean <- clean_raw_sequence(seq)

  all_positions <- get_stop_positions_raw_by_frame(
    seq_clean,
    frame = frame
  )

  if (all_positions == "") {
    return(NA_real_)
  }

  positions <- as.numeric(strsplit(all_positions, "\\|")[[1]])

  positions[1]
}

get_fragment_before_first_stop_length <- function(seq, frame = 1) {
  seq_clean <- clean_raw_sequence(seq)

  first_stop <- get_first_stop_position(seq_clean, frame = frame)

  if (is.na(first_stop)) {
    return(nchar(seq_clean))
  }

  nchar(substr(seq_clean, 1, first_stop - 1))
}

get_fragment_after_first_stop_length <- function(seq, frame = 1) {
  seq_clean <- clean_raw_sequence(seq)

  first_stop <- get_first_stop_position(seq_clean, frame = frame)

  if (is.na(first_stop)) {
    return(0)
  }

  nchar(substr(seq_clean, first_stop + 3, nchar(seq_clean)))
}

# 1.9 Stop-codon screening in three reading frames ----------------------------

stop_screen <- data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    stop_codons_frame1 = get_stop_codons_raw_by_frame(raw_sequence, frame = 1),
    stop_codons_frame2 = get_stop_codons_raw_by_frame(raw_sequence, frame = 2),
    stop_codons_frame3 = get_stop_codons_raw_by_frame(raw_sequence, frame = 3),

    stop_positions_frame1 = get_stop_positions_raw_by_frame(raw_sequence, frame = 1),
    stop_positions_frame2 = get_stop_positions_raw_by_frame(raw_sequence, frame = 2),
    stop_positions_frame3 = get_stop_positions_raw_by_frame(raw_sequence, frame = 3),

    n_stop_frame1 = count_stops(stop_codons_frame1),
    n_stop_frame2 = count_stops(stop_codons_frame2),
    n_stop_frame3 = count_stops(stop_codons_frame3),

    best_frame = which.min(c(
      n_stop_frame1,
      n_stop_frame2,
      n_stop_frame3
    )),

    best_frame_n_stops_raw = min(c(
      n_stop_frame1,
      n_stop_frame2,
      n_stop_frame3
    )),

    stop_codons_best_frame_raw = dplyr::case_when(
      best_frame == 1 ~ stop_codons_frame1,
      best_frame == 2 ~ stop_codons_frame2,
      best_frame == 3 ~ stop_codons_frame3,
      TRUE ~ ""
    ),

    stop_positions_best_frame_raw = dplyr::case_when(
      best_frame == 1 ~ stop_positions_frame1,
      best_frame == 2 ~ stop_positions_frame2,
      best_frame == 3 ~ stop_positions_frame3,
      TRUE ~ ""
    ),

    terminal_stop_in_best_frame = is_terminal_stop(
      raw_sequence,
      stop_positions_best_frame_raw
    ),

    stop_codons_best_frame = remove_terminal_stop_codons(
      raw_sequence,
      stop_positions_best_frame_raw,
      stop_codons_best_frame_raw
    ),

    stop_positions_best_frame = remove_terminal_stop_positions(
      raw_sequence,
      stop_positions_best_frame_raw
    ),

    best_frame_n_stops = count_stops(stop_codons_best_frame),

    has_stop_codon = best_frame_n_stops > 0,

    first_stop_position_best_frame_raw = get_first_stop_position(
      raw_sequence,
      frame = best_frame
    ),

    fragment_before_first_stop_length = get_fragment_before_first_stop_length(
      raw_sequence,
      frame = best_frame
    ),

    fragment_after_first_stop_length = get_fragment_after_first_stop_length(
      raw_sequence,
      frame = best_frame
    ),

    selected_fragment_after_stop_cut = get_fragment_selected_after_stop_cut(
      raw_sequence,
      frame = best_frame
    ),

    trimmed_sequence_longest_fragment = select_longest_fragment_around_first_stop(
      raw_sequence,
      frame = best_frame
    ),

    trimmed_seq_length = nchar(trimmed_sequence_longest_fragment),

    stop_codon_interpretation = dplyr::case_when(
      best_frame_n_stops == 0 &
        best_frame_n_stops_raw > 0 &
        terminal_stop_in_best_frame ~
        "no_internal_stop_terminal_stop_only_longest_fragment_selected",

      best_frame_n_stops == 0 ~
        "no_stop_in_best_frame",

      best_frame_n_stops > 0 &
        terminal_stop_in_best_frame ~
        "internal_stop_plus_terminal_stop_longest_fragment_selected",

      n_stop_frame1 > 0 &
        n_stop_frame2 > 0 &
        n_stop_frame3 > 0 ~
        "internal_stop_in_all_frames_possible_problem",

      TRUE ~
        "internal_stop_depends_on_frame_longest_fragment_selected"
    )
  ) %>%
  dplyr::ungroup()

# 1.10 Stop-codon screening output --------------------------------------------

stop_output <- stop_screen %>%
  dplyr::select(
    seq_id,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    pident,
    dplyr::all_of(sample_cols),
    SUM,
    raw_seq_length,
    raw_N_count,
    stop_codons_frame1,
    stop_codons_frame2,
    stop_codons_frame3,
    stop_positions_frame1,
    stop_positions_frame2,
    stop_positions_frame3,
    n_stop_frame1,
    n_stop_frame2,
    n_stop_frame3,
    best_frame,
    best_frame_n_stops_raw,
    stop_codons_best_frame_raw,
    stop_positions_best_frame_raw,
    terminal_stop_in_best_frame,
    best_frame_n_stops,
    stop_codons_best_frame,
    stop_positions_best_frame,
    has_stop_codon,
    first_stop_position_best_frame_raw,
    fragment_before_first_stop_length,
    fragment_after_first_stop_length,
    selected_fragment_after_stop_cut,
    stop_codon_interpretation,
    trimmed_seq_length,
    raw_sequence,
    trimmed_sequence_longest_fragment
  )

# 1.11 Trimmed sequences ------------------------------------------------------

trimmed_fasta <- Biostrings::DNAStringSet(
  stop_output$trimmed_sequence_longest_fragment
)

names(trimmed_fasta) <- stop_output$seq_id

# 1.12 Simplified stop-codon table --------------------------------------------

stop_output_simple <- stop_output %>%
  dplyr::select(
    seq_id,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    pident,
    dplyr::all_of(sample_cols),
    SUM,

    has_stop_codon,
    stop_codon_interpretation,

    first_stop_position_best_frame_raw,
    fragment_before_first_stop_length,
    fragment_after_first_stop_length,
    selected_fragment_after_stop_cut,

    raw_seq_length,
    raw_sequence,

    trimmed_seq_length,
    trimmed_sequence_longest_fragment
  )

# 1.13 Clean memory -----------------------------------------------------------

rm(
  list = setdiff(
    ls(),
    c(
      "stop_output_simple"
    )
  )
)

gc()

# Part 2. OTU pre-screening ===================================================
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(Biostrings)
library(DECIPHER)

options(scipen = 999)

# 2.1 Thresholds --------------------------------------------------------------

short_fragment_threshold <- 330
rare_threshold <- 0.01
very_rare_threshold <- 0.001
terminal_window <- 30

high_difference_threshold <- 20
middle_bias_threshold <- 0.5

# 2.2 Manual OTU groups -------------------------------------------------------

manual_groups <- tibble::tribble(
  ~seq_id,      ~manual_group,        ~manual_group_order,

  "Otu122046",  "Aequidens",          1,
  "Otu79958",   "Aequidens",          1,
  "Otu102086",  "Astyanax",           2,
  "Otu163370",  "Astyanax",           2,
  "Otu36201",   "Astyanax",           2,
  "Otu49654",   "Astyanax",           2,
  "Otu66060",   "Astyanax",           2,
  "Otu93479",   "Astyanax",           2,
  "Otu116104",  "Brachyhypopomus",    3,
  "Otu260877",  "Brachyhypopomus",    3,
  "Otu94629",   "Brachyhypopomus",    3,
  "Otu276732",  "Acanthodoras",       4,
  "Otu94079",   "Acanthodoras",       4,
  "Otu268636",  "Bryconops",          5,
  "Otu58549",   "Bryconops",          5,
  "Otu268267",  "Bryconops",          5,
  "Otu269429",  "Bryconops",          5,
  "Otu307567",  "Bryconops",          5,
  "Otu33065",   "Bryconops",          5,
  "Otu93748",   "Bryconops",          5,
  "Otu102682",  "Cichla",             6,
  "Otu154087",  "Cichla",             6,
  "Otu143770",  "Cichlasoma",         7,
  "Otu235776",  "Cichlasoma",         7,
  "Otu263649",  "Cichlasoma",         7,
  "Otu311766",  "Cichlasoma",         7,
  "Otu36459",   "Cichlasoma",         7,
  "Otu78631",   "Cichlasoma",         7,
  "Otu226096",  "Copella",            8,
  "Otu69290",   "Copella",            8,
  "Otu105641",  "Crenicichla",        9,
  "Otu6939",    "Crenicichla",        9,
  "Otu245364",  "Curimatopsis",       10,
  "Otu287554",  "Curimatopsis",       10,
  "Otu78645",   "Curimatopsis",       10,
  "Otu174316",  "Eigenmannia",        11,
  "Otu56611",   "Eigenmannia",        11,
  "Otu193643",  "Eigenmannia",        11,
  "Otu235230",  "Eigenmannia",        11,
  "Otu83393",   "Eigenmannia",        11,
  "Otu222626",  "Erythrinus",         12,
  "Otu264232",  "Erythrinus",         12,
  "Otu87916",   "Erythrinus",         12,
  "Otu123130",  "Geophagus",          13,
  "Otu255045",  "Geophagus",          13,
  "Otu295196",  "Geophagus",          13,
  "Otu93614",   "Geophagus",          13,
  "Otu192938",  "Geophagus",          13,
  "Otu140757",  "Gymnorhamphichthys", 14,
  "Otu184343",  "Gymnorhamphichthys", 14,
  "Otu93836",   "Gymnorhamphichthys", 14,
  "Otu102345",  "Gymnotus",           15,
  "Otu105065",  "Gymnotus",           15,
  "Otu286035",  "Gymnotus",           15,
  "Otu93835",   "Gymnotus",           15,
  "Otu94967",   "Gymnotus",           15,
  "Otu106524",  "Hoplerythrinus",     16,
  "Otu282914",  "Hoplerythrinus",     16,
  "Otu311137",  "Hoplerythrinus",     16,
  "Otu314662",  "Hoplerythrinus",     16,
  "Otu64406",   "Hoplerythrinus",     16,
  "Otu102116",  "Hoplias",            17,
  "Otu2235",    "Hoplias",            17,
  "Otu234067",  "Hoplias",            17,
  "Otu255052",  "Hoplias",            17,
  "Otu260774",  "Hoplias",            17,
  "Otu299260",  "Hoplias",            17,
  "Otu30735",   "Hoplias",            17,
  "Otu307471",  "Hoplias",            17,
  "Otu98035",   "Hoplias",            17,
  "Otu162714",  "Hyphessobrycon",     18,
  "Otu218869",  "Hyphessobrycon",     18,
  "Otu106862",  "Leporinus",          19,
  "Otu111220",  "Leporinus",          19,
  "Otu159448",  "Leporinus",          19,
  "Otu94639",   "Leporinus",          19,
  "Otu101709",  "Metynnis",           20,
  "Otu104952",  "Metynnis",           20,
  "Otu115058",  "Myloplus",           21,
  "Otu161004",  "Myloplus",           21,
  "Otu114699",  "Pimelodella",        22,
  "Otu94377",   "Pimelodella",        22,
  "Otu132078",  "Rhamdia",            23,
  "Otu93901",   "Rhamdia",            23,
  "Otu289706",  "Satanoperca",        24,
  "Otu84731",   "Satanoperca",        24,
  "Otu268728",  "Sternopygus",        25,
  "Otu79064",   "Sternopygus",        25,
  "Otu94207",   "Sternopygus",        25,
  "Otu101283",  "Trachelyopterus",    26,
  "Otu111027",  "Trachelyopterus",    26
)

# 2.3 Sample columns ----------------------------------------------------------

sample_cols <- names(stop_output_simple)[
  (which(names(stop_output_simple) == "pident") + 1):
    (which(names(stop_output_simple) == "SUM") - 1)
]

# 2.4 OTU metrics -------------------------------------------------------------

global_curation <- stop_output_simple %>%
  dplyr::left_join(manual_groups, by = "seq_id") %>%
  dplyr::mutate(
    in_manual_group = !is.na(manual_group),

    manual_group = dplyr::if_else(
      is.na(manual_group),
      "not_manually_grouped",
      manual_group
    ),

    manual_group_order = dplyr::if_else(
      is.na(manual_group_order),
      999L,
      manual_group_order
    ),

    RRA = SUM / sum(SUM, na.rm = TRUE),

    short_fragment = trimmed_seq_length <= short_fragment_threshold,

    very_low_abundance = RRA < very_rare_threshold | SUM == 1,

    low_abundance = RRA < rare_threshold,

    incomplete_taxonomy =
      is.na(Species) | Species == "" |
      is.na(Genus) | Genus == "",

    translation_problem = has_stop_codon %in% TRUE
  )

# 2.5 Group metrics -----------------------------------------------------------

manual_group_metrics <- global_curation %>%
  dplyr::filter(in_manual_group) %>%
  dplyr::group_by(manual_group) %>%
  dplyr::mutate(
    group_total_reads = sum(SUM, na.rm = TRUE),
    group_RRA = SUM / group_total_reads,
    group_dominant_SUM = max(SUM, na.rm = TRUE),
    dominant_within_group = SUM == group_dominant_SUM,
    abundance_ratio_to_group_dominant = SUM / group_dominant_SUM
  ) %>%
  dplyr::ungroup()

not_grouped_metrics <- global_curation %>%
  dplyr::filter(!in_manual_group) %>%
  dplyr::mutate(
    group_total_reads = NA_real_,
    group_RRA = NA_real_,
    group_dominant_SUM = NA_real_,
    dominant_within_group = NA,
    abundance_ratio_to_group_dominant = NA_real_
  )

global_curation <- dplyr::bind_rows(
  manual_group_metrics,
  not_grouped_metrics
)

# 2.6 Co-occurrence with the dominant OTU -------------------------------------

dominant_samples_by_group <- global_curation %>%
  dplyr::filter(in_manual_group, dominant_within_group) %>%
  dplyr::select(
    manual_group,
    dplyr::all_of(sample_cols)
  ) %>%
  dplyr::rename_with(
    ~ paste0("dominant_", .x),
    dplyr::all_of(sample_cols)
  )

global_curation <- global_curation %>%
  dplyr::left_join(dominant_samples_by_group, by = "manual_group") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    cooccurs_with_group_dominant = ifelse(
      !in_manual_group,
      NA,
      any(
        as.numeric(c(dplyr::across(dplyr::all_of(sample_cols)))) > 0 &
          as.numeric(c(dplyr::across(dplyr::all_of(paste0("dominant_", sample_cols))))) > 0,
        na.rm = TRUE
      )
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    -dplyr::all_of(paste0("dominant_", sample_cols))
  )

# 2.7 Sequence-difference function --------------------------------------------

difference_profile <- function(query_seq, reference_seq, window = 30) {
  q <- strsplit(query_seq, "")[[1]]
  r <- strsplit(reference_seq, "")[[1]]

  min_len <- min(length(q), length(r))
  q <- q[1:min_len]
  r <- r[1:min_len]

  valid <- !q %in% c("N", "-") & !r %in% c("N", "-")
  diff_pos <- which(valid & q != r)

  beginning_pos <- diff_pos[diff_pos <= window]
  end_pos <- diff_pos[diff_pos > (min_len - window)]
  terminal_pos <- union(beginning_pos, end_pos)
  middle_pos <- diff_pos[!(diff_pos %in% terminal_pos)]

  total_differences <- length(diff_pos)

  tibble::tibble(
    differences_to_group_dominant = total_differences,
    dominant_beginning_differences = length(beginning_pos),
    dominant_middle_differences = length(middle_pos),
    dominant_end_differences = length(end_pos),
    dominant_terminal_differences = length(terminal_pos),

    dominant_beginning_bias = ifelse(total_differences > 0, length(beginning_pos) / total_differences, 0),
    dominant_middle_bias = ifelse(total_differences > 0, length(middle_pos) / total_differences, 0),
    dominant_end_bias = ifelse(total_differences > 0, length(end_pos) / total_differences, 0),
    dominant_terminal_bias = ifelse(total_differences > 0, length(terminal_pos) / total_differences, 0),

    dominant_difference_positions = ifelse(
      total_differences > 0,
      paste(diff_pos, collapse = "|"),
      ""
    ),

    dominant_difference_pattern = dplyr::case_when(
      total_differences == 0 ~ "identical",
      length(terminal_pos) / total_differences >= 0.5 ~ "terminal_concentrated",
      length(middle_pos) / total_differences >= middle_bias_threshold ~ "middle_concentrated",
      TRUE ~ "mixed_distribution"
    )
  )
}

# 2.8 Alignment against the dominant OTU of each group ------------------------

get_group_dominant_metrics <- function(group_df, window = 30) {

  group_df <- group_df %>%
    dplyr::mutate(
      dominant_within_group = SUM == max(SUM, na.rm = TRUE)
    )

  seqs <- Biostrings::DNAStringSet(group_df$trimmed_sequence_longest_fragment)
  names(seqs) <- group_df$seq_id

  aligned <- DECIPHER::AlignSeqs(seqs, processors = 1)

  aligned_df <- tibble::tibble(
    seq_id = names(aligned),
    aligned_sequence = as.character(aligned),
    aligned_seq_length = Biostrings::width(aligned)
  )

  group_df2 <- group_df %>%
    dplyr::left_join(aligned_df, by = "seq_id")

  dominant_seq <- group_df2 %>%
    dplyr::filter(dominant_within_group %in% TRUE) %>%
    dplyr::arrange(dplyr::desc(SUM)) %>%
    dplyr::slice(1) %>%
    dplyr::pull(aligned_sequence)

  if (length(dominant_seq) == 0 || is.na(dominant_seq) || dominant_seq == "") {
    dominant_seq <- group_df2 %>%
      dplyr::arrange(dplyr::desc(SUM)) %>%
      dplyr::slice(1) %>%
      dplyr::pull(aligned_sequence)
  }

  dominant_profiles <- group_df2 %>%
    dplyr::rowwise() %>%
    dplyr::do({
      prof <- difference_profile(
        .$aligned_sequence,
        dominant_seq,
        window = window
      )

      dplyr::bind_cols(
        tibble::tibble(seq_id = .$seq_id),
        prof
      )
    }) %>%
    dplyr::ungroup()

  group_df2 %>%
    dplyr::left_join(dominant_profiles, by = "seq_id")
}

manual_dominant_metrics <- global_curation %>%
  dplyr::filter(in_manual_group) %>%
  dplyr::group_split(manual_group) %>%
  lapply(get_group_dominant_metrics, window = terminal_window) %>%
  dplyr::bind_rows()

not_grouped_dominant_metrics <- global_curation %>%
  dplyr::filter(!in_manual_group) %>%
  dplyr::mutate(
    aligned_sequence = NA_character_,
    aligned_seq_length = NA_integer_,
    differences_to_group_dominant = NA_integer_,
    dominant_beginning_differences = NA_integer_,
    dominant_middle_differences = NA_integer_,
    dominant_end_differences = NA_integer_,
    dominant_terminal_differences = NA_integer_,
    dominant_beginning_bias = NA_real_,
    dominant_middle_bias = NA_real_,
    dominant_end_bias = NA_real_,
    dominant_terminal_bias = NA_real_,
    dominant_difference_positions = NA_character_,
    dominant_difference_pattern = "not_evaluated_not_in_manual_group"
  )

global_curation <- dplyr::bind_rows(
  manual_dominant_metrics,
  not_grouped_dominant_metrics
)

# 2.9 Flags from the comparison with the dominant OTU -------------------------

global_curation <- global_curation %>%
  dplyr::mutate(
    terminal_concentrated =
      dominant_difference_pattern == "terminal_concentrated",

    middle_concentrated =
      dominant_difference_pattern == "middle_concentrated",

    identical_to_group_dominant =
      dominant_difference_pattern == "identical",

    high_difference_from_dominant =
      differences_to_group_dominant >= high_difference_threshold,

    divergent_middle_pattern =
      high_difference_from_dominant %in% TRUE &
        middle_concentrated %in% TRUE,

    block_auto_merge_due_divergence =
      divergent_middle_pattern %in% TRUE &
        !dominant_within_group %in% TRUE
  )

# 2.10 Pre-screening decision -------------------------------------------------

global_curation <- global_curation %>%
  dplyr::mutate(
    artefact_score =
      as.integer(very_low_abundance %in% TRUE) +
      as.integer(low_abundance %in% TRUE) +
      as.integer(short_fragment %in% TRUE) +
      as.integer(cooccurs_with_group_dominant %in% TRUE) +
      as.integer(incomplete_taxonomy %in% TRUE) +
      as.integer(terminal_concentrated %in% TRUE),

    curation_decision = dplyr::case_when(
      has_stop_codon %in% TRUE &
        dominant_within_group %in% TRUE ~
        "to_blast_revalidate_dominant_with_internal_stop_codon",

      has_stop_codon %in% TRUE ~
        "to_blast_validate_internal_stop_codon",

      !in_manual_group ~
        "outside_manual_group",

      dominant_within_group ~
        "keep_dominant",

      identical_to_group_dominant %in% TRUE ~
        "merge_dominant",

      block_auto_merge_due_divergence %in% TRUE ~
        "candidate_distinct_taxon",

      SUM <= 1 &
        terminal_concentrated %in% TRUE ~
        "exclude_artefact",

      terminal_concentrated %in% TRUE &
        low_abundance %in% TRUE ~
        "exclude_artefact",

      artefact_score >= 6 &
        !block_auto_merge_due_divergence %in% TRUE ~
        "exclude_artefact",

      low_abundance %in% TRUE &
        cooccurs_with_group_dominant %in% TRUE &
        !translation_problem %in% TRUE &
        !block_auto_merge_due_divergence %in% TRUE ~
        "merge_dominant",
      
      low_abundance %in% TRUE &
        middle_concentrated %in% TRUE &
        differences_to_group_dominant < high_difference_threshold &
        !translation_problem %in% TRUE &
        !block_auto_merge_due_divergence %in% TRUE ~
        "merge_dominant",

      low_abundance %in% TRUE &
        !translation_problem %in% TRUE ~
        "review_rare_variant",

      TRUE ~
        "review"
    ),

    interpretation = dplyr::case_when(
      curation_decision == "to_blast_revalidate_dominant_with_internal_stop_codon" ~
        "dominant_sequence_has_internal_stop_codon_validate_trimmed_taxonomy_and_reassign_reference",

      curation_decision == "to_blast_validate_internal_stop_codon" ~
        "internal_stop_codon_detected_validate_trimmed_taxonomy",

      curation_decision == "outside_manual_group" ~
        "outside_manual_group",

      curation_decision == "keep_dominant" ~
        "dominant_reference",

      curation_decision == "merge_dominant" &
        identical_to_group_dominant %in% TRUE ~
        "identical_to_dominant",

      curation_decision == "candidate_distinct_taxon" ~
        "divergent_internal_pattern",

      curation_decision == "exclude_artefact" &
        terminal_concentrated %in% TRUE ~
        "terminal_noise_pattern",

      curation_decision == "exclude_artefact" ~
        "technical_artefact",

      curation_decision == "merge_dominant" ~
        "likely_dominant_variant",

      curation_decision == "review_rare_variant" ~
        "rare_variant",

      TRUE ~
        "inconclusive"
    ),

    post_curation_action = dplyr::case_when(
      curation_decision == "to_blast_revalidate_dominant_with_internal_stop_codon" ~
        "blast_trimmed_sequence_then_reassign_reference",

      curation_decision == "to_blast_validate_internal_stop_codon" ~
        "blast_trimmed_sequence",

      curation_decision == "candidate_distinct_taxon" ~
        "validate_taxonomic_identity",

      TRUE ~
        "none"
    ),
     should_merge = curation_decision == "merge_dominant",

    should_keep_as_representative = curation_decision %in% c(
      "keep_dominant",
      "to_blast_revalidate_dominant_with_internal_stop_codon",
      "candidate_distinct_taxon",
      "review",
      "review_rare_variant",
      "outside_manual_group"
    ),

    merge_status = dplyr::case_when(
      curation_decision == "merge_dominant" ~ "merge_into_representative",
      curation_decision == "keep_dominant" ~ "representative",
      curation_decision == "candidate_distinct_taxon" ~ "do_not_merge_distinct_taxon",
      stringr::str_detect(curation_decision, "blast") ~ "do_not_merge_pending_blast_validation",
      curation_decision == "outside_manual_group" ~ "do_not_merge_outside_manual_group",
      TRUE ~ "do_not_merge_review"
    )
  )

# 2.11 Pre-screening output ---------------------------------------------------

global_curation_output <- global_curation %>%
  dplyr::arrange(
    manual_group_order,
    manual_group,
    dplyr::desc(in_manual_group),
    dplyr::desc(dominant_within_group),
    dplyr::desc(SUM)
  ) %>%
  dplyr::select(
    seq_id,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    in_manual_group,
    manual_group,
    manual_group_order,
    pident,
    dplyr::all_of(sample_cols),
    SUM,
    RRA,
    group_total_reads,
    group_RRA,
    group_dominant_SUM,
    abundance_ratio_to_group_dominant,
    dominant_within_group,
    cooccurs_with_group_dominant,
    has_stop_codon,
    stop_codon_interpretation,
    post_curation_action,
    raw_seq_length,
    trimmed_seq_length,
    short_fragment,
    very_low_abundance,
    low_abundance,
    incomplete_taxonomy,
    aligned_seq_length,

    differences_to_group_dominant,
    high_difference_from_dominant,
    dominant_beginning_differences,
    dominant_middle_differences,
    dominant_end_differences,
    dominant_terminal_differences,
    dominant_beginning_bias,
    dominant_middle_bias,
    dominant_end_bias,
    dominant_terminal_bias,
    dominant_difference_pattern,
    dominant_difference_positions,

    terminal_concentrated,
    middle_concentrated,
    divergent_middle_pattern,
    block_auto_merge_due_divergence,
    identical_to_group_dominant,
    artefact_score,
    curation_decision,
    interpretation,
    raw_sequence,
    trimmed_sequence_longest_fragment,
    aligned_sequence
  )

rm(
  list = setdiff(
    ls(),
    c(
      "stop_output_simple","global_curation_output"
    )
  )
)

gc()

# Part 3. BOLD and GenBank validation =========================================

# 3.1 Packages ----------------------------------------------------------------

library(Biostrings)
library(dplyr)
library(stringr)
library(httr)
library(xml2)
library(chromote)
library(tibble)

source("02_otu_prescreening_functions.R")

# 3.2 Helper: safe scalar character -------------------------------------------

safe_scalar_chr <- function(x, fallback) {
  
  if (is.null(x) || length(x) == 0) {
    return(fallback)
  }
  
  x <- as.character(x)[1]
  
  if (is.na(x) || x == "" || identical(x, "NA")) {
    return(fallback)
  }
  
  x
}

# 3.3 OTUs to validate --------------------------------------------------------

blast_rows <- which(
  global_curation_output$curation_decision %in% c(
    "to_blast_validate_internal_stop_codon",
    "to_blast_revalidate_dominant_with_internal_stop_codon"
  )
)

# 3.4 Working copy ------------------------------------------------------------

if (!"blast" %in% names(global_curation_output)) {
  global_curation_output$blast <- NA_character_
}

global_curation_output$blast <- as.character(
  global_curation_output$blast
)

global_curation_output <- global_curation_output %>%
  relocate(
    blast,
    .after = interpretation
  )

# 3.5 Audit table -------------------------------------------------------------

bold_genbank_audit <- list()

# 3.6 Browser session for the BOLD ID Engine ----------------------------------

b <- chromote::ChromoteSession$new()

# 3.7 Run validation ----------------------------------------------------------

for (i in blast_rows) {
  
  seq_id_i <- global_curation_output$seq_id[i]
  
  manual_group_i <- dplyr::coalesce(
    global_curation_output$manual_group[i],
    "not_manually_grouped"
  )
  
  curation_decision_i <- global_curation_output$curation_decision[i]
  
  query_id_i <- paste(
    seq_id_i,
    manual_group_i,
    curation_decision_i,
    sep = "|"
  )
  
  seq_i <- global_curation_output$trimmed_sequence_longest_fragment[i]
  seq_i <- safe_scalar_chr(seq_i, fallback = NA_character_)
  seq_i <- toupper(gsub("[^ACGTN]", "", seq_i))
  
  message("\n------------------------------")
  message("ROW: ", i)
  message("SEQ_ID: ", seq_id_i)
  message("QUERY_ID: ", query_id_i)
  message("SEQ_LENGTH: ", nchar(seq_i))
  
  bold_result_i <- NULL
  bold_hits_i <- NULL
  bold_hit_i <- NA_character_
  genbank_hit_i <- NA_character_
  fasta_path_i <- NA_character_
  fasta_deleted_i <- NA
  
  # Skip very short sequences
  
  if (is.na(seq_i) || nchar(seq_i) < 50) {
    
    bold_hit_i <- paste0("BOLD:sequence_too_short_", nchar(seq_i), "bp")
    genbank_hit_i <- paste0("GenBank:sequence_too_short_", nchar(seq_i), "bp")
    
  } else {
    
    # GenBank (NCBI BLAST) first; kept even if BOLD fails
    
    message("Submitting to GenBank BLAST: ", query_id_i)
    
    genbank_hit_i <- tryCatch(
      {
        out <- blast_genbank_top_taxon(
          sequence = seq_i
        )
        
        safe_scalar_chr(
          out,
          fallback = "GenBank:no_hit_or_unparsed_result"
        )
      },
      error = function(e) {
        message("GenBank BLAST error: ", e$message)
        "GenBank:error"
      }
    )
    
    genbank_hit_i <- safe_scalar_chr(
      genbank_hit_i,
      fallback = "GenBank:no_hit_or_unparsed_result"
    )
    
    # BOLD ID Engine
    
    fasta_i <- Biostrings::DNAStringSet(seq_i)
    names(fasta_i) <- query_id_i
    
    fasta_path_i <- tempfile(
      pattern = paste0("to_bold_id_engine_", seq_id_i, "_"),
      fileext = ".fasta"
    )
    
    Biostrings::writeXStringSet(
      fasta_i,
      filepath = fasta_path_i
    )
    
    message("Temporary FASTA created: ", fasta_path_i)
    message("Submitting to BOLD ID Engine: ", query_id_i)
    
    bold_result_i <- tryCatch(
      {
        run_bold_id_engine(
          b = b,
          fasta_path = fasta_path_i
        )
      },
      error = function(e) {
        message("BOLD ID Engine error: ", e$message)
        NULL
      }
    )
    
    bold_hits_i <- tryCatch(
      {
        if (is.null(bold_result_i) || is.null(bold_result_i$result_text)) {
          NULL
        } else {
          parse_bold_hits_v2(
            bold_result_i$result_text
          )
        }
      },
      error = function(e) {
        message("BOLD parse error: ", e$message)
        NULL
      }
    )
    
    bold_hit_i <- tryCatch(
      {
        if (
          is.null(bold_hits_i) ||
          !is.data.frame(bold_hits_i) ||
          nrow(bold_hits_i) == 0
        ) {
          "BOLD:no_hit_or_unparsed_result"
        } else {
          get_simple_bold_hit(bold_hits_i)
        }
      },
      error = function(e) {
        message("BOLD summary error: ", e$message)
        "BOLD:summary_error"
      }
    )
    
    bold_hit_i <- safe_scalar_chr(
      bold_hit_i,
      fallback = "BOLD:no_hit_or_unparsed_result"
    )
    
    # Remove temporary FASTA
    
    if (!is.na(fasta_path_i) && file.exists(fasta_path_i)) {
      unlink(fasta_path_i)
    }
    
    fasta_deleted_i <- !file.exists(fasta_path_i)
    
    Sys.sleep(3)
  }
  
  # Check before writing the blast column
  
  genbank_hit_i <- safe_scalar_chr(
    genbank_hit_i,
    fallback = "GenBank:no_hit_or_unparsed_result"
  )
  
  bold_hit_i <- safe_scalar_chr(
    bold_hit_i,
    fallback = "BOLD:no_hit_or_unparsed_result"
  )
  
  combined_hit_i <- paste(
    genbank_hit_i,
    bold_hit_i,
    sep = "|"
  )
  
  old_blast_i <- global_curation_output$blast[i]
  
  global_curation_output$blast[i] <- combined_hit_i
  
  # Store audit record
  
  bold_genbank_audit[[as.character(seq_id_i)]] <- list(
    row = i,
    seq_id = seq_id_i,
    manual_group = manual_group_i,
    curation_decision = curation_decision_i,
    query_id = query_id_i,
    sequence_length = nchar(seq_i),
    fasta_path = fasta_path_i,
    fasta_deleted = fasta_deleted_i,
    bold_result = bold_result_i,
    bold_hits = bold_hits_i,
    bold_hit_simple = bold_hit_i,
    genbank_hit_simple = genbank_hit_i,
    combined_hit = combined_hit_i,
    old_blast = old_blast_i,
    new_blast = global_curation_output$blast[i]
  )
  
  message("GenBank result: ", genbank_hit_i)
  message("BOLD result: ", bold_hit_i)
  message("Final blast: ", global_curation_output$blast[i])
}

# 3.8 Audit summary -----------------------------------------------------------

bold_genbank_summary <- dplyr::bind_rows(
  lapply(bold_genbank_audit, function(x) {
    tibble::tibble(
      row = x$row,
      seq_id = x$seq_id,
      manual_group = x$manual_group,
      curation_decision = x$curation_decision,
      query_id = x$query_id,
      sequence_length = x$sequence_length,
      fasta_path = x$fasta_path,
      fasta_deleted = x$fasta_deleted,
      bold_hit_simple = x$bold_hit_simple,
      genbank_hit_simple = x$genbank_hit_simple,
      combined_hit = x$combined_hit,
      old_blast = x$old_blast,
      new_blast = x$new_blast
    )
  })
)

# 3.9 Write validation results to the curation table --------------------------

blast_update <- bold_genbank_summary %>%
  dplyr::select(
    seq_id,
    blast_new = new_blast
  ) %>%
  dplyr::mutate(
    blast_new = as.character(blast_new)
  )

global_curation_output <- global_curation_output %>%
  dplyr::mutate(
    blast = as.character(blast)
  ) %>%
  dplyr::left_join(
    blast_update,
    by = "seq_id"
  ) %>%
  dplyr::mutate(
    blast = dplyr::if_else(
      !is.na(blast_new) & blast_new != "",
      blast_new,
      blast
    )
  ) %>%
  dplyr::select(
    -blast_new
  )

# 3.10 Final check ------------------------------------------------------------

global_curation_output %>%
  dplyr::filter(row_number() %in% blast_rows) %>%
  dplyr::select(
    seq_id,
    manual_group,
    curation_decision,
    interpretation,
    blast
  )

rm(
  list = setdiff(
    ls(),
    c(
      "stop_output_simple","global_curation_output"
    )
  )
)

gc()

selected_manual_groups <- global_curation_output %>%
  dplyr::group_by(manual_group) %>%
  dplyr::filter(
    all(
      interpretation %in% c(
        "dominant_reference",
        "divergent_internal_pattern",
        "dominant_sequence_has_internal_stop_codon_validate_trimmed_taxonomy_and_reassign_reference"
      )
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    manual_group_order,
    dplyr::desc(interpretation == "dominant_reference")
  )

# Part 4. Pairwise comparison within manual groups ============================

options(scipen = 999)

library(dplyr)
library(tidyr)
library(tibble)

terminal_window <- 30
middle_bias_threshold <- 0.5

pairwise_difference_profile <- function(query_seq, reference_seq, window = 30) {
  
  query_seq <- as.character(query_seq)
  reference_seq <- as.character(reference_seq)
  
  if (
    length(query_seq) == 0 ||
    length(reference_seq) == 0 ||
    is.na(query_seq) ||
    is.na(reference_seq) ||
    query_seq == "" ||
    reference_seq == ""
  ) {
    return(tibble::tibble(
      pairwise_differences = NA_integer_,
      pairwise_beginning_differences = NA_integer_,
      pairwise_middle_differences = NA_integer_,
      pairwise_end_differences = NA_integer_,
      pairwise_terminal_differences = NA_integer_,
      pairwise_beginning_bias = NA_real_,
      pairwise_middle_bias = NA_real_,
      pairwise_end_bias = NA_real_,
      pairwise_terminal_bias = NA_real_,
      pairwise_difference_positions = NA_character_,
      pairwise_difference_pattern = "not_evaluated_missing_sequence"
    ))
  }
  
  q <- strsplit(query_seq, "")[[1]]
  r <- strsplit(reference_seq, "")[[1]]
  
  min_len <- min(length(q), length(r))
  
  if (min_len == 0) {
    return(tibble::tibble(
      pairwise_differences = NA_integer_,
      pairwise_beginning_differences = NA_integer_,
      pairwise_middle_differences = NA_integer_,
      pairwise_end_differences = NA_integer_,
      pairwise_terminal_differences = NA_integer_,
      pairwise_beginning_bias = NA_real_,
      pairwise_middle_bias = NA_real_,
      pairwise_end_bias = NA_real_,
      pairwise_terminal_bias = NA_real_,
      pairwise_difference_positions = NA_character_,
      pairwise_difference_pattern = "not_evaluated_empty_alignment"
    ))
  }
  
  q <- q[1:min_len]
  r <- r[1:min_len]
  
  valid <- !q %in% c("N", "-") & !r %in% c("N", "-")
  diff_pos <- which(valid & q != r)
  
  beginning_pos <- diff_pos[diff_pos <= window]
  end_pos <- diff_pos[diff_pos > (min_len - window)]
  terminal_pos <- union(beginning_pos, end_pos)
  middle_pos <- diff_pos[!(diff_pos %in% terminal_pos)]
  
  total_differences <- length(diff_pos)
  
  tibble::tibble(
    pairwise_differences = total_differences,
    pairwise_beginning_differences = length(beginning_pos),
    pairwise_middle_differences = length(middle_pos),
    pairwise_end_differences = length(end_pos),
    pairwise_terminal_differences = length(terminal_pos),
    pairwise_beginning_bias = ifelse(total_differences > 0, length(beginning_pos) / total_differences, 0),
    pairwise_middle_bias = ifelse(total_differences > 0, length(middle_pos) / total_differences, 0),
    pairwise_end_bias = ifelse(total_differences > 0, length(end_pos) / total_differences, 0),
    pairwise_terminal_bias = ifelse(total_differences > 0, length(terminal_pos) / total_differences, 0),
    pairwise_difference_positions = ifelse(
      total_differences > 0,
      paste(diff_pos, collapse = "|"),
      ""
    ),
    pairwise_difference_pattern = dplyr::case_when(
      total_differences == 0 ~ "identical",
      length(terminal_pos) / total_differences >= 0.5 ~ "terminal_concentrated",
      length(middle_pos) / total_differences >= middle_bias_threshold ~ "middle_concentrated",
      TRUE ~ "mixed_distribution"
    )
  )
}

get_group_pairwise_metrics <- function(group_df, window = 30) {
  
  group_df <- group_df %>%
    dplyr::mutate(
      pairwise_sequence_used = as.character(aligned_sequence)
    )
  
  if (nrow(group_df) < 2) {
    group_df$pairwise_compared_seq_ids <- NA_character_
    group_df$pairwise_differences <- NA_character_
    group_df$pairwise_beginning_differences <- NA_character_
    group_df$pairwise_middle_differences <- NA_character_
    group_df$pairwise_end_differences <- NA_character_
    group_df$pairwise_terminal_differences <- NA_character_
    group_df$pairwise_beginning_bias <- NA_character_
    group_df$pairwise_middle_bias <- NA_character_
    group_df$pairwise_end_bias <- NA_character_
    group_df$pairwise_terminal_bias <- NA_character_
    group_df$pairwise_difference_pattern <- "not_evaluated_less_than_two_sequences"
    group_df$pairwise_difference_positions <- NA_character_
    return(group_df)
  }
  
  pairwise_pairs <- tidyr::expand_grid(
    query_seq_id = group_df$seq_id,
    comparison_seq_id = group_df$seq_id
  ) %>%
    dplyr::filter(query_seq_id != comparison_seq_id)
  
  pairwise_long <- dplyr::bind_rows(
    lapply(seq_len(nrow(pairwise_pairs)), function(k) {
      
      query_id_k <- pairwise_pairs$query_seq_id[k]
      comparison_id_k <- pairwise_pairs$comparison_seq_id[k]
      
      query_seq_k <- group_df$pairwise_sequence_used[
        match(query_id_k, group_df$seq_id)
      ]
      
      comparison_seq_k <- group_df$pairwise_sequence_used[
        match(comparison_id_k, group_df$seq_id)
      ]
      
      prof <- pairwise_difference_profile(
        query_seq = query_seq_k,
        reference_seq = comparison_seq_k,
        window = window
      )
      
      dplyr::bind_cols(
        tibble::tibble(
          query_seq_id = query_id_k,
          comparison_seq_id = comparison_id_k
        ),
        prof
      )
    })
  )
  
  pairwise_collapsed <- pairwise_long %>%
    dplyr::group_by(query_seq_id) %>%
    dplyr::summarise(
      pairwise_compared_seq_ids = paste(comparison_seq_id, collapse = " | "),
      pairwise_differences = paste(pairwise_differences, collapse = " | "),
      pairwise_beginning_differences = paste(pairwise_beginning_differences, collapse = " | "),
      pairwise_middle_differences = paste(pairwise_middle_differences, collapse = " | "),
      pairwise_end_differences = paste(pairwise_end_differences, collapse = " | "),
      pairwise_terminal_differences = paste(pairwise_terminal_differences, collapse = " | "),
      pairwise_beginning_bias = paste(round(pairwise_beginning_bias, 4), collapse = " | "),
      pairwise_middle_bias = paste(round(pairwise_middle_bias, 4), collapse = " | "),
      pairwise_end_bias = paste(round(pairwise_end_bias, 4), collapse = " | "),
      pairwise_terminal_bias = paste(round(pairwise_terminal_bias, 4), collapse = " | "),
      pairwise_difference_pattern = paste(pairwise_difference_pattern, collapse = " | "),
      pairwise_difference_positions = paste(pairwise_difference_positions, collapse = " || "),
      .groups = "drop"
    ) %>%
    dplyr::rename(seq_id = query_seq_id)
  
  group_df %>%
    dplyr::left_join(
      pairwise_collapsed,
      by = "seq_id"
    )
}

selected_manual_groups_pairwise <- selected_manual_groups %>%
  dplyr::group_split(manual_group) %>%
  lapply(
    get_group_pairwise_metrics,
    window = terminal_window
  ) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(
    manual_group_order,
    manual_group,
    dplyr::desc(interpretation == "dominant_reference"),
    dplyr::desc(SUM)
  )

# Part 5. Candidate taxon clusters within manual groups =======================

same_taxon_max_diff <- 5
distinct_taxon_min_diff <- 15

make_group_distance_matrix <- function(group_df) {
  
  ids <- as.character(group_df$seq_id)
  seqs <- as.character(group_df$aligned_sequence)
  
  seq_matrix <- do.call(
    rbind,
    strsplit(seqs, "")
  )
  
  rownames(seq_matrix) <- ids
  
  dist_mat <- matrix(
    0,
    nrow = nrow(seq_matrix),
    ncol = nrow(seq_matrix),
    dimnames = list(ids, ids)
  )
  
  for (i in seq_len(nrow(seq_matrix))) {
    for (j in seq_len(nrow(seq_matrix))) {
      
      x <- seq_matrix[i, ]
      y <- seq_matrix[j, ]
      
      valid <- !x %in% c("-", "N") &
        !y %in% c("-", "N")
      
      dist_mat[i, j] <- sum(x[valid] != y[valid])
    }
  }
  
  dist_mat
}

infer_candidate_taxon_clusters <- function(group_df) {
  
  if (nrow(group_df) == 1) {
    return(
      group_df %>%
        dplyr::mutate(
          candidate_taxon_cluster = paste0(manual_group, "_cluster_1"),
          cluster_divergence_status = "single_OTU_cluster",
          within_cluster_max_diff = 0,
          nearest_other_cluster_min_diff = NA_real_,
          cluster_decision_basis = "only_one_OTU_in_manual_group"
        )
    )
  }
  
  dist_mat <- make_group_distance_matrix(group_df)
  
  hc <- hclust(
    as.dist(dist_mat),
    method = "average"
  )
  
  cluster_id <- cutree(
    hc,
    h = same_taxon_max_diff
  )
  
  cluster_df <- tibble::tibble(
    seq_id = rownames(dist_mat),
    candidate_taxon_cluster = paste0(
      unique(group_df$manual_group),
      "_cluster_",
      as.integer(cluster_id[rownames(dist_mat)])
    )
  )
  
  group_df2 <- group_df %>%
    dplyr::left_join(cluster_df, by = "seq_id")
  
  cluster_stats <- group_df2 %>%
    dplyr::group_by(candidate_taxon_cluster) %>%
    dplyr::summarise(
      cluster_seq_ids = list(seq_id),
      n_otus_in_candidate_cluster = dplyr::n(),
      cluster_total_reads = sum(SUM, na.rm = TRUE),
      cluster_group_RRA = sum(group_RRA, na.rm = TRUE),
      cluster_blast_summary = paste(unique(blast), collapse = " | "),
      cluster_interpretation_summary = paste(unique(interpretation), collapse = " | "),
      cluster_has_stop_codon = any(has_stop_codon %in% TRUE, na.rm = TRUE),
      cluster_has_middle_concentrated = any(middle_concentrated %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    )
  
  cluster_stats <- cluster_stats %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      within_cluster_max_diff = {
        ids <- unlist(cluster_seq_ids)
        if (length(ids) <= 1) {
          0
        } else {
          max(dist_mat[ids, ids], na.rm = TRUE)
        }
      },
      nearest_other_cluster_min_diff = {
        ids <- unlist(cluster_seq_ids)
        other_ids <- setdiff(rownames(dist_mat), ids)
        if (length(other_ids) == 0) {
          NA_real_
        } else {
          min(dist_mat[ids, other_ids], na.rm = TRUE)
        }
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      cluster_divergence_status = dplyr::case_when(
        is.na(nearest_other_cluster_min_diff) ~
          "single_candidate_taxon_within_manual_group",
        
        within_cluster_max_diff <= same_taxon_max_diff &
          nearest_other_cluster_min_diff >= distinct_taxon_min_diff ~
          "candidate_taxon_diverges_from_other_clusters",
        
        within_cluster_max_diff <= same_taxon_max_diff &
          nearest_other_cluster_min_diff < distinct_taxon_min_diff ~
          "candidate_taxon_close_to_another_cluster_review",
        
        within_cluster_max_diff > same_taxon_max_diff ~
          "cluster_internally_heterogeneous_review",
        
        TRUE ~
          "manual_review"
      ),
      
      cluster_decision_basis = paste0(
        "within_max_diff=", within_cluster_max_diff,
        "; nearest_other_min_diff=", nearest_other_cluster_min_diff,
        "; n_otus=", n_otus_in_candidate_cluster,
        "; cluster_reads=", cluster_total_reads,
        "; middle_pattern=", cluster_has_middle_concentrated,
        "; stop_codon=", cluster_has_stop_codon
      )
    ) %>%
    dplyr::select(-cluster_seq_ids)
  
  group_df2 %>%
    dplyr::left_join(
      cluster_stats,
      by = "candidate_taxon_cluster"
    )
}

selected_manual_groups_pairwise <- selected_manual_groups_pairwise %>%
  dplyr::group_split(manual_group) %>%
  lapply(infer_candidate_taxon_clusters) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(
    manual_group_order,
    manual_group,
    candidate_taxon_cluster,
    dplyr::desc(SUM)
  )

# Part 6. Merge pairwise and cluster results ==================================

new_pairwise_cluster_columns <- setdiff(
  names(selected_manual_groups_pairwise),
  names(global_curation_output)
)

pairwise_cluster_results <- selected_manual_groups_pairwise %>%
  dplyr::select(
    seq_id,
    dplyr::all_of(new_pairwise_cluster_columns)
  )

global_curation_output_final <- global_curation_output %>%
  dplyr::left_join(
    pairwise_cluster_results,
    by = "seq_id"
  )

global_curation_output_final <- global_curation_output_final %>%
  dplyr::mutate(
    final_merge_group = dplyr::case_when(
      curation_decision %in% c("keep_dominant", "merge_dominant") &
        !is.na(candidate_taxon_cluster) ~ candidate_taxon_cluster,

      curation_decision %in% c("keep_dominant", "merge_dominant") &
        manual_group != "not_manually_grouped" ~ manual_group,

      TRUE ~ seq_id
    ),

    final_merge_flag = dplyr::case_when(
      curation_decision == "merge_dominant" ~ TRUE,
      TRUE ~ FALSE
    ),

    final_keep_representative = dplyr::case_when(
      curation_decision == "keep_dominant" ~ TRUE,
      curation_decision != "merge_dominant" ~ TRUE,
      TRUE ~ FALSE
    ),

    final_merge_status = dplyr::case_when(
      curation_decision == "merge_dominant" ~
        "merge_into_representative",

      curation_decision == "keep_dominant" ~
        "representative_of_merge_group",

      curation_decision == "candidate_distinct_taxon" ~
        "keep_distinct_candidate_taxon",

      curation_decision == "to_blast_revalidate_dominant_with_internal_stop_codon" ~
        "keep_pending_blast_revalidation",

      curation_decision == "to_blast_validate_internal_stop_codon" ~
        "keep_pending_blast_validation",

      curation_decision == "outside_manual_group" ~
        "keep_outside_manual_group",

      TRUE ~
        "keep_review"
    )
  )

write.table(
  global_curation_output_final,
  "global_curation_output_final_water_per_subsample_upd.txt",
  sep = "\t",
  dec = ",",
  quote = FALSE,
  row.names = FALSE
)

# Part 7. Presence/absence matrix after merging ===============================

global_curation_output_final <- read.table(
  "global_curation_output_final_water_per_subsample_upd.txt",
  header = TRUE,
  sep = "\t",
  dec = ",",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata_cols <- names(global_curation_output_final)[1:12]

sample_cols_final <- names(global_curation_output_final)[
  (which(names(global_curation_output_final) == "pident") + 1):
    (which(names(global_curation_output_final) == "SUM") - 1)
]

global_curation_output_final <- global_curation_output_final %>%
  dplyr::mutate(
    merge_action = dplyr::case_when(
      curation_decision == "merge_dominant" ~ "merge",
      TRUE ~ "keep"
    ),

    merge_group = dplyr::case_when(
      curation_decision %in% c("keep_dominant", "merge_dominant") &
        !is.na(candidate_taxon_cluster) ~ candidate_taxon_cluster,

      curation_decision %in% c("keep_dominant", "merge_dominant") &
        manual_group != "not_manually_grouped" ~ manual_group,

      TRUE ~ seq_id
    )
  )

presence_absence_matrix <- global_curation_output_final %>%
  dplyr::filter(curation_decision != "exclude_artefact") %>%
  dplyr::group_by(merge_group) %>%
  dplyr::arrange(
    dplyr::desc(curation_decision == "keep_dominant"),
    dplyr::desc(interpretation == "dominant_reference"),
    dplyr::desc(SUM),
    .by_group = TRUE
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(metadata_cols),
      dplyr::first
    ),
    dplyr::across(
      dplyr::all_of(sample_cols_final),
      ~ sum(as.numeric(.x), na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sample_cols_final),
      ~ dplyr::if_else(.x > 0, 1, 0)
    )
  ) %>%
  dplyr::select(
    dplyr::all_of(metadata_cols),
    dplyr::all_of(sample_cols_final)
  )

write.table(
  presence_absence_matrix,
  "presence_absence_matrix_curated_water_per_subsample_upd.txt",
  sep = "\t",
  dec = ",",
  quote = FALSE,
  row.names = FALSE
)

