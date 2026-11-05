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
# Script:  01_control_subtraction.R
# Purpose: Negative-control subtraction and PCR-replicate merging for the 137
#          fish OTUs (step 1).
# Input:   SuppInfo_S1.xlsx, sheet 1_reads_raw_PCR_replicates
# Output:  01_reads_after_control_subtraction.csv (per PCR replicate)
#          01_reads_per_water_sample.csv (identical to sheet
#          2_reads_after_controls)
# =============================================================================

# 1. Data ---------------------------------------------------------------------

workbook <- "SuppInfo_S1.xlsx"
raw <- as.data.frame(readxl::read_excel(workbook, sheet = "1_reads_raw_PCR_replicates"))

meta_cols <- names(raw)[1:12]  # qseqid ... length
read_cols <- names(raw)[-(1:12)]

# Negative controls: field (*CN*) and extraction/PCR blanks (CN#, NTC*)
control_cols <- read_cols[grepl("CN", read_cols) | grepl("^NTC", read_cols)]
sample_cols  <- setdiff(read_cols, control_cols)
cat("PCR replicates:", length(sample_cols),
    "| negative-control replicates:", length(control_cols), "\n")

reads <- as.matrix(raw[, read_cols])
reads[is.na(reads)] <- 0

# 2. Subtract the maximum control count per OTU -------------------------------

max_control <- apply(reads[, control_cols, drop = FALSE], 1, max)
cleaned <- reads[, sample_cols, drop = FALSE] - max_control  # row-wise
cleaned[cleaned < 0] <- 0

write.csv(data.frame(raw[, meta_cols], cleaned, check.names = FALSE),
          "01_reads_after_control_subtraction.csv", row.names = FALSE)

# 3. Sum PCR replicates per water sample --------------------------------------

sample_id <- sub("_[A-F]$", "", sample_cols)   # drop replicate suffix
sample_id <- sub("^X24REW_", "", sample_id)    # prefix of the 2024 run

ids <- sort(unique(sample_id))
merged <- sapply(ids, function(id)
  rowSums(cleaned[, sample_id == id, drop = FALSE]))

write.csv(data.frame(raw[, meta_cols], merged, check.names = FALSE),
          "01_reads_per_water_sample.csv", row.names = FALSE)

# 4. Summary ------------------------------------------------------------------

cat("Reads before subtraction:", sum(reads[, sample_cols]), "\n")
cat("Reads after subtraction:", sum(merged), "\n")
cat("Water samples:", length(ids), "\n")
