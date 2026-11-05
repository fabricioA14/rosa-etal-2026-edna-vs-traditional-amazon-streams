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
# Script:  02_otu_prescreening_functions.R
# Purpose: Helper functions for the BOLD/GenBank validation in
#          02_otu_prescreening.R.
# Input:   Sourced by 02_otu_prescreening.R
# Output:  None (functions only)
# =============================================================================

# 1. GenBank (NCBI BLAST) top hit ---------------------------------------------
blast_genbank_top_taxon <- function(sequence) {
  
  sequence <- toupper(sequence)
  sequence <- gsub("[^ACGTN]", "", sequence)
  
  if (is.na(sequence) || nchar(sequence) < 50) {
    return(paste0("GenBank:sequence_too_short_", nchar(sequence), "bp"))
  }
  
  put_res <- httr::POST(
    url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
    body = list(
      CMD = "Put",
      PROGRAM = "blastn",
      DATABASE = "nt",
      QUERY = sequence,
      HITLIST_SIZE = 1,
      MEGABLAST = "on",
      EXPECT = 10
    ),
    encode = "form"
  )
  
  put_txt <- httr::content(put_res, as = "text", encoding = "UTF-8")
  
  rid <- stringr::str_match(put_txt, "RID = ([A-Z0-9\\-]+)")[, 2]
  rtoe <- stringr::str_match(put_txt, "RTOE = ([0-9]+)")[, 2]
  
  if (is.na(rid)) {
    return("GenBank:no_RID_returned")
  }
  
  wait_time <- ifelse(is.na(rtoe), 20, as.numeric(rtoe) + 5)
  Sys.sleep(wait_time)
  
  max_tries <- 20
  
  for (j in seq_len(max_tries)) {
    
    get_res <- httr::GET(
      url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
      query = list(
        CMD = "Get",
        RID = rid,
        FORMAT_TYPE = "XML"
      )
    )
    
    get_txt <- httr::content(get_res, as = "text", encoding = "UTF-8")
    
    if (grepl("Status=WAITING", get_txt)) {
      Sys.sleep(15)
      next
    }
    
    if (grepl("Status=FAILED", get_txt)) {
      return("GenBank:failed")
    }
    
    if (grepl("Status=UNKNOWN", get_txt)) {
      return("GenBank:unknown_RID")
    }
    
    if (grepl("<BlastOutput", get_txt)) {
      
      xml <- xml2::read_xml(get_txt)
      
      hit_def <- xml2::xml_text(
        xml2::xml_find_first(xml, ".//Hit_def")
      )
      
      taxon <- extract_specific_taxon(hit_def)
      
      if (is.na(taxon) || taxon == "") {
        return("GenBank:no_hit")
      }
      
      return(paste0("GenBank:", taxon))
    }
    
    Sys.sleep(15)
  }
  
  "GenBank:timeout"
}

# 2. Taxon name from a GenBank hit description --------------------------------
extract_specific_taxon <- function(hit_def) {
  
  if (is.na(hit_def) || hit_def == "") {
    return(NA_character_)
  }
  
  hit_def <- stringr::str_squish(hit_def)
  
  # Remove common GenBank strings
  hit_def <- gsub("(?i),?\\s*complete\\s+mitochondrial\\s+genome.*$", "", hit_def, perl = TRUE)
  hit_def <- gsub("(?i),?\\s*mitochondrion,?\\s*complete\\s+genome.*$", "", hit_def, perl = TRUE)
  hit_def <- gsub("(?i),?\\s*mitochondrial\\s+DNA.*$", "", hit_def, perl = TRUE)
  hit_def <- gsub("(?i),?\\s*cytochrome\\s+c\\s+oxidase.*$", "", hit_def, perl = TRUE)
  hit_def <- gsub("(?i),?\\s*cytochrome\\s+oxidase.*$", "", hit_def, perl = TRUE)
  hit_def <- gsub("(?i),?\\s*COI.*$", "", hit_def, perl = TRUE)
  hit_def <- gsub("(?i),?\\s*COX1.*$", "", hit_def, perl = TRUE)
  hit_def <- gsub("(?i),?\\s*voucher.*$", "", hit_def, perl = TRUE)
  hit_def <- gsub("(?i),?\\s*isolate.*$", "", hit_def, perl = TRUE)
  hit_def <- gsub("(?i),?\\s*specimen.*$", "", hit_def, perl = TRUE)
  
  hit_def <- stringr::str_squish(hit_def)
  
  # Capture binomial: Genus species
  binomial <- stringr::str_extract(
    hit_def,
    "\\b[A-Z][a-z]+\\s+[a-z][a-z\\-]+\\b"
  )
  
  if (!is.na(binomial) && binomial != "") {
    return(binomial)
  }
  
  # Capture Genus sp.
  genus_sp <- stringr::str_extract(
    hit_def,
    "\\b[A-Z][a-z]+\\s+sp\\.?"
  )
  
  if (!is.na(genus_sp) && genus_sp != "") {
    return(genus_sp)
  }
  
  # Capture genus only
  genus <- stringr::str_extract(
    hit_def,
    "\\b[A-Z][a-z]+\\b"
  )
  
  if (!is.na(genus) && genus != "") {
    return(genus)
  }
  
  NA_character_
}

# 3. BOLD ID Engine query (chromote) ------------------------------------------
run_bold_id_engine <- function(
  b,
  fasta_path,
  url = "https://id.boldsystems.org/",
  wait_page = 5,
  wait_result = 10,
  max_attempts = 30
) {
  if (!file.exists(fasta_path)) {
    stop("FASTA file not found: ", fasta_path)
  }
  
  # 1. Open start page
  b$Page$navigate(url)
  Sys.sleep(wait_page)
  
  current_url <- b$Runtime$evaluate("document.location.href")$result$value
  message("Current URL: ", current_url)
  
  # 2. Check that the FASTA input exists
  input_check <- b$Runtime$evaluate("
    JSON.stringify({
      exists: document.getElementById('fasta-file') !== null,
      input_count: document.querySelectorAll('input').length,
      file_input_count: document.querySelectorAll('input[type=file]').length
    })
  ")$result$value
  
  message("Input check: ", input_check)
  
  has_fasta_input <- b$Runtime$evaluate("
    document.getElementById('fasta-file') !== null
  ")$result$value
  
  if (!isTRUE(has_fasta_input)) {
    stop("Could not find #fasta-file. The page may not be on the upload form.")
  }
  
  # 3. Get the actual input node
  doc <- b$DOM$getDocument()
  
  file_node <- b$DOM$querySelector(
    nodeId = doc$root$nodeId,
    selector = "#fasta-file"
  )
  
  if (is.null(file_node$nodeId) || file_node$nodeId == 0) {
    stop("DOM node for #fasta-file was not found.")
  }
  
  # 4. Insert FASTA into the input
  b$DOM$setFileInputFiles(
    nodeId = file_node$nodeId,
    files = list(normalizePath(fasta_path))
  )
  
  upload_check <- b$Runtime$evaluate("
    JSON.stringify({
      exists: document.getElementById('fasta-file') !== null,
      files_length: document.getElementById('fasta-file').files.length,
      file_name: document.getElementById('fasta-file').files.length > 0 ?
        document.getElementById('fasta-file').files[0].name : null
    })
  ")$result$value
  
  message("Upload check: ", upload_check)
  
  uploaded <- b$Runtime$evaluate("
    document.getElementById('fasta-file').files.length === 1
  ")$result$value
  
  if (!isTRUE(uploaded)) {
    stop("FASTA was not attached correctly.")
  }
  
  # 5. Click Identify
  clicked <- b$Runtime$evaluate("
    const btn = [...document.querySelectorAll('button,input[type=submit]')]
      .find(x =>
        x.innerText.trim() === 'Identify' ||
        x.value === 'Identify'
      );
    
    if (btn) {
      btn.click();
      true;
    } else {
      false;
    }
  ")$result$value
  
  if (!isTRUE(clicked)) {
    stop("Could not find/click the Identify button.")
  }
  
  # 6. Wait for results to load
  result_text <- NULL
  result_url <- NULL
  
  for (i in seq_len(max_attempts)) {
    Sys.sleep(wait_result)
    
    result_url <- b$Runtime$evaluate("document.location.href")$result$value
    result_text <- b$Runtime$evaluate("document.body.innerText")$result$value
    
    message("\n--- attempt ", i, " ---")
    message("URL: ", result_url)
    message(substr(result_text, 1, 500))
    
    if (
      grepl("Identification Engine - Result", result_text, ignore.case = TRUE) &&
      grepl("Classifications", result_text, ignore.case = TRUE)
    ) {
      message("BOLD ID results loaded.")
      break
    }
  }
  
  if (
    is.null(result_text) ||
    !grepl("Identification Engine - Result", result_text, ignore.case = TRUE)
  ) {
    warning("Result page was not clearly detected.")
  }
  
  return(list(
    result_url = result_url,
    result_text = result_text
  ))
}



# 4. Parse BOLD results -------------------------------------------------------
parse_bold_hits_v2 <- function(txt) {
  
  combined_block <- txt %>%
    str_split("Combined Hits", simplify = TRUE) %>%
    .[, 2] %>%
    str_split("Showing 1 to", simplify = TRUE) %>%
    .[, 1]
  
  lines <- combined_block %>%
    str_split("\n", simplify = FALSE) %>%
    unlist() %>%
    str_trim()
  
  lines <- lines[lines != ""]
  
  otu_idx <- which(str_detect(lines, "^Otu"))
  
  if (length(otu_idx) == 0) {
    return(tibble())
  }
  
  out <- lapply(seq_along(otu_idx), function(i) {
    
    start <- otu_idx[i]
    end <- if (i < length(otu_idx)) otu_idx[i + 1] - 1 else length(lines)
    
    block <- lines[start:end]
    
    query_id <- block[1]
    pid <- block[2]
    bin <- str_extract(block[3], "BOLD:[A-Z0-9]+")
    
    tax_line <- block[4]
    
    parts <- str_split(tax_line, "\\s+", simplify = TRUE)
    parts <- as.character(parts)
    
    id_percent <- parts[length(parts)]
    indels <- parts[length(parts) - 1]
    
    phylum <- parts[1]
    class <- parts[2]
    order <- parts[3]
    family <- parts[4]
    
    genus <- parts[length(parts) - 3]
    species <- paste(parts[(length(parts) - 3):(length(parts) - 2)], collapse = " ")
    
    tibble(
      query_id = query_id,
      pid = pid,
      bin = bin,
      phylum = phylum,
      class = class,
      order = order,
      family = family,
      subfamily = NA_character_,
      tribe = NA_character_,
      genus = genus,
      species = species,
      indels = as.integer(indels),
      id_percent = as.numeric(id_percent)
    )
  })
  
  bind_rows(out)
}

# 5. Best BOLD hit ------------------------------------------------------------
get_simple_bold_hit <- function(bold_hits_one) {
  
  if (is.null(bold_hits_one) || nrow(bold_hits_one) == 0) {
    return("BOLD:no_hit")
  }
  
  best <- bold_hits_one %>%
    dplyr::arrange(dplyr::desc(id_percent), indels) %>%
    dplyr::slice(1)
  
  if (
    !is.na(best$species[1]) &&
    best$species[1] != "" &&
    !stringr::str_detect(
      best$species[1],
      "\\bsp\\.?$|\\bsp\\.\\b|unidentified|uncultured"
    )
  ) {
    return(paste0("BOLD:", best$species[1]))
  }
  
  if (!is.na(best$genus[1]) && best$genus[1] != "") {
    return(paste0("BOLD:", best$genus[1]))
  }
  
  if (!is.na(best$family[1]) && best$family[1] != "") {
    return(paste0("BOLD:", best$family[1]))
  }
  
  if (!is.na(best$order[1]) && best$order[1] != "") {
    return(paste0("BOLD:", best$order[1]))
  }
  
  if (!is.na(best$class[1]) && best$class[1] != "") {
    return(paste0("BOLD:", best$class[1]))
  }
  
  if (!is.na(best$phylum[1]) && best$phylum[1] != "") {
    return(paste0("BOLD:", best$phylum[1]))
  }
  
  "BOLD:no_hit"
}

