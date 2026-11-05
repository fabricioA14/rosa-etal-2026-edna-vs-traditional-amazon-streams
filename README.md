<div align="center">

# Comparative Performance of Traditional Sampling and Water Environmental DNA for Detecting Fish Diversity in Amazonian Streams Across Hydrological Periods

**Fabricio dos Anjos Santa Rosa**<sup>1,2&#42;</sup>, **Derlan José Ferreira Silva**<sup>1&#42;</sup>, Cintia Oliveira Carvalho<sup>1,3</sup>, Silvia Britto Barreto<sup>2</sup>, Silvana Melo Sviggum<sup>4</sup>, Anna Karolina Oliveira de Queiroz<sup>1,2</sup>, Enki Gautier<sup>1</sup>, Marcia da Silva Anjos<sup>1</sup>, Jessica Conceição Nascimento Dergan<sup>1</sup>, Orlando Gabriel da Silva Solano Reis<sup>1</sup>, Tibério Cesar Tortola Burlamaqui<sup>5</sup>, Birgitte Lisbeth Graae Thorbek<sup>6</sup>, João Bráullio de Luna Sales<sup>1</sup>, Quentin Mauvisseau<sup>6</sup>, Hugo de Boer<sup>6</sup>, Jonathan Stuart Ready<sup>1</sup>

</div>

<sup>1</sup>Group for Integrated Biological Investigation, Center for Advanced Studies of Biodiversity, Federal University of Pará. Av. Perimetral da Ciência, km 01, PCT‐Guamá, Lot 11, Belém-PA. ZIP 66075‐750, Brazil<br>
<sup>2</sup>Vale Institute of Technology Sustainable Development (ITV DS), Belém-PA, Brazil<br>
<sup>3</sup>The School of Business, Innovation and Sustainability (FIH), Halmstad University, Halmstad, Sweden<br>
<sup>4</sup>Oslo University Hospital, Oslo, Norway<br>
<sup>5</sup>Applied Genetics Laboratory, Federal Rural University of the Amazon, Av. Tancredo Neves, Terra Firme, n 2501, Belém-PA. ZIP 66077-830<br>
<sup>6</sup>Natural History Museum of Oslo, University of Oslo, Norway

<sup>&#42;</sup>Joint first authors

Corresponding author: Fabricio dos Anjos Santa Rosa (fabriciodosanjos14@gmail.com)

## Study

Fish assemblages of 12 first- to third-order streams of the Potiritá River sub-basin (Capim River basin, eastern
Amazonia) were sampled with seine nets and hand sieves and with water eDNA (COI metabarcoding; nine 400 mL water
samples per stream) in four campaigns: August 2022 and August 2023 (dry season), February 2023 and February 2024
(rainy season).

## Contents

| File | Description |
|---|---|
| `01_control_subtraction.R` | Step 1. Negative-control subtraction and merging of PCR replicates (137 fish OTUs). |
| `02_otu_prescreening.R` | Step 2. Automated pre-screening of OTUs (stop codons, abundance and divergence metrics, BOLD/GenBank validation). |
| `02_otu_prescreening_functions.R` | Helper functions for the BOLD/GenBank validation (sourced by step 2). |
| `03_analyses.R` | Step 3. Taxon counts, heatmap by season and campaign (Figure 2), nMDS (Figure 3), ANOSIM/PERMANOVA/betadisper (Figure 4) and SIMPER (Figure 5). |

Curation path of the eDNA data: 137 OTUs → 122 (redundant OTUs merged) → 105 (technical artefacts removed) → 82 taxa
(further merges and one OTU below 92% identity excluded). Each step is a sheet of `SuppInfo_S1.xlsx`.

## How to run

This repository holds only the scripts. The dataset they use, `SuppInfo_S1.xlsx` (Supplementary Information S1 of
the article; README in its first sheet), will be made available upon acceptance of the associated article. To run the
scripts, place that file in this folder. Then run each script from this folder, in order:

```
Rscript 01_control_subtraction.R
Rscript 02_otu_prescreening.R
Rscript 03_analyses.R
```

`03_analyses.R` reads only the sheet `3_analysis_input` and can be run on its own. The order illustrations of
Figure 2 are optional: if a folder `fish_icons/` with one PNG per order (e.g., `Characiformes.png`) is present they are
drawn; otherwise the figure is produced without them. The BOLD/GenBank step of
`02_otu_prescreening.R` needs internet access and a Chrome browser (package `chromote`); its final, manually reviewed
decisions are already stored in the sheet `2_reads_curation_decisions`.

**Software:** R 4.6.1; vegan 2.7-5, ggplot2, ggh4x, ggrepel, patchwork, gtable, png, dplyr, tidyr, stringr, purrr,
readxl, Biostrings, DECIPHER, readr, tibble, httr, xml2, chromote. Permutation tests use a fixed seed
(`set.seed(123)`), so the statistics reported in the manuscript are reproduced exactly.

## Data availability

The dataset is not included in this repository and will be made available upon acceptance of the associated article:
raw sequencing data on Zenodo (https://doi.org/10.5281/zenodo.20977161) and processed data (read counts at each step,
curation decisions and final datasets) as Supplementary Information S1 of the article (`SuppInfo_S1.xlsx`). The reference database was built with refDBdelimiter
(https://github.com/fabricioA14/refDBdelimiter).

## License

The scripts are released under the MIT License (also in `LICENSE`). The license covers only the scripts; the data are
not part of this repository (see Data availability).

```
MIT License

Copyright (c) 2026 Rosa et al.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

"Rosa et al." refers to all authors listed above.

## Citation

Please cite the article (reference to be added upon publication) and the Zenodo dataset.
