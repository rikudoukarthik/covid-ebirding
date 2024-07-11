# Effects of COVID-19 pandemic on eBirding and birds in India

*Papers: [![DOI](https://zenodo.org/badge/DOI/10.1093/ornithapp/duae024.svg)](https://doi.org/10.1093/ornithapp/duae024)*

This repository contains the source code and data associated with a study on the effects of the COVID-19 pandemic on eBirding ([Thrikkadeeri & Viswanathan 2024](https://doi.org/10.1093/ornithapp/duae024)) and bird species reporting (Thrikkadeeri & Viswanathan, in prep.) in India. All analyses were performed in the R environment in RStudio, and with the addition of a couple of large (publicly available) data files that need to be procured externally (see below), the source code and files in this repo will allow the analyses to be fully reproduced.

The structure and workflow for the main analysis is outlined below.

## Data

**External**

The analysis is centred around eBird data from India, which is publicly available but is a very large file.

-   [*eBird Basic Dataset (EBD)*](https://ebird.org/data/download/ebd)*, relMay-2022*: The current version of this dataset is much more recent and contains data after May 2022, but filtering it should produce a fairly similar dataset to the ones used in this study, with only minor changes.
-   *eBird sensitive species data, relMay-2022*: Data on [sensitive species data](https://ebird.org/india/news/ebird-sensitive-species/) is not included in the public download, and was obtained after requesting separately from eBird. (Only required for bird reporting analyses, i.e., manuscript #2.)
-   Spatial polygon data of admin. unit boundaries

**Reproducible**

-   [*MODIS Land Cover Type Product MCD12Q1*](https://lpdaac.usgs.gov/products/mcd12q1v006/): The data used in the study is available in `00_data/in_LULC_MODIS/` & `00_data/rast_UNU.RData`, but it is recommended to delete those files before starting the analysis. This will then force `getmodisdata()` from `00_scripts/functions.R` to run and produce fresh LULC data files from MODIS.
-   *Timeline of COVID classification*: Available `00_data/covid_classification.csv`. Classification of the timeline of interest into various COVID categories.

### Reproduce results of Thrikkadeeri & Viswanathan 2024

Only the following data files, which are filtered, processed and prepared to be input into the analyses, are required to reproduce our results. These are mostly `.Rdata` files that can be readily imported into R (all analyses were performed in the R environment using RStudio).

-   The processed version of EBD ready for analyses is `00_data/data0_MY_d_slice.RData`
-   Shapefiles: `00_data/maps_sf.RData` contains admin unit polygons, and `00_data/grids_st_sf.RData` contains country-wide grids of multiple resolutions
-   The processed MODIS land cover data is available as `00_data/rast_UNU.RData`, but deleting this before starting the analysis will force fresh LULC data files from MODIS to be produced from the script
-   Classification of the timeline of interest into various COVID categories, available as `00_data/covid_classification.csv`

## Workflow

### Data preparation

`03_wrap.Rmd` is the place to start for the main analyses. This loads all necessary data (and processes it where required), packages and other 00_scripts/functions. It is essentially a setup file.

Lines 80--92 should be unhashed when running the analysis for the first time, as this will activate `data_qualfilt_prep()` from `00_scripts/functions.R`. In subsequent runs, the lines can be hashed again, since various `.RData` and other data files will have been generated in the original data filtering, processing and preparation step.

### Analysis

Then, the analyses can be run. The study is divided into two broad sections, one dealing with changes in birding and birder behaviour, and the other with changes in bird species reporting. Hence, the analysis is also split into two files, `04_birder.Rmd` and `05_bird.Rmd`. The contents of both these files are self-explanatory; they run various analyses, and generate and store outputs of the analyses such as figures (`03_wrap_figs/`) and summary data (`00_outputs/`), which can be directly utilised further downstream.

### Manuscript

Once the full analysis has been successfully run (and outputs generated and stored), the final manuscript can also be reproduced (and validated), since it is an `.Rmd` file and most summaries and results in it (textual or graphical) are obtained directly from the data and analysis outputs.

The full manuscript is divided into the main manuscript (`11_manuscript.Rmd`) and the supplementary material (`12_manuscript_supp.Rmd`). `index.Rmd` uses [bookdown](https://bookdown.org/) (see `_bookdown.yml`) to combine the two and produce the final document (`_book/_main.docx`). This building of the book is executed by `bookdown::render_book()`. This step **requires a template `.docx` file**, containing all necessary formatting rules and styles, to be present in the repository---in this case, named `rmd_word_template_journ_guide.docx`, which needs to be created separately. The bibliography for this study (`covid-ebirding.json`) and citation style language (CSL) for the journal of interest (`ibis.csl`) are also necessary.

A separate title page containing author affiliations, acknowledgements, etc. can be knitted directly from `13_titlepage.Rmd`, without having to build a book.

Some aspects or elements of the final manuscript are **non-reproducible**. They are detailed below.

**Model summary tables**

Due to the difficulty in converting the raw summary objects of LMMs/GLMMs into informative summary tables programmatically, the current heuristic is to take fixed effects information directly from the model summary output (`00_outputs/`) into a simple and neat table using `knitr::kable()`, while the random effects information as well as other elements of the table caption are altered manually each time. See last section of `12_manuscript_supp.Rmd` for details.

**Main text and supplementary material**

The reason for using `bookdown::render_book()` (instead of simply knitting the two documents separately) is to allow easy cross-referencing of sections, figures, tables, etc. between the main text and the supplementary material. However, one disadvantage is that [there is no easy way](https://stackoverflow.com/questions/50223141/rmarkdown-bookdown-separate-figure-numbering-for-supplemental-section) to restart numbering of sections, figures, tables, etc. from a given section onwards or to add a prefix to the numbering.

Therefore, the manual steps followed are:

1.  Create two copies of the final document (`_book/_main.docx`). In the first, manually delete all text from the supplementary section onwards. In the second, delete all text until the supplementary section.
2.  Adjust the numbering of supplementary sections.
    -   In the supplementary document, *find and replace* "7." with "". This removes the unnecessary tier of numbering.
    -   In the main document, *find and replace* "7." with "S". This adds the prefix "S" to all in-text references to supplementary sections.
3.  Manually adjust the numbering of *x* supplementary figures and tables one-by-one (for *x* from 1 to 46).
    -   In the supplementary document, *find and replace* "Figure x" or "Fig. x" with "Figure x-8" or "Fig. x-8" respectively (because there 8 figures in the main document). Similarly, *find and replace* "Table x" with "Table x-2" (because there 2 tables in the main document).
    -   In the main document, *find and replace* "Figure x" or "Fig. x" with "Figure Sx-8" or "Fig. Sx-8" respectively, and "Table x" with "Table Sx-2". This adds the prefix "S" to all in-text references to supplementary sections.

**Separate figure legends**

The journal requires all individual figure legends in the main manuscript to be listed together at the end of the document. For this, simply copy-paste the 8 figure legends one-by-one in the last section.

**Separate Conflict of Interest statement**

The journal requires the CoI statement to be provided separately as well. The steps are:

1.  Create a copy of the title page document generated from `13_titlepage.Rmd`.
2.  Retain only the title and CoI section.
3.  Add a closing section to sign off the document, including a digital version of the corresponding author's signature.
