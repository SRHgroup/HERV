Human endogenous retroviruses form a reservoir of T cell targets in
hematological cancers
================
Annie Borch
Sepetmber 23, 2020

  - [Project Overview](#project-overview)
      - [Abstract](#abstract)
  - [Description of folder contens](#description-of-folder-contens)
      - [Data](#data)
      - [Bin](#bin)
      - [Results](#results)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Project Overview

This directory contains data, scripts and results made in R from the
T-cell recognition and RNA expression analysis for this publication:
Human endogenous retroviruses form a reservoir of T cell targets in
hematological cancers If you would like to explore the scripts and data
from the Bayesian analysis see this page:
<https://github.com/SRHgroup/HERV_Bayesian_analyses>

## Abstract

Human endogenous retroviruses (HERV) form a substantial part of the
human genome, but mostly remain transcriptionally silent under strict
epigenetic regulation, yet can potentially be reactivated by malignant
transformation or epigenetic therapies. Here we evaluate the potential
for T cell recognition of HERV elements in myeloid malignancies by
mapping transcribed HERV genes and generating a library of 1,169
potential antigenic HERV-derived peptides predicted for presentation by
4 HLA class I molecules. Using DNA barcode-labeled MHC-I multimers, we
find CD8+ T cell populations recognizing 29 HERV-derived peptides
representing 18 different HERV loci, of which HERVH-5, HERVW-1, and
HERVE-3 have more profound responses; such HERV-specific T cells are
present in 17 of the 34 patients, but less frequently in healthy donors.
Transcriptomic analyses reveal enhanced transcription of the HERVs in
patients; meanwhile DNA-demethylating therapy causes a small and
heterogeneous enhancement in HERV transcription without altering T cell
recognition. Our study thus uncovers T cell recognition of HERVs in
myeloid malignancies, thereby implicating HERVs as potential targets for
immunotherapeutic therapies.

# Description of folder contens

## Data

### Raw data

#### Expression data

Kallisto is used to find expression from the raw RNA sequencing files.
Relevant gene sets are filtered from the raw files and saved in these
Rdata files. \#\#\#\# Clinical data All information about samples and
cell cycle

### Plot data

#### Barcode data

Output from the barcode screenings where the results from barracoda
(<https://services.healthtech.dtu.dk/service.php?Barracoda-1.8>) is
merge into one Rdata file.

#### Expression data

All extracted genes made from 1.Handle\_expression\_data.R

## Bin

### 1.Handle\_expressio\_data.R

This script will extract different gene signature and gene sets from the
RNA data.

### 2.HERV\_PaperPlots\_FinalRivision.R

This script is used for all plots according to the RNA expression
profile and T-cell recognition analysis.

## Results

### Paper Plots

All plots from made with R from the T-cell recognition analysis and the
RNA expression profile analysis.

### Tables

Tables including extra information (not used in the paper).

### Contact informations

  - Annie Borch <annbor@dtu.dk>
  - Anne-Mette Bjerregaard <anne.mette.bjerregaard@gmail.com>
  - Sunil Kumar Saini <sukusa@dtu.dk>
  - Sine Reker Hadrup <sirha@dtu.dk>
