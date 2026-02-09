# ZackPhen

This repository contains materials used to support the analyses for the
paper titled: “Inferences about phenological shifts in an Arctic
community vary with time-windows”. This study aims to examine how
signals of phenological change are manifested at different temporal
scales using data collected through the [BioBasis
program](https://data.g-e-m.dk/datasets) at the Zackenberg Valley from
1996 to present.

## Goals:

-   assess how trend directionality changes with time-series duration
    and choice of start years;

-   identify the minimum time-series length needed to obtain a high
    probability of agreement between shorter- and long-term trends.

## Folders:

-   data: contains raw and cleaned data used in the analysis

-   RScripts: contains analysis scripts for both plants and arthropods

-   diagnostics: contains taxon-specific files of model diagnostics

-   models: contains taxon-specific model objects

-   predictions: contains taxon-specific figures visualizing model
    parameters and predictions

-   Misc: contains exploratory analysis materials

## Relevant Scripts:

To reproduce this analysis, the following scripts in the RScripts folder
may be useful:

### Data manipulation

arthropod\_data\_cleaning.R: clean arthropod data, remove year-plot
combinations that do not form a bell curve

covariates\_data\_cleaning.R: clean air temperature and spring snow
cover data

plant\_data\_cleaning.R: clean plant data, remove year-plot combinations
that do not form a bell curve

### Analysis

#### Stan files

arth:phen2.stan: final Stan code for arthropod data

plant\_phen6b.stan: final Stan code for plant data

#### Model fitting files

arthropod\_functions.R: source code for fitting and examining arthropod
models

plant\_functions.R: source code for fitting and examining plant models

arthropod\_stanmod\_PD.R: code to run for interfacing with Stan in R,
and perform post-processing steps for arthropod models

plant\_stanmod\_PD.R: code to run for interfacing with Stan in R, and
perform post-processing steps for plant models

phen\_slide.R: script to extract slopes at different time-windows for
both plants and arthropods

signal\_strength.R: script to examine trend directionality and strength
for plants and arthropods

For inquiries, you may send an email to:
[patdumandan@gmail.com](patdumandan@gmail.com).
