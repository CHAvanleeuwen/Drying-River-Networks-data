**Drying River Networks Data**

This repository contains data-processing scripts and analyses associated with the Drying River Networks (DRYvER) project. The scripts are used to process aquatic invertebrate length and biomass data and to calculate secondary productivity across drying river networks.

**Overview**

Drying river networks experience changes in hydrological connectivity and environmental conditions as streams and rivers dry and rewet. Understanding how these changes affect aquatic communities and ecosystem processes requires comparable measurements across different river networks and sampling campaigns.

This repository contains the R scripts used to process biomass data from six drying river networks and to calculate secondary productivity from individual biomass, abundance, environmental, and campaign data.

**Study Areas**
The repository contains biomass-processing scripts for the following countries:

🇭🇷 Croatia
🇨🇿 Czech Republic
🇫🇮 Finland
🇫🇷 France
🇭🇺 Hungary
🇪🇸 Spain


**Repository Contents**

File	Description
Script biomass Croatia.R	Processing of biomass and length data for Croatia
Script biomass CzechR.R	Processing of biomass and length data for the Czech Republic
Script biomass Finland.R	Processing of biomass and length data for Finland
Script biomass France.R	Processing of biomass and length data for France
Script biomass Hungary.R	Processing of biomass and length data for Hungary
Script biomass Spain.R	Processing of biomass and length data for Spain
Secondary productivity calculations.R	Combines biomass, abundance and environmental data and calculates secondary productivity

The country-specific scripts follow a similar workflow for processing individual organism length measurements and estimating biomass. The secondary productivity script combines the resulting datasets across the six drying river networks.

Data Processing Workflow

The general workflow is:

Raw length measurements
        ↓
Country-specific biomass processing
        ↓
Individual length and biomass data
        ↓
Combine data from six drying river networks
        ↓
Merge with abundance data
        ↓
Merge with environmental data
        ↓
Estimate/interpolate missing length information
        ↓
Calculate secondary productivity

The secondary productivity workflow uses data from Croatia, France, Hungary, Finland, Spain, and the Czech Republic. It also incorporates environmental data and abundance information and calculates secondary productivity following Morin & Dumont (1994).

**Software**

The analyses are written in R.

The scripts use packages including:

ggplot2
dplyr
reshape2
readxl
vegan
GGally
leaflet
pals
stringr
lubridate

Install the required packages in R with:

install.packages(c(
  "ggplot2",
  "dplyr",
  "reshape2",
  "readxl",
  "vegan",
  "GGally",
  "leaflet",
  "pals",
  "stringr",
  "lubridate"
))

**Secondary Productivity**

The Secondary productivity calculations.R script:

Loads individual length and biomass data from the six drying river networks.
Combines the country-level datasets.
Removes incomplete records and selected records according to the analysis workflow.
Incorporates environmental data.
Calculates mean organism lengths and biomass values by taxonomic group, campaign, country, and stream type.
Incorporates abundance data.
Interpolates missing length information where required.
Calculates secondary productivity.

The secondary productivity calculations are based on Morin & Dumont (1994).


Repository: CHAvanleeuwen/Drying-River-Networks-data

