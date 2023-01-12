# Please ensure to read the README file before reproducing the analyses.
# This script is designed to set the environment by clearing the current session, setting the seed for reproducibility and installing packages that are not already installed. Then it loads all packages in the list, to be used in the script.

# Clear the environment
rm(list = ls())

# Set seed for reproducibility
set.seed(2022)

# Define list of packages to be used in script
packages <- 
  list(
    'gsheet', # package to interact with Google Sheets
    'urltools', # package for URL handling
    'RCurl', # package for connecting to URLs
    'stringr', # package for string manipulation
    'readxl', # package for reading excel files
    'stringr', # package for string manipulation
    'dplyr', # package for data manipulation
    'tidyr', # package for data tidying
    'parsedate', # package for date parsing
    'gt', # package for creating tables
    'gtsummary', # package for creating summary tables
    'gridExtra', # package for creating grid layouts
    'ggpubr', # package for creating ggplot2 based publication ready plots
    'glue', # package for string interpolation
    'lubridate', # package for working with dates
    'knitr', # package for creating dynamic documents
    'hrbrthemes', # package for creating ggplot2 themes
    'ggsci', # package for creating ggplot2 scientific plots
    'ggthemes', # package for creating ggplot2 themes
    'naniar', # package for handling missing data
    'flextable', # package for creating flexible tables
    'rlist', # package for handling lists
    'cowplot', # package for creating complex plots
    'psrwe', # package for creating progress bars
    'arm', # package for creating multivariate regression models
    'NPP', # package for nonparametric statistical methods
    'reshape2', # package for reshaping data
    'progress', # package for creating progress bars
    'logitnorm', # package for creating logit-normal models
    'gtools', # package for various tools
    'gifski', # package for creating gifs
    'psrwe', # package for creating progress bars
    'parallel', # package for parallel processing
    'progress' # package for creating progress bars
  )

# Check if packages are not installed, install if not
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  lapply(packages[!installed_packages], install.packages)
  install.packages(packages[!installed_packages])
}

# Load all packages in the list
invisible(lapply(packages, library, character.only = TRUE))
