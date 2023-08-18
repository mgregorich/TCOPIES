#additional system packages required on std ubuntu 22.04 installation
#sudo apt -y install cmake libcurl4-openssl-dev libfontconfig1-dev libxml2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libcairo2-dev libmagick++-dev

install.packages(c(
  "dplyr",
  "reshape2",
  "tidyr",
  "vegan",
  "ggplot2",
  "RColorBrewer",
  "stringdist",
  "tidyverse",
  "igraph",
  "sjPlot",
  "ineq",
  "BiodiversityR",
  "vegan",
  "SPECIES",
  "foreach",
  "doParallel",
  "doSNOW",
  "Matrix",
  "openxlsx",
  "pacman",
  "devtools",
  "sem",
  "rgl",
  "multcomp",
  "markdown",
  "lmtest",
  "leaps",
  "aplpack",
  "knitr",
  "flextable",
  "officer",
  "summarytools",
  "tableone",
  "stringr",
  "ggrepel",
  "gridExtra",
  "VennDiagram",
  "grid",
  "future.apply",
  "circlize",
  "RCircos",
  "ggpubr",
  "MASS",
  "questionr",
  "pander",
  "officedown"
  )
)

#install following packages no longer available from CRAN from archive
devtools::install_version("HDMD", version = "1.2")
devtools::install_version("tcR", version = "2.3.2")