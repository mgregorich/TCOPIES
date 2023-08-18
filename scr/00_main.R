###################################
# Author: MG
# Date:10/03/2020
# Info: 
###################################

setwd("scr")

rm(list=ls())

library(HDMD)             #Contains Atchley factor data
library(dplyr)
library(reshape2) 
library(tidyr)            #Packages for data cleaning and shaping
library(vegan)            #Contains similarity functions
library(ggplot2)          #Plotting package
library(RColorBrewer)    #Heatmap packages
library(stringdist)       # Calculate Levenshtein distance
library(tidyverse)
library(igraph)
library(sjPlot)
library(ineq)             # Gini index
library(tcR)              # Advanced Data Analysis of Immune Receptor Repertoires
library(BiodiversityR)
library(vegan)            #Contains similarity functions
library(SPECIES)
library(foreach)
library(doParallel)
library(doSNOW)
library(Matrix)           # Sparse matrix representation
library(openxlsx)

pacman::p_load(openxlsx)

CONVENTIONAL_ONLY <- T

# Data cleaning parameter
ambiguity.ratio=2
fold.change=NULL
frequency=NULL
nclones = 10000
layout="kk"


# Define paths
path <- "../Daten/All_corr/"
clean.path <- paste0(path, "clean_tcopies_fold", fold.change,"_freq",frequency,"_ar",ambiguity.ratio,"/")
path <- clean.path

foldername_sim <- paste0("sim_run_n",nclones,"_",layout)
out.path <- "../Output/"
out.path.final = file.path(out.path, "_final") #required by report rmd
if(dir.exists( paste0(out.path,foldername_sim,"/"))){
  do.call(file.remove, list(list.files(paste0(out.path,foldername_sim,"/"), full.names = TRUE)))
}else{
  dir.create(paste(out.path,foldername_sim, sep=""), recursive = T)
}
out.path <- paste0(out.path,foldername_sim,"/")

if(!dir.exists(out.path.final)){
  dir.create(out.path.final, recursive = T)
}

if(!dir.exists("../Daten/All_rerun_corr/clean_tcopies_fold5_freq1e-10_ar2")){ #ensure path used in data cleaning exists
  dir.create("../Daten/All_rerun_corr/clean_tcopies_fold5_freq1e-10_ar2", recursive = T)
}

# ------- Execute code

source("01_data_cleaning.R")

if(!CONVENTIONAL_ONLY){
  source("02_network_analysis.R")
}

rmarkdown::render("report_TCOPIES_publ.Rmd", envir = new.env(), params = list(conventional_only = CONVENTIONAL_ONLY))
