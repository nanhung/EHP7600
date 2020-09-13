# Instruction ------------------------------------------------------------------
# Set the working directory in R
#> setwd("C:/Users/nan_1/Documents/EHP7600")
# Select all R script (Ctrl + A) and run them (Ctrl + Enter)

# The file can run whole process and reproduce all results 
# CAUTION: Total time spent will be 9- 10 hours
# Lastest update: 9/12/20

# Install packages--------------------------------------------------------------
# Suggest installing all rstan dependencies
# install.packages("rstan", dependencies = T)

# The packages use in this study
# pkg <- c("gridGraphics", "knitr", "jpeg", "tiff", "tidyverse", 
#         "reshape2","scales", "ggpubr", "treemapify", "cowplot", "bayestestR",
#         "ggridges", "gridExtra", "PerformanceAnalytics", "rmarkdown")
# install.packages(pkg)

# RUN --------------------------------------------------------------------------
Str.time <- Sys.time()
cores <- 1 # use only 1 core to prevent the crash
chains <- 3

source("codes/0_mcmc_ind_chems.R")
Sys.time()
Sys.time() - Str.time 

source("codes/1_mcmc_mixtures.R")
Sys.time()
Sys.time() - Str.time 

source("codes/2_ec10_pred.R")
Sys.time()
Sys.time() - Str.time 

source("codes/3_conc_resp_pred.R")
Sys.time()
Sys.time() - Str.time 

source("codes/4_plot.R")
Sys.time()
Sys.time() - Str.time 

source("codes/5_plot_suppl.R")
Sys.time()
Sys.time() - Str.time 

# Create document (need to install pandoc)--------------------------------------
#rmarkdown::render('manuscript.Rmd')
#Sys.time()
#Sys.time() - Str.time 

rmarkdown::render('supplementary.Rmd')
Sys.time()
Sys.time() - Str.time

cat(paste0("\nStarting time: ", Str.time), "\n")
cat(paste0("\nEnding time:   ", Sys.time(), "\n"))
