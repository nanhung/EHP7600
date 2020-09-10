# Hint: Select all R script (Ctrl + A) and run them (Ctrl + Enter) 
# The file can run whole process and reproduce all results 
# CAUTION: Total time spent will be 7- 8 hours

library(rmarkdown)

Str.time <- Sys.time()

cores <- 1  
# Originally, we used 3 MCMC chains, but can change to 1 to save computational time 

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

#render('manuscript.Rmd')
#Sys.time()
#Sys.time() - Str.time 

render('supplementary.Rmd')
Sys.time()
Sys.time() - Str.time

cat(paste0("\nStarting time: ", Str.time), "\n\n")
cat(paste0("\nEnding time: ", Sys.time(), "\n\n"))
