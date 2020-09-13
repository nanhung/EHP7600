# Hint: Select all R script and run them ####
# ** CAUTION: Will take about 2 - 2.5 hours**

# Load packages
library(rstan)
library(dplyr)

if(!dir.exists("outputs")) dir.create("outputs")
if(!dir.exists("reports")) dir.create("reports")

# Location of model and data files 
hill_two <- "codes/hill_two.stan"
chem_data <- "datasets/chem_data.csv"

# Load data
x <- read.csv(chem_data)
chems <- unique(x$chemical) %>% as.character()     # Chemical name 
celltypes <- unique(x$celltype) %>% as.character() # cell type

# Create report
opt <- "reports/report_mcmc_ind_chems.log"
cat("", file = opt) 
sink(opt, append=TRUE)  
cat(paste0("Starting time: ", Sys.time()), "\n")
sink()

# Fit individual C-R  -----------------------------------------------------
Str.t <- Sys.time()

# Gloal setting of RStan
if (exists("cores")) cores <- cores else cores <- 3
rstan_options(auto_write = TRUE)
options(mc.cores = cores)

for(i in 1:5){ # 5 Human induced pluripotent stem cell (iPSC)-derived cells
  
  # Check number of phenotypes and its names
  resp <- x %>% filter(celltype == unique(x$celltype)[i]) %>%
    select(phenotype) %>% unique() %>% unlist() %>% as.character()
  n.response <- length(resp)
  ind.list <- list("ind")
  
  for (j in 1:n.response){ # Corresponding responses
    
    for (k in 1:42){ # 42 individual chemicals
      
      # Data wrangling
      X <- x %>% 
        filter(celltype == unique(x$celltype)[i]) %>%    
        filter(phenotype == resp[j]) %>%  
        filter(chemical == chems[k])
      X <- X[complete.cases(X),]
      
      # Prepared setting for simulation
      dat <- list(
        "len" = nrow(X),                     # No. of data points
        "y" = X$Response,               # Corresponding response
        "d" = X$Concentration,          # Given Concentration
        "p_b" = c(0.001, 10000),        # Prior for ec50
        "p_n" = c(1, 15),               # Prior for Hill coefficient
        "p_sig" = c(0.0, 0.1)           # Prior for error
      )
      
      # Run!
      fit <- stan(file = hill_two, data = dat, 
                  iter = 4000, chains = chains, seed = 42,  
                  control = list(adapt_delta = 0.95, max_treedepth = 20))
      
      cat("\n Cell type: ", i, "/", 5, "; ", celltypes[i], "\n", 
          "Response: ", j, "/", n.response, "; ", resp[j], "\n",
          "Chemical: ", k, "/", 42, "; ", chems[k], "\n\n")
      
      # Print ongoing process in console
      print(fit) 
      
      # Send Stand output to file
      sink(opt, append=TRUE)  
      cat("\n Cell type: ", i, "/", 5, "; ", celltypes[i], "\n", 
          "Response: ", j, "/", n.response, "; ", resp[j], "\n",
          "Chemical: ", k, "/", 42, "; ", chems[k], "\n\n")
      print(fit) 
      cat("---------------------------------------------------------------------------\n\n")
      sink()

      # Extract and save result
      mcmc <- fit %>% rstan::extract()
      m <- cbind(mcmc$b, mcmc$n)
      if (k==1) M <- tail(m, 500) else M <- cbind(M, tail(m, 500)) # only keep 500 iterations

    } # end 42 individual chemicals
    ind.list[[j]]<-M
    names(ind.list)[j] <- resp[j]
  }   # end responses 
  save(ind.list, file = paste0("outputs/", celltypes[i], "_42_chem_parms.rda"))
}     # end 5 iPSC-derived cells

# Report time spent (about 3.6 hours) 
message(paste0("Ending time: ", Sys.time()))
end.t <- Sys.time() - Str.t 
end.t
sink(opt, append=TRUE)  
end.t
cat(paste0("\nEnding time: ", Sys.time(), "\n\n"))
print(sessionInfo())
sink()
