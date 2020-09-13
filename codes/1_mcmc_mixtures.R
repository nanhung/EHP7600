# Hint: Select all R script and run them ####
# ** CAUTION: Will take about 1.5-2 hours**

# Load packages
library(rstan)
library(dplyr)

if(!dir.exists("outputs")) dir.create("outputs")
if(!dir.exists("reports")) dir.create("reports")

# Location of model and data files 
hill_two <- "codes/hill_two.stan"
mixture_data <- "datasets/mixture_data.csv"

# Load data
x <- read.csv(mixture_data)
mixtures <- c("AC50-L", "AC50-H", "Expo-L", "Expo-H", 
              "POD-L", "POD-H", "RFD-L", "RFD-H")      # Mixture name 
celltypes <- unique(x$celltype) %>% as.character()     # Cell type

# Create report
opt <- "reports/report_mcmc_mixtures.log"
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
  mix.list <- list("mix")
  
  for (j in 1:n.response){ # Corresponding responses
    
    for (k in 1:8){ # 8 mixtures
      
      # Data wrangling
      X <- x %>% 
        filter(celltype == unique(x$celltype)[i]) %>%    
        filter(phenotype == resp[j]) %>%  
        filter(mixture == mixtures[k])
      
      # Prepared setting for simulation
      dat <- list(
        "len" = nrow(X),  
        "y" = X$Response, 
        "d" = X$Dilution,
        "p_b" = c(0.00001, 100),
        "p_n" = c(0.1, 15.0),
        "p_sig" = c(0.0, 0.1) 
      )
      
      # Run!
      fit <- stan(file = hill_two, data = dat, 
                  iter = 8000, chains = chains, seed = 8,  
                  control = list(adapt_delta = 0.95, max_treedepth = 20))
      
      cat("\n Cell type: ", i, "/", 5, "; ", celltypes[i], "\n", 
          "Response: ", j, "/", n.response, "; ", resp[j], "\n",
          "Mix type: ", k, "/", 8, "; ", mixtures[k], "\n\n")
      
      # Print ongoing process in console
      print(fit) 
      
      # Send Stand output to file
      sink(opt, append=TRUE)  
      cat("\n Cell type: ", i, "/", 5, "; ", celltypes[i], "\n", 
          "Response: ", j, "/", n.response, "; ", resp[j], "\n",
          "Mix type: ", k, "/", 8, "; ", mixtures[k], "\n\n")
      print(fit) 
      cat("---------------------------------------------------------------------------\n\n")
      sink()
      
      # Extract and save result
      mcmc <- fit %>% rstan::extract()
      m <- cbind(mcmc$b, mcmc$n)
      if (k==1) M <- tail(m, 500) else M <- cbind(M, tail(m, 500)) # only keep 500 iterations
      
    } # end mixtures
    mix.list[[j]] <- M
    names(mix.list)[j] <- resp[j]
  }   # end responses 
  save(mix.list, file = paste0("outputs/", celltypes[i], "_mixtures_parms.rda"))
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
