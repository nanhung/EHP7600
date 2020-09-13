# Hint: Select all R script and run them ####
# ** CAUTION: Will take about 30 minutes**

# Load package ------------------------------------------------------------
library(dplyr)
library(reshape2)

# Define the variables ----------------------------------------------------
celltypes <- c("iCell Neurons", "HUVECs", "iCell Hepatocytes", "iCell Endothelial cells",  "iCell Cardiomyocytes")
mixture_info <- "datasets/mixture_info.csv"
mix.name <- read.csv(mixture_info)[2:9] %>% names()    # Mixture's name in data variable
mixtures <- c("AC50-L", "AC50-H", "Expo-L", "Expo-H", 
              "POD-L", "POD-H", "RFD-L", "RFD-H")      # mixture's label
#x <- read.csv(mixture_info)
#cum_conc_info <- read.csv(mixture_info) %>% 
#  `colnames<-`(c("Chemical", mixtures, "Class")) %>%
#  melt() %>% 
#  group_by(variable) %>% 
#  summarize(CumConc = sum(value))


# EC10 estimation ---------------------------------------------------------
# l: Mixtures (Total: 8)
# i: Celltype (Total: 5)
# j: Phenotype in each cell (Neurons: 10; HUVECs: 8; Hepatocytes: 5; Endothelial cells: 8; Cardiomyocytes: 16)
# k: Individual chemicals (Total: 42)
ec10.list <- list("ec10")
strt <- Sys.time()
# mixture loop
for (l in 1:8){ 
  conc <- read.csv(mixture_info) %>% select(mix.name[l])
  freq <- read.csv(mixture_info)[2:9] %>% melt() %>% 
    filter(variable == mix.name[l]) %>%
    mutate(freq = value/sum(value)) %>% select(freq)
  
  if(exists("ec10_n")) rm(ec10_n)
  
  # Celltype loop
  for (i in 1:5){
    
    load(file = paste0("outputs/", celltypes[i], "_mixtures_parms.rda")) 
    load(file = paste0("outputs/", celltypes[i], "_42_chem_parms.rda"))
    
    # Phenotype loop
    for (j in 1:length(names(mix.list))){
      
      resp <- names(mix.list)[j]
      M <- mix.list[[j]] # select parameter matrix for specific response from mixtures
      
      ec10_fit <- M[,l*2-1] * ((1/0.9)^(1/M[,l*2])-1) * sum(conc)
      
      M <- ind.list[[j]] # select parameter matrix for specific response from individual chemical
      
      # Independent action (refer to mixtox package)
      effPoints <- 0.9   # ec10
      if(exists("ec10_ia")) rm(ec10_ia)
      for(iter in 1:500){
        param <- matrix(M[iter,], nrow = 42, byrow = T) %>% 
          as.data.frame() %>% as.matrix()
        
        iaFun <- as.character(1)
        for (k in 1:42){
          iaFun <- paste0(iaFun, '*', '( 1 / (1 + (xx*', conc[k,],'/', param[k, 1], ')^', param[k, 2], '))')
        }
        fun <- paste(effPoints, '-',  iaFun, sep = '')
        f = function(xx) eval(parse(text = fun))
        root <- uniroot(f, c(0.001, 20),  extendInt = "yes", tol = 1e-10)$root
        if (!exists("ec10_ia")) ec10_ia <- root*sum(conc) else ec10_ia <- c(ec10_ia, root*sum(conc))
      }
      
      # Concentration addition
      for (k in 1:42){
        p <- freq[k,]
        ec10 <- M[,k*2-1] * ((1/effPoints)^(1/M[,k*2])-1)
        denominator10 <- p/ec10
        if (k == 1) {
          Denom10 <- denominator10
        } else{
          Denom10 <- cbind(Denom10, denominator10)
        }
      }
      ec10_ca <- 1 / apply(Denom10, 1, sum) 
      
      ec10_i <- data.frame(rep(c("CF", "CA", "IA"), each=500), c(ec10_fit, ec10_ca, ec10_ia), resp, celltypes[i]) %>% 
        `colnames<-`(c("Method", "EC10", "Phenotype", "Celltype"))
      if(!exists("ec10_n")){ 
        ec10_n <- ec10_i
      } else ec10_n <- rbind(ec10_n, ec10_i) 
      
      cat("\n Mixture: ", l, "/", 8,  "; ", mixtures[l], "\n",
          "Celltype: ", i, "/", 5,   "; ", celltypes[i], "\n", 
          "Response: ", j, "/", length(names(mix.list)), "; ", resp, "\n\n")
      print(Sys.time()-strt)
    }  # End phenotype loop
  }    # End celltype loop
  ec10.list[[l]] <- ec10_n
  names(ec10.list)[l] <- mixtures[l]
}      # End mixture loop

# Save outputs ------------------------------------------------------------
save(ec10.list, file = paste0("outputs/ec10_pred.rda"))
cat(paste0("Starting time: ", strt), "\n")
cat(paste0("Ending time: ", Sys.time()))
