# Hint: Select all R script and run them ####
# ** CAUTION: Will take about 4.5 hours**

# Load package ------------------------------------------------------------
library(dplyr)
library(reshape2)


# Define the variables ----------------------------------------------------
# Will be used in following analysis
celltypes <- c("iCell Neurons", "HUVECs", "iCell Hepatocytes", "iCell Endothelial cells",  "iCell Cardiomyocytes")
mixture_info <- "datasets/mixture_info.csv"
mix.name <- read.csv(mixture_info)[2:9] %>% names()    # Mixture's name in data variable
mixtures <- c("AC50-L", "AC50-H", "Expo-L", "Expo-H", 
              "POD-L", "POD-H", "RFD-L", "RFD-H")      # mixture's label 

if(exists("DF_CF")) rm("DF_CF")
if(exists("DF_CA")) rm("DF_CA")
if(exists("DF_IA")) rm("DF_IA")

# Exposure-response -------------------------------------------------------
# l: Mixtures (Total: 8)
# i: Celltype (Total: 5)
# j: Phenotype in each cell (Neurons: 10; HUVECs: 8; Hepatocytes: 5; Endothelial cells: 8; Cardiomyocytes: 16)
# k: Individual chemicals (Total: 42)
strt<-Sys.time()

# mixture loop
for (l in 1:8){
  
  conc <- read.csv(mixture_info) %>% select(mix.name[l])
  freq <- read.csv(mixture_info)[2:9] %>% melt() %>% 
    filter(variable == mix.name[l]) %>%
    mutate(freq = value/sum(value)) %>% select(freq)
  
  # Celltype loop
  for (i in 1:5){
    
    load(file = paste0("outputs/", celltypes[i], "_mixtures_parms.rda")) 
    load(file = paste0("outputs/", celltypes[i], "_42_chem_parms.rda"))
    
    # Phenotype loop
    for (j in 1:length(names(mix.list))){
      
      resp <- names(mix.list)[j]
      M <- mix.list[[j]] # select parameter matrix for specific response from mixtures
      
      d <- 1/10^seq(-1, 5, 0.1)
      parms.df <- M[,c(2*l-1, 2*l)] %>% as.data.frame()
      for(k in seq(length(d))){ # use for plot
        parms.df[,2+k] <- 1 / ((1 + d[k] / parms.df[,1] )^parms.df[,2])
      }
      df_CF <- apply(parms.df[,c(3:(length(d)+2))], 2, quantile, c(0.5, 0.025, 0.975)) %>% 
        t %>% `colnames<-`(c("y.med", "y.lcl95", "y.ucl95")) %>% as.data.frame() 
      df_CF$x.conc <- d * sum(conc)
      df_CF$celltype <- celltypes[i]
      df_CF$mixture <- mixtures[l]
      df_CF$phenotype <- resp
      rownames(df_CF) <- c()
      
      M <- ind.list[[j]] # select parameter matrix for specific response from individual chemical
      
      # Independent action
      for(iter in 1:500){
        param <- matrix(M[iter,], nrow = 42, byrow = T) %>% 
          as.data.frame() %>% as.matrix() # to be remove
        
        for(k in seq(length(d))){ 
          Conc <- conc*d[k]
          effect <- prod(( 1 / (1 + (Conc /  param[,1] )^ param[,2])))
          if(k==1) E <- effect else E <- c(E, effect)
        } 
        if(iter==1) total.E <- E else total.E <- c(total.E, E)
      }
      IA.matrix <- total.E %>% matrix(ncol = 61, byrow = T)
      
      df_IA <- apply(IA.matrix, 2, quantile, c(0.5, 0.025, 0.975)) %>% t %>%
        `colnames<-`(c("y.med", "y.lcl95", "y.ucl95")) %>% as.data.frame() 
      df_IA$x.conc <- d * sum(conc)
      df_IA$celltype <- celltypes[i]
      df_IA$mixture <- mixtures[l]
      df_IA$phenotype <- resp  

      
      # Concentration addition
      effv <- c(seq(0.01, 0.99, 0.01), 0.999)
      for (EC in seq(length(effv))){
        for (k in 1:42){
          p <- freq[k,]
          ED <- M[,k*2-1] * ((1 / effv[EC])^(1 / M[,k*2])-1)
          
          denominator <- p/ED
          if (k == 1) {
            Denom <- denominator
          } else{
            Denom <- cbind(Denom, denominator)
          }
        } 
        ED_mix <- 1 / apply(Denom, 1, sum)
        if (EC == 1) {
          ED_MIX <- ED_mix
        } else{
          ED_MIX <- cbind(ED_MIX, ED_mix)
        }
      } 
      df_CA <- apply(ED_MIX, 2, quantile,  c(0.5, 0.025, 0.975)) %>% t() %>%
        `colnames<-`(c("x.med", "x.lcl95", "x.ucl95")) %>% as.data.frame()
      df_CA$y.resp <- rev(1-effv)
      df_CA$celltype <- celltypes[i]
      df_CA$mixture <- mixtures[l]
      df_CA$phenotype <- resp
      rownames(df_CA) <- c()
      
      
      if(!exists("DF_CF")){
        DF_CF <- df_CF
      } else {
        DF_CF <- rbind(DF_CF, df_CF)
      }
      if(!exists("DF_CA")){
        DF_CA <- df_CA
      } else {
        DF_CA <- rbind(DF_CA, df_CA)
      }
      if(!exists("DF_IA")){
        DF_IA <- df_IA
      } else {
        DF_IA <- rbind(DF_IA, df_IA)
      }
      
      cat("\n Mixture: ", l, "/", 8,  "; ", mixtures[l], "\n",
          "Celltype: ", i, "/", 5,   "; ", celltypes[i], "\n", 
          "Response: ", j, "/", length(names(mix.list)), "; ", resp, "\n\n")
      print(Sys.time()-strt)
    } # End phenotype loop
  }   # End celltype loop
}     # End mixture loop

# Save outputs ------------------------------------------------------------
mix_CR.list <- list("mix_CR")
mix_CR.list[[1]] <- DF_CF
mix_CR.list[[2]] <- DF_CA
mix_CR.list[[3]] <- DF_IA
names(mix_CR.list) <- c("CF", "CA", "IA")
save(mix_CR.list, file = paste0("outputs/mix_CR.rda"))
cat(paste0("Starting time: ", strt), "\n")
cat(paste0("Ending time: ", Sys.time()))
