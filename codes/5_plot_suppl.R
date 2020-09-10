# Hint: Select all R script and run them ####
# ** CAUTION: Will take about 7 minutes**

library(dplyr)
library(reshape2)
library(bayestestR)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)

if(!dir.exists("plots/suppl")) dir.create("plots/suppl")

Str.t <- Sys.time()
# Plot C-R plot for individual chemicals
# cell: Celltype (Total: 5)
# resp: Phenotype in each cell (Neurons: 10; HUVECs: 8; Hepatocytes: 5; Endothelial cells: 8; Cardiomyocytes: 16)
plot.indCR <- function(cell = 1, resp = 1){
  
  # Load data
  chem_data <- "datasets/chem_data.csv"
  celltypes <- c("iCell Neurons", "HUVECs", "iCell Hepatocytes", "iCell Endothelial cells", "iCell Cardiomyocytes") 
  load(file = paste0("outputs/", celltypes[cell], "_42_chem_parms.rda"))
  n.response <- length(ind.list)
  name.resp <- names(ind.list)  
  m <- ind.list[[resp]]
  path <- paste0("Datasets/ind_chems/")
  x <- read.csv(chem_data) %>% filter(celltype == celltypes[cell] & phenotype == name.resp[resp])
  chem <- x$chemical %>% unique() %>% as.character()
  
  par(mfrow = c(7,6), mar = c(2,2,2,1), oma = c(2,2,2,0))
  for (k in 1:42){ # 42 individual chemicals
    
    # Extract data from specific chemical
    dat <- x %>% filter(chemical == chem[k])
    dat$conc <- rep(10^(c(-2:2)), each=2)
    
    # Estimate ec50 and its high density interval
    iter <- 100    # only pick 200 iters to plot
    m <- tail(m, iter) 
    ec50 <- m[,2*k-1] * ((1/0.9)^(1/m[,2*k])-1) 
    hdi50 <- hdi(ec50, ci =0.9) 
    
    # Generate corresponding concentration-response
    str <- 2*k-1
    M_2 <- m[,c(str,str+1)] %>% as.data.frame()
    d <- 100/10^seq(-0.5, 4.5, 0.1)
    for(l in seq(length(d))){ M_2[,2+l] <- 1 / ((1 + d[l] / M_2[,1] )^M_2[,2]) }
    
    # Plot
    plot(d, M_2[1, c(3:(length(d)+2))], log = "x", type="l", col="grey", ylim = c(0, 1.4), 
         xlab = "", ylab = "", main = chem[k], cex.main=0.8)
    for(iter in 2:iter){ lines(d, M_2[iter, c(3:(length(d)+2))], col="grey") }
    abline(0.5,0, lty=2, col= "grey20")
    lines(c(hdi50$CI_low, hdi50$CI_high), c(0.5,0.5), col = "red")
    abline(v=hdi50$CI_low, col="red", lty=2)
    abline(v=hdi50$CI_high, col="red", lty=2)
    points(dat$Concentration, dat$Response, col = "blue")
  } 
  main <- paste(celltypes[cell], "- ", name.resp[resp])
  mtext(main, NORTH<-3, line=0.5, outer=TRUE, cex = 1)
  mtext("Response", WEST<-2, line=0.5, outer=TRUE, cex = 1)
  mtext("Concentration, uM", SOUTH<-1, line=0.5, outer=TRUE, cex = 1)
}

# Plot C-R plot for mixtures
# cell: Celltype (Total: 5)
# resp: Phenotype in each cell (Neurons: 10; HUVECs: 8; Hepatocytes: 5; Endothelial cells: 8; Cardiomyocytes: 16)
plot.mixCR <- function(cell = 1, resp = 1){
  
  # Load data
  celltypes <- c("iCell Neurons", "HUVECs", "iCell Hepatocytes", "iCell Endothelial cells", "iCell Cardiomyocytes") 
  load(file = paste0("outputs/", celltypes[cell], "_mixtures_parms.rda"))
  n.response <- length(mix.list)
  name.resp <- names(mix.list)  
  m <- mix.list[[resp]]
  
  par(mfcol = c(2,4), mar = c(2,2,2,1), oma = c(2,2,2,0))
  mixture_info <- "datasets/mixture_info.csv"
  mixtype <- read.csv(mixture_info)[2:9] %>% names()
  
  mixtures <- c("AC50-L", "AC50-H", "Expo-L", "Expo-H", "POD-L", "POD-H", "RFD-L", "RFD-H")      
  cum_conc_info <- read.csv(mixture_info) %>% 
    `colnames<-`(c("Chemical", mixtures, "Class")) %>%
    melt() %>% 
    group_by(variable) %>% 
    summarize(CumConc = sum(value))
  
  mixture_data <- "datasets/mixture_data.csv"
  mix.dat <- read.csv(mixture_data)
  
  for (k in 1:8){      # 8 mixtures
    
    cum_conc <- as.numeric(cum_conc_info[k,2])
    # Extract data from specific mixture
    mix_conc <- mix.dat %>% 
      filter(celltype == celltypes[cell] & mixture == mixtures[k] & phenotype == name.resp[resp]) %>%
      mutate(Conc = Dilution * cum_conc) %>%
      select(Conc) %>% unlist()
    mix_resp <- mix.dat %>% 
      filter(celltype == celltypes[cell] & mixture == mixtures[k] & phenotype == name.resp[resp]) %>%
      select(Response) %>% unlist()
    
    # Estimate ec50 and its high density interval
    iter <- 100
    m <- tail(m, iter) # only pick 200 iters to plot
    ec10 <- m[,2*k-1] * ((1/0.9)^(1/m[,2*k])-1) * cum_conc
    
    hdi10 <- hdi(ec10, ci=0.9) 
    
    # Generate corresponding concentration-response
    str <- 2*k-1
    M_2 <- m[,c(str,str+1)] %>% as.data.frame()
    d <- 10^seq(-4.5, 0.5, 0.1)
    for(l in seq(length(d))){
      M_2[,2+l] <- 1 / ((1 + d[l] / M_2[,1] )^M_2[,2])
    }
    
    # Plot
    plot(d*cum_conc, M_2[1, c(3:(length(d)+2))], log = "x", type="l", col="grey", ylim = c(0, 1.4), 
         xlab = "", ylab = "", main = mixtype[k], cex.main=1)
    for(iter in 2:iter){
      lines(d*cum_conc, M_2[iter, c(3:(length(d)+2))], col="grey")
    }
    abline(0.9,0, lty=2, col= "grey20")
    lines(c(hdi10$CI_low, hdi10$CI_high), c(0.9,0.9), col = "red")
    abline(v=hdi10$CI_low, col="red", lty=2)
    abline(v=hdi10$CI_high, col="red", lty=2)
    points(mix_conc, mix_resp, col = "blue")
  } 
  main <- paste(celltypes[cell], "- ", name.resp[resp])
  mtext(main, NORTH<-3, line=0.5, outer=TRUE, cex = 1)
  mtext("Effect (propotion of control)", WEST<-2, line=0.5, outer=TRUE, cex = 1)
  mtext("Total concentration, uM", SOUTH<-1, line=0.5, outer=TRUE, cex = 1)
}

# Plot C-R profile to compare curve-fitting, independent action and concentration addition
mixplot <- function(cell = 1, pheno = 1, mix = "AC50-H"){
  
  load(file = "outputs/mix_CR.rda") 
  mixtures <- c("AC50-L", "AC50-H", "Expo-L", "Expo-H", 
                "POD-L", "POD-H", "RFD-L", "RFD-H")      # Mixture name 
  celltypes <- c("iCell Neurons", "HUVECs", "iCell Hepatocytes", "iCell Endothelial cells", "iCell Cardiomyocytes") 
  
  CumConc <- read.csv("datasets/mixture_info.csv") %>% 
    `colnames<-`(c("Chemical", mixtures, "Class")) %>%
    melt() %>% 
    group_by(variable) %>% 
    summarize(CumConc = sum(value)) %>%
    filter(variable == mix) %>% select(CumConc) %>% as.numeric()
  
  ct <- celltypes[cell]
  type <- mix_CR.list$CF %>% filter(mixture==mix & celltype == ct) %>% select(phenotype) %>% unique()
  CF <- mix_CR.list$CF %>% filter(mixture==mix & celltype == ct & phenotype == type[pheno,1])
  CA <- mix_CR.list$CA %>% filter(mixture==mix & celltype == ct & phenotype == type[pheno,1])
  IA <- mix_CR.list$IA %>% filter(mixture==mix & celltype == ct & phenotype == type[pheno,1])
  data <- read.csv("datasets/mixture_data.csv") %>% 
    filter(mixture==mix & celltype == ct & phenotype == type[pheno,1])
  
  plot(c(1e-1, 1e5), c(0,1.45), type = "n", log="x", 
       xaxt="n", 
       yaxt="n", 
       ylab="Response",
       xlab="Total concentration, uM")
  title(main = type[pheno,1], line = 1)
  
  at=seq(-1, 5, 2)
  at.minor <- log10(outer(1:9, 10^(min(at):max(at))))
  lab <- sapply(at, function(i) as.expression(bquote(10^ .(i))))
  axis(side=1, at=10^at, labels=lab)
  axis(side=2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.3, 1.4), 
       labels=c("0.0","0.2","0.4","0.6","0.8","1.0", "IA", "CA", "Fitting"), las = 2) 
  abline(h = 0.9, lty=3)
  lines(CA$x.med, CA$y.resp, col=rgb(0, 0, 1))
  lines(IA$x.conc, IA$y.med, col=rgb(0, 1, 0))
  lines(CF$x.conc, CF$y.med)
  polygon(c(rev(CA$x.lcl95), CA$x.ucl95), 
          c(rev(CA$y.resp), CA$y.resp), 
          col=rgb(0, 0, 1, 0.1), border = NA)
  polygon(c(rev(IA$x.conc), IA$x.conc), 
          c(rev(IA$y.lcl95), IA$y.ucl95), 
          col=rgb(0, 1, 0, 0.1), border = NA)
  polygon(c(rev(CF$x.conc), CF$x.conc), 
          c(rev(CF$y.lcl95), CF$y.ucl95), 
          col=rgb(0, 0, 0, 0.1), border = NA)
  points(data$Dilution * CumConc, data$Response, pch=19)
  
  load(file = "outputs/ec10_pred.rda")
  ec10.list[[mix]] %>% filter(Celltype == ct & Phenotype == type[pheno,1] & Method == "CF") %>% 
    select(EC10) %>% boxplot(horizontal=TRUE, add = T, axes=F, at = c(1.4), notch = F, 
                             outline = F, boxwex =0.15)
  ec10.list[[mix]] %>% filter(Celltype == ct & Phenotype == type[pheno,1] & Method == "CA") %>% 
    select(EC10) %>% boxplot(horizontal=TRUE, add = T, axes=F, at = c(1.3), notch = F, 
                             outline = F, boxwex =0.15, col = rgb(0, 0, 1, 0.1))
  ec10.list[[mix]] %>% filter(Celltype == ct & Phenotype == type[pheno,1] & Method == "IA") %>% 
    select(EC10) %>% boxplot(horizontal=TRUE, add = T, axes=F, at = c(1.2), notch = F, 
                             outline = F, boxwex =0.15, col = rgb(0, 1, 0, 0.1))
}

# The box plot that rank the individaul chemicals by margin of exposure
ind.MOE.plot <- function(cell = 1, pheno = 1, mix = "AC50-H", xlab) {
  
  load(file = "outputs/ec10_pred.rda") 
  mixtures <- c("AC50-L", "AC50-H", "Expo-L", "Expo-H", 
                "POD-L", "POD-H", "RFD-L", "RFD-H")      # Mixture name 
  celltypes <- c("iCell Neurons", "HUVECs", "iCell Hepatocytes", "iCell Endothelial cells", "iCell Cardiomyocytes") 
  
  indConc <- read.csv("datasets/mixture_info.csv") %>% 
    `colnames<-`(c("Chemical", mixtures, "Class")) %>% 
    select(mix)
  
  #
  chem_info <- read.csv("datasets/chem_info.csv")
  load(file = paste0("outputs/", celltypes[cell], "_42_chem_parms.rda")) # mixture results
  
  M <- ind.list[[pheno]]
  for(k in 1:42){
    EC10 <- M[,k*2-1] * ((1/0.9)^(1/M[,k*2])-1)
    y <- EC10/indConc[k,1]
    if (k==1) MOE <- y else MOE <- c(MOE, y)
  }
  df.MOE <- matrix(MOE, ncol=42) %>% as.data.frame() %>%`colnames<-`(chem_info$Chemical)
  rank.MOE <- apply(df.MOE, 2, median) %>% as.data.frame()
  rank.MOE$rank <- rank(rank.MOE$.)
  new.level <- rownames(rank.MOE[order(rank.MOE$rank),])
  df2.MOE <- melt(df.MOE) %>% `colnames<-`(c("Chemical", "MOE"))
  df2.MOE$Chemical <- factor(df2.MOE$Chemical, levels = new.level)  
  df2.MOE$Phenotype <- names(ind.list)[pheno]
  #
  load(file = "outputs/ec10_pred.rda")
  ct <- celltypes[cell]
  pt <- names(ind.list)[pheno]
  df.mixMOE <- ec10.list[[mix]] %>% filter(Method == "CF" & Celltype == ct & Phenotype == pt) %>% mutate(MOE = EC10/sum(indConc))
  med.CF.MOE <- df.mixMOE %>% select(MOE) %>% as.matrix() %>%  median()
  #
  level <- levels(df2.MOE[,1])
  ggplot(df2.MOE, aes(x=Chemical, y=MOE)) + 
    geom_boxplot(outlier.size = -1) +
    annotate("text", x = 1:42, y = 1e-12, label = level, hjust = 0, size = 2.5, col="grey20") +
    geom_hline(yintercept = 100, lty = 2, size=0.5, col="pink") + 
    geom_hline(yintercept = med.CF.MOE, lty = 2, size=0.5, col="red") + 
    scale_y_log10(lim = c(10^-12, 10^3),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    coord_flip() +
    theme_pubclean() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          panel.border=element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(colour = "black", size =0.5),
          title = element_text(size=8, face = "bold")) +
    labs(y = "", x = "", title = names(ind.list)[pheno])
}

# The box plot describe the margin of exposure for curve fitting, independent action, and concentration addition
mix.MOE.plot <- function(cell=1, pheno=1, mix = "AC50-H"){
  celltypes <- c("iCell Neurons", "HUVECs", "iCell Hepatocytes", "iCell Endothelial cells", "iCell Cardiomyocytes") 
  mixtures <- c("AC50-L", "AC50-H", "Expo-L", "Expo-H", "POD-L", "POD-H", "RFD-L", "RFD-H")      # Mixture name 
  load(file = paste0("outputs/", celltypes[cell], "_mixtures_parms.rda")) # mixture results
  load(file = "outputs/ec10_pred.rda")
  ct <- celltypes[cell]
  pt <- names(mix.list)[pheno]
  indConc <- read.csv("datasets/mixture_info.csv") %>% 
    `colnames<-`(c("Chemical", mixtures, "Class")) %>% 
    select(mix)
  df.mixMOE <- ec10.list[[mix]] %>% filter(Celltype == ct & Phenotype == pt) %>% 
    mutate(MOE = EC10/sum(indConc))
  df.mixMOE$Method <- factor(df.mixMOE$Method, levels = c("CF", "IA", "CA"))
  med.CF.MOE <- df.mixMOE %>% filter(Method == "CF") %>% select(MOE) %>% as.matrix() %>%  median()
  
  ggplot(df.mixMOE, aes(x=Method, y=MOE)) + 
    geom_boxplot(outlier.size = -1) +
    geom_hline(yintercept = 100, lty = 2, size=0.5, col="pink") + 
    geom_hline(yintercept = med.CF.MOE, lty = 2, size=0.5, col="red") + 
    scale_y_log10(lim = c(10^-12, 10^3),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    coord_flip() +
    theme_pubclean() +
    theme(axis.text.x=element_text(color = "black", face="bold", size =10),
          axis.text.y=element_blank(),
          panel.border=element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          axis.ticks.x = element_line(colour = "black", size =0.5)) +
    annotate("text", x = 1:3, y = 1e-12, 
             label = c("Fitting", "Independent action", "Concentration addition"),
             hjust = 0, size = 3) +
    labs(y = "Margin of Exposure", x ="")  
}


# Curve-fitting of single chemical concentration-respsonse ----------------

## iCell Neurons
for(i in 1:10){
  file = paste0("plots/suppl/FigS", i,".png")
  png(file = file, width=1800, height=1800, res=200)
  plot.indCR(1,i)
  dev.off()
  cat(paste("\n* Create", file, "\n\n"))
}
Sys.time() - Str.t 

## HUVECs
for(i in 1:8){
  file = paste0("plots/suppl/FigS", 10+i,".png")
  png(file = file, width=1800, height=1800, res=200)
  plot.indCR(2,i)
  dev.off()
  cat(paste("* Create", file, "\n\n"))
}
Sys.time() - Str.t 

## iCell Hepatocytes
for(i in 1:5){
  file = paste0("plots/suppl/FigS", 18+i,".png")
  png(file = file, width=1800, height=1800, res=200)
  plot.indCR(3,i)
  dev.off()
  cat(paste("* Create", file, "\n\n"))
}
Sys.time() - Str.t 

## iCell Endothelial cells
for(i in 1:8){
  file = paste0("plots/suppl/FigS", 23+i,".png")
  png(file = file, width=1800, height=1800, res=200)
  plot.indCR(4,i)
  dev.off()
  cat(paste("* Create", file, "\n\n"))
}
Sys.time() - Str.t 

## iCell Cardiomyocytes
for(i in 1:16){
  file = paste0("plots/suppl/FigS", 31+i,".png")
  png(file = file, width=1800, height=1800, res=200)
  plot.indCR(5,i)
  dev.off()
  cat(paste("* Create", file, "\n\n"))
}
Sys.time() - Str.t 

# Curve-fitting of mixture concentration-respsonse ------------------------

## iCell Neurons
for(i in 1:10){
  file = paste0("plots/suppl/FigS", 47+i,".png")
  png(file = file, width=2000, height=1000, res=200)
  plot.mixCR(1,i)  
  dev.off()
  cat(paste("* Create", file, "\n\n"))
}

## HUVECs
for(i in 1:8){
  file = paste0("plots/suppl/FigS", 57+i,".png")
  png(file = file, width=2000, height=1000, res=200)
  plot.mixCR(2,i)  
  dev.off()
  cat(paste("* Create", file, "\n\n"))
}

## iCell Hepatocytes
for(i in 1:5){
  file = paste0("plots/suppl/FigS", 65+i,".png")
  png(file = file, width=2000, height=1000, res=200)
  plot.mixCR(3,i)  
  dev.off()
  cat(paste("* Create", file, "\n\n"))
}

## iCell Endothelial cells
for(i in 1:8){
  file = paste0("plots/suppl/FigS", 70+i,".png")
  png(file = file, width=2000, height=1000, res=200)
  plot.mixCR(4,i)  
  dev.off()
  cat(paste("* Create", file, "\n\n"))
}

## iCell Cardiomyocytes
for(i in 1:16){
  file = paste0("plots/suppl/FigS", 78+i,".png")
  png(file = file, width=2000, height=1000, res=200)
  plot.mixCR(5,i)  
  dev.off()
  cat(paste("* Create", file, "\n\n"))
}
Sys.time() - Str.t 

# Curve-fitting and prediction of AC50-H concentration-response -----------

## iCell Neurons
file = paste0("plots/suppl/FigS95.png")
png(file = file, width=2200, height=1800, res=200)
par(mfrow = c(3,4))
for(i in 1:10){
  mixplot(1,i)
}
dev.off()

## iCell HUVECs
file = paste0("plots/suppl/FigS96.png")
png(file = file, width=2200, height=1200, res=200)
par(mfrow = c(2,4))
for(i in 1:8){
  mixplot(2,i)
}
dev.off()

## iCell Hepatocytes
file = paste0("plots/suppl/FigS97.png")
png(file = file, width=2200, height=1200, res=200)
par(mfrow = c(2,4))
for(i in 1:5){
  mixplot(3,i)
}
dev.off()

## iCell Endothelial cells
file = paste0("plots/suppl/FigS98.png")
png(file = file, width=2200, height=1200, res=200)
par(mfrow = c(2,4))
for(i in 1:8){
  mixplot(4,i)
}
dev.off()

## iCell Cardiomyocytes
file = paste0("plots/suppl/FigS99.png")
png(file = file, width=2200, height=2400, res=200)
par(mfrow = c(4,4))
for(i in 1:16){
  mixplot(5,i)
}
dev.off()
Sys.time() - Str.t 

# Estimation of the margin of exposure under AC50-H exposure --------------

## iCell Neurons
p1.1.1 <- ind.MOE.plot(cell = 1, pheno = 1)
p1.1.2 <- mix.MOE.plot(cell = 1, pheno = 1)
p1.2.1 <- ind.MOE.plot(cell = 1, pheno = 2)
p1.2.2 <- mix.MOE.plot(cell = 1, pheno = 2)
p1.3.1 <- ind.MOE.plot(cell = 1, pheno = 3)
p1.3.2 <- mix.MOE.plot(cell = 1, pheno = 3)
p1.4.1 <- ind.MOE.plot(cell = 1, pheno = 4)
p1.4.2 <- mix.MOE.plot(cell = 1, pheno = 4)
p1.5.1 <- ind.MOE.plot(cell = 1, pheno = 5)
p1.5.2 <- mix.MOE.plot(cell = 1, pheno = 5)
p1.6.1 <- ind.MOE.plot(cell = 1, pheno = 6)
p1.6.2 <- mix.MOE.plot(cell = 1, pheno = 6)
p1.7.1 <- ind.MOE.plot(cell = 1, pheno = 7)
p1.7.2 <- mix.MOE.plot(cell = 1, pheno = 7)
p1.8.1 <- ind.MOE.plot(cell = 1, pheno = 8)
p1.8.2 <- mix.MOE.plot(cell = 1, pheno = 8)
p1.9.1 <- ind.MOE.plot(cell = 1, pheno = 9)
p1.9.2 <- mix.MOE.plot(cell = 1, pheno = 9)
p1.10.1 <- ind.MOE.plot(cell = 1, pheno = 10)
p1.10.2 <- mix.MOE.plot(cell = 1, pheno = 10)

p1.1 <- arrangeGrob(p1.1.1, p1.1.2, nrow=2, heights = c(4/5,1/5))
p1.2 <- arrangeGrob(p1.2.1, p1.2.2, nrow=2, heights = c(4/5,1/5))
p1.3 <- arrangeGrob(p1.3.1, p1.3.2, nrow=2, heights = c(4/5,1/5))
p1.4 <- arrangeGrob(p1.4.1, p1.4.2, nrow=2, heights = c(4/5,1/5))
p1.5 <- arrangeGrob(p1.5.1, p1.5.2, nrow=2, heights = c(4/5,1/5))
p1.6 <- arrangeGrob(p1.6.1, p1.6.2, nrow=2, heights = c(4/5,1/5))
p1.7 <- arrangeGrob(p1.7.1, p1.7.2, nrow=2, heights = c(4/5,1/5))
p1.8 <- arrangeGrob(p1.8.1, p1.8.2, nrow=2, heights = c(4/5,1/5))
p1.9 <- arrangeGrob(p1.9.1, p1.9.2, nrow=2, heights = c(4/5,1/5))
p1.10 <- arrangeGrob(p1.10.1, p1.10.2, nrow=2, heights = c(4/5,1/5))

file = paste0("plots/suppl/FigS100.png")
png(file = file, width=2700, height=1800, res=150)
plot_grid(p1.1, p1.2, p1.3, p1.4, p1.5,
          p1.6, p1.7, p1.8, ncol = 4)
dev.off()

file = paste0("plots/suppl/FigS101.png")
png(file = file, width=2700, height=900, res=150)
plot_grid(p1.9, p1.10, ncol = 4)
dev.off()


## iCell HUVECs
p2.1.1 <- ind.MOE.plot(cell = 2, pheno = 1)
p2.1.2 <- mix.MOE.plot(cell = 2, pheno = 1)
p2.2.1 <- ind.MOE.plot(cell = 2, pheno = 2)
p2.2.2 <- mix.MOE.plot(cell = 2, pheno = 2)
p2.3.1 <- ind.MOE.plot(cell = 2, pheno = 3)
p2.3.2 <- mix.MOE.plot(cell = 2, pheno = 3)
p2.4.1 <- ind.MOE.plot(cell = 2, pheno = 4)
p2.4.2 <- mix.MOE.plot(cell = 2, pheno = 4)
p2.5.1 <- ind.MOE.plot(cell = 2, pheno = 5)
p2.5.2 <- mix.MOE.plot(cell = 2, pheno = 5)
p2.6.1 <- ind.MOE.plot(cell = 2, pheno = 6)
p2.6.2 <- mix.MOE.plot(cell = 2, pheno = 6)
p2.7.1 <- ind.MOE.plot(cell = 2, pheno = 7)
p2.7.2 <- mix.MOE.plot(cell = 2, pheno = 7)
p2.8.1 <- ind.MOE.plot(cell = 2, pheno = 8)
p2.8.2 <- mix.MOE.plot(cell = 2, pheno = 8)

p2.1 <- arrangeGrob(p2.1.1, p2.1.2, nrow=2, heights = c(4/5,1/5))
p2.2 <- arrangeGrob(p2.2.1, p2.2.2, nrow=2, heights = c(4/5,1/5))
p2.3 <- arrangeGrob(p2.3.1, p2.3.2, nrow=2, heights = c(4/5,1/5))
p2.4 <- arrangeGrob(p2.4.1, p2.4.2, nrow=2, heights = c(4/5,1/5))
p2.5 <- arrangeGrob(p2.5.1, p2.5.2, nrow=2, heights = c(4/5,1/5))
p2.6 <- arrangeGrob(p2.6.1, p2.6.2, nrow=2, heights = c(4/5,1/5))
p2.7 <- arrangeGrob(p2.7.1, p2.7.2, nrow=2, heights = c(4/5,1/5))
p2.8 <- arrangeGrob(p2.8.1, p2.8.2, nrow=2, heights = c(4/5,1/5))


file = paste0("plots/suppl/FigS102.png")
png(file = file, width=2700, height=1800, res=150)
plot_grid(p2.1, p2.2, p2.3, p2.4, p2.5, p2.6, p2.7, p2.8, nrow = 2)
dev.off()

## iCell Hepatocytes
p3.1.1 <- ind.MOE.plot(cell = 3, pheno = 1)
p3.1.2 <- mix.MOE.plot(cell = 3, pheno = 1)
p3.2.1 <- ind.MOE.plot(cell = 3, pheno = 2)
p3.2.2 <- mix.MOE.plot(cell = 3, pheno = 2)
p3.3.1 <- ind.MOE.plot(cell = 3, pheno = 3)
p3.3.2 <- mix.MOE.plot(cell = 3, pheno = 3)
p3.4.1 <- ind.MOE.plot(cell = 3, pheno = 4)
p3.4.2 <- mix.MOE.plot(cell = 3, pheno = 4)
p3.5.1 <- ind.MOE.plot(cell = 3, pheno = 5)
p3.5.2 <- mix.MOE.plot(cell = 3, pheno = 5)

p3.1 <- arrangeGrob(p3.1.1, p3.1.2, nrow=2, heights = c(4/5,1/5))
p3.2 <- arrangeGrob(p3.2.1, p3.2.2, nrow=2, heights = c(4/5,1/5))
p3.3 <- arrangeGrob(p3.3.1, p3.3.2, nrow=2, heights = c(4/5,1/5))
p3.4 <- arrangeGrob(p3.4.1, p3.4.2, nrow=2, heights = c(4/5,1/5))
p3.5 <- arrangeGrob(p3.5.1, p3.5.2, nrow=2, heights = c(4/5,1/5))


file = paste0("plots/suppl/FigS103.png")
png(file = file, width=2700, height=1800, res=150)
plot_grid(p3.1, p3.2, p3.3, p3.4, p3.5, ncol = 4)
dev.off()

## iCell Endothelial cells
p4.1.1 <- ind.MOE.plot(cell = 4, pheno = 1)
p4.1.2 <- mix.MOE.plot(cell = 4, pheno = 1)
p4.2.1 <- ind.MOE.plot(cell = 4, pheno = 2)
p4.2.2 <- mix.MOE.plot(cell = 4, pheno = 2)
p4.3.1 <- ind.MOE.plot(cell = 4, pheno = 3)
p4.3.2 <- mix.MOE.plot(cell = 4, pheno = 3)
p4.4.1 <- ind.MOE.plot(cell = 4, pheno = 4)
p4.4.2 <- mix.MOE.plot(cell = 4, pheno = 4)
p4.5.1 <- ind.MOE.plot(cell = 4, pheno = 5)
p4.5.2 <- mix.MOE.plot(cell = 4, pheno = 5)
p4.6.1 <- ind.MOE.plot(cell = 4, pheno = 6)
p4.6.2 <- mix.MOE.plot(cell = 4, pheno = 6)
p4.7.1 <- ind.MOE.plot(cell = 4, pheno = 7)
p4.7.2 <- mix.MOE.plot(cell = 4, pheno = 7)
p4.8.1 <- ind.MOE.plot(cell = 4, pheno = 8)
p4.8.2 <- mix.MOE.plot(cell = 4, pheno = 8)

p4.1 <- arrangeGrob(p4.1.1, p4.1.2, nrow=2, heights = c(4/5,1/5))
p4.2 <- arrangeGrob(p4.2.1, p4.2.2, nrow=2, heights = c(4/5,1/5))
p4.3 <- arrangeGrob(p4.3.1, p4.3.2, nrow=2, heights = c(4/5,1/5))
p4.4 <- arrangeGrob(p4.4.1, p4.4.2, nrow=2, heights = c(4/5,1/5))
p4.5 <- arrangeGrob(p4.5.1, p4.5.2, nrow=2, heights = c(4/5,1/5))
p4.6 <- arrangeGrob(p4.6.1, p4.6.2, nrow=2, heights = c(4/5,1/5))
p4.7 <- arrangeGrob(p4.7.1, p4.7.2, nrow=2, heights = c(4/5,1/5))
p4.8 <- arrangeGrob(p4.8.1, p4.8.2, nrow=2, heights = c(4/5,1/5))

file = paste0("plots/suppl/FigS104.png")
png(file = file, width=2700, height=1800, res=150)
plot_grid(p4.1, p4.2, p4.3, p4.4, p4.5, p4.6, p4.7, p4.8, nrow = 2)
dev.off()

## iCell Cardiomyocytes
p5.1.1 <- ind.MOE.plot(cell = 5, pheno = 1)
p5.1.2 <- mix.MOE.plot(cell = 5, pheno = 1)
p5.2.1 <- ind.MOE.plot(cell = 5, pheno = 2)
p5.2.2 <- mix.MOE.plot(cell = 5, pheno = 2)
p5.3.1 <- ind.MOE.plot(cell = 5, pheno = 3)
p5.3.2 <- mix.MOE.plot(cell = 5, pheno = 3)
p5.4.1 <- ind.MOE.plot(cell = 5, pheno = 4)
p5.4.2 <- mix.MOE.plot(cell = 5, pheno = 4)
p5.5.1 <- ind.MOE.plot(cell = 5, pheno = 5)
p5.5.2 <- mix.MOE.plot(cell = 5, pheno = 5)
p5.6.1 <- ind.MOE.plot(cell = 5, pheno = 6)
p5.6.2 <- mix.MOE.plot(cell = 5, pheno = 6)
p5.7.1 <- ind.MOE.plot(cell = 5, pheno = 7)
p5.7.2 <- mix.MOE.plot(cell = 5, pheno = 7)
p5.8.1 <- ind.MOE.plot(cell = 5, pheno = 8)
p5.8.2 <- mix.MOE.plot(cell = 5, pheno = 8)
p5.9.1 <- ind.MOE.plot(cell = 5, pheno = 9)
p5.9.2 <- mix.MOE.plot(cell = 5, pheno = 9)
p5.10.1 <- ind.MOE.plot(cell = 5, pheno = 10)
p5.10.2 <- mix.MOE.plot(cell = 5, pheno = 10)
p5.11.1 <- ind.MOE.plot(cell = 5, pheno = 11)
p5.11.2 <- mix.MOE.plot(cell = 5, pheno = 11)
p5.12.1 <- ind.MOE.plot(cell = 5, pheno = 12)
p5.12.2 <- mix.MOE.plot(cell = 5, pheno = 12)
p5.13.1 <- ind.MOE.plot(cell = 5, pheno = 13)
p5.13.2 <- mix.MOE.plot(cell = 5, pheno = 13)
p5.14.1 <- ind.MOE.plot(cell = 5, pheno = 14)
p5.14.2 <- mix.MOE.plot(cell = 5, pheno = 14)
p5.15.1 <- ind.MOE.plot(cell = 5, pheno = 15)
p5.15.2 <- mix.MOE.plot(cell = 5, pheno = 15)
p5.16.1 <- ind.MOE.plot(cell = 5, pheno = 16)
p5.16.2 <- mix.MOE.plot(cell = 5, pheno = 16)

p5.1 <- arrangeGrob(p5.1.1, p5.1.2, nrow=2, heights = c(4/5,1/5))
p5.2 <- arrangeGrob(p5.2.1, p5.2.2, nrow=2, heights = c(4/5,1/5))
p5.3 <- arrangeGrob(p5.3.1, p5.3.2, nrow=2, heights = c(4/5,1/5))
p5.4 <- arrangeGrob(p5.4.1, p5.4.2, nrow=2, heights = c(4/5,1/5))
p5.5 <- arrangeGrob(p5.5.1, p5.5.2, nrow=2, heights = c(4/5,1/5))
p5.6 <- arrangeGrob(p5.6.1, p5.6.2, nrow=2, heights = c(4/5,1/5))
p5.7 <- arrangeGrob(p5.7.1, p5.7.2, nrow=2, heights = c(4/5,1/5))
p5.8 <- arrangeGrob(p5.8.1, p5.8.2, nrow=2, heights = c(4/5,1/5))
p5.9 <- arrangeGrob(p5.9.1, p5.9.2, nrow=2, heights = c(4/5,1/5))
p5.10 <- arrangeGrob(p5.10.1, p5.10.2, nrow=2, heights = c(4/5,1/5))
p5.11 <- arrangeGrob(p5.11.1, p5.11.2, nrow=2, heights = c(4/5,1/5))
p5.12 <- arrangeGrob(p5.12.1, p5.12.2, nrow=2, heights = c(4/5,1/5))
p5.13 <- arrangeGrob(p5.13.1, p5.13.2, nrow=2, heights = c(4/5,1/5))
p5.14 <- arrangeGrob(p5.14.1, p5.14.2, nrow=2, heights = c(4/5,1/5))
p5.15 <- arrangeGrob(p5.15.1, p5.15.2, nrow=2, heights = c(4/5,1/5))
p5.16 <- arrangeGrob(p5.16.1, p5.16.2, nrow=2, heights = c(4/5,1/5))

file = paste0("plots/suppl/FigS105.png")
png(file = file, width=2700, height=1800, res=150)
plot_grid(p5.1, p5.2, p5.3, p5.4, p5.5, p5.6, p5.7, p5.8, nrow = 2)
dev.off()

file = paste0("plots/suppl/FigS106.png")
png(file = file, width=2700, height=1800, res=150)
plot_grid(p5.9, p5.10, p5.11, p5.12, p5.13, p5.14, p5.15, p5.16, nrow = 2)
dev.off()

cat(paste0("\nStarting time: ", Str.t), "\n\n")
cat(paste0("\nEnding time: ", Sys.time(), "\n\n"))
Sys.time() - Str.t 

