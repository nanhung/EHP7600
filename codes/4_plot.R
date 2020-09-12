# Hint: Select all R script and run them ####
rm(list=ls())

# Load packages
library(gridGraphics)
library(tidyverse)
library(magrittr)
library(reshape2)
library(scales)
library(ggpubr)
library(treemapify)
library(cowplot)
library(bayestestR)
library(ggridges)
library(gridExtra)
library(PerformanceAnalytics)
library(tiff)
library(jpeg)

# Required data info
chem_data <- "datasets/chem_data.csv"
mixture_data <- "datasets/mixture_data.csv"
mixture_info <- "datasets/mixture_info.csv"
mix.dat <- read.csv(mixture_data)
mixtures <- c("AC50-L", "AC50-H", "Expo-L", "Expo-H", 
              "POD-L", "POD-H", "RFD-L", "RFD-H")      # Mixture name 
celltypes <- read.csv(mixture_data) %>% 
  select(celltype) %>% 
  unique() %>% unlist() %>% as.character()             # Cell type
cum_conc_info <- read.csv(mixture_info) %>% 
  `colnames<-`(c("Chemical", mixtures, "Class")) %>%
  melt() %>% 
  group_by(variable) %>% 
  summarize(CumConc = sum(value))                      # Designed concentration

# Generate the plots folder to store result
if(!dir.exists("plots")) dir.create("plots")

#**************************************************************************
# mixtures properties (fig. 2.) -------------------------------------------
#**************************************************************************

names_mix_info <- c("Chemical", "AC50-L", "AC50-H", "Expo-L", "Expo-H",
                    "POD-L", "POD-H", "RFD-L", "RFD-H", "Class")

theme_set(theme_pubclean())
p2.1 <- read.csv("datasets/mixture_info.csv") %>% 
  mutate(`AC50-L` = AC50.L/sum(AC50.L)) %>%
  mutate(`AC50-H` = AC50.H/sum(AC50.H)) %>%
  mutate(`Expo-L` = Expo.L/sum(Expo.L)) %>%
  mutate(`Expo-H` = Expo.H/sum(Expo.H)) %>%
  mutate(`POD-L` = POD.L/sum(POD.L)) %>%
  mutate(`POD-H` = POD.H/sum(POD.H)) %>%
  mutate(`RFD-L` = RFD.L/sum(RFD.L)) %>%
  mutate(`RFD-H` = RFD.H/sum(RFD.H)) %>%
  select(names_mix_info) %>%
  `colnames<-`(c("Chemical", "paste(AC[50], '-L')", "paste(AC[50], '-H')", "'Expo-L'", "'Expo-H'",
                 "'POD-L'", "'POD-H'", "'RFD-L'", "'RFD-H'", "Class")) %>%
  melt() %>%
  ggplot(aes(area = value, fill = Class, label=Chemical)) +
  geom_treemap() +
  facet_wrap(~variable,nrow = 2, dir = "v", labeller =  label_parsed) +
  geom_treemap_text(fontface = "italic", colour = "white",
                    place = "centre", grow = F, reflow = F) + 
  scale_fill_viridis_d(begin = 0.8, end = 0.2, option = "inferno") +
  theme(legend.position = "top",
        legend.text = element_text(colour = "black", size = 12, face = "bold"),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(colour = "black", size = 14, face = "bold")) 
p2.1

p2.2 <- read.csv("datasets/mixture_info.csv") %>% 
  `colnames<-`(names_mix_info) %>%
  melt() %>% group_by(variable) %>% 
  summarize(cum_conc = sum(value)) %>%
  ggplot(aes(x = cum_conc, y=fct_rev(variable), label=round(cum_conc, d=1))) + 
  geom_segment(aes(x=10,xend=cum_conc,y=fct_rev(variable), yend=fct_rev(variable)), 
               size=5, color="grey80") + 
  geom_text(vjust=0) +
  labs(x = expression(paste("Total concentration, ",mu, "M")),
       y = "") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels =trans_format("log10", scales::math_format(10^.x))) +
  scale_y_discrete(labels = c("RFD-H", "RFD-L", "POD-H", "POD-L",
                              "Expo-H", "Expo-L", bquote(paste(AC[50], "-H")), bquote(paste(AC[50], "-L")))) +
  theme(axis.title=element_text(size=16, face = "bold", color = "black"),
        axis.text.y=element_text(size=12, face = "bold", color = "black"),
        axis.text.x=element_text(size=14, face = "bold", color = "black")) 
p2.2

#tiff("plots/fig2.tiff", res=600, compression = "lzw", height = 8, width = 11, units="in")
pdf("plots/fig2.pdf", height = 8, width = 11)
plot_grid(p2.1, p2.2, nrow = 2, rel_heights = c(2/3, 1/3), labels = c('A', 'B'), label_size = 22)
dev.off()

#**************************************************************************
# Concentration response (fig. 4.) ----------------------------------------
#**************************************************************************

i=1 # Neuron
j=9 # Total Out Growth
k=2 # AC50-H
load(file = paste0("outputs/", celltypes[i], "_mixtures_parms.rda"))
n.response <- length(mix.list)
name.resp <- names(mix.list)  
mix.m <- mix.list[[j]]

iter <- 100  # only pick 100 iters to plot
mix.m <- tail(mix.m, iter) 
str <- 2*k-1
X <- mix.m[,c(str,str+1)] %>% as.data.frame()
dil <- 10^seq(-4.5, 0.5, 0.1)
for(l in seq(length(dil))){
  X[,2+l] <- 1 / ((1 + dil[l] / X[,1] )^X[,2])
}

cum_conc <- as.numeric(cum_conc_info[2,2])
ec10 <- mix.m[,2*k-1] * ((1/0.9)^(1/mix.m[,2*k])-1) * cum_conc
hdi10 <- hdi(ec10, ci=0.9) 
ec10_med <- median(ec10)

# Extract mixture data
mix_conc <- mix.dat %>% filter(celltype == celltypes[i] & mixture == mixtures[k] & phenotype == name.resp[j]) %>%
  mutate(Conc = Dilution * cum_conc) %>%
  select(Conc) %>% unlist()
mix_resp <- mix.dat %>% filter(celltype == celltypes[i] & mixture == mixtures[k] & phenotype == name.resp[j]) %>%
  select(Response) %>% unlist()


# Extract individual chemical data
x <- read.csv(chem_data) %>% filter(celltype == "iCell Neurons" & phenotype == "Total Outgrowth")
chem <- x$chemical %>% unique()

# Plot
load(file = paste0("outputs/", celltypes[i], "_42_chem_parms.rda"))
m <- ind.list[[j]]

#tiff("plots/fig4.tiff", res=600, compression = "lzw", height = 10, width = 16, units="in")
pdf("plots/fig4.pdf", height = 10, width = 16)
par(mfrow = c(7,7), mar= c(2,0,2,0), oma=c(4,45,5,2))
for (k in 1:42){ # 42 individual chemicals
  
  dat <- x %>% filter(chemical == chem[k]) %>% select(Response)
  dat$conc <- rep(10^(c(-2:2)), each=2)
  
  iter <- 100 # only pick 100 iters to plot
  m <- tail(m, iter) 
  ec10 <- m[,2*k-1] * ((1/0.9)^(1/m[,2*k])-1)
  hdi_10 <- hdi(ec10, ci=0.9) 
  
  str <- 2*k-1
  M <- m[,c(str,str+1)] %>% as.data.frame()
  d <- 100/10^seq(-0.5, 4.5, 0.1)
  for(l in seq(length(d))){
    M[,2+l] <- 1 / ((1 + d[l] / M[,1] )^M[,2])
  }
  
  plot(d, M[1, c(3:(length(d)+2))], log = "x", type="l", col=alpha(rgb(0,0,0), 0.2), 
       ylim = c(0, 1.6), 
       xlab = "", ylab = "", main = chem[k], cex.main=0.9, 
       frame.plot = F, axes=F)
  axis(side=1, labels=FALSE)
  axis(side=2, labels=FALSE, tcl=0)
  for(iter in 2:iter){
    lines(d, M[iter, c(3:(length(d)+2))], col=alpha(rgb(0,0,0), 0.2))
  }
  abline(0.9,0, lty=3, col= "grey20")
  
  lines(c(hdi_10$CI_low, hdi_10$CI_high), c(0.9,0.9), col = "red")
  abline(v=hdi_10$CI_low, col="red", lty=2)
  abline(v=hdi_10$CI_high, col="red", lty=2)
  points(dat$conc, dat$Response, col = "darkgreen")
  
  y <- dat$Response
  d <- dat$conc
  param <- m[,c(str,str+1)]
  ypred <- matrix(ncol = length(y), nrow = iter) %>% as.data.frame()
  for(l in seq(length(d))){
    ypred[l] <- 1 / ((1 + d[l] / param[,1] )^param[,2])
  }
  text(0.01, 1.5, paste(round(median(ec10), digits = 2), "[", round(hdi_10$CI_low, digits = 2), "-", round(hdi_10$CI_high, digits = 2), "]"), adj = 0)
  
} # end chemicals

par(new = T, mfrow=c(1,1),mar=c(4,5,0,1), oma=c(0,0,3,0))
plot(dil*cum_conc, X[1, c(3:(length(dil)+2))],
     yaxt="n",
     log = "x", type="l", col=alpha(rgb(0,0,0), 0.1), ylim = c(0, 1.4), 
     xlab = "", 
     ylab = "",
     cex.axis=1.6)
axis(side=2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), 
     labels=c("0.0","0.2","0.4","0.6","0.8","1.0"), las = 2, cex.axis=1.6) 
mtext("Effect (Fraction of Vehicle)", side=2, line=3, cex=2)
mtext(expression(paste("Total concentration, ",mu, "M")), side=1, line=3, cex=2)

y <- mix_resp
d <- mix_conc / as.numeric(cum_conc_info[2,2])
k <- 2
str <- 2*k-1
param <- mix.m[,c(str,str+1)]
ypred <- matrix(ncol = length(y), nrow = iter) %>% as.data.frame()
for(l in seq(length(d))){
  ypred[l] <- 1 / ((1 + d[l] / param[,1] )^param[,2])
}

title <- expression(paste("iCell Neurons (Total Outgrowth) / AC"[50], "-H; EC"[10], ": ", "0.29 [0.23 - 0.33]"))

mtext(title, side=3, line=0, cex=2.5, adj=0.3, outer=TRUE)  
for(iter in 2:iter){
  lines(dil*cum_conc, X[iter, c(3:(length(dil)+2))], col=alpha(rgb(0,0,0), 0.2))
}

lines(c(hdi10$CI_low, hdi10$CI_high), c(0.9,0.9), col = "red")
abline(v=hdi10$CI_low, col="red", lty=2)
abline(v=hdi10$CI_high, col="red", lty=2)
points(mix_conc, mix_resp, col = "green", pch = 19)
dev.off()

#**************************************************************************
# Bayesian estimated EC10 (fig. 5.) ---------------------------------------
#**************************************************************************

for (i in 1:5){
  
  load(file = paste0("outputs/", celltypes[i], "_mixtures_parms.rda")) 
  
  for (j in 1:length(names(mix.list))){
    
    resp <- names(mix.list)[j]
    M <- mix.list[[j]] # select parameter matrix for specific response
    
    for (k in 1:8){ # 8 mixtures
      
      ec10 <- M[,k*2-1] * ((1/0.9)^(1/M[,k*2])-1) * as.numeric(cum_conc_info[k,2])
      
      if (k == 1) {
        EC10 <- ec10
      } else{
        EC10 <- c(EC10, ec10)
      }
    }
    EC10_mix <- matrix(EC10, ncol=8) %>% as.data.frame() %>% 
      `colnames<-`(mixtures) %>%
      reshape::melt() %>% mutate(response = names(mix.list)[j])
    
    if (j==1){
      EC10_resp <- EC10_mix
    } else {
      EC10_resp <- rbind(EC10_resp, EC10_mix)
    }
  }
  if (i==1){
    EC10_resp$celltype <- celltypes[i]
    EC10_cell <- EC10_resp
  } else {
    EC10_resp$celltype <- celltypes[i]
    EC10_cell <- rbind(EC10_cell, EC10_resp)
  }
}

EC10_df <- EC10_cell %>% group_by(response, variable, celltype) %>% 
  summarise(EC10.map = median(value),
            EC10.low = hdi(value, ci=0.9)$CI_low,
            EC10.high = hdi(value, ci=0.9)$CI_high)
EC10_df$label <- with(EC10_df, paste(response, variable, celltype, sep = "_"))

# Estimate positive rate
pos.rate <- merge(EC10_df, cum_conc_info, by = "variable", all.x = TRUE) %>% 
  mutate(ratio = EC10.map/CumConc) %>% mutate(pos = ifelse(ratio <= 1, 1,0)) %>% as.data.frame() %>%
  group_by(variable, CumConc) %>% summarize(sum(pos)/length(pos)*100) %>%
  `colnames<-`(c("variable", "CumConc", "positive_rate")) %>%
  mutate(label2 = paste0("'Phenotypes \"active\": ", round(positive_rate, 1), "%'"))

EC10_df2 <- merge(EC10_df, pos.rate, by = "variable", all.x = TRUE)

EC10_df2 %<>%
  mutate(variable2 = c(
    rep("paste(AC[50], '-H')", 47), 
    rep("paste(AC[50], '-L')", 47), 
    rep("'Expo-H'", 47), 
    rep("'Expo-L'", 47),
    rep("'POD-H'", 47),
    rep("'POD-L'", 47),
    rep("'RFD-H'", 47),
    rep("'RFD-L'", 47)  
  ))

var.order <- c("paste(AC[50], '-L')", "paste(AC[50], '-H')", "'Expo-L'", "'Expo-H'",
               "'POD-L'", "'POD-H'", "'RFD-L'", "'RFD-H'")
EC10_df2$variable2 <- factor(EC10_df2$variable2, 
                             levels = var.order)
pos.rate$variable2 <- var.order
EC10_df2$variable2 <- factor(EC10_df2$variable2, levels = var.order)
pos.rate$variable2 <- factor(pos.rate$variable2, levels = var.order)

theme_set(theme_bw())

p5.1 <- EC10_df2 %>%
  ggplot(aes(x=reorder(label, EC10.map), y=EC10.map, color=celltype, shape=celltype)) +
  facet_wrap(~variable2+label2, scales = 'free_y', ncol = 4, dir = "v", labeller =  label_parsed) +
  geom_hline(data = pos.rate, aes(yintercept = CumConc), col = "grey50", lty = 2) +
  coord_flip() +
  geom_errorbar(aes(ymin=EC10.low, ymax=EC10.high), size = 0.15,col="grey50", width=0) +
  geom_point(aes(), size = 1.5,fill="black") +
  scale_y_log10(lim = c(2*10^-3, 3*10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", scales::math_format(10^.x))) +
  scale_x_discrete(breaks = NULL) +
  scale_colour_grey(start = 0.0, end = 0.2)+ 
  scale_shape_manual(values=c(1, 15, 19, 17, 23)) +
  labs(y = "") +
  guides(colour=guide_legend(override.aes=list(size=2))) +
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(colour = "black", size = 14, face = "bold"),
        legend.title=element_blank(),
        legend.position = "top",
        legend.text=element_text(size=12, face = "bold"),
        axis.text=element_text(size=10),
        axis.ticks.x = element_line (colour = "black"), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=12, face = "bold"),
        legend.key.size = unit(1, "cm"),
        axis.text.y=element_blank())
p5.1

p5.2 <- EC10_cell %>%  
  ggplot(aes(x = value, y = fct_rev(variable), group = variable)) + 
  geom_density_ridges(stat = "binline", bins = 100, scale = 0.95, draw_baseline = FALSE, alpha = 0.7) +
  scale_x_log10(lim = c(2*10^-3, 3*10^5),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_discrete(labels = c("RFD-H", "RFD-L", "POD-H", "POD-L",
                              "Expo-H", "Expo-L", bquote(paste(AC[50], "-H")), bquote(paste(AC[50], "-L")))) +
  annotation_logticks(sides = "b") +
  labs(y = "", x= expression(paste("EC"[10]," (Total concentration, ",mu, "M)"))) +
  theme_bw()+
  theme(axis.title.x = element_text(hjust = 0.5), 
        legend.position = "none",
        axis.title=element_text(size=16, face = "bold", color = "black"),
        axis.text.y=element_text(size=14, face = "bold", color = "black"),
        axis.text.x=element_text(size=14, face = "bold", color = "black"))
p5.2

#tiff("plots/fig5.tiff", res=600, compression = "lzw", height = 12, width = 10, units="in")
pdf("plots/fig5.pdf", height = 12, width = 10)
plot_grid(p5.1, p5.2, nrow = 2, rel_heights = c(3/4, 1/4), labels = c('A', 'B'), label_size = 22)
dev.off()

#**************************************************************************
# IA & CA (Fig. 6) --------------------------------------------------------
#**************************************************************************

mixplot <- function(cell = 1, pheno = 1, mix = "AC50-H", main = "main"){
  
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
       ylab="",
       xlab="")
  title(main, line = 1, cex=1)
  
  at=seq(-1, 5, 2)
  at.minor <- log10(outer(1:9, 10^(min(at):max(at))))
  lab <- sapply(at, function(i) as.expression(bquote(10^ .(i))))
  axis(side=1, at=10^at, labels=lab, cex.axis=1.6)
  axis(side=2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.3, 1.4), 
       labels=c("0.0","0.2","0.4","0.6","0.8","1.0", "IA", "CA", "Fitting"), las = 2, cex.axis=1.6) 
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

load(file = "outputs/ec10_pred.rda")
CF <- ec10.list[["AC50-H"]] %>% filter(Method == "CF") %>%
  group_by(Phenotype, Celltype) %>% 
  summarise(CF.500 = median(EC10), 
            CF.025 = hdi(EC10, ci=0.9)$CI_low,
            CF.975 = hdi(EC10, ci=0.9)$CI_high) %>%
  as.data.frame()
CA <- ec10.list[["AC50-H"]] %>% filter(Method == "CA") %>%
  group_by(Phenotype, Celltype) %>% 
  summarise(CA.500 = median(EC10), 
            CA.025 = hdi(EC10, ci=0.9)$CI_low,
            CA.975 = hdi(EC10, ci=0.9)$CI_high) %>%
  as.data.frame()
IA <- ec10.list[["AC50-H"]] %>% filter(Method == "IA") %>%
  group_by(Phenotype, Celltype) %>% 
  summarise(IA.500 = median(EC10), 
            IA.025 = hdi(EC10, ci=0.9)$CI_low,
            IA.975 = hdi(EC10, ci=0.9)$CI_high) %>%
  as.data.frame()

ec10.df <- cbind(CF, IA[,c(3:5)], CA[,c(3:5)])
ec10.df$Celltype <- as.character(ec10.df$Celltype)

ec10.df$CF.rank <- rank(ec10.df$CF.500)
X <- data.frame(ec10.df$Phenotype, ec10.df$CF.rank)
ec10.rank.label <- X[order(-X$ec10.df.CF.rank),][,1]

dat_text <- data.frame(label = c("IA", "CA"), 
                       value = c(min(ec10.df$IA.500), min(ec10.df$CA.500)))

cum_conc <- as.numeric(cum_conc_info[2,2])
theme_set(theme_pubclean())
rng <- range(ec10.df[,c(3:11)])

p6.2 <- ec10.df %>% 
  unite("celltype_response", Celltype:Phenotype, remove = FALSE) %>%
  ggplot(aes(x=reorder(celltype_response, CF.500), y=CF.500)) +
  geom_hline(yintercept = cum_conc, col = "grey20", lty = 2) +
  geom_pointrange(aes(ymin=CF.025, ymax=CF.975, color=Celltype, shape=Celltype)) +
  scale_y_log10(lim = rng,
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_segment(aes(x=celltype_response, xend=celltype_response, y=IA.025, yend=IA.975),
               color = "green", size =2.5, alpha=0.1) +
  geom_segment(aes(x=celltype_response, xend=celltype_response, y=CA.025, yend=CA.975),
               color = "blue", size =2, alpha=0.1) +
  coord_flip() +
  geom_point(aes(x=celltype_response, y=IA.500), color = "green", shape = "l", size = 2) +
  geom_point(aes(x=celltype_response, y=CA.500), color = "blue", shape = "l", size = 2) +
  scale_colour_grey(start = 0.0, end = 0.2)+ 
  scale_shape_manual(values=c(1, 15, 19, 17, 23)) +
  theme(plot.margin = margin(0, 0.4, 0, 0, "cm"),
        legend.text=element_text(color = "black", size =13),
        axis.text.x=element_blank(),
        axis.text.y=element_text(color = "black", size =11),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "top",
        legend.title = element_blank()) +
  scale_x_discrete(labels = rev(ec10.rank.label)) +
  annotate("text", label = paste("Maximun (design/exposure) concentration"), x = 0.2, hjust = 1, y = 9000, angle = 270)

p6.3 <- ec10.df %>% select(CF.500, IA.500, CA.500, Phenotype, Celltype) %>%
  `colnames<-`(c("Fitting", "Independent \naction (IA)", "Concentration \naddition (CA)",
                 "Phenotype", "Celltype")) %>% melt() %>%
  ggplot(aes(x = variable, y = value)) +
  geom_boxplot() +
  scale_y_log10(lim = rng,
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept = cum_conc, col = "grey20", lty = 2) +
  coord_flip() +
  labs(x = "",
       y = expression(paste("EC"[10]," (Total concentration, ",mu, "M)"))) +
  theme(plot.margin = margin(0, 0.4, 0, 0, "cm"),
        axis.title.x = element_text(color = "black", size =16),
        axis.text.x=element_text(color = "black", size =14),
        axis.text.y=element_text(color = "black", size =14))

par(mfrow = c(3,1), mar=c(2,4,4,1), oma=c(3,3,0,1))  
mixplot(1, 9, main = "iCell Neurons (Total outgrowth)")
mixplot(5, 1, main = "iCell Cardiomyocytes (BPM, 15 min)")
mtext("Effect (Fraction of Vehicle)", side=2, line=5, cex=1.5)
mixplot(2, 3, main = "HUVECs (Mean Tube Length)")
mtext(expression(paste("Total concentration, ",mu, "M")), side=1, line=1, cex=1.5, outer=TRUE)  
p6.1 <- recordPlot()

#tiff("plots/fig6.tiff", res=600, compression = "lzw", height = 10, width = 16, units="in")
pdf("plots/fig6.pdf", height = 10, width = 16)
plot_grid(
  p6.1, 
  plot_grid(p6.2, p6.3, nrow = 2, rel_heights = c(3/4, 1/4), align = "v", axis = "l", labels = c('B', 'C'), label_size = 30),
  rel_widths = c(1/3, 2/3), labels = c('A', ''), label_size = 30)
dev.off()

#**************************************************************************
# MOE (Fig.7) --------------------------------------------------------------
#**************************************************************************

load(file = "outputs/ec10_pred.rda")
CF <- ec10.list[["AC50-H"]] %>% filter(Method == "CF") %>%
  group_by(Phenotype, Celltype) %>% 
  summarise(CF.500 = median(EC10), 
            CF.025 = hdi(EC10, ci=0.9)$CI_low,
            CF.975 = hdi(EC10, ci=0.9)$CI_high) %>%
  as.data.frame()
CA <- ec10.list[["AC50-H"]] %>% filter(Method == "CA") %>%
  group_by(Phenotype, Celltype) %>% 
  summarise(CA.500 = median(EC10), 
            CA.025 = hdi(EC10, ci=0.9)$CI_low,
            CA.975 = hdi(EC10, ci=0.9)$CI_high) %>%
  as.data.frame()
IA <- ec10.list[["AC50-H"]] %>% filter(Method == "IA") %>%
  group_by(Phenotype, Celltype) %>% 
  summarise(IA.500 = median(EC10), 
            IA.025 = hdi(EC10, ci=0.9)$CI_low,
            IA.975 = hdi(EC10, ci=0.9)$CI_high) %>%
  as.data.frame()

ec10.df <- cbind(CF, IA[,c(3:5)], CA[,c(3:5)])
ec10.df$Celltype <- as.character(ec10.df$Celltype)


theme_set(theme_pubclean())

ind.MOE.plot <- function(cell = 1, pheno = 1, mix = "AC50-H", title, xlab) {
  
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
                  breaks = c(10^-5, 10^-1, 10^3),
                  labels = c(expression(10^-5), expression(10^-1), expression(10^3))
                  #breaks = scales::trans_breaks("log10", function(x) 10^x),
                  #labels = scales::trans_format("log10", scales::math_format(10^.x))
    )+
    coord_flip() +
    theme_pubclean() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          panel.border=element_blank(),
          axis.title=element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(colour = "black", size =0.5),
          title = element_text(size=10, face = "bold")) +
    labs(y = "", x = xlab, title = title)
}

mix.MOE.plot <- function(cell=1, pheno=1, mix = "AC50-H"){
  load(file = paste0("outputs/", celltypes[cell], "_mixtures_parms.rda")) # mixture results
  load(file = "outputs/ec10_pred.rda")
  celltypes <- c("iCell Neurons", "HUVECs", "iCell Hepatocytes", "iCell Endothelial cells", "iCell Cardiomyocytes") 
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
                  breaks = c(10^-5, 10^-1, 10^3),
                  labels = c(expression(10^-5), expression(10^-1), expression(10^3))
                  #breaks = scales::trans_breaks("log10", function(x) 10^x),
                  #labels = scales::trans_format("log10", scales::math_format(10^.x))
    )+
    coord_flip() +
    theme_pubclean() +
    theme(axis.text.x=element_text(color = "black", face="bold", size =10),
          axis.text.y=element_blank(),
          panel.border=element_blank(),
          axis.title=element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none",
          axis.ticks.x = element_line(colour = "black", size =0.5)) +
    annotate("text", x = 1:3, y = 1e-12, 
             label = c("Fitting", "Independent action", "Concentration addition"),
             hjust = 0, size = 3) +
    ggpubr::color_palette("jco") +
    labs(y = "Margin of Exposure", x = xlab)  
}

p7.1 <- ind.MOE.plot(cell = 1, pheno = 9, title = "iCell Neurons (Total Outgrowth)", xlab = "")
p7.2 <- ind.MOE.plot(cell = 5, pheno = 1, title = "iCell Cardiomyocytes (BPM, 15 min)", xlab = "")
p7.3 <- ind.MOE.plot(cell = 2, pheno = 3, title = "HUVECs (Mean Tube Length)", xlab = "")
p7.4 <- mix.MOE.plot(cell = 1, pheno = 9)
p7.5 <- mix.MOE.plot(cell = 5, pheno = 1)
p7.6 <- mix.MOE.plot(cell = 2, pheno = 3)

p7 <- arrangeGrob(p7.1, p7.2, p7.3,  
                  p7.4, p7.5, p7.6, 
                  nrow=2, heights = c(4/5,1/5),
                  bottom="Margin of Exposure", top = "")

#tiff("plots/fig7.tiff", res=600, compression = "lzw", height = 6, width = 10, units="in")
pdf("plots/fig7.pdf", height = 7, width = 12)
plot_grid(p7)
dev.off()

#**************************************************************************
# MOE (Fig.8) -------------------------------------------------------------
#**************************************************************************

mixtures <- c("AC50-L", "AC50-H", "Expo-L", "Expo-H", "POD-L", "POD-H", "RFD-L", "RFD-H") 
indConc <- read.csv("datasets/mixture_info.csv") %>% 
  `colnames<-`(c("Chemical", mixtures, "Class")) %>% select(mixtures[2])

df.mixMOE <- ec10.list[[mixtures[2]]] %>% mutate(MOE =  EC10/sum(indConc))
df.mixMOE$Celltype <- as.character(df.mixMOE$Celltype)
df.mixMOE$Celltype[which(df.mixMOE$Celltype == "HUVECs")] <- "HUVECs (8)"
df.mixMOE$Celltype[which(df.mixMOE$Celltype == "iCell Neurons")] <- "iCell \nNeurons (10)"
df.mixMOE$Celltype[which(df.mixMOE$Celltype == "iCell Hepatocytes")] <- "iCell \nHepatocytes (5)"
df.mixMOE$Celltype[which(df.mixMOE$Celltype == "iCell Cardiomyocytes")] <- "iCell \nCardiomyocytes (16)"
df.mixMOE$Celltype[which(df.mixMOE$Celltype == "iCell Endothelial cells")] <- "iCell \nEndothelial cells (8)"
df.mixMOE$Celltype <- factor(df.mixMOE$Celltype, 
                             levels = c("iCell \nNeurons (10)",
                                        "iCell \nCardiomyocytes (16)",
                                        "HUVECs (8)",
                                        "iCell \nEndothelial cells (8)",
                                        "iCell \nHepatocytes (5)"))

df.mixMOE$Method <- as.character(df.mixMOE$Method)
df.mixMOE$Method[which(df.mixMOE$Method == "CF")] <- "Fitting"
df.mixMOE$Method[which(df.mixMOE$Method == "IA")] <- "Independent \naction (IA)"
df.mixMOE$Method[which(df.mixMOE$Method == "CA")] <- "Concentration \naddition (CA)"
df.mixMOE$Method <- factor(df.mixMOE$Method, 
                           levels = c("Fitting", "Independent \naction (IA)", "Concentration \naddition (CA)"))

p8.a <- df.mixMOE %>% 
  ggplot(aes(y = fct_rev(Celltype))) +
  geom_density_ridges(
    aes(x = MOE, fill = Method), scale = 0.8, alpha = .4, color = "grey") +
  scale_fill_brewer(palette = "Blues") + 
  scale_x_log10(lim = c(10^-6, 10^2),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_pubr() +
  theme(panel.background = element_rect (fill="#f3f3f3"),
        axis.text.y=element_text(color = "black", face="bold", size =10),
        legend.position="top",
        legend.text=element_text(color = "black", face="bold"),
        legend.title=element_blank(),
        axis.text.x=element_blank(),
        axis.title=element_blank())

p8.b <- df.mixMOE %>%
  ggplot(aes(y = Method)) +
  geom_density_ridges(
    aes(x = MOE, fill = Method), scale = 0.8, alpha = .4, color = "grey") +
  scale_x_log10(lim = c(10^-6, 10^2),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_pubr() +
  scale_fill_brewer(palette = "Blues") + 
  labs(x="Margin of Exposure")+ 
  theme(panel.background = element_rect (fill="#f3f3f3"),
        axis.text.y=element_text(color = "black", face="bold", size =10),
        axis.text.x=element_text(color = "black", face="bold", size =12),
        axis.title.y=element_blank(),
        axis.title.x=element_text(color = "black", face="bold", size =12),
        legend.position="none")

#tiff("plots/fig8.tiff", res=600, compression = "lzw", height = 7, width = 8, units="in")
pdf("plots/fig8.pdf", height = 7, width = 8)
plot_grid(p8.a, p8.b, nrow = 2, rel_heights = c(2/3, 1/3), labels = c("A", "B"), align = 'v', axis = "l", label_size = 30)
dev.off()


# Transfer tiff to jpeg file to use in word document ------------
#img <- readTIFF("plots/fig2.tiff", native=TRUE)
#writeJPEG(img, target = "plots/fig2.jpeg", quality = 0.5)
#img <- readTIFF("plots/fig4.tiff", native=TRUE)
#writeJPEG(img, target = "plots/fig4.jpeg", quality = 0.5)
#img <- readTIFF("plots/fig5.tiff", native=TRUE)
#writeJPEG(img, target = "plots/fig5.jpeg", quality = 0.5)
#img <- readTIFF("plots/fig6.tiff", native=TRUE)
#writeJPEG(img, target = "plots/fig6.jpeg", quality = 0.5)
#img <- readTIFF("plots/fig7.tiff", native=TRUE)
#writeJPEG(img, target = "plots/fig7.jpeg", quality = 0.5)
#img <- readTIFF("plots/fig8.tiff", native=TRUE)
#writeJPEG(img, target = "plots/fig8.jpeg", quality = 0.5)
