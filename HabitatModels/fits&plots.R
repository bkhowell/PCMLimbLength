# clear workspace 
rm(list=ls()) 

## load packages
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
library(ggplot2)
library(ggtree)
library(cowplot)
library(treeplyr)
library(phytools)
library(phangorn)
library(magrittr)
library(phytools)


##everything  fossorial&terrestrial as fossorial
fossallBM<-readRDS("foss_eco_fits_BM.RDS")
fossallOU<-readRDS("foss_eco_fits_OU.RDS")

##everything together except fossorial
fossBM<-readRDS("only_foss_eco_fits_BM.RDS")
fossOU<-readRDS("only_foss_eco_fits_OU.RDS")

##everything fossorial&terrestrial as fossorial&terrestrial
eco_fits_BM<-readRDS("eco_fits_BM.RDS")
eco_fits_OU<-readRDS("eco_fits_OU.RDS")

##everything fossorial&terrestrial as terrestrial
terrBM<-readRDS("terr_eco_fits_BM.RDS")
terrOU<-readRDS("terr_eco_fits_OU.RDS")

## set threshold for 0 eigenvalues
options(PCMBase.Threshold.EV = 0.00000001)

## find best model
min(unlist(lapply(eco_fits_BM, AIC)))
min(unlist(lapply(eco_fits_OU, AIC)))
min(unlist(lapply(fossBM, AIC))) 
min(unlist(lapply(fossOU, AIC)))#-863.0132
min(unlist(lapply(fossallBM, AIC))) #-845.4434
min(unlist(lapply(fossallOU, AIC)))
min(unlist(lapply(terrBM, AIC)))
min(unlist(lapply(terrOU, AIC)))

AIC(fossallBM[[96]])

## calculate slope and intercept
Vy_list <- vector(mode = "list", length = length(fossallBM[[96]])-1)
theta_list<- vector(mode= "list", length = length(fossallBM[[96]])-1)

StationaryVariance <- function(alpha,sigma){
  sigma <- sigma
  eig <- eigen(alpha)
  P <- eig$vectors
  invP <- solve(P)
  eigvalues <- eig$values
  p=dim(sigma)[1]
  Mat <- matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      Mat[i,j] <- 1/(eigvalues[i]+eigvalues[j])
    }
  }
  StVar <- P%*%(Mat*(invP%*%sigma%*%t(invP)))%*%t(P)
  return(StVar)
}

#for loop for getting Vy and thetas
for (i in 2:(length(fossallBM[[96]]))) {
  if (length(fossallBM[[96]][[i]])==3) {
    
    tempsig<-tcrossprod(fossallBM[[96]][[i]]$Sigma_x[1:2, 1:2,])
    Vy_list[[i-1]]<-StationaryVariance(alpha = PCMApplyTransformation(fossallBM[[96]][[i]])$H[1:2,1:2,], sigma = tempsig)
    theta_list[[i-1]]<-fossallBM[[96]][[i]]$Theta[1:2]
    
  } else {
    
    Vy_list[[i-1]]<-tcrossprod(fossallBM[[96]][[i]]$Sigma_x[1:2, 1:2,]) # using sig2 matrices for Vy for BM models
    # causes problems with ellipses and lines later on 
    theta_list[[i-1]]<- NaN 
  }
}

# trying to get away from for loops 
Vy_slope<- function(Vy) {
  covxy <- as.numeric(Vy[1,2])
  Vx<- as.numeric(Vy[1,1]) # I think I am grabbing the right number here because SVL from this run was trait 1 and predictor
  m<-covxy/Vx
  return(m)
}

slopeint<-matrix(data = NA, nrow = length(fossallBM[[96]])-1, ncol = 2)
rownames(slopeint)<-1:(length(fossallBM[[96]])-1)
colnames(slopeint)<-c("slope", "int")
slopeint[,"slope"]<-sapply(Vy_list, FUN = Vy_slope)
slopeint<-as.data.frame(slopeint)


thetas<-matrix(data = NA, nrow = length(fossallBM[[96]])-1, ncol = 2)
rownames(thetas)<-1:(length(fossallBM[[96]])-1)

for (i in 1:(length(fossallBM[[96]])-1)) {
  if(is.nan(theta_list[[i]][1])) {
    slopeint[i,2]<-fossallBM[[96]]$X0[2] - slopeint[i,1]*fossallBM[[96]]$X0[1]
  } else {
    thetax<-theta_list[[i]][1]
    thetay<-theta_list[[i]][2]
    slopeint[i,2]<-thetay-thetax*slopeint[i,1]
    thetas[i]<-thetay/thetax
  }
  
}

## plots of allometry
library(viridis)
cols<-viridis(5)
names(cols)<-1:5
eco_allometry_figs<-list()

for (i in 1:PCMNumRegimes(fossallBM[[96]])) {
  
  df<-as.data.frame(t(attributes(fossallBM[[96]])$X[,PCMTreeGetTipsInRegime(tree=attributes(fossallBM[[96]])$tree,regime = i)]))
  
  eco_allometry_figs[[i]] <- ggplot(data = df, aes(color = cols[i])) +
    geom_point(aes(x = SVL, y = TLL, color = cols[i]), color=cols[i]) + 
    theme_bw()+
    theme(panel.grid = element_blank())+
    coord_fixed()+
    xlab(element_blank())+
    ylab(element_blank()) +
    geom_abline(slope = slopeint[i,1], intercept = slopeint[i,2], color=cols[i], size=2) +
    xlim(3, 6) +
    ylim(0.5,5)
  
  
}

eco_combined_plot <- plot_grid(
  plotlist = eco_allometry_figs,  # the list of plots
  ncol = 3)+
  draw_label("snout-vent length (ln mm)", x= 0.5, y= 0, vjust=-0.5, angle= 0) +
  draw_label("hind limb length (ln mm)", x= 0, y= 0.5, vjust= 1.5, angle=90) 
  
  # Display the combined plot
print(eco_combined_plot)

## Plot tree
PCMTreePlot(attributes(fossallBM[[96]])$tree, palette = cols)

plot(summary(mtrees), colors=cols, ftype="off")
legend(legend = c("Terrestrial", "Arboreal", "Saxicolous", "Fossorial", "Semi-aquatic"), x=0, y= 275, col = cols, bty = 'n', pch = 19)




