# clear workspace 
rm(list=ls()) 
setwd("../HabitatModels")

library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(phytools)


td_eco<-readRDS("ext2_td.RDS")

## define different transition matrices
fiveState_q_restricted_FASA <- matrix(byrow = T, nrow = 5, ncol = 5, data = c(0,1,2,3,4, 5,0,6,0,0, 7,8,0,0,9, 10,0,0,0,0, 11,0,12,0,0))

fiveState_q_ARD <- matrix(byrow = T, nrow = 5, ncol = 5, data = c(0,1,2,3,4, 5,0,6,7,8, 9,10,0,11,12, 13,14,15,0,16, 17,18,19,20,0))

fiveState_q_SYM <- matrix(byrow = T, nrow = 5, ncol = 5, data = c(0,1,2,3,4, 1,0,5,6,7, 2,5,0,8,9, 3,6,8,0,10, 4,7,9,10,0))

fiveState_q_restricted <- matrix(byrow = T, nrow = 5, ncol = 5, data = c(0,1,2,3,4, 5,0,6,0,0, 7,8,0,0,0, 9,0,0,0,0, 10,0,0,0,0))

fiveState_q_restricted_restricted <- matrix(byrow = T, nrow = 5, ncol = 5, data = c(0,1,2,3,4, 5,0,0,0,0, 6,0,0,0,0, 7,0,0,0,0, 8,0,0,0,0))

## everything Fossorial&Terrestrial as Terrestrial

td_eco$dat$substrate2<-td_eco$dat$substrate
td_eco$dat$substrate2 <- dplyr::case_when(
  td_eco$dat$substrate2 == "Cryptic&Fossorial" ~ "Fossorial",
  td_eco$dat$substrate2 == "Cryptic&Terrestrial" ~ "Terrestrial",
  td_eco$dat$substrate2 == "Fossorial&Terrestrial" ~ "Terrestrial",
  TRUE ~ as.character(td_eco$dat$substrate2)
) 


td_eco$dat$substrate2 <- dplyr::case_when(
  td_eco$dat$substrate2 == "Terrestrial" ~ "1",
  td_eco$dat$substrate2 == "Arboreal" ~ "2",
  td_eco$dat$substrate2 == "Saxicolous" ~ "3",
  td_eco$dat$substrate2 == "Fossorial" ~ "4",
  td_eco$dat$substrate2 == "Semi_Aquatic" ~ "5",
  td_eco$dat$substrate2 == "Arboreal&Terrestrial" ~ "1&2",
  td_eco$dat$substrate2 == "Saxicolous&Terrestrial" ~ "1&3",
  td_eco$dat$substrate2 == "Fossorial&Terrestrial" ~ "1&4",
  td_eco$dat$substrate2 == "Arboreal&Saxicolous" ~ "2&3",
  td_eco$dat$substrate2 == "Arboreal&Saxicolous&Terrestrial" ~ "1&2&3",
  td_eco$dat$substrate2 == "Fossorial&Saxicolous&Terrestrial" ~ "1&3&4"
) 

## everything Fossorial&Terrestrial as Fossorial&Terrestrial

td_eco$dat$substrate3<-td_eco$dat$substrate
td_eco$dat$substrate3 <- dplyr::case_when(
  td_eco$dat$substrate3 == "Cryptic&Fossorial" ~ "Fossorial",
  td_eco$dat$substrate3 == "Cryptic&Terrestrial" ~ "Terrestrial",
  td_eco$dat$substrate3 == "Fossorial&Terrestrial" ~ "Fossorial&Terrestrial",
  TRUE ~ as.character(td_eco$dat$substrate3)
) 


td_eco$dat$substrate3 <- dplyr::case_when(
  td_eco$dat$substrate3 == "Terrestrial" ~ "1",
  td_eco$dat$substrate3 == "Arboreal" ~ "2",
  td_eco$dat$substrate3 == "Saxicolous" ~ "3",
  td_eco$dat$substrate3 == "Fossorial" ~ "4",
  td_eco$dat$substrate3 == "Semi_Aquatic" ~ "5",
  td_eco$dat$substrate3 == "Arboreal&Terrestrial" ~ "1&2",
  td_eco$dat$substrate3 == "Saxicolous&Terrestrial" ~ "1&3",
  td_eco$dat$substrate3 == "Fossorial&Terrestrial" ~ "1&4",
  td_eco$dat$substrate3 == "Arboreal&Saxicolous" ~ "2&3",
  td_eco$dat$substrate3 == "Arboreal&Saxicolous&Terrestrial" ~ "1&2&3",
  td_eco$dat$substrate3 == "Fossorial&Saxicolous&Terrestrial" ~ "1&3&4"
) 


## everything Fossorial&Terrestrial as Fossorial

td_eco$dat$substrate4<-td_eco$dat$substrate
td_eco$dat$substrate4 <- dplyr::case_when(
  td_eco$dat$substrate4 == "Cryptic&Fossorial" ~ "Fossorial",
  td_eco$dat$substrate4 == "Cryptic&Terrestrial" ~ "Terrestrial",
  td_eco$dat$substrate4 == "Fossorial&Terrestrial" ~ "Fossorial",
  TRUE ~ as.character(td_eco$dat$substrate4)
) 


td_eco$dat$substrate4 <- dplyr::case_when(
  td_eco$dat$substrate4 == "Terrestrial" ~ "1",
  td_eco$dat$substrate4 == "Arboreal" ~ "2",
  td_eco$dat$substrate4 == "Saxicolous" ~ "3",
  td_eco$dat$substrate4 == "Fossorial" ~ "4",
  td_eco$dat$substrate4 == "Semi_Aquatic" ~ "5",
  td_eco$dat$substrate4 == "Arboreal&Terrestrial" ~ "1&2",
  td_eco$dat$substrate4 == "Saxicolous&Terrestrial" ~ "1&3",
  td_eco$dat$substrate4 == "Fossorial&Terrestrial" ~ "1&4",
  td_eco$dat$substrate4 == "Arboreal&Saxicolous" ~ "2&3",
  td_eco$dat$substrate4 == "Arboreal&Saxicolous&Terrestrial" ~ "1&2&3",
  td_eco$dat$substrate4 == "Fossorial&Saxicolous&Terrestrial" ~ "1&3&4"
) 


## everything except Fossorial == state 1 and anything Fossorial == state 2 

td_eco$dat$substrate5<-td_eco$dat$substrate
td_eco$dat$substrate5 <- dplyr::case_when(
  td_eco$dat$substrate5 == "Cryptic&Fossorial" ~ "Fossorial",
  td_eco$dat$substrate5 == "Cryptic&Terrestrial" ~ "Terrestrial",
  td_eco$dat$substrate5 == "Fossorial&Terrestrial" ~ "Fossorial",
  TRUE ~ as.character(td_eco$dat$substrate5)
) 

td_eco$dat$substrate5 <- dplyr::case_when(
  td_eco$dat$substrate5 == "Terrestrial" ~ "1",
  td_eco$dat$substrate5 == "Arboreal" ~ "1",
  td_eco$dat$substrate5 == "Saxicolous" ~ "1",
  td_eco$dat$substrate5 == "Fossorial" ~ "2",
  td_eco$dat$substrate5 == "Semi_Aquatic" ~ "1",
  td_eco$dat$substrate5 == "Arboreal&Terrestrial" ~ "1",
  td_eco$dat$substrate5 == "Saxicolous&Terrestrial" ~ "1",
  td_eco$dat$substrate5 == "Fossorial&Terrestrial" ~ "2",
  td_eco$dat$substrate5 == "Arboreal&Saxicolous" ~ "1",
  td_eco$dat$substrate5 == "Arboreal&Saxicolous&Terrestrial" ~ "1",
  td_eco$dat$substrate5 == "Fossorial&Saxicolous&Terrestrial" ~ "2"
) 

## defining transition matrix for only two states
twoState_q_ARD <- matrix(byrow = T, nrow = 2, ncol = 2, data = c(0,1, 2,0))


library(castor)


tree<-td_eco$phy

make_tip_priors <- function(Nstates, tree, data) {
  bins2 <- Nstates
  tip_priors <- matrix(1e-08/(Nstates - 1), nrow=length(tree$tip.label), ncol=bins2)
  fitdat <- data
  
  for(i in 1:nrow(tip_priors)){
    if(class(fitdat)=="character"){
      
      .states <- as.numeric(strsplit(fitdat[i], "&")[[1]])
      
    } else{
      .states <- fitdat[i]
    }
    dens <- (1 - 1e-08)/length(.states)
    tip_priors[i, .states] <- dens
  }
  return(tip_priors)
}



andPriors_foss <- make_tip_priors(5, tree = tree, data = td_eco$dat$substrate4)

andPriors_terr <- make_tip_priors(5, tree = tree, data = td_eco$dat$substrate2)

andPriors <- make_tip_priors(5, tree = tree, data = td_eco$dat$substrate3)

andPriors_twostate <- make_tip_priors(2, tree = tree, data = td_eco$dat$substrate5)


castorAIC <- function(mat, loglik){
  uh <- length(unique(mat[which(!is.na(mat))]))
  yAIC = 2*uh - 2*loglik
  return(yAIC)
}


twostate_fit <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors_twostate, Nstates = 2, rate_model = twoState_q_ARD) 
twostate_fit$AIC

foss_fivestate_fit <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors_foss, Nstates = 5, rate_model = fiveState_q_restricted_FASA) 
foss_fivefit_and <- castor::hsp_mk_model(tree, tip_states = NULL, tip_priors = andPriors_foss, Nstates = 5, rate_model = fiveState_q_ARD) 

foss_fivestate_ARD <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors_foss, Nstates = 5, rate_model = fiveState_q_ARD)
foss_fivestate_SYM <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors_foss, Nstates = 5, rate_model = fiveState_q_SYM)

foss_fivestate_restricted <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors_foss, Nstates = 5, rate_model = fiveState_q_restricted)
foss_fivestate_restricted_restricted <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors_foss, Nstates = 5, rate_model = fiveState_q_restricted_restricted)


terr_fivestate_fit <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors_terr, Nstates = 5, rate_model = fiveState_q_restricted_FASA) 
terr_fivefit_and <- castor::hsp_mk_model(tree, tip_states = NULL, tip_priors = andPriors_terr, Nstates = 5, rate_model = fiveState_q_ARD) 

terr_fivestate_ARD <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors_terr, Nstates = 5, rate_model = fiveState_q_ARD)
terr_fivestate_SYM <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors_terr, Nstates = 5, rate_model = fiveState_q_SYM)

terr_fivestate_restricted <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors_terr, Nstates = 5, rate_model = fiveState_q_restricted)
terr_fivestate_restricted_restricted <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors_terr, Nstates = 5, rate_model = fiveState_q_restricted_restricted)


fivestate_fit <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors, Nstates = 5, rate_model = fiveState_q_restricted_FASA) 
fivefit_and <- castor::hsp_mk_model(tree, tip_states = NULL, tip_priors = andPriors, Nstates = 5, rate_model = fiveState_q_ARD) 

fivestate_ARD <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors, Nstates = 5, rate_model = fiveState_q_ARD)
fivestate_SYM <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors, Nstates = 5, rate_model = fiveState_q_SYM)

fivestate_restricted <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors, Nstates = 5, rate_model = fiveState_q_restricted)
fivestate_restricted_restricted <- castor::asr_mk_model(tree, tip_states = NULL, tip_priors = andPriors, Nstates = 5, rate_model = fiveState_q_restricted_restricted)



fivestate_SYM$AIC
fivestate_fit$AIC
fivestate_ARD$AIC
fivestate_restricted$AIC
fivestate_restricted_restricted$AIC
castorAIC(mat = fivefit_and$transition_matrix, loglik = fivefit_and$loglikelihood)


foss_fivestate_SYM$AIC
foss_fivestate_fit$AIC
foss_fivestate_ARD$AIC
foss_fivestate_restricted$AIC
foss_fivestate_restricted_restricted$AIC
castorAIC(mat = foss_fivefit_and$transition_matrix, loglik = foss_fivefit_and$loglikelihood)

terr_fivestate_SYM$AIC
terr_fivestate_fit$AIC
terr_fivestate_ARD$AIC
terr_fivestate_restricted$AIC
terr_fivestate_restricted_restricted$AIC
castorAIC(mat = terr_fivefit_and$transition_matrix, loglik = terr_fivefit_and$loglikelihood)

newq <- fivestate_ARD$transition_matrix
terr_newq <- terr_fivestate_ARD$transition_matrix
foss_newq <- foss_fivestate_ARD$transition_matrix
td_eco$dat$species<-td_eco$phy$tip.label
two_newq <- twostate_fit$transition_matrix




tree$tip.label<-1:length(td_eco$phy$tip.label)


names(td_eco$dat$substrate2)<-1:length(td_eco$phy$tip.label)
names(td_eco$dat$substrate3)<-1:length(td_eco$phy$tip.label)
names(td_eco$dat$substrate4)<-1:length(td_eco$phy$tip.label)
names(td_eco$dat$substrate5)<-1:length(td_eco$phy$tip.label)

colnames(andPriors)<-1:5
rownames(andPriors)<-1:length(td_eco$phy$tip.label)

colnames(andPriors_foss)<-1:5
rownames(andPriors_foss)<-1:length(td_eco$phy$tip.label)

colnames(andPriors_terr)<-1:5
rownames(andPriors_terr)<-1:length(td_eco$phy$tip.label)

colnames(andPriors_twostate)<-1:2
rownames(andPriors_twostate)<-1:length(td_eco$phy$tip.label)


rownames(newq)<-1:5
colnames(newq)<-1:5

rownames(foss_newq)<-1:5
colnames(foss_newq)<-1:5

rownames(terr_newq)<-1:5
colnames(terr_newq)<-1:5

rownames(two_newq)<-1:2
colnames(two_newq)<-1:2



cols<-viridis::viridis(5)
names(cols)<-1:5

mtrees<-make.simmap(tree, andPriors,  Q=newq, nsim = 100)
foss_mtrees<-make.simmap(tree, andPriors_foss,  Q=foss_newq, nsim = 100)
terr_mtrees<-make.simmap(tree, andPriors_terr,  Q=terr_newq, nsim = 100)
two_mtrees<-make.simmap(tree, andPriors_twostate,  Q=two_newq, nsim = 100)


library(doParallel)

n.cores <- parallel::detectCores() - 8

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK")

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()






eco_map_fits<-list()

eco_map_fits<-foreach (i=1:length(mtrees), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp', 'phytools'), .combine = 'c') %dopar% {
  

#index<-mtrees[[i]]$tip.label[which(td_eco$phy$tip.label %in% mean_list$Genus_Species)]
#phylo1<-phytools::drop.tip.simmap(mtrees[[i]], tip = setdiff(mtrees[[i]]$tip.label, index))
phylo1<-phytools::map.to.singleton(mtrees[[i]])


phylo1$tip.label<-1:length(phylo1$tip.label)
phylo1<-PCMBase::PCMTree(phylo1)


vector1<-phylo1$edge[,2]
present <- rep(FALSE, length(vector1))
for (num in 1:length(vector1)) {
  if (vector1[num] >= 1 && vector1[num] <= length(phylo1$tip.label)) {
    present[num] <- TRUE
  }
}
  rows<-which(present==TRUE)
  
  phylo1$tipstate<-names(phylo1$edge.length)[rows]
  
  phylo1$edge[rows,]
  
  phylo1$nodes<-phylo1$tipstate
  names(phylo1$nodes)<-phylo1$edge[,1][rows]
  phylo1$nodes
  
PCMBase::PCMTreeSetPartRegimes(
    phylo1, 
    part.regime = phylo1$nodes[unique(names(phylo1$nodes))], 
    setPartition = TRUE)
  
  
  phylo1


}

terr_eco_map_fits<-list()

terr_eco_map_fits<-foreach (i=1:length(terr_mtrees), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp', 'phytools'), .combine = 'c') %dopar% {
  
  
  #index<-mtrees[[i]]$tip.label[which(td_eco$phy$tip.label %in% mean_list$Genus_Species)]
  #phylo1<-phytools::drop.tip.simmap(mtrees[[i]], tip = setdiff(mtrees[[i]]$tip.label, index))
  phylo1<-phytools::map.to.singleton(terr_mtrees[[i]])
  
  
  phylo1$tip.label<-1:length(phylo1$tip.label)
  phylo1<-PCMBase::PCMTree(phylo1)
  
  
  vector1<-phylo1$edge[,2]
  present <- rep(FALSE, length(vector1))
  for (num in 1:length(vector1)) {
    if (vector1[num] >= 1 && vector1[num] <= length(phylo1$tip.label)) {
      present[num] <- TRUE
    }
  }
  rows<-which(present==TRUE)
  
  phylo1$tipstate<-names(phylo1$edge.length)[rows]
  
  phylo1$edge[rows,]
  
  phylo1$nodes<-phylo1$tipstate
  names(phylo1$nodes)<-phylo1$edge[,1][rows]
  phylo1$nodes
  
  PCMBase::PCMTreeSetPartRegimes(
    phylo1, 
    part.regime = phylo1$nodes[unique(names(phylo1$nodes))], 
    setPartition = TRUE)
  
  
  phylo1
  
  
}

foss_eco_map_fits<-list()

foss_eco_map_fits<-foreach (i=1:length(foss_mtrees), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp', 'phytools'), .combine = 'c') %dopar% {
  
  
  #index<-mtrees[[i]]$tip.label[which(td_eco$phy$tip.label %in% mean_list$Genus_Species)]
  #phylo1<-phytools::drop.tip.simmap(mtrees[[i]], tip = setdiff(mtrees[[i]]$tip.label, index))
  phylo1<-phytools::map.to.singleton(foss_mtrees[[i]])
  
  
  phylo1$tip.label<-1:length(phylo1$tip.label)
  phylo1<-PCMBase::PCMTree(phylo1)
  
  
  vector1<-phylo1$edge[,2]
  present <- rep(FALSE, length(vector1))
  for (num in 1:length(vector1)) {
    if (vector1[num] >= 1 && vector1[num] <= length(phylo1$tip.label)) {
      present[num] <- TRUE
    }
  }
  rows<-which(present==TRUE)
  
  phylo1$tipstate<-names(phylo1$edge.length)[rows]
  
  phylo1$edge[rows,]
  
  phylo1$nodes<-phylo1$tipstate
  names(phylo1$nodes)<-phylo1$edge[,1][rows]
  phylo1$nodes
  
  PCMBase::PCMTreeSetPartRegimes(
    phylo1, 
    part.regime = phylo1$nodes[unique(names(phylo1$nodes))], 
    setPartition = TRUE)
  
  
  phylo1
  
  
}


two_eco_map_fits<-list()

two_eco_map_fits<-foreach (i=1:length(two_mtrees), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp', 'phytools'), .combine = 'c') %dopar% {
  
  
  #index<-mtrees[[i]]$tip.label[which(td_eco$phy$tip.label %in% mean_list$Genus_Species)]
  #phylo1<-phytools::drop.tip.simmap(mtrees[[i]], tip = setdiff(mtrees[[i]]$tip.label, index))
  phylo1<-phytools::map.to.singleton(two_mtrees[[i]])
  
  
  phylo1$tip.label<-1:length(phylo1$tip.label)
  phylo1<-PCMBase::PCMTree(phylo1)
  
  
  vector1<-phylo1$edge[,2]
  present <- rep(FALSE, length(vector1))
  for (num in 1:length(vector1)) {
    if (vector1[num] >= 1 && vector1[num] <= length(phylo1$tip.label)) {
      present[num] <- TRUE
    }
  }
  rows<-which(present==TRUE)
  
  phylo1$tipstate<-names(phylo1$edge.length)[rows]
  
  phylo1$edge[rows,]
  
  phylo1$nodes<-phylo1$tipstate
  names(phylo1$nodes)<-phylo1$edge[,1][rows]
  phylo1$nodes
  
  PCMBase::PCMTreeSetPartRegimes(
    phylo1, 
    part.regime = phylo1$nodes[unique(names(phylo1$nodes))], 
    setPartition = TRUE)
  
  
  phylo1
  
  
}





X<-td_eco$dat[1:2]
rownames(X)<-td_eco$phy$tip.label


X<-t(X)


source('NewDefineParameterLimits.R', local=FALSE)
# threshold for 0 eigenvalues
options(PCMBase.Threshold.EV = 0.0000000001)
getOption("PCMBase.Threshold.EV")


# threshold for 0 eigenvalues
options(PCMBase.Threshold.EV = 0.000000001)
getOption("PCMBase.Threshold.EV")

pcmerror<-readRDS("pcmerror_ext.RDS")


modelTrueTypeMappingBM <- PCMBase::MixedGaussian(
  k = 2,
  modelTypes = PCMBase::MGPMDefaultModelTypes(),
  mapping = c(2, 2, 2, 2, 2),
  X0 = structure(
    0, class = c("VectorParameter",
                 "_Global"), description = "trait values at the root"), 
  Sigmae_x = structure(
    0, class = c("MatrixParameter", "_Omitted",
                 "_Global"),
    description =
      "Upper triangular Choleski factor of the non-phylogenetic variance-covariance"))

generatePCMModelsFunction <- function() {
  # make results reproducible
  # set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  PCMBase::PCMGenerateModelTypes()
  fileName <- 'NewDefineParameterLimits.R'
  codeDefineLimits <- readChar(fileName, file.info(fileName)$size)
  eval(parse(text = codeDefineLimits), .GlobalEnv)
}



eco_fits_BM<-foreach (i=1:length(eco_map_fits), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp'), .combine= 'c') %dopar% {
  

fitMGPMTrueTypeMapping <- PCMFit::PCMFit(
  model = modelTrueTypeMappingBM, 
  tree = eco_map_fits[[i]], 
  X = X, 
  SE = pcmerror, 
  metaI = PCMBaseCpp::PCMInfoCpp)


list(RetrieveBestModel(fitMGPMTrueTypeMapping))

}

terr_eco_fits_BM<-foreach (i=1:length(terr_eco_map_fits), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp'), .combine= 'c') %dopar% {
  
  
  fitMGPMTrueTypeMapping <- PCMFit::PCMFit(
    model = modelTrueTypeMappingBM, 
    tree = terr_eco_map_fits[[i]], 
    X = X, 
    SE = pcmerror, 
    metaI = PCMBaseCpp::PCMInfoCpp)
  
  
  list(RetrieveBestModel(fitMGPMTrueTypeMapping))
  
}

foss_eco_fits_BM<-foreach (i=1:length(foss_eco_map_fits), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp'), .combine= 'c') %dopar% {
  
  
  fitMGPMTrueTypeMapping <- PCMFit::PCMFit(
    model = modelTrueTypeMappingBM, 
    tree = foss_eco_map_fits[[i]], 
    X = X, 
    SE = pcmerror, 
    metaI = PCMBaseCpp::PCMInfoCpp)
  
  
  list(RetrieveBestModel(fitMGPMTrueTypeMapping))
  
}


modelTrueTypeMappingOU <- PCMBase::MixedGaussian(
  k = 2,
  modelTypes = PCMBase::MGPMDefaultModelTypes(),
  mapping = c(5, 5, 5, 5, 5),
  X0 = structure(
    0, class = c("VectorParameter",
                 "_Global"), description = "trait values at the root"), 
  Sigmae_x = structure(
    0, class = c("MatrixParameter", "_Omitted",
                 "_Global"),
    description =
      "Upper triangular Choleski factor of the non-phylogenetic variance-covariance"))

generatePCMModelsFunction <- function() {
  # make results reproducible
  # set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  PCMBase::PCMGenerateModelTypes()
  fileName <- 'NewDefineParameterLimits.R'
  codeDefineLimits <- readChar(fileName, file.info(fileName)$size)
  eval(parse(text = codeDefineLimits), .GlobalEnv)
}



eco_fits_OU<-foreach (i=1:length(eco_map_fits), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp'), .combine= 'c') %dopar% {
  
  
  fitMGPMTrueTypeMapping <- PCMFit::PCMFit(
    model = modelTrueTypeMappingOU, 
    tree = eco_map_fits[[i]], 
    X = X, 
    SE = pcmerror, 
    metaI = PCMBaseCpp::PCMInfoCpp)
  
  
  list(RetrieveBestModel(fitMGPMTrueTypeMapping))
  
}




terr_eco_fits_OU<-foreach (i=1:length(terr_eco_map_fits), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp'), .combine= 'c') %dopar% {
  
  
  fitMGPMTrueTypeMapping <- PCMFit::PCMFit(
    model = modelTrueTypeMappingOU, 
    tree = terr_eco_map_fits[[i]], 
    X = X, 
    SE = pcmerror, 
    metaI = PCMBaseCpp::PCMInfoCpp)
  
  
  list(RetrieveBestModel(fitMGPMTrueTypeMapping))
  
}




foss_eco_fits_OU<-foreach (i=1:length(foss_eco_map_fits), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp'), .combine= 'c') %dopar% {
  
  
  fitMGPMTrueTypeMapping <- PCMFit::PCMFit(
    model = modelTrueTypeMappingOU, 
    tree = foss_eco_map_fits[[i]], 
    X = X, 
    SE = pcmerror, 
    metaI = PCMBaseCpp::PCMInfoCpp)
  
  
  list(RetrieveBestModel(fitMGPMTrueTypeMapping))
  
}




modelTrueTypeMapping2BM <- PCMBase::MixedGaussian(
  k = 2,
  modelTypes = PCMBase::MGPMDefaultModelTypes(),
  mapping = c(2, 2),
  X0 = structure(
    0, class = c("VectorParameter",
                 "_Global"), description = "trait values at the root"), 
  Sigmae_x = structure(
    0, class = c("MatrixParameter", "_Omitted",
                 "_Global"),
    description =
      "Upper triangular Choleski factor of the non-phylogenetic variance-covariance"))

generatePCMModelsFunction <- function() {
  # make results reproducible
  # set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  PCMBase::PCMGenerateModelTypes()
  fileName <- 'NewDefineParameterLimits.R'
  codeDefineLimits <- readChar(fileName, file.info(fileName)$size)
  eval(parse(text = codeDefineLimits), .GlobalEnv)
}

two_eco_fits_BM<-foreach (i=1:length(two_eco_map_fits), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp'), .combine= 'c') %dopar% {
  
  
  
  fitMGPMTrueTypeMapping <- PCMFit::PCMFit(
    model = modelTrueTypeMapping2BM, 
    tree = two_eco_map_fits[[i]], 
    X = X, 
    SE = pcmerror, 
    metaI = PCMBaseCpp::PCMInfoCpp)
  
  
  list(RetrieveBestModel(fitMGPMTrueTypeMapping))
  
}

modelTrueTypeMapping2OU <- PCMBase::MixedGaussian(
  k = 2,
  modelTypes = PCMBase::MGPMDefaultModelTypes(),
  mapping = c(5, 5),
  X0 = structure(
    0, class = c("VectorParameter",
                 "_Global"), description = "trait values at the root"), 
  Sigmae_x = structure(
    0, class = c("MatrixParameter", "_Omitted",
                 "_Global"),
    description =
      "Upper triangular Choleski factor of the non-phylogenetic variance-covariance"))

generatePCMModelsFunction <- function() {
  # make results reproducible
  # set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  PCMBase::PCMGenerateModelTypes()
  fileName <- 'NewDefineParameterLimits.R'
  codeDefineLimits <- readChar(fileName, file.info(fileName)$size)
  eval(parse(text = codeDefineLimits), .GlobalEnv)
}


two_eco_fits_OU<-foreach (i=1:length(two_eco_map_fits), .packages = c('PCMFit', 'PCMBase', 'PCMBaseCpp'), .combine= 'c') %dopar% {
  
  
  fitMGPMTrueTypeMapping <- PCMFit::PCMFit(
    model = modelTrueTypeMapping2OU, 
    tree = two_eco_map_fits[[i]], 
    X = X, 
    SE = pcmerror, 
    metaI = PCMBaseCpp::PCMInfoCpp)
  
  
  list(RetrieveBestModel(fitMGPMTrueTypeMapping))
  
}


stopCluster(my.cluster)



##everything together except fossorial
saveRDS(two_eco_fits_BM, "only_foss_eco_fits_BM.RDS")
saveRDS(two_eco_fits_OU, "only_foss_eco_fits_OU.RDS")


##everything fossorial&terrestrial as fossorial&terrestrial
saveRDS(eco_fits_BM, "eco_fits_BM.RDS")
saveRDS(eco_fits_OU, "eco_fits_OU.RDS")


##everything fossorial&terrestrial as terrestrial
saveRDS(terr_eco_fits_BM, "terr_eco_fits_BM.RDS")
saveRDS(terr_eco_fits_OU, "terr_eco_fits_OU.RDS")


##everything  fossorial&terrestrial as fossorial
saveRDS(foss_eco_fits_BM, "foss_eco_fits_BM.RDS")
saveRDS(foss_eco_fits_OU, "foss_eco_fits_OU.RDS")


saveRDS(eco_map_fits, "eco_map_fits.RDS")
eco_map_fits<-readRDS("eco_map_fits.RDS")


library(viridis)
cols<-viridis(5)
names(cols)<-1:5
plot(mtrees[[47]], col=cols, ftype="off")
legend(legend = c("Terrestrial", "Arboreal", "Saxicolous", "Fossorial", "Semi-aquatic"), x=0, y= 275, col = cols, bty = 'n', lwd = 2)
plot(summary(mtrees), colors=cols, ftype="off", cex=0.25, type="fan")
legend(legend = c("Terrestrial", "Arboreal", "Saxicolous", "Fossorial", "Semi-aquatic"), x=0, y= 275, col = cols, bty = 'n', pch = 19)

col.morph<-setNames(viridis::viridis(5),c("Terrestrial", "Arboreal", "Saxicolous", "Fossorial", "Semi-aquatic"))

plotTree(td_eco$phy,ftype="off",lwd=1, type="fan")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)

pp<-summary(foss_mtrees)$ace[1:td_eco$phy$Nnode,]
ii<-order(rowSums(pp^2),decreasing=TRUE)
h<-max(nodeHeights(td_eco$phy))
for(i in ii){
  plotrix::floating.pie(obj$xx[i+Ntip(tree)],
                        obj$yy[i+Ntip(tree)],
                        radius=if(max(pp[i,])<=0.90) 0.025*h else 0.0*h,
                        x=pp[i,],col=col.morph,border="black")
}
legend(legend = c("Terrestrial", "Arboreal", "Saxicolous", "Fossorial", "Semi-aquatic"), x=0, y= 275, col = cols, bty = 'n', pch = 19)


rownames(foss_newq)<-c("T", "A", "S", "F", "SA")
colnames(foss_newq)<-c("T", "A", "S", "F", "SA")
phytools::plot.Qmatrix(foss_newq, width=TRUE)


