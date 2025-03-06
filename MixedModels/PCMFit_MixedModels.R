##### Preparing working directory #####
setwd("MixedModels")
# clear workspace 
rm(list=ls()) 


# load tree and data
ext_td<-readRDS("ext2_td.RDS") # new tree and updated data
mm<-readRDS("final_mm.RDS") # individual data for creating error objects

# load packages
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

# threshold for 0 eigenvalues
options(PCMBase.Threshold.EV = 0.00000001)
getOption("PCMBase.Threshold.EV")




spname<-unique(ext_td$dat$names_match)


listcov<-
  plyr::llply(spname, .inform = TRUE, function(i){
    subsetdat<-subset(mm,Genus_Species==i)
    n<-dim(subsetdat)[1]
    if (n>2) {
      dplyr::select(subsetdat, SVL, TLL) %>% var
    } else {
      
      matrix(data = NA, nrow = 2, ncol = 2)
    }
  })


x<-ext_td$dat$n
x<-as.numeric(x)
x<-case_when(x==1 ~ NA_real_,
             x==2 ~ NA_real_,
             TRUE ~ x)


x<-x[!is.na(x)]

cov<-plyr::llply(listcov, .inform = TRUE, function(m){
  if (all(is.na(m))==TRUE) {
    
  } else {
    m
  }
}  
)

cov<-cov[!sapply(cov,is.null)]

cov<-plyr::llply(cov, .inform = TRUE, function(m){
  if (any(is.na(m))==TRUE) {
    
  } else {
    m
  }
}  
)

x <- x[-c(which(sapply(cov, function(q){
  is.null(q)
})==TRUE))]



cov<-cov[!sapply(cov,is.null)]
cov<-abind::abind(cov, along = 3)

count<-0
sum(sapply(cov, function(p){
  if (is.null(p)==TRUE) {
    count<-count+1
  } else {
    count
  }
}
))


covwa<-apply(cov, 1:2, FUN=weighted.mean, w=x)


VCV<-plyr::llply(listcov, .inform = TRUE, function(o){
  if (any(is.na(o))==TRUE) {
    covwa
  } else {
    o
  }
}  
)

VCVn<-vector(mode = "list", length = length(VCV))
pcmerror<-vector(mode = "list", length = length(VCV))

for (o in 1:length(VCV)) {
  
  VCVn[[o]]<-VCV[[o]]/ext_td$dat$n[o]
  pcmerror[[o]]<-UpperTriFactor(VCVn[[o]])
  
}

VCVn<-abind::abind(VCVn, along = 3)
pcmerror<-abind::abind(pcmerror, along = 3)



tree<-ext_td$phy
X<-ext_td$dat[,1:2]
rownames(X)<-tree$tip.label


tree<-PCMTree(tree)
class(tree)

X<-t(X)


# set cores for running in parallel
doParallel::registerDoParallel(cores = 21)
parallel::mcaffinity(24:44)



##### fitting BM and OU models #####


# set model types 
modelBM <- PCM(
  PCMDefaultModelTypes()["B"], modelTypes = PCMDefaultModelTypes(), k = 2)

modelOU <- PCM(
  PCMDefaultModelTypes()["F"], modelTypes = PCMDefaultModelTypes(), k = 2)

# fit models
pcmerrorfitBM <- PCMFit(model = modelBM, tree = tree, X = X, SE = pcmerror,
                        metaI = PCMBaseCpp::PCMInfoCpp)

pcmerrorfitOU <- PCMFit(model = modelOU, tree = tree, X = X, SE = pcmerror,
                        metaI = PCMBaseCpp::PCMInfoCpp)







##### for loop for running PCMFitMixed #####

# generate model types
generatePCMModelsFunction <- function() {
  # make results reproducible
  #set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  PCMGenerateModelTypes()
  fileName <- 'NewDefineParameterLimits.R'
  codeDefineLimits <- readChar(fileName, file.info(fileName)$size)
  eval(parse(text = codeDefineLimits), .GlobalEnv)
}

# set number of runs (1 for each tree)
nruns<-10

# create list to contain each fit
fitMIXED<-vector("list",nruns)

# for loop to run PCMFitMixed 
for(j in 1:length(fitMIXED)) {
  
  #Fitting unknown shifts-----
  
  prefixFiles = paste0("mixed_",j)
  
  currentResultFile <- paste0("Current_",j,"_", prefixFiles, ".RData")
  if(file.exists(currentResultFile)) {
    load(currentResultFile)
    tableFitsPrev <- listResults$tableFits
  } else {
    tableFitsPrev <- NULL
  }
  
  fitMIXED[[j]] <- PCMFitMixed(modelTypes = MGPMDefaultModelTypes(), SE = pcmerror,
                               X = X, tree = tree, metaIFun = PCMInfoCpp,
                               generatePCMModelsFun = generatePCMModelsFunction, 
                               maxNumRoundRobins = , maxNumPartitionsInRoundRobins = 2,
                               tableFitsPrev = tableFitsPrev,
                               prefixFiles = prefixFiles,
                               doParallel = TRUE, minCladeSizes = 5)
  
  
  #saveRDS(fitMIXED, file = "fitMixed.RDS")
  #save.image(paste0("fitmixed", j, ".RData"))
  
} 

listModels <- list(
  RetrieveBestModel(pcmerrorfitBM), 
  RetrieveBestModel(pcmerrorfitOU))

dtSummary <- data.table(
  model = 1:2,
  p = sapply(listModels, PCMParamCount),
  logLik = sapply(listModels, logLik), 
  AIC = sapply(listModels, AIC))
knitr::kable(dtSummary)
