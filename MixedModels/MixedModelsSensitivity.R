##### Preparing working directory #####

# clear workspace 
rm(list=ls()) 

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

# load tree and data
ext_td<-readRDS("ext2_td.RDS") # new tree and updated data
ext_td$dat$species<-ext_td$phy$tip.label # for combining data with other trees later
dat<-ext_td$dat # making data independent of trees
mm<-readRDS("final_mm.RDS") # individual data for creating error objects
trees<-read.tree("pseudoposterior.100.trees") # pseudo posterior trees from Title et al 2024
set.seed(20) # set seed for reproducibility 
rand<-sample(1:100, size = 10) # random sample of size 10 
tree<-trees[rand] # 10 trees randomly sampled from 100 posterior trees



# threshold for 0 eigenvalues
options(PCMBase.Threshold.EV = 0.00000001)
getOption("PCMBase.Threshold.EV")

# make lists to contain error data and trees
error<-list()
X<-list()



##### for loop for creating error objects for 10 tree topologies #####
for (t in 1:10) {
  
  #is.ultrametric(tree[[t]])
  #tree[[i]]<-nnls.tree(cophenetic(tree[[t]]),tree[[t]],
  #                rooted=TRUE,trace=0)
  #is.ultrametric(tree[[t]])
  
  # make a treedata object which combines data with tree on each for loop iteration
  td<-make.treedata(tree[[t]], dat)
  
  # make an index to connect tip labels with individual data names
  spname<-unique(td$dat$names_match)
  
  # subset individual data by spname, calculating variance, or assign NA if n <= to 2 
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
  
  # create vector of sample sizes
  x<-td$dat$n
  x<-as.numeric(x)
  # assign NAs when n <= 2
  x<-case_when(x==1 ~ NA_real_,
               x==2 ~ NA_real_,
               TRUE ~ x)
  
  # remove NA values from vector
  x<-x[!is.na(x)]
  
  # remove NA values from list of matrices 
  cov<-plyr::llply(listcov, .inform = TRUE, function(m){
    if (all(is.na(m))==TRUE) {
      
    } else {
      m
    }
  }  
  )
  
  # remove null values 
  cov<-cov[!sapply(cov,is.null)]
  
  # remove matrices with missing values 
  cov<-plyr::llply(cov, .inform = TRUE, function(m){
    if (any(is.na(m))==TRUE) {
      
    } else {
      m
    }
  }  
  )
  
  # remove values from vector which had missing values in list 
  x <- x[-c(which(sapply(cov, function(q){
    is.null(q)
  })==TRUE))]
  
  # remove null values again
  cov<-cov[!sapply(cov,is.null)]
  # condense
  cov<-abind::abind(cov, along = 3)
  
  # check for null values
  count<-0
  sum(sapply(cov, function(p){
    if (is.null(p)==TRUE) {
      count<-count+1
    } else {
      count
    }
  }
  ))
  
  # calculate weighted average 
  covwa<-apply(cov, 1:2, FUN=weighted.mean, w=x)
  
  # add weighted average in places with missing values 
  VCV<-plyr::llply(listcov, .inform = TRUE, function(o){
    if (any(is.na(o))==TRUE) {
      covwa
    } else {
      o
    }
  }  
  )
  
  # create lists for VCV and errors
  VCVn<-vector(mode = "list", length = length(VCV))
  pcmerror<-vector(mode = "list", length = length(VCV))
  
  # for loop to create VCV matrices and format them for PCMFitMixed
  for (o in 1:length(VCV)) {
    
    VCVn[[o]]<-VCV[[o]]/td$dat$n[o]
    pcmerror[[o]]<-UpperTriFactor(VCVn[[o]])
    
  }
  
  # condense to arrays 
  VCVn<-abind::abind(VCVn, along = 3)
  pcmerror<-abind::abind(pcmerror, along = 3)
  
  # make list for each tree iteration
  error[[t]]<-pcmerror
  
  # format trees and data for PCMFitMixed
  tree[[t]]<-td$phy
  X[[t]]<-td$dat[,1:2]
  rownames(X[[t]])<-td$phy$tip.label
  tree[[t]]<-PCMTree(tree[[t]])
  X[[t]]<-t(X[[t]])
  
}


# set cores for running in parallel
doParallel::registerDoParallel(cores = 21)
parallel::mcaffinity(24:44)


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

  fitMIXED[[j]] <- PCMFitMixed(modelTypes = MGPMDefaultModelTypes(), SE = error[[j]],
  X = X[[j]], tree = tree[[j]], metaIFun = PCMInfoCpp,
  generatePCMModelsFun = generatePCMModelsFunction, 
  maxNumRoundRobins = , maxNumPartitionsInRoundRobins = 2,
  tableFitsPrev = tableFitsPrev,
  prefixFiles = prefixFiles,
  doParallel = TRUE, minCladeSizes = 5)
  
  
  #saveRDS(fitMIXED, file = "fitMixed.RDS")
  #save.image(paste0("fitmixed", j, ".RData"))
  
} 



#Error in mapping[[which.min(score)]] : 
#attempt to select less than one element in get1index

# investigate tree[[i]]
# try try/while


