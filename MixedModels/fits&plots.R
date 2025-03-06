# clear workspace 
rm(list=ls()) 

## load data and packages
library(PCMFit)
library(PCMBase)
library(phytools)
library(data.table)

## set threshold for 0 eigenvalues
options(PCMBase.Threshold.EV = 0.00000001)


## find best model
ext_fitMIXED<-readRDS("ext_fitMIXED.RDS")

listModels_ext <- list(
  RetrieveBestFitScore(ext_fitMIXED[[1]])$inferredModel, 
  RetrieveBestFitScore(ext_fitMIXED[[2]])$inferredModel,
  RetrieveBestFitScore(ext_fitMIXED[[3]])$inferredModel,
  RetrieveBestFitScore(ext_fitMIXED[[4]])$inferredModel,
  RetrieveBestFitScore(ext_fitMIXED[[5]])$inferredModel,
  RetrieveBestFitScore(ext_fitMIXED[[6]])$inferredModel, 
  RetrieveBestFitScore(ext_fitMIXED[[7]])$inferredModel,
  RetrieveBestFitScore(ext_fitMIXED[[8]])$inferredModel,
  RetrieveBestFitScore(ext_fitMIXED[[9]])$inferredModel,
  RetrieveBestFitScore(ext_fitMIXED[[10]])$inferredModel)

ext_dtSummary <- data.table(
  model = 1:10,
  p = sapply(listModels_ext, PCMParamCount),
  logLik = sapply(listModels_ext, logLik), 
  AIC = sapply(listModels_ext, AIC))
knitr::kable(ext_dtSummary)

ext_8<-RetrieveBestFitScore(ext_fitMIXED[[8]])

## calculate slope and intercept
Vy_list <- vector(mode = "list", length = length(ext_8$inferredModel)-1)
theta_list<- vector(mode= "list", length = length(ext_8$inferredModel)-1)

## stationary variance function from MVMorph
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
for (i in 2:(length(ext_8$inferredModel))) {
  if (length(ext_8$inferredModel[[i]])==3) {
    
    tempsig<-tcrossprod(ext_8$inferredModel[[i]]$Sigma_x[1:2, 1:2,])
    Vy_list[[i-1]]<-StationaryVariance(alpha = PCMApplyTransformation(ext_8$inferredModel[[i]])$H[1:2,1:2,], sigma = tempsig)
    theta_list[[i-1]]<-ext_8$inferredModel[[i]]$Theta[1:2]
    
  } else {
    
    Vy_list[[i-1]]<-tcrossprod(ext_8$inferredModel[[i]]$Sigma_x[1:2, 1:2,]) # using sig2 matrices for Vy for BM models
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

slopeint<-matrix(data = NA, nrow = length(ext_8$inferredModel)-1, ncol = 2)
rownames(slopeint)<-1:(length(ext_8$inferredModel)-1)
colnames(slopeint)<-c("slope", "int")
slopeint[,"slope"]<-sapply(Vy_list, FUN = Vy_slope)
slopeint<-as.data.frame(slopeint)


thetas<-matrix(data = NA, nrow = length(ext_8$inferredModel)-1, ncol = 2)
rownames(thetas)<-1:(length(ext_8$inferredModel)-1)

for (i in 1:(length(ext_8$inferredModel)-1)) {
  if(is.nan(theta_list[[i]][1])) {
    slopeint[i,2]<-ext_8$inferredModel$X0[2] - slopeint[i,1]*ext_8$inferredModel$X0[1]
  } else {
    thetax<-theta_list[[i]][1]
    thetay<-theta_list[[i]][2]
    slopeint[i,2]<-thetay-thetax*slopeint[i,1]
    thetas[i]<-thetay/thetax
  }
  
}

## Plot tree
PCMTreePlot(ext_8$tree, layout="fan")



## Plot allometry
library(ggplot2)

slope<-slopeint[,"slope"]
names(slope)<-1:7
int<-slopeint[,"int"]
names(int)<-1:7

cols<-PCMColorPalette(7, 1:7)

xplot<-PCMPlotTraitData2D(
  ext_8$X, 
  ext_8$tree, 
  scaleSizeWithTime = FALSE,
  numTimeFacets = 1) +
  geom_text(
    aes(x = x, y = y, label = NA, color = regime), 
    size=2, 
    hjust=2) +
  theme_bw() +
  theme(legend.position = "bottom")+
  theme(panel.grid = element_blank())+
  coord_fixed()+
  xlab(expression(paste("snout vent length (", italic("ln"), " mm)")))+
  ylab(expression(paste("hind limb length (", italic("ln"), " mm)"))) +
  geom_abline(slope = slope[1], intercept = int[1], color = cols[1], linewidth=2) +
  geom_abline(slope = slope[2], intercept = int[2], color = cols[2], linewidth=2) +
  geom_abline(slope = slope[3], intercept = int[3], color = cols[3], linewidth=2) +
  geom_abline(slope = slope[4], intercept = int[4], color = cols[4], linewidth=2) +
  geom_abline(slope = slope[5], intercept = int[5], color = cols[5], linewidth=2) +
  geom_abline(slope = slope[6], intercept = int[6], color = cols[6], linewidth=2) +
  geom_abline(slope = slope[7], intercept = int[7], color = cols[7], linewidth=2) # +
  


allometry_figs<-list()

for (i in 1:PCMNumRegimes(ext_8$inferredModel)) {
  
  df<-as.data.frame(t(ext_8$X[,PCMTreeGetTipsInRegime(tree=ext_8$tree,regime = i)]))
  
  allometry_figs[[i]] <- ggplot(data = df, aes(color = PCMColorPalette(7, 1:7)[i])) +
    geom_point(aes(x = SVL, y = TLL, color = PCMColorPalette(7, 1:7)[i]), color=PCMColorPalette(7, 1:7)[i]) + 
    theme_bw()+
    theme(panel.grid = element_blank())+
    coord_fixed()+
    xlab(element_blank())+
    ylab(element_blank()) +
    geom_abline(slope = slope[i], intercept = int[i], color=PCMColorPalette(7, 1:7)[i], size=2) +
    xlim(2.5, 6.5) +
    ylim(-0.5,6)
  
  
}



library(cowplot)


combined_plot <- plot_grid(
  plotlist = allometry_figs,  # the list of plots
  ncol = 4)+
  draw_label("snout-vent length (ln mm)", x= 0.5, y= 0, vjust=-0.5, angle= 0) +
  draw_label("hind limb length (ln mm)", x= 0, y= 0.5, vjust= 1.5, angle=90)

# Display the combined plot
print(combined_plot)




## PGLS 

df<-as.data.frame(t(ext_8$X))

vector1<-ext_8$tree$edge[,2]
present <- rep(FALSE, length(vector1))
for (num in 1:length(vector1)) {
  if (vector1[num] >= 1 && vector1[num] <= length(ext_8$tree$tip.label)) {
    present[num] <- TRUE
  }
}

rows<-which(present==TRUE)

df$regimes<-ext_8$tree$edge.part[rows]

df$regimes<-as.factor(df$regimes)

rownames(df)<-ext_8$tree$tip.label
df$Names<-ext_8$tree$tip.label


library(phylolm)

mod1<-phylolm(df$TLL~df$SVL*df$regimes,phy=ext_8$tree, data = df)


estimates <- summary(mod1)$coefficients[, "Estimate"]
std_errors <- summary(mod1)$coefficients[, "StdErr"]

lower_CI <- estimates - (1.96 * std_errors)
upper_CI <- estimates + (1.96 * std_errors)

results <- data.frame(
  Term = rownames(summary(mod1)$coefficients),
  Estimate = estimates,
  Std_Error = std_errors,
  Lower_95_CI = lower_CI,
  Upper_95_CI = upper_CI
)

intercept <- results$Estimate[1]  # First term is usually the intercept
base_slope_term <- grep(".*:.*", results$Term, invert = TRUE, value = TRUE)[2]  # First non-intercept term
base_slope <- results$Estimate[results$Term == base_slope_term]

# Initialize new columns for adjusted values
results$Adjusted_Estimate <- results$Estimate
results$Adjusted_Lower_95_CI <- results$Lower_95_CI
results$Adjusted_Upper_95_CI <- results$Upper_95_CI

# Adjust intercepts for categorical groups
group_rows <- grep("^[^:]+$", results$Term)  # Terms without ":"
for (i in group_rows[-c(1,2)]) {  # Skip the first row (intercept)
  results$Adjusted_Estimate[i] <- results$Estimate[i] + intercept
  results$Adjusted_Lower_95_CI[i] <- results$Lower_95_CI[i] + intercept
  results$Adjusted_Upper_95_CI[i] <- results$Upper_95_CI[i] + intercept
}

# Adjust slopes for interaction terms
interaction_rows <- grep(":", results$Term)  # Find interaction terms
for (i in interaction_rows) {
  # Adjust slopes by adding base slope to interaction effect
  results$Adjusted_Estimate[i] <- results$Estimate[i] + base_slope
  results$Adjusted_Lower_95_CI[i] <- results$Lower_95_CI[i] + base_slope
  results$Adjusted_Upper_95_CI[i] <- results$Upper_95_CI[i] + base_slope
}


library(xtable)
print(xtable(cbind(results[,c("Adjusted_Estimate", "Std_Error")], summary(mod1)$coefficients[,"t.value"], summary(mod1)$coefficients[,"p.value"]), type='latex'))








X<-readRDS("X.RDS")
fitMIXED<-readRDS("fitMIXED.RDS")
fitMIXED<-lapply(fitMIXED, RetrieveBestFitScore)

tips_list<-list()
index_list<-list()

for (x in 1:length(fitMIXED)) {
  
  tips_list[[x]]<-fitMIXED[[x]]$tree$edge.part[which(fitMIXED[[x]]$tree$edge[,2] %in% 1:length(fitMIXED[[x]]$tree$tip.label)==TRUE)]
  names(tips_list[[x]])<-colnames(X[[x]])
  index_list[[x]]<-fitMIXED[[x]]$tree$part.regime
  
}




# Function to replace values based on the index
replace_values <- function(df, index) {
  rp<-as.numeric(index[as.character(df)])
  names(rp)<-names(df)
  return(rp)
}

tip_maps <- Map(replace_values, tips_list, index_list)




# Function to calculate mean across named vectors
mean_across_vectors <- function(vector_list) {
  # Identify common names across all vectors
  common_names <- Reduce(intersect, lapply(vector_list, names))
  
  # Extract only the values for the common names from each vector
  filtered_vectors <- lapply(vector_list, function(vec) vec[common_names])
  
  # Combine the filtered vectors into a matrix
  combined_matrix <- do.call(cbind, filtered_vectors)
  
  # Take the mean across the columns (i.e., across vectors)
  mean_vector <- rowMeans(combined_matrix)
  
  # Assign the names back to the mean vector
  names(mean_vector) <- common_names
  
  return(mean_vector)
}

# Apply the function to calculate the mean across vectors
mean_vector <- mean_across_vectors(tip_maps)

# View the result
mean_vector

# Function to calculate median across named vectors
median_across_vectors <- function(vector_list) {
  # Identify common names across all vectors
  common_names <- Reduce(intersect, lapply(vector_list, names))
  
  # Extract only the values for the common names from each vector
  filtered_vectors <- lapply(vector_list, function(vec) vec[common_names])
  
  # Combine the filtered vectors into a matrix
  combined_matrix <- do.call(cbind, filtered_vectors)
  
  # Take the median across the columns (i.e., across vectors)
  median_vector <- apply(combined_matrix, 1, median)
  
  # Assign the names back to the median vector
  names(median_vector) <- common_names
  
  return(median_vector)
}

# Apply the function to calculate the median across vectors
median_vector <- median_across_vectors(tip_maps)

# View the result
median_vector

# Function to calculate mode and proportion of occurrence across named vectors
mode_and_proportion_across_vectors <- function(vector_list) {
  # Identify common names across all vectors
  common_names <- Reduce(intersect, lapply(vector_list, names))
  
  # Extract only the values for the common names from each vector
  filtered_vectors <- lapply(vector_list, function(vec) vec[common_names])
  
  # Combine the filtered vectors into a matrix
  combined_matrix <- do.call(cbind, filtered_vectors)
  
  # Function to calculate mode and proportion for each row (name)
  calculate_mode_and_proportion <- function(row_values) {
    # Calculate the frequency of each unique value
    value_counts <- table(row_values)
    
    # Identify the mode (value with the highest frequency)
    mode_value <- as.numeric(names(value_counts)[which.max(value_counts)])
    
    # Calculate the proportion of times the mode appears
    mode_proportion <- max(value_counts) / length(row_values)
    
    return(c(mode = mode_value, proportion = mode_proportion))
  }
  
  # Apply the function to each row of the matrix
  mode_proportion_matrix <- t(apply(combined_matrix, 1, calculate_mode_and_proportion))
  
  # Convert the result to a dataframe and assign row names
  result_df <- as.data.frame(mode_proportion_matrix)
  rownames(result_df) <- common_names
  
  return(result_df)
}

# Apply the function to calculate mode and proportion across vectors
result_df <- mode_and_proportion_across_vectors(tip_maps)



# Function to calculate mode, proportions, and additional info for values
mode_proportion_details <- function(vector_list) {
  # Identify common names across all vectors
  common_names <- Reduce(intersect, lapply(vector_list, names))
  
  # Extract only the values for the common names from each vector
  filtered_vectors <- lapply(vector_list, function(vec) vec[common_names])
  
  # Combine the filtered vectors into a matrix
  combined_matrix <- do.call(cbind, filtered_vectors)
  
  # Function to calculate mode, proportions, and other counts
  calculate_mode_proportion <- function(row_values) {
    # Calculate the frequency of each unique value
    value_counts <- table(row_values)
    
    # Identify the mode (value with the highest frequency)
    mode_value <- as.numeric(names(value_counts)[which.max(value_counts)])
    
    # Calculate the proportion of the mode
    mode_proportion <- max(value_counts) / length(row_values)
    
    # If mode proportion is less than 0.8, capture other values and their proportions
    if (mode_proportion < 0.8) {
      other_proportions <- as.numeric(value_counts) / length(row_values)
      return(list(mode = mode_value, 
                  mode_proportion = mode_proportion, 
                  other_values = names(value_counts), 
                  other_proportions = other_proportions))
    } else {
      return(list(mode = mode_value, 
                  mode_proportion = mode_proportion))
    }
  }
  
  # Apply the function to each row of the matrix
  result_list <- apply(combined_matrix, 1, calculate_mode_proportion)
  
  return(result_list)
}

# Apply the function to calculate mode, proportion, and other counts
result_details <- mode_proportion_details(tip_maps)

# View the result
result_details

# Function to calculate mode, accounting for ties and returning the largest
mode_with_ties_handled <- function(vector_list) {
  # Identify common names across all vectors
  common_names <- Reduce(intersect, lapply(vector_list, names))
  
  # Extract only the values for the common names from each vector
  filtered_vectors <- lapply(vector_list, function(vec) vec[common_names])
  
  # Combine the filtered vectors into a matrix
  combined_matrix <- do.call(cbind, filtered_vectors)
  
  # Function to calculate mode, handle ties, and calculate proportions
  calculate_mode_and_proportion <- function(row_values) {
    # Calculate the frequency of each unique value
    value_counts <- table(row_values)
    
    # Find the maximum frequency (mode frequency)
    max_count <- max(value_counts)
    
    # Get all values with the maximum frequency
    tied_values <- as.numeric(names(value_counts)[value_counts == max_count])
    
    # If there are ties, select the largest value
    mode_value <- max(tied_values)
    
    # Calculate the proportion of the mode
    mode_proportion <- max_count / length(row_values)
    
    # Return mode and proportion
    return(c(mode = mode_value, proportion = mode_proportion))
  }
  
  # Apply the function to each row of the matrix
  mode_proportion_matrix <- t(apply(combined_matrix, 1, calculate_mode_and_proportion))
  
  # Convert the result to a dataframe and assign row names
  result_df <- as.data.frame(mode_proportion_matrix)
  rownames(result_df) <- common_names
  
  return(result_df)
}

# Apply the function to calculate mode and proportion with tie-handling
result_df <- mode_with_ties_handled(tip_maps)


selected_vectors <- result_details[sapply(result_details, function(x) length(x) == 4)]

library(treeplyr)
tree<-readRDS("ext2_td.RDS")

mmap<-make.simmap(tree$phy, median_vector)

mode_vector<-result_df[,1]
names(mode_vector)<-rownames(result_df)

library(viridis)


values <- c(mean(result_df$proportion[which(result_df[,1]==1)]), mean(result_df$proportion[which(result_df[,1]==2)]), mean(result_df$proportion[which(result_df[,1]==3)]), mean(result_df$proportion[which(result_df[,1]==4)]))

# Generate a magma color palette
# The 'magma' palette with a range of 256 colors (default)
palette <- viridis(256, option = "plasma", direction = -1)


colors <- palette[scales::rescale(values, to = c(1, 256), from = c(0, 1))]

# View the colors
names(colors)<-1:4
colors[3]<-"black"

smap<-phytools::make.simmap(tree$phy, mode_vector)
plot(smap, col=colors, ftype="off")





