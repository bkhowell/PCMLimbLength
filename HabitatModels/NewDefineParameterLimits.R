# lower limits for models A and B
library(PCMBase)
PCMParamLowerLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$Sigma_x)) {
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- -.0
    }
  } else {
    if(!is.Diagonal(o$Sigma_x)) {
      for(r in seq_len(R)) {
        o$Sigma_x[1, 2, r] <- -.0
      }
    }
  }
  o
}

# upper limits for models A and B
PCMParamUpperLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- 1.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- 1.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- 1.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- 1.0
      }
    }
  }
  o
}

# lower limits for models C, ..., F.
PCMParamLowerLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$Theta)) {
    o$Theta[1] <- 0.0
    o$Theta[2] <- -1.2
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- -0.5
      o$Theta[2, r] <- -0.5
    }
  }
  if(is.Global(o$Sigma_x)) {
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- -.0
    }
  } else {
    if(!is.Diagonal(o$Sigma_x)) {
      for(r in seq_len(R)) {
        o$Sigma_x[1, 2, r] <- -.0
      }
    }
  }
  if(is.Global(o$H)) {
    o$H[1,1]<-o$H[2,2]<- 0.001
    if(!is.Diagonal(o$H)) {
      o$H[1, 2] <- -10
    }
  } else {
        for(r in seq_len(R)) {
        o$H[1, 1, r] <- o$H[2, 2, r] <- .001
        if(!is.Diagonal(o$H)) {
          o$H[1, 2, r] <- -10
          
        }
    }
  }
  o
}

# upper limits for models C, ..., F.
PCMParamUpperLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))
  
  if(is.Global(o$Theta)) {
    o$Theta[1] <- 7.8
    o$Theta[2] <- 6.8
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- 7.8
      o$Theta[2, r] <- 6.8
    }
  }
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- 1.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- 1.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- 1.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- 1.0
      }
    }
  }
  if(is.Global(o$H)) {
    o$H[1,1]<-o$H[2,2]<- 10
    if(!is.Diagonal(o$H)) {
      o$H[1, 2] <- 10
    }
  } else {
    for(r in seq_len(R)) {
      o$H[1, 1, r] <- o$H[2, 2, r] <- 10
      if(!is.Diagonal(o$H)) {
        o$H[1, 2, r] <- 10
        
      }
    }
  }
  o
}