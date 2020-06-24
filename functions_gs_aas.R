




####### CULTURAL CONTINUITIES AND DISCONTINUITIES AT THE NEOLITHIC TRANSITION IN EASTERN
####### IBERIA: AN ANLYSIS OF THE MORPHOMETRY OF THE GEOMETRIC MICROLITHS 

#### Archaeological and Anthropological Sciences

#### Article authors: A. Cortell-Nicolau, O. García-Puchol, S. Shennan

#### Code author: A. Cortell-Nicolau

#### London, 2018-2019 

#### Reviewed: València, 2020

#### License: Permission is granted to use and adapt this code. Please do not distribute
#### without permission of the author

#### Warranty: No warranty is expressed or implied. Use at your own risk

#### About: This script provides five functions used during the Geomeasure analysis script

#### FUNCTIONS USED IN GEOMEASURE ANALYSIS

## FUNCTION 1. Selects only full L-measures
select_full_lengths <- function(x){ ## Specific for the data as it is ordered in this article
  x <- x[,c(25:134)] 
}


## FUNCTION 2. Selects only half L-measures
select_half_lengths <- function(x){ ## Specific for the data as it is ordered in this article
  x <- x[,c(137:236)] 
}

## FUNCTION 3. Computes Huber M-estimate. If MAD==0, computes mean
## Modified from huber() function, in MASS package (Ripley, 2020)
huber_modified <- function (y, k = 1.5, tol = 1e-06) 
{
  y <- y[!is.na(y)]
  n <- length(y)
  mu <- median(y)
  s <- mad(y)
  if (s == 0) {
    mu = mean(y)
  } else
    
    repeat {
      yy <- pmin(pmax(mu - k * s, y), mu + k * s)
      mu1 <- sum(yy)/n
      if (abs(mu - mu1) < tol * s) 
        break
      mu <- mu1
    }
  list(mu = mu, s = s)
}

## FUNCTION 4.1 Substitutes na with means (for each column in the list)
na_to_mean <- function(x){
  for (i in 1:length(x)){
    for (j in 1:ncol(x[[i]])){
      mean_wo_na <- x[[i]][,j][!is.na(x[[i]][,j])] 
      subs_hub <- huber_modified(mean_wo_na,2)
      subs_value <- subs_hub$mu
      x[[i]][,j][is.na(x[[i]][,j])] <- subs_value
      x[[i]][,j] <- round(x[[i]][,j],2)
    }  
  }
  return(x)
}

## FUNCTION 4.2 Substitutes na with means (for vectors)
na_to_mean_vec <- function(x){
  vect <- x[!is.na(x)]
  subs_hub <- huber_modified(vect,2)
  subs_value <- round(subs_hub$mu,2)
  x[is.na(x)] <- subs_value
  return(x)
}

## FUNCTION 5. Only variables with five or more measures taken are accepted. 0 values means no presence
## of the geometric
five_or_more <- function(x){
  for (j in 1:length(x)){
    for (i in ncol(x[[j]]):1){
      dimen <- length(x[[j]][,i][which(x[[j]][,i]>0)])  
      if (dimen < 5){
        x[[j]][,i] <- NULL
      }
    }
  }
  return(x)
}

## FUNCTION 6. Intervariance function
# x = data (must be a list)
# n = sample size
# m = average of group
# y = overall average
# g = number of groups

## Using Huber M-estimator for robust means

intervariance <- function(x){
  
  ## Compute the averages (huber) for each group/phase 
  nvar <- length(x[[1]]) ## All objects in list have the same length
  var_hub <- c()
  for (i in 1:nvar){
    measures <- c()
    for (j in 1:length(x)){
      measures <- append(measures,x[[j]][,i])
      measures <- measures[measures!=0] ## Remove 0s (0 == geometric not present on line)
    }
    sing_hub <- huber(measures,2)
    sing_hub <- sing_hub$mu
    var_hub <- append(var_hub,sing_hub)
  }
  names(var_hub) <- colnames(x[[1]])
  
  ## Compute intervariance
  results <- c()
  for (j in 1:length(var_hub)){ ## Will work for each variable
    numerator <- 0
    for (i in 1:length(x)){
      n <- nrow(x[[i]]) ## Extracts sample size (per group)
      for_means <- unlist(x[[i]][j])
      for_means <- for_means[for_means!=0] ## Remove 0s (means unexistent)
      m <- huber(for_means,2) ## Computes the mean (per group) winsorizes at a = 2
      m <- m$mu
      var <- n*(m-var_hub[j])^2 ## Computes elements for numerator individually
      numerator <- numerator+var ## Sums the elements of the numerator
    }
    loc_results <- numerator/(length(x)-1) ## Intervariance product  
    results <- append(results,loc_results)
  }
  return(results)
}




