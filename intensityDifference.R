#' Intensity difference test using pre-calculated differences and average
#' 
#' This is a modification of Simon's intensity difference code which, instead of two sets of values, 
#' takes pre-calculated raw differences and an average. For example, where the differences to assess
#' are from ratios of expression values (stability, translation efficiency), where variability will 
#' scale with the variability in the expression values from which they were calculated. So the 
#' averages provided would be based on the expression values (one of them, or combined in some way), 
#' and the differences on the ratios.
#' 
#' Also now added option to provide values or averages/differences as a list of numeric vectors. 
#' This would allow you to calculate z scores across multiple comparisons, with joint modelling - 
#' so theoretically this might make them more comparable. Although for current manuscript not needed.
#' 
#' This function calculates an intensity corrected z-score
#' statistic for a matched pair of numeric vectors. For each
#' point a local sub-sampling of the data is performed where
#' a set of points with the most similar mean values across
#' the two datasets are collected.  The distribution of
#' differences for these selected points is then modelled with
#' a normal distribution and a z-score value for the point
#' being tested is then calculated, along with the p.value
#' corresponding to that z-score.  The absolute difference
#' and FDR corrected p-values are also reported.
#'
#' @param values.1 A numeric vector of values
#' @param values.2 A numeric vector of values
#' @param averages A numeric vector of averages based eg on expression level
#' @param differences A numeric vector of raw differences
#' @param window.proportion What proportion of the data will be used to create each local model - Default 0.01
#' @export
#' @examples
#' intensity.difference()

intensity.difference <- function (
  values.1 = NULL,
  values.2 = NULL,
  averages = NULL,
  differences = NULL,
  window.proportion=0.01
) {
  
  return_as_list <- FALSE
  if(length(values.1) > 0 & length(values.2) > 0) {
    if(is.list(values.1)) {
      if(!length(values.1) == length(values.2)) {
        stop("To provide as a list, values.1 and values.2 must both be lists of equal length")
      }
      list_index <- integer()
      for(l in 1:length(values.1)) {
        print(paste("processing list",l))
        if(!is.numeric(values.1[[l]])) {
          stop("The data in values.1 is not numeric")
        }
        if(!is.numeric(values.2[[l]])) {
          stop("The data in values.2 is not numeric")
        }
        if(!length(values.1[[l]]) == length(values.2[[l]])) {
          stop("Length of values.1 must be equal to length of values.2")
        }
        list_index <- c(list_index, rep(l,length(values.1[[l]])))
      }
      values.1 <- unlist(values.1)
      values.2 <- unlist(values.2)
      return_as_list <- TRUE
    }
    
    if (!is.numeric(values.1)) {
      stop("The data in values.1 was not numeric")
    }
    
    if (!is.numeric(values.2)) {
      stop("The data in values.2 was not numeric")
    }
    
    
    if (length(values.1) != length(values.2)) {
      stop("The two vectors passed to intensity.difference must be the same length")
    }
    
    return.frame <- data.frame(values.1,values.2,index=1:length(values.1))
    
    (values.1+values.2)/2 -> return.frame$average.intensity
    
    values.1-values.2 -> return.frame$difference
  } else {
    if(!(length(averages) > 0 & length(differences) > 0)){
      stop("Corresponding vectors either for values.1 and values.2, or for averages and differences, must be provided")
    }
    if(is.list(averages)) {
      if(!length(averages) == length(differences)) {
        stop("To provide as a list, averages and differences must both be lists of equal length")
      }
      list_index <- integer()
      for(l in 1:length(averages)) {
        print(paste("processing list",l))
        if(!is.numeric(averages[[l]])) {
          stop("The data in averages is not numeric")
        }
        if(!is.numeric(differences[[l]])) {
          stop("The data in differences is not numeric")
        }
        if(!length(averages[[l]]) == length(differences[[l]])) {
          stop("Length of averages must be equal to length of differences")
        }
        list_index <- c(list_index, rep(l,length(averages[[l]])))
      }
      averages <- unlist(averages)
      differences <- unlist(differences)
      return_as_list <- TRUE
    }
    if (!is.numeric(averages)) {
      stop("The data in averages was not numeric")
    }
    
    if (!is.numeric(differences)) {
      stop("The data in differences was not numeric")
    }
    
    
    if (length(values.1) != length(values.2)) {
      stop("The two vectors passed to intensity.difference must be the same length")
    }
    
    return.frame <- data.frame(averages, differences, index = 1:length(averages))
    colnames(return.frame)[1:2] <- c("average.intensity", "difference")
  }
  if(exists("list_index") == TRUE) { return.frame$list_index <- list_index }
  
  return.frame[order(return.frame$average.intensity),] -> return.frame
  
  return.frame$p.value <- 1
  return.frame$local.sd <- 0
  return.frame$z.score <- 0
  
  window.size <- as.integer(nrow(return.frame)*window.proportion)
  half.window.size <- as.integer(window.size/2)
  
  if (half.window.size < 1) {
    stop(paste("Sample size is too small when using window.proportion of",window.proportion))
  }
  
  
  sapply(1:nrow(return.frame), function(x) {
    start <- x-half.window.size
    if (start < 0) start <- 0
    end <- start+window.size
    if (end > nrow(return.frame)) {
      end <- nrow(return.frame)
      start <- end-window.size
    }
    
    local.diffs <- return.frame$difference[start:end]
    
    # We assume a mean of 0 and calculate the sd
    sqrt(mean(local.diffs*local.diffs)) -> local.sd
    
    # Now we work out the p.value for the value we're actually
    # looking at in the context of this distibution
    pnorm(return.frame$difference[x],mean=0,sd=local.sd) -> local.p
    
    if (local.p >0.5) local.p <- (1 - local.p)
    
    return.frame$z.score[x] <<- return.frame$difference[x] / local.sd
    
    return.frame$local.sd[x] <<- local.sd
    return.frame$p.value[x] <<- local.p
    
  }
  )
  
  return.frame$fdr.value <- p.adjust(return.frame$p.value,method="fdr")
  
  return.frame[order(return.frame$index),] -> return.frame
  
  if(return_as_list == TRUE) {
    return_frame_list <- list()
    for(l in unique(list_index)){
      return_frame_list[[l]] <- return.frame[return.frame$list_index == l,]
    }
    return.frame <- return_frame_list
  }
  
  return(return.frame)
  
}
