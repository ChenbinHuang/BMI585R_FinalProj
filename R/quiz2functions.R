#' Compute Sensitivity based on prediction and reference truth value
#' 
#' Compute Sensitivity based on binary prediction and reference truth value.
#' Data should be boolean, or 1 for positive and 0 for negative
#' @export
#'
#' @param pred prediction result
#' @param truth reference(truth) result
#' @return sensitivity
#' @examples
#' sensitivity(c(1,0,1,1,0),c(1,0,1,1,1))
#' sensitivity(c(T,F,T,T,F),c(F,T,T,F,F))
sensitivity = function(pred,truth){
  pred = as.logical(pred)
  truth = as.logical(truth)
  ntp = sum(pred[truth])
  nap = sum(truth)
  return(ntp/nap)
}

#' Compute specificity based on prediction and reference truth value
#' 
#' Compute specificity based on prediction and reference truth value
#' Data should be boolean, or 1 for positive and 0 for negative
#' @export
#'
#' @param pred prediction result
#' @param truth reference(truth) result
#' @return specificity
#' @examples
#' specificity(c(1,0,1,1,0),c(1,0,1,1,1))
#' specificity(c(T,F,T,T,F),c(F,T,T,F,F))
specificity = function(pred,truth){
  pred = as.logical(pred)
  truth = as.logical(truth)
  ntn = sum(!pred[!truth])
  nan = sum(!truth)
  return(ntn/nan)
}

#' Compute accuracy based on prediction and reference truth value
#' 
#' Compute accuracy based on prediction and reference truth value
#' Data should be boolean, or 1 for positive and 0 for negative
#' @export
#'
#' @param pred prediction result
#' @param truth reference(truth) result
#' @return accuracy
#' @examples
#' accuracy(c(1,0,1,1,0),c(1,0,1,1,1))
#' accuracy(c(T,F,T,T,F),c(F,T,T,F,F))
accuracy = function(pred,truth){
  pred = as.logical(pred)
  truth = as.logical(truth)
  # ntp = sum(pred[truth])
  # ntn = sum(!pred[!truth])
  nt = sum(pred == truth)
  return(nt/length(truth))
}

#' Compute precision based on prediction and reference truth value
#' 
#' Compute precision based on prediction and reference truth value
#' Data should be boolean, or 1 for positive and 0 for negative
#' @export
#'
#' @param pred prediction result
#' @param truth reference(truth) result
#' @return precision
#' @examples
#' ppv(c(1,0,1,1,0),c(1,0,1,1,1))
#' ppv(c(T,F,T,T,F),c(F,T,T,F,F))
ppv = function(pred,truth){
  pred = as.logical(pred)
  truth = as.logical(truth)
  ntp = sum(pred[truth])
  return(ntp/sum(pred))
}

#' Compute f1 score based on prediction and reference truth value
#' 
#' Compute f1 score based on prediction and reference truth value
#' Data should be boolean, or 1 for positive and 0 for negative
#' @export
#'
#' @param pred prediction result
#' @param truth reference(truth) result
#' @return f1 score
#' @examples
#' f1(c(1,0,1,1,0),c(1,0,1,1,1))
#' f1(c(T,F,T,T,F),c(F,T,T,F,F))
f1 = function(pred,truth){
  recall = sensitivity(pred,truth)
  precision = ppv(pred,truth)
  2*(precision * recall) / (precision + recall)
}
