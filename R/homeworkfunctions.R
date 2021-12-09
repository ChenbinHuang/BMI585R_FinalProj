#' Simulate Normal Data with Box Muller Tranformation
#' @export
#'
#' @param n number of data to simulate
#' @return simulated result, 1D vector with length n
#' @examples
#' sim = boxMuller(500)
boxMuller = function(n = 1){
  u1 = runif(n)
  u2 = runif(n)
  sqrt(-2*log(u1))*cos(2*pi*u2)
}

#' two sided t test
#' @export
#'
#' @param t t-statistic
#' @param df degree of freedom
#' @return p-value as decimal
#' @examples
#' twoSidedT(2,5)
twoSidedT = function(t,df){
  t = abs(t)
  1-(pt(t,df)-pt(-t,df))
} 

#' two sided Z test
#' @export
#'
#' @param z z-statistic
#' @return p-value as decimal
#' @examples
#' twoSidedZ(2)
twoSidedZ = function(z){
  z = abs(z)
  1-(pnorm(z)-pnorm(-z))
} 

#' effect size from data and group
#' @export
#'
#' @param x numeric data
#' @param g group allocation of data
#' @return effect size(cohen's d)
#' @examples
#' effectSize(runif(20), sample(0:1,20, replace=T))
effectSize = function(x,g){
  sd = sd(x,na.rm=T)
  dplyr::tibble(x,g) |> dplyr::group_by(g) |> 
    dplyr::summarize(x = mean(x,na.rm=T)) |> dplyr::pull(x) |> 
    diff() |> abs() / sd
}

#' welch t test
#' @export
#'
#' @param x dataset1
#' @param y dataset2 to compare
#' @return class Welch.t.test. includes p.value, t-statistic, degree of freedom
#' mean of dataset1, mean of dataset2
#' @examples
#' welchT(c(1:10), c(7:20))
welchT = function(x,y){
  # exclude empty values
  x = x[!is.na(x)]
  y = y[!is.na(y)]
  
  # calculate group summary statistics
  n.x = length(x)
  n.y = length(y)
  mn.x = mean(x)
  mn.y = mean(y)
  sd.x = sd(x)
  sd.y = sd(y)
  
  # calculate the test statistic
  t = abs(mn.x - mn.y)/sqrt(sd.x^2/n.x + sd.y^2/n.y)
  
  # calculate the degrees of freedom
  r = (sd.x^2*n.x^-1 + sd.y^2*n.y^-1)^2 / (sd.x^4*n.x^-2*(n.x-1)^-1 + sd.y^4*n.y^-2*(n.y-1)^-1)
  
  p.value = 1-(pt(t,r)-pt(-t,r))
  
  # return the results as a tibble
  out = dplyr::tibble(p.value,t,df = r, mn.x, mn.y)
  
  # create a class statement for a custom print method
  class(out) <- "welch.t.test"
  out
}

#' minimum sample size n with power=0.8
#' @export
#'
#' @param d effect size (cohen's d)
#' @return minimum sample size n
#' @examples
#' minimumN(0.8)
minimumN = function(d) {
  ptt = stats::power.t.test(power=0.8, delta =d)
  ceiling(ptt$n)
}

#' chi-square test for count data in a tibble or similar object
#' @export
#'
#' @param tib count data in a tibble or similar object
#' @return chisquare test result include chisq statistic, p-value, degree of freedom
#' @examples
#' data = data.frame(x1=c(15,34), x2 =c(65,46))
#' chiSquareCounts(data)
chiSquareCounts = function(tib) {
  
  ratio = colSums(tib) / sum(colSums(tib))
  data_row_exp = rowSums(tib)%*%t(ratio)
  
  df = (nrow(tib) -1 )*(ncol(tib) -1 )
  if (df == 1){
    cat("Corrected with Yates's\n")
    chisq_v = sum((abs(tib - data_row_exp)-0.5)^2 / data_row_exp)
  } else {
    cat("Uncorrected\n")
    chisq_v = sum((tib - data_row_exp)^2 / data_row_exp)
  }
  
  p_val=1-pchisq(chisq_v,df=df)
  
  list("chisq"=chisq_v, "pValue"=p_val, "df"=df)
}

#' estimate of post-hoc power based on 1000 simultions
#' @export
#'
#' @param d effect size (cohen's d)
#' @param n1 size of group 1
#' @param n2 size of group 2
#' @return estimated power
#' @examples
#' postHocPower(1,12,13)
postHocPower = function(d, n1, n2){
  sim = function(d, n1, n2){
    essd = 1
    est_M = essd*d
    gp1 <- rnorm(n1,mean=-est_M/2,sd=essd)
    gp2 <- rnorm(n2,mean=est_M/2,sd=essd)
    tt = stats::t.test(gp1,gp2)
    tt$p.value
  }
  data = replicate(1000, sim(d, n1, n2))
  power = sum(data > 0.05)
  1- power/1000
}

#' Bonferroni-Holm adjusted p-value
#' @export
#'
#' @param p p-values
#' @return boolean vector of each p-value is significant or not
#' @examples
#' test.data = seq(0.0025,0.0250,0.0025)
#' test.shuffle = sample(test.data)
#' bhAdjust(test.data)
#' bhAdjust(test.shuffle)
bhAdjust = function(p) {
  alpha = 0.05
  adj.order = order(p)
  alpha_levels = alpha / c(length(p):1)
  res = p < alpha_levels[order(adj.order)]
  resend = FALSE
  for (i in adj.order) {
    if (resend){
      res[i] = FALSE
    } else {
      if (res[i]) {
        next
      } else {
        resend = TRUE
      }
    }
  }
  res
}

#' FDR adjusted p-value
#' @export
#'
#' @param p p-values
#' @return boolean vector of each p-value is significant or not
#' @examples
#' test.data = seq(0.0025,0.0250,0.0025)
#' test.shuffle = sample(test.data)
#' fdrAdjust(test.data)
#' fdrAdjust(test.shuffle)
fdrAdjust = function(p){
  alpha = 0.05
  adj.order = order(p, decreasing=T)
  res = p < 0.05*order(order(p)) / length(p)
  findtrue = FALSE
  for (i in adj.order) {
    if (findtrue){
      res[i] = TRUE
    } else {
      if (res[i]==FALSE) {
        next
      } else {
        findtrue = TRUE
      }
    }
  }
  res
}

#' Compute R square
#' @export
#'
#' @param pred prediction result
#' @param truth reference(truth) result
#' @return r square score
#' @examples
#' r2(1:5,2:6)
r2= function(pred,truth){
  residuals = truth-pred
  rss = sum(residuals^2)
  tss = sum((truth-mean(truth))^2)
  1-(rss/tss)
}

#' Compute Adjusted R square
#' @export
#'
#' @param pred prediction result
#' @param truth reference(truth) result
#' @param n the number of predictors
#' @return adjusted r square score
#' @examples
#' adjR2(1:5,2:6,1)
adjR2= function(pred,truth,n){
  residuals = truth-pred
  rss = sum(residuals^2)
  tss = sum((truth-mean(truth))^2)
  n.data = length(truth)
  1-(rss/(n.data-n-1))/(tss/(n.data-1))
}
