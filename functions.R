## https://github.com/timothyslau/ICC.merMod
# Gaussian ICC lme4 (numerator can also be a list of values)
# function based off of http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1540459/
ICC.GAU <- function(model, numerator){
  require(lme4)
  mout <- data.frame(VarCorr(model)) # random intercept model variances
  sigma_a2 <- sum(mout[mout$grp %in% numerator, "vcov"]) # random effect(s) in numerator
  sigma_2 <- sum(mout["vcov"]) # sum of random effects variance in denominator
  icc <- sigma_a2 / sigma_2
  return(icc)
}

## Function to calculate chi-square for a one-way random effects ANOVA in lme4
## (e.g. the unconditional model)
rANOVA_chisq = function(x) {
  require(dplyr)
  require(magrittr)
  
  # get name of grouping and response variables
  levname = names(ranef(x))
  resp_name = as.character(x@call[[2]][[2]])
  
  ## Reconstruct dataset from model object
  dat = x@frame[,c(levname, resp_name)]
  names(dat) = c("grouping", "response")
  
  ## Grand mean
  fixed = fixef(x)
  
  # get OLS level-2 residuals
  ols_lev2 = dat %>%
    group_by(grouping) %>%
    summarize(resid = mean(response) - fixed) %>%
    extract2("resid")
  
  # Get n per group
  nj = as.integer(table(dat$grouping))
  
  # Get numbers of groups
  J = length(nj)
  
  # Get residual variance
  sig2 = summary(x)$sigma^2
  
  ## Test statistic for variance components
  chisq = sum(nj * ols_lev2^2)/sig2
  
  pval_TS = pchisq(chisq, df = J - 1, lower.tail = FALSE)
  
  return(c(chisq = chisq, df = J - 1,  p = pval_TS))
}
