#===============================================================================
# allelicimbalance.R
#===============================================================================

# Allelic imbalance utilities


# Imports ======================================================================

#' @import VGAM




# Functions ====================================================================

#' @title p-value for a beta-binomial variable
#'
#' @description Hypothesis test against a beta-binomial null.
#' @param q vector of quantiles
#' @param size vector of sample size parameters
#' @param shape1,shape2 the two (positive) shape parameters of the standard
#'   beta distribution. See the documentation for \code{Betabinom} in the
#'   \code{VGAM} package.
#' @return numeric. The p-value.
#' @export
beta_binom_pval <- function(q, size, shape1, shape2) {
  lower_tail_area <- pbetabinom.ab(q, size, shape1, shape2)
  upper_tail_area <- 1 - lower_tail_area + dbetabinom.ab(
    q,
    size,
    shape1,
    shape2
  )
  min(1, 2 * pmin(lower_tail_area, upper_tail_area))
}

#' @title log posterior allelic fold change
#'
#' @description Bayesian effect size for imbalance data
#' @param count1,count2 allele counts
#' @param shape1,shape2 shape parameters
#' @param level confidence level for credible interval
#' @return \describe{
#'   \item{lpafc}{the effect size.}
#'   \item{lower}{lower bound for the credible interval.}
#'   \item{upper}{upper bound for the credible interval.}
#' }
#' @export
log_posterior_allelic_fold_change <- function(
  count1,
  count2,
  shape1,
  shape2,
  level = 0.99
) {
  posterior_mean <- (shape1 + count1) / (shape1 + count1 + shape2 + count2)
  posterior_mean_quantile <- pbeta(
    posterior_mean,
    shape1 + count1,
    shape2 + count2
  )
  lower_prob <- posterior_mean_quantile - level / 2
  upper_prob <- posterior_mean_quantile + level / 2
  if (lower_prob > 0 && upper_prob < 1) {
    lower <- qbeta(lower_prob, shape1 + count1, shape2 + count2)
    upper <- qbeta(upper_prob, shape1 + count1, shape2 + count2)
  } else if (lower_prob > 0) {
    lower <- qbeta(1 - level, shape1 + count1, shape2 + count2)
    upper <- 1
  } else if (upper_prob < 1) {
    lower <- 0
    upper <- qbeta(level, shape1 + count1, shape2 + count2)
  }
  shift_term <- log2(shape1) - log2(shape2)
  list(
    lpafc = log2(1 - posterior_mean) - log2(posterior_mean) + shift_term,
    lower = log2(1 - lower) - log2(lower) + shift_term,
    upper = log2(1 - upper) - log2(upper) + shift_term
  )
}
