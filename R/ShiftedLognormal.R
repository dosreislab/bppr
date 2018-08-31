# #############################################################################
# THE SHIFTED LOG-NORML DISTRIBUTION
# #############################################################################

#' The Shifted Log-normal Distribution
#'
#' Density, distribution and quantile functions, and random number generation
#' for the shifted log-normal distribution.
#'
#' @param x,q vector of quantiles.
#' @param shift vector of shifts.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'   to be the number required.
#' @param meanlog,sdlog mean and standard deviation of the distribution on the
#'   log scale with default values of 0 and 1 respectively.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X ≤
#'   x]}, otherwise, \eqn{P[X > x]}.
#'
#' @details Let \eqn{Y} have a log-normal distribution with parameters \eqn{\mu}
#'   (\code{meanlog}) and \eqn{\sigma} (\code{sdlog}). Then \eqn{X = Y + s, (s >
#'   0)} has a shifted log-normal distribution with shift \eqn{s}
#'   (\code{shift}), mean \eqn{E(X) = exp(μ + 1/2 σ^2) + s} and variance
#'   \eqn{Var(X) = exp(2*μ + σ^2)*(exp(σ^2) - 1)}.
#'
#'   Note \code{[dpqr]slnorm} are wrappers for the corresponding
#'   \code{[dpqr]lnorm} functions.
#'
#' @return \code{dslnorm} gives the density, \code{pslnorm} gives the distribution
#'   function, \code{qslnorm} gives the quantile function, and \code{rslnorm}
#'   generates random deviates.
#'
#'   The length of the result is determined by \code{n} for \code{rlnorm}, and
#'   is the maximum of the lengths of the numerical arguments for the other
#'   functions.
#'
#'   The numerical arguments other than \code{n} are recycled to the length of
#'   the result. Only the first elements of the logical arguments are used.
#'
#' @seealso \code{\link{Lognormal}}
#'
#' @examples
#' curve(dslnorm(x, shift=6.5), from=0, to=15, n=1e3)
#' rr <- rslnorm(1e3, shift=6.5)
#' lines(density(rr, adj=.1), lty=2)
#'
#' all.equal (qslnorm(c(.025, .9), shift=6.5) - 6.5, qlnorm(c(.025, .9)))
#' all.equal (pslnorm(10, shift=6.5), plnorm(10 - 6.5))
#'
#' @name ShiftedLognormal
NULL

#' @rdname ShiftedLognormal
#' @export
dslnorm <- function (x, shift, meanlog = 0, sdlog = 1, log = FALSE) {
  return( dlnorm(x - shift, meanlog, sdlog, log) )
}

#' @rdname ShiftedLognormal
#' @export
pslnorm <- function (q, shift, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE) {
  return( plnorm(q - shift, meanlog, sdlog, lower.tail, log.p) )
}

#' @rdname ShiftedLognormal
#' @export
qslnorm <- function (p, shift, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE) {
  return( qlnorm(p, meanlog, sdlog, lower.tail, log.p) + shift)
}

#' @rdname ShiftedLognormal
#' @export
rslnorm <- function (n, shift, meanlog = 0, sdlog = 1) {
  return( rlnorm(n, meanlog, sdlog) + shift )
}
