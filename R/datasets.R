#' A BPP A00 MCMC sample for a mouse lemur phylogeny
#'
#' @description This dataset contains the results from the BPP A00 analysis of mouse
#' lemur evolution in Madagascar from Yoder et al. (2016).
#'
#' @format \code{microcebus} is a list with elements \code{mcmc}, a dataframe
#'   with 20,000 rows and 12 columns, and \code{tree}, an object of class
#'   \code{phylo} from the \code{ape} package.
#'
#'   \code{mcmc} is a posterior sample from a BPP A00 MCMC analysis containing
#'   the relative divergence times (tau's) and nucleotide diversities (theta's)
#'   for the six species mouse lemur (\emph{Microcebus} spp) phylogeny.
#'
#'   \code{tree} contains the phylogeny with node ages given as the posterior
#'   means of the tau's in \code{mcmc}.
#'
#' @source A. D. Yoder, C. R. Campbell, M. B. Blanco, M. dos Reis, J. U.
#'   Ganzhorn, S. M. Goodman, K. E. Hunnicutt, P. A. Larsen, P. M. Kappeler, R.
#'   M. Rasoloarison, J. M. Ralison, D. L. Swofford, and D. W. Weisrock. (2016)
#'   \emph{Geogenetic patterns in mouse lemurs (genus Microcebus) reveal the
#'   ghosts of Madagascar's forests past.} Proc. Nat. Acad. Sci. USA., 113:
#'   8049--8056.
#'
#' @seealso \code{\link{hominids}}
#'
"microcebus"

#' A BPP A00 MCMC sample for an hominid phylogeny
#'
#' @description This dataset contains the results from the BPP A00 analysis of
#'   hominid evolution from Angelis and dos Reis (2015).
#'
#' @format \code{hominids} is a list with elements \code{mcmc}, a dataframe with
#'   20,000 rows and 8 columns, and \code{tree}, an object of class
#'   \code{phylo} from the \code{ape} package.
#'
#'   \code{mcmc} is a posterior sample from a BPP A00 MCMC analysis containing
#'   the relative divergence times (tau's) and nucleotide diversities (theta's)
#'   for the four species ape (hominid) phylogeny.
#'
#'   \code{tree} contains the phylogeny with node ages given as the posterior
#'   means of the tau's in \code{mcmc}.
#'
#' @source K. Angelis and M. dos Reis (2015) \emph{The impact of ancestral
#'   population size and incomplete lineage sorting on Bayesian estimation of
#'   species divergence times.} Curr. Zool., 61: 874--885.
#'
#' @seealso \code{\link{microcebus}}
"hominids"
