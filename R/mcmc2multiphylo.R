#' Convert an MCMC sample from BPP or MCMCTree to a list of trees
#'
#' @param tree an object of class phylo
#' @param mcmc data frame with an MCMC sample from MCMCTree or a BPP A00
#'   analysis
#' @param time.name character vector of length one
#' @param thin numeric, the fraction of MCMC samples to keep
#'
#' @details \code{tree} must be rooted and strictly bifurcating, and it must
#'   match the tree used by BPP or MCMCTree to obtain the MCMC sample. The
#'   function uses the node ages in \code{mcmc} to calculate branch lengths and
#'   generate a list of trees (with the same topology as \code{tree}), one tree
#'   per (thinned) MCMC sample. The tips of the phylogeny are assumed to have
#'   age zero.
#'
#' @return An object of class multiPhylo (i.e., a list of trees).
#'
#' @examples
#' data(microcebus)
#' # convert a BPP A00 MCMC sample of Microcebus spp. to a list of trees
#' mtts <- mcmc2multiphylo(microcebus$tree, microcebus$mcmc, "tau_", thin=0.01)
#' length(mtts)
#'
#' data(hominids)
#' # Calibrate the hominid phylogeny with a uniform fossil calibration of
#' # between 6.5 to 10 Ma for the human-chimp divergence.
#' calmsc <- msc2time.t(mcmc=hominids$mcmc, node="7humanchimp", calf=runif,
#'                     min=6.5, max=10)
#' # convert the time-calibrated MCMC sample to a list of trees
#' htts <- mcmc2multiphylo(hominids$tree, calmsc, "t_", thin=0.01)
#' htts[[1]]
#'
#' \dontrun{
#' # If you have the ape package installed, you can output the trees in Newick
#' ape::write.tree(htts[1:5])
#'
#' # The trees are suitable for plotting with the phangorn package
#' # Relative node ages (tau's):
#' mcon <- microcebus$tree$tip.label
#' phangorn::densiTree(mtts, col="blue", alpha=0.04, cons=mcon, label.offset=.01)
#'
#' # Absolute node ages (in millions of years):
#' hcon <- hominids$tree$tip.label
#' phangorn::densiTree(htts, col="blue", alpha=0.04, cons=hcon, label.offset=.01)
#' }
#'
#' @author Mario dos Reis
#'
#' @export
mcmc2multiphylo <- function(tree, mcmc, time.name, thin) {
  tts <- list()
  N <- dim(mcmc)[1]
  ti <- grep(time.name, names(mcmc))
  ii <- floor(seq(from=1, to=N, length.out = N * thin))
  n <- length(ii)

  for (i in 1:n) {
    tts[[i]] <- .mcmc2multiphylo.treef(tree, unlist(mcmc[ii[i],ti]))
  }
  class(tts) <- "multiPhylo"

  return (tts)
}

# This helper function does pre-order tree traversal to calculate branch lengths
# from node ages.
.mcmc2multiphylo.treef <- function(tree, ages) { # TODO: Tip ages
  # tree is assummed rooted and perfectly bifurcating
  ns <- tree$Nnode + 1
  root <- ns + 1
  tip.ages <- rep(0, ns) # assumes tips are extant
  all.ages <- c(tip.ages, ages)
  blens <- numeric(2*ns - 2)

  porder <- function(node) { # parent is NULL if root
    #print(paste("i'm", node))
    desc <- which(tree$edge[,1] == node) # my descendant branches

    if (node != root) {
      parent <- tree$edge[which(tree$edge[,2] == node)]
      mybranch <- which(tree$edge[,2] == node)
      blens[mybranch] <<- all.ages[parent] - all.ages[node]
      #print(blens)
    }

    # If I'm a tip:
    if (!length(desc)) { return(0) }
    # Otherwise I'm an internal node:
    else {
      porder(tree$edge[desc[1], 2]) # visit my left daughter
      porder(tree$edge[desc[2], 2]) # visit my right daughter
    }

    return(0)
  }

  porder(root)

  tt <- tree
  tt$edge.length <- blens

  return(tt)
}
