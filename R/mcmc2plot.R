#' Plot a densi-tree from an MCMC sample
#'
#' Plot a densi-tree from an MCMC sample from a BPP or MCMCTree analysis
#'
#' @param tree an object of class phylo.
#' @param mcmc data frame with an MCMC sample from MCMCTree or a BPP A00
#'   analysis.
#' @param time.name character vector of length one.
#' @param thin numeric, the fraction of MCMC samples to keep.
#' @param col character, the color for branches.
#' @param alpha numeric, between 0 and 1, the branch color transparency.
#' @param y.offset numeric, the vertical offset for plotting the tree.
#' @param pfract, numeric, how much of the plotting space to used for plotting
#'   the tip labels. If \code{pfrac = 1}, the same amount of space is used for
#'   the tree and the labels. Use large values if your tip labels are long.
#' @param plot.labels, logical, whether to plot the tip labels. Ignored if
#'   \code{add = TRUE}.
#' @param axis, logical, whether to plot the x axis.
#' @param add logical, if TRUE add the trees to an existing plot, otherwise
#'   create a new plot.
#' @param tip.ages numeric, the ages of the tips, with the most recent tip
#'   having age zero, and the oldest tip having the largest age. If \code{NULL},
#'   tips are assumed to have all age zero.
#'
#' @details The function will reduce the MCMC sample to \code{dim(mcmc)[1] *
#'   thin} observations. Then the node ages in each observarion are used to plot
#'   each tree in the sample. For a tree with \code{s} species. The y
#'   coordinates of the tips are given by \code{0:(s - 1) + y.offset}.
#'
#'   The \code{tree} must be rooted, strictly bifurcating, and be the same tree
#'   used to genarate the BPP (A00) or MCMCTree MCMC samples.
#'
#' @author Mario dos Reis
#'
#' @examples
#' data(microcebus)
#' mcmc2densitree(microcebus$tree, microcebus$mcmc, time.name="tau_", thin=0.05,
#'  alpha=0.01, col="blue")
#'  title(xlab="Distance (substitutions per site)")
#'
#' data(hominids)
#' # Calibrate the hominid phylogeny with a uniform fossil calibration of
#' # between 6.5 to 10 Ma for the human-chimp divergence, and plot the
#' # calibrated sample
#' calmsc <- msc2time.t(mcmc=hominids$mcmc, node="7humanchimp", calf=runif,
#'   min=6.5, max=10)
#' mcmc2densitree(hominids$tree, calmsc, "t_", thin=0.05, alpha=0.01)
#' title(xlab="Divergence time (Ma)")
#'
#' @export
# FIXME: I'm very slow
mcmc2densitree <- function(tree, mcmc, time.name, thin, col="blue", alpha=1, y.offset=0,
                           pfract = 0.1, plot.labels = TRUE, axis=TRUE, add=FALSE, tip.ages=NULL) {
  ns <- tree$Nnode + 1
  if (is.null(tip.ages)) tip.ages <- rep(0, ns)
  ti <- grep(time.name, names(mcmc))
  n <- dim(mcmc)[1]; N <- n * thin
  ii <- floor(seq(from=1, to=n, length.out = N))

  # set up plotting area:
  ycoo <- .get.node.y(tree) + y.offset
  if (!add) {
    xend <- max(mcmc[ii,ti])
    xmin <- - xend * pfract
    plot(rep(0, 2*ns - 1), ycoo, xlim=c(xend, xmin), ylim=c(0, max(ycoo)), ty='n', axes=FALSE, xlab=NA, ylab=NA)
    if (plot.labels) text(x=0, y=ycoo[1:ns], tree$tip.label, pos=4)
    if (axis) axis(1)
  }

  # set colors:
  clv <- col2rgb(col) / 255
  clr <- rgb(clv[1,], clv[2,], clv[3,], alpha=alpha, maxColorValue = 1)

  # plot branch segments:
  for (i in ii) {
    ages <- c(tip.ages, unlist(mcmc[i,ti]))
    coo <- .get.segment.coo(tree, ages, ycoo)
    segments(coo$x0, coo$y0, coo$x1, coo$y1, col=clr)
  }
}

# calculate branch segment coordinates for plotting
# ages <- c(tip.ages, node.ages)
.get.segment.coo <- function(tree, ages, ycoo) {

  x0 <- ages[tree$edge[,1]]
  x1 <- ages[tree$edge[,2]]
  y0 <- ycoo[tree$edge[,1]]
  y1 <- ycoo[tree$edge[,2]]

  return(list(x0=x0, x1=x1, y0=y0, y1=y1))
}

# caculate y coordinates for nodes
.get.node.y <- function(tree) {
  ns <- tree$Nnode + 1
  y.node <- numeric(2 * ns - 1)

  # y coordinates for tips
  y.node[1:ns] <- 0:(ns - 1)

  # post-order traversal to calculate y coordinates
  # for internal nodes
  porder <- function(node) {
    desc <- which(tree$edge[,1] == node) # my descendant branches
    # left and right daughter nodes:
    ld <- tree$edge[desc[1], 2]
    rd <- tree$edge[desc[2], 2]

    # If I'm a tip:
    if (!length(desc)) { return(0) }

    # Otherwise I'm an internal node:
    else {
      porder(ld) # visit my left daughter
      porder(rd) # visit my right daughter
    }

    y.node[node] <<- (y.node[ld] + y.node[rd]) / 2

    return(0)
  }
  porder(ns + 1)

  return (y.node)
}

# Functions for plotting trees from MCMC samples of node ages
# mcmc.sum <- mcmc.summary(hominids$mcmc[,-1])
# mcmc2citree(hominids$tree, mcmc.sum, time.name="tau_")
# mcmc.sum <- mcmc.summary(microcebus$mcmc[,-1])
# mcmc2citree(microcebus$tree, mcmc.sum, time.name="tau_")

mcmc2citree <- function(tree, mcmc.sum, time.name) {
  ns <- tree$Nnode + 1
  tip.ages <- rep(0, ns)
  ns <- tree$Nnode + 1 # number of species
  ti <- grep(time.name, names(mcmc.sum$means))
  ycoo <- .get.node.y(tree)
  node.ages <- mcmc.sum$means[ti]
  ages <- c(tip.ages, node.ages)

  ycoo <- .get.node.y(tree)
  coo <- .get.segment.coo(tree, ages, ycoo)

  plot(ages, ycoo, pch=19, xlim=rev(range(ages)))
  segments(coo$x0, coo$y0, coo$x1, coo$y1)
}
