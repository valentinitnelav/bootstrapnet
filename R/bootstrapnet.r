## This file documents the package.

#' @title
#' Bootstrap network indices
#'
#' @description
#' R package for resampling network metrics/indices (bootstrapping without
#' replacement). It wraps around the [bipartite][bip] package functions
#' \code{\link[bipartite]{networklevel}} and
#' \code{\link[bipartite]{specieslevel}}.
#'
#' Assuming a network/web like `Safariland` from the `bipartite` package, one
#' can sample interactions without replacement from the network until all
#' interactions have been used. The sampling procedure starts with a small
#' sample size to which interactions are added until all are consumed. Every
#' time we sample interactions, a smaller version of the entire network can be
#' built and a network index/metric can be computed. The sampling procedure can
#' be repeated as many times as needed, giving the possibility to compute mean
#' values with quantile-based confidence intervals. The mean values across
#' sample sizes can be plotted and indices of different networks can be visually
#' compared. See examples at
#' [README](https://github.com/valentinitnelav/plotbiomes/blob/develop/README.md)
#'
#' [bip]: https://cran.r-project.org/web/packages/bipartite/index.html
#'
#' @section Help:
#' Examples at: [README](https://github.com/valentinitnelav/plotbiomes/blob/develop/README.md)
#'
#' GitHub page [here](https://github.com/valentinitnelav/bootstrapnet)
#'
#' @author
#' Valentin Stefan, Tiffany Knight
#'
#' @md
#'
#' @docType package
#'
#' @name plotbiomes
NULL
