#' @title
#'  Bootstrap network level metric multiple times.
#'
#' @description
#'  Bootstrap a single network of interactions n times in parallel and collects
#'  the network level metrics. Starts with a small sample size of interactions
#'  (e.g. `start = 30`), builds the corresponding web matrix/network, computes
#'  its metric (e.g. `index = "nestedness"`) using
#'  \code{\link[bipartite]{networklevel}}, then adds new interactions (e.g.
#'  `step = 20`) until all interactions are sampled. The last sample is actually
#'  the entire network. Repeats n times (as given in `n_boot`) these steps in
#'  parallel on multiple CPUs.
#'
#' @param data
#'  Data frame of interactions from which to build and sample web matrices. Each
#'  interaction (row in the data frame) must be repeated as many times as it was
#'  observed. E.g. if the interaction species_1 x species_2 was observed 5
#'  times, then repeat that row 5 times within the data frame. See examples
#'  below.
#'
#' @param col_lower
#'  Quoted column name in `data` for lower trophic level species (plants).
#'
#' @param col_higher
#'  Quoted column name in `data` for higher trophic level species (insects).
#'
#' @param index
#'  Passed to \code{\link[bipartite]{networklevel}}. See
#'  `?bipartite::networklevel` for details.
#'
#' @param level
#'  Passed to \code{\link[bipartite]{networklevel}}. See
#'  `?bipartite::networklevel` for details. For which level should the
#'  level-specific indices be computed: 'both' (default), 'lower' or 'higher'?
#'
#' @param start
#'  Integer. The sample size (number of interactions) to start the bootstrap
#'  with. If the start sample size is small (e.g. 5 or 10), then first
#'  iterations might results in NaN-s and warning messages are displayed.
#'
#' @param step
#'  Integer. Sample size (number of interactions) used to increase gradually the
#'  sampled network until all interactions are sampled. If `step` is too small
#'  (e.g. 1) then the computation time is very long depending on your total
#'  number of interactions from which samples are taken.
#'
#' @param n_boot
#'  Number of desired bootstraps (50 or 100 can be enough).
#'
#' @param n_cpu
#'  Number of CPU-s to use for parallel processing.
#'
#' @param ...
#'  Other arguments passed to \code{\link[bipartite]{networklevel}} like
#'  `logbase`, etc.
#'
#' @return
#'  Returns a list of a single matrix or a list of two matrices depending on the
#'  provided `index` metric and `level` value. For example if
#'  `index='nestedness'` and `level='both'`, then it returns a list of a single
#'  matrix. But in the case of `index='niche overlap'`, it returns a list of two
#'  matrices, first matrix (`niche.overlap.HL`) corresponding to the higher level
#'  (e.g. insects), and the second (`niche.overlap.LL`) for the lower level (e.g.
#'  plants). The number of rows of a matrix indicates how many iterations took
#'  place. This is decided internally based on the given values to `start`,
#'  `step` and the total number of rows (interactions) in `data`.The row names
#'  give the sample size at each iteration. The last iteration (last row name) is
#'  always the entire network (total number of interactions in `data`). The
#'  number of columns corresponds to `n_boot` (number of bootstraps).
#'
#' @examples
#'
#' library(bootstrapnet)
#' library(magrittr)
#' library(bipartite)
#' data(Safariland)
#'
#' Safariland %>%
#'   web_matrix_to_df() %>%
#'   boot_networklevel_n(col_lower = "lower", # column name for plants
#'                       col_higher = "higher", # column name for insects
#'                       index = "niche overlap",
#'                       level = "both",
#'                       start = 100,
#'                       step = 100,
#'                       n_boot = 10,
#'                       n_cpu = 2)
#'
#' @import data.table
#'
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom parallel splitIndices makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators iter
#'
#' @export
#'
#' @md
boot_networklevel_n <- function(data,
                                col_lower,
                                col_higher,
                                index,
                                level,
                                start,
                                step,
                                n_boot,
                                n_cpu,
                                ...){
  cls_data <- class(data)
  if (! "data.table" %in% cls_data) data.table::setDT(data)

  if (any(c("", "NA", "na") %in% unique(data[[col_lower]])))
    stop("You have undefined/empty species names. Check the lower level species names for NA-s or empty strings.")

  if (any(c("", "NA", "na") %in% unique(data[[col_higher]])))
    stop("You have undefined/empty species names. Check the higher level species names for NA-s or empty strings.")


  { # Start parallel processing
    chunks <- parallel::splitIndices(n_boot, n_cpu)
    cl <- parallel::makeCluster(n_cpu)
    doParallel::registerDoParallel(cl)

    i <- NULL # to avoid 'Undefined global functions or variables: i' in R CMD check

    boot_lst <-
      foreach::foreach(i = iterators::iter(chunks),
                       .packages = c("magrittr",
                                     "data.table",
                                     "bipartite"),
                       .export = c("boot_networklevel_once",
                                   "split_in_chunks",
                                   "sample_indices")) %dopar% {
                                     lapply(i, FUN = function(x) # note the lapply!
                                       boot_networklevel_once(data = data,
                                                              col_lower = col_lower,
                                                              col_higher = col_higher,
                                                              index = index,
                                                              level = level,
                                                              start = start,
                                                              step = step,
                                                              seed = x,
                                                              ...)
                                     )
                                   }
    parallel::stopCluster(cl)
    remove(cl)
    foreach::registerDoSEQ()

    } # End of parallel processing

  # Row names of the resulting bootstrap matrices will carry information about
  # the sample size at each iteration/bootstrap step.
  iter_spl_size <- sample_indices(data = data, start = start, step = step, seed = 42) %>%
    sapply(length)

  boot_lst <- unlist(boot_lst, recursive = FALSE)

  metric_names <- names(boot_lst[[1]])
  n <- length(metric_names)

  if (n == 1) {
    results <- get_list_of_matrices(boot_lst, metric_names, n, iter_spl_size)
    # Convert data back to data frame if applicable
    if (! "data.table" %in% cls_data) data.table::setDF(data)
    return(results)

  } else if (n == 2) {
    results <- get_list_of_matrices(boot_lst, metric_names, n, iter_spl_size)
    if (! "data.table" %in% cls_data) data.table::setDF(data)
    return(results)
  }
}


#' @title
#'  Bootstrap network level metric once.
#'
#' @description
#'  You will rarely use this function alone. It was designed to be executed in
#'  parallel by \code{\link[bootstrapnet]{boot_networklevel_n}}. See more
#'  details there.
#'
#' @param data
#'  See \code{\link[bootstrapnet]{boot_networklevel_n}}.
#'
#' @param col_lower
#'  See \code{\link[bootstrapnet]{boot_networklevel_n}}.
#'
#' @param col_higher
#'  See \code{\link[bootstrapnet]{boot_networklevel_n}}.
#'
#' @param index
#'  See \code{\link[bootstrapnet]{boot_networklevel_n}}.
#'
#' @param level
#'  See \code{\link[bootstrapnet]{boot_networklevel_n}}.
#'
#' @param start
#'  See \code{\link[bootstrapnet]{boot_networklevel_n}}.
#'
#' @param step
#'  See \code{\link[bootstrapnet]{boot_networklevel_n}}.
#'
#' @param seed
#'  Set seed to get reproducible random results. Passed to `set.seed()`.
#'
#' @param ...
#'  See \code{\link[bootstrapnet]{boot_networklevel_n}}.
#'
#' @return
#'  Returns a data frame of one or two columns depending on the provided `index`
#'  metric and `level` value. For example if `index='nestedness'` and
#'  `level='both'`, then it returns a one column data frame. But in case of
#'  `index='niche overlap'`, it returns a two columns data frame, first column
#'  ('niche.overlap.HL') corresponding to the higher level (e.g. insects), and
#'  the second column ('niche.overlap.LL') for the lower level (e.g. plants).
#'
#' @references
#'  This function is a wrapper of \code{\link[bipartite]{networklevel}}.
#'
#' @import data.table
#'
#' @importFrom magrittr %>%
#' @importFrom bipartite networklevel
#'
#' @export
#'
#' @md
boot_networklevel_once <- function(data,
                                   col_lower,
                                   col_higher,
                                   index,
                                   level,
                                   start,
                                   step,
                                   seed,
                                   ...){

  if (! "data.table" %in% class(data)) data.table::setDT(data)

  # List of sampled indices as vectors for bootstrapping. Each vector of indices
  # is used to build a web matrix/network from data.
  ids_lst <- sample_indices(data = data, start = start, step = step, seed = seed)

  # Allocate memory for the objects that will grow during iterations/bootstrapping.
  metric_lst <- vector(mode = "list", length = length(ids_lst))

  # Sample interactions with the indices build above and compute the network
  # level metric until all interactions are sampled. In case of errors (most
  # probably related to small sample size of a network at first iterations),
  # then `try()` returns NaN-s for the metric so it doesn't break the
  # bootstrapping process.
  for (i in 1:length(ids_lst)){
    metric_lst[[i]] <- try({
      # Some metrics, like nestedness have random processes in their
      # computation, so set seed here also for reproducibility.
      web <- data[ids_lst[[i]], table(get(col_lower), get(col_higher))]
      set.seed(42)
      bipartite::networklevel(web = web, index = index, level = level, ...)
    })
  }

  # Prepare results as data frame.
  df <- metric_lst %>%
    lapply(rbind) %>%
    lapply(as.data.frame) %>%
    data.table::rbindlist()
  setDF(df)

  return(df)
}
