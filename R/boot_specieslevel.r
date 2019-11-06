#' @title
#'  Prepare for plotting the bootstrapped metrics at the species level of one or
#'  multiple networks.
#'
#' @description
#'  Takes a list of network interactions (each interaction being repeated as
#'  many times as it was observed), and bootstraps the given species level
#'  metric (`index`) for each network. Runs `boot_specieslevel_n()` for a list
#'  of network interactions and prepares the data for plotting with `ggplot`.
#'  The output list can be passed to
#'  \code{\link[bootstrapnet]{gg_specieslevel_compare_webs}} or to
#'  \code{\link[bootstrapnet]{get_stats_multi}} and then its output to
#'  \code{\link[bootstrapnet]{gg_specieslevel_web_by_web}}. See examples below.
#'
#' @param lst
#'  A list of one or multiple data frames of interactions from which to build
#'  and sample web matrices. Each interaction (row in the data frame) must be
#'  repeated as many times as it was observed. E.g. if the interaction species_1
#'  x species_2 was observed 5 times, then repeat that row 5 times within the
#'  data frame.
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
#'  Consider to set `start` to maybe 10\% of your total unique interactions.
#'
#' @param step
#'  Integer. Sample size (number of interactions) used to increase gradually the
#'  sampled network until all interactions are sampled. If `step` is too small
#'  (e.g. 1) then the computation time is very long depending on your total
#'  number of interactions from which samples are taken. Consider to set `step`
#'  to maybe 5-10\% of your total unique interactions.
#'
#' @param n_boot
#'  Number of desired bootstraps (50 or 100 can be enough).
#'
#' @param n_cpu
#'  Number of CPU-s to use for parallel processing.
#'
#' @param probs
#'  A numeric vector of two probabilities in `[0, 1]`. Passed to
#'  \code{\link[matrixStats]{rowQuantiles}} and used for building the lower and
#'  upper bounds of the confidence intervals. Defaults to `c(0.025, 0.975)`,
#'  which corresponds to a 95\% confidence interval.
#'
#' @param ...
#'  Other arguments passed to \code{\link[bipartite]{specieslevel}}.
#'
#' @return
#'  A 3 levels list of the same length as `lst`. The list can be passed to
#'  \code{\link[bootstrapnet]{gg_specieslevel_compare_webs}} or to
#'  \code{\link[bootstrapnet]{get_stats_multi}} and then its output to
#'  \code{\link[bootstrapnet]{gg_specieslevel_web_by_web}}. Each element of the
#'  list is a sub-list corresponding to each network. Each sub-list contains
#'  further sub-lists with the bootstrapped metric values at 'lower', 'higher'
#'  or 'both' levels. Each sub-sub-lists contains two data frames: `stats_df`
#'  and `lines_df`, which can be used by the `ggplot2::geom_line()` function.
#'  See the return section of \code{\link[bootstrapnet]{get_stats_single}} for
#'  more details about `stats_df` and `lines_df` data frames.
#'
#' @examples
#'
#' library(bootstrapnet)
#' library(bipartite)
#' library(magrittr)
#' data(Safariland)
#'
#' set.seed(321)
#' Safariland_1 <- Safariland[, sort(sample.int(ncol(Safariland), 20))]
#' set.seed(123)
#' Safariland_2 <- Safariland[, sort(sample.int(ncol(Safariland), 20))]
#'
#' lst <- list(s1 = Safariland_1, s2 = Safariland_2) %>%
#'   lapply(web_matrix_to_df) %>%
#'   boot_specieslevel(col_lower = "lower", # column name for plants
#'                     col_higher = "higher", # column name for insects
#'                     index = "betweenness",
#'                     level = "both",
#'                     start = 50,
#'                     step = 20,
#'                     n_boot = 10,
#'                     n_cpu = 2)
#'
#' lst %>%
#'   get_stats_multi() %>%
#'   gg_specieslevel_compare_webs(sp_lower = "Alstroemeria aurea",
#'                                sp_higher = "Allograpta.Toxomerus")
#'
#' lst %>%
#'   gg_specieslevel_web_by_web(sp_lower = c("Alstroemeria aurea",
#'                                           "Aristotelia chilensis"))
#'
#' @export
#'
#' @md
boot_specieslevel <- function(lst,
                              col_lower,
                              col_higher,
                              index,
                              level,
                              start,
                              step,
                              n_boot,
                              n_cpu,
                              probs = c(0.025, 0.975),
                              ...){

  webs_stats <- vector(mode = "list", length = length(lst))
  names(webs_stats) <- names(lst)

  for (i in 1:length(webs_stats)){
    webs_stats[[i]] <- lst[[i]] %>%
      boot_specieslevel_n(col_lower = col_lower,
                          col_higher = col_higher,
                          index = index,
                          level = level,
                          start = start,
                          step = step,
                          n_boot = n_boot,
                          n_cpu = n_cpu,
                          ...) %>%
      lapply(FUN = get_stats_single, probs = probs)
  }

  return(webs_stats)
}


#' @title
#'  Bootstrap metric at the species level multiple times.
#'
#' @description
#'  Bootstrap a single network of interactions n times in parallel. This
#'  function is used for network metrics at the species level. Starts with a
#'  small sample size of interactions (e.g. `start = 30`), builds the
#'  corresponding web matrix/network, computes its metric (e.g. `index = "d"`)
#'  using \code{\link[bipartite]{specieslevel}}, then adds new interactions
#'  (e.g. `step = 20`) until all interactions are sampled. The last sample is
#'  actually the entire network. Repeats n times (as given in `n_boot`) these
#'  steps in parallel on multiple CPUs.
#'
#' @param data
#'  Data frame of interactions from which to build and sample web matrices. Each
#'  interaction (row in the data frame) must be repeated as many times as it was
#'  observed. E.g. if species_1 x species_2 was observed 5 times, then repeat
#'  that row 5 times within the data frame. See examples below.
#'
#' @param col_lower
#'  Quoted column name in `data` for lower trophic level species (plants).
#'
#' @param col_higher
#'  Quoted column name in `data` for higher trophic level species (insects).
#'
#' @param index
#'  Passed to \code{\link[bipartite]{specieslevel}}. See
#'  `?bipartite::specieslevel` for details.
#'
#' @param level
#'  Passed to \code{\link[bipartite]{specieslevel}}. See
#'  `?bipartite::specieslevel` for details. For which level should the
#'  level-specific indices be computed: 'both' (default), 'lower' or 'higher'?
#'
#' @param start
#'  Integer. The sample size (number of interactions) to start the bootstrap
#'  with. If the start sample size is small (e.g. 5 or 10), then first
#'  iterations might results in NaN-s and warning messages are displayed.
#'  Consider to set `start` to maybe 10\% of your total unique interactions.
#'
#' @param step
#'  Integer. Sample size (number of interactions) used to increase gradually the
#'  sampled network until all interactions are sampled. If `step` is too small
#'  (e.g. 1) then the computation time is very long depending on your total
#'  number of interactions from which samples are taken. Consider to set `step`
#'  to maybe 5-10\% of your total unique interactions.
#'
#' @param n_boot
#'  Number of desired bootstraps (50 or 100 can be enough).
#'
#' @param n_cpu
#'  Number of CPU-s to use for parallel processing.
#'
#' @param ...
#'  Other arguments passed to \code{\link[bipartite]{specieslevel}} like
#'  `logbase`, etc.
#'
#' @return
#'  Returns a list of 1, 2 or 4 arrays of matrices. The species names are stored
#'  as the row names of each matrix. The number of columns of a matrix indicates
#'  how many iterations took place. This is decided internally based on the given
#'  values to `start`, `step` and the total number of rows (interactions) in
#'  `data`. The column names give the sample size at each iteration. The last
#'  iteration (last column name) is always the entire network (total number of
#'  interactions in `data`). The 3rd dimension (number of matrices in the array)
#'  corresponds to `n_boot` (number of bootstraps).
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
#'   boot_specieslevel_n(col_lower = "lower", # column name for plants
#'                       col_higher = "higher", # column name for insects
#'                       index = "d",
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
boot_specieslevel_n <- function(data,
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
  if (! "data.table" %in% class(data)) data.table::setDT(data)

  if (any(c("", "NA", "na") %in% unique(data[[col_lower]])))
    stop("You have undefined/empty species names. Check the lower level species names for NA-s or empty strings.")

  if (any(c("", "NA", "na") %in% unique(data[[col_higher]])))
    stop("You have undefined/empty species names. Check the higher level species names for NA-s or empty strings.")

  # Get sample sizes. Column names of the resulting bootstrap matrices will
  # carry information about the sample size at each iteration/bootstrap step.
  # This is run also before the parallel processing initiation because it can
  # throw error messages if the start and step values are not adequate. No need
  # to initiate parallel processing for something that will fail.
  iter_spl_size <- sample_indices(data = data, start = start, step = step, seed = 42) %>%
    sapply(length)

  { # Start parallel processing
    chunks <- parallel::splitIndices(n_boot, n_cpu)
    cl <- parallel::makeCluster(n_cpu)
    doParallel::registerDoParallel(cl)
    i <- NULL # to avoid 'Undefined global functions or variables: i' in R CMD check

    boot_lst <- foreach::foreach(i = iterators::iter(chunks),
                                 .errorhandling = 'pass',
                                 .packages = c("data.table",
                                               "bipartite",
                                               "magrittr"),
                                 .export = c("boot_specieslevel_once",
                                             "split_in_chunks",
                                             "sample_indices",
                                             "boot_specieslevel_lower_or_higher",
                                             "boot_specieslevel_both_levels")) %dopar%
                                             {
                                               lapply(i, FUN = function(x)
                                                 boot_specieslevel_once(data = data,
                                                                        col_lower = col_lower,
                                                                        col_higher = col_higher,
                                                                        index = index,
                                                                        level = level,
                                                                        start = start,
                                                                        step = step,
                                                                        seed = x,
                                                                        ...) )
                                             }
    parallel::stopCluster(cl)
    remove(cl)
    foreach::registerDoSEQ()
    } # End of parallel processing

  boot_lst <- unlist(boot_lst, recursive = FALSE)

  metric_names <- names(boot_lst[[1]])
  n <- length(metric_names)

  if (n == 1) {
    results <- get_list_of_arrays(boot_lst, metric_names, n, iter_spl_size)
    # Convert data back to data frame if applicable
    if (! "data.table" %in% cls_data) data.table::setDF(data)
    return(results)

  } else if (n == 2) {
    results <- get_list_of_arrays(boot_lst, metric_names, n, iter_spl_size)
    if (! "data.table" %in% cls_data) data.table::setDF(data)
    return(results)

  } else if (n == 4) {
    results <- get_list_of_arrays(boot_lst, metric_names, n, iter_spl_size)
    if (! "data.table" %in% cls_data) data.table::setDF(data)
    return(results)
  }
}


#' @title
#'  Bootstrap metric at the species level once.
#'
#' @description
#'  You will rarely use this function alone. It was designed to be executed in
#'  parallel by \code{\link[bootstrapnet]{boot_specieslevel_n}}. See more
#'  details there.
#'
#' @param data
#'  See \code{\link[bootstrapnet]{boot_specieslevel_n}}.
#'
#' @param col_lower
#'  See \code{\link[bootstrapnet]{boot_specieslevel_n}}.
#'
#' @param col_higher
#'  See \code{\link[bootstrapnet]{boot_specieslevel_n}}.
#'
#' @param index
#'  See \code{\link[bootstrapnet]{boot_specieslevel_n}}.
#'
#' @param level
#'  See \code{\link[bootstrapnet]{boot_specieslevel_n}}.
#'
#' @param start
#'  See \code{\link[bootstrapnet]{boot_specieslevel_n}}.
#'
#' @param step
#'  See \code{\link[bootstrapnet]{boot_specieslevel_n}}.
#'
#' @param seed
#'  Set seed to get reproducible random results. Passed to `set.seed()`.
#'
#' @param ...
#'  See \code{\link[bootstrapnet]{boot_specieslevel_n}}.
#'
#' @return
#'  A list of 1, 2 or 4 matrices. This depends on the provided `index` metric
#'  and `level` value.
#'
#' @references
#'  This function is a wrapper of \code{\link[bipartite]{specieslevel}}.
#'
#' @import data.table
#'
#' @importFrom bipartite specieslevel
#'
#' @export
#'
#' @md
boot_specieslevel_once <- function(data,
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

  # Allocate memory for a list of outputs (data frames) from
  # bipartite::specieslevel()
  metric_lst <- vector(mode = "list", length = length(ids_lst))

  # Sample interactions with the indices build above and compute the species
  # level metric until all interactions are sampled. In case of warnings or
  # errors (most probably related to small sample size of a network at first
  # iterations), then return NA-s for the metric so it doesn't break the
  # simulation.
  for (i in 1:length(ids_lst)){
    metric_lst[[i]] <- try({
      # Some metrics may have random processes in their computation, so set seed
      # here also for reproducibility.
      web <- data[ids_lst[[i]], table(get(col_lower), get(col_higher))]
      set.seed(42)
      bipartite::specieslevel(web = web, index = index, level = level, ...)
    })
  }

  if (level == "lower") {
    boot_results <- boot_specieslevel_lower_or_higher(data = data,
                                                      col = col_lower,
                                                      metric_lst = metric_lst,
                                                      level)
    return(boot_results)

  } else if (level == "higher") {
    boot_results <- boot_specieslevel_lower_or_higher(data = data,
                                                      col = col_higher,
                                                      metric_lst = metric_lst,
                                                      level)
    return(boot_results)

  } else if (level == "both") {
    boot_results <- boot_specieslevel_both_levels(data = data,
                                                  col_lower = col_lower,
                                                  col_higher = col_higher,
                                                  metric_lst = metric_lst)
    return(boot_results)
  }
}
