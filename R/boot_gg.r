#' @title
#'  Prepare for plotting the bootstrapped network level metric results of one or
#'  multiple networks.
#'
#' @description
#'  Takes a list of network interactions (each interaction being repeated as
#'  many times as it was observed), and bootstraps the network level metric
#'  (`index`) for each network. Runs `boot_networklevel_n()` for a list of
#'  network interactions and prepares the data for plotting with `ggplot`. The
#'  output list can be passed to \code{\link[bootstrapnet]{gg_networklevel}}.
#'  See examples below.
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
#'  The name of the network level metric. Passed to
#'  \code{\link[bipartite]{networklevel}}. See `?bipartite::networklevel` for
#'  details.
#'
#' @param level
#'  For which level should the level-specific indices be computed: 'both'
#'  (default), 'lower' or 'higher'? Passed to
#'  \code{\link[bipartite]{networklevel}}. See `?bipartite::networklevel` for
#'  details.
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
#' @param probs
#'  A numeric vector of two probabilities in `[0, 1]`. Passed to
#'  \code{\link[matrixStats]{rowQuantiles}} and used for building the lower and
#'  upper bounds of the confidence intervals. Defaults to `c(0.025, 0.975)`,
#'  which corresponds to a 95\% confidence interval.
#'
#' @param ...
#'  Other arguments passed to \code{\link[bipartite]{networklevel}} like
#'  `logbase`, etc.
#'
#' @return
#'  A list of one (when `level = 'lower'` or `level = 'higher'`) or two
#'  sub-lists (when `level = 'both'`). The list can be passed to
#'  \code{\link[bootstrapnet]{gg_networklevel}}. Each sub-list, contains two
#'  data frames: `stats_df` and `lines_df`, which can be used by the
#'  `ggplot2::geom_line()` function. See the return section of
#'  \code{\link[bootstrapnet]{get_stats_single}} for more details about
#'  `stats_df` and `lines_df` data frames.
#'
#' @examples
#'
#' library(bootstrapnet)
#' library(bipartite)
#' library(magrittr)
#' data(Safariland)
#'
#' set.seed(321)
#' Safariland_1 <- Safariland[, sort(sample.int(ncol(Safariland), 10))]
#' set.seed(123)
#' Safariland_2 <- Safariland[, sort(sample.int(ncol(Safariland), 10))]
#'
#' lst <- list(s1 = Safariland_1, s2 = Safariland_2) %>%
#'   lapply(web_matrix_to_df) %>%
#'   boot_networklevel(col_lower = "lower", # column name for plants
#'                     col_higher = "higher", # column name for insects
#'                     index = "nestedness",
#'                     level = "both",
#'                     start = 10,
#'                     step = 10,
#'                     n_boot = 10,
#'                     n_cpu = 2)
#' gg_networklevel(lst)
#'
#' @export
#'
#' @md
boot_networklevel <- function(lst,
                              col_lower,
                              col_higher,
                              index,
                              level,
                              start,
                              step,
                              n_boot,
                              n_cpu,
                              probs = c(0.025, 0.975),
                              ...) {

  webs_stats <- vector(mode = "list", length = length(lst))
  names(webs_stats) <- names(lst)

  for (i in 1:length(webs_stats)){
    webs_stats[[i]] <- lst[[i]] %>%
      boot_networklevel_n(col_lower = col_lower,
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

  webs_stats %>%
    get_stats_multi() %>%
    return()
}


#' @title
#'  Plot bootstrapped network level metrics for one or multiple networks.
#'
#' @description
#'  This function takes as input the output of
#'  \code{\link[bootstrapnet]{boot_networklevel}}. Can be used to plot the
#'  bootstrapped network level metrics of a single network (not that useful) or
#'  for multiple networks (to check for overlaps in their confidence intervals).
#'
#' @param metrics_stats
#'  The output of \code{\link[bootstrapnet]{boot_networklevel}}.
#'
#' @param size_boot_lines
#'  Size of the line connecting the bootstrapped values of the metric. Default
#'  is 0.1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @param alpha_boot_lines
#'  Alpha parameter of the line connecting the bootstrapped values of the
#'  metric. Default is 0.1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @param size_ci
#'  Size of the line used for the confidence intervals. Default is 1. Passed to
#'  \code{\link[ggplot2]{geom_line}}.
#'
#' @param linetype_ci
#'  Type of the line used for the confidence intervals. Default is 2. Passed to
#'  \code{\link[ggplot2]{geom_line}}.
#'
#' @param size_mean
#'  Size of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @param linetype_mean
#'  Type of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @return
#'  A list of one or more ggplots.
#'
#' @examples
#' # See example section of ?boot_networklevel
#'
#' @import ggplot2
#'
#' @export
#'
#' @md
gg_networklevel <- function(metrics_stats,
                            size_boot_lines = 0.1,
                            alpha_boot_lines = 0.1,
                            size_ci = 1,
                            linetype_ci = 2,
                            size_mean = 1,
                            linetype_mean = 1){

  n <- length(metrics_stats)
  plots_gg <- vector(mode = "list", length = n)
  names(plots_gg) <- names(metrics_stats)

  for (i in 1:n){
    stats_df <- metrics_stats[[i]]$stats_df
    lines_df <- metrics_stats[[i]]$lines_df

    plots_gg[[i]] <- ggplot() +
      # Bootstrap lines
      geom_line(data = lines_df,
                aes(x = spl_size,
                    y = value,
                    group = interaction(web, simulation_id),
                    color = web),
                size = size_boot_lines,
                alpha = alpha_boot_lines) +
      # Lower confidence interval bound
      geom_line(data = stats_df,
                aes(x = spl_size,
                    y = ci_low,
                    group = web,
                    color = web),
                size = size_ci,
                linetype = linetype_ci) +
      # Upper confidence interval bound
      geom_line(data = stats_df,
                aes(x = spl_size,
                    y = ci_up,
                    group = web,
                    color = web),
                size = size_ci,
                linetype = linetype_ci) +
      # Average line
      geom_line(data = stats_df,
                aes(x = spl_size,
                    y = mean,
                    group = web,
                    color = web),
                size = size_mean,
                linetype = linetype_mean) +
      ylab(names(plots_gg[i]))
  }

  return(plots_gg)
}


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
#'  Plot bootstrapped species level metrics for one or more webs/networks.
#'
#' @description
#'  To be used in combination with \code{\link[bootstrapnet]{boot_specieslevel}}
#'  and \code{\link[bootstrapnet]{get_stats_multi}}. See the example section of
#'  \code{\link[bootstrapnet]{boot_specieslevel}}.
#'
#' @param sp_lower
#'  A character vector of one or more lower level species names (plants) shared
#'  among the given networks. Only the bootstrapped metric values for these
#'  species will be plotted. It can be helpful when there are too many species
#'  in a network and the graphs could get easily cluttered. If `sp_lower` is
#'  left empty, then all species will be considered for plotting. If the species
#'  names are not shared across all given networks, then a warning will be given
#'  for the specific network from where they are missing and the plot will be
#'  empty or show only existing shared species names.
#'
#' @param sp_higher
#'  Similar to `sp_lower`, but valid for higher level species names (insects).
#'
#' @param metrics_stats
#'  The output of \code{\link[bootstrapnet]{get_stats_multi}}.
#'
#' @param size_boot_lines
#'  Size of the line connecting the bootstrapped values of the metric. Default
#'  is 0.1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @param alpha_boot_lines
#'  Alpha parameter of the line connecting the bootstrapped values of the
#'  metric. Default is 0.1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @param size_ci
#'  Size of the line used for the confidence intervals. Default is 1. Passed to
#'  \code{\link[ggplot2]{geom_line}}.
#'
#' @param linetype_ci
#'  Type of the line used for the confidence intervals. Default is 2. Passed to
#'  \code{\link[ggplot2]{geom_line}}.
#'
#' @param size_mean
#'  Size of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @param linetype_mean
#'  Type of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @return
#'  A list of one or more ggplots.
#'
#' @examples
#'  # See example section of ?boot_specieslevel
#'
#' @import ggplot2
#'
#' @export
#'
#' @md
gg_specieslevel_web_by_web <- function(metrics_stats,
                                       sp_lower = NULL,
                                       sp_higher = NULL,
                                       size_boot_lines = 0.1,
                                       alpha_boot_lines = 0.1,
                                       size_ci = 1,
                                       linetype_ci = 2,
                                       size_mean = 1,
                                       linetype_mean = 1){

  metrics_stats <- unlist(metrics_stats, recursive = FALSE)
  n <- length(metrics_stats)

  plots_gg <- vector(mode = "list", length = n)
  names(plots_gg) <- names_web_metrics <- names(metrics_stats)

  for (i in 1:n){
    lst_df <- filter_species(lst = metrics_stats[[i]],
                             sp_lower = sp_lower,
                             sp_higher = sp_higher,
                             name_metric = names_web_metrics[i],
                             name_plot = names(plots_gg[i]))
    stats_df <- lst_df[["stats_df"]]
    lines_df <- lst_df[["lines_df"]]

    plots_gg[[i]] <- ggplot() +
      # Bootstrap lines
      geom_line(data = lines_df,
                aes(x = spl_size,
                    y = value,
                    group = interaction(sp, simulation_id),
                    color = sp),
                size = size_boot_lines,
                alpha = alpha_boot_lines) +
      # Lower confidence interval bound
      geom_line(data = stats_df,
                aes(x = spl_size,
                    y = ci_low,
                    color = sp),
                size = size_ci,
                linetype = linetype_ci) +
      # Upper confidence interval bound
      geom_line(data = stats_df,
                aes(x = spl_size,
                    y = ci_up,
                    color = sp),
                size = size_ci,
                linetype = linetype_ci) +
      # Average line
      geom_line(data = stats_df,
                aes(x = spl_size,
                    y = mean,
                    color = sp),
                size = size_mean,
                linetype = linetype_mean) +
      ylab(names(plots_gg[i]))
  }

  return(plots_gg)
}


#' @title
#'  Plot bootstrapped species level metrics for a common species in two or more
#'  webs.
#'
#' @description
#'  Allows the visual comparison of the bootstrapped metric values of a selected
#'  species (lower and higher level) across two or more networks. To be used in
#'  combination with \code{\link[bootstrapnet]{boot_specieslevel}}. See the
#'  example section of \code{\link[bootstrapnet]{boot_specieslevel}}.
#'
#' @param sp_lower
#'  A character vector of a single lower level species name (plant) shared among
#'  the given networks. Only the bootstrapped metric values for this single
#'  species will be plotted.
#'
#' @param sp_higher
#'  Similar to `sp_lower`, but valid for higher level species names (insects).
#'
#' @param metrics_stats
#'  The output of \code{\link[bootstrapnet]{get_stats_multi}}.
#'
#' @param size_boot_lines
#'  Size of the line connecting the bootstrapped values of the metric. Default
#'  is 0.1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @param alpha_boot_lines
#'  Alpha parameter of the line connecting the bootstrapped values of the
#'  metric. Default is 0.1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @param size_ci
#'  Size of the line used for the confidence intervals. Default is 1. Passed to
#'  \code{\link[ggplot2]{geom_line}}.
#'
#' @param linetype_ci
#'  Type of the line used for the confidence intervals. Default is 2. Passed to
#'  \code{\link[ggplot2]{geom_line}}.
#'
#' @param size_mean
#'  Size of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @param linetype_mean
#'  Type of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2]{geom_line}}.
#'
#' @return
#'  A list of one or more ggplots.
#'
#' @examples
#'  # See example section of ?boot_specieslevel
#'
#' @import ggplot2
#'
#' @export
#'
#' @md
gg_specieslevel_compare_webs <- function(metrics_stats,
                                         sp_lower = NULL,
                                         sp_higher = NULL,
                                         size_boot_lines = 0.1,
                                         alpha_boot_lines = 0.1,
                                         size_ci = 1,
                                         linetype_ci = 2,
                                         size_mean = 1,
                                         linetype_mean = 1){


  n <- length(metrics_stats)
  plots_gg <- vector(mode = "list", length = n)
  names(plots_gg) <- names_web_metrics <- names(metrics_stats)

  test_species_filter_compare_webs(metric_names = names_web_metrics,
                                   sp_lower = sp_lower,
                                   sp_higher = sp_higher)

  for (i in 1:n){
    lst_df <- filter_species(lst = metrics_stats[[i]],
                             sp_lower = sp_lower,
                             sp_higher = sp_higher,
                             name_metric = names_web_metrics[i],
                             name_plot = names(plots_gg[i]))
    stats_df <- lst_df[["stats_df"]]
    lines_df <- lst_df[["lines_df"]]

    plots_gg[[i]] <- ggplot() +
      # Bootstrap lines
      geom_line(data = lines_df,
                aes(x = spl_size,
                    y = value,
                    group = interaction(web, simulation_id),
                    color = web),
                size = size_boot_lines,
                alpha = alpha_boot_lines) +
      # Lower confidence interval bound
      geom_line(data = stats_df,
                aes(x = spl_size,
                    y = ci_low,
                    color = web),
                size = size_ci,
                linetype = linetype_ci) +
      # Upper confidence interval bound
      geom_line(data = stats_df,
                aes(x = spl_size,
                    y = ci_up,
                    color = web),
                size = size_ci,
                linetype = linetype_ci) +
      # Average line
      geom_line(data = stats_df,
                aes(x = spl_size,
                    y = mean,
                    color = web),
                size = size_mean,
                linetype = linetype_mean) +
      ylab(paste(names(plots_gg[i]), unique(stats_df[["sp"]]), sep = ": "))
  }

  return(plots_gg)
}


#' @title
#' Prepare bootstrapped results of a single network for `ggplot` friendly data
#' format.
#'
#' @description
#'  Helper to prepare the bootstrapped metric output of a single web/network
#'  from \code{\link[bootstrapnet]{boot_networklevel_n}} or
#'  \code{\link[bootstrapnet]{boot_specieslevel_n}} in a format ready for
#'  \code{\link[ggplot2]{ggplot}}.
#'
#' @param x
#'  A bootstrap matrix or an array of matrices from
#'  \code{\link[bootstrapnet]{boot_networklevel_n}} or
#'  \code{\link[bootstrapnet]{boot_specieslevel_n}}.
#'
#' @param probs
#'  A numeric vector of two probabilities in `[0, 1]`. Passed to
#'  \code{\link[matrixStats]{rowQuantiles}} and used for building the lower and
#'  upper bounds of the confidence intervals. Defaults to `c(0.025, 0.975)`,
#'  which corresponds to a 95\% confidence interval.
#'
#' @return
#'  A list of two data frames to be used with \code{\link[ggplot2]{ggplot}}.
#'  First one, `stats_df` is a 4 columns data frame, containing the sample size
#'  (`spl_size`) at each iteration and the corresponding mean metric (`mean`)
#'  together with its bootstrap confidence interval limits (`ci_low`, `ci_up`).
#'  The second one, `lines_df`, contains all the bootstrapped values (`value`) of
#'  a given network metric at each iteration (`simulation_id`) and its
#'  corresponding sample size (`spl_size`). Can be used for enhancing visual
#'  effect when plotting the mean bootstrap values.
#'
#' @examples
#' library(bootstrapnet)
#' library(bipartite)
#' library(magrittr)
#' data(Safariland)
#'
#' stats_d_lower <- Safariland %>%
#'   web_matrix_to_df() %>%
#'   boot_specieslevel_n(col_lower = "lower", # column name for plants
#'                         col_higher = "higher", # column name for insects
#'                         index = "d",
#'                         level = "both",
#'                         start = 100,
#'                         step = 100,
#'                         n_boot = 10,
#'                         n_cpu = 2) %>%
#'   .[["lower_level"]] %>%
#'   get_stats_single()
#'
#' @importFrom magrittr %>%
#' @importFrom data.table melt rbindlist
#' @importFrom dplyr mutate inner_join
#' @importFrom tibble rownames_to_column
#' @importFrom purrr reduce
#' @importFrom matrixStats rowMeans2 rowQuantiles
#'
#' @export
#'
#' @md
get_stats_single <- function(x, probs = c(0.025, 0.975)){
  cls <- class(x)
  if (! cls %in% c("matrix", "array")) stop("Expecting a matrix or a 3 dimensions array")

  if (cls == "matrix") {
    stats_df <- data.frame(spl_size = as.integer(rownames(x))) %>%
      dplyr::mutate(mean   = matrixStats::rowMeans2(x, na.rm = TRUE),
                    ci_low = matrixStats::rowQuantiles(x, probs = probs[1], na.rm = TRUE) %>% unname,
                    ci_up  = matrixStats::rowQuantiles(x, probs = probs[2], na.rm = TRUE) %>% unname)

    lines_df <- x %>%
      as.data.frame %>%
      tibble::rownames_to_column(var = "spl_size") %>%
      data.table::melt(id.vars = "spl_size",
                       variable.name = "simulation_id") %>%
      dplyr::mutate(spl_size = as.integer(spl_size),
                    simulation_id = as.integer(simulation_id))

    return(list(stats_df = stats_df,
                lines_df = lines_df))

  } else if (cls == "array") {
    means <- x %>%
      apply(MARGIN = 1:2, FUN = mean, na.rm = TRUE) %>%
      as.data.frame %>%
      tibble::rownames_to_column(var = "sp") %>%
      data.table::melt(id.vars = "sp", variable.name = "spl_size", value.name = "mean")

    quantiles_low <- x %>%
      apply(MARGIN = 1:2, FUN = quantile, na.rm = TRUE, probs = probs[1]) %>%
      as.data.frame %>%
      tibble::rownames_to_column(var = "sp") %>%
      data.table::melt(id.vars = "sp", variable.name = "spl_size", value.name = "ci_low")

    quantiles_up <- x %>%
      apply(MARGIN = 1:2, FUN = quantile, na.rm = TRUE, probs = probs[2]) %>%
      as.data.frame %>%
      tibble::rownames_to_column(var = "sp") %>%
      data.table::melt(id.vars = "sp", variable.name = "spl_size", value.name = "ci_up")

    stats_df <- list(means, quantiles_low, quantiles_up) %>%
      purrr::reduce(dplyr::inner_join, by = c("sp", "spl_size")) %>%
      dplyr::mutate(spl_size = as.integer(as.character(spl_size)))


    lines_df <- x %>%
      apply(MARGIN = 3, FUN = as.data.frame) %>%
      lapply(rownames_to_column, var = "sp") %>%
      data.table::rbindlist(idcol = "simulation_id") %>%
      data.table::melt(id.vars = c("sp", "simulation_id"), variable.name = "spl_size") %>%
      dplyr::mutate(spl_size = as.integer(as.character(spl_size)))

    return(list(stats_df = stats_df,
                lines_df = lines_df))
  }
}

#' @title
#'  Prepare bootstrapped results of one or more networks for `ggplot2` friendly
#'  data format.
#'
#' @description
#'  Helper to prepare the bootstrapped metric output of one or more
#'  webs/networks from \code{\link[bootstrapnet]{boot_specieslevel}},
#'  \code{\link[bootstrapnet]{boot_networklevel_n}} or
#'  \code{\link[bootstrapnet]{boot_specieslevel_n}} in a format ready for
#'  \code{\link[ggplot2]{ggplot}}. See the examples section of
#'  \code{\link[bootstrapnet]{boot_specieslevel}}.
#'
#' @param webs_stats
#'  A list of bootstrap results from
#'  \code{\link[bootstrapnet]{boot_specieslevel}}, or from multiple runs of
#'  \code{\link[bootstrapnet]{boot_networklevel_n}} or
#'  \code{\link[bootstrapnet]{boot_specieslevel_n}}.
#'
#' @return
#'  A list of two data frames to be used with \code{\link[ggplot2]{ggplot}}.
#'  First one, `stats_df` is a 5 columns data frame, containing for each network
#'  (`web`) the sample size (`spl_size`) at each iteration and the corresponding
#'  mean metric (`mean`), together with its bootstrap confidence interval limits
#'  (`ci_low`, `ci_up`). The second one, `lines_df`, contains all the
#'  bootstrapped values (`value`) of a given network metric at each iteration
#'  (`simulation_id`) and its corresponding sample size (`spl_size`). Can be
#'  used for enhancing visual effect when plotting the mean bootstrap values.
#'
#' @examples
#'  # See example section of ?boot_specieslevel
#'
#' @export
#'
#' @md
get_stats_multi <- function(webs_stats){

  metric_names <- names(webs_stats[[1]])
  n <- length(metric_names)

  metrics_stats <- vector(mode = "list", length = n)
  names(metrics_stats) <- metric_names

  for (i in 1:n){
    temp <- webs_stats %>%
      lapply("[[", metric_names[i])

    temp_stats_df <- temp %>%
      lapply("[[", "stats_df") %>%
      data.table::rbindlist(idcol = "web") %>%
      data.table::setDF()

    temp_lines_df <- temp %>%
      lapply("[[", "lines_df") %>%
      data.table::rbindlist(idcol = "web") %>%
      data.table::setDF()

    metrics_stats[[i]] <- list(stats_df = temp_stats_df,
                               lines_df = temp_lines_df)

    rm(temp, temp_stats_df, temp_lines_df)
  }

  return(metrics_stats)
}
