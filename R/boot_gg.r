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
#'  is 0.1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param alpha_boot_lines
#'  Alpha parameter of the line connecting the bootstrapped values of the
#'  metric. Default is 0.1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param size_ci
#'  Size of the line used for the confidence intervals. Default is 1. Passed to
#'  \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param linetype_ci
#'  Type of the line used for the confidence intervals. Default is 2. Passed to
#'  \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param size_mean
#'  Size of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param linetype_mean
#'  Type of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
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
#'  is 0.1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param alpha_boot_lines
#'  Alpha parameter of the line connecting the bootstrapped values of the
#'  metric. Default is 0.1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param size_ci
#'  Size of the line used for the confidence intervals. Default is 1. Passed to
#'  \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param linetype_ci
#'  Type of the line used for the confidence intervals. Default is 2. Passed to
#'  \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param size_mean
#'  Size of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param linetype_mean
#'  Type of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
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
#'  is 0.1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param alpha_boot_lines
#'  Alpha parameter of the line connecting the bootstrapped values of the
#'  metric. Default is 0.1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param size_ci
#'  Size of the line used for the confidence intervals. Default is 1. Passed to
#'  \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param linetype_ci
#'  Type of the line used for the confidence intervals. Default is 2. Passed to
#'  \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param size_mean
#'  Size of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
#'
#' @param linetype_mean
#'  Type of the line connecting the bootstrapped mean values of the metric.
#'  Default is 1. Passed to \code{\link[ggplot2:geom_path]{geom_line}}.
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
