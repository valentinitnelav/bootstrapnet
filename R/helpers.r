#' @title
#' Transform web matrix to expanded data frame
#'
#' @description
#'  Transform a web matrix to an expanded data frame which is needed as input
#'  for the bootstrap functions. Each row (interaction) is repeated as many
#'  times as it was observed. This is important for the sampling procedure. Rows
#'  are shuffled.
#'
#' @param web
#'  A web matrix, e.g. `Safariland` from \code{\link{bipartite}}.
#'
#' @param seed
#'  Passed to `set.seed()`. Used for shuffling rows in a reproducible way.
#'
#' @return
#'  Returns an expanded data frame needed as input for the bootstrap functions.
#'  Rows are shuffled.
#'
#' @examples
#'
#' library(bootstrapnet)
#' library(bipartite)
#' data(Safariland)
#'
#' Safariland_df <- web_matrix_to_df(Safariland)
#'
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom dplyr filter slice n sample_frac
#'
#' @export
#'
#' @md
web_matrix_to_df <- function(web, seed = 42){
  if (class(web) != "matrix") stop("`web` must be a matrix, e.g. `bipartite::Safariland`")

  set.seed(seed)
  web %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "lower") %>%
    tidyr::gather(key = "higher", value = "counts", -"lower") %>%
    dplyr::filter(counts != 0) %>%
    dplyr::slice(rep(1:n(), times = counts)) %>%
    dplyr::sample_frac()
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
#'  which corresponds to a 95\\% confidence interval.
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
      setDT() %>%
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
      setDT() %>%
      data.table::melt(id.vars = "sp", variable.name = "spl_size", value.name = "mean")

    quantiles_low <- x %>%
      apply(MARGIN = 1:2, FUN = quantile, na.rm = TRUE, probs = probs[1]) %>%
      as.data.frame %>%
      tibble::rownames_to_column(var = "sp") %>%
      setDT() %>%
      data.table::melt(id.vars = "sp", variable.name = "spl_size", value.name = "ci_low")

    quantiles_up <- x %>%
      apply(MARGIN = 1:2, FUN = quantile, na.rm = TRUE, probs = probs[2]) %>%
      as.data.frame %>%
      tibble::rownames_to_column(var = "sp") %>%
      setDT() %>%
      data.table::melt(id.vars = "sp", variable.name = "spl_size", value.name = "ci_up")

    stats_df <- list(means, quantiles_low, quantiles_up) %>%
      purrr::reduce(dplyr::inner_join, by = c("sp", "spl_size")) %>%
      dplyr::mutate(spl_size = as.integer(as.character(spl_size)))


    lines_df <- x %>%
      apply(MARGIN = 3, FUN = as.data.frame) %>%
      lapply(rownames_to_column, var = "sp") %>%
      data.table::rbindlist(idcol = "simulation_id") %>%
      setDT() %>%
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


# Hidden helper functions -------------------------------------------------


#' @title
#' Splits a vector of indices `x` in `n` chunks.
#'
#' @description
#'  Splits a vector of indices `x` in `n` chunks. Inspired from
#'  https://stackoverflow.com/a/16275428/5193830
#'
#' @param x
#'  A vector of indices.
#'
#' @param n
#'  Number of chunks.
#'
#' @return
#'  A list of `n` vectors of indices.
#'
#' @noRd
split_in_chunks <- function(x, n) {
  split(x = x,
        f = cut(x = seq_along(x),
                breaks = n,
                labels = FALSE))
}


#' @title
#' Prepare a list of interaction indices for bootstrapping.
#'
#' @description
#'  Prepares a list of row indices from `data`. The first vector of indices (the
#'  first element of the list) has the length given by `start`. This is the
#'  smallest network of interactions sampled from `data`. Then the next vector
#'  of indices grows by `step`. The final vector contains all the interactions.
#'
#' @param data
#'  `data.table` of interactions from which to build and sample web matrices.
#'
#' @param start
#'  Integer. The sample size (number of interactions) of the first sampled
#'  network.
#'
#' @param step
#'  Integer. Sample size (number of interactions) used to increase gradually the
#'  sampled network until all interactions are sampled.
#'
#' @param seed
#'  Set seed to get reproducible random results. Passed to `set.seed()`.
#'
#' @return
#'  A list of vectors containing interaction indices.
#'
#' @noRd
sample_indices <- function(data, start, step, seed){
  set.seed(seed)
  # Shuffle the ids to sample from
  ids_shuffle <- sample(1:nrow(data))
  # Draw the start sample; set.seed() from above has effect here as well
  start_sample <- tryCatch(
    {
      sample(ids_shuffle, start)
    },
    error = function(cond){
      message(cond, "\n")
      stop("\nPossibly because your `start` value, ", start,
              ", is bigger than the number of interactions in data, ",
              unique(data) %>% dplyr::select(counts) %>% sum(),
              ".\nIf this is the case, then consider to reduce `start` to maybe 10 % of your interactions.")
    }
  )
  # Remove the already sampled ids used for the start sample
  ids_remaining <- ids_shuffle[! ids_shuffle %in% start_sample]

  # Build the list of growing sampled indices. The first element (first vector
  # of indices) was already sampled above (the start sample). Then add gradually
  # a new chunk of indices until all indices are used.
  ids_remaining_chunks <- tryCatch(
    {
      split_in_chunks(ids_remaining, length(ids_remaining)/step)
    },
    error = function(cond){
      message(cond, "\n")
      stop("\nPossibly because your `step` value, ", step,
           ", is bigger than the number of interactions in data, ",
           unique(data) %>% dplyr::select(counts) %>% sum(),
           ".\nIf this is the case, then consider to reduce `step` to maybe 5-10 % of your interactions.")
    }
  )
  ids_lst <- vector(mode = "list", length = (length(ids_remaining_chunks) + 1))
  ids_lst[[1]] <- start_sample
  for (i in 2:length(ids_lst)){
    ids_lst[[i]] <- c(ids_lst[[i-1]], ids_remaining_chunks[[i-1]])
  }

  return(ids_lst)
}


#' @title
#'  Helper function used in `boot_specieslevel_once()`. Gathers results in
#'  matrices.
#'
#' @description
#'  Prepares a list of matrices with the bootstrapped network metrics at species
#'  level.
#'
#' @param data
#'  `data.table` of interactions.
#'
#' @param col
#'  Quoted column name in `data` for lower trophic level species (plants) or for
#'  higher trophic level species (insects). See `boot_specieslevel_n()`.
#'
#' @param metric_lst
#'  A list of bootstrap outputs (data frames) from `bipartite::specieslevel`.
#'
#' @param level
#'  For which level should the level-specific indices be computed: 'both'
#'  (default), 'lower' or 'higher'?
#'
#' @return
#'  A list of matrices.
#'
#' @noRd
boot_specieslevel_lower_or_higher <- function(data, col, metric_lst, level){
  # Get species names from the web
  sp <- data[[col]] %>% unique %>% sort
  # Allocate memory for the matrix that will store bootstrap results. Rows
  # correspond to each species and columns correspond to each bootstrap
  # iteration. The last column will store metric results for the entire web.
  # Allocate memory also for a matrix that can contain the weighted version of
  # the metric if such a case is applicable (e.g. for index = "betweenness"
  # the output returned by bipartite::specieslevel is a two columns data
  # frame)
  n <- length(metric_lst)
  if (dim(metric_lst[[n]])[2] == 1) {
    boot_mat <- matrix(nrow = length(sp), ncol = length(metric_lst))
    rownames(boot_mat) <- sp
    # Because the first iterations most probably do not include all species
    # from the web, then must match species names in the matrix when storing
    # results.
    for (j in 1:n){
      indices_logical <- rownames(boot_mat) %in% rownames(metric_lst[[j]])
      boot_mat[indices_logical, j] <- metric_lst[[j]][, 1]
    }

    if (level == "lower") {
      return(list(lower_level = boot_mat))
    } else if (level == "higher") {
      return(list(higher_level = boot_mat))
    }

    # In case of weighted indices:
  } else if (dim(metric_lst[[n]])[2] == 2) {
    boot_mat <- boot_mat_weighted <- matrix(nrow = length(sp), ncol = length(metric_lst))
    rownames(boot_mat) <- rownames(boot_mat_weighted) <- sp
    for (j in 1:n){
      indices_logical <- rownames(boot_mat) %in% rownames(metric_lst[[j]])
      boot_mat[indices_logical, j] <- metric_lst[[j]][, 1]
      boot_mat_weighted[indices_logical, j] <- metric_lst[[j]][, 2]
    }

    if (level == "lower") {
      return(list(lower_level = boot_mat,
                  lower_level_weighted = boot_mat_weighted))
    } else if (level == "higher") {
      return(list(higher_level = boot_mat,
                  higher_level_weighted = boot_mat_weighted))
    }
  }
}

#' @title
#'  Helper function used in `boot_specieslevel_once()`. Gathers results in
#'  matrices.
#'
#' @description
#'  Similar to `boot_specieslevel_lower_or_higher` above, but with the extra
#'  layer of complexity of dealing with both levels at once.
#'
#' @param data
#'  `data.table` of interactions.
#'
#' @param col_lower
#'  Quoted column name in `data` for lower trophic level species (plants).
#'
#' @param col_higher
#'  Quoted column name in `data` for higher trophic level species (insects).
#'
#' @param metric_lst
#'  A list of bootstrap outputs (data frames) from `bipartite::specieslevel`.
#'
#' @return
#'  A list of matrices.
#'
#' @noRd
boot_specieslevel_both_levels <- function(data, col_lower, col_higher, metric_lst){
  # Get species names from the web
  sp_low <- data[[col_lower]] %>% unique %>% sort
  sp_high <- data[[col_higher]] %>% unique %>% sort

  n <- length(metric_lst)
  if (dim(metric_lst[[n]][[1]])[2] == 1) {
    boot_mat_lower_level <- matrix(nrow = length(sp_low), ncol = length(metric_lst))
    rownames(boot_mat_lower_level) <- sp_low

    boot_mat_higher_level <- matrix(nrow = length(sp_high), ncol = length(metric_lst))
    rownames(boot_mat_higher_level) <- sp_high

    for (j in 1:n){
      low_indices_logical <- rownames(boot_mat_lower_level) %in% rownames(metric_lst[[j]]$`lower level`)
      boot_mat_lower_level[low_indices_logical, j] <- metric_lst[[j]]$`lower level`[, 1]

      high_indices_logical <- rownames(boot_mat_higher_level) %in% rownames(metric_lst[[j]]$`higher level`)
      boot_mat_higher_level[high_indices_logical, j] <- metric_lst[[j]]$`higher level`[, 1]
    }

    return(list(higher_level = boot_mat_higher_level,
                lower_level = boot_mat_lower_level))

    # In case of weighted indices:
  } else if (dim(metric_lst[[n]][[1]])[2] == 2) {
    boot_mat_lower_level <- boot_mat_lower_level_weighted <- matrix(nrow = length(sp_low), ncol = length(metric_lst))
    rownames(boot_mat_lower_level) <- rownames(boot_mat_lower_level_weighted) <- sp_low

    boot_mat_higher_level <- boot_mat_higher_level_weighted <- matrix(nrow = length(sp_high), ncol = length(metric_lst))
    rownames(boot_mat_higher_level) <-  rownames(boot_mat_higher_level_weighted) <- sp_high

    for (j in 1:n){
      low_indices_logical <- rownames(boot_mat_lower_level) %in% rownames(metric_lst[[j]]$`lower level`)
      boot_mat_lower_level[low_indices_logical, j] <- metric_lst[[j]]$`lower level`[, 1]
      boot_mat_lower_level_weighted[low_indices_logical, j] <- metric_lst[[j]]$`lower level`[, 2]

      high_indices_logical <- rownames(boot_mat_higher_level) %in% rownames(metric_lst[[j]]$`higher level`)
      boot_mat_higher_level[high_indices_logical, j] <- metric_lst[[j]]$`higher level`[, 1]
      boot_mat_higher_level_weighted[high_indices_logical, j] <- metric_lst[[j]]$`higher level`[, 2]
    }

    return(list(higher_level = boot_mat_higher_level,
                higher_level_weighted = boot_mat_higher_level_weighted,
                lower_level = boot_mat_lower_level,
                lower_level_weighted = boot_mat_lower_level_weighted))
  }
}


#' @title
#'  Helper function used in `boot_specieslevel_n()`. Gathers results in an array
#'  of matrices.
#'
#' @description
#'  Simplifies the results from the parallel processing into a list of arrays.
#'
#' @param boot_lst
#'  A list of matrices. Is the output of the `foreach` loop from
#'  `boot_specieslevel_n`.
#'
#' @param metric_names
#'  The name(s) of the network metric. There can be multiple names, depending on
#'  the `index` and `level` values.
#'
#' @param n
#'  Number of metric names.
#'
#' @param iter_spl_size
#'  Sample size at each iteration/bootstrap step.
#'
#' @return
#'  A list of arrays of matrices.
#'
#' @noRd
get_list_of_arrays <- function(boot_lst, metric_names, n, iter_spl_size){
  array_lst <- vector(mode = "list", length = n)
  for (i in 1:n) {
    array_lst[[i]] <- boot_lst %>%
      lapply("[[", metric_names[i]) %>%
      simplify2array
    colnames(array_lst[[i]]) <- iter_spl_size
  }
  names(array_lst) <- metric_names
  return(array_lst)
}


#' @title
#'  Helper function used in `boot_networklevel_n()`. Gathers results in a list
#'  of matrices.
#'
#' @description
#'  Simplifies the results from the parallel processing into a list of matrices.
#'
#' @param boot_lst
#'  A list of matrices. Is the output of the `foreach` loop from
#'  `boot_networklevel_n`.
#'
#' @param metric_names
#'  The name(s) of the network metric. There can be multiple names, depending on
#'  the `index` and `level` values.
#'
#' @param n
#'  Number of metric names.
#'
#' @param iter_spl_size
#'  Sample size at each iteration/bootstrap step.
#'
#' @return
#'  A list of matrices.
#'
#' @noRd
get_list_of_matrices <- function(boot_lst, metric_names, n, iter_spl_size){
  mat_lst <- vector(mode = "list", length = n)
  for (i in 1:n) {
    mat_lst[[i]] <- boot_lst %>%
      lapply("[[", metric_names[i]) %>%
      do.call(what = "cbind")
    rownames(mat_lst[[i]]) <- iter_spl_size
  }
  names(mat_lst) <- metric_names
  return(mat_lst)
}


#' @title
#'  Helper used in `gg_specieslevel_web_by_web` and
#'  `gg_specieslevel_compare_webs`
#'
#' @description
#'  Filters `stats_df` and `lines_df` by the proper lower and higher level
#'  species names.
#'
#' @noRd
filter_species <- function(lst, sp_lower, sp_higher, name_metric, name_plot) {

  is_lower <- grepl(pattern = "lower", x = name_metric)
  is_higher <- grepl(pattern = "higher", x = name_metric)

  if (is_lower & is.null(sp_lower)) {
    warning("`sp_lower` not provided, all lower level species will be used!")
    sp_filter <- lst[["stats_df"]][["sp"]] %>% unique()
    stats_df <- lst[["stats_df"]] %>% dplyr::filter(sp %in% sp_filter)
    lines_df <- lst[["lines_df"]] %>% dplyr::filter(sp %in% sp_filter)

  } else if (is_lower & !is.null(sp_lower)) {
    sps <- lst[["stats_df"]][["sp"]] %>% unique()
    sp_test <- ! sp_lower %in% sps
    if (any(sp_test))
      warning("At case ", name_plot, ", ",
              "these species are not in the lower level species of the network:\n",
              paste(sp_lower[sp_test], collapse = "\n"))
    stats_df <- lst[["stats_df"]] %>% dplyr::filter(sp %in% sp_lower)
    lines_df <- lst[["lines_df"]] %>% dplyr::filter(sp %in% sp_lower)

  } else if (is_higher & is.null(sp_higher)) {
    warning("`sp_higher` not provided, all higher level species will be used!")
    sp_filter <- lst[["stats_df"]][["sp"]] %>% unique()
    stats_df <- lst[["stats_df"]] %>% dplyr::filter(sp %in% sp_filter)
    lines_df <- lst[["lines_df"]] %>% dplyr::filter(sp %in% sp_filter)

  } else if (is_higher & !is.null(sp_higher)) {
    sps <- lst[["stats_df"]][["sp"]] %>% unique()
    sp_test <- ! sp_higher %in% sps
    if (any(sp_test))
      warning("At case ", name_plot, ", ",
              "these species are not in the higher level species of the network:\n",
              paste(sp_higher[sp_test], collapse = "\n"))
    stats_df <- lst[["stats_df"]] %>% dplyr::filter(sp %in% sp_higher)
    lines_df <- lst[["lines_df"]] %>% dplyr::filter(sp %in% sp_higher)
  }

  return(list(stats_df = stats_df,
              lines_df = lines_df))
}


#' @title
#'  Helper used in `gg_specieslevel_compare_webs`.
#'
#' @description
#'  Checks if the values for `sp_lower` and `sp_higher` parameters are ok and
#'  throws errors otherwise.
#'
#' @noRd
test_species_filter_compare_webs <- function(metric_names, sp_lower, sp_higher){
  if (is.null(sp_lower) & is.null(sp_higher))
    stop("Please provide a species name to compare in `sp_lower` and/or `sp_higher`")

  if (!is.null(sp_lower) & length(sp_lower) != 1)
    stop("Please provide a single common species name in `sp_lower` to compare across the networks")

  if (!is.null(sp_higher) & length(sp_higher) != 1)
    stop("Please provide a single common species name in `sp_higher` to compare across the networks")

  has_lower <- grepl(pattern = "lower", x = metric_names) %>% any()
  has_higher <- grepl(pattern = "higher", x = metric_names) %>% any()

  if (has_lower & is.null(sp_lower))
    stop("Please provide a lower species name in `sp_lower`")

  if (has_higher & is.null(sp_higher))
    stop("Please provide a higher species name in `sp_higher`")
}


#' @title
#'  Helper used in `boot_networklevel_n()` and `boot_specieslevel_n()`.
#'
#' @description
#'  Checks if the row and column values of `data` are ok and throws errors
#'  otherwise.
#'
#' @noRd
test_data_species_names <- function(data, col_lower, col_higher) {
  if (any(c("", "NA", "na", NA) %in% unique(data[[col_lower]])))
    stop("You have undefined/empty species names. Check the lower level species names for NA-s or empty strings.")

  if (any(c("", "NA", "na", NA) %in% unique(data[[col_higher]])))
    stop("You have undefined/empty species names. Check the higher level species names for NA-s or empty strings.")
}


#' @title
#'  Helper used in `boot_networklevel_n()` and `boot_specieslevel_n()`.
#'
#' @description
#'  Checks if:
#'  - the `data` argument is a non-empty `data.frame` or `data.table` class,
#'  - data should have at least 3 columns (e.g.: plants, insects, counts),
#'  - the columns names should match the user input for `col_lower` and `col_higher`,
#'  - each interaction (row) was repeated as many times as it was observed
#'
#' @noRd
test_data <- function(data, col_lower, col_higher) {
  if ( ! inherits(data, what = c("data.frame", "data.table")) )
    strwrap('`data` must be a data.frame or data.table.
             Maybe check out the function web_matrix_to_df() for examples.',
            prefix = " ", initial = "") %>% stop()

  if ( dim(data)[1] == 0 )
    strwrap('Your `data` has no rows.
             Maybe check out the function web_matrix_to_df() for examples.',
            prefix = " ", initial = "") %>% stop()

  if ( dim(data)[2] < 3 )
    strwrap('`data` must have at least 3 columns, e.g: plants, insects, counts.
             Maybe check out the function web_matrix_to_df() for examples.',
            prefix = " ", initial = "") %>% stop()

  if ( ! col_lower %in% colnames(data) )
    paste0('Cannot find the column name provided in `col_lower` as a column of `data`.
            Maybe a typo? You gave: <', col_lower, '>, but columns in data are: ',
           paste(colnames(data), collapse = ", ")) %>%
    strwrap(prefix = " ", initial = "") %>% stop()

  if ( ! col_higher %in% colnames(data) )
    paste0('Cannot find the column name provided in `col_higher` as a column of `data`.
            Maybe a typo? You gave: <', col_higher, '>, but columns in data are: ',
           paste(colnames(data), collapse = ", ")) %>%
    strwrap(prefix = " ", initial = "") %>% stop()

  if ( nrow(data) == unique(data) %>% nrow() )
    strwrap('Each interaction (row) should be repeated as many times as it was observed.
             If you indeed observed each interaction just once, then carry on.
             Maybe check out the function web_matrix_to_df() for examples.',
            prefix = " ", initial = "") %>% warning()
}


#' @title
#'  Helper used in `boot_networklevel_n()` and `boot_specieslevel_n()`.
#'
#' @description
#'  Checks if the user gave correct values for `level`.
#'
#' @noRd
test_level_value <- function(level) {
  if (length(level) != 1)
    stop('Provide a single value for `level`: "both", "lower" or "higher"')

  if (! level %in% c("both", "lower", "higher"))
    stop('The value provided for `level` must be one of the following: "both", "lower" or "higher"')
}


#' @title
#'  Helper used in `boot_networklevel_n()`.
#'
#' @description
#'  Checks if the user gave correct values for `index`.
#'
#' @noRd
test_index_networklevel <- function(index) {
  if (length(index) != 1) stop('Provide a single value for `index`')

  # As given in https://github.com/biometry/bipartite/blob/master/bipartite/R/networklevel.R at 2019-11-07
  allindex <- c(
    # descriptive:
    "number of species",
    "connectance",
    "web asymmetry",
    # binary based and related:
    "links per species",
    "number of compartments",
    "compartment diversity",
    "cluster coefficient",
    "degree distribution",
    "mean number of shared partners",
    "togetherness",
    "C score",
    "V ratio",
    "discrepancy",
    "nestedness",
    "NODF",
    "weighted nestedness",
    # miscelleneous:
    "ISA",
    "SA",
    "extinction slope",
    "robustness",
    "niche overlap",
    # quantitative series:
    "weighted cluster coefficient",
    "weighted NODF",
    "partner diversity",
    "generality",
    "vulnerability",
    "linkage density",
    "weighted connectance",
    "Fisher alpha",
    "interaction evenness",
    "Alatalo interaction evenness",
    "Shannon diversity",
    "functional complementarity",
    "H2"
  )

  if (! index %in% allindex)
    stop('Your index is not recognised. Typo? Did you mean one of these?:\n',
         agrep(pattern = index,
               x = allindex,
               ignore.case = TRUE,
               value = TRUE,
               fixed = FALSE) %>%
           paste(collapse = "\n"))
}


#' @title
#'  Helper used in `boot_specieslevel_n()`.
#'
#' @description
#'  Checks if the user gave correct values for `index`.
#'
#' @noRd
test_index_specieslevel <- function(index) {
  if (length(index) != 1) stop('Provide a single value for `index`')

  # As given in https://github.com/biometry/bipartite/blob/master/bipartite/R/specieslevel.R at 2019-11-07
  allindex <- c(
    "degree",
    "normalised degree",
    "species strength",
    "nestedrank",
    "interaction push pull",
    "PDI",
    "resource range",
    "species specificity",
    "PSI",
    "NSI",
    "betweenness",
    "closeness",
    "Fisher alpha",
    "partner diversity",
    "effective partners",
    "d",
    "dependence",
    "proportional generality",
    "proportional similarity"
  )

  if (! index %in% allindex)
    stop('Your index is not recognised. Typo? Did you mean one of these?:\n',
         agrep(pattern = index,
               x = allindex,
               ignore.case = TRUE,
               value = TRUE,
               fixed = FALSE) %>%
           paste(collapse = "\n"))
}
