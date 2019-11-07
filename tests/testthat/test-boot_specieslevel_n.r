context("Species level indices - test boot_specieslevel_n")

library(magrittr)
library(bipartite)
data(Safariland)


# Helper function to get the length of unique values from the last columns of
# each matrix of the given array.
get_length_unq <- function(test_array) {
  test_array %>%
    .[, dim(.)[2],] %>%
    apply(MARGIN = 1, FUN = function(x) length(unique(x))) %>%
    unique(.)
}


# Test d ------------------------------------------------------------------

# Note that "d" doesn’t have an weighted version.

index <- "d"

# For full reproducibly is safer to reorder the Safariland matrix
# alphabetically.
sf_sort <- Safariland %>%
  .[, order(colnames(.))] %>%
  .[order(rownames(.)), ]

set.seed(42)
test_both <- specieslevel(sf_sort, index = index, level = "both")


# ~ both levels -----------------------------------------------------------

test_that('"d" index output is ok for both levels', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_specieslevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "both",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 2) # one array for each level

  expect_match(names(sf_boot)[1], "higher_level", fixed = TRUE)
  expect_match(names(sf_boot)[2], "lower_level", fixed = TRUE)

  expect_true(inherits(sf_boot[["higher_level"]], 'array'))
  expect_true(inherits(sf_boot[["lower_level"]], 'array'))

  expect_true(length(dim(sf_boot[["higher_level"]])) == 3)
  expect_true(length(dim(sf_boot[["lower_level"]])) == 3)

  expect_true(dim(sf_boot[["higher_level"]])[3] == 6) # because n_boot = 6
  expect_true(dim(sf_boot[["lower_level"]])[3] == 6) # because n_boot = 6

  expect_true(typeof(sf_boot[["higher_level"]]) == "double")
  expect_true(typeof(sf_boot[["lower_level"]]) == "double")

  # Are all last columns of the matrices identical?
  unq_length_hl <- get_length_unq(sf_boot[["higher_level"]])
  expect_equal(unq_length_hl, 1)

  unq_length_ll <- get_length_unq(sf_boot[["lower_level"]])
  expect_equal(unq_length_ll, 1)

  # Is the last bootstrap value correct? As per bipartite::specieslevel ?
  final_boot_val_hl <- sf_boot[["higher_level"]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val_hl, test_both[["higher level"]])

  final_boot_val_ll <- sf_boot[["lower_level"]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val_ll, test_both[["lower level"]])
})


# ~ higher level ----------------------------------------------------------

test_that('"d" index output is ok for higher level', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_specieslevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "higher",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 1)

  expect_match(names(sf_boot), "higher_level", fixed = TRUE)

  expect_true(inherits(sf_boot[[1]], 'array'))

  expect_true(length(dim(sf_boot[[1]])) == 3)

  expect_true(dim(sf_boot[[1]])[3] == 6) # because n_boot = 6

  expect_true(typeof(sf_boot[[1]]) == "double")

  # Are all last columns of the matrices identical?
  unq_length <- get_length_unq(sf_boot[[1]])
  expect_equal(unq_length, 1)

  # Is the last bootstrap value correct? As per bipartite::specieslevel ?
  final_boot_val <- sf_boot[[1]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val, test_both[["higher level"]])
})


# ~ lower level -----------------------------------------------------------

test_that('"d" index output is ok for lower level', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_specieslevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "lower",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 1)

  expect_match(names(sf_boot), "lower_level", fixed = TRUE)

  expect_true(inherits(sf_boot[[1]], 'array'))

  expect_true(length(dim(sf_boot[[1]])) == 3)

  expect_true(dim(sf_boot[[1]])[3] == 6) # because n_boot = 6

  expect_true(typeof(sf_boot[[1]]) == "double")

  # Are all last columns of the matrices identical?
  unq_length <- get_length_unq(sf_boot[[1]])
  expect_equal(unq_length, 1)

  # Is the last bootstrap value correct? As per bipartite::specieslevel ?
  final_boot_val <- sf_boot[[1]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val, test_both[["lower level"]])
})


# Test betweenness --------------------------------------------------------

# Note that "betweenness" has an weighted version. That is,
# bipartite::specieslevel() to return a two column data frame - "betweenness" &
# "weighted.betweenness" (see below).

index <- "betweenness"

# For full reproducibly is safer to reorder the Safariland matrix
# alphabetically.
sf_sort <- Safariland %>%
  .[, order(colnames(.))] %>%
  .[order(rownames(.)), ]

set.seed(42)
test_both <- specieslevel(sf_sort, index = index, level = "both")


# ~ both levels -----------------------------------------------------------

test_that('"betweenness" index output is ok for both levels', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_specieslevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "both",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 4) # one array for each level

  expect_match(names(sf_boot)[1], "higher_level", fixed = TRUE)
  expect_match(names(sf_boot)[2], "higher_level_weighted", fixed = TRUE)
  expect_match(names(sf_boot)[3], "lower_level", fixed = TRUE)
  expect_match(names(sf_boot)[4], "lower_level_weighted", fixed = TRUE)

  expect_true(inherits(sf_boot[["higher_level"]], 'array'))
  expect_true(inherits(sf_boot[["higher_level_weighted"]], 'array'))
  expect_true(inherits(sf_boot[["lower_level"]], 'array'))
  expect_true(inherits(sf_boot[["lower_level_weighted"]], 'array'))

  expect_true(length(dim(sf_boot[["higher_level"]])) == 3)
  expect_true(length(dim(sf_boot[["higher_level_weighted"]])) == 3)
  expect_true(length(dim(sf_boot[["lower_level"]])) == 3)
  expect_true(length(dim(sf_boot[["lower_level_weighted"]])) == 3)

  expect_true(dim(sf_boot[["higher_level"]])[3] == 6) # because n_boot = 6
  expect_true(dim(sf_boot[["higher_level_weighted"]])[3] == 6)
  expect_true(dim(sf_boot[["lower_level"]])[3] == 6)
  expect_true(dim(sf_boot[["lower_level_weighted"]])[3] == 6)

  expect_true(typeof(sf_boot[["higher_level"]]) == "double")
  expect_true(typeof(sf_boot[["higher_level_weighted"]]) == "double")
  expect_true(typeof(sf_boot[["lower_level"]]) == "double")
  expect_true(typeof(sf_boot[["lower_level_weighted"]]) == "double")

  # Are all lastcolumns of the matrices identical?
  unq_length_hl <- get_length_unq(sf_boot[["higher_level"]])
  expect_equal(unq_length_hl, 1)

  unq_length_hlw <-  get_length_unq(sf_boot[["higher_level_weighted"]])
  expect_equal(unq_length_hlw, 1)

  unq_length_ll <-  get_length_unq(sf_boot[["lower_level"]])
  expect_equal(unq_length_ll, 1)

  unq_length_llw <-  get_length_unq(sf_boot[["lower_level_weighted"]])
  expect_equal(unq_length_llw, 1)

  # Is the last bootstrap value correct? As per bipartite::specieslevel ?
  final_boot_val_hl <- sf_boot[["higher_level"]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val_hl, test_both[["higher level"]][, "betweenness", drop = FALSE])

  final_boot_val_hlw <- sf_boot[["higher_level_weighted"]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val_hlw, test_both[["higher level"]][, "weighted.betweenness", drop = FALSE])

  final_boot_val_ll <- sf_boot[["lower_level"]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val_ll, test_both[["lower level"]][, "betweenness", drop = FALSE])

  final_boot_val_llw <- sf_boot[["lower_level_weighted"]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val_llw, test_both[["lower level"]][, "weighted.betweenness", drop = FALSE])
})


# ~ higher level ----------------------------------------------------------

test_that('"betweenness" index output is ok for higher level', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_specieslevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "higher",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 2) # one array for each level

  expect_match(names(sf_boot)[1], "higher_level", fixed = TRUE)
  expect_match(names(sf_boot)[2], "higher_level_weighted", fixed = TRUE)

  expect_true(inherits(sf_boot[["higher_level"]], 'array'))
  expect_true(inherits(sf_boot[["higher_level_weighted"]], 'array'))

  expect_true(length(dim(sf_boot[["higher_level"]])) == 3)
  expect_true(length(dim(sf_boot[["higher_level_weighted"]])) == 3)

  expect_true(dim(sf_boot[["higher_level"]])[3] == 6) # because n_boot = 6
  expect_true(dim(sf_boot[["higher_level_weighted"]])[3] == 6)

  expect_true(typeof(sf_boot[["higher_level"]]) == "double")
  expect_true(typeof(sf_boot[["higher_level_weighted"]]) == "double")

  # Are all lastcolumns of the matrices identical?
  unq_length_hl <- get_length_unq(sf_boot[["higher_level"]])
  expect_equal(unq_length_hl, 1)

  unq_length_hlw <-  get_length_unq(sf_boot[["higher_level_weighted"]])
  expect_equal(unq_length_hlw, 1)

  # Is the last bootstrap value correct? As per bipartite::specieslevel ?
  final_boot_val_hl <- sf_boot[["higher_level"]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val_hl, test_both[["higher level"]][, "betweenness", drop = FALSE])

  final_boot_val_hlw <- sf_boot[["higher_level_weighted"]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val_hlw, test_both[["higher level"]][, "weighted.betweenness", drop = FALSE])
})


# ~ lower level -----------------------------------------------------------

test_that('"betweenness" index output is ok for lower level', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_specieslevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "lower",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 2) # one array for each level

  expect_match(names(sf_boot)[1], "lower_level", fixed = TRUE)
  expect_match(names(sf_boot)[2], "lower_level_weighted", fixed = TRUE)

  expect_true(inherits(sf_boot[["lower_level"]], 'array'))
  expect_true(inherits(sf_boot[["lower_level_weighted"]], 'array'))

  expect_true(length(dim(sf_boot[["lower_level"]])) == 3)
  expect_true(length(dim(sf_boot[["lower_level_weighted"]])) == 3)

  expect_true(dim(sf_boot[["lower_level"]])[3] == 6)
  expect_true(dim(sf_boot[["lower_level_weighted"]])[3] == 6)

  expect_true(typeof(sf_boot[["lower_level"]]) == "double")
  expect_true(typeof(sf_boot[["lower_level_weighted"]]) == "double")

  # Are all lastcolumns of the matrices identical?
  unq_length_ll <-  get_length_unq(sf_boot[["lower_level"]])
  expect_equal(unq_length_ll, 1)

  unq_length_llw <-  get_length_unq(sf_boot[["lower_level_weighted"]])
  expect_equal(unq_length_llw, 1)

  # Is the last bootstrap value correct? As per bipartite::specieslevel ?
  final_boot_val_ll <- sf_boot[["lower_level"]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val_ll, test_both[["lower level"]][, "betweenness", drop = FALSE])

  final_boot_val_llw <- sf_boot[["lower_level_weighted"]] %>% .[, ncol(.), 1] %>% as.data.frame
  expect_equivalent(final_boot_val_llw, test_both[["lower level"]][, "weighted.betweenness", drop = FALSE])
})
