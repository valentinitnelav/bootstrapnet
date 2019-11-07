context("Network level indices - test boot_networklevel_n")

library(magrittr)
library(bipartite)
data(Safariland)


# Helper function to get the length of unique values from the last row of the
# given matrix.
get_length_unq <- function(test_matrix) {
  test_matrix %>%
    .[nrow(.), ] %>%
    unique() %>%
    length()
}


# Test nestedness ---------------------------------------------------------

# Note that, "nestedness" is not be affected by `level` as opposed to "niche
# overlap".

index <- "nestedness"

# For full reproducibly is safer to reorder the Safariland matrix
# alphabetically.
sf_sort <- Safariland %>%
  .[, order(colnames(.))] %>%
  .[order(rownames(.)), ]

set.seed(42)
test_both <- networklevel(sf_sort, index = index, level = "both")
set.seed(42)
test_higher <- networklevel(sf_sort, index = index, level = "higher")
set.seed(42)
test_lower <- networklevel(sf_sort, index = index, level = "lower")


# ~ both levels -----------------------------------------------------------

test_that('"nestedness" index output is ok for both levels', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_networklevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "both",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  # Is there just a single matrix in the output list?
  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 1)
  expect_match(names(sf_boot), "nestedness", fixed = TRUE)
  expect_true(inherits(sf_boot[["nestedness"]], 'matrix'))
  expect_true(dim(sf_boot[["nestedness"]])[2] == 6) # because n_boot = 6

  # Are all final bootstraped network values equal? That is, for each column,
  # last row values must stay the same. For the code idea, see the comment of
  # https://stackoverflow.com/a/35899767/5193830
  unq_length <- get_length_unq(sf_boot[["nestedness"]])
  expect_equal(unq_length, 1)

  # Is the last bootstrap value correct? As per bipartite::networklevel ?
  final_boot_val <- sf_boot[["nestedness"]] %>% .[nrow(.), 1]
  expect_equivalent(final_boot_val, test_both)
})


# ~ higher level ----------------------------------------------------------

test_that('"nestedness" index output is ok for higher level', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_networklevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "higher",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  # Is there just a single matrix in the output list?
  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 1)
  expect_match(names(sf_boot), "nestedness", fixed = TRUE)
  expect_true(inherits(sf_boot[["nestedness"]], 'matrix'))
  expect_true(dim(sf_boot[["nestedness"]])[2] == 6) # because n_boot = 6

  # Are all final bootstraped network values equal? That is, for each column,
  # last row values must stay the same. For the code idea, see the comment of
  # https://stackoverflow.com/a/35899767/5193830
  unq_length <- get_length_unq(sf_boot[["nestedness"]])
  expect_equal(unq_length, 1)

  # Is the last bootstrap value correct? As per bipartite::networklevel ?
  final_boot_val <- sf_boot[["nestedness"]] %>% .[nrow(.), 1]
  expect_equivalent(final_boot_val, test_higher)
})


# ~ lower level -----------------------------------------------------------

test_that('"nestedness" index output is ok for lower level', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_networklevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "lower",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  # Is there just a single matrix in the output list?
  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 1)
  expect_match(names(sf_boot), "nestedness", fixed = TRUE)
  expect_true(inherits(sf_boot[["nestedness"]], 'matrix'))
  expect_true(dim(sf_boot[["nestedness"]])[2] == 6) # because n_boot = 6

  # Are all final bootstraped network values equal? That is, for each column,
  # last row values must stay the same. For the code idea, see the comment of
  # https://stackoverflow.com/a/35899767/5193830
  unq_length <- get_length_unq(sf_boot[["nestedness"]])
  expect_equal(unq_length, 1)

  # Is the last bootstrap value correct? As per bipartite::networklevel ?
  final_boot_val <- sf_boot[["nestedness"]] %>% .[nrow(.), 1]
  expect_equivalent(final_boot_val, test_lower)
})


# Test niche overlap ------------------------------------------------------

# Note that, "niche overlap" is affected by `level`.

index <- "niche overlap"

sf_sort <- Safariland %>%
  .[, order(colnames(.))] %>%
  .[order(rownames(.)), ]

set.seed(42)
test_both <- networklevel(sf_sort, index = index, level = "both")
set.seed(42)
test_higher <- networklevel(sf_sort, index = index, level = "higher")
set.seed(42)
test_lower <- networklevel(sf_sort, index = index, level = "lower")


# ~ both levels -----------------------------------------------------------

test_that('"niche overlap" index output is ok for both levels', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_networklevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "both",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  # Are there two matrices in the output list?
  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 2)

  expect_match(names(sf_boot)[1], "niche.overlap.HL", fixed = TRUE)
  expect_match(names(sf_boot)[2], "niche.overlap.LL", fixed = TRUE)

  expect_true(inherits(sf_boot[["niche.overlap.HL"]], 'matrix'))
  expect_true(inherits(sf_boot[["niche.overlap.LL"]], 'matrix'))

  expect_true(dim(sf_boot[["niche.overlap.HL"]])[2] == 6) # because n_boot = 6
  expect_true(dim(sf_boot[["niche.overlap.LL"]])[2] == 6) # because n_boot = 6

  # Are all final bootstraped network values equal? That is, for each column,
  # last row values must stay the same. For the code idea, see the comment of
  # https://stackoverflow.com/a/35899767/5193830
  unq_length_hl <- get_length_unq(sf_boot[["niche.overlap.HL"]])
  expect_equal(unq_length_hl, 1)

  unq_length_ll <- get_length_unq(sf_boot[["niche.overlap.LL"]])
  expect_equal(unq_length_ll, 1)

  # Is the last bootstrap value correct? As per bipartite::networklevel ?
  final_boot_val_hl <- sf_boot[["niche.overlap.HL"]] %>% .[nrow(.), 1]
  expect_equivalent(final_boot_val_hl, test_both["niche.overlap.HL"])

  final_boot_val_ll <- sf_boot[["niche.overlap.LL"]] %>% .[nrow(.), 1]
  expect_equivalent(final_boot_val_ll, test_both["niche.overlap.LL"])
})


# ~ higher level ----------------------------------------------------------

test_that('"niche overlap" index output is ok for higher level', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_networklevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "higher",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  # Is there just a single matrix in the output list?
  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 1)
  expect_match(names(sf_boot), "niche.overlap.HL", fixed = TRUE)
  expect_true(inherits(sf_boot[["niche.overlap.HL"]], 'matrix'))
  expect_true(dim(sf_boot[["niche.overlap.HL"]])[2] == 6) # because n_boot = 6

  # Are all final bootstraped network values equal? That is, for each column,
  # last row values must stay the same. For the code idea, see the comment of
  # https://stackoverflow.com/a/35899767/5193830
  unq_length <- get_length_unq(sf_boot[["niche.overlap.HL"]])
  expect_equal(unq_length, 1)

  # Is the last bootstrap value correct? As per bipartite::networklevel ?
  final_boot_val <- sf_boot[["niche.overlap.HL"]] %>% .[nrow(.), 1]
  expect_equivalent(final_boot_val, test_higher)
})


# ~ lower level -----------------------------------------------------------

test_that('"niche overlap" index output is ok for lower level', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_networklevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "lower",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  # Is there just a single matrix in the output list?
  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 1)
  expect_match(names(sf_boot), "niche.overlap.LL", fixed = TRUE)
  expect_true(inherits(sf_boot[["niche.overlap.LL"]], 'matrix'))
  expect_true(dim(sf_boot[["niche.overlap.LL"]])[2] == 6) # because n_boot = 6

  # Are all final bootstraped network values equal? That is, for each column,
  # last row values must stay the same. For the code idea, see the comment of
  # https://stackoverflow.com/a/35899767/5193830
  unq_length <- get_length_unq(sf_boot[["niche.overlap.LL"]])
  expect_equal(unq_length, 1)

  # Is the last bootstrap value correct? As per bipartite::networklevel ?
  final_boot_val <- sf_boot[["niche.overlap.LL"]] %>% .[nrow(.), 1]
  expect_equivalent(final_boot_val, test_lower)
})


# Test generality ---------------------------------------------------------

# Note that "generality" is affected by `level`, but outputs as results
# "generality" for higher level and "vulnerability" for lower level.

index <- "generality"

sf_sort <- Safariland %>%
  .[, order(colnames(.))] %>%
  .[order(rownames(.)), ]

set.seed(42)
test_both <- networklevel(sf_sort, index = index, level = "both")
set.seed(42)
test_higher <- networklevel(sf_sort, index = index, level = "higher")
set.seed(42)
test_lower <- expect_warning(networklevel(sf_sort, index = index, level = "lower"),
                             "You requested 'generality' for the lower level.*")
# Expected warning: You requested 'generality' for the lower level, although it
# is not a lower level index! You will get 'vulnerability' instead (same thing,
# really).

# ~ both levels -----------------------------------------------------------

test_that('"generality" index output is ok for both levels', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_networklevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "both",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  # Are there two matrices in the output list?
  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 2)

  expect_match(names(sf_boot)[1], "generality.HL", fixed = TRUE)
  expect_match(names(sf_boot)[2], "vulnerability.LL", fixed = TRUE)

  expect_true(inherits(sf_boot[["generality.HL"]], 'matrix'))
  expect_true(inherits(sf_boot[["vulnerability.LL"]], 'matrix'))

  expect_true(dim(sf_boot[["generality.HL"]])[2] == 6) # because n_boot = 6
  expect_true(dim(sf_boot[["vulnerability.LL"]])[2] == 6) # because n_boot = 6

  # Are all final bootstraped network values equal? That is, for each column,
  # last row values must stay the same. For the code idea, see the comment of
  # https://stackoverflow.com/a/35899767/5193830
  unq_length_hl <- get_length_unq(sf_boot[["generality.HL"]])
  expect_equal(unq_length_hl, 1)

  unq_length_ll <- get_length_unq(sf_boot[["vulnerability.LL"]])
  expect_equal(unq_length_ll, 1)

  # Is the last bootstrap value correct? As per bipartite::networklevel ?
  final_boot_val_hl <- sf_boot[["generality.HL"]] %>% .[nrow(.), 1]
  expect_equivalent(final_boot_val_hl, test_both["generality.HL"])

  final_boot_val_ll <- sf_boot[["vulnerability.LL"]] %>% .[nrow(.), 1]
  expect_equivalent(final_boot_val_ll, test_both["vulnerability.LL"])
})


# ~ higher level ----------------------------------------------------------

test_that('"generality" index output is ok for higher level', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_networklevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "higher",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  # Is there just a single matrix in the output list?
  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 1)
  expect_match(names(sf_boot), "generality.HL", fixed = TRUE)
  expect_true(inherits(sf_boot[["generality.HL"]], 'matrix'))
  expect_true(dim(sf_boot[["generality.HL"]])[2] == 6) # because n_boot = 6

  # Are all final bootstraped network values equal? That is, for each column,
  # last row values must stay the same. For the code idea, see the comment of
  # https://stackoverflow.com/a/35899767/5193830
  unq_length <- get_length_unq(sf_boot[["generality.HL"]])
  expect_equal(unq_length, 1)

  # Is the last bootstrap value correct? As per bipartite::networklevel ?
  final_boot_val <- sf_boot[["generality.HL"]] %>% .[nrow(.), 1]
  expect_equivalent(final_boot_val, test_higher)
})


# ~ lower level -----------------------------------------------------------

test_that('"generality" index output is ok for lower level', {
  sf_boot <- Safariland %>%
    web_matrix_to_df() %>%
    boot_networklevel_n(col_lower = "lower",
                        col_higher = "higher",
                        index = index,
                        level = "lower",
                        start = 100,
                        step = 100,
                        n_boot = 6,
                        n_cpu = 2)

  # Is there just a single matrix in the output list?
  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 1)
  expect_match(names(sf_boot), "vulnerability.LL", fixed = TRUE)
  expect_true(inherits(sf_boot[["vulnerability.LL"]], 'matrix'))
  expect_true(dim(sf_boot[["vulnerability.LL"]])[2] == 6) # because n_boot = 6

  # Are all final bootstraped network values equal? That is, for each column,
  # last row values must stay the same. For the code idea, see the comment of
  # https://stackoverflow.com/a/35899767/5193830
  unq_length <- get_length_unq(sf_boot[["vulnerability.LL"]])
  expect_equal(unq_length, 1)

  # Is the last bootstrap value correct? As per bipartite::networklevel ?
  final_boot_val <- sf_boot[["vulnerability.LL"]] %>% .[nrow(.), 1]
  expect_equivalent(final_boot_val, test_lower)
})
