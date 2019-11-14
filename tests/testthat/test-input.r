context("Wrong inputs should trigger errors")

library(dplyr)
library(bipartite)

data(Safariland)

# Generate a fictive smaller network.
set.seed(321)
saf_df <- Safariland[, sort(sample.int(ncol(Safariland), 10))] %>%
  web_matrix_to_df()


# Too big start and step arguments ----------------------------------------

test_that('error is triggered when the arguments start and step are bigger than the number of interactions', {
  # boot_networklevel - wrong start
  expect_error(
    boot_networklevel(lst = list(s1 = saf_df),
                      col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "nestedness",
                      level = "both",
                      start = 50,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    regexp = "Possibly because your `start` value"
  )

  # boot_networklevel - wrong step
  expect_error(
    boot_networklevel(lst = list(s1 = saf_df),
                      col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "nestedness",
                      level = "both",
                      start = 10,
                      step = 50,
                      n_boot = 6,
                      n_cpu = 2),
    regexp = "Possibly because your `step` value"
  )

  # boot_specieslevel - wrong start
  expect_error(
    boot_specieslevel(lst = list(s1 = saf_df),
                      col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "betweenness",
                      level = "both",
                      start = 50,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    regexp = "Possibly because your `start` value"
  )

  # boot_specieslevel - wrong step
  expect_error(
    boot_specieslevel(lst = list(s1 = saf_df),
                      col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "betweenness",
                      level = "both",
                      start = 10,
                      step = 50,
                      n_boot = 6,
                      n_cpu = 2),
    regexp = "Possibly because your `step` value"
  )
})


# Test species names ------------------------------------------------------

test_that('error is triggered when species names are NA or empty', {
  # boot_networklevel_n - higher level
  expect_error(
    boot_networklevel_n(data = saf_df %>% mutate(higher = NA),
                        col_lower = "lower", # column name for plants
                        col_higher = "higher", # column name for insects
                        index = "nestedness",
                        level = "both",
                        start = 10,
                        step = 10,
                        n_boot = 6,
                        n_cpu = 2),
    "You have undefined/empty species names. Check the higher level species names for NA-s or empty strings."
  )
  # boot_networklevel_n - lower level
  expect_error(
    boot_networklevel_n(data = saf_df %>% mutate(lower = NA),
                        col_lower = "lower", # column name for plants
                        col_higher = "higher", # column name for insects
                        index = "nestedness",
                        level = "both",
                        start = 10,
                        step = 10,
                        n_boot = 6,
                        n_cpu = 2),
    "You have undefined/empty species names. Check the lower level species names for NA-s or empty strings."
  )

  # boot_specieslevel_n - higher level
  expect_error(
    boot_specieslevel_n(data = saf_df %>% mutate(higher = NA),
                        col_lower = "lower", # column name for plants
                        col_higher = "higher", # column name for insects
                        index = "betweenness",
                        level = "both",
                        start = 10,
                        step = 10,
                        n_boot = 6,
                        n_cpu = 2),
    "You have undefined/empty species names. Check the higher level species names for NA-s or empty strings."
  )
  # boot_specieslevel_n - lower level
  expect_error(
    boot_specieslevel_n(data = saf_df %>% mutate(lower = NA),
                        col_lower = "lower", # column name for plants
                        col_higher = "higher", # column name for insects
                        index = "betweenness",
                        level = "both",
                        start = 10,
                        step = 10,
                        n_boot = 6,
                        n_cpu = 2),
    "You have undefined/empty species names. Check the lower level species names for NA-s or empty strings."
  )
})


# Test typo in level ------------------------------------------------------

test_that('error is triggered when typo in `level`', {
  # for networklevel
  expect_error(
    boot_networklevel(lst = list(s1 = saf_df),
                      col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "nestedness",
                      level = "nice_typo",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    'The value provided for `level` must be one of the following: "both", "lower" or "higher"'
  )

  # for specieslevel
  expect_error(
    boot_specieslevel(lst = list(s1 = saf_df),
                      col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "betweenness",
                      level = "nice_typo",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    'The value provided for `level` must be one of the following: "both", "lower" or "higher"'
  )
})


# Test when more than one value in `level` --------------------------------

test_that('error is triggered when more than one value in `level`', {
  # for networklevel
  expect_error(
    boot_networklevel(lst = list(s1 = saf_df),
                      col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "nestedness",
                      level = c("higher", "lower"),
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    'Provide a single value for `level`: "both", "lower" or "higher"'
  )

  # for specieslevel
  expect_error(
    boot_specieslevel(lst = list(s1 = saf_df),
                      col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "betweenness",
                      level = c("higher", "lower"),
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    'Provide a single value for `level`: "both", "lower" or "higher"'
  )
})


# Test typo in index ------------------------------------------------------

test_that('error is triggered when typo in index', {
  # for networklevel
  expect_error(
    boot_networklevel(lst = list(s1 = saf_df),
                      col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "all", # or some typo
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "Your index is not recognised.*"
  )

  # for specieslevel
  expect_error(
    boot_specieslevel(lst = list(s1 = saf_df),
                      col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "all", # or some typo
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "Your index is not recognised.*"
  )
})


# Test data ---------------------------------------------------------------

test_that('error is triggered when `data` is of different class', {
  # for networklevel
  expect_error(
    boot_networklevel(lst = list(s1 = Safariland), # Safariland is a matrix
                      col_lower = "lower",   # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "nestedness",
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "`data` must be a data.frame or data.table.*"
  )

  # for specieslevel
  expect_error(
    boot_specieslevel(lst = list(s1 = Safariland), # Safariland is a matrix
                      col_lower = "lower",   # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "betweenness",
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "`data` must be a data.frame or data.table.*"
  )
})


test_that('error is triggered when `data` is empty', {
  empty_data <- saf_df %>% filter(lower == "bla-bla")

  # for networklevel
  expect_error(
    boot_networklevel(lst = list(s1 = empty_data),
                      col_lower = "lower",   # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "nestedness",
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "Your `data` has no rows.*"
  )

  # for specieslevel
  expect_error(
    boot_specieslevel(lst = list(s1 = empty_data),
                      col_lower = "lower",   # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "betweenness",
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "Your `data` has no rows.*"
  )
})


test_that('error is triggered when `data` has less than 3 columns', {
  wrong_data <- saf_df %>% select(- counts)

  # for networklevel
  expect_error(
    boot_networklevel(lst = list(s1 = wrong_data),
                      col_lower = "lower",   # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "nestedness",
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "`data` must have at least 3 columns.*"
  )

  # for specieslevel
  expect_error(
    boot_specieslevel(lst = list(s1 = wrong_data),
                      col_lower = "lower",   # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "betweenness",
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "`data` must have at least 3 columns.*"
  )
})


test_that('error is triggered when typo in col_lower or col_higher', {
  # for networklevel
  expect_error(
    boot_networklevel(lst = list(s1 = saf_df),
                      col_lower = "typo",   # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "nestedness",
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "Cannot find the column name provided in `col_lower`.*"
  )

  expect_error(
    boot_networklevel(lst = list(s1 = saf_df),
                      col_lower = "lower",   # column name for plants
                      col_higher = "typo", # column name for insects
                      index = "nestedness",
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "Cannot find the column name provided in `col_higher`.*"
  )

  # for specieslevel
  expect_error(
    boot_specieslevel(lst = list(s1 = saf_df),
                      col_lower = "typo",   # column name for plants
                      col_higher = "higher", # column name for insects
                      index = "betweenness",
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "Cannot find the column name provided in `col_lower`.*"
  )

  expect_error(
    boot_specieslevel(lst = list(s1 = saf_df),
                      col_lower = "lower",   # column name for plants
                      col_higher = "typo", # column name for insects
                      index = "nestedness",
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2),
    "Cannot find the column name provided in `col_higher`.*"
  )
})


test_that('warning is triggered when interactions are unique', {
  doubtful_data <- web_matrix_to_df(Safariland) %>% unique()

  # for networklevel
  expect_warning(
    out <- boot_networklevel(lst = list(s1 = doubtful_data),
                             col_lower = "lower",   # column name for plants
                             col_higher = "higher", # column name for insects
                             index = "nestedness",
                             level = "both",
                             start = 10,
                             step = 10,
                             n_boot = 4,
                             n_cpu = 2),
    "Each interaction .* should be repeated.*"
  )

  # for specieslevel
  expect_warning(
    out <- boot_specieslevel(lst = list(s1 = doubtful_data),
                             col_lower = "lower",   # column name for plants
                             col_higher = "higher", # column name for insects
                             index = "betweenness",
                             level = "both",
                             start = 10,
                             step = 10,
                             n_boot = 4,
                             n_cpu = 2),
    "Each interaction .* should be repeated.*"
  )
})
