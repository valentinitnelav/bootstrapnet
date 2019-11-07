context("Wrong inputs should trigger errors")

library(dplyr)
library(bipartite)

data(Safariland)

# Generate two fictive networks to compare. For full reproducibly is safer to
# reorder the Safariland matrices alphabetically.
set.seed(321)
Safariland_1 <- Safariland[, sort(sample.int(ncol(Safariland), 10))] %>%
  .[, order(colnames(.))] %>% .[order(rownames(.)), ]
set.seed(123)
Safariland_2 <- Safariland[, sort(sample.int(ncol(Safariland), 10))] %>%
  .[, order(colnames(.))] %>% .[order(rownames(.)), ]


# Too big start and step arguments ----------------------------------------

test_that('error is triggered when the arguments start and step are bigger than the number of interactions', {
  # boot_networklevel - wrong start
  expect_error(
    boot_networklevel(lst = list(s1 = web_matrix_to_df(Safariland_1)),
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
    boot_networklevel(lst = list(s1 = web_matrix_to_df(Safariland_1)),
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
    boot_specieslevel(lst = list(s1 = web_matrix_to_df(Safariland_1)),
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
    boot_specieslevel(lst = list(s1 = web_matrix_to_df(Safariland_1)),
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
    boot_networklevel_n(data = web_matrix_to_df(Safariland_1) %>% mutate(higher = NA),
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
    boot_networklevel_n(data = web_matrix_to_df(Safariland_1) %>% mutate(lower = NA),
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
    boot_specieslevel_n(data = web_matrix_to_df(Safariland_1) %>% mutate(higher = NA),
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
    boot_specieslevel_n(data = web_matrix_to_df(Safariland_1) %>% mutate(lower = NA),
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
    boot_networklevel(lst = list(s1 = web_matrix_to_df(Safariland_1)),
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
    boot_specieslevel(lst = list(s1 = web_matrix_to_df(Safariland_1)),
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
    boot_networklevel(lst = list(s1 = web_matrix_to_df(Safariland_1)),
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
    boot_specieslevel(lst = list(s1 = web_matrix_to_df(Safariland_1)),
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
    boot_networklevel(lst = list(s1 = web_matrix_to_df(Safariland_1)),
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
    boot_specieslevel(lst = list(s1 = web_matrix_to_df(Safariland_1)),
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
