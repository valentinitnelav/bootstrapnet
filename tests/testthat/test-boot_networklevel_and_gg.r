context("Network level indices - test boot_networklevel & gg_networklevel")

library(dplyr)
library(bipartite)
data(Safariland)

# Helper function to get the length of unique values of the last bootstrap.
get_length_unq <- function(df) {
  df %>%
    select(mean, ci_low, ci_up) %>%
    as.numeric() %>%
    round(digits = 3) %>%
    unique() %>%
    length()
}

# Generate two fictive networks to compare. For full reproducibly is safer to
# reorder the Safariland matrices alphabetically.
set.seed(321)
Safariland_1 <- Safariland[, sort(sample.int(ncol(Safariland), 10))] %>%
  .[, order(colnames(.))] %>% .[order(rownames(.)), ]
set.seed(123)
Safariland_2 <- Safariland[, sort(sample.int(ncol(Safariland), 10))] %>%
  .[, order(colnames(.))] %>% .[order(rownames(.)), ]


# Test nestedness ---------------------------------------------------------

# Note that, "nestedness" is not be affected by `level` as opposed to "niche
# overlap".

index <- "nestedness"

set.seed(42)
test_both_s1 <- networklevel(Safariland_1, index = index, level = "both")
set.seed(42)
test_both_s2 <- networklevel(Safariland_2, index = index, level = "both")


test_that('"nestedness" index output is ok for both levels', {
  sf_boot <- list(s1 = Safariland_1, s2 = Safariland_2) %>%
    lapply(web_matrix_to_df) %>%
    boot_networklevel(col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = index,
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2)

  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 1)
  expect_match(names(sf_boot), "nestedness", fixed = TRUE)

  expect_true(inherits(sf_boot[["nestedness"]], 'list'))
  expect_true(length(sf_boot[["nestedness"]]) == 2)

  expect_match(names(sf_boot[["nestedness"]])[1], "stats_df", fixed = TRUE)
  expect_match(names(sf_boot[["nestedness"]])[2], "lines_df", fixed = TRUE)

  expect_true(inherits(sf_boot[["nestedness"]][["stats_df"]], 'data.frame'))
  expect_true(inherits(sf_boot[["nestedness"]][["lines_df"]], 'data.frame'))

  expect_true(ncol(sf_boot[["nestedness"]][["stats_df"]]) == 5)
  expect_true(ncol(sf_boot[["nestedness"]][["lines_df"]]) == 4)

  # Get final boot values from the stats_df data frame
  last_boot_s1 <- sf_boot[["nestedness"]][["stats_df"]] %>%
    filter(web == "s1") %>%
    filter(row_number() == n())

  last_boot_s2 <- sf_boot[["nestedness"]][["stats_df"]] %>%
    filter(web == "s2") %>%
    filter(row_number() == n())

  # Are all final boostrap values converging (are all equal) for each web?
  expect_equal(get_length_unq(last_boot_s1), 1)
  expect_equal(get_length_unq(last_boot_s2), 1)

  # Is the last bootstrap value correct? As per bipartite::networklevel ?
  expect_equivalent(last_boot_s1 %>% select(mean),
                    test_both_s1)
  expect_equivalent(last_boot_s2 %>% select(mean),
                    test_both_s2)


  # Test graph function
  sf_gg_lst <- gg_networklevel(sf_boot)
  expect_true(inherits(sf_gg_lst, 'list'))
  expect_true(length(sf_gg_lst) == 1)
  expect_true(inherits(sf_gg_lst[["nestedness"]], c("gg", "ggplot")))
})


# Test niche overlap ------------------------------------------------------

# Note that, "niche overlap" is affected by `level`.

index <- "niche overlap"

set.seed(42)
test_both_s1 <- networklevel(Safariland_1, index = index, level = "both")
set.seed(42)
test_both_s2 <- networklevel(Safariland_2, index = index, level = "both")


test_that('"niche overlap" index output is ok for both levels', {
  sf_boot <- list(s1 = Safariland_1, s2 = Safariland_2) %>%
    lapply(web_matrix_to_df) %>%
    boot_networklevel(col_lower = "lower", # column name for plants
                      col_higher = "higher", # column name for insects
                      index = index,
                      level = "both",
                      start = 10,
                      step = 10,
                      n_boot = 6,
                      n_cpu = 2)

  expect_true(inherits(sf_boot, 'list'))
  expect_true(length(sf_boot) == 2)

  expect_match(names(sf_boot)[1], "niche.overlap.HL", fixed = TRUE)
  expect_match(names(sf_boot)[2], "niche.overlap.LL", fixed = TRUE)

  expect_true(inherits(sf_boot[["niche.overlap.HL"]], 'list'))
  expect_true(inherits(sf_boot[["niche.overlap.LL"]], 'list'))

  expect_true(length(sf_boot[["niche.overlap.HL"]]) == 2)
  expect_true(length(sf_boot[["niche.overlap.LL"]]) == 2)

  expect_match(names(sf_boot[["niche.overlap.HL"]])[1], "stats_df", fixed = TRUE)
  expect_match(names(sf_boot[["niche.overlap.HL"]])[2], "lines_df", fixed = TRUE)

  expect_match(names(sf_boot[["niche.overlap.LL"]])[1], "stats_df", fixed = TRUE)
  expect_match(names(sf_boot[["niche.overlap.LL"]])[2], "lines_df", fixed = TRUE)

  expect_true(inherits(sf_boot[["niche.overlap.HL"]][["stats_df"]], 'data.frame'))
  expect_true(inherits(sf_boot[["niche.overlap.HL"]][["lines_df"]], 'data.frame'))

  expect_true(inherits(sf_boot[["niche.overlap.LL"]][["stats_df"]], 'data.frame'))
  expect_true(inherits(sf_boot[["niche.overlap.LL"]][["lines_df"]], 'data.frame'))

  expect_true(ncol(sf_boot[["niche.overlap.HL"]][["stats_df"]]) == 5)
  expect_true(ncol(sf_boot[["niche.overlap.HL"]][["lines_df"]]) == 4)

  expect_true(ncol(sf_boot[["niche.overlap.LL"]][["stats_df"]]) == 5)
  expect_true(ncol(sf_boot[["niche.overlap.LL"]][["lines_df"]]) == 4)

  # Get final boot values from the stats_df data frame
  # - high level
  last_boot_s1_hl <- sf_boot[["niche.overlap.HL"]][["stats_df"]] %>%
    filter(web == "s1") %>%
    filter(row_number() == n())

  last_boot_s2_hl <- sf_boot[["niche.overlap.HL"]][["stats_df"]] %>%
    filter(web == "s2") %>%
    filter(row_number() == n())

  # - low level
  last_boot_s1_ll <- sf_boot[["niche.overlap.LL"]][["stats_df"]] %>%
    filter(web == "s1") %>%
    filter(row_number() == n())

  last_boot_s2_ll <- sf_boot[["niche.overlap.LL"]][["stats_df"]] %>%
    filter(web == "s2") %>%
    filter(row_number() == n())

  # Are all final boostrap values converging (are all equal) for each web?
  # - high level
  expect_equal(get_length_unq(last_boot_s1_hl), 1)
  expect_equal(get_length_unq(last_boot_s2_hl), 1)
  # - low level
  expect_equal(get_length_unq(last_boot_s1_ll), 1)
  expect_equal(get_length_unq(last_boot_s2_ll), 1)

  # Is the last bootstrap value correct? As per bipartite::networklevel ?
  # - high level
  expect_equivalent(last_boot_s1_hl %>% select(mean),
                    test_both_s1["niche.overlap.HL"])
  expect_equivalent(last_boot_s2_hl %>% select(mean),
                    test_both_s2["niche.overlap.HL"])
  # - low level
  expect_equivalent(last_boot_s1_ll %>% select(mean),
                    test_both_s1["niche.overlap.LL"])
  expect_equivalent(last_boot_s2_ll %>% select(mean),
                    test_both_s2["niche.overlap.LL"])

  # Test graph function
  sf_gg_lst <- gg_networklevel(sf_boot)
  expect_true(inherits(sf_gg_lst, 'list'))
  expect_true(length(sf_gg_lst) == 2)
  expect_true(inherits(sf_gg_lst[["niche.overlap.HL"]], c("gg", "ggplot")))
  expect_true(inherits(sf_gg_lst[["niche.overlap.LL"]], c("gg", "ggplot")))
})
