# Declare global variables and native symbol objects ----------------------

# Doing so, avoids the note from devtools::check():
# "no visible binding for global variable".
# See https://stackoverflow.com/a/12429344/5193830
# or https://stackoverflow.com/a/17807914/5193830

.onLoad <- function(...) {
  if (getRversion() >= "2.15.1")
    utils::globalVariables(
      c(
        '.',
        ':=',
        '%>%',
        '%dopar%',
        'ggplot',
        'geom_line',
        'aes',
        'type',
        'simulation_id',
        'ci_low',
        'ci_up',
        'n_inter',
        'value',
        'sp',
        'quantile',
        'counts',
        'spl_size',
        'web'
      )
    )
}


# Package startup message -------------------------------------------------

.onAttach <- function(...) {
  packageStartupMessage(strwrap("Bootstrap network metrics", indent = 5))
}
