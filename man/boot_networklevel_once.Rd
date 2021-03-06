% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_networklevel.r
\name{boot_networklevel_once}
\alias{boot_networklevel_once}
\title{Bootstrap network level metric once.}
\usage{
boot_networklevel_once(
  data,
  col_lower,
  col_higher,
  index,
  level,
  start,
  step,
  seed,
  ...
)
}
\arguments{
\item{data}{See \code{\link[bootstrapnet]{boot_networklevel_n}}.}

\item{col_lower}{See \code{\link[bootstrapnet]{boot_networklevel_n}}.}

\item{col_higher}{See \code{\link[bootstrapnet]{boot_networklevel_n}}.}

\item{index}{See \code{\link[bootstrapnet]{boot_networklevel_n}}.}

\item{level}{See \code{\link[bootstrapnet]{boot_networklevel_n}}.}

\item{start}{See \code{\link[bootstrapnet]{boot_networklevel_n}}.}

\item{step}{See \code{\link[bootstrapnet]{boot_networklevel_n}}.}

\item{seed}{Set seed to get reproducible random results. Passed to \code{set.seed()}.}

\item{...}{See \code{\link[bootstrapnet]{boot_networklevel_n}}.}
}
\value{
Returns a data frame of one or two columns depending on the provided \code{index}
metric and \code{level} value. For example if \code{index='nestedness'} and
\code{level='both'}, then it returns a one column data frame. But in case of
\code{index='niche overlap'}, it returns a two columns data frame, first column
('niche.overlap.HL') corresponding to the higher level (e.g. insects), and
the second column ('niche.overlap.LL') for the lower level (e.g. plants).
The number of rows in the returned data frame corresponds to the number of
splits given by \code{bootstrapnet:::sample_indices} (is the length/number of
elements of the list returned by this function). The last value(s) in the
data frame correspond to the results of
\code{\link[bipartite]{networklevel}} on the entire network - see in
examples.
}
\description{
You will rarely use this function alone. It was designed to be executed in
parallel by \code{\link[bootstrapnet]{boot_networklevel_n}}. See more
details there.
}
\examples{

library(bootstrapnet)
library(bipartite)
library(data.table)

data(Safariland)

df <- web_matrix_to_df(Safariland)
setDT(df) # df must be a data.table and not data.frame

# Example with "nestedness":

boot_networklevel_once(data = df,
                       col_lower = "lower", # column name for plants
                       col_higher = "higher", # column name for insects
                       index = "nestedness",
                       level = "both",
                       start = 100,
                       step = 100,
                       seed = 2020-11-5)
# Returns a one column data frame. The last value must equal the result of:
set.seed(42) # this is the same seed used in the boot_networklevel_once() function
bipartite::networklevel(table(df$lower, df$higher),
                       index = "nestedness")
# which is also the equivalent of:
Safariland_sorted <- Safariland[order(rownames(Safariland)), order(colnames(Safariland))]
set.seed(42)
bipartite::networklevel(Safariland_sorted, index = "nestedness")


# Example with "niche overlap":

boot_networklevel_once(data = df,
                       col_lower = "lower", # column name for plants
                       col_higher = "higher", # column name for insects
                       index = "niche overlap",
                       level = "both",
                       start = 100,
                       step = 100,
                       seed = 2020-11-5)
# Returns a two column data frame

}
\references{
This function is a wrapper of \code{\link[bipartite]{networklevel}}.
}
