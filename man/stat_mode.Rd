% Generated by roxygen2 (4.0.2): do not edit by hand
\name{stat_mode}
\alias{stat_mode}
\title{Find the mode from a set of samples}
\usage{
stat_mode(x)
}
\arguments{
\item{x}{vector}
}
\value{
a numeric
}
\description{
This function finds the mode from a set of samples on a discrete space, i.e. it returns the value which is sampled the most often.
}
\examples{
x <- rpois(lambda = 10,n=1000)
stat_mode(x)
}

