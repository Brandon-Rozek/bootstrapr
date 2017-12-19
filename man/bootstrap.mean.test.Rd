\name{bootstrap.mean.test}
\alias{bootstrap.mean.test}
\title{Bootstrap Test for Means}
\description{Performs one and two sample bootstrap tests for means on vectors of data}
\usage{
bootstrap.mean.test(x, ...)

## Default method:
bootstrap.mean.test(x, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, conf.level = 0.95, bootstrap_samples = 1e+06)
}
\arguments{
  \item{x}{a (non-empty) numeric vector of data values.}
  \item{y}{an optional (non-empty) numeric vector of data values.}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.}
  \item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test).}
  \item{conf.level}{confidence level of the interval.}
  \item{bootstrap_samples}{the number of bootstrap replicants to generate. Influences the accuracy of the p-value and confidence interval}
}
\details{
alternative = "greater" is the alternative that x has a larger mean than y.

If the input data are effectively constant (compared to the larger of the two means) an error is generated.
}
\value{
A list with class \code{"htest"} containing the following components:
  \item{p.value}{the p-value for the test.}
  \item{conf.int}{a confidence interval for the mean appropriate to the
    specified alternative hypothesis.}
  \item{estimate}{the estimated mean or difference in means depending on
    whether it was a one-sample test or a two-sample test.}
  \item{null.value}{the specified hypothesized value of the mean or mean
    difference depending on whether it was a one-sample test or a
    two-sample test.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a character string indicating what type of t-test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{
Gould, Rob. Bootstrap Hypothesis Test. PDF.
http://www.stat.ucla.edu/~rgould/110as02/bshypothesis.pdf

Myung, Jay. Bootstrap Hypothesis Testing. PDF.
http://faculty.psy.ohio-state.edu/myung/personal/course/826/bootstrap_hypo.pdf

}
\author{
Brandon Rozek (brozek@mail.umw.edu)
}

\seealso{
  \code{\link{t.test}}
}
\examples{
x = rnorm(200, 2, 5)
y = rnorm(200, 2, 5)
# One sample test
bootstrap.mean.test(x, mu = 1)

# Two sample test
bootstrap.mean.test(x, y)
}

\keyword{htest}
