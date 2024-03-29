\name{binomCI}
\alias{binomCI}
\title{ Confidence Intervals for Binomial Proportions }
\description{
  This function can be used to compute confidence intervals for binomial proportions.
}
\usage{
binomCI(x, n, conf.level = 0.95, method = "wilson", rand = 123, 
        R = 9999, bootci.type = "all", 
        alternative = c("two.sided", "less", "greater"), ...)
}
\arguments{
  \item{x}{ number of successes }
  \item{n}{ number of trials }
  \item{conf.level}{ confidence level }
  \item{method}{ character string specifing which method to use; see details. }
  \item{rand}{ seed for random number generator; see details. }
  \item{R}{ number of bootstrap replicates. }
  \item{bootci.type}{ type of bootstrap interval; see \code{\link[boot]{boot.ci}}. }
  \item{alternative}{ a character string specifying one- or two-sided confidence intervals. 
  Must be one of "two.sided" (default), "greater" or "less" (one-sided intervals). 
  You can specify just the initial letter. }
  \item{\dots}{ further arguments passed to function \code{\link[boot]{boot}},
  e.g. for parallel computing. }
}
\details{
  The Wald interval is obtained by inverting the acceptance region of the Wald
  large-sample normal test. There is also a Wald interval with continuity 
  correction (\code{"wald-cc"}).

  The Wilson interval, which is the default, was introduced by Wilson (1927) and is
  the inversion of the CLT approximation to the family of equal tail tests of p = p0.
  The Wilson interval is recommended by Agresti and Coull (1998) as well as by
  Brown et al (2001).

  The Agresti-Coull interval was proposed by Agresti and Coull (1998) and is a slight
  modification of the Wilson interval. The Agresti-Coull intervals are never shorter
  than the Wilson intervals; cf. Brown et al (2001).

  The Jeffreys interval is an implementation of the equal-tailed Jeffreys prior interval
  as given in Brown et al (2001).

  The modified Wilson interval is a modification of the Wilson interval for x close to 0
  or n as proposed by Brown et al (2001).

  The modified Jeffreys interval is a modification of the Jeffreys interval for
  \code{x == 0 | x == 1} and \code{x == n-1 | x == n} as proposed by
  Brown et al (2001).

  The Clopper-Pearson interval is based on quantiles of corresponding beta
  distributions. This is sometimes also called exact interval.

  The arcsine interval is based on the variance stabilizing distribution for the binomial
  distribution.

  The logit interval is obtained by inverting the Wald type interval for the log odds.

  The Witting interval (cf. Beispiel 2.106 in Witting (1985)) uses randomization to
  obtain uniformly optimal lower and upper confidence bounds (cf. Satz 2.105 in
  Witting (1985)) for binomial proportions.
  
  The bootstrap interval is calculated by using function \code{\link[boot]{boot.ci}}.

  For more details we refer to Brown et al (2001) as well as Witting (1985).
}
\value{
  A list with class \code{"confint"} containing the following components:
  \item{estimate}{ the estimated probability of success. }
  \item{conf.int}{ a confidence interval for the probability of success. }
}
\references{
  A. Agresti and B.A. Coull (1998). Approximate is better than "exact" for interval
  estimation of binomial proportions.
  \emph{American Statistician}, \bold{52}, 119-126.

  L.D. Brown, T.T. Cai and A. Dasgupta (2001). Interval estimation for a binomial
  proportion. \emph{Statistical Science}, \bold{16}(2), 101-133.

  H. Witting (1985). \emph{Mathematische Statistik I}. Stuttgart: Teubner.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\note{
A first version of this function appeared in R package SLmisc.
}
\seealso{\code{\link[stats]{binom.test}}}
\examples{
binomCI(x = 42, n = 43, method = "wald")
binomCI(x = 42, n = 43, method = "wald-cc")
binomCI(x = 42, n = 43, method = "wilson")
binomCI(x = 42, n = 43, method = "agresti-coull")
binomCI(x = 42, n = 43, method = "jeffreys")
binomCI(x = 42, n = 43, method = "modified wilson")
binomCI(x = 42, n = 43, method = "modified jeffreys")
binomCI(x = 42, n = 43, method = "clopper-pearson")
binomCI(x = 42, n = 43, method = "arcsine")
binomCI(x = 42, n = 43, method = "logit")
binomCI(x = 42, n = 43, method = "witting")
## bootstrap intervals (R = 999 to reduce computation time for R checks)
binomCI(x = 42, n = 43, method = "boot", R = 999) # may generate values > 1!

## the confidence interval computed by binom.test
## corresponds to the Clopper-Pearson interval
binomCI(x = 42, n = 43, method = "clopper-pearson")$conf.int
binom.test(x = 42, n = 43)$conf.int

## one-sided intervals
binomCI(x = 10, n = 43, alternative = "less")
binomCI(x = 10, n = 43, alternative = "less", method = "boot", 
        bootci.type = "bca", R = 999)
binomCI(x = 10, n = 43, alternative = "greater", method = "boot", 
        bootci.type = "perc", R = 999)
        
\donttest{
## parallel computing for bootstrap
binomCI(x = 10, n = 43, method = "boot", R = 9999, 
        parallel = "multicore", ncpus = 2)
}
}
\keyword{univar}
