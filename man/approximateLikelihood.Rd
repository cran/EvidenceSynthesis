% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LikelihoodApproximation.R
\name{approximateLikelihood}
\alias{approximateLikelihood}
\title{Approximate a likelihood function}
\usage{
approximateLikelihood(
  cyclopsFit,
  parameter = 1,
  approximation = "custom",
  bounds = c(log(0.1), log(10))
)
}
\arguments{
\item{cyclopsFit}{A model fitted using the \code{\link[Cyclops:fitCyclopsModel]{Cyclops::fitCyclopsModel()}} function.}

\item{parameter}{The parameter in the \code{cyclopsFit} object to profile.}

\item{approximation}{The type of approximation. Valid options are \code{'normal'}, \code{'skew normal'},
\code{'custom'}, \code{'grid'}, or \code{'adaptive grid'}.}

\item{bounds}{The bounds on the effect size used to fit the approximation.}
}
\value{
A vector of parameters of the likelihood approximation.
}
\description{
Approximate the likelihood function using a parametric (normal, skew-normal, or custom parametric),
or grid approximation. The approximation does not reveal person-level information, and can
therefore be shared among data sites. When counts are low, a normal approximation might not be
appropriate.
}
\examples{
# Simulate some data for this example:
populations <- simulatePopulations()

cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
  data = populations[[1]],
  modelType = "cox"
)
cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
approximation <- approximateLikelihood(cyclopsFit, "x")
approximation

# (Estimates in this example will vary due to the random simulation)

}
\seealso{
\link{computeConfidenceInterval}, \link{computeFixedEffectMetaAnalysis}, \link{computeBayesianMetaAnalysis}
}
