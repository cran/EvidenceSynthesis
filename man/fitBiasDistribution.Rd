% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesianBiasCorrection.R
\name{fitBiasDistribution}
\alias{fitBiasDistribution}
\title{Fit Bias Distribution}
\usage{
fitBiasDistribution(
  likelihoodProfiles,
  priorSds = c(2, 0.5),
  numsamps = 10000,
  thin = 10,
  minNCs = 5,
  robust = FALSE,
  df = 4,
  seed = 1
)
}
\arguments{
\item{likelihoodProfiles}{A list of grid profile likelihoods regarding negative controls.}

\item{priorSds}{A two-dimensional vector with the standard deviation of the prior for
the average bias and the sd/scale parameter, respectively.}

\item{numsamps}{Total number of MCMC samples needed.}

\item{thin}{Thinning frequency: how many iterations before another sample is obtained?}

\item{minNCs}{Minimum number of negative controls needed to fit a bias distribution;
default (also recommended): 5.}

\item{robust}{Whether or not to use a t-distribution model; default: FALSE.}

\item{df}{Degrees of freedom for the t-model, only used if robust is TRUE.}

\item{seed}{Seed for the random number generator.}
}
\value{
A dataframe with three columns and \code{numsamps} number of rows.
Column \code{mean} includes MCMC samples for the average bias,
\code{scale} for the sd/scale parameter,
and \code{bias} for predictive samples of the bias.
}
\description{
Learn an empirical distribution on estimation bias by simultaneously analyzing
a large set of negative control outcomes by a Bayesian hierarchical model through MCMC.
Analysis is based on a list of extracted likelihood profiles.
}
\examples{
# load example data
data("ncLikelihoods")

# fit a bias distributions by analyzing a set of negative control outcomes
# for example, for the 5th analysis period, and using the t model
# NOT RUN
# biasDistribution = fitBiasDistribution(ncLikelihoods[[5]], robust = TRUE)

}
\seealso{
\link{computeBayesianMetaAnalysis}
}
