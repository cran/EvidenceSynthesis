## ----eval=TRUE,echo=FALSE-----------------------------------------------------
set.seed(1)

## ----eval=TRUE, warning=FALSE, message=FALSE----------------------------------
library(EvidenceSynthesis)
simulationSettings <- createSimulationSettings(nSites = 5,
                                               n = 10000,
                                               treatedFraction = 0.75,
                                               nStrata = 5,
                                               hazardRatio = 2,
                                               randomEffectSd = 0.5)
populations <- simulatePopulations(simulationSettings)

## ----eval=TRUE----------------------------------------------------------------
library(Cyclops)
# Assume we are at site 1:
population <- populations[[1]]

cyclopsData <- createCyclopsData(Surv(time, y) ~ x + strata(stratumId), 
                                 data = population, 
                                 modelType = "cox")
cyclopsFit <- fitCyclopsModel(cyclopsData)

## ----eval=TRUE----------------------------------------------------------------
# Hazard ratio:
exp(coef(cyclopsFit))

# 95% confidence interval:
exp(confint(cyclopsFit, parm = "x")[2:3])

## ----eval=TRUE----------------------------------------------------------------
approximation <-  approximateLikelihood(cyclopsFit, parameter = "x", approximation = "custom")
approximation

## ----eval=TRUE----------------------------------------------------------------
plotLikelihoodFit(approximation = approximation,
                  cyclopsFit = cyclopsFit,
                  parameter = "x")

## ----eval=TRUE----------------------------------------------------------------
fitModelInDatabase <- function(population) {
  cyclopsData <- createCyclopsData(Surv(time, y) ~ x + strata(stratumId), 
                                   data = population, 
                                   modelType = "cox")
  cyclopsFit <- fitCyclopsModel(cyclopsData)
  approximation <-  approximateLikelihood(cyclopsFit, 
                                          parameter = "x",
                                          approximation = "custom")
  return(approximation)
}
approximations <- lapply(populations, fitModelInDatabase)
approximations

## ----eval=TRUE, message=FALSE-------------------------------------------------
# Combine the various approximations into a single data frame, one row per database:
approximationsTable <- do.call(rbind, approximations)

# Do meta-analysis
estimate <- computeFixedEffectMetaAnalysis(approximationsTable)
estimate

## ----eval=TRUE, message=FALSE-------------------------------------------------
# Combine the various approximations into a single data frame, one row per database:
approximationsTable <- do.call(rbind, approximations)

# Do meta-analysis
estimate <- computeBayesianMetaAnalysis(approximationsTable)
exp(estimate[1:3])

## ----eval=TRUE, message=FALSE-------------------------------------------------
plotPosterior(estimate)

## ----eval=TRUE, message=FALSE-------------------------------------------------
plotMcmcTrace(estimate)

## ----eval=TRUE, message=FALSE-------------------------------------------------
estimate2 <- computeBayesianMetaAnalysis(approximationsTable, priorSd = c(2, 0.1))
exp(estimate2[1:3])

## ----eval=TRUE, message=FALSE-------------------------------------------------
# Make up some data site labels:
labels <- paste("Data site", LETTERS[1:length(populations)])

plotMetaAnalysisForest(approximationsTable, labels, estimate)

