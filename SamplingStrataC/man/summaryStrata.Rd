\name{summaryStrata}
\alias{summaryStrata}
\title{Information on strata structure}
\description{
  Information on strata (population, allocation, sampling rate 
  and X variables ranges)
}
\usage{
summaryStrata(x,outstrata,progress,writeFiles)
}
\arguments{
  \item{x}{the sampling frame}
  \item{outstrata}{the optimized strata}
  \item{progress}{progress bar}
  \item{writeFiles}{csv output of the strata structure}
 }
\value{
A formatted output containing information on the strata in the given domain
}
\examples{
\dontrun{
library(SamplingStrata)
data("swissmunicipalities")
data("errors")
errors$CV1 <- 0.1
errors$CV2 <- 0.1
errors <- errors[rep(row.names(errors),7),]
errors$domainvalue <- c(1:7)
errors
swissmunicipalities$id <- c(1:nrow(swissmunicipalities))
swissmunicipalities$domain = 1
frame <- buildFrameDF(swissmunicipalities,
                      id = "id",
                      domainvalue = "REG",
                      X = c("Surfacesbois","Surfacescult"),
                      Y = c("Pop020", "Pop2040")
)
solution <- optimizeStrata2 (
  errors, 
  frame,
  nStrata = 5,
  iter = 10, 
  pops = 10, 
  writeFiles = FALSE, 
  showPlot = TRUE, 
  parallel = FALSE) 
strataStructure <- summaryStrata(solution$framenew, solution$aggr_strata)
strataStructure
}
}
\keyword{survey}
