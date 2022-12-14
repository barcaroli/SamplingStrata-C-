\name{optimizeStrata}
\Rdversion{1.3}
\alias{optimizeStrata}
\title{
Best stratification of a sampling frame for multipurpose surveys
}
\description{
This function runs a set of other functions to optimise the stratification of a sampling frame 
}
\usage{
	optimizeStrata(
	errors , 
	strata , 
	cens = NULL, 
	strcens = FALSE,
	alldomains = TRUE,
	dom = NULL,	
	initialStrata = NA, 
	addStrataFactor = 0.0, 
	minnumstr = 2, 
	iter = 50, 
	pops = 20, 
	mut_chance = NA, 
	elitism_rate = 0.2,
	highvalue = 1e+08, 
	suggestions = NULL,
	realAllocation = TRUE,
	writeFiles = FALSE,
	showPlot = TRUE, 
	parallel = TRUE,
	cores
	)
}
\arguments{
  \item{errors}{
  This is the (mandatory) dataframe containing the precision levels expressed in terms of maximum expected value of the Coefficients of Variation
  related to target variables of the survey.
}
  \item{strata}{
  This is the (mandatory) dataframe containing the information related to "atomic" strata, i.e. the strata obtained by 
  the Cartesian product of all auxiliary variables X's. Information concerns the identifiability of strata
  (values of X's) and variability of Y's (for each Y, mean and standard error in strata).
}
  \item{cens}{
  This the (optional) dataframe containing the takeall strata, those strata whose units must be selected in whatever sample. It has same structure than "strata" dataframe.
}
  \item{strcens}{
  Flag (TRUE/FALSE) to indicate if takeall strata do exist or not. Default is FALSE.
}
  \item{alldomains}{
  Flag (TRUE/FALSE) to indicate if the optimization must be carried out on all domains (default is TRUE). If it is set to FALSE,
  then a value must be given to parameter 'dom'.
}
  \item{dom}{
  Indicates the domain on which the optimization must be carried. It is an integer value that has to be internal to the interval 
  (1 <--> number of domains). If 'alldomains' is set to TRUE, it is ignored.
}
  \item{initialStrata}{
  This is the initial limit on the number of strata in the different domains for each solution. Default is NA, and in this case it is set equal to the
  number of atomic strata in each domain.
}
  \item{addStrataFactor}{
  This parameter indicates the probability that at each mutation the number of strata may increase with 
  respect to the current value. Default is 0.0.
}
  \item{minnumstr}{
  Indicates the minimum number of units that must be allocated in each stratum. Default is 2.
}
  \item{iter}{
  Indicated the maximum number of iterations (= generations) of the genetic algorithm. Default is 50.
}
  \item{pops}{
  The dimension of each generations in terms of individuals. Default is 20.
}
  \item{mut_chance}{
  Mutation chance: for each new individual, the probability to change each single chromosome, i.e. one bit
  of the solution vector. High values of this parameter allow a deeper exploration of the solution space,
  but a slower convergence, while low values permit a faster convergence, but the final solution can be
  distant from the optimal one. Default is NA, in correspondence of which it is computed as 1/(vars+1) 
  where vars is the length of elements in the solution.
}
  \item{elitism_rate}{
  This  parameter indicates the rate of better solutions that must be preserved from one generation
	to another. Default is 0.2 (20%).
}
  \item{highvalue}{
  Parameter for genetic algorithm. In should not be changed
}
  \item{suggestions}{
  Optional parameter for genetic algorithm that indicates a suggested solution to be introduced in the initial population. The most convenient is the one found by the function "KmeanSolution". Default is NULL. 
}
  \item{realAllocation}{
  If FALSE, the allocation is based on INTEGER values; if TRUE, the allocation is based on REAL values.
  Default is TRUE.
}  
  \item{writeFiles}{
  Indicates if the various dataframes and plots produced during the execution have to be written in the working
  directory.
  Default is FALSE.
}
  \item{showPlot}{
  Indicates if the plot showing the trend in the value of the objective function has to be shown or not. In parallel = TRUE, this defaults to FALSE
  Default is TRUE.
}
  \item{parallel}{
  Should the analysis be run in parallel.
  Default is TRUE.
}
  \item{cores}{
  If the analysis is run in parallel, how many cores should be used.
  If not specified n-1 of total available cores are used OR if number of domains < (n-1) cores, then number of cores equal to number of domains are used.
}
}

\value{
A list containing (1) the vector of the solution and (2) the optimal aggregated strata
}

\author{
Giulio Barcaroli
}

\examples{
\dontrun{
library(SamplingStrata)
############################
# Example of "atomic" method
############################
data(swissmunicipalities)
swissmunicipalities$id <- c(1:nrow(swissmunicipalities))
frame <- buildFrameDF(df = swissmunicipalities,
                      id = "id",
                      domainvalue = "REG",
                      X = c("POPTOT","HApoly"),
                      Y = c("Surfacesbois", "Airind"))
ndom <- length(unique(frame$domainvalue))
cv <- as.data.frame(list(DOM = rep("DOM1",ndom),
                         CV1 = rep(0.1,ndom),
                         CV2 = rep(0.1,ndom),
                         domainvalue = c(1:ndom)))
strata <- buildStrataDF(frame)
kmean <- KmeansSolution(strata,cv,maxclusters=30)
nstrat <- tapply(kmean$suggestions, kmean$domainvalue,
                 FUN=function(x) length(unique(x)))
solution <- optimizeStrata(strata = strata,
                        errors = cv,
                        initialStrata = nstrat,
                        suggestions = kmean,
                        iter = 50,
                        pops = 10)
outstrata <- solution$aggr_strata
newstrata <- updateStrata(strata,solution)
framenew <- updateFrame(frame, newstrata)
s <- selectSample(framenew, outstrata)
}
}
\keyword{ survey }

