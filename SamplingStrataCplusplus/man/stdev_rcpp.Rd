\name{stdev_Rcpp}
\Rdversion{1.3}
\alias{stdev_Rcpp}
\title{
Standard deviations of a given Y target variables in sampling strata
}
\description{
Calculates standard deviations of Y target variables in sampling strata
also taking into account spatial correlation
}
\usage{
buildStrataDFSpatial(dataset,
                     i,
                     fitting=c(1),
                     range=c(0),
                     kappa=3)
}
\arguments{
  \item{dataset}{
  This is the name of the dataframe containing the sample data, or the frame data. 
  It is strictly required that auxiliary information is organised in variables named 
  as X1, X2, ... , Xm (there should be at least one of them) and the target variables 
  are denoted by Y1, Y2, ... , Yn.
  In addition, in case of sample data, a variable named 'WEIGHT' must be present in the dataframe, 
  containing   the weigths associated to each sampling unit
	}
\item{i}{
	Integer indicating which Y target variable.
}
\item{fitting}{
	Fitting of the spatial model(s).
	Default is 1.
}
\item{range}{
	Maximum range for spatial autocorrelation
}	
\item{kappa}{
	Factor used in evaluating spatial autocorrelation. 
	Default is 3.
}	

}

\value{
A dataframe containing strata
}
\author{
Giulio Barcaroli
}

\examples{
\dontrun{
stdev_RCpp(dataset=frame,i=1,fitting=1,range=800,kappa=3)
}
}
\keyword{ survey }

