\name{finner.prod.iter}
\alias{finner.prod.iter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A fast implementation of the CUSUM statistic

}
\description{
The function returns the CUSUM statistic of a time series by calling the C function fast_inner_prod
}
\usage{
finner.prod.iter(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A time series
}
}

\references{
K. Korkas and P. Fryzlewicz (2015), Multiple change-point detection for non-stationary time series using Wild Binary Segmentation. (http://stats.lse.ac.uk/fryzlewicz/WBS_LSW/WBS_LSW.pdf)
}
\author{
K. Korkas and P. Fryzlewicz
}
\examples{ 
cps=seq(from=1000,to=2000,by=200)
y=sim.pw.arma(N =3000,sd_u = c(1,1.5,1,1.5,1,1.5,1),
b.slope=rep(0.99,7),b.slope2 = rep(0.,7), mac = rep(0.,7),br.loc = cps)[[2]]
z=ews.trans(y,scales=c(11,9,8,7,6))
ts.plot(abs(finner.prod.iter(z[10:2990,2])))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ piecewise ARMA }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
