% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05_kf_LL.R
\name{kfLLH}
\alias{kfLLH}
\title{Computes the data distribution}
\usage{
kfLLH(yObs, wReg, xtt1, Ptt1, C, D, R, dimX, dimY, TT, LOG = TRUE)
}
\arguments{
\item{yObs}{A matrix or vector of measurements (observations):
 \itemize{
 \item{rows: multivariate dimension}
 \item{columns: time series dimension \code{T}}
 }
If \code{Y} is a univariate process, \code{yObs} can be passed as a vector
of length \code{T}. If \code{nrow(yObs) = 1}, then \code{yObs} becomes a
vector of length \code{T}.}

\item{wReg}{Matrix (vector) of regressors for the measurement process of
dimension \code{ncol(D) x T}. For a single regressors \code{wReg} is a
vector of length \code{T}.}

\item{xtt1}{predictive means as produced by \code{\link{kfMFPD}} (if
\code{PDSTORE = TRUE})}

\item{Ptt1}{predictive variances as produced by \code{\link{kfMFPD}} (if
\code{PDSTORE = TRUE})}

\item{C}{Parameter (or system) matrix of dimension \code{dimY x dimX}.}

\item{D}{Parameter (or system) matrix of dimension \code{dimY x numW}.}

\item{R}{Error VCM of measurement process of dimension \code{dimY x dimY}}

\item{dimX}{integer giving the dimension of the latent state process}

\item{dimY}{integer giving the dimension of the measurement process}

\item{TT}{integer giving the length of the time series}

\item{LOG}{logical; if \code{TRUE}, then the logarithm of the likelihood is
returned}
}
\value{
the (logarithmic, if \code{LOG=TRUE}) value of the data likelihood
}
\description{
The data distribution, i.e. the (log-)likelihood, given all parameters
\eqn{\theta=A, B, C, D, P, Q}, as described in \emph{Data distribution -
observed likelihood computation} of the \code{Details} section from
\code{\link{kfLGSSM}}.
}
