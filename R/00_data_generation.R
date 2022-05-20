#' Generates data from a linear Gaussian state space model.
#'
#' The data is simulated from a multivariate Linear Gaussian state space model
#' (LGSSM), and is returned as a named list for convenient access.
#'
#' The specification of the LGSSM is of the following form:
#' \deqn{x_t = Ax_{t-1} + Bu_t + v_t \\}
#' \deqn{y_t = Cx_t + Dw_t + \varepsilon_t\;,}
#' with serially uncorrelated state innovations \eqn{v_t\sim\mathcal{N}(0, Q)}
#' and measurement errors \eqn{\varepsilon_t\sim\mathcal{N}(0,R)}, for
#' \eqn{t=1,\ldots,T}. The argument names match the above equations; if an
#' argument is missing,the corresponding component is dropped possibly returning
#' \code{NULL} for missing elements in the output list.
#'
#' There must be an initial state (vector) value passed via \code{initX} and an
#' initial regressor (vector) value \code{initU} upon which the time series of
#' state and measureent (observation) processes, running for
#' \eqn{t=1,\ldots,TT}, is initialized. If the initial state process is not
#' specified, then one sample is drawn internally from the prior distribution
#' \deqn{x_0\sim\mathcal{N}\left(Bu_0(I-A)^{-1},\left(I-A\right)^{-1}Q
#' \left[\left(I-A\right)^{-1}\right]^{\top}\right)\;,}
#' where
#' \eqn{I} is the identity matrix, and \eqn{u_0} is \code{initU}. If
#' \code{initU} is missing, then \eqn{u_0\sim\mathcal{N}(0,I)} of dimension
#' \code{dimX} is used internally. In the univariate case (\code{dimX = 1}), the
#' prior reduces to
#' \eqn{x_0\sim\mathcal{N}\left(\frac{B\cdot u_0}{1-A},
#' \frac{Q}{(1-A^2)}\right)},
#' and \eqn{B\cdot u_0} is the "dot" or scalar product if the corresponding
#' number of \code{u}-type regressors \code{>1}.
#'
#' @param TT Integer giving the time series length.
#' @param dimX Integer giving the state process dimension; defaults to
#'   \code{dimX = 1}.
#' @param dimY Integer giving the measurement process dimension; defaults to
#'   \code{dimX = 1}.
#' @param numU Integer giving the number of regressors attached to the
#'   state process; defaults to \code{numU = 0}.
#' @param numW Integer giving the number of regressors attached to the
#'   measurement process; defaults to \code{numW = 0}.
#' @param A Parameter (or system) matrix of dimension \code{dimX x dimX}.
#' @param B Parameter (or system) matrix of dimension \code{dimX x numU}.
#' @param C Parameter (or system) matrix of dimension \code{dimY x dimX}.
#' @param D Parameter (or system) matrix of dimension \code{dimY x numW}.
#' @param Q Error VCM of state process of dimension \code{dimX x dimX}.
#' @param R Error VCM of measurement process of dimension \code{dimY x dimY}
#' @param initX Initial value for state process. Think of \eqn{X_{t=0}} as a
#'   starting/initial condition for a time series running for \code{t=1,...,T}.
#'   If not specified, then \code{initX} is sampled from the prior under the
#'   stationary distribution of the latent state process (see \code{Details}).
#' @param initU Initial values of regressors of the state process specification
#'   (must have length equal to \code{numU}).
#'
#' @return A named list with elements "\code{data}" and "\code{init}" each being
#'   a (named) list of the following form:
#'   \itemize{
#'     \item{\code{data:}}
#'     \itemize{
#'       \item{\code{xStates:} matrix/vector of the simulated state process of
#'         dimension \code{dimX x TT}.}
#'       \item{\code{yMeasurements:} matrix/vector of simulated measurement
#'         process of dimension \code{dimY x TT}.}
#'       \item{\code{uRegs:} matrix/vector of simulated regressors attached to
#'         the state process of dimension \code{numU x TT}.}
#'       \item{\code{wRegs:} matrix/vector of simulated regressors attached to
#'         the measurement process of dimension \code{numW x TT}.}
#'     }
#'     \item{\code{init:}}
#'     \itemize{
#'       \item \code{initX}, same as argument, see \code{Argument} list
#'       \item \code{initU}, same as argument, see \code{Argument} list
#'       \item \code{A}, same as argument, see \code{Argument} list
#'       \item \code{B}, same as argument, see \code{Argument} list
#'       \item \code{C}, same as argument, see \code{Argument} list
#'       \item \code{D}, same as argument, see \code{Argument} list
#'       \item \code{Q}, same as argument, see \code{Argument} list
#'       \item \code{R}, same as argument, see \code{Argument} list
#'     }
#'
#'   }
#'
#' @seealso To run Kalman filtering, smoothing, prediction and likelihood
#'    computations, refer to \code{\link{kfLGSSM}} and the sub-functions
#'    therein.
#'
#' @export
dataGenLGSSM <- function(TT, dimX = 1, dimY = 1, numU = 0, numW = 0,
                         A = NULL, B = NULL, C = NULL, D = NULL,
                         Q = NULL, R = NULL,
                         initX = NULL, initU = NULL) {
  setDefaultsIfNULL(dimX, dimY, numU, numW, A, B, C, D, Q, R)
  # Call helper function to check dimensions of user input args
  checkDimMatch(initX, initU, dimX, dimY, numU, numW, A, B, C, D, Q, R)
  # Simulate u-type regressors (those attached to the latent state process)
  uRegOut <- simulateRegsU(initU, numU = numU, TT = TT)
  u0   <- uRegOut[["initU"]]
  uReg <- uRegOut[["regsU"]]

  # Simulate w-type regressors (those attached to the measurement process)
  wReg <- sampleRegs(numRegs = numW, TT = TT)

  # Simulate latent state process
  xStatesOut <- simulateStatesX(TT, dimX, numU,
                                A, B, Q,
                                x0 = initX, u0, uReg)
  x0      <- xStatesOut[["xInit"]]
  xStates <- xStatesOut[["xStates"]]

  # Simulate observation process
  yObs <- simulateObsY(dimY, TT, x0,
                       xStates, wReg,
                       C, D, R)

  return(list(data = list(xStates = xStates, yMeasurements = yObs,
                          uRegs = uReg, wRegs = wReg),
              init = list(initX = x0, initU = u0,
                          A = A, B = B, C = C, D = D, Q = Q, R = R))
  )
}
