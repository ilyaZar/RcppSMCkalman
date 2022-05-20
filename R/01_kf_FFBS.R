#' Kalman forward filtering, prediction, smoothing and likelihood evaluation
#'
#' Under a Linear Gaussian State Space Model (LGSSM) specification, the function
#' runs an (exact) Kalman forward (marginal) filter, computes the prediction,
#' marginal and joint smoothing densities (the former sometimes referred to as
#' the RST-smoother), and generates an unbiased estimate of the
#' (log-)likelihood.
#'
#' The amount of computation can be controlled via the \code{computeXXX}-flags
#' (e.g. only a forward filter and log-likelihood computations, but skipping
#' prediction the smoothers). If required the \code{nSimXX}-type arguments set
#' the number of simulations from above densities; the defaults are set to 1,
#' but if set to 0, then no simulation is performed and only the means and
#' variances of the Gaussian distributions of the \code{computeXXX}-quantities
#' are returned.
#'
#' @details
#' \subsection{Model specification}{
#'
#' The specification of the LGSSM is of the following form:
#' \deqn{x_t = Ax_{t-1} + Bu_t + v_t \\}
#' \deqn{y_t = Cx_t + Dw_t + \varepsilon_t\;,}
#' with serially uncorrelated state innovations \eqn{v_t\sim\mathcal{N}(0, Q)}
#' and measurement errors \eqn{\varepsilon_t\sim\mathcal{N}(0,R)}, for
#' \eqn{t=1,\ldots,T}. The argument names match the previous equations. If an
#' argument is missing, the corresponding component is dropped (i.e. set to
#' zero), and forward filtering, prediction and backward smoothing are run under
#' this specification as usual.
#'
#' There can be an initial state (vector) value passed via \code{initX} as well
#' as an initial "prediction" or VCM matrix at period \eqn{t=0} via \code{initP}
#' for initialization of the algorithm. If these are not provided, the
#' stationary prior is used as a default initializer having the following form:
#' \deqn{x_0\sim\mathcal{N}\left(Bu_0(I-A)^{-1},\left(I-A\right)^{-1}Q
#' \left[\left(I-A\right)^{-1}\right]^{\top}\right)\;,} where \eqn{I} is the
#' identity matrix, and \eqn{u_0} is \code{initU}. In the univariate case
#' (\code{length(initX) = 1}), the prior reduces to
#' \eqn{x_0\sim\mathcal{N}\left(\frac{B\cdot u_0}{1-A},
#' \frac{Q}{(1-A^2)}\right)}, and \eqn{B\cdot u_0} is the "dot" or scalar
#' product if the corresponding number of \code{u}-type regressors \code{>1}.
#'
#' Note that above specification requires an initial regressor (vector) value
#' \code{initU} if \code{initX} is missing, but can be dropped if \code{initX}
#' is provided.
#' }
#'
#' \subsection{(Marginal) Filtering and prediction - implemented recursions}{
#' \enumerate{
#'   \item {\code{For }\eqn{t=0:} \code{compute} \eqn{\hat{x}_{0|0}}
#'   and \eqn{\hat{x}_{0|0}} \code{according to:}
#'     \deqn{\hat{x}_{0|0} = Bu_0\left(I-A\right)^{-1}}
#'     \deqn{\hat{P}_{0|0} =\left(I-A\right)^{-1}Q
#'     \left[\left(I-A\right)^{-1}\right]^{\top}}
#'   }
#'   \item {\code{For }\eqn{t=1,\ldots,T:} \code{compute:}
#'     \itemize{
#'       \item \eqn{K_t} and \eqn{L_t}
#'       \deqn{K_t = \hat{P}_{t|t-1} C^{\top} }
#      # = \left(A\hat{P}_{t-1|t-1}
#     #  A^{\top} + Q\right) C^{\top}}
#'       \deqn{L_t = \left[C\hat{P}_{t|t-1} C^{\top} + R\right]^{-1}}
#'      with \eqn{\hat{x}_{t|t-1}} and \eqn{\hat{P}_{t|t-1}} as specified
#'      below under \code{prediction}.
#      =
#       \left[C \left(A\hat{P}_{t-1|t-1} A^{\top} + Q\right) C^{\top} +
#      R\right]^{-1}}
#'       \item \code{for marginal filtering:}
#'        \deqn{\hat{x}_{t|t} = \hat{x}_{t|t-1} + K_tL_t
#'        \left(y_t - C\hat{x}_{t|t-1}-Dw_t\right)}
#         = A\hat{x}_{t-1|t-1} +
#        Bu_{t-1} + K_tL_t\left(y_t -C\left[A\hat{x}_{t-1|t-1} + Bu_{t-1}
#        \right]-Dw_t\right)}
#'        \deqn{\hat{P}_{t|t} = \hat{P}_{t|t-1} - K_t L_t K_t^{\top}}
#'        with \eqn{\hat{x}_{t|t-1}} and \eqn{\hat{P}_{t|t-1}} as specified
#'         below under \code{prediction}.
#        =
#        A\hat{P}_{t-1|t-1} A^{\top} + Q - K_t L_t K_t^{\top}}
#'      \item \code{for prediction:}
#'       \deqn{\hat{P}_{t|t-1} = A\hat{P}_{t-1|t-1} A^{\top} + Q}
#'       \deqn{\hat{x}_{t|t-1} = A\hat{x}_{t-1|t-1} + Bu_{t-1}}
#'       }
#'   }
#' }
#' These recursions are implemented in \code{\link{kfMFPD}}.
#' }
#' \subsection{Marginal smoothing - implemented recursions}{
#'  To obtain the marginal smoothing densities
#' \eqn{p\left(x_t|y_{1:T}\right)\forall t=1,\ldots,T} a complete Kalman forward
#' filtering and prediction pass is needed. The marginal smoothing density for
#' period \eqn{T} equals the marginal filtering density for the same period,
#' i.e. after the Kalman filtering and prediction pass \eqn{x_{T|T}} and
#' \eqn{P_{T|T}} for \eqn{p\left(x_T|y_{1:T}\right)} are already given. To
#' compute the other marginals, the following recursion is run backwards in
#' time:
#'
#' \code{For }\eqn{t=T-1, \ldots,1}\code{ compute:}
#'  \enumerate{
#'    \item \eqn{J_t=P_{t|t}A^{\top}P_{t|t-1}^{-1}}
#'    \item \eqn{\hat{x}_{t|T}=\hat{x}_{t|t} + J_t \left(\hat{x}_{t+1|T}-
#'    \hat{x}_{t+1|t}\right)}
#'    \item \eqn{\hat{P}_{t|T}=\hat{P}_{t|t} + J_t\left(\hat{x}_{t+1|T} -
#'    \hat{x}_{t+1|t}\right)J_t^{\top}}
#'  }
#'
#' Each marginal smoothing density is then given via
#' \eqn{p\left(x_t|y_{t:T}\right)=\mathcal{N}\left(x_t|\hat{x}_{t|T},
#' \hat{P}_{t|T}\right)}.
#'
#' These recursions are implemented in \code{\link{kfMSD}}.
#'
#' }
#' \subsection{Joint smoothing density- implemented backward simulation
#' recursions}{
#'   Backward simulation is a strategy for generating realizations of a whole
#'   trajectory \eqn{\tilde{x}_{0:T}} from the joint smoothing density
#'   \eqn{p\left(x_{0:T}|y_{1:}\right)}:
#'   \itemize{
#'   \item \code{Draw nSimJSD samples from} \deqn{\tilde{x}_T\sim
#'   p\left(x_T|y_{1:T}\right)=\mathcal{N}\left(x_T|\hat{x}_{T|T}|P_{T:T}\right)}
#'   i.e. from the last iteration of the forward filter.
#'   \item \code{Draw nSimJSD samples from backwards in time for}
#'   \eqn{t=T-1,\ldots,0}:
#'   \deqn{\tilde{x}_t\sim p\left(x_t|\tilde{x}_{t+1},y_{1:t}\right)
#'   =\mathcal{N}\left(x_t|\hat{\mu}_{t},L_t\right)}}
#'   with means and variances obtained
#'   via the following recursions:
#'   \itemize{
#'     \item \eqn{J_t=\hat{P}_{t|t}A^{\top}\left(A\hat{P}_{t|t}A^{\top} +
#'     Q\right)^{-1}}
#'     \item \eqn{\mu_t=\hat{x}_{t|t} + J_t \left(x_{t+1} - Bu_t -A\hat{x}_{t|t}
#'     \right)}
#'     \item \eqn{L_t=\hat{P}_{t|t} - J_t A\hat{P}_{t|t}}
#'   }
#'
#'  After a complete backward sweep, the \code{nSimJSD} backward trajectories
#'  \eqn{\tilde{x}_{0:T}^{nSimJS}} are (by construction) a realization from the
#'   above joint smoothing density. These recursions are implemented in
#'   \code{\link{kfJSD}}.
#'
#' }
#' \subsection{Data distribution - observed likelihood computation}{
#'   The data likelihood (if logarithmized also sometimes simply referred to as
#'   the log-likelihood), with \eqn{\theta=A, B, C, D, P, Q}, is defined as
#'   \deqn{p\left(y_{1:T}|\theta\right)\prod_{t=1}^{T} \int p\left(y_t|
#'   \theta,x_t\right)p\left(x_t|y_{t:(t-1)}\right)dx_t=\mathcal{N}
#'   \left(y_t|C\hat{x}_{t|t-1}+Dw_t, C\hat{P}_{t|t-1}C^{\top} + R\right)}
#'   with the logarithmic version given as
#'   \deqn{\log p\left(y_{1:T}\right) =
#'   \sum_{t=1}^T\log\left(\mathcal{N}
#'   \left(y_t|C\hat{x}_{t|t-1}+Dw_t, C\hat{P}_{t|t-1}C^{\top} + R\right)\right)
#'   =
#'   -\frac{Tn_y}{2}\log\left(2\pi\right) + \sum_{t=1}^T
#'   \left(-\frac{1}{2} \log\det\lbrace\Lambda_t\rbrace -
#'   \frac{1}{2} \left(y_t-C\hat{x}_{t|t-1}-Dw_t\right)
#'   \Lambda_t^{-1}
#'   \left(y_t-C\hat{x}_{t|t-1}-Dw_t\right)^{\top}
#'    \right)
#'   }
#' with \eqn{\Lambda_t = C\hat{P}_{t|t-1}C^{\top} + R}. These computations are
#' implemented in \code{\link{kfLLH}}.
#' }
#' @param yObs A matrix or vector of measurements (observations):
#'   \itemize{
#'   \item{rows: multivariate dimension}
#'   \item{columns: time series dimension \code{T}}
#'   }
#'  If \code{Y} is a univariate process, \code{yObs} can be passed as a vector
#'  of length \code{T}. If \code{nrow(yObs) = 1}, then \code{yObs} becomes a
#'  vector of length \code{T}.
#' @param uReg Matrix (vector) of regressors for the latent state process of
#'   dimension \code{ncol(B) x T}. For a single regressors \code{uReg} is a
#'   vector of length \code{T}.
#' @param wReg Matrix (vector) of regressors for the measurement process of
#'   dimension \code{ncol(D) x T}. For a single regressors \code{wReg} is a
#'   vector of length \code{T}.
#' @inheritParams dataGenLGSSM
#' @param initP VCM matrix initialization for state process; if not specified,
#'   then the prior covariance-matrix under the stationary distribution of the
#'   latent states is used (see \code{Details}).
#' @param initU Initial values of regressors of the state process specification
#'   (cannot be missing).
#' @param nSimMF Number of forward filtering runs; defaults to \code{nSimMF=1}.
#' @param nSimPD Number of prediction runs; defaults to \code{nSimPD=1}.
#' @param nSimMS Number of marginal smoothing runs; defaults to \code{nSimMS=1}.
#' @param nSimJS Number of joint smoothing (backward simulation) runs; defaults
#'   to \code{nSimJS=1}.
#' @param computeMFD Logical: if \code{TRUE}, then the marginal filtering
#'   density (means and variances) and \code{nSimMF} forward
#'   filtering runs are computed and returned.
#' @param computePDD Logical: if \code{TRUE}, then the prediction density
#'   (means and variances) and \code{nSimPD} predictions are
#    are computed and returned.
#' @param computeMSD Logical: if \code{TRUE}, then the marginal smoothing
#'   density (means and variances) and \code{nSimMS} marginal
#'   smoothing runs are computed and returned.
#' @param computeJSD Logical: if \code{TRUE}, then the joint smoothing density
#'   (means and variances) and \code{nSimJS} joint smoothing (i.e
#'   backward simulation) runs are computed and returned.
#' @param computeLLH Logical: if \code{TRUE}, then the log-likelihood (data
#'   density) is returned.
#'
#' @return A named list of 4:
#'   \itemize{
#'     \item{\code{kfMFDout:}} if \code{computeMFD=TRUE} the output produced by
#'        \code{\link{kfMFPD}} with \code{PDSTORE=FALSE} (or \code{NULL} if
#'        \code{computeMFD=FALSE})
#'     \item{\code{kfMPDout:}} if \code{computePDD=TRUE} the output produced by
#'        \code{\link{kfMFPD}} with \code{PDSTORE=TRUE} (or \code{NULL} if
#'        \code{computePDD=FALSE})
#'     \item{\code{kfMSDout:}} if \code{computeMSD=TRUE} the output produced by
#'        \code{\link{kfMSD}} (or \code{NULL} if \code{computeMSD=FALSE})
#'     \item{\code{kfJSDout:}} if \code{computeJSD=TRUE} the output produced by
#'        \code{\link{kfJSD}} (or \code{NULL} if \code{computeJSD=FALSE})
#'     \item {\code{kfLLHout:}} if \code{computeMFD=TRUE} the output produced by
#'        \code{\link{kfLLH}} (or \code{NULL} if \code{computeMFD=FALSE})
#'   }
#' @seealso If simulated data is required, the function
#'   \code{\link{dataGenLGSSM}} can be used to simulate from a LGSSM.
#'
#'    To use real data that comes with this package: \code{data(IBMdata)}, see
#'    \code{\link{IBMdata}}.
#'
#' @export
kfLGSSM <- function(yObs,
                    uReg = NULL, wReg = NULL,
                    A = NULL, B = NULL,
                    C = NULL, D = NULL,
                    Q = NULL, R = NULL,
                    initX = NULL,
                    initP = NULL,
                    initU = NULL,
                    nSimMF = 1, nSimPD = 1, nSimMS = 1, nSimJS = 1,
                    computeMFD = TRUE,
                    computePDD = FALSE, computeMSD = FALSE,
                    computeJSD = FALSE, computeLLH = FALSE) {
  # 0. Housekeeping and argument checks
  if(missing(yObs)) stop("Measurements/Observations must be provided via yObs.")
  if(is.vector(yObs)) yObs <- matrix(yObs, nrow = 1)
  checkIsMatrix(list("yObs" = yObs, "A" = A, "B" = B, "C" = C, "D" = D,
                     "Q" = Q, "R" = R))
  out = list()
  # 0.1 infer dimension/length from system matrices and data
  TT   <- ifelse(is.matrix(yObs), ncol(yObs), length(yObs))
  dimY <- ifelse(is.matrix(yObs), nrow(yObs), 1)
  dimX <- getDimX(initX, A, B, Q)
  numU <- getNumReg(uReg)
  numW <- getNumReg(wReg)
  # 0.2 check for consistency of arguments: system matrices vs data dimensions
  checkArgumentInputs(TT, dimY, dimX, numU, numW, # true reference for dimension
                      yObs, uReg, wReg, # check all data for correct dimension
                      A, B, C, D, Q, R, # check all system matrices for the same
                      initX, initP, initU) # check all initial vals for the same
  # 0.3 infer model dimension case from Y and X to speed up computation later
  dimCase <- setDimCase(dimX, dimY, silent = FALSE)
  # 0.4 decide upon flags which computations to perform (smooth, filter, etc...)
  cmpCase <- setCmpCase(computeMFD, computePDD,
                        computeMSD, computeJSD,
                        computeLLH)
  # 1. Initialization at t=0:
  u00 <-initU # just pass over and rename for internal naming consistency
  x00 <-initializeX00(initX, initU, A, B, dimX)
  P00 <-initializeP00(initP, A, Q, dimX)

  if (1 %in% cmpCase) {
    out$kfMFDout <- kfMFPD(yObs, uReg, wReg, dimX, dimY, TT,
                           x00, P00, u00,
                           A, B, C, D, Q, R,
                           PDSTORE = FALSE)
  } else {
    out$kfMFDout <- NULL
  }
  if (2 %in% cmpCase) {
    out$kfMPDout <- kfMFPD(yObs, uReg, wReg, dimX, dimY, TT,
                            x00, P00, u00,
                            A, B, C, D, Q, R,
                            PDSTORE = TRUE)
  } else {
    out$kfMPDout <- NULL
  }
  if (3 %in% cmpCase) {
    out$kfMSDout <- kfMSD(dimX, TT,
                          xtt  = out$kfMPDout$mfdEXP,
                          Ptt  = out$kfMPDout$mfdVAR,
                          xtt1 = out$kfMPDout$pddEXP,
                          Ptt1 = out$kfMPDout$pddVAR,
                          A, Q)
  } else {
    out$kfMSDout <- NULL
  }
  if (4 %in% cmpCase) {
    out$kfJSDout  <- kfJSD(uReg = uReg,
                           xtt  = out$kfMFDout$mfdEXP,
                           Ptt  = out$kfMFDout$mfdVAR,
                           A, B, Q, dimX, TT, nSimJS = nSimJS)
  } else {
    out$kfJSDout <- NULL
  }
  if (5 %in% cmpCase) {
    out$kfLLHout <- kfLLH(yObs, wReg,
                          out$kfMPDout$pddEXP,
                          out$kfMPDout$pddVAR,
                          C, D, R,
                          dimX, dimY, TT)
  } else {
    out$kfLLHout <- NULL
  }
  return(out)
}

