#' Kalman forward filter and backward smoother.
#'
#' The function runs a Kalman (exact) forward (marginal) filter, obtains
#' predictions, and runs a full backward smoother; flags can control the amount
#' of computation (e.g. only obtaining a forward filter, skipping the backward
#' smoother).
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
#' as an initial "prediction" or VCM matrix at period \eqn{t=0} via
#' \code{initP}, while an initial regressor (vector) value \code{initU} must be
#' provided for initialization of the algorithm. If the former two are not
#' provided, the stationary prior is used as a default initializer having the
#' following form:
#' \deqn{x_0\sim\mathcal{N}\left(Bu_0(I-A)^{-1},\left(I-A\right)^{-1}Q
#' \left[\left(I-A\right)^{-1}\right]^{\top}\right)\;,} where \eqn{I} is the
#' identity matrix, and \eqn{u_0} is \code{initU}. In the univariate case
#' (\code{length(initX) = 1}), the prior reduces to
#' \eqn{x_0\sim\mathcal{N}\left(\frac{B\cdot u_0}{1-A},
#' \frac{Q}{(1-A^2)}\right)}, and \eqn{B\cdot u_0} is the "dot" or scalar
#' product if the corresponding number of \code{u}-type regressors \code{>1}.}
#'
#' \subsection{(Marginal) Filtering and prediciton - implemented recursions}{
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
#'       \deqn{K_t = \hat{P}_{t|t-1} C^{\top} = \left(A\hat{P}_{t-1|t-1}
#'        A^{\top} + Q\right) C^{\top}}
#'       \deqn{L_t = \left[C\hat{P}_{t|t-1} C^{\top} + R\right]^{-1} =
#'       \left[C \left(A\hat{P}_{t-1|t-1} A^{\top} + Q\right) C^{\top} +
#'       R\right]^{-1}}
#'       \item \code{for marginal filtering:}
#'        \deqn{\hat{x}_{t|t} = \hat{x}_{t|t-1} + K_tL_t
#'        \left(y_t - C\hat{x}_{t|t-1}-Dw_t\right) = A\hat{x}_{t-1|t-1} +
#'        Bu_{t-1} + K_tL_t\left(y_t -C\left[A\hat{x}_{t-1|t-1} + Bu_{t-1}
#'        \right]-Dw_t\right)}
#'        \deqn{\hat{P}_{t|t} = \hat{P}_{t|t-1} - K_t L_t K_t^{\top} =
#'        A\hat{P}_{t-1|t-1} A^{\top} + Q - K_t L_t K_t^{\top}}
#'       \item \code{for prediction:}
#'       \deqn{\hat{P}_{t|t-1} = A\hat{P}_{t-1|t-1} A^{\top} + Q}
#'       \deqn{\hat{x}_{t|t-1} = A\hat{x}_{t-1|t-1} + Bu_{t-1}}
#'       }
#'   }
#' }
#' }
# '
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
#' @param nSImMfd Number of forward filtering runs; defaults to
#'   \code{nSImMfd=1}.
#' @param nSimJsd Number of backward smoothing runs; defaults to
#'   \code{nSimJsd=1}.
#' @param computeMFD Logical: if \code{TRUE}, then the marginal filtering
#'   density is computed.
#' @param computeJSD Logical: if \code{TRUE}, then the joint smoothing density
#'   is computed.
#'
#' @return A named list of 4:
#'   \itemize{
#'     \item{\code{kfMarginalFilteringDensity:} list of length \code{nSImMfd}
#'     with each element being one complete marginal filtering series.}
#'     \item{\code{kfVCMmfd:} the marginal filtering VCM estimates i.e. the
#'     \eqn{P_{t|t}} matrices, see \code{Details}.}
#'     \item{\code{kfJointSmoothingDensity:} list of length \code{nSimJsd} with
#'     each element being one complete joint smoothing density pass.}
#'     \item{\code{kfVCMjsd:} the joint smoothing VCM estimates i.e. the
#'     \eqn{J_{t|t}} matrices, i.e. see \code{Details}.}
#'   }
#' @export
kfLGSSM <- function(yObs, uReg, wReg,
                    A, B, C, D, Q, R,
                    initX, initP, initU,
                    nSImMfd = 1,
                    nSimJsd = 1,
                    computeMFD = TRUE,
                    computeJSD = TRUE) {
  TT   <- ncol(yObs)
  dimX <- nrow(yObs)
  A <- as.matrix(diag(A))
  B <- as.matrix(diag(B))
  C <- as.matrix(diag(C))
  D <- as.matrix(diag(D))

  xtt <- matrix(0, nrow = dimX, ncol = TT)
  Ptt <- rep(list(list()), times = TT)

  Ptt1    <- tcrossprod(A, tcrossprod(A, initP)) + Q
  Lt      <- solve(tcrossprod(C, tcrossprod(C, Ptt1)) + R)
  Kt      <- Ptt1 %*% t(C) %*% Lt

  Ptt[[1]]  <- Ptt1 - Kt %*% Lt %*% t(Kt)
  if (!matrixcalc::is.positive.definite(Ptt[[1]])) {
    stop(paste0("matrix is no longer p.d. at iteration number: ", 1))
  }
  xtt[, 1]  <- A %*% initX + B %*% initU + Kt %*% (yObs[, 1] - C %*% (A %*% initX + B %*% initU) - D %*% wReg[, 1])
  for (t in 2:TT) {
    Ptt1     <- tcrossprod(A, tcrossprod(A, Ptt[[t - 1]])) + Q
    Lt       <- solve(tcrossprod(C, tcrossprod(C, Ptt1)) + R)
    Kt       <- Ptt1 %*% t(C) %*% Lt

    Ptt[[t]] <- Ptt1 - Kt %*% Lt %*% t(Kt)
    if (!matrixcalc::is.positive.definite(Ptt[[t]])) {
      stop(paste0("matrix is no longer p.d. at iteration number: ", t))
    }
    xtt[, t] <- A %*% xtt[, t - 1] + B %*% uReg[, t - 1] + Kt %*% (yObs[, t] - C %*% (A %*% xtt[, t - 1] + B %*% uReg[, t - 1]) - D %*% wReg[, t])
  }
  if (computeMFD) {
    kfMFDres <- rep(list(list()), times = nSImMfd)
    for (n in 1:nSImMfd) {
      kfMFD <- matrix(0, ncol = TT, nrow = dimX)
      for (t in 1:TT) {
        kfMFD[, t] <- rmvnorm(1, mean = xtt[, t], sigma = Ptt[[t]])
      }
      kfMFDres[[n]] <-kfMFD
    }
  }
  if (computeJSD) {
    kfJSDout <- matrix(0, nrow = dimX*nSimJsd, ncol = TT)
    for (n in 1:nSimJsd) {
      Jtt <- rep(list(list()), times = TT)
      xttJSD <- matrix(0, nrow = dimX, ncol = TT)
      kfJSD  <- matrix(0, nrow = dimX, ncol = TT)

      kfJSD[, TT] <- rmvnorm(1, mean = xtt[, TT], sigma = Ptt[[TT]])

      for (t in (TT - 1):1) {
        Jtt[[t]] <- Ptt[[t]] %*% t(A) %*% solve(tcrossprod(A, tcrossprod(A, Ptt[[t]])) + Q)
        xttJSD[, t]  <- xtt[, t] + Jtt[[t]]  %*% (kfJSD[, t + 1] - B %*% uReg[, t] - A %*% xtt[, t])
        kfJSD[, t] <- rmvnorm(1, mean = xttJSD[, t], sigma = Jtt[[t]])
      }
      # browser()
      print(n)
      rowID <- (dimX*(n - 1) + 1):(dimX*n)
      kfJSDout[rowID, ] <-  kfJSD
    }
  }
  return(list(kfMarginalFilteringDensity = kfMFDres, kfVCMmfd = Ptt,
              kfJointSmoothingDensity = kfJSDout, kfVCMjsd = Jtt))
}
