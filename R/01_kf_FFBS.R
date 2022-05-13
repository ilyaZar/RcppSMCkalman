#' Kalman forward filtering, prediction and smoothing.
#'
#' Under a Linear Gaussian State Space Model (LGSSM) specification, the function
#' runs an (exact) Kalman forward (marginal) filter, obtains predictions,
#' computes marginal and joint smoothing densities (the former sometimes
#' referred to as the RST-smoother), and generates an unbiased estimate of the
#' (log-)likelihood. \code{computeXXX}-flags control the amount of computation
#' (e.g. only a forward filter and log-likelihood computations, but skipping
#' prediction the smoothers), while \code{nSimXX} arguments set the number of
#' simulations from above densities (if required).
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
#'  vector of length \code{T}. This argument can not be missing.
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
#' @param computePRD Logical: if \code{TRUE}, then the prediction density
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
#'     \item{\code{kfMarginalFilteringDensity:} list of length \code{nSimMF}
#'     with each element being one complete marginal filtering series.}
#'     \item{\code{kfVCMmfd:} the marginal filtering VCM estimates i.e. the
#'     \eqn{P_{t|t}} matrices, see \code{Details}.}
#'     \item{\code{kfJointSmoothingDensity:} list of length \code{nSimJS} with
#'     each element being one complete joint smoothing density pass.}
#'     \item{\code{kfVCMjsd:} the joint smoothing VCM estimates i.e. the
#'     \eqn{J_{t|t}} matrices, i.e. see \code{Details}.}
#'   }
#' @export
kfLGSSM <- function(yObs, uReg = NULL, wReg = NULL,
                    A = NULL, B = NULL, C = NULL, D = NULL, Q = NULL, R = NULL,
                    initX = NULL, initP = NULL, initU = NULL,
                    nSimMF = 1, nSimPD = 1, nSimMS = 1, nSimJS = 1,
                    computeMFD = TRUE, computePRD = TRUE,
                    computeMSD = TRUE, computeJSD = TRUE,
                    computeLLH = TRUE) {
  # 0. Housekeeping and argument checks
  if(missing(yObs)) stop("Measurements/Observations must be provided via yObs.")
  checkIsMatrix(list("A" = A, "B" = B, "C" = C, "D" = D, "Q" = Q, "R" = R))
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
  cmpCase <- setCmpCase(computeMFD, computePRD,
                        computeMSD, computeJSD,
                        computeLLH)
  # 1. Initialization
  A <- as.matrix(diag(A))
  B <- as.matrix(diag(B))
  C <- as.matrix(diag(C))
  D <- as.matrix(diag(D))

  xtt <- matrix(0, nrow = dimX, ncol = TT)
  Ptt <- rep(list(list()), times = TT)

  Ptt1    <- tcrossprod(A, tcrossprod(A, initP)) + Q
  Lt      <- solve(tcrossprod(C, tcrossprod(C, Ptt1)) + R)
  Kt      <- Ptt1 %*% t(C) %*% Lt

  # 1. Initialization
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
    kfMFDres <- rep(list(list()), times = nSimMF)
    for (n in 1:nSimMF) {
      kfMFD <- matrix(0, ncol = TT, nrow = dimX)
      for (t in 1:TT) {
        kfMFD[, t] <- rmvnorm(1, mean = xtt[, t], sigma = Ptt[[t]])
      }
      kfMFDres[[n]] <-kfMFD
    }
  }
  if (computeJSD) {
    kfJSDout <- matrix(0, nrow = dimX*nSimJS, ncol = TT)
    for (n in 1:nSimJS) {
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
getDimX <- function(initX, A, B, Q) {
  if(!is.null(initX)) {
    dimX <- length(initX)
  } else if(!is.null(A)) {
    dimX <- ifelse(is.matrix(A), nrow(A), 1)
  } else if(!is.null(Q)) {
    dimX <- ifelse(is.matrix(Q), nrow(Q), 1)
  } else if(!is.null(B)) {
    dimX <- ifelse(is.matrix(B), nrow(B), 1)
  } else {
    stop("Model miss-specified: can not infer state process dimension.")
  }
  return(dimX)
}
getNumReg <- function(reg) {
  ifelse(test = is.matrix(reg),
         yes  = nrow(reg),
         no   = ifelse(is.null(reg), 0, 1))
}
checkArgumentInputs <- function(TT, dimY, dimX, numU, numW,
                                yObs, uReg, wReg,
                                A, B, C, D, Q, R,
                                initX, initP, initU) {
  ### check matching dimension/length for T
  ## skip check for yObs since TT is inferred from measurement length/dim
  ## check uReg series length
  checkDim(TT, uReg, "col", "T", "U-type regressors")
  ## check wReg series length
  checkDim(TT, wReg, "col", "T", "W-type regressors")
  ### check matching dimension/length for dimY
  ## skip check for yObs since dimY is inferred from measurement length/dim
  checkDim(dimY, wReg, "row", "dimY", "W-type regressors")
  checkDim(dimY, C, "row", "dimY", "'C' matrix")
  checkDim(dimY, D, "row", "dimY", "'D' matrix")
  checkDim(dimY, R, "row", "dimY", "'R' matrix")
  checkSym(R, "'R'")
  ### check matching dimension/length for dimX
  ## skip check for initX since dimX is inferred from state length/dim
  checkDim(dimX, uReg, "row", "dimX", "U-type regressors")
  checkDim(dimX, A, "row", "dimX", "'A' matrix")
  checkSym(A, "'A'")
  checkDim(dimX, B, "col", "dimX", "'B' matrix")
  checkDim(dimX, Q, "col", "dimX", "'Q' matrix")
  checkSym(Q, "'Q'")
  checkDim(dimX, initP, "col", "dimX", "'initP' matrix")
  ### check matching dimension/length for numU
  checkDim(numU, uReg, "row", "numU", "U-type regressors")
  checkDim(numU, B, "col", "numU", "'B' matrix")
  msg <- paste("Initial values for U-type regressors 'initU' do not match",
               "number of regressors derived from matrix B")
  if(length(initU) != numU) stop(msg)
  ### check matching dimension/length for numW
  checkDim(numW, wReg, "row", "numW", "U-type regressors")
  checkDim(numW, D, "col", "numW", "'D' matrix")
}
checkDim <- function(dim, mat, type = "row", nameDim, nameMat) {
  if (!is.null(mat)) {
    msg <- paste("Wrong dimension for", nameMat,
                 ": does not match", nameDim, ".")
    if (type == "row") {
      checkMe <- ifelse(is.matrix(mat), nrow(mat), 1)
    } else if (type == "col") {
      checkMe <- ifelse(is.matrix(mat), ncol(mat), 1)
    } else {
      stop("Unknown check type: use only 'col' or 'row'.")
    }
    if (checkMe != dim) stop(msg)
  }
}
checkSym <- function(mat, matName) {
  msg <- paste(matName, "matrix must be symmetric.")
  checkMe <- ((length(mat) %in% c(0, 1)) || nrow(mat) == ncol(mat))
  if(!checkMe) stop(msg)
}
checkIsMatrix <- function(listOfMatrices) {
  numMat <- length(listOfMatrices)
  namMat <- names(listOfMatrices)
  for (i in 1:numMat) {
    if(isTRUE(is.null(listOfMatrices[[i]]))) next
    msg <- paste0("Matrix '", namMat[i],
                  "' is neither a matrix nor default NULL.")
    if(isFALSE(is.matrix(listOfMatrices[[i]]))) stop(msg)
  }
}
setDimCase <- function(dimX, dimY, silent = FALSE) {
  if (dimX == 1 && dimY == 1) dimCase <- c("dimX=1 && dimY=1" = 1)
  if (dimX  > 1 && dimY == 1) dimCase <- c("dimX>1 && dimY=1" = 2)
  if (dimX == 1 && dimY >  1) dimCase <- c("dimX=1 && dimY>1" = 3)
  if (dimX >  1 && dimY >  1) dimCase <- c("dimX>1 && dimY>1" = 4)
  if(isTRUE(silent)) return(dimCase)
  print(paste0("The following model dimensions are inferred: ",
               names(dimCase), "."))
  return(dimCase)
}
setCmpCase <- function(computeMFD, computePRD,
                       computeMSD, computeJSD,
                       computeLLH, silent = FALSE) {
  if (isTRUE(computeMFD)) cmpCase <- ("Plain MFD computation" = 1)
  if (isTRUE(computePRD)) cmpCase <- ("Plain PRD computation" = 2)
  if (isTRUE(computeMSD)) cmpCase <- ("Plain MSD computation" = 3)
  if (isTRUE(computeJSD)) cmpCase <- ("Plain JSD computation" = 4)
  if (isTRUE(computeLLH)) cmpCase <- ("Plain LLH computation" = 5)
  if(isTRUE(silent)) return(cmpCase)
  print(paste0("The following computational task is inferred: ",
               names(cmpCase), "."))
  return(cmpCase)
}
