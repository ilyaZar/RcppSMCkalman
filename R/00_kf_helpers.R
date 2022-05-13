#' Generates data from a Linear Gaussian state space model.
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
#' @param Q VCM of state process of dimension \code{dimX x dimX}.
#' @param R VCM of measurement process of dimension \code{dimY x dimY}
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
#'       \item{\code{measurements:} matrix/vector of simulated measurement
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
#' @export
dataGenLGSSM <- function(TT, dimX = 1, dimY = 1, numU = 0, numW = 0,
                         A = NULL, B = NULL, C = NULL, D = NULL,
                         Q = NULL, R = NULL,
                         initX = NULL, initU = NULL) {
  # browser()
  # to only deal with one regressor case even if numU/numW set to defaults
  numU <- max(c(1, numU))
  numW <- max(c(1, numW))

  setDefaultsIfNULL(dimX, dimY, numU, numW, A, B, C, D, Q, R)
  # Call helper function to check dimensions of user input args
  checkDimMatch(initX, initU, dimX, dimY, numU, numW, A, B, C, D, Q, R)
  # browser()
  # Simulate u-type regressors (those attached to the latent state process)
  uRegOut <- simulateRegsU(initU, numU = numU, TT = TT)
  u0   <- uRegOut[["initU"]]
  uReg <- uRegOut[["regsU"]]

  # Simulate w-type regressors (those attached to the measurement process)
  # browser()
  wReg <- sampleRegs(numRegs = numW, TT = TT)

  # Simulate latent state process
  xStatesOut <- simulateStatesX(TT, dimX, numU,
                                A, B, Q,
                                x0 = initX, u0, uReg)
  x0      <- xStatesOut[["xInit"]]
  xStates <- xStatesOut[["xStates"]]

  # browser()
  # Simulate observation process
  yObs <- simulateObsY(dimY, TT, x0,
                       xStates, wReg,
                       C, D, R)

  return(list(data = list(xStates = xStates, measurements = yObs,
                          uRegs = uReg, wRegs = wReg),
              init = list(initX = x0, initU = u0,
                          A = A, B = B, C = C, D = D, Q = Q, R = R))
  )
}
sampleError <- function(TT, vcm) {
  dimProc <- ifelse(is.null(nrow(vcm)), 1, nrow(vcm))
  if (dimProc == 1) {
    errOut <- rnorm(n = TT, mean = 0, sd = sqrt(vcm))
  } else {
    errOut <- mvtnorm::rmvnorm(TT, mean = rep(0, times = dimProc), sigma = vcm)
  }
  return(errOut)
}
initializeRegsU <- function(u0, numU) {
  if (is.null(u0)) {
    if (numU == 1) {
      u0 <- rnorm(n = 1, mean = 1, sd = sqrt(1))
    } else {
      u0 <- as.vector(mvtnorm::rmvnorm(1, mean = rep(1, times = numU),
                                       sigma = diag(numU)))
    }
  }
  return(u0)
}
sampleRegs <- function(numRegs, TT,
                       regMean = rep(0, times = numRegs),
                       regVar  = diag(numRegs)) {
  if (numRegs == 1) {
    regs <- rnorm(n = TT, mean = regMean, sd = sqrt(regVar))
  } else {
    regs <- t(mvtnorm::rmvnorm(TT, mean = regMean, sigma = regVar))
  }
  return(regs)
}
simulateRegsU <- function(initU, numU, TT) {
  u0    <- initializeRegsU(initU, numU = numU)
  regsU <- sampleRegs(numRegs = numU, TT = TT)
  return(list(initU = u0,
              regsU = regsU))
}
initializeStatesX <- function(x0, u0, numU, A, B, Q) {
  browser()
  if (is.null(x0)) {
    dimX <- ifelse(is.null(nrow(A)), 1, nrow(A))
    if (dimX == 1) {
      x0 <- rnorm(n = 1, mean = sum(B*u0)/(1 - A), sd = sqrt(Q/(1 - A^2)))
    } else {
      invIminA <- solve(diag(dimX) - A)
      BuReg <- computeBDtimesRegs(mat = B, regs = as.matrix(u0))
      x0 <- as.vector(mvtnorm::rmvnorm(1, mean = as.vector(invIminA %*% BuReg),
                                       sigma = invIminA %*% Q %*% t(invIminA)))
    }
  }
  return(x0)
}
sampleX <- function(x0, dimX, TT, A, B, Q, uReg) {
  browser()
  nuError <- sampleError(TT, Q)
  if (dimX == 1) {
    xStatesOut    <- rep(0, times = TT)
    xStatesOut[1] <- A*x0 + sum(B*uReg[1]) + nuError[1]
    for (t in 2:TT) {
      xStatesOut[t] <- A*xStatesOut[t - 1] + sum(B*uReg[t]) + nuError[t]
    }
  } else {
    xStatesOut <- matrix(0, nrow = dimX, ncol = TT)
    BuReg <- computeBDtimesRegs(mat = B, regs = uReg)

    xStatesOut[, 1] <- A %*% x0 + BuReg[, 1] + nuError[1, ]
    for (t in 2:TT) {
      xStatesOut[, t] <- A %*% xStatesOut[, t - 1] + BuReg[, t] + nuError[t, ]
    }
  }
  return(xStatesOut)
}
simulateStatesX <- function(TT, dimX, numU,
                            A, B, Q,
                            x0 = NULL, u0, uReg) {
  x0 <- initializeStatesX(x0, u0, numU, A, B, Q)
  xStates <- sampleX(x0, dimX, TT,
                     A, B, Q,
                     uReg)
  return(list(xInit = x0,
              xStates = xStates))
}
simulateObsY <- function(dimY, TT, x0,
                         xStates, wReg,
                         C, D, R) {
  # browser()
  epsError <- sampleError(TT, R)
  if (dimY == 1) {
    yObsOut  <- numeric(TT)

    yObsOut[1] <- sum(C * x0) + sum(D * wReg[1]) + epsError[t]
    for(t in 1:TT) {
      yObsOut[t] <- sum(C * xStates[t]) + sum(D * wReg[t]) + epsError[t]
    }
  } else {
    yObsOut <- matrix(0, nrow = dimY, ncol = TT)
    DwReg    <- computeBDtimesRegs(mat = D, regs = wReg)
    CxStates <- computeBDtimesRegs(mat = C, regs = xStates)

    # yObsOut[, 1] <- CxStates[, 1] + DwReg[, 1] + epsError[1, ]
    # for(t in 1:TT) {
    #   yObsOut[, t] <- CxStates[, t] + DwReg[, t] + epsError[t, ]
    # }
    yObsOut <- CxStates + DwReg + t(epsError)
  }
  return(yObsOut)
}
computeBDtimesRegs <- function (mat, regs) {
  if (is.matrix(mat) && is.matrix(regs)) {
    matRegs <- mat %*% regs
  } else if (!is.matrix(mat) && !is.matrix(regs)) {
    regs  <- t(as.matrix(regs))
    mat     <- as.matrix(mat)
    matRegs <- mat %*%regs
  } else if (!is.matrix(mat)) {
    mat     <- as.matrix(mat)
    matRegs <- mat %*%regs
  } else {
    regs    <- t(as.matrix(regs))
    matRegs <- mat %*%regs
  }
  return(matRegs)
}
setDefaultsIfNULL <- function(dimX, dimY, numU, numW, A, B, C, D, Q, R,
                              varQ = 0.5, varR = 0.25) {
  if (is.null(A)) {
    if (dimX == 1) {
      A <- 1
    } else if (dimX > 1) {
      A <- diag(dimX)
    }
  }
  if (is.null(B)) {
    if (numU == 1) {
      B <- rep(0, times = dimX)
    } else {
      B <- matrix(0, nrow = dimX, ncol = numU)
    }

  }
  if (is.null(C)) {
    if (dimY == 1) {
      C <- 1
    } else if (dimY > 1) {
      C <- diag(dimY)
    }
  }
  if (is.null(D)) {
    if (numW == 1) {
      D <- rep(0, times = dimY)
    } else {
      D <- matrix(0, nrow = dimY, ncol = numW)
    }
  }
  if (is.null(Q)) Q <- diag(dimX) * varQ
  if (is.null(R)) R <- diag(dimY) * varR

  assign("A", A, envir = parent.frame())
  assign("B", B, envir = parent.frame())
  assign("C", C, envir = parent.frame())
  assign("D", D, envir = parent.frame())
  assign("Q", Q, envir = parent.frame())
  assign("R", R, envir = parent.frame())
}
checkDimMatch <- function(initX, initU, dimX, dimY, numU, numW,
                          A, B, C, D, Q, R) {
  dimX <- as.integer(dimX)
  dimY <- as.integer(dimY)
  numU <- as.integer(numU)
  numW <- as.integer(numW)
  if (!is.null(initX)) {
    stopifnot(`Initial state process unequal to dimX` = length(initX) == dimX)
  }
  if (!is.null(initU)) {
    stopifnot(`Initial u_00 regressors unequal to dimX` = length(initU) == dimX)
  }

  # Check dimension missmatches for matrix A
  matchAX <- identical(dimX, ncol(A)) || identical(dimX, length(A))
  stopifnot(`'dimX' dimension does not match 'A' dimension` = matchAX)

  # Check dimension missmatches for matrix B
  if (is.null(dim(B))) {
    checkLength <- length(B) # since B is a vector either of the following
    if (dimX == 1) { # case where dimX=1 implies numU and length(B) must match
      matchXREG <- identical(numU, checkLength)
      stopifnot(`For dimX=1, 'numU' does not match 'B' length` = matchXREG)
    } else { # case where dimX>1 must mean numU=1 (otherwise is.null(dim(B))
      # evaluates to FALSE), s.th. dimX and length(B) must match
      matchXREG <- identical(dimX, checkLength)
      stopifnot(`'dimX' does not match 'B' length` = matchXREG)
      matchXREG <- matchXREG && (numU == 1)
      stopifnot(`'dimX' matches length 'B' but 'numU>1' ambiguous` = matchXREG)
    }
  } else {
    matchXREG <- identical(dimX, nrow(B))
    stopifnot(`'dimX' does not match row-dimension of matrix 'B'` = matchXREG)
    matchXREG <- identical(numU, ncol(B))
    stopifnot(`'numU' does not match col-dimension of matrix 'B'` = matchXREG)
  }

  # Check dimension missmatches for matrix C
  if (is.null(dim(C))) {
    checkLength <- length(C) # since C is a vector either of the following
    if (dimX == 1) { # case where dimX=1 implies dimY and length(C) must match
      matchYCX <- identical(dimY, checkLength)
      stopifnot(`'Y' (measurement) dimension does not match 'C'` = matchYCX)
    } else { # case where dimX>1 implies that dimY=1 (otherwise is.null(dim(C))
      # must evaluate to FALSE), s.th. dimX and length(C) must match
      matchYCX <- identical(dimX, checkLength)
      stopifnot(`Dimension of states 'X' does not match 'C'` = matchYCX)
    }
  } else {
    matchYCX <- identical(dimX, ncol(C)) && identical(dimY, nrow(C))
    stopifnot(`Dimension missmatch of 'C' vs. 'X' or 'Y'` = matchYCX)
  }

  # Check dimension missmatches for matrix D
  if (is.null(dim(D))) {
    checkLength <- length(D) # since D is a vector either of the following
    if (dimY == 1) {
      matchYREG <- identical(numW, checkLength)
      stopifnot(`For dimY=1, 'numW' does not match 'D' length` = matchYREG)
    } else {# case where dimY>1 must mean numW=1 (otherwise is.null(dim(B))
      # evaluates to FALSE), s.th. dimY and length(D) must match
      matchYREG <- identical(dimY, checkLength)
      stopifnot(`'dimY' does not match 'D' length` = matchYREG)
      matchYREG <- matchYREG && (numW == 1)
      stopifnot(`'dimY' matches length 'D' but 'numW>1' ambiguous` = matchYREG)
    }
  } else {
    matchYREG <- identical(dimX, nrow(D))
    stopifnot(`Dimension missmatch of 'D' vs. 'Y'` = matchYREG)
    matchYREG <- identical(numW, ncol(D))
    stopifnot(`Dimension missmatch of 'D' vs. 'W'` = matchYREG)
  }

  # Check dimension missmatch for matrix Q
  matchQ <- identical(unique(dim(Q)), dimX) || identical(length(Q), dimX)
  stopifnot(`Error VCM 'Q' does not match state ('X') dimension` = matchQ)
  # Check dimension missmatch for matrix R
  matchR <- identical(unique(dim(R)), dimY) || identical(length(Q), dimY)
  stopifnot(`Error VCM 'R' does not match 'Y' (measurement) dimension` = matchR)
}
# kflg <- function(yObs, A, C, Q, R, P00, x00, n_sim = 1) {
#   TT    <- length(yObs)
#   xStates  <- numeric(TT)
#   Ptt  <- numeric(TT)
#   # Ptt1    <- A %*% P00 %*% A  + Q
#   Ptt1    <- tcrossprod(A, tcrossprod(A, P00))
#   Lt <-
#   Kt      <- Ptt1*C*(((C^2)*Ptt1 + R)^(-1))
#   Ptt[1]  <- Ptt1 - Kt*C*Ptt1
#   xStates[1]  <- A*x00 + Kt*(yObs[1] - C*A*x00)
#   for (t in 2:TT) {
#     Ptt1     <- (A^2)*Ptt[t - 1] + Q
#     Kt     <- Ptt1*C*(((C^2)*Pttm + R)^(-1))
#     Ptt[t] <- Ptt1 - Kt*C*Ptt1
#     xStates[t] <- A*xStates[t - 1] + Kt*(yObs[t] - C*A*xStates[t - 1])
#   }
#   # KFapprox <- xStates
#   res_sim <- matrix(0, ncol = TT, nrow = n_sim)
#   for (n in 1:n_sim) {
#     res_sim[n, ] <- rnorm(TT, mean = xStates, sd = sqrt(Ptt))
#   }
#   if (n_sim == 1) {
#     res_sim <- as.vector(res_sim)
#   }
#   return(list(KFx = res_sim, KFvar = Ptt))
# }
