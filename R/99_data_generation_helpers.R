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
  if (numRegs == 0) {
    regs <- rep(0, times = TT)
  } else if (numRegs == 1) {
    regs <- rnorm(n = TT, mean = regMean, sd = sqrt(regVar))
  } else  if (numRegs > 1) {
    regs <- t(mvtnorm::rmvnorm(TT, mean = regMean, sigma = regVar))
  }
  return(regs)
}
simulateRegsU <- function(initU, numU, TT) {
  if (is.null(initU) && numU == 0) {
    u0    <- 0
    regsU <- rep(0, times = TT)
  } else {
    u0    <- initializeRegsU(initU, numU = numU)
    regsU <- sampleRegs(numRegs = numU, TT = TT)
  }
  return(list(initU = u0,
              regsU = regsU))
}
initializeStatesX <- function(x0, u0, numU, A, B, Q) {
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
  epsError <- sampleError(TT, R)
  if (dimY == 1) {
    yObsOut  <- numeric(TT)

    yObsOut[1] <- sum(C * x0) + sum(D * wReg[1]) + epsError[1]
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
  numU <- max(c(1, numU))
  numW <- max(c(1, numW))
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
  numU <- max(c(1, numU))
  numW <- max(c(1, numW))

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
