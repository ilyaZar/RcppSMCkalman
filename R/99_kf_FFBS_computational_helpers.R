initializeX00 <- function(initX, initU, A, B, dimX) {
  if(!is.null(initX)) return(initX)
  flagZero <- FALSE
  if (is.null(B) && is.null(initU)) {
    flagZero <- TRUE
  } else if (is.null(B) || is.null(initU)) {
    msg <- paste("Error when initializing states:, ",
                 "'B' and 'initU' either both NULL or real values!")
    stop(msg)
  }
  if (dimX == 1) {
    BuReg <- ifelse(flagZero, 0, as.numeric(B %*% initU))
    x00 <- BuReg/(1 - A)
  } else {
    BuReg <- ifelse(flagZero,
                    matrix(0, nrow = 1, ncol = dimX),
                    B %*% initU)
    I   <- diag(dimX)
    x00 <- BuReg %*% solve(I - A)
  } <-
  return(x00)
}
initializeP00 <- function(initP, A, Q, dimX) {
  if(!is.null(initP)) return(initP)
  if (dimX == 1) {
    P00 <- Q/(1 - A^2)
  } else {
    I   <- diag(dimX)
    IA  <- solve(I - A)
    P00 <- IA %*% Q %*% t(IA)
  }
  return(P00)
}
computeMatReg <- function(mat, reg, dim, lenT) {
  if (is.null(mat) || is.null(reg)) return(matrix(0, nrow = dim, ncol = lenT))
  mat %*% reg
}
computeXtt1 <- function(A, xtt, BuReg, dimX) {
  if (dimX == 1) return(A*xtt + BuReg)
  A %*% xtt + BuReg
}
computePtt1 <- function(A, Ptt, Q) {
  A %*% tcrossprod(Ptt, A) + Q
}
computeLt <- function(C, Ptt1, R) {
  solve(C %*% tcrossprod(Ptt1, C) + R)
}
computeKt <-function(Ptt1, C) {
  tcrossprod(Ptt1, C)
}
computekG <- function(yObs, C, xtt1, DwReg) {
  yObs - C %*% xtt1 - DwReg
}
computeXtt <- function(xtt1, Kt, Lt, kGain) {
  xtt1 + Kt %*% Lt %*% kGain
}
computePtt <- function(Ptt1, Kt, Lt, t) {
  Ptt <- Ptt1 - Kt %*% tcrossprod(Lt, Kt)
  if (!matrixcalc::is.positive.definite(Ptt)) {
    stop(paste0("matrix is no longer p.d. at iteration number: ", t))
  }
  return(Ptt)
}
computeJtt2 <- function(Ptt, A , Q) {
  tcrossprod(Ptt, A) %*% solve(A %*% tcrossprod(Ptt, A) + Q)
}
