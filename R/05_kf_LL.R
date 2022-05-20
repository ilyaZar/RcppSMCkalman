#' Computes the data distribution
#'
#' The data distribution, i.e. the (log-)likelihood, given all parameters
#' \eqn{\theta=A, B, C, D, P, Q}, as described in \emph{Data distribution -
#' observed likelihood computation} of the \code{Details} section from
#' \code{\link{kfLGSSM}}.
#'
#' @inheritParams kfMFPD
#' @inheritParams kfMSD
#' @inheritParams kfLGSSM
#' @param dimY integer specifying the dimension of the measurement process
#' @param LOG logical; if \code{TRUE}, then the logarithm of the likelihood is
#'   returned
#'
#' @return the (logarithmic, if \code{LOG=TRUE}) value of the data likelihood
#' @export
kfLLH <- function(yObs, wReg, xtt1, Ptt1, C, D, R, dimX, dimY, TT, LOG = TRUE) {
  DwReg <- computeMatReg(mat = D, reg = wReg, dim = dimX, lenT = TT)
  part1 <- -TT*dimY/2 * log(2*pi)
  part2 <- 0

  for (t in 1:TT) {
    meanY      <- C %*% xtt1[, t] + DwReg[, t]
    VarY       <- C %*% tcrossprod(Ptt1[, , t], C) + R
    logDetVarY <- determinant(VarY, logarithm = TRUE)$modulus[1]

    part2 <- part2 + logDetVarY
    part2 <- part2 + (yObs[, t] - meanY) %*% tcrossprod(solve(VarY),
                                                        (yObs[, t] - meanY))
  }
  llOUT <- part1 - 0.5 * part2
  if(isFALSE(LOG)) {
    llOUT <- exp(llOUT)
  } else if(!isTRUE(LOG) || isFALSE(LOG)) {
    stop("Wrong argument type for 'LOG': must be logical either TRUE or FALSE.")
  }
  return(llOUT)
}
