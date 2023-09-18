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
#' @param LOG logical; if \code{TRUE}, then the logarithm of the likelihood is
#'   returned
#'
#' @return the (logarithmic, if \code{LOG=TRUE}) value of the data likelihood
#' @export
kfLLH <- function(yObs, wReg, xtt1, Ptt1, C, D, R, dimX, dimY, TT, LOG = TRUE) {
  DwReg <- computeMatReg(mat = D, reg = wReg, dim = dimY, lenT = TT)
  part2 <- 0

  AllObsMissing <- FALSE
  MissingsObsTotal <- sum(is.na(yObs))

  for (t in 1:TT) {
    CAdj <- C
    RAdj <- R[, , t]
    yObsAdj <- yObs[, t]
    DwRegAdj <- DwReg[, t]
    MissingObs <- which(is.na(yObsAdj))

    if(length(MissingObs) != 0) {
      if (length(MissingObs) == dimY) {
        AllObsMissing <- TRUE
      } else {
        W <- diag(dimY)
        W <- W[-MissingObs, ]

        yObsAdj[MissingObs] <- 0
        CAdj[MissingObs, ] <- 0
        RAdj[MissingObs, ] <- 0
        RAdj[ ,MissingObs] <- 0

        CAdj <- W %*% CAdj
        RAdj <- W %*% RAdj %*% t(W)
        yObsAdj <- W %*% yObsAdj
        DwRegAdj <- W %*% DwRegAdj
      }
    }

    if(!AllObsMissing) {
      meanY      <- CAdj %*% xtt1[, t] + DwRegAdj
      VarY       <- CAdj %*% tcrossprod(Ptt1[, , t], CAdj) + RAdj
      logDetVarY <- determinant(VarY, logarithm = TRUE)$modulus[1]

      part2 <- part2 + logDetVarY
      part2 <- part2 + t(yObsAdj - meanY) %*% solve(VarY) %*% (yObsAdj - meanY)
    }
  }

  part1 <- -(TT*dimY - MissingsObsTotal) /2 * log(2*pi)

  llOUT <- part1 - 0.5 * part2
  if(isFALSE(LOG)) {
    llOUT <- exp(llOUT)
  } else if(!isTRUE(LOG) || isFALSE(LOG)) {
    stop("Wrong argument type for 'LOG': must be logical either TRUE or FALSE.")
  }
  return(llOUT)
}
