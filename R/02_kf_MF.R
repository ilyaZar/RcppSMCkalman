#' Kalman marginal (forward) filtering computation and prediction
#'
#' Runs (Marginal) filtering and prediction recursions, as described in the
#' corresponding \emph{"(Marginal) Filtering and prediction - implemented
#' recursions"} subsection of the \code{Details} section in
#' \code{\link{kfLGSSM}}.
#'
#' @inheritParams kfLGSSM
#' @param dimX integer giving the dimension of the latent state process
#' @param dimY integer giving the dimension of the measurement process
#' @param TT integer giving the length of the time series
#' @param x00 see \code{initX} from \code{\link{kfLGSSM}}
#' @param P00 see \code{initP} from \code{\link{kfLGSSM}}
#' @param u00 see \code{initU} from \code{\link{kfLGSSM}}
#' @param PDSTORE logical; if \code{TRUE} prediction quantities
#'   \eqn{\hat{P}_{t|t-1}} and \eqn{\hat{x}_{t|t-1}} are stored
#'   \eqn{\forall t=1,\ldots,T}; otherwise, only stored for current
#'   iteration \eqn{t} and not returned.
#'
#' @return if \code{PDSTORE = FALSE} a named list of two
#'   \itemize{
#'   \item \code{mfdEXP} a matrix of dimension \code{dimX x TT} with each
#'      column being the corresponding \eqn{\hat{x}_{t|t}} (see the
#'      \code{Details} section in \code{\link{kfLGSSM}}).
#'   \item \code{mfdVAR} an array of dimension \code{dimX x dimX x TT} with
#'   matrices \eqn{\hat{P}_{t|t}} of dimension \code{dimX x dimX} \eqn{\forall
#'   t = 1,\ldots,TT} (see the \code{Details} section in \code{\link{kfLGSSM}}).
#'   }
#'   if \code{PDSTORE = TRUE} a named list of four:
#'   \itemize{
#'   \item \code{mfdEXP} as above
#'   \item \code{mfdVAR} as above
#'   \item \code{pddEXP} a matrix of dimension \code{dimX x (TT + 1)} with each
#'      column being the corresponding \eqn{\hat{x}_{t+1|t}}, starting from
#'      \eqn{\hat{x}_{1|0}} and running to \eqn{\hat{x}_{T+1|T}} (see the
#'      \code{Details} section in \code{\link{kfLGSSM}}).
#'   \item \code{pddVAR} an array of dimension \code{dimX x dimX x (TT + 1)}
#'   with matrices \eqn{\hat{P}_{t+1|t}} of dimension \code{dimX x dimX},
#'   \eqn{\forall t = 1, \ldots,TT + 1}, starting with \eqn{\hat{P}_{1|0}} and
#'   running to \eqn{\hat{P}_{T+1|T}} (see the \code{Details} section in
#'   \code{\link{kfLGSSM}}).
#'   }
#'
#' @export
kfMFPD <- function(yObs, uReg, wReg,
                   dimX, dimY, TT,
                   x00, P00, u00,
                   A, B, C, D, Q, R,
                   PDSTORE = FALSE) {
  xtt <- matrix(0, nrow = dimX, ncol = TT)
  Ptt <- array(0, dim = c(dimX, dimX, TT),
               dimnames = list(NULL, NULL, as.character(1:TT)))
  if (isTRUE(PDSTORE)) {
    xtt1STORE <- matrix(0, nrow = dimX, ncol = TT + 1)
    Ptt1STORE <- array(0, dim = c(dimX, dimX, TT + 1),
                       dimnames = list(NULL, NULL, as.character(1:(TT + 1))))
  }
  BuRegInit <- computeMatReg(mat = B, reg = u00, dim = dimX, lenT = 1)
  BuReg     <- computeMatReg(mat = B, reg = uReg, dim = dimX, lenT = TT)
  DwReg     <- computeMatReg(mat = D, reg = wReg, dim = dimY, lenT = TT)

  xtt1 <- computeXtt1(A[, ,1], x00, BuRegInit[, 1], dimX)
  Ptt1 <- computePtt1(A[, ,1], P00, Q[, ,1])
  if (isTRUE(PDSTORE)) {
    xtt1STORE[, 1]   <- xtt1
    Ptt1STORE[, , 1] <- Ptt1
  }
  Lt   <- computeLt(C, Ptt1, R)
  Kt   <- computeKt(Ptt1, C)
  for(t in 1:TT) {
    # period t quantities for current iteration
    kGain      <- computekG(yObs[, t], C, xtt1, DwReg[, t])
    xtt[, t]   <- computeXtt(xtt1, Kt, Lt, kGain)
    Ptt[, , t] <- computePtt(Ptt1, Kt, Lt, t)

    xtt1 <- computeXtt1(A[, ,t+1], xtt[, t], BuReg[, t], dimX)
    Ptt1 <- computePtt1(A[, ,t+1], Ptt[, , t], Q[, ,t+1])
    # period t+1 quantities for next iteration
    if (isTRUE(PDSTORE)) {
      xtt1STORE[, t + 1]   <- xtt1
      Ptt1STORE[, , t + 1] <- Ptt1
    }
    Lt   <- computeLt(C, Ptt1, R)
    Kt   <- computeKt(Ptt1, C)
  }
  if (PDSTORE) {
    out <- list(mfdEXP = xtt, mfdVAR = Ptt,
                pddEXP = xtt1STORE, pddVAR = Ptt1STORE)
  } else {
    out <- list(mfdEXP = xtt, mfdVAR = Ptt)
  }
  return(out)
}
