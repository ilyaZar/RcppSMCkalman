#' Kalman marginal smoothing
#'
#' Performs marginal smoothiong i.e. computation of means and variances for
#' \eqn{p(x_t|y_{1:T})\forall t=1,\ldots,TT}, as described in the corresponding
#' \emph{"Marginal smoothing - implemented recursions"} subsection of the
#' \code{Details} section in \code{\link{kfLGSSM}}.
#'
#' @inheritParams kfMFPD
#' @param xtt forward filtering means as produced by \code{\link{kfMFPD}}
#' @param Ptt forward filtering variances as produced by \code{\link{kfMFPD}}
#' @param xtt1 predictive means as produced by \code{\link{kfMFPD}} (if
#' \code{PDSTORE = TRUE})
#' @param Ptt1 predictive variances as produced by \code{\link{kfMFPD}} (if
#' \code{PDSTORE = TRUE})
#'
#' @return a named list of two:
#'   \itemize{
#'   \item\code{msdEXP:} predicitve means \eqn{x_{t+1|t}}
#'      (see \emph{"Marginal smoothing - implemented recursions"} in
#'      \code{\link{kfLGSSM}})
#'   \item\code{msdEXP:} predicitve variances \eqn{P_{t+1|t}}
#'      (see \emph{"Marginal smoothing - implemented recursions"} in
#'      \code{\link{kfLGSSM}})
#'   }
#'
#' @export
kfMSD <- function(dimX, TT,
                  xtt, Ptt, xtt1, Ptt1,
                  A, Q) {
  # 0. Housekeeping of output containers
  xtT <- matrix(0, nrow = dimX, ncol = TT)
  PtT <- array(0, dim = c(dimX, dimX, TT),
               dimnames = list(NULL, NULL, as.character(1:TT)))
  Jtt <- array(0, dim = c(dimX, dimX, TT),
               dimnames = list(NULL, NULL, as.character(1:TT)))
  # 1. Initialization: t = T
  xtT[, TT]   <- xtt[, TT]
  PtT[, , TT] <- Ptt[, , TT]
  Jtt[, , TT] <- matrix(0, dimX, dimX)
  # 2. Iteration: t = T-1,...,1
  for (t in (TT - 1):1) {
    Jtt[, , t] <- tcrossprod(Ptt[, , t], A) %*% solve(PtT[, , t + 1])
    matGain    <- PtT[, , t + 1] - Ptt1[, , t + 1]
    PtT[, , t] <- Ptt[, , t] + Jtt[, , t] %*% tcrossprod(matGain, Jtt[, , t])
    xtT[, t]   <- xtt[, t] + Jtt[, , t] %*% xtT[, t + 1] - xtt1[, t + 1]
  }
  return(list(msdEXP = xtT, msdVAR = PtT))
}
