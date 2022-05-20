#' Kalman joint smoothing
#'
#' Performs joint smoothiong via backward simulation i.e. generating
#' \code{nSimJS} realizations from \eqn{p(x_{1:T}|y_{1:T})}, as described in
#' the corresponding \emph{"Joint smoothing density- implemented backward
#' simulation recursions"} subsection of the \code{Details} section in
#' \code{\link{kfLGSSM}}.
#'
#' @inheritParams kfMFPD
#' @inheritParams kfMSD
#' @inheritParams kfLGSSM
#'
#' @return a named list of three:
#'   \itemize{
#'   \item\code{jsdEXP:} an array of dimension \code{dimX x TT x nSimJS} giving
#'      the smoothing means \eqn{\mu_t} of dimension \code{dimX} for all
#'      \eqn{t=TT,\ldots,1} per \eqn{1,\ldots,}\code{nSimJS} simulation run
#'   \item\code{jsdVAR:} an array of dimension \code{dimX x dimX  TT} giving
#'      the smoothing variances \eqn{L_t} of dimension \code{dimX x dimX} for
#'      all \eqn{t=TT,\ldots,1} (note: these do not change with the number of
#'      simulation runs \code{nSimJS}, in contrast to previous means)
#'   \item\code{jsdTRJ:} an array of dimension \code{dimX x TT x nSimJS} giving
#'      the sampled trajectories \eqn{\tilde{x}_{1:T}} (\code{dimX x TT}) from
#'      \eqn{p(x_{1:T}|y_{1:T})} for each simulation \code{nSimJS}
#'   }
#' @export
kfJSD <- function(uReg, xtt, Ptt, A, B, Q, dimX, TT, nSimJS = 1) {
  # 0. Housekeeping of output containers
  mut <- array(0, dim = c(dimX, TT, nSimJS),
               dimnames = list(NULL,
                               paste0("t_", as.character(1:TT)),
                               paste0("sim_", as.character(1:nSimJS))))
  Lt2 <- array(0, dim = c(dimX, dimX, TT),
               dimnames = list(NULL, NULL,
                               paste0("t_", as.character(1:TT))))
  xJSD <- array(0, dim = c(dimX, TT, nSimJS),
                dimnames = list(NULL,
                                paste0("t_", as.character(1:TT)),
                                paste0("sim_", as.character(1:nSimJS))))
  Jtt2 <- array(0, dim = c(dimX, dimX, TT),
                dimnames = list(NULL, NULL,
                                paste0("t_", as.character(1:TT))))
  # 1.A Initialization: t = T
  for (i in 1:nSimJS) {
    xJSD[, TT, i] <- mvtnorm::rmvnorm(n = 1, mean = xtt[, TT],
                                      sigma = matrix(Ptt[, , TT]))
  }
  # 1.B Pre-compute necessary quantities
  BuReg <- computeMatReg(B, uReg, dimX, TT)
  # 2. Iteration
  for (t in (TT - 1):1) {
    Jtt2[, , t] <- computeJtt2(Ptt[, , t], A, Q)
    Lt2[, , t]  <- Ptt[, , t] - Jtt2[, , t] %*% A %*% Ptt[, , t]
    for (i in 1:nSimJS) {
      jsdG          <- xJSD[, t + 1, i] - BuReg[, t] - A %*% xtt[, t]
      mut[, t, i]   <- xtt[, t] + Jtt2[, , t] %*% jsdG
      xJSD[, t, i] <- mvtnorm::rmvnorm(n = 1, mean = mut[, t, i],
                                       sigma = matrix(Lt2[, , t]))
    }
  }
  return(list(jsdEXP = mut,
              jsdVAR = Lt2,
              jsdTRJ = xJSD))
}
