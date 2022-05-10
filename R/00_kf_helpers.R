#' Generates simulated data from a multivariate LGSSM.
#'
#' @param TT integer giving the time series length
#' @param dim_x integer giving the state process dimension; defaults to
#'   \code{dim_x = 1}
#' @param dim_y integer giving the measurement process dimension; defaults to
#'   \code{dim_x = 1}
#' @param A state process matrix1
#' @param B state process matrix2
#' @param C measurement process matrix1
#' @param D measurement process matrix2
#' @param Q VCM state process
#' @param R VCM measurement process
#' @param initX Initial value of state process (must have same length as
#'   \code{dim_x}); can be thought of X_{t=0}, where time series runs from
#'   \code{t=1,...,TT}
#'
#' @return a named list of 7:
#'   \itemize{
#'     \item matrix or vector of simulated state process: \code{dim_x x TT}
#'     \item matrix or vector of simulated measurements: \code{dim_y x TT}
#'     \item matrix or vector of simulated regressors for states:
#'       \code{ncol(B) x TT}
#'     \item matrix or vector of simulated regressors for measurements:
#'       \code{ncol(D) x TT}
#'     \item initial state regressor values: \code{ncol(B) x 1}
#'     \item initial measurement regressor values: \code{ncol(D) x 1}
#'     \item \code{initX}, see argument to function call
#'   }
#' @export
data_gen_lgssm <- function(TT, dim_x = 1, dim_y = 1,
                           A, B, C, D, Q, R,
                           initX) {

  stopifnot(`Initial state process unequal to dim_x` = length(initX) == dim_x)
  x0 <- initX
  if (dim_x == 1) {
    xt <- numeric(TT)
    yt <- numeric(TT)
    zt <- rnorm(1, 0, 1)
    ut <- rnorm(1, 0, 1)

    xt[1] <- A*x0 + B*zt[1] + rnorm(1, 0, sd = Q^2)
    for (t in 2:TT) {
      xt[t] <- A*xt[t - 1] + B*zt[t] + rnorm(1, 0, sd = Q^2)
    }
    yt <- C*xt + D*ut + rnorm(TT, 0, sd = R)
    return(list(xt, yt))
  } else if (dim_x >= 2) {
    xt <- matrix(0, nrow = dim_x, ncol = TT)
    yt <- matrix(0, nrow = dim_y, ncol = TT)
    zt <- matrix(rnorm(dim_x*TT, 0, 1), nrow = dim_x, ncol = TT)
    ut <- matrix(rnorm(dim_y*TT, 0, 1), nrow = dim_y, ncol = TT)
    z0 <- rnorm(dim_x, 0, 1)

    xt[, 1] <- A*x0 + B*z0 + rmvnorm(n = 1, mean = rep(0, times = dim_x), sigma = Q)
    for (t in 2:TT) {
      xt[, t] <- A*xt[, t - 1] + B*zt[, t] + rmvnorm(n = 1, mean = rep(0, times = dim_x), sigma = Q)
      yt[, t] <- C*xt[, t] + D*ut[, t] +  rmvnorm(n = 1, mean = rep(0, times = dim_y), sigma = R)
    }
    return(list(states = xt,
                measurements = yt,
                x_reg = zt, y_reg = ut,
                z_init = z0,
                x_init = x0))
  }
}

# kflg <- function(yt, A, C, Q, R, P00, x00, n_sim = 1) {
#   TT    <- length(yt)
#   xtt  <- numeric(TT)
#   Ptt  <- numeric(TT)
#   # Ptt1    <- A %*% P00 %*% A  + Q
#   Ptt1    <- tcrossprod(A, tcrossprod(A, P00))
#   Lt <-
#   Kt      <- Ptt1*C*(((C^2)*Ptt1 + R)^(-1))
#   Ptt[1]  <- Ptt1 - Kt*C*Ptt1
#   xtt[1]  <- A*x00 + Kt*(yt[1] - C*A*x00)
#   for (t in 2:TT) {
#     Ptt1     <- (A^2)*Ptt[t - 1] + Q
#     Kt     <- Ptt1*C*(((C^2)*Pttm + R)^(-1))
#     Ptt[t] <- Ptt1 - Kt*C*Ptt1
#     xtt[t] <- A*xtt[t - 1] + Kt*(yt[t] - C*A*xtt[t - 1])
#   }
#   # KFapprox <- xtt
#   res_sim <- matrix(0, ncol = TT, nrow = n_sim)
#   for (n in 1:n_sim) {
#     res_sim[n, ] <- rnorm(TT, mean = xtt, sd = sqrt(Ptt))
#   }
#   if (n_sim == 1) {
#     res_sim <- as.vector(res_sim)
#   }
#   return(list(KFx = res_sim, KFvar = Ptt))
# }
