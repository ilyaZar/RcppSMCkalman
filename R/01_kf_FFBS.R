#' Runs a Kalman forward filter and backward smoother
#'
#' @param yt a matrix or vector of measurements (observations)
#' @param zt a matrix or vector of linear covariates of the latent state process
#' @param ut a matrix or vector of linear covariates of the measurement process
#' @param A autoregressive state process matrix (if univariate: autocorrelation)
#' @param B linear covariate effects from \code{zt} on the measurements
#' @param C linear effects from latent state process on the measurements
#' @param D linear covariate effects from \code{ut} on the measurements
#' @param Q VCM of the latent state process
#' @param R VCM of the measurement process variance
#' @param P00 initial XXX matrix
#' @param x00 initial/starting value for latent states (\code{t = 0})
#' @param z00 initial value of regressors (\code{t = 0})
#' @param n_sim number of estimation/simulation runs for forward filtering
#' @param n_sim_jsd number of estimation/simulation runs for backward smoothing
#' @param MFD logical; if \code{TRUE}, then the marginal filtering density is
#'   computed
#' @param JSD logical; if \code{TRUE}, then the joint smoothing density is
#'   computed
#'
#' @return a named list of 4:
#'   \itemize{
#'     \item{kf_marginal_filtering_density: a list of length \code{n_sim} each
#'        containing the marginal filtering density for the whole time series}
#'     \item{kf_VCM_mfd: the intermediate \code{P_{tt}} matrices, i.e. the
#'        marginal filtering VCM estimates, see \code{Details} }
#'     \item{kf_joint_smoothing_density: a list of length \code{n_sim_jsd} each
#'        containing the joint smoothing density for the whole time series}
#'     \item{kf_VCM_jsd: the intermediate \code{J_{tt}} matrices, i.e. the
#'        joint smoothing VCM estimates see \code{Details} }
#'   }
#' @export
kf_lgssm <- function(yt, zt, ut, A, B, C, D, Q, R,
                     P00, x00, z00,
                     n_sim = 1,
                     n_sim_jsd = 1,
                     MFD = TRUE,
                     JSD = TRUE) {
  TT     <- ncol(yt)
  dim_x <- nrow(yt)
  A <- as.matrix(diag(A))
  B <- as.matrix(diag(B))
  C <- as.matrix(diag(C))
  D <- as.matrix(diag(D))

  xtt <- matrix(0, nrow = dim_x, ncol = TT)
  Ptt <- rep(list(list()), times = TT)

  Ptt1    <- tcrossprod(A, tcrossprod(A, P00)) + Q
  Lt      <- solve(tcrossprod(C, tcrossprod(C, Ptt1)) + R)
  Kt      <- Ptt1 %*% t(C) %*% Lt

  Ptt[[1]]  <- Ptt1 - Kt %*% Lt %*% t(Kt)
  if (!matrixcalc::is.positive.definite(Ptt[[1]])) {
    stop(paste0("matrix is no longer p.d. at iteration number: ", 1))
  }
  xtt[, 1]  <- A %*% x00 + B %*% z00 + Kt %*% (yt[, 1] - C %*% (A %*% x00 + B %*% z00) - D %*% ut[, 1])
  for (t in 2:TT) {
    Ptt1     <- tcrossprod(A, tcrossprod(A, Ptt[[t - 1]])) + Q
    Lt       <- solve(tcrossprod(C, tcrossprod(C, Ptt1)) + R)
    Kt       <- Ptt1 %*% t(C) %*% Lt

    Ptt[[t]] <- Ptt1 - Kt %*% Lt %*% t(Kt)
    if (!matrixcalc::is.positive.definite(Ptt[[t]])) {
      stop(paste0("matrix is no longer p.d. at iteration number: ", t))
    }
    xtt[, t] <- A %*% xtt[, t - 1] + B %*% zt[, t - 1] + Kt %*% (yt[, t] - C %*% (A %*% xtt[, t - 1] + B %*% zt[, t - 1]) - D %*% ut[, t])
  }
  if (MFD) {
    kf_MFD_res <- rep(list(list()), times = n_sim)
    for (n in 1:n_sim) {
      kf_MFD <- matrix(0, ncol = TT, nrow = dim_x)
      for (t in 1:TT) {
        kf_MFD[, t] <- rmvnorm(1, mean = xtt[, t], sigma = Ptt[[t]])
      }
      kf_MFD_res[[n]] <- kf_MFD
    }
  }
  if (JSD) {
    kf_JSD_all <- matrix(0, nrow = dim_x*n_sim_jsd, ncol = TT)
    for (n in 1:n_sim_jsd) {
      Jtt <- rep(list(list()), times = TT)
      xtt_jsd <- matrix(0, nrow = dim_x, ncol = TT)
      kf_JSD  <- matrix(0, nrow = dim_x, ncol = TT)

      kf_JSD[, TT] <- rmvnorm(1, mean = xtt[, TT], sigma = Ptt[[TT]])

      for (t in (TT - 1):1) {
        Jtt[[t]] <- Ptt[[t]] %*% t(A) %*% solve(tcrossprod(A, tcrossprod(A, Ptt[[t]])) + Q)
        xtt_jsd[, t]  <- xtt[, t] + Jtt[[t]]  %*% (kf_JSD[, t + 1] - B %*% zt[, t] - A %*% xtt[, t])
        kf_JSD[, t] <- rmvnorm(1, mean = xtt_jsd[, t], sigma = Jtt[[t]])
      }
      # browser()
      print(n)
      row_ID <- (dim_x*(n - 1) + 1):(dim_x*n)
      kf_JSD_all[row_ID, ] <-  kf_JSD
    }
  }
  return(list(kf_marginal_filtering_density = kf_MFD_res, kf_VCM_mfd = Ptt,
              kf_joint_smoothing_density = kf_JSD_all, kf_VCM_jsd = Jtt))
}
