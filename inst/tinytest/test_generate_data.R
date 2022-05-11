# The state vector X_t=(s_t^x, u_t^x, s_t^y, u_t^y) contains position and
# <- <- <- <- <- <- <- <- <-  velocity of an object moving in a plane (this is known as the almost constant
# velocity model from the tracking literature). The position, s_t^x and s_t^y,
# can be observed imperfectly (with noise) while the velocities, u_t^x and u_t^y
# are latent.
# A <- diag(4)*2
# C <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 2)
# Q <- diag(c(0.2, 0.001, 0.2, 0.001))
# R <- diag(4) * 10
# TT <- 100
# RcppSMCkalman::dataGenLGSSM(TT, dimX = 4, dimY = 2,
#                             A = A, C = C, Q = Q, R = R)
