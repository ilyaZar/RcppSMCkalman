#### Test dimension missmatches of user argument inputs for system matrices,
#### innovations variances and measurment error variances
### Testing state dimension against missmatches
## State dimension = 1
# missmatch A
A <- diag(4)*0.9
B <- c(3, 4)
C <- c(1, 2)
D <- c(4, 5)
Q <- 0.2
R <- diag(2) * 10
TT <- 100
expect_error(dataGenLGSSM(TT, dimX = 1, dimY = 2, numU = 2,
                          A = A, B = B, C = C, D = D,
                          Q = Q, R = R),
             pattern = "'dimX' dimension does not match 'A' dimension")
A <- 0.9
expect_silent(dataGenLGSSM(TT, dimX = 1, dimY = 2, numU = 2,
                           A = A, B = B, C = C, D = D,
                           Q = Q, R = R))
# missmatch B
A <- 0.9
B <- matrix(3, nrow = 2, ncol = 2) # Wrong B-matrix since dimX = 1
C <- c(1, 2)
D <- c(4, 5)
Q <- 0.2
R <- diag(2) * 10
TT <- 100
expect_error(dataGenLGSSM(TT, dimX = 1, dimY = 2,
                          A = A, B = B, C = C, D = D,
                          Q = Q, R = R),
             pattern = "'dimX' does not match row-dimension of matrix 'B'")
B <- 3
expect_silent(dataGenLGSSM(TT, dimX = 1, dimY = 2,
                           A = A, B = B, C = C, D = D,
                           Q = Q, R = R))
B <- c(3, 1)
expect_error(dataGenLGSSM(TT, dimX = 1, dimY = 2,
                          A = A, B = B, C = C, D = D,
                          Q = Q, R = R),
             pattern = "For dimX=1, 'numU' does not match 'B' length")
B <- c(3, 1)
expect_silent(dataGenLGSSM(TT, dimX = 1, dimY = 2, numU = 2,
                           A = A, B = B, C = C, D = D,
                           Q = Q, R = R))
# missmatch C

# missmatch Q

### State dimension > 1
# missmatch A

# missmatch B
A <- matrix(0.9, nrow = 3, ncol = 3)
B <- matrix(0.3, nrow = 4, ncol = 2) # Wrong B-matrix since dimX = 3
C <- matrix(c(0.1, 0.2), nrow = 2, ncol = 3)
D <- c(0.4, 0.5)
Q <- diag(3) * 0.2
R <- diag(2) * 10
TT <- 100
expect_error(dataGenLGSSM(TT, dimX = 3, dimY = 2, numU = 2,
                          A = A, B = B, C = C, D = D,
                          Q = Q, R = R),
             pattern = "'dimX' does not match row-dimension of matrix 'B'")
B <- matrix(0.3, nrow = 3, ncol = 2) # Correct matrix dimension since dimX=3
expect_silent(dataGenLGSSM(TT, dimX = 3, dimY = 2, numU = 2,
                           A = A, B = B, C = C, D = D,
                           Q = Q, R = R))
B <- matrix(0.3, nrow = 3, ncol = 3)
expect_error(dataGenLGSSM(TT, dimX = 3, dimY = 2, numU = 2,
                          A = A, B = B, C = C, D = D,
                          Q = Q, R = R),
             pattern = "'numU' does not match col-dimension of matrix 'B'")
B <- matrix(0.3, nrow = 3, ncol = 2)
expect_silent(dataGenLGSSM(TT, dimX = 3, dimY = 2, numU = 2,
                           A = A, B = B, C = C, D = D,
                           Q = Q, R = R))
B <- matrix(0.3, nrow = 3, ncol = 1)
expect_silent(dataGenLGSSM(TT, dimX = 3, dimY = 2,
                           A = A, B = B, C = C, D = D,
                           Q = Q, R = R))
# missmatch C

# missmatch Q

### Observation dimension = 1
# missmatch C

# missmatch D

# missmatch R

### Observation dimension > 1
# missmatch C

# missmatch D

# missmatch R

#
# A <- diag(4)*2
# C <- matrix(c(2, 0, 0, 0, 0, 1, 0, 0), nrow = 2)
# Q <- diag(c(0.2, 0.001, 0.2, 0.001))
# R <- diag(4) * 10
# TT <- 100
# RcppSMCkalman::dataGenLGSSM(TT, dimX = 4, dimY = 2,
#                             A = A, C = C, Q = Q, R = R)
