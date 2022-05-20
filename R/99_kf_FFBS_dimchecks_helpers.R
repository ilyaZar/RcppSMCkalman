getDimX <- function(initX, A, B, Q) {
  if(!is.null(initX)) {
    dimX <- length(initX)
  } else if(!is.null(A)) {
    dimX <- ifelse(is.matrix(A), nrow(A), 1)
  } else if(!is.null(Q)) {
    dimX <- ifelse(is.matrix(Q), nrow(Q), 1)
  } else if(!is.null(B)) {
    dimX <- ifelse(is.matrix(B), nrow(B), 1)
  } else {
    stop("Model miss-specified: can not infer state process dimension.")
  }
  return(dimX)
}
getNumReg <- function(reg) {
  ifelse(test = is.matrix(reg),
         yes  = nrow(reg),
         no   = ifelse(is.null(reg), 0, 1))
}
checkArgumentInputs <- function(TT, dimY, dimX, numU, numW,
                                yObs, uReg, wReg,
                                A, B, C, D, Q, R,
                                initX, initP, initU) {
  ### check matching dimension/length for T
  ## skip check for yObs since TT is inferred from measurement length/dim
  ## check uReg series length
  checkDim(TT, uReg, "col", "T", "U-type regressors")
  ## check wReg series length
  checkDim(TT, wReg, "col", "T", "W-type regressors")
  ### check matching dimension/length for dimY
  ## skip check for yObs since dimY is inferred from measurement length/dim
  checkDim(dimY, wReg, "row", "dimY", "W-type regressors")
  checkDim(dimY, C, "row", "dimY", "'C' matrix")
  checkDim(dimY, D, "row", "dimY", "'D' matrix")
  checkDim(dimY, R, "row", "dimY", "'R' matrix")
  checkSym(R, "'R'")
  ### check matching dimension/length for dimX
  ## skip check for initX since dimX is inferred from state length/dim
  checkDim(dimX, uReg, "row", "dimX", "U-type regressors")
  checkDim(dimX, A, "row", "dimX", "'A' matrix")
  checkSym(A, "'A'")
  checkDim(dimX, B, "col", "dimX", "'B' matrix")
  checkDim(dimX, Q, "col", "dimX", "'Q' matrix")
  checkSym(Q, "'Q'")
  checkDim(dimX, initP, "col", "dimX", "'initP' matrix")
  ### check matching dimension/length for numU
  checkDim(numU, uReg, "row", "numU", "U-type regressors")
  checkDim(numU, B, "col", "numU", "'B' matrix")
  msg <- paste("Initial values for U-type regressors 'initU' do not match",
               "number of regressors derived from matrix B")
  if(length(initU) != numU) stop(msg)
  ### check matching dimension/length for numW
  checkDim(numW, wReg, "row", "numW", "W-type regressors")
  checkDim(numW, D, "col", "numW", "'D' matrix")
}
checkDim <- function(dim, mat, type = "row", nameDim, nameMat) {
  if (!is.null(mat)) {
    msg <- paste("Wrong dimension for", nameMat,
                 ": does not match", nameDim, ".")
    if (type == "row") {
      checkMe <- ifelse(is.matrix(mat), nrow(mat), 1)
    } else if (type == "col") {
      checkMe <- ifelse(is.matrix(mat), ncol(mat), 1)
    } else {
      stop("Unknown check type: use only 'col' or 'row'.")
    }
    if (checkMe != dim) stop(msg)
  }
}
checkSym <- function(mat, matName) {
  msg <- paste(matName, "matrix must be symmetric.")
  checkMe <- ((length(mat) %in% c(0, 1)) || nrow(mat) == ncol(mat))
  if(!checkMe) stop(msg)
}
checkIsMatrix <- function(listOfMatrices) {
  numMat <- length(listOfMatrices)
  namMat <- names(listOfMatrices)
  for (i in 1:numMat) {
    if(isTRUE(is.null(listOfMatrices[[i]]))) next
    if (isTRUE(is.null(dim(listOfMatrices[[i]])))) {
      if (length(listOfMatrices[[i]]) ==  1L) {
        msg <- paste0("No 0 scalar values permitted: leave arg '", namMat[[i]],
                      "' as default NULL,",
                      " if systemt matrix should be dropped from the model.")
        if(listOfMatrices[[i]] == 0) stop(msg)
        next
      } else {
        msg <- paste0("Matrix '", namMat[i],
                      "' appears to be passed as vector argument (convert to
                      column or row matrix).")
        stop(msg)
      }
    }
    msg <- paste0("Matrix '", namMat[i],
                  "' is neither a matrix nor default NULL.")
    if(isFALSE(is.matrix(listOfMatrices[[i]]))) stop(msg)
  }
}
setDimCase <- function(dimX, dimY, silent = FALSE) {
  if (dimX == 1 && dimY == 1) dimCase <- c("dimX=1 and dimY=1" = 1)
  if (dimX  > 1 && dimY == 1) dimCase <- c("dimX>1 and dimY=1" = 2)
  if (dimX == 1 && dimY >  1) dimCase <- c("dimX=1 and dimY>1" = 3)
  if (dimX >  1 && dimY >  1) dimCase <- c("dimX>1 and dimY>1" = 4)
  if(isTRUE(silent)) return(dimCase)
  cat(paste0(crayon::green("The following model dimensions are inferred: \n"),
             crayon::red(paste0(names(dimCase), ".\n"))))
  return(dimCase)
}
setCmpCase <- function(computeMFD, computePDD,
                       computeMSD, computeJSD,
                       computeLLH, silent = FALSE) {
  updateCases(computeMFD, computePDD,
              computeMSD, computeJSD,
              computeLLH, silent = FALSE)
  if (isTRUE(computeMFD)){
    cmpCase1 <- c("Marginal filtering density computations" = 1)
  } else {
    cmpCase1 <- NULL
  }
  if (isTRUE(computePDD)) {
    cmpCase2 <- c("Marginal filtering and predictive density computations" = 2)
  } else {
    cmpCase2 <- NULL
  }
  if (isTRUE(computeMSD)) {
    cmpCase3 <- c("Marginal smoothing density computations" = 3)
  } else {
    cmpCase3 <- NULL
  }
  if (isTRUE(computeJSD)) {
    cmpCase4 <- c("Joint smoothing density computations" = 4)
  } else {
    cmpCase4 <- NULL
  }
  if (isTRUE(computeLLH)) {
    cmpCase5 <- c("Observed data likelihood computations" = 5)
  } else {
    cmpCase5 <- NULL
  }
  cmpCase <- c(cmpCase1, cmpCase2, cmpCase3, cmpCase4, cmpCase5)
  cmpCasePrint <- paste0(cmpCase, ". ", names(cmpCase), "\n")
  if(isTRUE(silent)) return(cmpCase)
  cat(crayon::green("The following computational tasks are inferred: \n"))
  cat(crayon::cyan(paste(cmpCasePrint, collapse =  "")))
  return(cmpCase)
}
updateCases <- function(computeMFD, computePDD,
                        computeMSD, computeJSD,
                        computeLLH, silent = FALSE) {
  if (isTRUE(computeMSD)) {
    if (isFALSE(computePDD)) {
      msg <- paste0("Setting 'computePDD=TRUE' internally.\n",
                    "For marginal smoothing density computations, ",
                    "the marginal filtering and prediction densities are ",
                    "required and computed.\n",
                    "Results are also added to the function output!")
      warning(msg)
      assign("computePDD", TRUE, envir = parent.frame())
    }
    assign("computeMFD", FALSE, envir = parent.frame())
  }

  if (isTRUE(computeJSD)) {
    if (isFALSE(computePDD)) {
      msg <- paste0("Setting 'computePDD=TRUE' internally. ",
                    "For joint smoothing density computations, ",
                    "the marginal filtering density is ",
                    "required and computed.\n",
                    "Results are also added to the function output!")
      warning(msg)
      assign("computeMFD", TRUE, envir = parent.frame())
    }
  }

  if (isTRUE(computeLLH)) {
    if (isFALSE(computePDD)) {
      msg <- paste0("Setting 'computePDD=TRUE' internally. ",
                    "For computing the observed data likelihood, ",
                    "the prediction density is required and computed.\n",
                    "Results are also added to the function output.")
      warning(msg)
      assign("computePDD", TRUE, envir = parent.frame())
    }
    assign("computeMFD", FALSE, envir = parent.frame())
  }
}
