# This file contains developer level functions mainly adapted from SPM12 for use in fitting image models

# To Do:
# check: .c2tsp, maybe .res
# I can probably use 'identical()' instead of these other ones
# .C1inC2
# .notunique (~unique)

## from spm_sp.m----
# set the filtered and whitened design matrix
# the `nu` and `nv` arguments seem to make the difference between the matlab svd and R svd
.setx <- function(x) {
  out <- list()
  if (missing(x)) {
    out <- list(X = c(), v = c(), u = c(), d = c(), tol = c(), rk = c())
  } else {
    out$X <- x
    if (nrow(x) < ncol(x))
      svdx <- svd(t(x), nu = nrow(x))
    else
      svdx <- svd(x, nu = nrow(x))
    
    out$v <- svdx$v
    out$u <- svdx$u
    out$d <- svdx$d
    out$tol <- max(dim(out$X)) * max(abs(out$d)) * .Machine$double.eps
    out$rk <- sum(out$d > out$tol)
  }
  return(out)
}

# orthogonal projectors on space of X
.op <- function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- tcrossprod(x$u[, seq_len(x$rk)], x$u[, seq_len(x$rk)])
  else  
    out <- matrix(0, nrow(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out < x$tol)] <- 0
  return(out)
}

# orthogonal projectors on space of t(X)
.opp <- function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- crossprod(x$v[, seq_len(x$rk)], x$v[, seq_len(x$rk)])
  else
    out <- matrix(0, ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out < x$tol)] <- 0
  return(out)
}

# pseudo inverse of X
.pinvx <- function(x, y, check = FALSE) {
  if (x$rk > 0 ) {
    out <- x$v[, seq_len(x$rk)] %*% 
      tcrossprod(diag(rep(1, x$rk) / x$d[seq_len(x$rk)]),
                 x$u[, seq_len(x$rk)])
  } else
    out <- matrix(0, ncol(x$X), nrow(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

# pseudo-inverse of t(X)
.pinvxp <- function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- x$u[, seq_len(x$rk)] %*% 
      tcrossprod(diag(rep(1, x$rk) / x$d[seq_len(x$rk)]),
                 x$v[, seq_len(x$rk)])
  else
    out <- matrix(0, nrow(x$X), ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out < x$tol)] <- 0
  return(out)
}

# coordinates of pseudo-inverse of t(X) in the base of uk
.cukxp <- function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- tcrossprod(diag(rep(1, x$rk) / x$d[seq_len(x$rk)]), x$v[, seq_len(x$rk)])
  else
    out <- rep(0, ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

# coordinates of X in the base of uk
.cukx <- function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- tcrossprod(diag(x$d[seq_len(x$rk)]), x$v[, seq_len(x$rk)])
  else
    out <- rep(0, ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[out < x$tol] <- 0
  return(out)
}

# orthonormal basis sets for X
.ox <- function(x, y, check = FALSE) {
  if (x$rk > 0) {
    out <- x$u[,seq_len(x$rk)]
    if (!missing(y))
      out <- out %*% y
    if (check)
      out[abs(out) < x$tol] <- 0
  } else
    out <- matrix(0, nrow(x$X), 1)
  return(out)
}

# orthonormal basis sets for t(X)
.oxp <- function(x, y, check = FALSE) {
  if (x$rk > 0)
    x$v[, seq_len(x$rk)]
  else
    matrix(0, ncol(x$X), 1)
}

# pseudo-inverse of crossprod(X)'
.pinvxpx <- function(x, y, check = FALSE) {
  if (x$rk > 0)
    out <- x$v[, seq_len(x$rk)] %*%
      tcrossprod(diag(rep(1, x$rk) / x$d[seq_len(x$rk)]), x$v[, seq_len(x$rk)])
  else
    out <- rep(0, ncol(x$X))
  if (!missing(y))
    out %*% y
  if (check)
    out[abs(out < x$tol)] <- 0
  return(out)
}

# computation of crossprod(X)
.xpx <- function(x, y, check = FALSE) {
  if (x$rk > 0) {
    out <- x$v[, seq_len(x$rk)] %*%
      tcrossprod(diag(x$d[seq_len(x$rk)]^2),
                 x$v[, seq_len(x$rk)])
  } else
    out <- rep(0, ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

# coordinates in the basis of X to basis u
.cx2cu <- function(x, y, check = FALSE) {
  if (mussing(y))
    out <- tcrossprod(diag(x$d), x$v)
  else
    out <- tcrossprod(diag(x$d), x$v) %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

# pseudo-inverse of tcrossprod(X)'
.pinvxxp <- function(x, y, check = FALSE) {
  if (x$rk > 0) {
    out <- x$u[, seq_len(x$rk)] %*%
      tcrossprod(diag(x$d[seq_len(x$rk)]^(-2)),
                 x$u[, seq_len(x$rk)])
  } else
    out <- rep(0, nrow(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

# computation of tcrossprod(X)
.xxp <- function(x, y, check = FALSE) {
  if (x$rk > 0) {
    if(missing(y)) {
      out <- x$u[, seq_len(x$rk)] %*%
        tcrossprod(diag(x$d[seq_len(x$rk)]^2),
                   x$u[, seq_len(x$rk)])
    } else
      x$u[, seq_len(x$rk)] %*% diag(x$d^2) %*%
      crossprod(x$u[, seq_len(x$rk)], y)
  } else
    out <-rep(0, nrow(x$X))
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

# computes u * (diag(x^n)) * t(u)
.power <- function(x, n, y, check = FALSE) {
  if (missing(n))
    n <- 1
  if (x$rk > 0) {
    if (missing(y))
      out <- x$u[, seq_len(x$rk)] %*% tcrossprod(diag(x$d[seq_len(x$rk)] ^ n), x$u[, seq_len(x$rk)])
    else
      out <- x$u[, seq_len(x$fk)] %*% diag(x$d[seq_len(x$fk)] ^ n) %*% crossprod(x$u[, seq_len(x$fk)], y)
  } else {
    if (missing(y))
      out <- rep(0, nrow(x$X))
    else
      matrix(0, nrow(x$X), ncol(y))
  }
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out) < x$tol] <- 0
  return(out)
}

# computes v * (diag(x^n)) * t(v)
.powerp <- function(x, n, y, check = FALSE) {
  if (missing(n))
    n <- 1
  if (x$rk > 0)
    out <- x$v[, seq_len(x$rk)] %*% tcrossprod(diag(x$d[seq_len(x$rk)]^(n)), x$v[, seq_len(x$rk)])
  else
    out <- rep(0, ncol(x$X))
  if (!missing(y))
    out <- out %*% y
  if (check)
    out[abs(out < x$tol)] <- 0
  
  return(out)
}

# null space of t(X)
.np <- function(x) {
  out  <- list()
  out$X <- t(x$X)
  out$v <- x$u
  out$u <- x$v
  
  out$oP <- x$opp
  out$oPp <- x$op
  
  dimX <- dim(out$X)
  if (x$rk > 0) {
    if (dimX[1] >= dimX[2]) {
      if (out$rk == dimX[2])
        n <- matrix(0, dimX[2], 1)
      else
        n <- out$v[, (x$rk):dimX[2]]
    } else
      n <- MASS::Null(out$X)
  } else
    n <- diag(dimX[2])
  return(n)
}

# project into null space X
.nop <- function(x, y, check = FALSE) {
  dimX <- dim(x$X)
  if (x$rk > 0) {
    if (dimX[1] >= dimX[2]) {
      if (x$rk == dimX[2])
        n <- matrix(0, dimX[2], 1)
      else
        n <- x$v[, (x$rk):dimX[2]]
    } else
      n <- MASS::Null(x$X)
  } else
    n <- diag(dimX[2])
  if (missing(y))
    out <- tcrossprod(n)
  else
    out <- n %*% crossprod(n, y)
  if (check)
    out[abs(out < x$tol)] <- 0
  return(out)
}

# projection into null space of t(X)'
.nopp <- function(x, y, check = FALSE) {
  out <- list()
  out$X <- t(x$X)
  out$v <- x$u
  out$u <- x$v
  out$tol <- x$tol
  out$rk <- x$rk
  if (missing(y))
    .nop(out, check = check)
  else
    .nop(out, y, check = check)
}

# return residual forming matrix or set residuals
.res <- function(x, y, check = FALSE) {
  if (missing(y)) {
    out <- diag(nrow(x$X)) - .op(x)
  } else {
    # sf_rY(x, y)
    dimx <- dim(x$X)
    if (x$rk > 0) {
      if (x$rk < (dimx[1] - x$rk))
        out <- y - x$u[, seq_len(x$rk)] %*% crossprod(x$u[, seq_len(x$rk)], y)
      else {
        if (ncol(y) < (5 * dimx[1]))
          out <- (diag(nrow(x$X)) - .op(x)) %*% y
        else {
          n <- .np(x)
          out <- n %*% crossprod(n, y)
        }
      }
    }
  }
  if (check)
    out[out < x$tol] <- 0
  return(out)
}

# check whether vectors are in row/column space of x
.isinsp <- function(x, c, tol) {
  'is in space or is in dual space'
  if (missing(tol))
    tol <- x$tol
  if (nrow(x$X) != nrow(c)) {
    warning('Vector dim dont match col dim : not in space')
    return(0)
  } else {
    out <- all(abs(.op(x) %*% c - c) <= tol)
    return(colSums(out) > 0)
  }
}

# check whether vectors are in row/column space of x
.isinspp <- function(x, c, tol) {
  'is in space or is in dual space'
  if (missing(tol))
    tol <- x$tol
  if (ncol(x$X) != nrow(c)) {
    warning('Vector dim dont match row dim : not in space')
    return(0)
  } else {
    out <- all(abs(.opp(x) %*% c - c) <= tol)
    return(sum(out) > 0)
  }
}

# each column of c in space or in dual space
.eachinsp <- function(x, c, tol) {
  if (missing(tol))
    tol <- x$tol
  if (nrow(x$X) != nrow(c)) {
    warning('Vector dim dont match col. dim : not in space')
    return(0)
  } else
    return(all(abs(.op(x) %*% c - c) <= tol))
}

# each column of c in space or in dual space
.eachinspp <- function(x, c, tol) {
  if (missing(tol))
    tol <- x$tol
  if (ncol(x$X) != nrow(c)) {
    warning('Vector dim dont match row. dim : not in space')
    return(0)
  } else
    return(all(abs(.opp(x) %*% c - c) <= tol))
}

# test wether two spaces are the same
.equal <- function(x, X2) {
  x2 <- .setx(X2)
  maxtol <- max(x$tol, x2$tol)
  return(all(.isinsp, x, X2, maxtol) && all(.isinsp(x2, x$X, maxtol)))
}

# space structure check
.isspc <- function(x, oneout = TRUE) {
  if (!is.list(x))
    out <- 0
  else {
    b <- 1
    fnames <- names(x)
    tmp <- .setx()
    for (i in fnames) {
      b <- b & any(names(tmp) == i)
      if (!b)
        break
    }
    
    if (!oneout) {
      if (b)
        out <- list(b, "ok")
      else
        out <- list(b, "not a space")
    } else
      out <- b
    return(out)
  }
}

# is it a completed space structure?
.issetspc <- function(x) {
  !(is.null(x$X) | is.null(x$u) | is.null(x$v) | is.null(x$ds) |
      is.null(x$tol) | is.null(x$rk))
}

## from spm_SpUtil.m----
#test whether weight vectors specify contrast
.iscon <- function(x, c) {
  if (!.isspc(x))
    x <- .setx(X)
  .eachinspp(x, c)
}

.allcon <- function() {
  if (!.isspc(x))
    x <- .setx(X)
  .isinspp(x, c)
}

.conR <- function(x, c) {
  if (!.isspc(x))
    x <- .setx(X)
  if (nrow(c) != ncol(x$X))
    stop("The contrast size doesn't match the number of columns in the design matrix.")
  r <- crossprod(c, .xpx(x)) %*% c
  r <- r / crossprod(matrix(sqrt(diag(r))), matrix(sqrt(diag(r))))
  r[abs(r) < x$tol] <- 0
  return(r)
}

.conO <- function(x, c) {
  if (!.isspc(x))
    x <- .setx(X)
  return(abs(crossprod(c, .xpx(x)) %*% c) < x$tol)
}

# check iX0
.iX0check <- function(i0, sL) {
  if (all(any(i0 == 1) && any(i0 == 0)) && (length(i0) == sL))
    i0c <- matrix(seq_len(sL)[i0 != 0])
  else if (all(dim(i0) > 0) && any(floor(i0) != i0) || any(i0 < 1) || any(i0 > sL))
    stop('logical mask or vector of column indices required')
  else
    i0c <- matrix(i0)
  return(i0c)
}

# get the estimable parts of C0 and C1'
.i02c <- function(x, i0) {
  # argument checks
  if (!.isspc(x))
    x <- .setx(x)
  if (x$rk == 0)
    stop("i02c null rank x == 0")
  sL <- ncol(x$X)
  i0 <- .iX0check(i0, sL)
  
  # function
  c0 <- diag(sL)
  c0 <- matrix(c0[, i0])
  c1 <- diag(sL)
  c1 <- matrix(c1[, sort(setdiff(seq_len(sL), i0))])
  
  if (!.isinspp(x, c0))
    c0 <- .opp(x, c0)
  if (!.isinspp(x, c1))
    c1 <- .opp(x, c1)
  
  if (!(length(c1) == 0)) {
    if (!(length(c0 == 0)))
      out <- .r(.setx(c0), c1, check = TRUE)
    else
      out <- .xpx(x)
  } else
    out <- c()
  return(out)
}

# orthonormal partitioning implied by F-contrast
.c2tsp <- function(x, c, oneout = FALSE, plus = FALSE) {
  # check arguments
  if (!.isspc(x))
    x <- .setx(x)
  
  if (oneout) {
    if (length(c) > 0 && !is.null(c[[1]])) {  # check argument here
      if (plus)
        out <- .cukxp(x, c, check = TRUE)
      else
        out <- .pinvxp(x, c, check = TRUE)
    } else if (length(c) == 0)  # check argument here
      out <- list()
    else {
      if (plus)
        out <- .cukx(x, c)
      else
        out <- x$X %*% c
    }
  } else {
    if (length(c) > 0 && !is.null(c[[1]])) {
      if (plus)
        out <- list(.cukxp(x, c, check = TRUE),  # X1
                    .cukx(x, .res(.setx(c))))  # X0
      else {
        out1 <- .pinvxp(x, c, check = TRUE)  # X1
        out1[out1 < x$tol] <- 0
        out2 <- x$X %*% .r(.setx(c))  # X0
        out2[out2 < x$tol] <- 0
        out <- list(out1, out2)
      }
    } else {
      if (length(c) == 0) {
        if (plus)
          out <- c(list(), .cukx(x))
        else
          out <- c(list(), x$X)
      } else {
        if (plus) {
          out1 <- .cukx(x, c)
          out1[out1 < x$tol] <- 0
          out2 <- .cukx(x)
          out2[out2 < x$tol] <- 0
          out <- list(out1, out2)
        } else {
          out <- list(x$X %*% c, x$X)
        }
      }
    }
  }
  return(out)
}

# space tested while keeping size of X$i0'
.i02X1 = function(x, c, plus = FALSE) {
  c <- .i02C(x, c)
  if (plus)
    out <- .c2tsp(x, c, plus = TRUE)
  else
    out <- .c2tsp(x, c)
  return(out)
}

# get the orthogonal compliment and project onto the estimable space'
.X02c <- function(x, y, plus = FALSE) {
  if (plus) {
    cukX0 <- y
    if (!is.list(cukX0))
      X0 <- c()
    else
      X0 <- .ox(x) %*% cukX0
  } else
    X0 <- y
  
  if (!.isspc(x))
    x <- .setx(x)
  if (x$rk == 0)
    stop("null rank x == 0")
  
  if (plus) {
    if ((is.null(cukX0) | length(cukX0) == 0) && x$rk != nrow(cukX0))
      stop("cukX0 of wrong size")
  } else {
    if ((is.null(X0) | length(X0) == 0) && nrow(x$X) != nrow(X0))
      stop("X0 of wrong size")
  }
  
  if (is.null(X0) | length(X0) == 0) {
    sc0 <- .setx(.pinvx(x, X0))
    if (sc0$rk)
      c <- .opp(x, .r(sc0))
    else
      c <- .opp(x)
    
    # dont know if this is equivalent to matlab command c  = c(:,any(c))
    c <- c[, colSums(abs(c)) > 0]
    sL <- ncol(x$X)
    
    if ((ncol(c) == 0) && (ncol(X0) != sL))
      c <- rep(0, sL)
  } else
    c <- .xpx(x)
  return(c)
}

# effective F degrees of freedom dof(idf, rdf)'
.i02edf <- function(x, i0, V) {
  if (.isspc(x))
    x <- .setx(x)
  i0 <- .iX0check(i0, ncol(x$X))
  if (misisng(V))
    V <- diag(nrow(x$X))
  
  r <- .trMV(x, V)
  m <- .trMV(.i02X1(x, i0), V)
  
  return(c(m$trMV^2 / m$trMVMV, r$trMV^2 / r$trMVMV))
}

# traces for effective df calculation. If oneout = TRUE then returns trRV.'
.trRV <- function(KWX, XV, oneout = FALSE) {
  if (!.isspc(KWX))
    KWX <- .setx(KWX)
  
  rk <- KWX$rk
  sL <- nrow(KWX$X)
  if (missing(XV)) {
    if (oneout) {
      if (is.null(rk))
        out <- sL
      else
        out <- sL - rk
    } else{
      if (is.null(rk))
        out <- list(trRV = sL, trRVRV = sL)
      else
        out <- list(trRV = sL - rk, trRVRV = sL)
    }

  } else {
    u <- KWX$u[, seq_len(rk)]
    if (oneout) {
      if (rk == 0)
        return(0)
      else
        trmv <- sum(u * (XV %*% u))
      return(sum(diag(XV)) - trmv)
    } else {
      if (rk == 0) {
        trmv <- 0
        tmp <- norm(XV, "f")^2
        trv <- sum(diag(XV))
      } else {
        Vu <- XV %*% u
        trv <- sum(diag(XV))
        tmp <- norm(XV, "F")^2
        tmp <- tmp - 2 * norm(Vu, "f")^2
        trRVRV <- tmp + norm(crossprod(u, Vu), "f")^2
        trmv <- sum(u * Vu)
      }
      trRV <- trv - trmv
    }
    rdf <- (trRV^2) / trRVRV
    out <- list(trRV = trRV, trRVRV = trRVRV, rdf = rdf)
  }
  return(out)
}

# compute the traceof MV, MVMV, and find the degrees of interest'
.trMV <- function(x, v, oneout = FALSE) {
  if (!.isspc(x))
    x <- .setx(x)
  rk <- x$rk
  if (length(rk) == 0) {
    warning("Rank is empty.")
    if (oneout)
      return(c())
    else
      return(list(trMV = c(), trMVMV = c(), idf = c()))
  } else {
    if (missing(v)) {
      if (oneout)
        return(rk)
      else
        return(list(trMV = rk, trMVMV = rk, idf = rk - 1))
    } else {
      u <- x$u[, seq_len(rk)]
      if (oneout)
        return(sum(t(u) * crossprod(u, v)))
      else {
        vu <- v %*% u
        trmv <- sum(u * vu)
        trmvmv <- norm(crossprod(u, vu), "f")^2
        idf <- (trmv^2) / trmvmv
        return(list(trMV = trmv, trMVMV = trmvmv, idf = idf))
      }
    }
  }
}

## from sp_FcUtil.m----
.fconfields <- function() {
  list(name = "",
       fieldType = "",
       c = matrix(0, 0, 0),
       X0 = matrix(0, 0, 0),
       X1 = matrix(0, 0, 0),
       iX0 = "",
       idf = 0)
}

# set contrast fields'
.setcon <- function(name, fieldType, action, c, x) {
  if (!is.character(name))
    stop("Name must be a character.")
  if (!(fieldType == "F" | fieldType == "T" | fieldType == "P"))
    stop("fieldType must be F, T, or P.")
  if (!(action == "c" | action == "c+" | action == "X0" | action == "ukX0" | action == "iX0"))
    stop("Action must be c, c+, X0, ukX0, iX0.")
  
  Fc <- .fconfields()
  Fc$name <- name
  Fc$fieldType <- fieldType
  if (Fc$fieldType == "T" && !any(action == c("c", "c+")))
    warning("enter T stat with contrast - here no check rank == 1")
  sC <- nrow(x$X)
  sL <- ncol(x$X)
  
  if (action == "c" | action == "c+") {
    Fc$iX0 <- action
    c[abs(c) < x$tol] <- 0
    if (length(c) == 0) {
      out <- .c2tsp(x, c(), plus = TRUE)
      suppressWarnings(Fc$X1$ukX1 <- out[[1]])
      suppressWarnings(Fc$X0$ukX0 <- out[[2]])
      Fc$c <- c
    } else if (nrow(c) != sL)
      stop("not contrast dim")
    else {
      if (action == "c+") {
        if (!.isinspp(x, c))
          c <- .opp(x, c)
      }
      if (Fc$fieldType == "T" && !.isT(x, c))
        stop("trying to define a T that looks like an F")
      Fc$c <- c
      out <- .c2tsp(x, c, plus = TRUE)
      suppressWarnings(Fc$X1$ukX1 <- out[[1]])
      suppressWarnings(Fc$X0$ukX0 <- out[[2]])
    }
  } else if (action == "X0") {
    Fc$iX0 <- action
    X0 <- c
    X0[X0 < x$tol] <- 0
    if (length(X0) == 0) {
      Fc$c <- .xpx(x)
      suppressWarnings(Fc$X1$ukX1 <- .cukx(x))
      suppressWarnings(Fc$X0$ukX0 <- matrix(0, 0, 0))
    } else if (nrow(X0) != sC) {
      stop("dimesion of X0 wrong in set")
    } else {
      Fc$c <- .X02c(x, X0)
      suppressWarnings(Fc$X0$ukX0 <- crossprod(.ox(x), X0))
      suppressWarnings(Fc$X1$ukX1 <- .c2tsp(x, Fc$c, plus = TRUE))
    }
  } else if (action == "ukX0") {
    Fc$iX0 <- action
    if (length(ukX0) == 0) {
      Fc$c <- .xpx(x)
      Fc$X1$ukX1 <- .cukx(x)
      suppressWarnings(Fc$X0$ukX0 <- matrix(0, 0, 0))
    } else if (nrow(X0) != sC)
      stop("dimension of X0 wrong in set")
    else {
      Fc$c <- .X02c(x, X0)
      suppressWarnings(Fc$X0$ukX0 <- crossprod(.ox(x), X0))
      suppressWarnings(Fc$X1$ukX1 <- .c2tsp(x, Fc$c, plus = TRUE))
    }
  } else {
    iX0 <- .iX0check(c, sL)
    Fc$iX0 <- iX0
    suppressWarnings(Fc$X0$ukX0 <- crossprod(.ox(x), matrix(x$X[, iX0])))
    if (length(iX0) == 0) {
      FC$c <- .xpx(x)
      suppressWarnings(Fc$X1$ukX1 <- .cukx(x))
    } else {
      Fc$c <- .i02c(x, iX0)
      suppressWarnings(Fc$X1$ukX1 <- .c2tsp(x, Fc$c, oneout = TRUE, plus = TRUE))
    }
  }
  return(Fc)
}

.X0 <- function(Fc, x) {
  if (!.isfcon(Fc))
    stop("argument is not a contrast structure")
  if (!.isspc(x))
    x <- .setx(x)
  
  if (is.list(Fc$X0))
    .ox(x) %*% Fc$X0$ukX0
  else
    Fc$X0
}

.X1 <- function(Fc, x) {
  if (!.isfcon(Fc))
    stop("Argument is not a contrast structure.")
  if (!.isspc(x))
    x <- .setx(x)
  
  if (is.list(Fc$X0))
    .ox(x) %*% Fc$X1$ukX1
  else
    Fc$X1
}

.isfcon <- function(Fc) {
  if (!is.list(Fc))
    return(0)
  b <- 1
  minnames <- names(.minFc())
  FCnames <- names(Fc)
  for (i in 1:length(minnames)) {
    b <- (b && any(minnames[i] == FCnames))
    if (!b)
      break
  }
  return(b)
}

.fconedf <- function(Fc, x, V, oneout = FALSE) {
  if (!.isfcon(Fc))
    stop("Fc must be Fcon")
  if (!.isspc(x))
    x <- .setx(x)
  
  if (!.isemptyX1(Fc)) {
    trmv <- .trMV(.X1(Fc, x), V)
  } else {
    trmv <- c(0, 0)
  }
  
  if (!trmv[2]) {
    edf_tsp <- 0
    warning("edf_tsp = 0")
  }
  
  if (oneout) {
    return(edf_tsp)
  } else {
    out <- .trRV(x, V)
    if (!out$trMVMV) {
      edf_Xsp <- 0
      warning("edf_Xsp = 0")
    } else
      edf_Xsp <- out$trRV^2 / out$trRVRV
    return(edf_tsp, edf_Xsp)
  }
}

# extra sum of squares matrix for betas from contrast'
.hsqr <- function(Fc, x) {
  if (!.isfcon(Fc))
    stop("Fc must be F-contrast")
  if (!.isset(Fc))
    stop("F-contrast must be set")
  if (!.isspc(x))
    x <- .setx(x)
  
  if (.isemptyX1(Fc)) {
    if (!.isemptyX0(Fc))
      return(rep(0, ncol(x$X)))
    else
      stop("Fc must be set")
  } else {
    if (is.list(Fc$X0))
      return(crossprod(.ox(.setx(Fc$X1$ukX1)), .cukx(x)))
    else
      return(crossprod(.ox(.setx(Fc$X1)), x$X))
  }
}

# extra sum of squares matrix for betas from contrast'
.H <- function(Fc, sX) {
  if (!.isfcon(Fc))
    stop("Fc must be F-contrast")
  if (.isset(Fc))
    stop("Fcon must be set")
  if (!.isspc(x))
    x <- .setx(x)
  if (.isemptyX1(Fc)) {
    if (!.isemptyX0(Fc))
      return(rep(0, ncol(x$X)))
    else
      stop("Fc must be set")
  } else {
    if (is.list(Fc$X0))
      return(crossprod(.hsqr(Fc, sX)))
    else
      return(Fc$c %*% tcrossprod(MASS::ginv(crossprod(Fc$X1)), Fc$c))
  }
}

# fitted data corrected for confounds defined by Fc'
.yc <- function(Fc, x, b) {
  if (!.isfcon(Fc))
    stop("Fc must be F-contrast")
  if (!.isset(Fc))
    stop("Fcon must be set")
  if (!.isspc(x))
    x <- .setx(x)
  if (ncol(x$X) != nrow(b))
    stop("ncol(x$X) must equal nrow(b)")
  
  if (.isemptyX1(Fc)) {
    if (!.isemptyX0(Fc))
      return(matrix(0, nrow(x$X), ncol(b)))
    else
      stop("Fc must be set")
  } else
    return(x$X %*% .pinvxpx(x) %*% .H(Fc, x) %*% b)
}

# fitted data corrected for confounds defined by Fc'
.y0 <- function(Fc, x, b) {
  if (!.isfcon(Fc))
    stop("Fc must be F-contrast")
  if (!.isset(Fc))
    stop("F-contrast must be set")
  if (!.isspc(x))
    x <- .setx(x)
  if (ncol(x$X) != nrow(b))
    stop("ncol(x$X) must equal nrow(b)")
  
  if (.isemptyX1(Fc)) {
    if (!.isemptyX0(Fc))
      return(x$X %*% b)
    else
      stop("Fc must be set")
  } else
    return(x$X %*% (diag(ncol(x$X)) - .xpx(sX)) %*% .H(Fc, x) %*% b)
}

# Fc orthogonalisation'
.fcortho <- function(Fc1, x, Fc2) {
  L1 <- length(Fc1)
  
  if (~L1) {
    warning("no contrast provided")
    return(c())
  }
  for (i in seq_len(L1)) {
    if (!.isfcon(Fc1[[i]]))
      stop(paste("Fc1[[", i, "]] must be a contrast", sep = ""))
  }
  L2 <- length(Fc2)
  if (!L2)
    stop("must have at least one contrast in Fc2")
  for (i in seq_len(L2)) {
    if (!.isfcon(Fc2[[i]]))
      stop(paste("Fc2[[", i, "]] must be a contrast", sep = ""))
  }
  if (!.isspc(x))
    x <- .setx(x)
  
  strnm <- paste(Fc2[[1]]$name, sep = "")
  tmpc <- Fc2[[1]]$c
  for (i in 2:L2) {
    strnm <- paste(strnm, " ", Fc2[[i]]$name, sep = "")
    tmpc <- cbind(tmpc, Fc2[[i]]$c)
  }
  Fc2 <- .setcon(strnm, "F", "c+", tmpc, x)
  
  if (.isemptyX1(Fc2) || .isnull(Fc2, x))
    return(Fc1)
  else {
    out <- list()
    for (i in seq_len(L1)) {
      if (.isemptyX1(Fc1[[i]]) || .isnull(Fc1[[i]], x))
        out[[i]] <- Fc1[[i]]
      else {
        c1_2  = Fc1[[i]]$c - .H(Fc2, x) %*% .xpx(x, Fc1[[i]]$c)
        out  = .setcon(paste("(", Fc1[[i]]$name, "|_ (", Fc2$name, "))", sep = ""), Fc1$fieldType, "c+", c1_2, x)
      }
    }
    return(out)
  }
}

# are contrasts orthogonals?
.Rortho <- function (Fc1, sX, Fc2) {
  if (length(Fc1) == 0 | length(Fc1[[1]]) == 0)
    stop("must provide at least one non empty contrast")
  if (!.isspc(sX))
    sX <- .setxt(sX)
  
  L1 <- length(Fc1)
  for (i in seq_len(L1)) {
    if (!.isfcon(Fc1[[i]]))
      stop(paste("Fc1[[", i, "]] must be a contrast", sep = ""))
  }
  L2 <- length(Fc2)
  for (i in seq_len(L2)) {
    if (!.isfcon(Fc2[[i]]))
      stop(paste("Fc2[[", i, "]] must be a contrast", sep = ""))
  }
  
  if (length(Fc2) == 0) {
    if (length(Fc1) <= 1)
      b <- 0
    else {
      c1 <- Fc1[[1]]$c
      for (i in 2:length(Fc1))
        c1 <- cbind(c1, Fc1[[i]]$c)
      b <- !any(abs(upper.tri(crossprod(c1, .xpx(sX, c1)))) > sX$tol)
    }
  } else {
    c1 <- Fc1[[1]]$c
    for (i in 2:length(Fc1))
      c1 <- cbind(c1, Fc1[[i]]$c)
    c2 <- Fc2[[1]]$c
    for (i in 2:length(Fc2))
      c2 <- cbind(c1, Fc1[[i]]$c)
    b <- !any(abs(crossprod(c1, .pinvxpx(x, c2)) > sX$tol))
  }
  return(b)
}

# Fc1 is in list of contrasts Fc2'
.C1inC2 <- function(Fc1, x, Fc2) {
  l1 <- length(Fc1)
  l2 <- length(Fc2)
  for (j in seq_len(l1)) {
    if (!.isinspp(x, Fc1[[j]]$c))
      c1 <- .opp(x, Fc1[[j]]$c)
    else
      c1 <- Fc1[[j]]$c
    sc1 <- .setx(c1)
    S <- Fc1[[j]]$fieldType
    
    boul = 0
    idxFc1 <- c()
    idxFc2 <- c()
    for (i in seq_len(l2)) {
      if (Fc2[[i]]$fieldType == S) {
        boul <- .equal(sc1, .opp(x, Fc2[[i]]$c))
        
        if (boul && S == "T") {
          atmp <- .X1(Fc1[[j]], x)
          btmp <- .X1(Fc2[[i]], x)
          boul <- !any(crossprod(atmp, btmp) < 0)
        }
        if (boul) {
          idxFc1 <- c(idxFc1, j)
          idxFc2 <- c(idxFc2, i)
        }
      }
    }
  }
  return(list(idxFc1, idxFc2))
}

# Fc list unique
.listnotunique <- function(Fc, x) {
  L1 <- length(Fc)
  if (!L1) {
    warning("no contrast provided")
    return(c())
  } else {
    for (i in seq_len(L1)) {
      if (!.isfcon(Fc[[i]]))
        stop(paste("Fc[[", i, "]] must be a contrast"))
    }
    if (!.isspc(x))
      x <- .setx(x)
    return(.notunique(Fc, x))
  }
}

.isnull <- function(Fc, x) {
  any(colSums(abs(.opp(x, Fc$c))) > 0)
}

.isT <- function(x, c) {
  boul <- 1
  if (!.isinspp(x, c))
    c <- .opp(x, c, check = TRUE)
  # SPM12 uses "rank(c) > 1" instead of "sum(abs(c))) > 1"
  if ((sum(abs(c)) > 1) || any(crossprod(c) < 0))  
    boul <- 0
  return(boul)
}

.minFc <- function() {
  out <- list()
  out$name <- ""
  out$fieldType <- ""
  out$c <- matrix(0, 0, 0)
  out$X0 <- matrix(0, 0, 0)
  out$X1 <- matrix(0, 0, 0)
  return(out)
}

.isset <- function(Fc) {
  !.isemptyX0(Fc) | !.isemptyX1(Fc)
}

.isemptyX1 <- function(Fc) {
  if (is.list(Fc$X0)) {
    b <- is.null(Fc$X1$ukX1)
    if (!(b == is.null(Fc$c) | (b == (length(Fc$c) == 0))))
      stop("Contrast is internally inconsistent")
  } else {
    b <- is.null(Fc$X1) | (length(Fc$X1) == 0)
    if (!(b == is.null(Fc$c) | (b == (length(Fc$c) == 0))))
      stop("Contrast is internally inconsistent")
  }
  return(b)
}

.isemptyX0 <- function(Fc) {
  if (is.list(Fc$X0))    
    is.null(Fc$X0$ukX0)
  else
    is.null(Fc$X0)
}

.in <- function(Fc1, x, Fc2) {
  L1 <- length(Fc1)
  L2 <- length(Fc2)
  
  idxFc1 <- c()
  idxFc2 <- c()
  for (j in seq_len(L1)) {
    if (!.isinspp(x, Fc1[[j]]$c))
      c1 <- .opp(x, Fc1[[j]]$c)
    else
      c1 <- Fc1[[j]]$c
    sc1 <- .setx(c1)
    S <- Fc1[[j]]$fieldType
    
    boul <- 0
    for (i in seq_len(L2)) {
      if (Fc2[[i]]$fieldType == S) {
        boul <- .equal(sc1, .opp(x, Fc2[[i]]$c))
        if (boul && S == "T") {
          atmp <- .X1(Fc1[[j]], x)
          btmp <- .X1(Fc2[[i]], x)
          boul <- !any(crossprod(atmp, btmp) < 0)
        }
        if (boul) {
          idxFc1 <- c(idxFc1, j)
          idxFc2 <- c(idxFc2, i)
        }
      }
    }
    
  }
  return(boul)
}

# Fc list unique
.notunique <- function(Fc, x) {
  l <- length(Fc)
  if (l == 1)
    out <- c()
  else {
    out <- list(1 + .in(Fc[1], x, Fc[2:l]), 1 + .notnunique(Fc[2:l], x))
    for (i in length(out))
      out[[i]] <- out[[i]] + 1
  }
  return(out)
}

.cX1 <- function(x, c, fieldType) {
  if (missing(x))
    x <- KWX
  
  out <- list()
  out$name <- name
  out$fieldType <- fieldType
  out$c <- matrix(c, ncol = 1)
  
  tmp <- .c2tsp(x, c, plus = TRUE)
  out$X1$ukX1 <- tmp[[1]]
  out$X0$ukX0 <- tmp[[2]]
  
  X1 <- .ox() %*% out$X1$ukX1
}

# general utilities----
# Discrete cosine transform
# 
# Creates a matrix for the first few basis functions of a one dimensional discrete cosine transform.
# 
# @param N dimension
# @param K order
# @param n optional points to sample
# @param f type of transform
.dctmtx <- function(N, K, n, f) {
  d <- 0
  
  if (nargs() == 1)
    K <- N
  if (args() == 1 | args() == 2)
    n <- seq_len(N - 1) - 1
  else if (nargs() > 2) {
    if (f == "diff")
      d <- 1
    else if (f == "diff2")
      d <- 2
    if (missing(n))
      n <- seq_len(N - 1) - 1
  }
  
  C <- matrix(0, lenght(n), K)
  
  if (d == 0) {
    C[, 1] <- rep(1, dim(n)[1]) / sqrt(N)
    for (k in 2:K)
      C[, k] <- sqrt(2 / N) * cos(pi * (2 * n + 1) * (k - 1) / (2 * N))
  } else if (d == 1) {
    for (k in 2:K)
      C[, k] = -2^(1 / 2) * (1 / N)^(1 / 2) * sin(1 / 2 * pi * (2 * n * k - 2 * n + k - 1) / N) * pi * (k - 1) / N;
  } else if (d == 2) {
    for (k in 2:K)
      C[, k] = -2^(1 / 2) * (1 / N)^(1 / 2) * cos(1 / 2 * pi * (2 * n + 1) * (k - 1) / N) * pi^2 * (k - 1)^2 / N^2;
  }
  return(C)
}

# @param dfdx = df/dx
# @param f = dx/dt
# @param t = integration time: (default t = Inf);
# if t is a cell (i.e., {t}) then t is set to:
# exp(t - log(diag(-dfdx))
#
.dx <- function(dfdx, f, t = Inf) {
  # t is a regulariser
  if (length(t) == 1)
    t <- exp(t - log(diag(-dfdx)))
  
  if (min(t) > exp(16)) {
    dx = - MASS::ginv(dfdx) %*% as.matrix(f, ncol = 1)
    dx =  as.array(dx, dim = dim(f));
  } else {
    # ensure t is a scalar or matrix
    if (length(t) > 1)
      t = diag(t)
    
    q <- matrix(0, nrow = max(dim(dfdx)) + 1, ncol = 1)
    q[1, 1] <- 1
    
    # augment Jacobian and take matrix exponential
    Jx <- rbind(0, cbind(t %*% f, t %*% dfdx))
    dx <- Matrix::expm(Jx) %*% q
    dx <- dx[2:nrow(dx), ]
  }
  dx
}

# Filter 
#
# @param Y data matrix
# @param K filter matrix
# @param K struct array containing partition specific specifications
# @param K[[s]]$RT observation interval in seconds
# @param K[[s]]$row row of Y constituting block/partitions
# @param K[[s]]$HParams cut-off period in seconds
# @param K[[s]]$X0 low frequencies to be removed (DCT)
# @return filtered data
# 
# out <- iModelMake(X = z$X, y = z$y[i], data = z$iData)
#
# @export iFilter
.filter <- function(K, Y) {
  if (missing(Y) && is.data.frame(K)) {
    # set K$X0
    out <- list()
    for (s in seq_len(length(K))) {
      # create filter index from filter data.frame
      subout <- list()
      subout$row <- seq_len(nrow(K))[K$Filter == paste("F", s, sep = "")]
      subout$RT <- K$RT[subout$row[1]]
      subout$HParam <- K$HParam[subout$row[1]]
      
      # determine low frequencies to be removed
      nk <- length(subout$row)
      n <- floor(2 * (nk * subout$RT) / subout$HParam + 1)
      X0 <- .dctmtx(nk, n)
      subout$X0 <- X0[, 2:ncol(X0)]
      out[[i]] <- subout
    }
    return(out)
  } else {
    if (is.data.frame(K)) {
      # ensure requisite fields are present
      if (is.null(K[[1]]$X0))
        K <- iFilter(K)
      
      # apply high pass filter
      if (length(K) == 1 && length(K[[1]]$row == ncol(Y)))
        Y <- Y - K[[1]]$X0 %*% crossprod(K[[1]]$X0, Y)
      else {
        for (i in seq_len(length(K))) {
          
          # select data
          y <- Y[K[[i]]$row, ]
          
          # apply high pass filter
          y <- y - K[[i]]$X0 %*% crossprod(K[[i]]$X0, y)
          
          # reset filtered data in Y
          Y[K[[i]]$row, ] <- y
        }
      }
    } else {
      Y <- K * Y
    }
    return(Y)
  }
}