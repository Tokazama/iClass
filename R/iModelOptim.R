#' @export 
#' @docType methods
#' @details \strong{iModelOptim} Optimize iModel objects.
#' @rdname iModel-method
iModelOptim <- function(object, its) {
  # from spm_est_non_sphericity
  if (!is.null(object@xVi$Fcontrast))
    con <- .setcon("User-specified contrast", "F", "c", object@xVi$Fcontrast, object@X$KWX)
  else {
    iX0 <- rep(TRUE, ncol(object@X$X))
    iX0[colnames(object@X$X) == "(Intercept)"] <- FALSE
    con <- .setcon("Effects of interest", "F", "iX0", iX0, object@X$KWX)
  }
  
  if ((!is.null(con$c)) | length(con$c) == 0) {
    X1 <- .X1(con, object@X$KWX)
    hsqr <- .hsqr(con, object@X$KWX)
    trmv <- .trMV(X1, oneout = TRUE)
  } else {
    trmv <- 1
    hsqr <- Inf
  }
  trrv <- .trRV(object@X$KWX, oneout = TRUE)
  
  # threshold for voxels entering non-sphericity estimates
  object@xVi$UFp <- object@control$cf
  UF <- qf(1 - object@control$cf, trmv, trrv)
  
  
  chunksize <- object@Y@chunksize[2]
  chunkseq <- seq_len(chunksize)
  nchunks <- floor(object@dims$nvox / chunksize)
  Cy <- matrix(0, object@dims$nimg, object@dims$nimg)
  s <- 0
  for (j in seq_len(nchunks)) {
    KWY <- .filter(object@X$K, object@X$W %*% object@Y[, chunkseq])
    object@beta[, chunkseq] <- object@X$pKX %*% KWY
    object@res[, chunkseq] <- .res(object@X$KWX, KWY)
    object@mrss[, chunkseq] <- colSums(object@res[, chunkseq]^2) / object@X$trRV
    if (object@control$scr)
      object@res[, chunkseq] <- t(t(object@res[, chunkseq]) * (1 / as.numeric(object@mrss[, chunkseq])))
    
    good <- chunkseq[(colSums((hsqr %*% object@beta[, chunkseq])^2) / trmv) > (UF * object@mrss[, chunkseq])]
    
    if (length(good) > 0) {
      q <- as.numeric(sqrt(1 / object@mrss[, good]))
      q <- t(t(object@Y[, good]) * q)
      Cy <- tcrossprod(q)
      s <- s + length(good)
    }
    chunkseq <- chunkseq + chunksize
  }
  if (chunkseq[1] < object@dims$nvox) {
    chunkseq <- chunkseq[1]:object@dims$nvox
    KWY <- .filter(object@X$K, object@X$W %*% object@Y[, chunkseq])
    object@beta[, chunkseq] <- object@X$pKX %*% KWY
    object@res[, chunkseq] <- .res(object@X$KWX, KWY)
    object@mrss[, chunkseq] <- colSums(object@res[, chunkseq]^2) / object@X$trRV
    if (object@control$scr)
      object@res[, chunkseq] <- t(t(object@res[, chunkseq]) * (1 / as.numeric(object@mrss[, chunkseq])))
    
    good <- chunkseq[(colSums((hsqr %*% object@beta[, chunkseq])^2) / trmv) > (UF * object@mrss[, chunkseq])]
    
    if (length(good) > 0) {
      q <- as.numeric(sqrt(1 / object@mrss[, good]))
      q <- t(t(object@Y[, good]) * q)
      Cy <- tcrossprod(q)
      s <- s + length(good)
    }
  }
  
  if (s > 0)
    Cy <- Cy / s
  else
    warning("No voxels appear to be significant.")
  
  if (is.list(object@X$K)) {
    m <- length(object@xVi$Vi)
    h <- rep(0, m)
    V <- matrix(0, n, n)
    for (i in seq_len(length(object@X$K))) {
      # extract blocks from bases
      q <- object@X$K[[i]]$row
      p <- c()
      QP <- list()
      for (j in seq_len(m)) {
        if (any(object@xVi$Vi[[j]][q, q] != 0)) {
          Qp <- lappend(Qp, object@xVi$Vi[[j]][q, q])
          p <- c(p, j)
        }
      }
      
      # design space for ReML (with confounds in filter)
      Xp <- X[q, ]
      Xp <- lappend(Xp, object@X$K[[i]]$X0)
      
      # ReML
      reml <- .reml(Cy[q, q], Xp, Qp, its = its)
      V[q, q] <- V[q, q] + reml$Vp
      h[p] <- reml$hp
    }
  } else {
    reml <- .reml(Cy, object@X$X, object@xVi$Vi, its = its)
    V <- reml$V
    h <- reml$h
  }
  
  V <- V * n / sum(diag(V))
  
  object@xVi$h <- h
  object@xVi$V <- V
  object@xVi$Cy <- Cy
  return(x)
}


.reml <- function(YY, X, Q, N, D, t, hE, hP, its) {
  if (missing(N))
    N <- 1
  if (missing(D))
    D <- 0
  if (missing(t))
    t <- 0
  if (missing(hE))
    hE <- 0
  if (missing(hP))
    hP <- 1e-16
  
  # ortho-normalise X----
  if (missing(X))
    stop("Please specify X")
  else
    X <- svd(X)$u
  if (!is.list(Q))
    Q <- list(Q)
  
  # dimensions----
  n <- nrow(Q[[1]])
  m <- length(Q)
  
  # catch NaNs----
  W <- Q
  #q <- is.finite(YY)
  #YY <- YY[q, q]
  #for (i in 1:m)
  #  Q[[i]] <- Q[[i]][q, q]
  
  # initialise h and specify hyperpriors----
  h <- matrix(0, m, 1)
  for (i in 1:m)
    h[i, 1] <- sum(diag(Q[[i]]))
  hE <- matrix(0, m, 1) + hE
  hP <- matrix(0, m, m) * hP
  dF <- Inf
  dh <- matrix(0, m, 1)
  dFdh <- matrix(0, m, 1)
  dFdhh <- matrix(0, m, m)
  PQ <- list()
  
  # ReML (EM/VB)----
  for (it in seq_len(its)) {
    # compute current estimate of covariance----
    #C <- matrix(0, n, n)
    C <- 0
    for (i in 1:m)
      C <- C + Q[[i]] * h[i]
    
    # positive [semi]-definite check----
    # might be able to use nlme::pdMat()
    for (i in 1:D) {
      if (min(eigen(C)$values) < 0) {
        t <- t - 1
        h <- h - dh
        dh <- .dx(dFdhh, dFdh, t)
        h <- h + dh
        C <- matrix(0, n, n)
        for (j in 1:m)
          C <- C + Q[[i]] * h[i]
      }
    }
    
    # E-step: conditional covariance cov(B|y) {Cq}----
    iC <- solve(C)
    iCX <- iC %*% X
    Cq <- solve(crossprod(X, iCX))
    # if (!any(X != 0))
    #   Cq <- solve(crossprod(X, iCX))
    # else
    #   Cq <- 0
    
    # M-step: ReML estimate of hyperparameters----
    P <- iC - iCX %*% Cq %*% t(iCX) # P = iV - iV %*% X %*% solve(crossprod(X, iV) %*% X) %*% crossprod(X, iV)
    U <- diag(n) - P %*% YY / N
    
    # dF/dh
    for (i in 1:m) {
      PQ[[i]] <- P %*% Q[[i]]
      dFdh[i] <- -sum(diag(PQ[[i]] %*% U)) * N / 2
    }
    
    # expected curvature E{dF / dhhh}
    for (i in 1:m) {
      for (j in 1:m) {
        # dF/dhh
        dFdhh[i, j] <- -sum(diag(PQ[[i]] %*% PQ[[j]])) * N / 2
        dFdhh[j, i] <- dFdhh[i, j]
      }
    }
    
    # add hyperpriors
    e <- h - hE
    dFdh <- dFdh - hP %*% e
    dFdhh <- dFdhh - hP
    
    # fisher scoring: update dh = -inv(ddF/dhh) * dF / dh
    dh <- .dx(dFdhh, dFdh, {t})
    h <- h + dh
    
    # predicted change in F - increase regularisation if increasing
    pF <- crossprod(dFdh, dh)
    if (pF > dF)
      t <- t - 1
    else
      t <- t + 1/4
    
    # final estimate of covarience (with missing data points)
    if (dF < 1e-1)
      break
  }
  
  # rebuild predicted covariance----
  V <- 0
  for (i in 1:m) {
    V <- V + W[[i]] * h[i]
  }
  
  # check V is positive semi-definite
  if (!D) {
    if (min(eigen(V)$values) < 0)
      .reml(YY, X, Q, N, 1, 2, hE[1], hP[1])
  }
  
  # log evidence = ln p(y|X, Q) = ReML objective = F = trace(R' *iC * R * YY) / 2 ...
  Ph <- - dFdhh
  
  
  if (nargs() > 4) {
    # tr(hP * inv(Ph)) - nh + tr...
    Ft <- sum(diag(hP %*% MASS::ginv(Ph))) - length(Ph) - length(Cq)
    
    # complexity
    Fc <- Ft / 2 +
      crossprod(e, hP) %*% e/2 +
      determinant(Ph %*% MASS::ginv(hP), logarithm = TRUE)$modulus / 2
    
    # accuracy - lnp(Y|h)
    Fa = Ft / 2 -
      sum(diag(C * P * YY * P)) / 2 -
      N * n * log(2 * pi) / 2 -
      N * determinant(C, logarithm = TRUE)$modulus / 2
    
    # free-energy
    FE <- Fa - Fc
    return(list(V = V, h = h, Ph = Ph, FE = FE, Fa = Fa, Fc = Fc))
  } else {
    return(list(V = V, h = h, Ph = Ph))
  }
}