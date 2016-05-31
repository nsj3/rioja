randomPTF <- function (spec, env, fun, ncol=1, nVar, nTF=500, verbose=TRUE, do.parallel=FALSE, ...) {
  if (do.parallel) {
     if (!requireNamespace("foreach", quietly = TRUE)) 
       stop("Package foreach is required for parallel operation. Please install and try again.")
  }
  do.TF <- function(y, x, nsp, nsam, fun, ncol, nVar, ...) {
    res2 <- vector("numeric", length=nsp)
    sel.sp <- sample.int(nsp, nVar)
    # bootstrap and test set sample indices      
    boot <- sample.int(nsam, replace=TRUE)
    test <- setdiff(1:nsam, boot)
    x.b <- x[boot]
    y.b <- y[boot, sel.sp]
    x.t <- x[test]
    y.t <- y[test, sel.sp]
    # remove taxa with no occur in bootstrap sample
    sel <- apply(y.b, 2, sum) > 0.001
    y.b <- y.b[, sel]
    y.t <- y.t[, sel]
    sel.sp <- sel.sp[sel]
    # remove samples with no taxa
    sel <- apply(y.b, 1, sum) > 0.01
    y.b <- y.b[sel, ]
    x.b <- x.b[sel]
    call.fit <- paste0(fun, ".fit")
    mod <- do.call(call.fit, args=list(y.b, x.b, lean=TRUE, ...))
    call.pred <- paste0("predict.internal.", fun)
    pred.oob <- do.call(call.pred, args=list(mod, y.t, lean=TRUE))
#    mod <- WA.fit(y.b, x.b, lean=TRUE)
#    pred.oob <- rioja:::predict.internal.WA(mod, y.t, lean=TRUE)
    mseOOB <- mean((pred.oob[, ncol] - x.t)^2)
    nv <- ncol(y.t)
    nsam.t <- nrow(y.t)
    for (j in 1:nv) {
      y.t2 <- y.t
      y.t2[, j] <- y.t[sample.int(nsam.t), j]
      pred.oob_h <- do.call(call.pred, args=list(mod, y.t2, lean=TRUE))
#      pred.oob_h <- rioja:::predict.internal.WA(mod, y.t2, lean=TRUE)
      mseOOB_H <- mean((pred.oob_h[, ncol] - x.t)^2)
      res2[sel.sp[j]] <- mseOOB_H - mseOOB
    }
    res2
  }
  fun <- substitute(fun)
  y <- as.matrix(spec)
  x <- as.matrix(env)
  nsam <- nrow(y)
  nsp <- ncol(spec)
  if (missing(nVar))   
    nVar <- max(as.integer(nsp/3), 1)
  if (do.parallel) {

    # this is some magic copied from plyr
    .paropts <- NULL
    .paropts$.combine <- function(...) NULL
    i <- seq_len(nTF)
    fe_call <- as.call(c(list(quote(foreach::foreach), i = i), .paropts))
    fe <- eval(fe_call)
    res <- foreach::`%dopar%`(fe, do.TF(y, x, nsp, nsam, fun, ncol, nVar, ...))
    
#     res <- foreach::foreach(1:nTF, .packages=c('rioja')) %dopar% { do.TF(y, x, nsp, nsam, fun, ncol, nvar, ...) }
     res <- t(sapply(res, "["))
  } else {
    if (verbose) {
      writeLines("Running:")
      pb <- txtProgressBar(min = 0, max = 1, style = 3)
      on.exit(close(pb))
    }
    res <- matrix(ncol=nsp, nrow=nTF)
     for (i in 1:nTF) {
       if (verbose) {
          setTxtProgressBar(pb, i/nTF)
       }        
       res[i, ] <- do.TF(y, x, nsp, nsam, fun, ncol, nVar, ...)
     }
    #      res <- foreach(1:nTF, .packages=c('rioja')) %do% do.TF(y, x, nsp, nsam, fun, ncol, nvar, ...)
  }
  colnames(res) <- colnames(y)
  SumErr <- colSums(res, na.rm=TRUE)
  nTree <- colSums(!is.na(res)) 
  VI <- SumErr / nTree
  names(VI) <- colnames(y)
  VI <- data.frame(VI=sort(VI, decreasing=TRUE))
  res <- list(VI=VI, spec=y, env=x)
  class(res) <- "randomPTF"
  res
}

plot.randomPTF <- function(x, use.pointLabel=TRUE, ...) {
  mt <- match(rownames(x$VI), colnames(x$spec))
  mx <- apply(x$spec, 2, max)
  plot(x$VI[, 1], mx[mt], xlab="VI", ylab="Maximum abundance", ...)
  if (use.pointLabel) {
    if (!requireNamespace("maptools")) {
      text(x$VI[, 1], mx[mt], labels=rownames(x$VI), ...)
    } else {
      maptools::pointLabel(x$VI[, 1], mx[mt], labels=rownames(x$VI), ...)
    }
  }
}

print.randomPTF <- function(x, ...) {
   print(x$VI, ...)
}