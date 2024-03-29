sample.perm <- function(x, k = NULL, R = 1000, replace = FALSE){
  if(is.null(k)) k <- length(x)
  max.R <- try(npermutations(x, k = k, replace = replace), silent = TRUE)
  if(inherits(max.R, "try-error")) max.R <- Inf
  if(max.R < R){
    warning("The requested number of permutations is larger than ",
            "the total number of possible permutations.\n",
            "Hence all possible permutations are computed.")
    res <- permutations(x, k = k, replace = replace)
  }else{
    res <- NULL
    iter <- R
    repeat{
      res <- rbind(res, 
                   t(replicate(iter, sample(x, size = k, 
                                            replace = replace))))
      if(all(!duplicated(res))) break
      iter <- R - sum(!duplicated(res))
      res <- res[!duplicated(res),]
    }
  }
  res
}
perm.t.test <- function (x, ...){ 
  UseMethod("perm.t.test")
}
perm.t.test.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
                        mu = 0, paired = FALSE, var.equal = FALSE, 
                        conf.level = 0.95, R = 9999, symmetric = TRUE, ...){
  alternative <- match.arg(alternative)
  if(!missing(mu) && (length(mu) != 1 || is.na(mu))) 
    stop("'mu' must be a single number")
  if(!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
                               conf.level < 0 || conf.level > 1)) 
    stop("'conf.level' must be a single number between 0 and 1")
  if(!is.null(y)){
    dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
    if (paired) 
      xok <- yok <- complete.cases(x, y)
    else{
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  }else{
    dname <- deparse1(substitute(x))
    if (paired) 
      stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]
  if(paired){
    x <- x - y
    y <- NULL
  }
  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)
  if (is.null(y)) {
    if (nx < 2) 
      stop("not enough 'x' observations")
    df <- nx - 1
    stderr <- sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx)) 
      stop("data are essentially constant")
    tstat <- (mx - mu)/stderr
    method <- if (paired) "Permutation Paired t-test" else "Permutation One Sample t-test"
    estimate <- setNames(mx, if (paired) "mean of the differences" else "mean of x")
    x.cent <- x - mx
    X <- abs(x.cent)*sample.perm(c(-1,1), k = nx, R = R, replace = TRUE)
    MX <- rowMeans(X)
    VX <- rowSums((X-MX)^2)/(nx-1)
    STDERR <- sqrt(VX/nx)
    perm.stderr <- mean(STDERR)
    TSTAT <- MX/STDERR
    EFF <- MX+mx
    perm.estimate <- mean(EFF)
    if(paired){
      names(perm.estimate) <- "permutation mean of the differences" 
    }else{
      names(perm.estimate) <- "permutation mean of x"
    } 
  }else{
    ny <- length(y)
    if(nx < 1 || (!var.equal && nx < 2)) 
      stop("not enough 'x' observations")
    if(ny < 1 || (!var.equal && ny < 2)) 
      stop("not enough 'y' observations")
    if(var.equal && nx + ny < 3) 
      stop("not enough observations")
    my <- mean(y)
    vy <- var(y)
    method <- paste("Permutation", paste(if (!var.equal) "Welch", "Two Sample t-test"))
    estimate <- c(mx, my)
    names(estimate) <- c("mean of x", "mean of y")
    z <- c(x, y)
    Z <- sample.perm(z, k = nx+ny, R = R)
    X <- Z[,1:nx]
    Y <- Z[,(nx+1):(nx+ny)]
    MX <- rowMeans(X)
    MY <- rowMeans(Y)
    EFF <- (MX+mx) - (MY+my)
    if(var.equal){
      df <- nx + ny - 2
      v <- 0
      if (nx > 1) 
        v <- v + (nx - 1) * vx
      if (ny > 1) 
        v <- v + (ny - 1) * vy
      v <- v/df
      stderr <- sqrt(v * (1/nx + 1/ny))
      V <- (rowSums((X-MX)^2) + rowSums((Y-MY)^2))/df
      STDERR <- sqrt(V*(1/nx + 1/ny))
    }else{
      stderrx <- sqrt(vx/nx)
      stderry <- sqrt(vy/ny)
      stderr <- sqrt(stderrx^2 + stderry^2)
      df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1))
      VX <- rowSums((X-MX)^2)/(nx-1)
      VY <- rowSums((Y-MY)^2)/(ny-1)
      STDERR <- sqrt(VX/nx + VY/ny)
    }
    perm.stderr <- mean(STDERR)
    perm.estimate <- mean(EFF) 
    names(perm.estimate) <- "permutation difference of means"
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))) 
      stop("data are essentially constant")
    tstat <- (mx - my - mu)/stderr
    TSTAT <- (MX - MY)/STDERR
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    perm.pval <- mean(TSTAT < tstat)
    cint <- c(-Inf, tstat + qt(conf.level, df))
    perm.cint <- c(-Inf, quantile(EFF, conf.level))
  }else if(alternative == "greater") {
    perm.pval <- mean(TSTAT > tstat)
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
    perm.cint <- c(quantile(EFF, 1-conf.level), Inf)
  }else{
    pval <- 2 * pt(-abs(tstat), df)
    if(symmetric)
      perm.pval <- mean(abs(TSTAT) > abs(tstat))
    else
      perm.pval <- 2*min(mean(TSTAT <= tstat), mean(TSTAT > tstat))
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
    perm.cint <- quantile(EFF, c(alpha/2, 1-alpha/2))
  }
  cint <- mu + cint * stderr
  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- if (paired || !is.null(y)) "difference in means" else "mean"
  attr(cint, "conf.level") <- conf.level
  attr(perm.cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = df, p.value = pval, 
               perm.p.value = perm.pval,
               conf.int = cint, perm.conf.int = perm.cint,
               estimate = estimate, perm.estimate = perm.estimate, 
               null.value = mu, stderr = stderr, perm.stderr = perm.stderr,
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- c("perm.htest", "htest")
  rval
}
perm.t.test.formula <- function (formula, data, subset, na.action, ...){
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]), 
                                                                  "term.labels")) != 1L)) 
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L) 
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("perm.t.test", c(DATA, list(...)))
  y$data.name <- DNAME
  if (length(y$estimate) == 2L) 
    names(y$estimate) <- paste("mean in group", levels(g))
  y
}
print.perm.htest <- function (x, digits = getOption("digits"), prefix = "\t", ...) {
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n", sep = "")
  out <- character()
  if (!is.null(x$perm.p.value)) {
    bfp <- format.pval(x$perm.p.value, digits = max(1L, digits - 3L))
    cat("(Monte-Carlo) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("=", bfp), "\n")
  }
  if (!is.null(x$perm.estimate)) {
    cat(paste(names(x$perm.estimate), "(SE) =", 
              format(x$perm.estimate, digits = digits),
              paste("(", format(x$perm.stderr, digits = digits), ")", sep = "")), 
        "\n")
  }
  if (!is.null(x$perm.conf.int)) {
    cat(format(100 * attr(x$perm.conf.int, "conf.level")), 
        " percent (Monte-Carlo) permutation percentile confidence interval:\n", 
        " ", paste(format(x$perm.conf.int[1:2], digits = digits), 
                   collapse = " "), "\n", sep = "")
  }
  cat("\nResults without permutation:\n")
  if (!is.null(x$statistic)) 
    out <- c(out, paste(names(x$statistic), "=", format(x$statistic, 
                                                        digits = max(1L, digits - 2L))))
  if (!is.null(x$parameter)) 
    out <- c(out, paste(names(x$parameter), "=", format(x$parameter, 
                                                        digits = max(1L, digits - 2L))))
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = max(1L, digits - 
                                                3L))
    out <- c(out, paste("p-value", 
                        if (substr(fp, 1L, 1L) == "<") fp else paste("=", fp)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(x$alternative)) {
    cat("alternative hypothesis: ")
    if (!is.null(x$null.value)) {
      if (length(x$null.value) == 1L) {
        alt.char <- switch(x$alternative, two.sided = "not equal to", 
                           less = "less than", greater = "greater than")
        cat("true ", names(x$null.value), " is ", alt.char, 
            " ", x$null.value, "\n", sep = "")
      }
      else {
        cat(x$alternative, "\nnull values:\n", sep = "")
        print(x$null.value, digits = digits, ...)
      }
    }
    else cat(x$alternative, "\n", sep = "")
  }
  if (!is.null(x$conf.int)) {
    cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval:\n", 
        " ", paste(format(x$conf.int[1:2], digits = digits), 
                   collapse = " "), "\n", sep = "")
  }
  if (!is.null(x$estimate)) {
    cat("sample estimates:\n")
    print(x$estimate, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}