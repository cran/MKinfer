sample.perm <- function(x, k, nx = NULL, R, replace = FALSE, useCombn = FALSE){
  if(is.null(k)) k <- length(x)
  
  if(useCombn){
    max.R <- choose(k, nx)
  }else{
    max.R <- try(npermutations(x, k = k, replace = replace), silent = TRUE)
  }
  if(inherits(max.R, "try-error")) max.R <- Inf
  if(max.R < R){
    if(useCombn){
      message("The requested number of combinations (", R, ") is larger than ",
              "the total number of possible combinations (", max.R, ").\n",
              "Hence all possible combinations are computed.")
    }else{
      message("The requested number of permutations (", R, ") is larger than ",
              "the total number of possible permutations (", max.R, ").\n",
              "Hence all possible permutations are computed.")
    }
    if(useCombn){
      res1 <- combinations(x, k = nx, replace = replace)
      res2 <- combinations(x, k = k-nx, replace = replace)
      res <- cbind(res1, res2[nrow(res2):1,])
    }else{
      res <- permutations(x, k = k, replace = replace)
    }
  }else{
    res <- permutations(x, k = k, replace = replace, nsample = R)
  }
  res
}
perm.t.test <- function (x, ...){ 
  UseMethod("perm.t.test")
}
perm.t.test.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
                        mu = 0, paired = FALSE, var.equal = FALSE, 
                        conf.type = "pivot", conf.level = 0.95, R = 9999, 
                        symmetric = TRUE, permStat = FALSE, useCombn = FALSE, ...){
  alternative <- match.arg(alternative)
  if(!missing(mu) && (length(mu) != 1 || is.na(mu))) 
    stop("'mu' must be a single number")
  if(conf.type %notin% c("pivot", "exact", "stud", "perc", "all")){
    stop("'conf.type' must be one of 'pivot', 'exact', 'stud', 'perc', 'all'")
  }
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
    stddev <- sqrt(vx)
    if (stderr < 10 * .Machine$double.eps * abs(mx)) 
      stop("data are essentially constant")
    tstat <- (mx - mu)/stderr
    method <- if (paired) "Permutation Paired t-test" else "Permutation One Sample t-test"
    estimate <- setNames(mx, if (paired) "mean of the differences" else "mean of x")
    eff <- mx
    x.cent <- x - mx
    X <- abs(x.cent)*sample.perm(c(-1,1), k = nx, R = R, replace = TRUE)
    R.true <- nrow(X)
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
    eff <- mx-my
    names(estimate) <- c("mean of x", "mean of y")
    z <- c(x, y)
    Z.perm <- sample.perm(seq_along(z), k = nx+ny, nx = nx, R = R, useCombn = useCombn)
    R.true <- nrow(Z.perm)
    Z <- matrix(z[Z.perm], nrow = R.true, ncol = nx + ny)
    R.true <- nrow(Z)
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
    stddev <- sqrt(vx + vy)
    perm.stderr <- mean(STDERR)
    perm.estimate <- mean(EFF) 
    names(perm.estimate) <- "permutation difference of means"
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))) 
      stop("data are essentially constant")
    tstat <- (mx - my - mu)/stderr
    TSTAT <- (MX - MY)/STDERR
  }
  ## pivot inversion
  get.p.pivot <- function(mu.cand, target.type) {
    t.obs.s <- (eff - mu.cand) / stderr
    if (target.type == "less") {
      return((sum(TSTAT <= t.obs.s) + 1) / (R.true + 1))
    } else if (target.type == "greater") {
      return((sum(TSTAT >= t.obs.s) + 1) / (R.true + 1))
    } else {
      return((sum(abs(TSTAT) >= abs(t.obs.s)) + 1) / (R.true + 1))
    }
  }
  get.p.exact <- function(mu.cand, target.type) {
    if (is.null(y)) {
      t.obs.s <- (eff - mu.cand) / stderr
      t.perm.s <- (MX + (mx - mu.cand)) / STDERR
    } else {
      t.obs.s <- (eff - mu.cand) / stderr
      z.s <- c(x, y + mu.cand)
      Z.s <- matrix(z.s[Z.perm], nrow = R.true, ncol = nx + ny)
      X.s <- Z.s[, 1:nx]
      Y.s <- Z.s[, (nx+1):(nx+ny)]
      MX.s <- rowMeans(X.s)
      MY.s <- rowMeans(Y.s)
      
      if (var.equal) {
        V.s <- (rowSums((X.s - MX.s)^2) + rowSums((Y.s - MY.s)^2)) / df
        STDERR.s <- sqrt(V.s * (1/nx + 1/ny))
      } else {
        VX.s <- rowSums((X.s - MX.s)^2) / (nx - 1)
        VY.s <- rowSums((Y.s - MY.s)^2) / (ny - 1)
        STDERR.s <- sqrt(VX.s/nx + VY.s/ny)
      }
      t.perm.s <- (MX.s - MY.s) / STDERR.s
    }
    
    if (target.type == "less") {
      return((sum(t.perm.s <= t.obs.s) + 1) / (R.true + 1))
    } else if (target.type == "greater") {
      return((sum(t.perm.s >= t.obs.s) + 1) / (R.true + 1))
    } else {
      return((sum(abs(t.perm.s) >= abs(t.obs.s)) + 1) / (R.true + 1))
    }
  }
  ## robust inversion with grid search and bisection
  find.pinv.bound <- function(target.type, target.p, p.func, side = c("lower", "upper"), tol = 1e-8) {
    side <- match.arg(side)
    
    mult <- 15
    repeat {
      grid.vals <- seq(eff - mult * stderr, eff + mult * stderr, length.out = 150)
      p.vals <- sapply(grid.vals, function(m) p.func(m, target.type))
      in.ci <- grid.vals[p.vals > target.p]
      
      if (length(in.ci) > 0 || mult > 500) break
      mult <- mult * 3
    }
    
    if (length(in.ci) == 0) return(if(side == "lower") -Inf else Inf)
    
    if (side == "lower") {
      cand <- min(in.ci)
      step <- (grid.vals[2] - grid.vals[1])
      a <- cand - step
      b <- cand + step
    } else {
      cand <- max(in.ci)
      step <- (grid.vals[2] - grid.vals[1])
      a <- cand - step
      b <- cand + step
    }
    
    mid <- (a + b) / 2
    p.mid <- p.func(mid, target.type)
    while ((b - a) > tol && abs(p.mid - target.p) > 1/R.true) {
      if (side == "lower") {
        if (p.mid <= target.p) a <- mid else b <- mid
      } else {
        if (p.mid > target.p) a <- mid else b <- mid
      }
      mid <- (a + b) / 2
      p.mid <- p.func(mid, target.type)
    }
    return((a + b) / 2)
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    perm.pval <- max(mean(TSTAT <= tstat), 1/R.true)
    cint <- c(-Inf, tstat + qt(conf.level, df))
    ## confidence interval
    if(conf.type == "pivot"){
      ## pivot inversion
      u.bound <- find.pinv.bound("less", 1 - conf.level, get.p.pivot, side = "upper")
      perm.cint <- c(-Inf, u.bound)
    }
    if(conf.type == "exact"){
      ## exact inversion
      u.bound.exact <- find.pinv.bound("less", 1 - conf.level, get.p.exact, 
                                       side = "upper")
      perm.cint <- c(-Inf, u.bound.exact)
    }
    if(conf.type == "stud"){
      ## studentized
      perm.cint <- c(-Inf, eff - quantile(TSTAT, 1-conf.level)*stderr)
    }
    if(conf.type == "perc"){
      ## percentile
      perm.cint <- c(-Inf, quantile(EFF, conf.level))
    }
    if(conf.type == "all"){
      ## pivot inversion
      u.bound <- find.pinv.bound("less", 1 - conf.level, get.p.pivot, side = "upper")
      perm.cint.pivot <- c(-Inf, u.bound)
      u.bound.exact <- find.pinv.bound("less", 1 - conf.level, get.p.exact, 
                                       side = "upper")
      perm.cint.exact <- c(-Inf, u.bound.exact)
      ## studentized
      perm.cint.stud <- c(-Inf, eff - quantile(TSTAT, 1-conf.level)*stderr)
      ## percentile
      perm.cint.perc <- c(-Inf, quantile(EFF, conf.level))
      perm.cint <- rbind(perm.cint.pivot, perm.cint.exact, 
                         perm.cint.stud, perm.cint.perc)
      rownames(perm.cint) <- c("pivot", "exact", "stud", "perc")
    }
  }else if(alternative == "greater") {
    perm.pval <- max(mean(TSTAT >= tstat), 1/R.true)
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
    ## confidence interval
    if(conf.type == "pivot"){
      ## pivot inversion
      l.bound <- find.pinv.bound("greater", 1 - conf.level, get.p.pivot, side = "lower")
      perm.cint <- c(l.bound, Inf)
    }
    if(conf.type == "exact"){
      ## exact inversion
      l.bound.exact <- find.pinv.bound("greater", 1 - conf.level, get.p.exact, 
                                       side = "lower")
      perm.cint <- c(l.bound.exact, Inf)
    }
    if(conf.type == "stud"){
      ## studentized
      perm.cint <- c(eff - quantile(TSTAT, conf.level)*stderr, Inf)
    }
    if(conf.type == "perc"){
      ## percentile
      perm.cint <- c(quantile(EFF, 1-conf.level), Inf)
    }
    if(conf.type == "all"){
      ## pivot inversion
      l.bound <- find.pinv.bound("greater", 1 - conf.level, get.p.pivot, side = "lower")
      perm.cint.pivot <- c(l.bound, Inf)
      ## exact inversion
      l.bound.exact <- find.pinv.bound("greater", 1 - conf.level, get.p.exact, 
                                       side = "lower")
      perm.cint.exact <- c(l.bound.exact, Inf)
      ## studentized
      perm.cint.stud <- c(eff - quantile(TSTAT, conf.level)*stderr, Inf)
      ## percentile
      perm.cint.perc <- c(quantile(EFF, 1-conf.level), Inf)
      perm.cint <- rbind(perm.cint.pivot, perm.cint.exact, 
                         perm.cint.stud, perm.cint.perc)
      rownames(perm.cint) <- c("pivot", "exact", "stud", "perc")
    }
  }else{
    pval <- 2 * pt(-abs(tstat), df)
    if(symmetric)
      perm.pval <- max(mean(abs(TSTAT) >= abs(tstat)), 1/R.true)
    else
      perm.pval <- max(2*min(mean(TSTAT <= tstat), mean(TSTAT > tstat)), 1/R.true)
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
    ## confidence interval
    if(conf.type == "pivot"){
      ## pivot inversion
      l.bound <- find.pinv.bound("two.sided", alpha, get.p.pivot, side = "lower")
      u.bound <- find.pinv.bound("two.sided", alpha, get.p.pivot, side = "upper")
      perm.cint <- c(l.bound, u.bound)
    }
    if(conf.type == "exact"){
      ## exact inversion
      l.bound.exact <- find.pinv.bound("two.sided", alpha, get.p.exact, side = "lower")
      u.bound.exact <- find.pinv.bound("two.sided", alpha, get.p.exact, side = "upper")
      perm.cint <- c(l.bound.exact, u.bound.exact)
    }
    if(conf.type == "stud"){
      ## studentized
      perm.cint <- eff - quantile(TSTAT, c(1-alpha/2, alpha/2))*stderr
    }
    if(conf.type == "perc"){
      ## percentile
      perm.cint <- quantile(EFF, c(alpha/2, 1-alpha/2))
    }
    if(conf.type == "all"){
      ## pivot inversion
      l.bound <- find.pinv.bound("two.sided", alpha, get.p.pivot, side = "lower")
      u.bound <- find.pinv.bound("two.sided", alpha, get.p.pivot, side = "upper")
      perm.cint.pivot <- c(l.bound, u.bound)
      ## exact inversion
      l.bound.exact <- find.pinv.bound("two.sided", alpha, get.p.exact, side = "lower")
      u.bound.exact <- find.pinv.bound("two.sided", alpha, get.p.exact, side = "upper")
      perm.cint.exact <- c(l.bound.exact, u.bound.exact)
      ## studentized
      perm.cint.stud <- eff - quantile(TSTAT, c(1-alpha/2, alpha/2))*stderr
      ## percentile
      perm.cint.perc <- quantile(EFF, c(alpha/2, 1-alpha/2))
      perm.cint <- rbind(perm.cint.pivot, perm.cint.exact, 
                         perm.cint.stud, perm.cint.perc)
      rownames(perm.cint) <- c("pivot", "exact", "stud", "perc")
    }
  }
  cint <- mu + cint * stderr
  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- if (paired || !is.null(y)) "difference in means" else "mean"
  attr(cint, "conf.level") <- conf.level
  attr(perm.cint, "conf.level") <- conf.level
  if(permStat){ 
    perm.statistic <- TSTAT
  }else{
    perm.statistic <- NULL
  }
  rval <- list(statistic = tstat, parameter = df, p.value = pval, 
               perm.p.value = perm.pval, R = R, R.true = R.true, 
               p.min = perm.pval == 1/R.true,
               conf.int = cint, conf.type = conf.type, perm.conf.int = perm.cint,
               estimate = estimate, perm.estimate = perm.estimate, 
               null.value = mu, stderr = stderr, perm.stderr = perm.stderr,
               alternative = alternative, method = method, data.name = dname,
               perm.statistic = perm.statistic)
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
  cat("number of permutations:  ", x$R.true, "\n", sep = "")
  out <- character()
  if (!is.null(x$perm.p.value)) {
    bfp <- format.pval(x$perm.p.value, digits = max(1L, digits - 3L))
    if(x$R.true < x$R){
      if(x$p.min){
        cat("(Exact) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("<", bfp), "\n")
      }else{
        cat("(Exact) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("=", bfp), "\n")
      }
    }else{
      if(x$p.min){
        cat("(Monte-Carlo) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("<", bfp), "\n")
      }else{
        cat("(Monte-Carlo) permutation p-value", if (substr(bfp, 1L, 1L) == "<") bfp else paste("=", bfp), "\n")
      }
    }
  }
  if (!is.null(x$perm.estimate)) {
    cat(paste(names(x$perm.estimate), "(SE) =", 
              format(x$perm.estimate, digits = digits),
              paste("(", format(x$perm.stderr, digits = digits), ")", sep = "")), 
        "\n")
  }
  if (!is.null(x$perm.conf.int)) {
    if(x$R.true < x$R){
      if(x$conf.type == "all"){
        cat(format(100 * attr(x$perm.conf.int, "conf.level")), 
            " percent (exact) permutation confidence interval:\n", 
            "pivot:\t", paste(format(x$perm.conf.int[1,1:2], digits = digits), 
                                      collapse = " "), "\n", 
            "exact:\t", paste(format(x$perm.conf.int[2,1:2], digits = digits), 
                              collapse = " "), "\n", 
            "stud:\t", paste(format(x$perm.conf.int[3,1:2], digits = digits), 
                              collapse = " "), "\n", 
            "perc:\t", paste(format(x$perm.conf.int[4, 1:2], digits = digits), 
                             collapse = " "), "\n", sep = "")
      }else{
        cat(format(100 * attr(x$perm.conf.int, "conf.level")), 
            " percent (exact) permutation confidence interval:\n", 
            x$conf.type, ":\t", paste(format(x$perm.conf.int[1:2], digits = digits), 
                                      collapse = " "), "\n", sep = "")
      }
    }else{
      if(x$conf.type == "all"){
        cat(format(100 * attr(x$perm.conf.int, "conf.level")), 
            " percent (Monte-Carlo) permutation confidence interval:\n", 
            "pivot:\t", paste(format(x$perm.conf.int[1,1:2], digits = digits), 
                              collapse = " "), "\n", 
            "exact:\t", paste(format(x$perm.conf.int[2,1:2], digits = digits), 
                              collapse = " "), "\n", 
            "stud:\t", paste(format(x$perm.conf.int[3,1:2], digits = digits), 
                             collapse = " "), "\n", 
            "perc:\t", paste(format(x$perm.conf.int[4, 1:2], digits = digits), 
                             collapse = " "), "\n", sep = "")
      }else{
        cat(format(100 * attr(x$perm.conf.int, "conf.level")), 
            " percent (Monte-Carlo) permutation confidence interval:\n", 
            x$conf.type, ":\t", paste(format(x$perm.conf.int[1:2], digits = digits), 
                                      collapse = " "), "\n", sep = "")
      }
    }
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
