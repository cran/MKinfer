## -----------------------------------------------------------------------------
library(MKinfer)

## -----------------------------------------------------------------------------
## default: "wilson"
binomCI(x = 12, n = 50)
## Clopper-Pearson interval
binomCI(x = 12, n = 50, method = "clopper-pearson")
## identical to 
binom.test(x = 12, n = 50)$conf.int

## -----------------------------------------------------------------------------
## default: wilson 
binomDiffCI(a = 5, b = 0, c = 51, d = 29)
## default: wilson with continuity correction
binomDiffCI(a = 212, b = 144, c = 256, d = 707, paired = TRUE)

## -----------------------------------------------------------------------------
x <- rnorm(50, mean = 2, sd = 3)
## mean and SD unknown
normCI(x)
meanCI(x)
sdCI(x)
## SD known
normCI(x, sd = 3)
## mean known
normCI(x, mean = 2, alternative = "less")
## bootstrap
meanCI(x, boot = TRUE)

## -----------------------------------------------------------------------------
x <- rnorm(20)
y <- rnorm(20, sd = 2)
## paired
normDiffCI(x, y, paired = TRUE)
## compare
normCI(x-y)
## bootstrap
normDiffCI(x, y, paired = TRUE, boot = TRUE)

## unpaired
y <- rnorm(10, mean = 1, sd = 2)
## classical
normDiffCI(x, y, method = "classical")
## Welch (default as in case of function t.test)
normDiffCI(x, y, method = "welch")
## Hsu
normDiffCI(x, y, method = "hsu")
## Xiao
normDiffCI(x, y, method = "xiao")
## bootstrap: assuming equal variances
normDiffCI(x, y, method = "classical", boot = TRUE, bootci.type = "bca")
## bootstrap: assuming unequal variances
normDiffCI(x, y, method = "welch", boot = TRUE, bootci.type = "bca")

## -----------------------------------------------------------------------------
M <- 100
CIxiao <- CIhsu <- CIwelch <- CIclass <- matrix(NA, nrow = M, ncol = 2)
for(i in 1:M){
  x <- rnorm(10)
  y <- rnorm(30, sd = 0.1)
  CIclass[i,] <- normDiffCI(x, y, method = "classical")$conf.int
  CIwelch[i,] <- normDiffCI(x, y, method = "welch")$conf.int
  CIhsu[i,] <- normDiffCI(x, y, method = "hsu")$conf.int
  CIxiao[i,] <- normDiffCI(x, y, method = "xiao")$conf.int
}
## coverage probabilies
## classical
sum(CIclass[,1] < 0 & 0 < CIclass[,2])/M
## Welch
sum(CIwelch[,1] < 0 & 0 < CIwelch[,2])/M
## Hsu
sum(CIhsu[,1] < 0 & 0 < CIhsu[,2])/M
## Xiao
sum(CIxiao[,1] < 0 & 0 < CIxiao[,2])/M

## -----------------------------------------------------------------------------
x <- rnorm(100, mean = 10, sd = 2) # CV = 0.2
## default: "miller"
cvCI(x)
## Gulhar et al. (2012)
cvCI(x, method = "gulhar")
## bootstrap
cvCI(x, method = "boot")

## -----------------------------------------------------------------------------
x <- rexp(100, rate = 0.5)
## exact
quantileCI(x = x, prob = 0.95)
## asymptotic
quantileCI(x = x, prob = 0.95, method = "asymptotic")
## boot
quantileCI(x = x, prob = 0.95, method = "boot")

## -----------------------------------------------------------------------------
## exact
medianCI(x = x)
## asymptotic
medianCI(x = x, method = "asymptotic")
## boot
medianCI(x = x, method = "boot")

## -----------------------------------------------------------------------------
## exact
madCI(x = x)
## aysymptotic
madCI(x = x, method = "asymptotic")
## boot
madCI(x = x, method = "boot")
## unstandardized
madCI(x = x, constant = 1)

## -----------------------------------------------------------------------------
t.test(1:10, y = c(7:20))      # P = .00001855
t.test(1:10, y = c(7:20, 200)) # P = .1245    -- NOT significant anymore
hsu.t.test(1:10, y = c(7:20))
hsu.t.test(1:10, y = c(7:20, 200))

## Traditional interface
with(mtcars, t.test(mpg[am == 0], mpg[am == 1]))
with(mtcars, hsu.t.test(mpg[am == 0], mpg[am == 1]))
## Formula interface
t.test(mpg ~ am, data = mtcars)
hsu.t.test(mpg ~ am, data = mtcars)

## ----fig.width=7, fig.height=7------------------------------------------------
h0plot(t.test(mpg ~ am, data = mtcars))
h0plot(hsu.t.test(mpg ~ am, data = mtcars))

## ----fig.height = 7, fig.width = 7--------------------------------------------
## equal variances
dgt(0, n1 = 5, n2 = 10, v1tov2 = 1)
dt(0, df = 5+10-2)

pgt(2, n1 = 5, n2 = 10, v1tov2 = 1)
pt(2, df = 5+10-2)

qgt(0.975, n1 = 5, n2 = 10, v1tov2 = 1)
qt(0.975, df = 5+10-2)

library(ggplot2)
ggplot(data = data.frame(x = c(-4, 4)), aes(x)) +
  stat_function(fun = function(x) dgt(x, n1 = 5, n2 = 10, v1tov2 = 1), 
                n = 101, col = "darkblue", lwd = 1) + 
  stat_function(fun = function(x) dt(x, df = 5+10-2), 
                n = 101, col = "darkgreen", lwd = 1) + ylab("density") +
  ggtitle("Generalized Central (blue) t vs. Student t Distribution (green)")

## unequal variances
nu <- (1^2/5+5^2/10)^2/(1^4/(5^2*(5-1)) + 5^4/(10^2*(10-1)))
dgt(0, n1 = 5, n2 = 10, v1tov2 = 1^2/5^2)
dt(0, df = nu)

pgt(2, n1 = 5, n2 = 10, v1tov2 = 1^2/5^2)
pt(2, df = nu)

qgt(0.975, n1 = 5, n2 = 10, v1tov2 = 1^2/5^2)
qt(0.975, df = nu)

ggplot(data = data.frame(x = c(-4, 4)), aes(x)) +
  stat_function(fun = function(x) dgt(x, n1 = 5, n2 = 10, v1tov2 = 1^2/5^2), 
                n = 101, col = "darkblue", lwd = 1) + 
  stat_function(fun = function(x) dt(x, df = nu), 
                n = 101, col = "darkred", linetype = "dashed", lwd = 1) + 
  ylab("density") +
  ggtitle("Generalized Central t (blue) vs. Welch t Distribution (red)")

## -----------------------------------------------------------------------------
## Examples taken and adapted from function t.test
t.test(1:10, y = c(7:20))      # P = .00001855
t.test(1:10, y = c(7:20, 200)) # P = .1245    -- NOT significant anymore
xiao.t.test(1:10, y = c(7:20))
xiao.t.test(1:10, y = c(7:20, 200))

## Traditional interface
with(mtcars, t.test(mpg[am == 0], mpg[am == 1]))
with(mtcars, xiao.t.test(mpg[am == 0], mpg[am == 1]))
## Formula interface
t.test(mpg ~ am, data = mtcars)
xiao.t.test(mpg ~ am, data = mtcars)

## Example 6.1 in Xiao (2018)
x <- c(134, 146, 104, 119, 124, 161, 107, 83, 113, 129, 97, 123)
y <- c(70, 118, 101, 85, 107, 132, 94)
xiao.t.test(x, y, alternative = "greater")
t.test(x, y, var.equal = TRUE, alternative = "greater")
t.test(x, y, alternative = "greater")

## ----fig.width=7, fig.height=7------------------------------------------------
h0plot(t.test(mpg ~ am, data = mtcars))
h0plot(xiao.t.test(mpg ~ am, data = mtcars))

## -----------------------------------------------------------------------------
boot.t.test(1:10, y = c(7:20)) # without bootstrap: P = .00001855
boot.t.test(1:10, y = c(7:20, 200)) # without bootstrap: P = .1245

## Traditional interface
with(mtcars, boot.t.test(mpg[am == 0], mpg[am == 1]))
## Formula interface
boot.t.test(mpg ~ am, data = mtcars)

## ----fig.width=7, fig.height=7------------------------------------------------
h0plot(boot.t.test(mpg ~ am, data = mtcars, bootStat = TRUE))

## -----------------------------------------------------------------------------
perm.t.test(1:10, y = c(7:20), conf.type = "all") 
perm.t.test(1:10, y = c(7:20, 200), conf.type = "pivot") 
perm.t.test(1:10, y = c(7:20, 200), conf.type = "exact") 
perm.t.test(1:10, y = c(7:20, 200), conf.type = "stud") 
## contradiction between p value and confidence interval!
perm.t.test(1:10, y = c(7:20, 200), conf.type = "perc") 

## Traditional interface
with(mtcars, perm.t.test(mpg[am == 0], mpg[am == 1]))
## Formula interface
perm.t.test(mpg ~ am, data = mtcars)

## ----fig.width=7, fig.height=7------------------------------------------------
h0plot(perm.t.test(mpg ~ am, data = mtcars, permStat = TRUE))

## -----------------------------------------------------------------------------
## all permutations vs combinations
perm.t.test(1:3, y = 2:5, useCombn = TRUE)
perm.t.test(1:3, y = 2:5, useCombn = FALSE)

## -----------------------------------------------------------------------------
res <- perm.t.test(extra ~ group, data = sleep)
p2ses(res$p.value)

## -----------------------------------------------------------------------------
library(mvtnorm)
## effect size
delta <- c(0.25, 0.5)
## covariance matrix
Sigma <- matrix(c(1, 0.75, 0.75, 1), ncol = 2)
## sample size
n <- 50
## generate random data from multivariate normal distributions
X <- rmvnorm(n=n, mean = delta, sigma = Sigma)
Y <- rmvnorm(n=n, mean = rep(0, length(delta)), sigma = Sigma)
## perform multivariate z-test
mpe.z.test(X = X, Y = Y, Sigma = Sigma)
## perform multivariate t-test
mpe.t.test(X = X, Y = Y)

## -----------------------------------------------------------------------------
## Generate some data
set.seed(123)
x <- rnorm(25, mean = 1)
x[sample(1:25, 5)] <- NA
y <- rnorm(20, mean = -1)
y[sample(1:20, 4)] <- NA
pair <- c(rnorm(25, mean = 1), rnorm(20, mean = -1))
g <- factor(c(rep("yes", 25), rep("no", 20)))
D <- data.frame(ID = 1:45, response = c(x, y), pair = pair, group = g)

## Use Amelia to impute missing values
library(Amelia)
res <- amelia(D, m = 10, p2s = 0, idvars = "ID", noms = "group")

## Per protocol analysis (Welch two-sample t-test)
t.test(response ~ group, data = D)
## Intention to treat analysis (Multiple Imputation Welch two-sample t-test)
mi.t.test(res, x = "response", y = "group")

## Per protocol analysis (Two-sample t-test)
t.test(response ~ group, data = D, var.equal = TRUE)
## Intention to treat analysis (Multiple Imputation two-sample t-test)
mi.t.test(res, x = "response", y = "group", var.equal = TRUE)

## Specifying alternatives
mi.t.test(res, x = "response", y = "group", alternative = "less")
mi.t.test(res, x = "response", y = "group", alternative = "greater")

## One sample test
t.test(D$response[D$group == "yes"])
mi.t.test(res, x = "response", subset = D$group == "yes")
mi.t.test(res, x = "response", mu = -1, subset = D$group == "yes",
          alternative = "less")
mi.t.test(res, x = "response", mu = -1, subset = D$group == "yes",
          alternative = "greater")

## paired test
t.test(D$response, D$pair, paired = TRUE)
mi.t.test(res, x = "response", y = "pair", paired = TRUE)

## -----------------------------------------------------------------------------
library(mice)
res.mice <- mice(D, m = 10, print = FALSE)
mi.t.test(res.mice, x = "response", y = "group")

## -----------------------------------------------------------------------------
## Per protocol analysis (Exact Wilcoxon rank sum test)
library(exactRankTests)
wilcox.exact(response ~ group, data = D, conf.int = TRUE)
## Intention to treat analysis (Multiple Imputation Exact Wilcoxon rank sum test)
mi.wilcox.test(res, x = "response", y = "group")

## Specifying alternatives
mi.wilcox.test(res, x = "response", y = "group", alternative = "less")
mi.wilcox.test(res, x = "response", y = "group", alternative = "greater")

## One sample test
wilcox.exact(D$response[D$group == "yes"], conf.int = TRUE)
mi.wilcox.test(res, x = "response", subset = D$group == "yes")
mi.wilcox.test(res, x = "response", mu = -1, subset = D$group == "yes",
               alternative = "less")
mi.wilcox.test(res, x = "response", mu = -1, subset = D$group == "yes",
               alternative = "greater")

## paired test
wilcox.exact(D$response, D$pair, paired = TRUE, conf.int = TRUE)
mi.wilcox.test(res, x = "response", y = "pair", paired = TRUE)

## -----------------------------------------------------------------------------
mi.wilcox.test(res.mice, x = "response", y = "group")

## ----fig.width=7, fig.height=7------------------------------------------------
## small effect size
mdplot(delta = 0.2)
## medium effect size
mdplot(delta = 0.5)
## large effect size
mdplot(delta = 0.8)
## z-factor = 0  (z-factor = 1 - 3*2*sd/delta)
mdplot(delta = 6)
## z-factor = 0.5  (z-factor = 1 - 3*2*sd/delta)
mdplot(delta = 12)

## unequal variances
mdplot(delta = 0.8, sd1 = 1, sd2 = 2)
mdplot(delta = 0.8, sd1 = 2, sd2 = 1)

## ----fig.width=7, fig.height=7------------------------------------------------
library(ggplot2)
## (standardized) mean difference to sensitivity/specificity
## equal variances
delta <- seq(from = 0.0, to = 6, by = 0.05)
res <- sapply(delta, md2sens)
DF <- data.frame(SMD = delta, sensitivity = res[1,], 
                 specificity = res[2,])
ggplot(DF, aes(x = SMD, y = sensitivity)) +
  geom_line() + ylim(0.5, 1.0) + xlab("(standardized) mean difference") +
  ylab("sensitivity = specificity") + ggtitle("SD1 = SD2 = 1")

## unequal variances
delta <- seq(from = 0.0, to = 6, by = 0.05)
res <- sapply(delta, md2sens, sd1 = 1, sd2 = 2)
DF <- data.frame(MD = delta, performance = c(res[1,], res[2,]),
                 measure = c(rep("sensitivity", length(delta)),
                             rep("specificity", length(delta))))
ggplot(DF, aes(x = MD, y = performance, color = measure)) +
  geom_line() + ylim(0, 1.0) + xlab("mean difference") +
  scale_color_manual(values = c("darkblue", "darkred")) +
  ggtitle("SD1 = 1, SD2 = 2")

## ----fig.width=7, fig.height=7------------------------------------------------
## (standardized) mean difference to sensitivity/specificity
## equal variances
library(ggplot2)
delta <- seq(from = 2, to = 18, by = 0.05)
res <- sapply(delta, md2zfactor)
DF <- data.frame(SMD = delta, zfactor = res)
ggplot(DF, aes(x = SMD, y = zfactor)) +
  geom_line() + xlab("(standardized) mean difference") +
  ylab("z-factor") + ggtitle("SD1 = SD2 = 1") + 
  geom_hline(yintercept = 1, linetype = "dotted")

## unequal variances
delta <- seq(from = 2.5, to = 20, by = 0.05)
res <- sapply(delta, md2zfactor, sd1 = 1, sd2 = 2)
DF <- data.frame(MD = delta, zfactor = res)
ggplot(DF, aes(x = MD, y = zfactor)) +
  geom_line() + xlab("mean difference") +
  ylab("z-factor") + ggtitle("SD1 = 1, SD2 = 2") +
  geom_hline(yintercept = 1, linetype = "dotted")

## -----------------------------------------------------------------------------
set.seed(123)
outcome <- c(rnorm(10), rnorm(10, mean = 1.5), rnorm(10, mean = 1))
timepoints <- factor(rep(1:3, each = 10))
patients <- factor(rep(1:10, times = 3))
rm.oneway.test(outcome, timepoints, patients)
rm.oneway.test(outcome, timepoints, patients, method = "lme")
rm.oneway.test(outcome, timepoints, patients, method = "quade")
rm.oneway.test(outcome, timepoints, patients, method = "friedman")

## ----fig.width=7, fig.height=7------------------------------------------------
h0plot(rm.oneway.test(outcome, timepoints, patients))
h0plot(rm.oneway.test(outcome, timepoints, patients, method = "lme"))
h0plot(rm.oneway.test(outcome, timepoints, patients, method = "quade"))
h0plot(rm.oneway.test(outcome, timepoints, patients, method = "friedman"),
       qtail = 1e-4)

## ----fig.width=7, fig.height=7------------------------------------------------
## Generate some data
x <- matrix(rnorm(1000, mean = 10), nrow = 10)
g1 <- rep("control", 10)
y1 <- matrix(rnorm(500, mean = 11.75), nrow = 10)
y2 <- matrix(rnorm(500, mean = 9.75, sd = 3), nrow = 10)
g2 <- rep("treatment", 10)
group <- factor(c(g1, g2))
Data <- rbind(x, cbind(y1, y2))
## compute Xiao t-test
pvals <- apply(Data, 2, function(x, group) xiao.t.test(x ~ group)$p.value,
               group = group)
## compute log-fold change
logfc <- function(x, group){
  res <- tapply(x, group, mean)
  log2(res[1]/res[2])
}
lfcs <- apply(Data, 2, logfc, group = group)
volcano(lfcs, p.adjust(pvals, method = "fdr"), 
        effect.low = -0.25, effect.high = 0.25, 
        xlab = "log-fold change", ylab = "-log10(adj. p value)")

## ----fig.width=7, fig.height=7------------------------------------------------
data("fingsys")
baplot(fingsys$fingsys, fingsys$armsys, 
       title = "Approximative Confidence Intervals", 
       type = "parametric", ci.type = "approximate")
baplot(fingsys$fingsys, fingsys$armsys, 
       title = "Exact Confidence Intervals", 
       type = "parametric", ci.type = "exact")
baplot(fingsys$fingsys, fingsys$armsys, 
       title = "Bootstrap Confidence Intervals", 
       type = "parametric", ci.type = "boot", R = 999)

## ----fig.width=7, fig.height=7------------------------------------------------
data("fingsys")
baplot(fingsys$fingsys, fingsys$armsys, 
       title = "Approximative Confidence Intervals", 
       type = "nonparametric", ci.type = "approximate")
baplot(fingsys$fingsys, fingsys$armsys, 
       title = "Exact Confidence Intervals", 
       type = "nonparametric", ci.type = "exact")
baplot(fingsys$fingsys, fingsys$armsys, 
       title = "Bootstrap Confidence Intervals", 
       type = "nonparametric", ci.type = "boot", R = 999)

## -----------------------------------------------------------------------------
SD1 <- c(0.149, 0.022, 0.036, 0.085, 0.125, NA, 0.139, 0.124, 0.038)
SD2 <- c(NA, 0.039, 0.038, 0.087, 0.125, NA, 0.135, 0.126, 0.038)
SDchange <- c(NA, NA, NA, 0.026, 0.058, NA, NA, NA, NA)
imputeSD(SD1, SD2, SDchange)

## -----------------------------------------------------------------------------
SDchange2 <- rep(NA, 9)
imputeSD(SD1, SD2, SDchange2, corr = c(0.85, 0.9, 0.95))

## -----------------------------------------------------------------------------
pairwise.wilcox.test(airquality$Ozone, airquality$Month, 
                     p.adjust.method = "none")
## To avoid the warnings
pairwise.fun(airquality$Ozone, airquality$Month, 
             fun = function(x, y) wilcox.exact(x, y)$p.value)

## -----------------------------------------------------------------------------
pairwise.wilcox.exact(airquality$Ozone, airquality$Month)

## -----------------------------------------------------------------------------
pairwise.t.test(airquality$Ozone, airquality$Month, pool.sd = FALSE)
pairwise.ext.t.test(airquality$Ozone, airquality$Month)
pairwise.ext.t.test(airquality$Ozone, airquality$Month,
                    method = "hsu.t.test")
pairwise.ext.t.test(airquality$Ozone, airquality$Month,
                    method = "xiao.t.test")
pairwise.ext.t.test(airquality$Ozone, airquality$Month, 
                    method = "boot.t.test")
pairwise.ext.t.test(airquality$Ozone, airquality$Month, 
                    method = "perm.t.test")

## -----------------------------------------------------------------------------
sessionInfo()

