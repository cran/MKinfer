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
## bootstrap: assuming equal variances
normDiffCI(x, y, method = "classical", boot = TRUE, bootci.type = "bca")
## bootstrap: assuming unequal variances
normDiffCI(x, y, method = "welch", boot = TRUE, bootci.type = "bca")

## -----------------------------------------------------------------------------
M <- 100
CIhsu <- CIwelch <- CIclass <- matrix(NA, nrow = M, ncol = 2)
for(i in 1:M){
  x <- rnorm(10)
  y <- rnorm(30, sd = 0.1)
  CIclass[i,] <- normDiffCI(x, y, method = "classical")$conf.int
  CIwelch[i,] <- normDiffCI(x, y, method = "welch")$conf.int
  CIhsu[i,] <- normDiffCI(x, y, method = "hsu")$conf.int
}
## coverage probabilies
## classical
sum(CIclass[,1] < 0 & 0 < CIclass[,2])/M
## Welch
sum(CIwelch[,1] < 0 & 0 < CIwelch[,2])/M
## Hsu
sum(CIhsu[,1] < 0 & 0 < CIhsu[,2])/M

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
with(sleep, t.test(extra[group == 1], extra[group == 2]))
with(sleep, hsu.t.test(extra[group == 1], extra[group == 2]))
## Formula interface
t.test(extra ~ group, data = sleep)
hsu.t.test(extra ~ group, data = sleep)

## -----------------------------------------------------------------------------
boot.t.test(1:10, y = c(7:20)) # without bootstrap: P = .00001855
boot.t.test(1:10, y = c(7:20, 200)) # without bootstrap: P = .1245

## Traditional interface
with(sleep, boot.t.test(extra[group == 1], extra[group == 2]))
## Formula interface
boot.t.test(extra ~ group, data = sleep)

## -----------------------------------------------------------------------------
perm.t.test(1:10, y = c(7:20)) # without permutation: P = .00001855
## permutation confidence interval sensitive to outlier!
perm.t.test(1:10, y = c(7:20, 200)) # without permutation: P = .1245

## Traditional interface
with(sleep, perm.t.test(extra[group == 1], extra[group == 2]))
## Formula interface
res <- perm.t.test(extra ~ group, data = sleep)
res

## -----------------------------------------------------------------------------
p2ses(res$p.value)

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

## -----------------------------------------------------------------------------
set.seed(123)
outcome <- c(rnorm(10), rnorm(10, mean = 1.5), rnorm(10, mean = 1))
timepoints <- factor(rep(1:3, each = 10))
patients <- factor(rep(1:10, times = 3))
rm.oneway.test(outcome, timepoints, patients)
rm.oneway.test(outcome, timepoints, patients, method = "lme")
rm.oneway.test(outcome, timepoints, patients, method = "friedman")
rm.oneway.test(outcome, timepoints, patients, method = "quade")

## ---- fig.width=7, fig.height=7-----------------------------------------------
## Generate some data
x <- matrix(rnorm(1000, mean = 10), nrow = 10)
g1 <- rep("control", 10)
y1 <- matrix(rnorm(500, mean = 11.75), nrow = 10)
y2 <- matrix(rnorm(500, mean = 9.75, sd = 3), nrow = 10)
g2 <- rep("treatment", 10)
group <- factor(c(g1, g2))
Data <- rbind(x, cbind(y1, y2))
## compute Hsu t-test
pvals <- apply(Data, 2, function(x, group) hsu.t.test(x ~ group)$p.value,
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

## ---- fig.width=7, fig.height=7-----------------------------------------------
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

## ---- fig.width=7, fig.height=7-----------------------------------------------
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
                    method = "boot.t.test")
pairwise.ext.t.test(airquality$Ozone, airquality$Month, 
                    method = "perm.t.test")

## -----------------------------------------------------------------------------
sessionInfo()

