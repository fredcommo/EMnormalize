####################################
```r
# Example
op <- par(no.readonly=TRUE)

require("rGithubClient")
git <- getRepo('fredcommo/EMnormazlize')
sourceRepoFile(git, "EMnormalize.R")

set.seed(112335)
N <- 4
p.init <- runif(N-1, .1, .3)
p.init <- c(p.init, 1- sum(p.init))
x <- lapply(1:N, function(ii) rnorm(round(1e4*p.init[ii]), ii, .5) )
x <- do.call(c, x)

# Original EM
model <- Mclust(x)
m <- model$parameters$mean
s <- model$parameters$variance$sigmasq
p <- model$parameters$pro
D <- computeDensities(length(x), m, s, p)
bestPeak <- chooseBestPeak(D$peaks, m, peakThresh=.95)

par(mfrow=c(3, 1))
plotEMmodel(x, m, s, p, bestPeak, xlim=range(0, 6), main="Original Mclust")

# Resampling method
xnorm <- EMnormalize(x, xlim=range(0, 6), main="Resampling method")

# True groups
plotEMmodel(x, 1:4, rep(0.5, 4), p.init, which.max(p.init), xlim=range(0, 6), main="True distribution")
par(op)
```
