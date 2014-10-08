##############################
# Helper functions
require(mclust)
.smooth <- function(x, cut, K){
  '
	Called by EMnormalize(object)
	Smoothing the LR vector to improve the EM classification.
	'
  if(any(is.na(x)))
    x <- x[!is.na(x)]
  runx <- runmed(x, k = K)	
  q1 <- cut[1]; q2 <- cut[2]
  index <- which(runx>=q1 & runx<=q2)
  return (runx[index])
}
buildEMmodel <- function(x, G, useN){
  '
	Called by EMnormalize(object)
	Model the distribution as a gaussian mixture.
	'
  if(is.na(useN))
    useN <- floor(length(x)*.25)
  model <- Mclust(x[sample(1:length(x), size=useN)], G=G)
  nG <- model$G
  p <- model$parameters$pro
  m <- model$parameters$mean
  s <- model$parameters$variance$sigmasq
  if(length(s)<length(m)) s <- rep(s, length(m))
  p <- p[order(m)]
  s <- s[order(m)]
  m <- m[order(m)]
  return(list(nG = nG, m = m, p = p, s =s))
}
computeDensities <- function(n, m, s, p){
  '
	Called in EMnormalize(object)
	Simulates the mixture model according to the returned EM paramaters.
	'
  if(length(s)<length(m)) s <- rep(s, length(m))
  D <- lapply(1:length(m), function(ii){
    tmp <- rnorm(n, m[ii], sqrt(s[ii]))
    tmpD <- density(tmp, na.rm = T)
    tmpD$y = tmpD$y *p[ii]
    return(tmpD)
  })
  peaks <- sapply(D, function(d) max(d$y, na.rm=TRUE))
  return(list(D=D, peaks=peaks))
}
chooseBestPeak <- function(peaks, m, peakThresh){
  '
	Called in EMnormalize(object)
  Estimates what peak as to be used as the centralization value.
  '
  best <- which(peaks>=max(peaks)*peakThresh)
  if(length(best)>0){
    cat(length(best), 'peak(s) above', sprintf("%s%s", peakThresh*100, "%"), 'of max peak.\n')
  }
  else{
    cat('No peak above threshold:', peakThresh, '\n')
  }
  bestPeak <- best[1]
  cat('Peak at', m[bestPeak], 'has been chosen.\n')
  return(bestPeak)
}
plotEMmodel <- function(x, m, s, p, bestPeak, ...){
  
'
	Called in EMnormalize(object)
	Visualization of the mixture model
	'
  if(length(s)<length(m)) s <- rep(s, length(m))
  dx <- density(x, na.rm=TRUE)
  plot(dx, ...)
  polygon(dx$x, dx$y, col="grey90")

  N <- length(m)
  maxy <- max(dx$y)

  lapply(1:N, function(ii){
    x <- sort(rnorm(1000, m[ii], sqrt(s[ii])))
    d <- dnorm(x, m[ii], sqrt(s[ii]))
    d <- d/max(d)*p[ii]
    lines(x, d, col = rgb(ii/N, 0.2, (N-ii)/N, 0.75))
    alpha=.1; font=1; Cex=1
    if(ii==bestPeak){
      alpha=.5; font=2; Cex=1.25
    } 
    polygon(x, d, col = rgb(ii/N, 0.2, (N-ii)/N, alpha))
    text(m[ii], min(maxy, max(d)*1.25), round(m[ii], 3), font=font, cex=Cex)
  })
}

# End helper functions
##############################
##############################
# Main function
EMnormalize <- function(x, cut=c(.05, .95), G=3:6, B=100, peakThresh=0.95, ksmooth=11, useN=1e3, Plot=TRUE, ...){
  cuts <- as.numeric(quantile(x, c(cut[1], cut[2])))
  runx <- .smooth(x, cuts, K=ksmooth)
  
  cat("Searching for parameters...\n")
  EMmodels <- lapply(1:B, function(b){cat(b, "\t"); buildEMmodel(runx, G, useN)})										# helper function
  nG <- do.call(c, lapply(EMmodels, function(m) m$nG))
  models <- EMmodels[nG==median(nG)]
  m <- do.call(rbind, lapply(models, function(m) m$m))
  p <- do.call(rbind, lapply(models, function(m) m$p))
  s <- do.call(rbind, lapply(models, function(m) m$s))
  cat('\nDone.\n')
  
  m <- apply(m, 2, median, na.rm=TRUE)
  s <- apply(s, 2, median, na.rm=TRUE)
  p <- apply(p, 2, median, na.rm=TRUE)
  
  # compute densities
  cat('Computing densities...')
  computeD <- computeDensities(length(runx), m, s, p)								# helper function
  bestPeak <- chooseBestPeak(computeD$peaks, m, peakThresh)
  correct <- m[bestPeak]

  cat('\tDone.\n')
  cat('Gaussian mixture:')
  cat("\nn.peaks = ", median(nG, na.rm=TRUE), '\n') 
  cat("\n\nEM correction factor = ", correct, "\n\n")

  if(Plot){
    plotEMmodel(runx, m, s, p, bestPeak, ...)
  }
return(list(original=x, corrected=x-correct, means=m, sd=s, props=p, correction=correct))
}


