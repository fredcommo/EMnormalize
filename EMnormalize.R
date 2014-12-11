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
  idx <- sample(1:length(x), size=useN)
  model <- Mclust(x[idx], G=G)
  nG <- model$G
  k <- as.factor(model$classification)
  p <- model$parameters$pro
  m <- model$parameters$mean
  s2 <- model$parameters$variance$sigmasq
  if(length(s2)<length(m)) s2 <- rep(s2, length(m))
  levels(k) <- order(m)
  p <- p[order(m)]
  s2 <- s2[order(m)]
  m <- m[order(m)]
  return(list(nG = nG, m = m, p = p, s2 = s2, idx=idx, k=as.numeric(as.vector(k))))
}
computeDensities <- function(n, m, s2, p){
  '
	Called in EMnormalize(object)
	Simulates the mixture model according to the returned EM paramaters.
	'
  if(length(s2)<length(m)) s2 <- rep(s2, length(m))
  D <- lapply(1:length(m), function(ii){
    tmp <- rnorm(n, m[ii], sqrt(s2[ii]))
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
plotEMmodel <- function(x, m, s2, p, bestPeak, ...){
  
  '
	Visualization of the mixture model
	Called in EMnormalize(object)
	'
  
  if(length(s2)<length(m)) s2 <- rep(s2, length(m))
  dx <- density(x, na.rm=TRUE)
  bw <- dx$bw/10
  h <- hist(x, nclass=floor(length(x)*bw), freq=FALSE, border="lightblue",...)
  lines(dx, lty=2)
  polygon(dx$x, dx$y, col=rgb(.5, .5, .5, .05), border=NA)

  N <- length(m)
  maxy <- max(dx$y)*1.25

  lapply(1:N, function(ii){
    x <- sort(rnorm(1e4, m[ii], sqrt(s2[ii])))
    d <- dnorm(x, m[ii], sqrt(s2[ii]))
    d <- d*p[ii]
    lines(x, d, col = rgb(ii/N, 0.2, (N-ii)/N, 0.75))
    alpha=.1; font=1; Cex=1
    if(ii==bestPeak){
      alpha=.5; font=2; Cex=1.25
    } 
    polygon(x, d, col = rgb(ii/N, 0.2, (N-ii)/N, alpha))
    yrange <- par("yaxp")[1:2]
    text(m[ii], max(d)+yrange[2]*.1, round(m[ii], 3), font=font, cex=Cex)
    
#    text(m[ii], maxy*p[ii], round(m[ii], 3), font=font, cex=Cex)
    #    text(m[ii], min(maxy, max(d)*1.35), round(m[ii], 3), font=font, cex=Cex)
  })
}

# End helper functions
##############################
##############################
# Main function
EMnormalize <- function(x, cut=c(0, 1), G=3:6, B=100, peakThresh=0.95, Ksmooth=11, useN=1e3, Plot=TRUE, ...){

  cuts <- as.numeric(quantile(x, c(cut[1], cut[2]), na.rm=TRUE))

  if(Ksmooth>0)
    xprim <- .smooth(x, cuts, K=Ksmooth)
  else
    xprim <- x
  
  if(useN>length(xprim)) useN <- floor(length(xprim)*.3)
  
  cat("Searching for parameters...\n")
  EMmodels <- lapply(1:B, function(b){cat(b, "\t"); buildEMmodel(xprim, G, useN)})										# helper function
  nG <- do.call(c, lapply(EMmodels, function(m) m$nG))
  models <- EMmodels[nG==median(nG)]
  m <- do.call(rbind, lapply(models, function(m) m$m))
  p <- do.call(rbind, lapply(models, function(m) m$p))
  s2 <- do.call(rbind, lapply(models, function(m) m$s2))
  idx <- do.call(c, lapply(models, function(m) m$idx))
  k <- do.call(c, lapply(models, function(m) m$k))
  tk <- table(idx, k)
  k <- apply(tk, 1, function(x) colnames(tk)[which.max(x)])
  cat('\nDone.\n')
  
  m <- apply(m, 2, median, na.rm=TRUE)
  s2 <- apply(s2, 2, median, na.rm=TRUE)
  if(Ksmooth>0)
    s2 <- s2*Ksmooth
  p <- apply(p, 2, median, na.rm=TRUE)
  
  # compute densities
  cat('Computing densities...')
  computeD <- computeDensities(length(xprim), m, s2, p)								# helper function
  bestPeak <- Inf
  correct <- 0
  if(!is.na(peakThresh)){
    bestPeak <- chooseBestPeak(computeD$peaks, m, peakThresh)
    correct <- m[bestPeak]
  }
  
  cat('\tDone.\n')
  cat('Gaussian mixture:')
  cat("\nn.peaks = ", median(nG, na.rm=TRUE), '\n') 
  cat("\n\nEM correction factor = ", correct, "\n\n")

  if(Plot){
    plotEMmodel(x, m, s2, p, bestPeak, ...)
  }
return(list(original=x, corrected=x-correct, means=m, sq=s2, props=p, correction=correct, k=as.numeric(k)))
}


