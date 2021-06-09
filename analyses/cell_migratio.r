# devtools::install_github( "ingewortel/celltrackR" )
# setwd("~/Uni/4/NC/evolutionary-potts/analyses/")
library(celltrackR)
library(VGAM)

levywalk <- function (nsteps = 100, dim = 3, mean = .3, s = 1) 
{
  if (dim < 1) {
    stop("1 or more spatial dimensions required!")
  }
  d <- pmin(rlevy(nsteps, location = mean, scale = s/100), replicate(nsteps, 400))
  
  a <- runif(nsteps, 0, 2*pi)
  m <- cbind(0:nsteps, stats::diffinv(d * cos(a)), stats::diffinv(d * sin(a)))
  if (dim <= 3) {
    colnames(m) <- c("t", "x", "y", "z")[1:(dim + 1)]
  }
  else {
    colnames(m) <- c("t", paste0("x", 1:dim))
  }
  m
}

fitlevy <- function(track) {
  x <- sapply(subtracks(track, 1), FUN = trackLength)
  x <- data.frame(x)
  fit <- vglm(x ~ 1, levy(location = 0), data=x)
  return(fit)
}

stats <- function(tracks) {
  dat <- aggregate(tracks, squareDisplacement, FUN="mean.se", max.overlap = 0, 
                   subtrack.length = seq(1, (maxTrackLength(tracks) / 2)))
  with( dat ,{
    plot( mean ~ i, xlab="step size",
          ylab="mean square displacement", type="l", main="Mean square displacement")
    segments( i[seq(1, length(i), 5)], lower[seq(1, length(i), 5)], y1=upper[seq(1, length(i), 5)] )
  } )
  
  dat <- aggregate(tracks, overallNormDot, subtrack.length = seq(2, (maxTrackLength(tracks) - 1)))
  plot(dat, type='l', main="autocorrelation for different step sizes", xlab="step size", ylab="autocorrelation")
  
  dat <- aggregate(tracks, overallAngle, subtrack.length = seq(2, (maxTrackLength(tracks) - 1)))
  plot(dat, type='l', main="mean turning angle for different step sizes", xlab="step size", ylab="mean turning angle")
  
  fit <- fitlevy(tracks[["cell"]])
  hist(sapply(subtracks(tracks, 1), FUN = trackLength), breaks = 50, 
       main = paste("loglikelihood:", as.character(logLik(fit)), sep=" "),
       xlab="average track length")
}

data <- read.csv("data1.csv")
data[] <- lapply(data, as.numeric)

par(mfcol=c(5,3))

tracks <- as.tracks(list(cell = as.matrix(data)))
plot(tracks, main="Evolutionary cell")

stats(tracks)

br <- as.tracks(list(cell = brownianTrack(nsteps = nrow(data), dim = 2, sd = trackLength(tracks[["cell"]])/nrow(data))))
plot(br, main="brownian motion")

stats(br)

lv <- as.tracks(list(cell = levywalk(nsteps = nrow(data), dim = 2, s = -0.855579)))
plot(lv, main="levy walk")

stats(lv)