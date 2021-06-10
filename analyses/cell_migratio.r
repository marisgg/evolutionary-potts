# devtools::install_github( "ingewortel/celltrackR" )
setwd("~/Uni/4/NC/evolutionary-potts/analyses/")
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
  x <- sapply(subtracks(track, 1), FUN = displacement)
  x <- data.frame(x)
  fit <- vglm(x ~ 1, levy(location = 0), data=x)
  return(fit)
}

tracklhist <- function(tracks, l) {
  hist(sapply(subtracks(tracks, l), FUN = trackLength), breaks = 50, 
       main =  paste("tracks of length", as.character(l), sep=" "),
       xlab="average track length")
}

stats <- function(tracks) {
  dat <- aggregate(tracks, squareDisplacement, FUN="mean.se", 
                   subtrack.length = seq(1, (maxTrackLength(tracks) / 2)))
  with( dat ,{
    plot( mean ~ i, xlab="step size",
          ylab="mean square displacement", type="l", main="Mean square displacement")
    segments( i[seq(1, length(i), 5)], lower[seq(1, length(i), 5)], y1=upper[seq(1, length(i), 5)] )
  } )
  
  # dat <- aggregate(tracks, overallNormDot, subtrack.length = seq(2, (maxTrackLength(tracks) - 1)))
  # plot(dat, type='l', main="autocorrelation for different step sizes", 
  #      xlab="step size", ylab="autocorrelation", ylim=c(-1,1))
  # 
  # dat <- aggregate(tracks, overallAngle, subtrack.length = seq(2, (maxTrackLength(tracks) - 1)))
  # plot(dat, type='l', main="mean turning angle for different step sizes", 
  #      xlab="step size", ylab="mean turning angle", ylim=c(0, pi))
  
  tracklhist(tracks, 1)
  tracklhist(tracks, 10)
  tracklhist(tracks, 30)
  tracklhist(tracks, 50)
}

others <- function() {
  par(mfcol=c(5,2))
  
  br <- as.tracks(list(cell = brownianTrack(nsteps = nrow(data), dim = 2, sd = trackLength(tracks[["cell"]])/nrow(data))))
  plot(br, main="brownian motion")

  stats(br)

  lv <- as.tracks(list(cell = levywalk(nsteps = nrow(data), dim = 2, s = -0.855579)))
  plot(lv, main="levy walk")

  stats(lv)
}

analysetype <- function(prefix, names) {
  data <- list()
  for (n in names) {
    print(paste(prefix, "/", n, ".csv", sep=""))
    d <- read.csv(paste(prefix, "/", n, ".csv", sep=""))
    data[[n]] <- as.matrix(d)
  }
  
  par(mfrow=c(3,2))
  
  tracks <- as.tracks(data)
  plot(tracks, main="Evolutionary cell", xlim=c(0,400), ylim=c(0,400), cex = 0, pch.start=NULL, pch.end=1)
  
  stats(tracks)
  return(tracks)
}

# analysetype("none", c("600_best", "600_0", "700_best", "700_1", "800_best", "800_0"))
# analysetype("medium", c("300_0", "300_45", "400_0", "400_2", "500_0", "500_46"))
# analysetype("high", c("10_0", "10_1", "100_0", "100_44", "200_17", "200_0"))
analysetype("high", c("10_0"))
analysetype("medium", c("300_45"))
analysetype("none", c("600_0"))
