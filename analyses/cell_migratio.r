# devtools::install_github( "ingewortel/celltrackR" )
library(celltrackR)

stats <- function(tracks) {
  dat <- aggregate(tracks, squareDisplacement, FUN="mean.se")
  with( dat ,{
    plot( mean ~ i, xlab="time step",
          ylab="mean square displacement", type="l" )
    segments( i, lower, y1=upper )
  } )
  
  dat <- aggregate(tracks, overallNormDot)
  plot(dat, type='l')
  
  hist(sapply(subtracks(tracks, 1), FUN = trackLength), breaks = 50, main = "")
}

data <- read.csv("data.csv")
data[] <- lapply(data, as.numeric)

par(mfcol=c(4,2), mar=c(2, 3, 2, 2))

tracks <- as.tracks(list(cell = as.matrix(data)))
plot(tracks)

stats(tracks)

br <- as.tracks(list(cell = brownianTrack(nsteps = nrow(data), dim = 2, sd = trackLength(tracks[["cell"]])/nrow(data))))
plot(br)

stats(br)
