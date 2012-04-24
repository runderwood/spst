library(maptools)
library(spdep)

args <- commandArgs(trailingOnly=TRUE)

infile <- args[1]
infile

pts <- readShapePoints(infile)
pts.nb <- dnearneigh(pts, 0, 1000000000)
dists <- nbdists(pts.nb, pts)
idw <- lapply(dists, function(x) 1/x)
idw.list <- nb2listw(pts.nb, glist=idw, style="W")
moran.test(pts$sample, listw=idw.list)
