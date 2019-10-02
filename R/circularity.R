## Circularity - this takes a long time, so don't run again unless changes to data.

data$circularity <- NA
outlines <- dir("data/outlines")

for (i in 1:nrow(data)) {

  fil <- outlines[grepl(paste0(data$species_code[i], sprintf("%04d", data$image_tag.x[i]), "_", data$fieldtrip_id.x[i]), outlines)][1]
  ol <- read.delim(paste0("data/outlines/", fil), header=FALSE)
  sps <- SpatialPolygons(list(Polygons(list(Polygon(cbind(ol[,1], ol[,2]))),1)))

  plot(sps)
  sps2 <- gConvexHull(sps) #*(dim(img)[1]/2000)
  lines(sps2, col="blue", lwd=2)
  circularity(area(sps2), polyPerimeter(sps2))
  circularity(area(sps), polyPerimeter(sps))

  data$circularity[i] <- circularity(area(sps2), polyPerimeter(sps2))
}
