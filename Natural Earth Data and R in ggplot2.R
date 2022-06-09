# Natural Earth Data and R in ggplot2

library(raster)
library(rgdal)
library(ggplot2)
library(reshape2)
library(plyr)

#  Assuming you have a path 'Maps' that you store your spatial files in.  This
#  is all downloaded from <a href="http://www.naturalearthdata.com/downloads/">http://www.naturalearthdata.com/downloads/</a> using the
#  1:50m "Medium" scale data.

nat.earth <- stack('G:/GIS/natural_earth/NE2_HR_LC_SR_W/NE2_HR_LC_SR_W.tif')
ne_lakes <- readOGR('G:/GIS/natural_earth//10m_physical/ne_10m_lakes.shp', 'ne_10m_lakes')
ne_rivers <- readOGR('G:/GIS/natural_earth/10m_physical/ne_10m_rivers_lake_centerlines.shp', 'ne_10m_rivers_lake_centerlines')
ne_coast <- readOGR('G:/GIS/natural_earth/10m_physical/ne_10m_coastline.shp', 'ne_10m_coastline')

#  I have a domain I'm interested in, but there's no reason you can't define something else:
quick.subset <- function(x, longlat){
  
  # longlat should be a vector of four values: c(xmin, xmax, ymin, ymax)
  x@data$id <- rownames(x@data)
  
  x.f = fortify(x, region="id")
  x.join = join(x.f, x@data, by="id")
  
  x.subset <- subset(x.join, x.join$long > longlat[1] & x.join$long < longlat[2] &
                       x.join$lat > longlat[3] & x.join$lat < longlat[4])
  
  x.subset
}

domain <- c(-158.5, -157.5, 21, 22)
lakes.subset <- quick.subset(ne_lakes, domain)
river.subset <- quick.subset(ne_rivers, domain)
coast.subset <- quick.subset(ne_coast, domain)

nat.crop <- crop(nat.earth, y=extent(domain))

rast.table <- data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)),
                         getValues(nat.crop/255))

rast.table$rgb <- with(rast.table, rgb(NE2_HR_LC_SR_W.1,
                                       NE2_HR_LC_SR_W.2,
                                       NE2_HR_LC_SR_W.3,
                                       1))
# et voila!

ggplot(data = rast.table, aes(x = x, y = y)) +
  geom_tile(fill = rast.table$rgb) +
  geom_polygon(data=lakes.subset, aes(x = long, y = lat, group = group), fill = '#ADD8E6') +
  scale_alpha_discrete(range=c(1,0)) +
  geom_path(data=lakes.subset, aes(x = long, y = lat, group = group), color = 'blue') +
  geom_path(data=river.subset, aes(x = long, y = lat, group = group), color = 'blue') +
  geom_path(data=coast.subset, aes(x = long, y = lat, group = group), color = 'blue') +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_map(projection = "mercator") + 
  xlab('') + ylab('')

