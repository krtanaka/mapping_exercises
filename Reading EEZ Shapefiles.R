# Palacios, D.M. 2016. Marine Boundaries in R: Reading EEZ Shapefiles. 
# RPubs. Accessed dd/mm/yyyy. https://rpubs.com/danielequs/marine_boundaries

rm(list = ls())

library("rgdal") # for `ogrInfo()` and `readOGR()`
library("tools") # for `file_path_sans_ext()`
library("dplyr") # for `inner_join()`, `filter()`, `summarise()`, and the pipe operator (%>%)
library("ggplot2") # for `fortify()` and for plotting
library("sp") # for `point.in.polygon()` and `spDists()`
library("tidyr") # for `gather()`
library("readr") # for `write_tsv()`
library("mapproj") # for ggplot:coord_map()

# Provide the function fortify.shape(), which puts the shapefile data in the object class data.frame, so that it can be used by ggplot2. (These steps follow Hadley Wickham's post "Plotting polygon shapefiles").

fortify.shape <- function(x){
  x@data$id <- rownames(x@data)
  x.f <- fortify(x, region = "id")
  x.join <- inner_join(x.f, x@data, by = "id")
}

# Provide the function subset.shape(), which we will use to extract portions of the data (from the fortified data.frame object) for a smaller domain, since shapefiles often contain global data.

subset.shape <- function(x, domain){
  x.subset <- filter(x, 
                     long > domain[1] & 
                       long < domain[2] & 
                       lat > domain[3] & 
                       lat < domain[4])
  x.subset
}

###########################################################
### Plotting the coastline and some animal observations ###
###########################################################
path.ne.coast <- ("G:/GIS/nm/ne_10m_coastline")
fnam.ne.coast <- "ne_10m_coastline.shp"
dat.coast <- readOGR(dsn = path.ne.coast, layer = file_path_sans_ext(fnam.ne.coast))

# Fortify the shapefile data using `fortify.shape()`:
dat.coast <- fortify.shape(dat.coast) # a 410951x8 dataframe

# Specify the desired domain (the West Coast of the USA):
domain <- c(-128, -115, 29, 50)

# Extract the coastline data for the desired domain using `subset.shape()`:
dat.coast.wc <- subset.shape(dat.coast, domain) # a 4871x8 dataframe

xlims <- c(-135, -116)
ylims <- c(29.5, 49.5)

# Generate a base map with the coastline:
(p0 <- ggplot() + 
    geom_path(data = dat.coast.wc, 
              aes(x = long, y = lat, group = group), 
              color = "black", size = 0.25) + 
    coord_map(projection = "mercator") + 
    scale_x_continuous(limits = xlims, expand = c(0, 0)) + 
    scale_y_continuous(limits = ylims, expand = c(0, 0)) + 
    labs(list(title = "", x = "Longitude", y = "Latitude")))

# Simulate 500 whale observations:
set.seed(250)
coord.ctr <- c(-124.25, 40)
long <- runif(n = 500, min = coord.ctr[1]-10, max = coord.ctr[1])
lat <- runif(n = 500, min = coord.ctr[2]-10, max = coord.ctr[2]+10)

# Compute all pair-wise distances between observations:
dist.obs <- spDists(cbind(long, lat), longlat = TRUE) # a matrix, in km
diag(dist.obs) <- NA # convert 0s to NAs along the diagonal

# For each observation obtain the distance to the nearest observation:
obs.dist.near <- as.data.frame(dist.obs) %>% 
  gather(key = obs_pair, value = pair_dist) %>% 
  group_by(obs_pair) %>% 
  summarise(near_dist = min(pair_dist, na.rm = TRUE))

obs.dist.near <- transform(obs.dist.near, obs_pair = as.factor(obs_pair))
range(obs.dist.near$near_dist) # 1.668178 132.716169 km

# Simulate group size for each observation:
set.seed(123)
n.whales <- rpois(500, lambda = sqrt(obs.dist.near$near_dist)/4)+1
range(n.whales) # 1-8

# Combine all vectors into a data frame:
sim.obs <- data.frame(long = long, lat = lat, n.whales = n.whales) # 500x3

ggplot(data = sim.obs, aes(x = n.whales, ..count..)) + 
  geom_density(adjust = 2, colour = "blue", fill = "blue", alpha = 0.25)

(p1 <- p0 + 
    geom_point(data = sim.obs, 
               aes(x = long, y = lat, 
                   colour = as.factor(n.whales)), 
               shape = 16, size = 1) + 
    scale_colour_brewer(palette = "YlOrBr", 
                        guide_legend(title = "Number \nof whales")))

########################################
### Reading NOAA's USA EEZ shapefile ###
########################################

path.eez.usa <- ("G:/GIS/eez/USMaritimeLimitsAndBoundariesSHP")
fnam.eez.usa <- "USMaritimeLimitsNBoundaries.shp"
eez.usa <- readOGR(dsn = path.eez.usa, layer = file_path_sans_ext(fnam.eez.usa))
# eez.usa has 259 features and 16 fields
# A Large SpatialLinesDataFrame object with 259 features and 16 fields (3.3 Mb)

# Fortify the shapefile data using `fortify.shape()`:
dat.eez.usa1 <- fortify.shape(eez.usa) # a 180400x22 dataframe

(p.eez.noaa <- p0 +  
    geom_path(data = filter(dat.eez.usa1, REGION == "Pacific Coast"), 
              aes(x = long, y = lat, group = group), 
              colour = "red", size = 0.75))

(p.eez.noaa2 <- p0 + 
    geom_path(data = filter(dat.eez.usa1, REGION == "Pacific Coast" & TS == 1), 
              aes(x = long, y = lat, group = group), 
              colour = "tomato2", size = 0.75) + 
    geom_path(data = filter(dat.eez.usa1, REGION == "Pacific Coast" & CZ == 1), 
              aes(x = long, y = lat, group = group), 
              colour = "purple2", size = 0.75) + 
    geom_path(data = filter(dat.eez.usa1, REGION == "Pacific Coast" & EEZ == 1), 
              aes(x = long, y = lat, group = group), 
              colour = "royalblue2", size = 0.75))

###############################################################
### Reading the global EEZ shapefile from Marineregions.org ###
###############################################################

#path.eez.world.v8 <- ("/Users/danielpalacios/Documents/DMP/dmp_data/World_EEZ/World_EEZ_v8_20140228_LR")
path.eez.world.v9 <- ("G:/GIS/eez/World_EEZ_v9_20161021")
#fnam.eez.world.v8 <- "World_EEZ_v8_2014.shp"
fnam.eez.world.v9 <- "eez"
#eez.world.v8 <- readOGR(dsn = path.eez.world.v8, 
#                        layer = file_path_sans_ext(fnam.eez.world.v8))
eez.world.v9 <- readOGR(dsn = path.eez.world.v9, 
                        layer = file_path_sans_ext(fnam.eez.world.v9))
# A Large SpatialLinesDataFrame object with 281 features and 23 fields (18.9 Mb)

# Extract the EEZ for the USA:
#dat.eez.usa2 <- eez.world.v8[eez.world.v8@data$Country == "United States", ]
# For v. 9 use $Territory1 instead of $Country:
dat.eez.usa2 <- eez.world.v9[eez.world.v9@data$Territory1 == "United States", ]
# A Formal class Large SpatialPolygonsDataFrame

# Fortify the shapefile data:
#dat.eez.usa2 <- fortify(dat.eez.usa2)
# `fortify.shape()` did not work for v. 8 so we had to use `fortify()`.
# message: Regions defined for each Polygons
dat.eez.usa2 <- fortify.shape(dat.eez.usa2) # a 10298x30 dataframe

(p.eez.vliz <- p0 + 
    geom_path(data = filter(dat.eez.usa2, piece == 2 & PolygonID == 221), 
              aes(x = long, y = lat, group = group), 
              colour = "blue", size = 0.75))

##############################################################
### Extracting simulated whale observations inside the EEZ ###
##############################################################

# Let's extract the Pacific EEZ polygon into a dataframe for further use:
dat.eez.pac <- droplevels(filter(dat.eez.usa2, piece == 2 & PolygonID == 221)) # 2031x30

dat.eez.exp <- dat.eez.pac %>% select(long, lat) # 2031x2
path.eez.exp <- ("/Users/kisei.tanaka/Desktop/")
fnam.eez.exp <- "eez_usa_pac_lr_v9.txt"
write_tsv(dat.eez.exp, path = paste(path.eez.exp, fnam.eez.exp, sep = "/"))

inside.eez <- point.in.polygon(sim.obs$long, 
                               sim.obs$lat, 
                               dat.eez.pac$long, 
                               dat.eez.pac$lat) # 500x1

# Add a column to sim.obs with this information:
sim.obs$eez <- inside.eez # 500x4

# Extract the observations in sim.obs that occur inside the EEZ:
sim.obs.eez <- sim.obs %>% 
  filter(eez == 1)  %>% 
  select(-eez) # 126x3

(p2 <- p0 + 
    geom_point(data = sim.obs, aes(x = long, y = lat), 
               colour = "gray75", shape = 16, size = 1) + 
    geom_point(data = sim.obs.eez, 
               aes(x = long, y = lat, colour = as.factor(n.whales)), 
               shape = 16, size = 1) + 
    scale_colour_brewer(palette = "YlOrBr", 
                        guide_legend(title = "Number \nof whales")) + 
    geom_path(data = dat.eez.pac, aes(x = long, y = lat, group = group), 
              colour = "blue", size = 0.75))

(p3 <- p0 + 
    stat_density_2d(data = sim.obs.eez, geom = "polygon", aes(x = long, y = lat, fill = ..level..)) + 
    scale_fill_viridis_c(guide_legend(title = "Probability \ndensity")) + 
    geom_point(data = sim.obs, aes(x = long, y = lat), colour = "gray75", shape = 16, size = 1) + 
    geom_point(data = sim.obs.eez, aes(x = long, y = lat), colour = "black", shape = 16, size = 1) + 
    geom_path(data = dat.eez.pac, aes(x = long, y = lat, group = group), colour = "blue", size = 0.75) + 
    labs(list(title = "Simulated observations in \nthe Pacific EEZ of the USA", x = "Longitude", y = "Latitude")))

