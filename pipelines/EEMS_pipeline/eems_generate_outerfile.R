args <- commandArgs(TRUE)

if (any(is.na(args[1:3]))) { stop("Provide arguments.") }

sdfile <- args[1]
regionsfile <- args[2]
outputfile <- args[3]
plot <- ifelse(is.na(args[4]), TRUE, args[4])

# Load library
library('sf')
library(raster)
library(rgdal)
library(plyr)
library(dplyr)
library(sp)
library(ecospat)

# Load shapefile
cat("\nLoading shape files...")
sdm <- read_sf(sdfile)
regions <- read_sf(regionsfile)
regions <- st_transform(regions, "+init=epsg:4326") # transform to WGS84

cat("\nProcessing...")
# Convert sf object to SpatialPolygons
sp_polygons <- st_as_sf(regions)
sp_polygons <- as(sp_polygons, "Spatial")

# Extract the attributes from the sf object
sp_data <- data.frame(regions)

# Create a SpatialPolygonsDataFrame
spdf <- SpatialPolygonsDataFrame(sp_polygons, data = sp_data)

# extract coordinates from sps distribution
coords <- as.data.frame(st_coordinates(sdm)[,1:2])
coords[coords[,1] > 180,1] <- coords[coords[,1] > 180,1] - 360
coordinates(coords) <- ~X+Y
crs(coords) <- "+init=epsg:4326"

# overlay to get polygons of interest
sp_idx <- names(na.exclude(over(spdf, coords)))
spdf <- spdf[sp_idx,]
spdf <- st_as_sf(spdf)
spdf <- st_shift_longitude(spdf)

res <- st_bbox(extent(spdf))
buffer_distance <- oce::geodDist(as.numeric(res[1]), as.numeric(res[2]),
                      as.numeric(res[3]), as.numeric(res[4])) * 10

cat("\nBuffering...")
# Create a buffer around the points
buffered_polygons<- st_buffer(spdf, dist = buffer_distance)
# Join the buff ered polygons into a single polygon
unioned_polygon <- st_union(buffered_polygons)
unioned_polygon <- st_shift_longitude(unioned_polygon)
# Create a convex hull around the unioned polygon to connect the buffers
convex_hull <- st_convex_hull(unioned_polygon)

coords <- st_coordinates(convex_hull)[,1:2]
cat("Writting file...")
# write .outer
write.table(coords, quote = F, col.names = F, row.names = F , file = paste0(outputfile, ".outer"))

# plot
if (plot) {
  cat("\nPloting area...")
  tiff(filename = paste0(outputfile, "_area.tiff"), width = 1200, height = 800, units = "px")
  par(mfrow = c(1,1))
  maps::map("world", wrap = c(0,360), 
            xlim = c(min(res[1], 20), max(res[3], 250)), 
            ylim = c(min(res[2], -40), max(res[4], 40)), 
            col = "grey95", border = "grey85", fill=TRUE, bg = "#79B0D3")
  plot(sdm, border = "orange", lwd = 2, add = T, col = "orange")
  #plot(unioned_polygon, border = "red", lwd = 2, add = T, col = scales::alpha("pink", 0.5))
  plot(convex_hull, add = TRUE, border = "black", lwd = 2, lty = 2)
  dev.off()
}



