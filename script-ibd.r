#' This R script
#'
#' @author  Laura Benestan
#'
#' @date 2023/07/14
#'
#' --------------------------------------------------------------------------   @libraries

### Download libraries
library(marmap)
library(dplyr)
library(ade4)
library(adespatial)
library(ggmap)
library(ggplot2)
library(reshape2)
library(mapdata)

############## IMPORT SPATIAL COORDINATES #############

### Download geographic coordinates
sites.all <- read.table("geo_urchin_modified.txt", header=TRUE,dec=".",sep="\t")

### Order the site object
sites.subset <- sites.all %>% dplyr::arrange(CODE)

### Check the sites positions
summary(sites.subset)

### Keep only latitude and longitude info
sites <- dplyr::select(sites.subset,Longitude,Latitude)
sites

############## GET BATHYMETRIC DATA #############

### Get bathymetric dataset
bathydata <- marmap::getNOAA.bathy(lon1= 94,
                                   lon2= 129,
                                   lat1= -9,
                                   lat2= 6,
                                   resolution = 1)

### Summarizing the data
summary(bathydata)

############## EXTRACT BATHYMETRIC DATA #############

### Plot map according to different levels of depth
# set colors for each level
blues <- colorRampPalette(c("lightblue", "cadetblue2", "cadetblue1", "white"))
blues <- c("lightsteelblue4", "lightsteelblue3","lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

### Plot a map with color
plot(bathydata, image = TRUE, land = TRUE, n=1, 
     bpal = list(c(0, max(bathydata), greys), 
                 c(min(bathydata), 0, blues)))

# add sampling points, and add text to the plot:
points(sites$Longitude, sites$Latitude, pch = 21, col = "black", 
       bg = "yellow", cex = 1.3)
text(sites$Longitude, sites$Latitude,sites$CODE, pos = 2)

### Plot a map without color
pdf("Marmap_seaurchins.pdf")
plot(bathydata, lwd = c(0.3, 1), lty = c(1, 1),
     deep = c(-4500, 0), shallow = c(-50, 0), 
     step = c(500, 0),
     col = c("grey", "black"), drawlabels = c(FALSE, FALSE))
scaleBathy(bathydata, deg = 3, x = "bottomleft", inset = 5)
points(sites$Longitude, sites$Latitude, pch = 21, col = "black", bg = "red", cex = 1)
text(sites$Longitude, sites$Latitude,sites$Code, pos = 1,cex = 0.5)
dev.off()

############## GET BATHYMETRIC INFO #############

### Get the bathymetric value of each sampling point using the interactive map
indv_depth <- marmap::get.depth(bathydata, sites, locator=FALSE)
indv_depth

############## CALCULATE GEOGRAPHIC DISTANCES #############
#### IN-WATER DISTANCES
### Calculate the shortest by-water path
#first, define the constraint for calculating a least-cost-path. make transition object.
t <- marmap::trans.mat(bathydata, min.depth = -0.1) #travel cant be shallower than 0.1 metres depth. It can take a while to calculate. This is a "transition object"

### Get km distance matrix
leastDist.km <- marmap::lc.dist(t, sites, res = "dist") #use "dist" instead of path to get a kilometres distance matrix between all sites.
leastDist.km

### Get the minimum distance value
min_dist <-min(leastDist.km)
min_dist

### Get the minimum distance value
max_dist <-max(leastDist.km)
max_dist

### Check the path
leastDist.path <- marmap::lc.dist(t, sites, res = "path") #use "path" to get visualize the path 
lapply(leastDist.path, lines, col = "red", lwd = 2, lty = 1) 
leastDist.path

#### EUCLIDEAN DISTANCES
### Create an euclidean matrix
euclidean_distances <- dist(sites$Latitude,sites$Longitude, method="euclidean")
geo_distances <- as.dist(euclidean_distances)

############## MANTEL TEST WITH ALL SAMPLES #############

### Download fst values
load("fst_urchins.Rdata")

### Download FST data
fst_distances <- as.dist(urchin_fst)

### Perform a Mantel test on euclidean and in-water distances
ade4::mantel.rtest(fst_distances, euclidean_distances, nrepet = 9999)
ade4::mantel.rtest(fst_distances, leastDist.km, nrepet = 9999)

############## MANTEL TEST WITHOUT SAB SAMPLES #############

### Download FST data
fst_distances.data <- as.matrix(urchin_fst)
fst.subset <- fst_distances.data[,colnames(fst_distances.data)!="SAB"]
fst.subset.no.sab<- fst_distances.data[,rownames(fst_distances.data)!="SAB"]

### Create an euclidean matrix
leastDist.km.no.sab <- as.matrix(leastDist.km)[,-10]

### Perform a Mantel test on euclidean and in-water distances
ade4::mantel.rtest(as.dist(fst.subset.no.sab), as.dist(leastDist.km.no.sab), nrepet = 9999)

############## MANTEL TEST WITHOUT SAB, PAR AND BAL SAMPLES #############

### Download FST data
drop.cols <- c("SAB","BAT","PAR")
fst.subset.columns <- fst_distances.data[, setdiff(colnames(fst_distances.data), drop.cols)]
fst.subset.no.sab.par.bat <-fst.subset.columns[-c(4,9,10),]

### Create an euclidean matrix
leastDist.km.no.sab.par.bat <- as.matrix(leastDist.km)[-c(4,9,10),]

### Perform a Mantel test on euclidean and in-water distances
ade4::mantel.rtest(as.dist(fst.subset.no.sab.par.bat), as.dist(leastDist.km.no.sab.par.bat), nrepet = 9999)

############## NICE MANTEL TEST VISUALIZATION #############

fst_distances.data.melted <- melt(fst_distances.data)
in.water.distances.melted <- melt(leastDist.km)
colnames(leastDist.km) <- colnames(fst_distances.data)

cbind(fst_distances.data.melted, in.water.distances.melted)
