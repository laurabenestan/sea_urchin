rm(list=ls())
#' --------------------------------------------------------------------------   @Header
#'
#' @title Figure 1 - Sampling design of sea urchin
#'
#' @description
#' This R script...
#'
#' @author   Laura Benestan, \email{lmbenestan@gmail.com}
#'
#' @date 2023/04/25
#'
#' --------------------------------------------------------------------------   @VVVVVV

#'  -------------------------------------------------------------------------   @library

### Download libraries
install_github("ericpante/marmap")
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
sites.all <- read.table("../../../00-Structure/04-spatial/geo_urchin_modified.txt", header=TRUE,dec=".",sep="\t")

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
                                   lat1= -15,
                                   lat2= 8,
                                   resolution = 1)

### Check in marmap the plot
autoplot.bathy(bathydata, geom=c("r", "c"), colour="white", size=0.1) + 
  scale_fill_etopo()+
  geom_point(data = sites, aes(x = Longitude, y = Latitude),
             fill = 'red', size = 3, alpha = 1, shape = 21) +
  labs(y = "Latitude", x = "Longitude", fill = "Elevation") +
  coord_cartesian(expand = 0)+
  theme_classic()+ 
  geom_label_repel(data=sites.all,aes(x=Longitude, Latitude, label = Code), alpha=1,size=3, fill="white")+
  theme(legend.position = "bottom")+
  theme(legend.key.size = unit(1, 'cm'))
ggsave("Sampling_map.pdf", width=12, height=7)

### Plot a map without color
pdf("Marmap_seaurchins.pdf")
plot(bathydata, lwd = c(0.1, 1), lty = c(1, 1),
     deep = c(-4500, 0), shallow = c(-50, 0), 
     step = c(500, 0),
     col = c("white", "black"))
#scaleBathy(bathydata, deg = 3, x = "bottomleft", inset = 5)
points(sites$Longitude, sites$Latitude, pch = 21, col = "black", bg = "grey", cex = 1)
text(sites$Longitude, sites$Latitude,sites$Code, pos = 1,cex = 0.5,col = "black")
dev.off()
