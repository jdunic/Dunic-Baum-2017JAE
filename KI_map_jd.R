#  KI_maps.R
#  Maps of sites at Kiritimati
#  
#  Created by Rowan Trebilco on 2011-06-14.
#  Copyright 2011 __MyCompanyName__. All rights reserved.
# 

# note: I originally tried reading in the whole line islands shapefile, but 
#couldn't work out how to subset to just Kiritimati

# ultimately I used QGIS to do the subset using: Vector/Data Management tools
#/Split Vector layer

#source("//Users//RowanT//Dropbox//PhD/RT_Kiritimati//r\ code//kiritimati_analysis_functions.R")

library(maptools)
library(RColorBrewer)
library(ggplot2)

library(maps)
library(mapdata)
library(mapproj)

# Location of Kiritimati in the Pacific
data(worldHiresMapEnv)

dev.new(width=2.5, height=2.5)

par(plt = c(0.75, 1, 0.75, 1), new = TRUE)
grid.rect()

print(
plot.map("worldHires", center=180, col="grey70", bg="white",
   fill=TRUE, ylim=c(-65, 65), mar=c(0,0,0,0), xlimits = c(-40, 65), 
   resolution = 0, lwd = 0.1))
points(x = 23, y = 1.88, pch = 0, cex = 3)

png("/Users/jillian/R_projects/Allometry/Pacific_map.png", units = 'in',
  width = 2.75, height = 2.75, res = 0)
plot.map("worldHires", center=180, col="grey70", bg="white",
   fill=TRUE, ylim=c(-65, 65), mar=c(0,0,0,0), xlimits = c(-40, 65), 
   resolution = 0, lwd = 0.1)
points(x = 23, y = 1.88, pch = 0, cex = 4)

dev.copy2png(device=quartz, file = "/Users/jillian/R_projects/Allometry/Pacific_map.pdf")
dev.copy(device = quartz, "/Users/jillian/R_projects/Allometry/Pacific_map.png")
quartz.save("/Users/jillian/R_projects/Allometry/Pacific_map.png", type = "png", 
  device = dev.cur(), dpi = 200)


# Kiritimati map
kiM <- readShapePoly("KI_maps/Line_v3_Island__Kiritimati", 
        proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
       )
sites <- read.csv(file = "KI_maps/Sites.csv", header = TRUE, 
                  stringsAsFactors = FALSE) %>% dplyr::as_data_frame()

jd_sites <- subset(sites, sites$SiteName %in% unique(pento$Site))

# Must tell R what the coordinate values are. This does not need to be set more
# than once.
coordinates(jd_sites)<-c("Longitude","Latitude")
#proj4string(sites09)<- CRS("+proj=longlat+ellps=WGS84")

dev.new(width=6, height=6)
grid.newpage()
KI_map <- plot(kiM,
               col = ifelse(kiM$RB_DEPTH_A == "land", "light grey", "dark grey"),
               axes = FALSE,
               xlim = c(-157.6, -157.1), ylim = c(1.6, 2.15)
               )
par(mar = c(5.3, 5.1, 4.1, 5.1))
#box(bty = 'L')
points(jd_sites, pch = 21, col = "black", bg = "white", cex = 1.2)
points(jd_sites, pch = 21, col = "black", bg = "black", cex = 0.3)
#points(jd_sites, pch = 21, cex = 1.25, bg = as.character(jd_sites$colour))
axis(side = 2, cex.axis = 0.8, pos = -157.6) #mgp = c(2, 0.5, 0))
axis(side = 1, cex.axis = 0.8, pos = 1.6)
title(xlab = "Longitude", line = 2, ylab = "Latitude")
#mtext("Latitude", side = 2, line = 1)
# legend(x = -157.583, y = 1.76, 
#       legend = fp_colours$reg,
#       pch = 21,
#       pt.bg = as.character(fp_colours$colour),
#       bty = 'n',
#       cex = 0.7,
#       pt.cex = 1.2,
#       y.intersp = 1.1
#       )
map.scale(x = -157.45, y = 1.68, relwidth = 0.2, cex = 0.7, ratio = FALSE)
x <- c(-157.56, -157.573, -157.56, -157.547, -157.56) + 0.04
y <- c(2.08, 2.04, 2.05, 2.04, 2.08) - 0.38
polygon(x = x, y = y, border = "black", col = "black")
text(x = x[3], y = 1.72, "N")


dev.copy2pdf(device=quartz, file = "/Users/jillian/R_projects/Allometry/KI_map.pdf")
