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

grid.newpage()


a <- viewport(width = 5, height = 5, default.units = "inches")  # sets a viewport
pushViewport(a)  # opens the active viewport (not necessary in RStudio)
grid.rect()

b <- viewport(x = 0.8, y = 0.8, height = 0.35, width = 0.35)
    # height and width 0.7 make viewport b 70% of the size of a (the active
    # viewport)
pushViewport(b)
grid.rect()  # draws a rectangle in viewport b
grid.xaxis() # draws x-axes in viewport b
grid.yaxis() # draws y-axes in viewport b

x <- runif(100, 0, 1) 
y <- runif(100, 0, 1)

grid.points(x,y)

c <- viewport(x = 0.25, y = 0.75, height = 0.2, width = 0.2) 
    # will make a viewport 'c' in the topleft corner of viewport b
    # this viewport will be 20% the size of viewport b
pushViewport(c) # don't see anything after doing this
grid.rect()  # grid.rect() draws a rectangle the size of the viewport by default
grid.xaxis() # makes x-axis in viewport c 
grid.yaxis()
grid.points(x,y)

print(xyplot(decrease ~ treatment, OrchardSprays, groups = rowpos, type = "a"), newpage=FALSE)

data(worldHiresMapEnv)

dev.new(width=2.5, height=2.5)

par(plt = c(0.75, 1, 0.75, 1), new = TRUE)
grid.rect()

print(
plot.map("worldHires", center=180, col="grey70", bg="white",
   fill=TRUE, ylim=c(-65, 65), mar=c(0,0,0,0), xlimits = c(-40, 65), 
   resolution = 0, lwd = 0.1)
points(x = 23, y = 1.88, pch = 0, cex = 3)

png("/Users/jillian/R_projects/Allometry/Pacific_map.png", units = 'in',
  width = 2.75, height = 2.75, res = 0)
plot.map("worldHires", center=180, col="grey70", bg="white",
   fill=TRUE, ylim=c(-65, 65), mar=c(0,0,0,0), xlimits = c(-40, 65), 
   resolution = 0, lwd = 0.1)
points(x = 23, y = 1.88, pch = 0, cex = 4)
dev.off()


dev.copy2png(device=quartz, file = "/Users/jillian/R_projects/Allometry/Pacific_map.pdf")
dev.copy(device = quartz, "/Users/jillian/R_projects/Allometry/Pacific_map.png")
quartz.save("/Users/jillian/R_projects/Allometry/Pacific_map.png", type = "png", 
  device = dev.cur(), dpi = 200)

# See note about this function in /Allometry/03_func.r
plot.map <- function(database,center, xlimits, ...){
    Obj <- map(database,...,plot=F)
    coord <- cbind(Obj[[1]],Obj[[2]])
    # split up the coordinates
    id <- rle(!is.na(coord[,1]))
    id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
    polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})
    # split up polygons that differ too much
    polygons <- lapply(polygons,function(x){
        x[,1] <- x[,1] + center
        x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
        if(sum(diff(x[,1])>300,na.rm=T) >0){
          id <- x[,1] < 0
          x <- rbind(x[id,],c(NA,NA),x[!id,])
       }
       x
    })
    # reconstruct the object
    polygons <- do.call(rbind,polygons)
    Obj[[1]] <- polygons[,1]
    Obj[[2]] <- polygons[,2]
    map(Obj,..., xlim=xlimits)
}

kiM <- readShapePoly("KI_maps/Line_v3_Island__Kiritimati", 
        proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
       )
sites <- read.csv(file = "KI_maps/Sites.csv", header = TRUE, 
                  stringsAsFactors = FALSE)


# renaming CombinedFishingProductivity
colnames(sites)[2] <- "fp"

jd_sites <- subset(sites, sites$SiteName %in% unique(pento$Site))

# set colour codes for plotting productivity and fishing
# Want to use Set 1 from the Brewer colour pallet
regions <- c("LP.MF", "LP.LF", "HP.MF", "HP.LF", "HP.HF")

# Colour blind friendly palette:
# Black,   Orange,   Light blue,   Green,    Yellow,   Dark Blue,  Red,      Pink 
"#000000", "#E6"9"F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"

cb <- c("#D55E00", "#56B4E9", "#009E73", "#0072B2", "#CC79A7")
cb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
colour <- c("green", "blue", "dark green", "orange", "red")

fp_colours <- data.frame(colour = cb, reg = regions)

jd_sites <- merge(jd_sites, fp_colours, by.x = "fp", by.y = "reg", all.x = TRUE)


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

line<- readShapePoly("Line_v3", proj4string=CRS("+proj=longlat 
                                                +ellps=WGS84 
                                                +datum=WGS84 
                                                +no_defs"))

line_map <- plot(line, 
               col = ifelse(line$RB_DEPTH_A == "land", "light grey", "dark grey")
               axes = FALSE)



ggmap(
    get_googlemap(center = "Kiritimati",
                    zoom = 4,
                  colour = "bw"))

ggmap(
    get_googlemap("Kiritimati",
                  zoom = 100,
                  maptype = "terrain",
                  colour = "bw"
)
        "Kiritimati"))


get_googlemap(center = c(lon = -95.3632715, lat = 29.7632836),
              zoom = 10, size = c(640, 640), scale = 2,
              format = c("png8", "gif", "jpg", "jpg-baseline", "png32"),
              maptype = c("terrain", "satellite", "roadmap", "hybrid"),
              language = "en-EN", region, markers, path, visible,
              style, sensor = FALSE, messaging = FALSE,
              urlonly = FALSE, filename = "ggmapTemp",
              color = c("color", "bw"), ...)

ggmap(
    get_googlemap(
        center=c(-3.17486, 55.92284), #Long/lat of centre, or "Edinburgh"
        zoom=14, 
        maptype='satellite', #also hybrid/terrain/roadmap
        scale = 2), #resolution scaling, 1 (low) or 2 (high)
    size = c(600, 600), #size of the image to grab
    extent='device', #can also be "normal" etc
    darken = 0) #you can dim the map when plotting on top

ggsave ("/Users/s0679701/Desktop/map.png", dpi = 200) #this saves the output to a file


plotKImap<- function(sitedat=sites,basemap=kiM, categories="smw", to.pdf=F, 
                     pdf.name=paste("kiritimati sites map", categories, "categories.pdf"), 
                     site.lab.exp=0.6, suppress.sites=F){
    if(to.pdf) pdf(file=pdf.name)
    plot(basemap, col=ifelse(kiM$RB_DEPTH_A=="land","light grey","dark grey"), axes=T)
    switch(categories,
           "smw.all"={
               points(sitedat, pch=21, bg=as.character(sitedat$prod.col), cex=2)
               points(sitedat, pch=24, bg=as.character(sitedat$f.col.HL), cex=1)
               legend("topright", 
                      legend=c(levels(sitedat$chlA), levels(sitedat$fishingHL)), 
                      pch=c(21,21,24,24), 
                      pt.bg=c("green", "blue", "red", "dark green"), 
                      bty="n")
           } ,
           "six_levels"={
               points(sitedat, 
                      pch=21, 
                      bg=as.character(sitedat$f.col6), cex=2)
               legend("topright", legend=rev(levels(sitedat$f.pressure)), 
                      fill=as.character(fcols$col), bty="n")	
           } ,
           "smw.prod"={
               points(sitedat, pch=21, bg=as.character(sitedat$prod.col), cex=2)
               legend("topright",legend=c(levels(sitedat$chlA)), 
                      pch=c(21,21), pt.bg=c("green", "blue"), bty="n")
           },
           "smw.fp" ={
               points(sitedat, pch=24, bg=as.character(sitedat$f.col.HL), cex=1)
               legend("topright",legend=c(levels(sitedat$fishingHL)), 
                      pch=c(24,24), pt.bg=c("red", "dark green"), bty="n")
           } ,
           "new.cat" ={
               points(sitedat[sitedat$new.cat!="EXTRA",], pch=21, 
                      bg=as.character(sitedat[sitedat$new.cat!="EXTRA",]$prod.col), cex=2)
               points(sitedat[sitedat$new.cat=="EXTRA",], pch=21, col="gray",cex=2)
               points(sitedat[sitedat$new.cat!="EXTRA",], pch=24, 
                      bg=as.character(sitedat[sitedat$new.cat!="EXTRA",]$f.col3), cex=1)
               legend("topright",
                      legend=c(levels(sitedat$chlA),
                               "Low fishing","Medium fishing","High fishing","extra sites"), 
                      pch=c(21,21,24,24,24,21), 
                      pt.bg=c("green", "blue", "dark green","orange","red","white"), 
                      col=c("black","black","black","black","black","grey"), bty="n")
           } ,
           "none"={
               points(sitedat, pch=19, bg=as.character(sitedat$f.col.HL), 
                      pt.bg="gray30",cex=1)
           }
           
    )	
    cexp<-site.lab.exp
    if(!suppress.sites){
        text(sitedat[sitedat$smw.region=="LON",], 
             labels=sitedat[sitedat$smw.region=="LON",]$site, pos=2, cex=cexp)
        text(sitedat[sitedat$smw.region%in%c("NOR","NWA"),], 
             labels=sitedat[sitedat$smw.region%in%c("NOR","NWA"),]$site, pos=3,cex=cexp)
        text(sitedat[sitedat$smw.region=="BOW",], 
             labels=sitedat[sitedat$smw.region=="BOW",]$site, pos=4,cex=cexp)
        text(sitedat[sitedat$smw.region=="POL",], 
             labels=sitedat[sitedat$smw.region=="POL",]$site, pos=2,cex=cexp)
        text(sitedat[sitedat$smw.region=="VAS",], 
             labels=sitedat[sitedat$smw.region=="VAS",]$site, pos=1,cex=cexp)
        text(sitedat[sitedat$smw.region=="KWR",], 
             labels=sitedat[sitedat$smw.region=="KWR",]$site, pos=2,cex=cexp)
    }
    if(!suppress.sites&categories=="new.cat"){
        ext.sites<-sitedat[sitedat$new.cat=="EXTRA",]
        text(ext.sites[ext.sites$smw.region=="NOR",], labels=ext.sites[ext.sites$smw.region=="NOR",]$site, pos=3,cex=cexp, col="grey")
        #text(ext.sites[ext.sites$smw.region=="POL",], labels=ext.sites[ext.sites$smw.region=="POL",]$site, pos=2,cex=cexp, col="grey")
        text(ext.sites[ext.sites$smw.region=="VAS",], labels=ext.sites[ext.sites$smw.region=="VAS",]$site, pos=1,cex=cexp, col="grey")
        text(ext.sites[ext.sites$smw.region=="KWR",], labels=ext.sites[ext.sites$smw.region=="KWR",]$site, pos=2,cex=cexp,col="grey")
    }
    
    if(to.pdf) dev.off()
}


plotKImap(categories="smw.all",to.pdf=F)
plotKImap(categories="six_levels",to.pdf=F)

## region labels
cexp1<-0.6
text(sites[sites$smw.region=="LON",], labels=sites[sites$smw.region=="LON",]$smw.region, pos=2, cex=cexp1)
text(sites[sites$smw.region%in%c("NOR","NWA"),], labels=sites[sites$smw.region%in%c("NOR","NWA"),]$smw.region, pos=3,offset=1,cex=cexp1)
text(sites[sites$smw.region=="BOW",], labels=sites[sites$smw.region=="BOW",]$smw.region, pos=4,offset=1,cex=cexp1)
text(sites[sites$smw.region=="POL",], labels=sites[sites$smw.region=="POL",]$smw.region, pos=2,offset=1,cex=cexp1)
text(sites[sites$smw.region=="VAS",], labels=sites[sites$smw.region=="VAS",]$smw.region, pos=1,offset=1,cex=cexp1)
text(sites[sites$smw.region=="KWR",], labels=sites[sites$smw.region=="KWR",]$smw.region, pos=2,offset=1,cex=cexp1)

## new category labels
text(sites, labels=sites$new.cat, cex=0.6, pos=2, col="black")

# make a new data frame called area labels and another called area boundaries
# in area labels, put labels for the groups of sites in about the middle of each region
# in area boundaries, put a line centred on the coast at the boundaries between the groups of sites
# get these line start and endpoints from qgis

# break lines at:
# SW corner between 13 and 7
# W face between 34 and 9
# W face between 40 and 32
# NW corner, right through site or just to the S of it
# N coast E of 2
# top of bay of wrecks
# end of the wang

# labels at:
# middle of bay of wrecks: Low productivity/Low fishing 

segments(x=c(-157.4, -157.3), y=c(1.7,1.7))



## for reading from GPS
readGPS(i = "garmin", f = "usb:", type="w", invisible=TRUE, ...)

