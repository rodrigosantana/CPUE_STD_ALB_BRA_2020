########################################################################
## Description: Standardizing CPUE of Albacore on brazilian pelagic
## longline fishery...
##
## Maintainer: UNIVALI/EMCT/LEMA / ProTUNA Project
## Author: Rodrigo Sant'Ana
## Created: Sat Apr 25 03:40:53 2020 (-0300)
## Version: 0.0.1
##
## URL:
## Doc URL:
##
## Database info:
##
### Commentary:
##
### Code:
########################################################################

########################################################################
######@> Loading R packages...

######@> Package list...
library(dplyr)
library(tidyr)
library(lubridate)
library(janitor)
library(ggplot2)
library(patchwork)
library(extrafont)
library(maps)
library(maptools)
library(tmap)
library(tmaptools)
library(rgeos)
library(rgdal)
library(sf)
library(viridis)
library(INLA)
library(Matrix)
library(spdep)
library(lunar)
library(suncalc)
library(oce)

########################################################################
######@> Setup R...

######@> Loading and registering new fonts...
font_import()
loadfonts(device = "pdf")

######@> DatenKraft theme...
seta <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open")
my_theme <- function(base_size = 16, base_family = "Trebuchet MS") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(axis.ticks = element_blank(),
              axis.line = element_line(arrow = seta, color = "#373538"),
              legend.background = element_blank(),
              legend.key = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              strip.background = element_rect(fill = "#555353"),
              plot.background = element_blank(),
              complete = TRUE)
}

######@> Standardizing the number of decimal places...
options(scipen = 10)

########################################################################
######@> Loading datasets...

######@> Pelagic longline fishery - BNDA Dataset...
db01 <- read.table("data/BNDA_1978-2018_07-08-2019.csv",
                   header = TRUE, sep = ";", dec = ".",
                   fileEncoding = "ISO-8859-1")

######@> World shapefile...
data(World)

####@> converting CRS...
World <- st_transform(World, crs = 4326)

#####@> World base for simple maps...
mm <- map_data("world")

########################################################################
######@> Cleaning and tidying data...

######@> Selecting only variables that will be used in analysis...
db02 <- select(db01, boat, setmonth, setyear, hpb, effort, set.latb,
               set.lonb, ALB.n, ALB.kg)

######@> Cleaning latitude data...
unique(db02$set.latb)
db02$lat <- trunc(as.numeric(gsub(",", ".", db02$set.latb)))

######@> Cleaning longitude data...
unique(db02$set.lonb)
db02$lon <- trunc(as.numeric(gsub(",", ".", db02$set.lonb)))

#####@> Total effort - Pelagic longline...
total <- db02 %>%
    filter(!is.na(set.latb)) %>%
    group_by(lat, lon) %>%
    summarise(hooks = sum(effort, na.rm = TRUE)) %>%
    as.data.frame()

####@> Transforming db06 into a spatial object...
total.sp <- st_as_sf(total, coords = c("lon", "lat"))

####@> converting CRS...
st_crs(total.sp) <- st_crs(World)

####@> Map of total effort...
p00 <- tm_shape(World, xlim = c(-60, 20), ylim = c(-60, 10)) +
    tm_polygons(col = "#555353", alpha = 0.7) +
    tm_grid(x = seq(-60, 20, by = 10), y = seq(-60, 10, by = 10)) +
    tm_shape(total.sp) +
    tm_dots(col = "hooks", palette = "viridis", size = 0.2,
            title = "Hooks",
            shape = 22, alpha = 0.9,
            style = "quantile", border.col = "#555353") +
    tm_xlab(text = "Longitude", size = 1) +
    tm_ylab(text = "Latitude", size = 1) +
    tm_layout(scale = 2, legend.position = c("right", "bottom"),
              bg.color = "white", legend.bg.color = "white",
              legend.frame = "black") +
    tm_compass(type = "4star", position = c("left", "bottom")) +
    tm_scale_bar(width = 0.15, position = c("left", "bottom")) +
    tm_credits("Source: BNDA - 2020")
p00

####@> exporting figure...
dir.create("Figs")

png("Figs/Map_total_effort_ver00.png", units = "cm", res = 300,
    width = 30, height = 30)
print(p00)
dev.off()

#####@> Creating two new variables that could be used as spatial id for
#####@> each square spatial block...

####@> Text id...
db02$id.sp01 <- factor(paste(db02$lat, db02$lon, sep = " "))

####@> Numeric id...
db02$id.sp02 <- as.numeric(db02$id.sp01)

#####@> Creating a new variable to represent the seasonal effect...
db02$season <- ifelse(db02$setmonth < 4, 1,
               ifelse(db02$setmonth < 7, 2,
               ifelse(db02$setmonth < 10, 3, 4)))

########################################################################
######@> Preparing the spatial structure for the models...

#####@> Building a spatial model matrix...
loc <- db02 %>%
    group_by(id.sp02, lat, lon) %>%
    summarise(n = n()) %>%
    select(id.sp02, lon, lat) %>%
    as.data.frame()

####@> spatial matrix...
loc01 <- as.matrix(loc[, 2:3])
colnames(loc01) <- c("X", "Y")

#####@> Neighborhoods matrix...
neibor01.a <- dnearneigh(loc01, 0, 1, longlat = NULL)
neibor01.b <- dnearneigh(loc01, 0, 2, longlat = NULL)

####@> Transforming neiborhoods objects into spatial objects...
sp.neibor01.a <- st_as_sf(nb2lines(neibor01.a, coords = loc01))
sp.neibor01.b <- st_as_sf(nb2lines(neibor01.b, coords = loc01))

####@> converting CRS...
st_crs(sp.neibor01.a) <- st_crs(World)
st_crs(sp.neibor01.b) <- st_crs(World)

####@> Map Neighborhoods 01.a...
p01 <- tm_shape(World, xlim = c(-60, 20), ylim = c(-60, 10)) +
    tm_polygons(col = "#555353", alpha = 0.7) +
    tm_grid(x = seq(-60, 20, by = 10), y = seq(-60, 10, by = 10)) +
    tm_shape(sp.neibor01.a) +
    tm_lines() +
    tm_xlab(text = "Longitude", size = 1) +
    tm_ylab(text = "Latitude", size = 1) +
    tm_compass(type = "4star", position = c("left", "bottom")) +
    tm_scale_bar(width = 0.15, position = c("left", "bottom"))
p01

####@> Map Neighborhoods 01.b...
p02 <- tm_shape(World, xlim = c(-60, 20), ylim = c(-60, 10)) +
    tm_polygons(col = "#555353", alpha = 0.7) +
    tm_grid(x = seq(-60, 20, by = 10), y = seq(-60, 10, by = 10)) +
    tm_shape(sp.neibor01.b) +
    tm_lines() +
    tm_xlab(text = "Longitude", size = 1) +
    tm_ylab(text = "Latitude", size = 1) +
    tm_compass(type = "4star", position = c("left", "bottom")) +
    tm_scale_bar(width = 0.15, position = c("left", "bottom"))
p02

###@> Exporting maps...
png("Figs/Spatial_Discrete_Distribution_01.png", units = "cm", w = 20, h = 15,
    res = 100)
print(p01)
dev.off()

png("Figs/Spatial_Discrete_Distribution_02.png", units = "cm", w = 20, h = 15,
    res = 100)
print(p02)
dev.off()

####@> converting neiborhoods matrix in INLA format...
nb2INLA("neibor01a.adj", neibor01.a)
nb2INLA("neibor01b.adj", neibor01.b)

######@> Preparing temporal matrix...

#####@> creating an index to time series (years)...
db02$t <- db02$setyear - 1977

########################################################################
######@> Models with INLA - Spatial Discrete Models - iCar...

######@> Clean dataset to run the models...
tmp <- select(db02, boat, setmonth, season, setyear, hpb, effort, lat,
              lon, id.sp02, t, ALB.n)

######@> converting data to factor...
tmp$setyear <- factor(tmp$setyear)
tmp$season <- factor(tmp$season)
tmp$boat <- factor(tmp$boat)

######@> Testing distinct likelihood distributions to response
######@> variable...

#####@> Likelihood distributions tested...
fams <- c("poisson", "nbinomial", "zeroinflatedpoisson0",
          "zeroinflatedpoisson1", "zeroinflatednbinomial0",
          "zeroinflatednbinomial1")

#####@> Model 0 - without spatial structure...
form0 <- ALB.n ~ setyear + boat + season + I(hpb)

#####@> fitting...
info00 <- sapply(fams, function(fam) {
    cat(fam, "\n")
    res <- inla(form0, family = fam,
                data = tmp,
                offset = log(effort),
                control.fixed = list(expand.factor.strategy = "inla"),
                control.predictor = list(compute = TRUE, link = 1),
                ## control.inla = list(strategy = "laplace"),
                control.compute =
                    list(dic = TRUE, waic = TRUE, cpo = TRUE))
    c(dic = res$dic$dic, waic = res$waic$waic, cpu = res$cpu[4])
})

#####@> Model 1 - only with spatial structure...

####@> Scenario A - Besag...
form1.a <- birds ~ -1 + f(ind, model = "besag", graph = g01.a)

###@> fitting...
info01.a <- sapply(fams, function(fam) {
    cat(fam, "\n")
    res <- inla(form1.a, family = fam,
                data = tmp,
                offset = log(hooks),
                control.fixed = list(expand.factor.strategy = "inla"),
                control.predictor = list(compute = TRUE, link = 1),
                ## control.inla = list(strategy = "laplace"),
                control.compute =
                    list(dic = TRUE, waic = TRUE, cpo = TRUE))
    c(dic = res$dic$dic, waic = res$waic$waic, cpu = res$cpu[4])
})

####@> Scenario B - Besag...
form1.b <- birds ~ -1 + f(ind, model = "besag", graph = g01.b)

###@> fitting...
info01.b <- sapply(fams, function(fam) {
    cat(fam, "\n")
    res <- inla(form1.b, family = fam,
                data = tmp,
                offset = log(hooks),
                control.fixed = list(expand.factor.strategy = "inla"),
                control.predictor = list(compute = TRUE, link = 1),
                ## control.inla = list(strategy = "laplace"),
                control.compute =
                    list(dic = TRUE, waic = TRUE, cpo = TRUE))
    c(dic = res$dic$dic, waic = res$waic$waic, cpu = res$cpu[4])
})

####@> Scenario A - Bym...
form1.a2 <- birds ~ -1 + f(ind, model = "bym", graph = g01.a)

###@> fitting...
info01.a2 <- sapply(fams, function(fam) {
    cat(fam, "\n")
    res <- inla(form1.a2, family = fam,
                data = tmp,
                offset = log(hooks),
                control.fixed = list(expand.factor.strategy = "inla"),
                control.predictor = list(compute = TRUE, link = 1),
                ## control.inla = list(strategy = "laplace"),
                control.compute =
                    list(dic = TRUE, waic = TRUE, cpo = TRUE))
    c(dic = res$dic$dic, waic = res$waic$waic, cpu = res$cpu[4])
})

####@> Scenario B - Bym...
form1.b2 <- birds ~ -1 + f(ind, model = "bym", graph = g01.b)

###@> fitting...
info01.b2 <- sapply(fams, function(fam) {
    cat(fam, "\n")
    res <- inla(form1.b2, family = fam,
                data = tmp,
                offset = log(hooks),
                control.fixed = list(expand.factor.strategy = "inla"),
                control.predictor = list(compute = TRUE, link = 1),
                ## control.inla = list(strategy = "laplace"),
                control.compute =
                    list(dic = TRUE, waic = TRUE, cpo = TRUE))
    c(dic = res$dic$dic, waic = res$waic$waic, cpu = res$cpu[4])
})

########################################################################
##
##                  Creative Commons License 4.0
##                       (CC BY-NC-SA 4.0)
##
##  This is a humam-readable summary of (and not a substitute for) the
##  license (https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
##
##  You are free to:
##
##  Share - copy and redistribute the material in any medium or format.
##
##  The licensor cannot revoke these freedoms as long as you follow the
##  license terms.
##
##  Under the following terms:
##
##  Attribution - You must give appropriate credit, provide a link to
##  license, and indicate if changes were made. You may do so in any
##  reasonable manner, but not in any way that suggests the licensor
##  endorses you or your use.
##
##  NonCommercial - You may not use the material for commercial
##  purposes.
##
##  ShareAlike - If you remix, transform, or build upon the material,
##  you must distributive your contributions under the same license
##  as the  original.
##
##  No additional restrictions — You may not apply legal terms or
##  technological measures that legally restrict others from doing
##  anything the license permits.
##
########################################################################
