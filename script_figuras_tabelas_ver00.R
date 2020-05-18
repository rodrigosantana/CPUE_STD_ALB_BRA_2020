########################################################################
## Description: Preparing all figures to article...
##
## Maintainer: ICCAT / Brazilian SCC for Tunas and Tunas like
## Author: Rodrigo Sant'Ana
## Created: dom mai 17 23:58:45 2020 (-0300)
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
######@> Loading R packages for the exploratory analyses...

######@> Packages list
## library(plyr)
library(dplyr)
## library(dtplyr)
library(ggplot2)
## library(lme4)
## library(doBy)
## library(gridExtra)
## library(MuMIn)
## library(lsmeans)
## library(sjstats)
## ##library(lmerTest)
## ##library(pbkrtest)
## library(corrplot)
library(cluster)
library(reshape2)
## library(grImport)
library(date)
library(splines)
library(maps)
library(mapdata)
library(maptools)
## library(data.table)
## library(lunar)
library(lubridate)
library(readr)
## library(tm)
## library(readxl)
library(rpart)
library(randomForest)
## library(mgcv)
library(influ)
library(nFactors)
library(boot)
library(beanplot)
library(tmap)
library(tmaptools)
library(sf)

########################################################################
######@> Installing cpue.rfmo package...

#####@> Installing package builded in my pc...
## install.packages("/home/rodrigo/Github/cpue.rfmo_0.1.0.tar.gz",
##                  repo = NULL)

######@> Loading cpue.rfmo package...
library(cpue.rfmo)

######@> Package Description...
packageDescription("cpue.rfmo")

########################################################################
######@> Setup R...

######@> Customization for ggplot2 theme...
seta <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open")
my_theme <- function(base_size = 12, base_family = "Arial") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(axis.ticks = element_blank(),
              axis.line = element_line(arrow = seta),
              legend.background = element_blank(),
              legend.key = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              ## strip.background = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(hjust = 1),
              complete = TRUE)
}

########################################################################
######@> Loading spatial datasets...

######@> Base map...
mm <- map_data("world")

#####@> Loading a base map from R...
data(World)

####@> converting CRS...
World <- st_transform(World, crs = 4326)

########################################################################
######@> Preparing maps...

######@> Mapping areas using in the models...
load("ICCAT/2020_Albacore/BR/analyses/BRdat.RData")

#####@> Choosing default colors...
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
col <- gg_color_hue(5)

#####@> Filtering dataset - areas 2 and 3 only...
tmp <- prepdat3 %>%
    filter(regY1 %in% 2:3) %>%
    group_by(lon5, lat5, regY1) %>%
    summarise(n = n()) %>%
    as.data.frame()
tmp$regY1 <- factor(tmp$regY1)

#####@> Converting tmp into spatial...
tmp.sp <- st_as_sf(tmp, coords = c("lon5", "lat5"))

#####@> Setting the crs...
st_crs(tmp.sp) <- st_crs(World)

p00 <- tm_shape(World, xlim = c(-60, 20), ylim = c(-50, 20)) +
    tm_polygons() +
    tm_grid(x = seq(-60, 20, by = 10), y = seq(-50, 20, by = 10)) +
    tm_shape(tmp.sp) +
    tm_dots(col = "regY1", palette = col[4:5], size = 3.2,
            title = "Areas", shape = 15, alpha = 0.8) +
    tm_xlab(text = "Longitude", size = 1) +
    tm_ylab(text = "Latitude", size = 1) +
    tm_layout(scale = 2, legend.position = c("left", "top"),
              bg.color = "white", legend.bg.color = "white",
              legend.frame = "black") +
    tm_scale_bar(width = 0.15, position = c("left", "bottom")) +
    tm_compass(position = c("right", "top"), type = "8star",
               color.dark = "gray30") +
    tm_credits("Fonte: Brazilian National Tuna Database")
p00

png("ICCAT/2020_Albacore/BR/figures/Map_areas.png", res = 200,
    units = "cm", w = 25, h = 25)
print(p00)
dev.off()

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
