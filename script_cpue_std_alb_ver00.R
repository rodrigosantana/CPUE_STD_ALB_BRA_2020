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
library(readxl)
library(dplyr)
library(tidyr)
library(janitor)
library(ggplot2)
library(patchwork)
library(extrafont)
library(tmap)
library(tmaptools)
library(rgeos)
library(rgdal)
library(sf)
library(viridis)
library(INLA)

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

######@> Cerco SC...
db01 <- read.table("data/BNDA_1978-2018_07-08-2019.csv",
                   header = TRUE, sep = ";", dec = ".",
                   fileEncoding = "ISO-8859-1")

########################################################################
######@> Cleaning and tidying data...

######@> Selecting only variables that will be used in analysis...
db02 <- select(db01, hpb, effort, boat, setmonth, setyear, set.latb,
               set.lonb, ALB.n, ALB.kg)

######@> Cleaning latitude data...
unique(db02$set.latb)
db02$lat <- as.numeric(gsub(",", ".", db02$set.latb))

######@> Cleaning longitude data...
unique(db02$set.lonb)
db02$lon <- as.numeric(gsub(",", ".", db02$set.lonb))

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
