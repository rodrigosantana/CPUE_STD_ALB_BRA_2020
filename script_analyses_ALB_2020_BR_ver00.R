########################################################################
## Description: Albacore tuna CPUE Standardization - Longline Brazil
##
## Maintainer: ICCAT / Brazilian SCC for Tunas and Tunas like
## Author: Rodrigo Sant'Ana
## Created: seg mai 04 10:01:43 2020 (+0200)
##
## URL:
## Doc URL:
##
## Database info:
## LOA - Boat size (m)
## HBF - Number of hooks per basket
## TYPE - (1) Pelagic LL; (2) Demersal LL
## MOON - Moon illumination (%)
## BAT - Bathymetry (m)
## SST - Sea Surface Temperature (C)
## MLD - Mixture Layer Depth
## TTD - Top Thermocline Depth (m)
## CH - Clorophyll (mg)
## REGB e REGB1 - Areas for bigeye tuna
## TOTAL2 - The sum of YFT, ALB and BET
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
######@> Loading datasets...

######@> Rawdata from Brazil...
db <- read.table("data/BNDA_1978-2018_07-08-2019.csv", header = TRUE,
                 sep = ";", dec = ",")

######@> Base map...
mm <- map_data("world")

########################################################################
######@> Adapting some cpue.rfmo functions...

######@> Dataprep Brazil - descriptive function - it was builded below...
dataprep_BR <- function(dat, splist) {
    dat <- dat %>%
        mutate(dmy = ymd(paste(setyear, setmonth, set.day, sep = "/"))) %>%
        mutate(op_yr = setyear) %>%
        mutate(op_mon = setmonth) %>%
        mutate(op_day = set.day) %>%
        mutate(hbf = hpb) %>%
        mutate(hooks = effort) %>%
        mutate(floats = hooks/hbf)
    dat$lon <- dat$set.lonb
    dat$lat <- dat$set.latb
    dat$qtr <- ceiling(as.numeric(dat$op_mon) / 3)
    dat$yrqtr <- dat$op_yr + floor((dat$op_mon - 1)/3)/4 + 0.125
    dat$latlong <- paste(dat$lat5, dat$lon5, sep = "_")
    dat$vessid <- as.factor(as.numeric(dat$boat))
    dat$tripidmon <- paste(dat$vessid, dat$op_yr, dat$op_mon)
    dat$Total <- apply(dat[, splist], 1, sum, na.rm = TRUE)
    dat$Total2 <- apply(dat[, c("bet", "yft", "alb")], 1, sum,
                        na.rm = TRUE)
    dat <- dat %>%
        select(vessid, flag, tripidmon, dmy, op_yr, op_mon, op_day,
               qtr, yrqtr, hooks, hbf, floats, lat, lon, lat5, lon5,
               latlong, Total, Total2, bft, yft, alb, bet, blf, skj,
               oth, swo, sai, whm, bum, spf, spg, bsh, spx, bth, sma,
               ocs, fal, ccs) %>%
        as.data.frame()
    return(dat)
}

######@> Modifying dataclean function for Brazil - same here...
dataclean_BR <- function(dat, yearlim = 2019, splist) {
    for (sp in splist) {
        dat[, sp] <- as.numeric(dat[, sp])
        if (sum(is.na(dat[, sp])) > 0)
            dat[is.na(dat[, sp]), sp] <- 0
    }
    dat <- dat[!is.na(dat$hooks), ]
    dat <- dat[dat$hooks < 5000, ]
    dat <- dat[dat$hooks >= 500, ]
    dat <- dat[is.na(dat$hbf) == FALSE, ]
    dat <- dat[dat$op_yr > 1976, ]
    dat <- dat[dat$yrqtr < yearlim, ]
    dat <- dat[dat$hbf >= 3, ]
    return(dat)
}

########################################################################
######@> Preparing the structure of folders to receive the outputs...

######@> Directories...
projdir <- "ICCAT/2020_Albacore/"
brdir <- paste0(projdir, "BR/")
datadir <- paste0(brdir, "data/")
bralysis_dir <- paste0(brdir, "analyses/")
brfigs <- paste0(brdir, "figures/")
Rdir <- paste0(projdir, "Rfiles/")
## dir.create(bralysis_dir, recursive = TRUE)
## dir.create(brfigs)
setwd(bralysis_dir)

########################################################################
######@> Cleaning dataset...

######@> Selecting variables...
prepdat <- db %>%
    dplyr::select(boat, flag, hpb, effort, setyear, setmonth, set.day, Arte,
                  set.latb, set.lonb,
                  SBF.n, BFT.n, YFT.n, ALB.n, BET.n, BLF.n, SKJ.n,
                  OTH.n, SWO.n, SAI.n, WHM.n, BUM.n, SPF.n, SPG.n,
                  BSH.n, SPX.n, BTH.n, SMA.n, OCS.n, FAL.n, CCS.n) %>%
    mutate(flag = as.character(flag)) %>%
    as.data.frame()

######@> Correcting flag...
prepdat$flag <- ifelse(is.na(prepdat$flag), "BRA", prepdat$flag)

######@> Correcting hpb data...

#####@> Filtering 2018 boats...
boats.2018 <- data.frame(
    boat = sort(unique(prepdat$boat[prepdat$setyear == 2018])))

####@> Estimating the hpb pattern for 2018 boats...
hpb <- prepdat %>%
    filter(boat %in% boats.2018$boat) %>%
    group_by(boat) %>%
    summarise(hpb.mea = round(mean(hpb, na.rm = TRUE))) %>%
    mutate(hpb.mea = replace(hpb.mea, is.na(hpb.mea),
                             round(mean(hpb.mea, na.rm = TRUE)))) %>%
    as.data.frame()

####@> Imputing the average for each boat in 2018...
for(i in 1:nrow(prepdat)) {
    if(is.na(prepdat$hpb[i])) {
        prepdat$hpb[i] <- hpb$hpb.mea[which(hpb$boat == prepdat$boat[i])]
    } else {
        prepdat$hpb[i] <- prepdat$hpb[i]
    }
}

######@> Correcting the number of floats...
prepdat$floats <- with(prepdat, effort / hpb)

######@> Correcting the lat5 and lon5...
prepdat$lat5 <- 5 * floor(prepdat$set.latb / 5) + 2.5
prepdat$lon5 <- 5 * floor(prepdat$set.lonb / 5) + 2.5
prepdat$latlong <- paste(prepdat$lat5, prepdat$lon5, sep = "_")

######@> Correcting species names...
names(prepdat)[11:31] <- tolower(gsub(".n", "", names(prepdat)[11:31]))

######@> Species list...
splist <- c("yft", "alb", "bet", "swo", "sai", "whm", "bum", "bsh",
            "spx", "bth", "sma", "ocs", "fal", "ccs")

######@> Cleaning some data...

#####@> Dataprep Brazil function...
dataprep_BR <- function(dat, splist) {
    dat <- dat %>%
        mutate(dmy = ymd(paste(setyear, setmonth, set.day, sep = "/"))) %>%
        mutate(op_yr = setyear) %>%
        mutate(op_mon = setmonth) %>%
        mutate(op_day = set.day) %>%
        mutate(hbf = hpb) %>%
        mutate(hooks = effort) %>%
        mutate(floats = hooks/hbf)
    dat$lon <- dat$set.lonb
    dat$lat <- dat$set.latb
    dat$qtr <- ceiling(as.numeric(dat$op_mon) / 3)
    dat$yrqtr <- dat$op_yr + floor((dat$op_mon - 1)/3)/4 + 0.125
    dat$latlong <- paste(dat$lat5, dat$lon5, sep = "_")
    dat$vessid <- as.factor(as.numeric(dat$boat))
    dat$tripidmon <- paste(dat$vessid, dat$op_yr, dat$op_mon)
    dat$Total <- apply(dat[, splist], 1, sum, na.rm = TRUE)
    dat$Total2 <- apply(dat[, c("bet", "yft", "alb")], 1, sum,
                        na.rm = TRUE)
    dat <- dat %>%
        dplyr::select(vessid, flag, tripidmon, dmy, op_yr, op_mon,
                      op_day, qtr, yrqtr, hooks, hbf, floats, lat, lon,
                      lat5, lon5, latlong, Total, Total2, bft, yft, alb,
                      bet, blf, skj, oth, swo, sai, whm, bum, spf, spg,
                      bsh, spx, bth, sma, ocs, fal, ccs) %>%
        as.data.frame()
    return(dat)
}

#####@> Data preparation 2...
prepdat2 <- dataprep_BR(prepdat, splist)

######@> Setup Regions...
prepdat3 <- setup_AO_regions(dat = prepdat2, regB = TRUE, regB1 = TRUE,
                             regY = TRUE, regY1 = TRUE, regY2 = TRUE)

#####@> Calling regions 0 by 1...
prepdat3$regB <- ifelse(prepdat3$regB == 0, 1, prepdat3$regB)
prepdat3$regB1 <- ifelse(prepdat3$regB1 == 0, 1, prepdat3$regB1)
prepdat3$regY <- ifelse(prepdat3$regY == 0, 1, prepdat3$regY)
prepdat3$regY1 <- ifelse(prepdat3$regY1 == 0, 1, prepdat3$regY1)
prepdat3$regY2 <- ifelse(prepdat3$regY2 == 0, 1, prepdat3$regY2)

######@> Choosing default colors...
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
col <- gg_color_hue(5)

#####@> Visualizing spatial distribution of sets...
p00 <- ggplot(data = prepdat3, aes(x = lon5, y = lat5)) +
    geom_tile(aes(fill = factor(regB)), colour = "black") +
    geom_polygon(data = mm, aes(x = long, y = lat, group = group)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    scale_fill_manual(name = "Areas",
                      values = col[1:3]) +
    coord_fixed(xlim = c(-55, 0), ylim = c(-50, 20)) +
    xlab(expression(paste("Longitude ", "(", degree, ")"))) +
    ylab(expression(paste("Latitude ", "(", degree, ")"))) +
    my_theme()
p00

p01 <- ggplot(data = prepdat3, aes(x = lon5, y = lat5)) +
    geom_tile(aes(fill = factor(regB1)), colour = "black") +
    geom_polygon(data = mm, aes(x = long, y = lat, group = group)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    scale_fill_manual(name = "Areas",
                      values = col[1:3]) +
    coord_fixed(xlim = c(-55, 0), ylim = c(-50, 20)) +
    xlab(expression(paste("Longitude ", "(", degree, ")"))) +
    ylab(expression(paste("Latitude ", "(", degree, ")"))) +
    my_theme()
p01

p02 <- ggplot(data = prepdat3, aes(x = lon5, y = lat5)) +
    geom_tile(aes(fill = factor(regY)), colour = "black") +
    geom_polygon(data = mm, aes(x = long, y = lat, group = group)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    scale_fill_manual(name = "Areas",
                      values = col[1:3]) +
    coord_fixed(xlim = c(-55, 0), ylim = c(-50, 20)) +
    xlab(expression(paste("Longitude ", "(", degree, ")"))) +
    ylab(expression(paste("Latitude ", "(", degree, ")"))) +
    my_theme()
p02

p03 <- ggplot(data = filter(prepdat3, lat5 <= 10),
              aes(x = lon5, y = lat5)) +
    geom_tile(aes(fill = factor(regY1)), colour = "black") +
    geom_polygon(data = mm, aes(x = long, y = lat, group = group)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    scale_fill_manual(name = "Areas",
                      values = col[1:3]) +
    coord_fixed(xlim = c(-55, 0), ylim = c(-50, 20)) +
    xlab(expression(paste("Longitude ", "(", degree, ")"))) +
    ylab(expression(paste("Latitude ", "(", degree, ")"))) +
    my_theme()
p03

p04 <- ggplot(data = prepdat3, aes(x = lon5, y = lat5)) +
    geom_tile(aes(fill = factor(regY2)), colour = "black") +
    geom_polygon(data = mm, aes(x = long, y = lat, group = group)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    scale_fill_manual(name = "Areas",
                      values = col[1:5]) +
    coord_fixed(xlim = c(-55, 0), ylim = c(-50, 20)) +
    xlab(expression(paste("Longitude ", "(", degree, ")"))) +
    ylab(expression(paste("Latitude ", "(", degree, ")"))) +
    my_theme()
p04

######@> Cleaning data set...

#####@> Function to do that
dataclean_BR <- function(dat, yearlim = 2019, splist) {
    for (sp in splist) {
        dat[, sp] <- as.numeric(dat[, sp])
        if (sum(is.na(dat[, sp])) > 0)
            dat[is.na(dat[, sp]), sp] <- 0
    }
    dat <- dat[!is.na(dat$hooks), ]
    dat <- dat[dat$hooks < 5000, ]
    dat <- dat[dat$hooks >= 500, ]
    dat <- dat[is.na(dat$hbf) == FALSE, ]
    dat <- dat[dat$op_yr > 1976, ]
    dat <- dat[dat$yrqtr < yearlim, ]
    dat <- dat[dat$hbf >= 3, ]
    return(dat)
}

#####@> Cleaning set...
dat <- dataclean_BR(dat = prepdat3, splist = splist)

######@> Export final data set...
save(prepdat3, dat, file = "BRdat.RData")

########################################################################
######@> Exploring dataset...

######@> Looking for some data...

#####@> Some exploratory tables - Vessels per year | Sets without
#####@> Vessels | Years without hooks...
str(dat)
summary(dat)
table(dat$vessid, dat$op_yr)
table(dat$op_yr, dat$vessid)
table(dat$op_yr, is.na(dat$vessid))
table(dat$op_yr, (dat$hooks > 0))

#####@> number of vessels per year...
xfun <- function(x) sum(x > 0)
a <- table(dat$vessid, dat$op_yr)
apply(a, 2, xfun)
apply(a, 2, lu)

#####@> Sets per day...
dev.new(width = 15, height = 9)
hist(dat$dmy, breaks = "days", freq = T, xlab = "Date",
     main = "Sets per day")
savePlot(filename = "sets_per_day.png", type = "png")
dev.off()

#####@> Plot grid squares with sets by region, for each regional
#####@> structure...
a <- unique(paste(dat$lat, dat$lon))
a0 <- dat[match(a, paste(dat$lat, dat$lon)),
          c("lat","lon","regB", "regB1", "regY", "regY1", "regY2")]

for (fld in c("regB", "regB1", "regY", "regY1", "regY2")) {
    dev.new(widath = 10, height = 10)
    reg <- with(a0, get(fld))
    plot(a0$lon, a0$lat, type = "n", xlab = "Longitude",
         ylab = "Latitude", main = fld)
    text(a0$lon, a0$lat, labels = reg, cex = 0.8, col = reg + 1)
    map(add = T)
    savePlot(paste0("map_", fld), type = "png")
    dev.off()
}

#####@> Map of hook distribution, all time...
a <- aggregate(dat$hooks, list(dat$lat5, dat$lon5), sum, na.rm = T)
dev.new(width = 11,height = 9)
symbols(x = a[,2], y = a[,1], circles = .0002 * sqrt(a[,3]),
        inches = F, bg = 2, fg = 2, xlab = "Longitude", ylab = "Latitude",
        ylim = c(-50, 20), xlim = c(-60, -10))
map(add = T, interior = F, fill = T)
savePlot(filename = "map_hooks.png", type = "png")
dev.off()

#####@> Histogram of hooks per set
table(dat$hooks[dat$hooks > 2000])
hist(dat$hooks, main = "", breaks = seq(0, 5000, 100),
     include.lowest = TRUE, right = FALSE, ylim = c(0, 15000),
     xlab = "Hooks per set", xlim = c(0, 5000))
savePlot("Hook_histogram.png", type = "png")
dev.off()

#####@> Check catch distribtions for outliers. Probably no need to remove.
table(dat$yft)
table(dat$alb)
table(dat$bet)
table(dat$swo)
table(dat$sai)
table(dat$whm)
table(dat$bum)
table(dat$bsh)
table(dat$spx)
table(dat$bth)
table(dat$sma)
table(dat$ocs)
table(dat$fal)
table(dat$ccs)

table(dat$hbf, useNA = "always")
table(dat$hbf, dat$op_yr, useNA = "always")
table(dat$op_yr, is.na(dat$hbf))

a <- table(dat$op_yr, round(dat$hbf, 0), useNA = "always")
write.csv(a, "table_hbf_by_year.csv")

#####@> Set density map by 5 degree cell
a <- log(table(dat$lon5, dat$lat5))
dev.new(width = 13, height = 10)
image(as.numeric(dimnames(a)[[1]]), as.numeric(dimnames(a)[[2]]), a,
      xlab = "Longitude", ylab = "Latitude",
      col  =  rev(heat.colors(12)))
contour(as.numeric(dimnames(a)[[1]]), as.numeric(dimnames(a)[[2]]), a,
        xlab = "Longitude", ylab = "Latitude", add  =  TRUE)
map("world", add = T, interior = T, fill = T)
savePlot("setmap_logscale.png", type = "png")
dev.off()

#####@> Mean fishing location  by yearqtr
dev.new(width = 15, height = 10)
par(mfrow = c(1, 2))
ax <- tapply(dat$yrqtr, dat$yrqtr, mean)
ay <- tapply(dat$lat5, dat$yrqtr, mean)
plot(ax, ay, xlab = "yr", ylab = "Mean latitude", type = "n")
a <- 4 * (.125+dat$yrqtr-floor(dat$yrqtr))
a <- tapply(a, dat$yrqtr, mean)
text(ax, ay, a, cex = 0.7)
ax <- tapply(dat$lon5, dat$yrqtr, mean)
ay <- tapply(dat$yrqtr, dat$yrqtr, mean)
plot(ax, ay, ylab = "yr", xlab = "Mean longitude", type = "n")
text(ax, ay, a, cex = 0.7)
savePlot("mean_fishing_location1.png",type = "png")
dev.off()

#####@> Mean fishing location by year
dev.new(width = 15, height = 10)
par(mfrow = c(1, 2))
plot(tapply(dat$op_yr, dat$op_yr, mean), tapply(dat$lat5, dat$op_yr, mean),
     xlab = "yr", ylab = "Mean latitude")
plot(tapply(dat$lon5, dat$op_yr,mean), tapply(dat$op_yr, dat$op_yr,mean),
     ylab = "yr",xlab = "Mean longitude")
savePlot("mean_fishing_location2.png", type = "png")
dev.off()

#####@> Plot hbf... Change spatial selection criteria for AO.
dev.new(20, 14)
par(mfrow = c(3, 3), mar = c(2, 2, 2, 2))
for (y in seq(1975, 2015, 5)) {
    a <- dat[dat$op_yr %in% y:(y+4),]
    ## max.hbf <- max(a$hbf, na.rm = TRUE)
    a <- tapply(a$hbf, list(a$lon5, a$lat5), mean, na.rm = T)
    image(as.numeric(dimnames(a)[[1]]), as.numeric(dimnames(a)[[2]]) ,
          a, main = y, zlim = c(0, 30), col = rev(heat.colors(30)),
          xlab = "Lon", ylab = "Lat", ylim = c(-50, 20),
          xlim = c(-60, -20))
    contour(as.numeric(dimnames(a)[[1]]), as.numeric(dimnames(a)[[2]]),
            a, add = T, levels = seq(0, 7, 1))
    map("world", add = T, interior = F, fill = T)
}
savePlot("mean_HBF.png", type = "png")
dev.off()

#####@> Plot hbf per quarter...
qqs <- c(0.125, 0.375, 0.625, 0.875)
for (qq in 1:4) {
    dev.new(20, 14)
    par(mfrow = c(3, 3), mar = c(2, 2, 2, 2), oma = c(0, 0, 1, 0))
    for (y in seq(1975, 2015, 5)) {
        a <- dat[dat$yrqtr %in% (qqs[qq]+y:(y+4)), ]
        a <- tapply(a$hbf, list(a$lon5, a$lat5), mean, na.rm = T)
        image(as.numeric(dimnames(a)[[1]]),
              as.numeric(dimnames(a)[[2]]),
              a, main = y, zlim = c(0, 30),
              col = rev(heat.colors(30)), xlab = "Lon", ylab = "Lat",
              ylim = c(20, -50), xlim = c(-60, -20))
        contour(as.numeric(dimnames(a)[[1]]),
                as.numeric(dimnames(a)[[2]]), a, add = T,
                levels = seq(0, 15, 1), col = "blue")
        map("world", add = T, interior = F, fill = T)
    }
    title(paste("Quarter",qq), outer = T, line = 0)
    savePlot(paste0("mean_HBF_q",qq,".png"), type = "png")
    dev.off()
}

write.csv(table(dat$lat5, dat$lon5), file = "ops_by_lat-long.csv")
write.csv(table(dat$lat5, dat$lon5, 5 * floor(dat$yrqtr/5)),
          file = "ops_by_lat-long-5yr.csv")

#####@> Exploration Regressions tree...

####@> preparing data...
a <- dat
head(dat)
dim(a)

###@> nominal cpue...
a$yftcpue <- a$yft/a$hooks
a$albcpue <- a$alb/a$hooks
a$betcpue <- a$bet/a$hooks
a$swocpue <- a$swo/a$hooks
a$saicpue <- a$sai/a$hooks
a$whmcpue <- a$whm/a$hooks
a$bumcpue <- a$bum/a$hooks
a$bshcpue <- a$bsh/a$hooks
a$spxcpue <- a$spx/a$hooks
a$bthcpue <- a$bth/a$hooks
a$smacpue <- a$sma/a$hooks
a$ocscpue <- a$ocs/a$hooks
a$falcpue <- a$fal/a$hooks
a$ccscpue <- a$ccs/a$hooks

###@> simple model 01...
simplemod01 <- rpart(a$albcpue ~ a$lon + a$lat + a$yrqtr + a$yftcpue +
                         a$betcpue + a$swocpue + a$saicpue + a$whmcpue +
                         a$bumcpue + a$bshcpue + a$spxcpue + a$bthcpue +
                         a$smacpue + a$ocscpue + a$falcpue + a$ccscpue)

dev.new(width = 11, height = 7)
plot(simplemod01)
text(simplemod01)
savePlot("Rpart_alb_cpue_full", type = "png")
dev.off()

simplemod02 <- rpart(a$albcpue ~ a$lon + a$lat + a$yrqtr + a$bshcpue +
                         a$yftcpue + a$swocpue + a$betcpue)

dev.new(width = 11, height = 7)
plot(simplemod02)
text(simplemod02)
savePlot("Rpart_alb_cpue_subset", type = "png")
dev.off()

#####@> Exploration with Random Forest...

####@> These take a long time and use a lot of memory, but are useful...
system.time(
    simplefor <- randomForest(albcpue ~ lon + lat + yrqtr + hbf + bshcpue +
                                  yftcpue + swocpue + betcpue, data = a)
)

print(simplefor)

dev.new(width = 11, height = 7)
plot(simplefor)

varImpPlot(simplefor)
savePlot("Rforest_alb_cpue", type = "png")
dev.off()

partialPlot(simplefor, pred.data = a, x.var = "hbf")
savePlot("Rforest_alb_cpue_partial", type = "png")
dev.off()

########################################################################
######@> Cluster Analyses...

######@> Change the directory...
clustdir <- "../clustering/"
## dir.create(clustdir)
setwd(clustdir)

######@> Preparing the data...
gc(reset = TRUE)
br_splist <- c("yft", "alb", "bet", "swo", "sai", "whm", "bum", "bsh",
               "spx", "bth", "sma", "ocs", "fal", "ccs")
use_splist <- c("yft", "alb", "bet", "swo", "sai", "bum", "bsh", "whm",
                "sma")
allabs <- c("vessid", "yrqtr", "latlong", "op_yr", "hbf", "hooks",
            "tripidmon", use_splist, "Total", "lat", "lon", "lat5",
            "lon5", "regB", "regB1", "regY", "regY1", "regY2")
dat <- data.frame(dat)

#####@> Number of albacore clusters. Will need to be adjusted for each
#####@> fleet...
nclB <- c(4, 4, 4)
flag <- "BR"
cvn <- c("yrqtr", "latlong", "hooks", "hbf", "vessid", "Total", "lat",
         "lon", "lat5", "lon5", "op_yr", "tripidmon", "regB", "regB1",
         "regY", "regY1", "regY2")
allabs <- c("vessid", "yrqtr", "latlong", "op_yr", "hbf", "hooks",
            "tripidmon", use_splist, "Total", "lat", "lon", "lat5",
            "lon5", "regB", "regB1", "regY", "regY1", "regY2")

#####@> Looping for the analysis of the species - regB region...
for (r in unique(dat$regB)) {
    dev.new(15, 12)
    par(mfrow = c(5, 3), mar = c(3, 2, 2, 1), oma = c(0, 0, 2, 0))
    a <- dat[dat$regB == r, ]
    for (sp in br_splist) {
        plot(sort(unique(a$yrqtr)), tapply(a[, sp], a$yrqtr, mean), main = sp)
        title(paste("Region", r ), outer = TRUE)
        savePlot(filename = paste("freq", flag, "regB_Region", r, sep = "_"),
                 type = "png")
    }
    dev.off()
}

#####@> Looping for the analysis of the species - regY1 region...
for (r in unique(dat$regY1)) {
    dev.new(15, 12)
    par(mfrow = c(5, 3), mar = c(3, 2, 2, 1), oma = c(0, 0, 2, 0))
    a <- dat[dat$regY1 == r, ]
    for (sp in br_splist) {
        plot(sort(unique(a$yrqtr)), tapply(a[, sp], a$yrqtr, mean), main = sp)
        title(paste("Region", r ), outer = TRUE)
        savePlot(filename = paste("freq", flag, "regY1_Region", r, sep = "_"),
                 type = "png")
    }
    dev.off()
}

#####@> Looping for the analysis of the species - regY2 region...
for (r in unique(dat$regY2)) {
    dev.new(15, 12)
    par(mfrow = c(5, 3), mar = c(3, 2, 2, 1), oma = c(0, 0, 2, 0))
    a <- dat[dat$regY2 == r, ]
    for (sp in br_splist) {
        plot(sort(unique(a$yrqtr)), tapply(a[, sp], a$yrqtr, mean), main = sp)
        title(paste("Region", r ), outer = TRUE)
        savePlot(filename = paste("freq", flag, "regY2_Region", r, sep = "_"),
                 type = "png")
    }
    dev.off()
}

#####@> Cluster analyses for regB region...
regtype <- "regB"
for (r in unique(dat$regB)) {
    fnh <- paste(flag, regtype, r, sep = "_")
    dataset <- clust_PCA_run(r = r, ddd = dat, allsp = use_splist,
                             allabs = allabs, regtype = regtype,
                             ncl = nclB[r], plotPCA = FALSE,
                             clustid = "tripidmon",
                             allclust = FALSE, ll5 = TRUE, flag = flag,
                             fnhead = fnh, covarnames = cvn)
    save(dataset, file = paste0(fnh, ".RData"))
}

#####@> Cluster analyses for regB1 region...
regtype <- "regB1"
for (r in unique(dat$regB1)) {
    fnh <- paste(flag, regtype, r, sep = "_")
    dataset <- clust_PCA_run(r = r, ddd = dat, allsp = use_splist,
                             allabs = allabs, regtype = regtype,
                             ncl = nclB[r], plotPCA = FALSE,
                             clustid = "tripidmon",
                             allclust = FALSE, ll5 = TRUE, flag = flag,
                             fnhead = fnh, covarnames = cvn)
    save(dataset, file = paste0(fnh, ".RData"))
}

#####@> Cluster analyses for regY1 region...
regtype <- "regY1"
for (r in unique(dat$regY1)) {
    fnh <- paste(flag, regtype, r, sep = "_")
    dataset <- clust_PCA_run(r = r, ddd = dat, allsp = use_splist,
                             allabs = allabs, regtype = regtype,
                             ncl = nclB[r], plotPCA = FALSE,
                             clustid = "tripidmon",
                             allclust = FALSE, ll5 = TRUE, flag = flag,
                             fnhead = fnh, covarnames = cvn)
    save(dataset, file = paste0(fnh, ".RData"))
}

#####@> Cluster analyses for all regions integrated...
dat2 <- dat
dat2$reg <- 1
nclB <- 4
cvn <- c("yrqtr", "latlong", "hooks", "hbf", "vessid", "Total", "lat",
         "lon", "lat5", "lon5", "op_yr", "tripidmon", "regB", "regB1",
         "regY", "regY1", "reg")
allabs <- c("vessid", "yrqtr", "latlong", "op_yr", "hbf", "hooks",
            "tripidmon", use_splist, "Total", "lat", "lon", "lat5",
            "lon5", "regB", "regB1", "regY", "regY1", "reg")

regtype <- "reg"
for (r in unique(dat2$reg)) {
    fnh <- paste(flag, regtype, "All_Regions_Integrated", r, sep = "_")
    dataset <- clust_PCA_run(r = r, ddd = dat2, allsp = use_splist,
                             allabs = allabs, regtype = regtype,
                             ncl = nclB[r], plotPCA = FALSE,
                             clustid = "tripidmon",
                             allclust = FALSE, ll5 = TRUE, flag = flag,
                             fnhead = fnh, covarnames = cvn)
    save(dataset, file = paste0(fnh, "All_Regions_Integrated", ".RData"))
}

########################################################################
######@> Standardization CPUE...

######@> Defining the new folder for the next analyses...
resdir <- "../analyses/std_cl_BRonly_hbf"
dir.create(resdir)
setwd(resdir)

####@> defining the projdir...
projdir <- "~/Github/CPUE_STD_ALB_BRA_2020/ICCAT/2020_Albacore/"

######@> Defining the species list...
splist <- c("alb", "bet", "yft", "swo", "bsh", "whm", "sai", "sma",
            "bum")

######@> Variables to standardization...
stdlabs <- c("vessid", "yrqtr", "latlong", "op_yr", "op_mon", "hbf",
             "hooks", splist, "hcltrp", "Total", "lat", "lon", "lat5",
             "lon5", "regB", "regY1", "flag")

######@> Defining the clusters for the standardization...

#####@> Joint standardization...
## clkeepJP_Y <- list("alb" = list(c(1, 2, 4), c(1, 2, 3, 4), c(1, 2, 3)))
## clkeepKR_Y <- list("alb" = list(c(0), c(1, 2, 3, 4), c(1, 2, 3)))
clkeepBR_A <- list("alb" = list(c(1, 2, 3, 4), c(1, 2, 3, 4), c(1, 2, 3, 4)))
## clkeepTW_Y <- list("alb" = list(c(4), c(2, 3), c(0)))
## clkeepUS_Y <- list("alb" = list(c(2, 3), c(1, 3), c(0)))
## clk_Y <- list(JP = clkeepJP_Y, KR = clkeepKR_Y, BR = clkeepBR_Y,
##               TW = clkeepTW_Y, US = clkeepUS_Y)

######@> Brazil only...
clk_A <- list(BR = clkeepBR_A)

######@> Defining the parameters...
runpars <- list()
runpars[["alb"]] <- list(regtype = "regY1", regtype2 = "Y1", clk = clk_A,
                         doregs = 2:3, addcl = TRUE, dohbf = TRUE,
                         cltype = "hcltrp")
runreg = 1; runsp = "alb"
keepd = TRUE; maxyr = 2018; maxqtrs = 200; minqtrs_byreg = c(5, 5, 5)

######@> Running standardization by hand...

#####@> Region Y1 02...

####@> Loading datasets from cluster folder...
load("../../clustering/BR_regY1_2.RData")
jdat02 <- dataset; rm(dataset)

####@> watching the dataset for the proportions of catches of alb per
####@> vessel...
tab <- jdat02 %>%
    group_by(vessid) %>%
    summarise(alb = sum(alb, na.rm = TRUE),
              total = sum(Total, na.rm = TRUE)) %>%
    mutate(prop = alb/total) %>%
    arrange(-prop) %>%
    as.data.frame()
table(tab$prop == 0)

####@> id for boats that never caught a alb in life...
id <- tab$vessid[tab$prop == 0]

####@> removing the boats that never catched a alb in life...
jdat02.cut <- filter(jdat02, !vessid %in% id)

####@> Defining data for glm...
glmdat <- select_data_JointIO(jdat02.cut,
                              runreg = 2,
                              clk = clk_A,
                              minqtrs = 2,
                              runsp = "alb",
                              mt = "deltabin",
                              vars = c("vessid", "hooks", "yrqtr",
                                       "latlong", "hbf"),
                              oneflag = "BR",
                              maxqtrs = 500,
                              minvess = 5,
                              minll = 5,
                              minyrqtr = 50,
                              yrlims = c(1998, 2019),
                              addcl = TRUE,
                              cltype = "hcltrp",
                              addpca = NA,
                              samp = NA,
                              strsmp = NA)

#####@> Coverage in region 02...

####@> Prepdat3...
tmp01 <- prepdat3 %>%
    filter(regY1 == 2) %>%
    filter(op_yr %in% 1998:2018) %>%
    group_by(op_yr) %>%
    summarise(Total = n())

####@> Jdat02.cut...
tmp02 <- jdat02.cut %>%
    filter(regY1 == 2) %>%
    filter(op_yr %in% 1998:2018) %>%
    group_by(op_yr) %>%
    summarise(Total2 = n())

####@> Modelo...
tmp03 <- glmdat %>%
    mutate(op_yr = as.numeric(substr(yrqtr, 1, 4))) %>%
    filter(op_yr %in% 1998:2018) %>%
    group_by(op_yr) %>%
    summarise(Amostra = n())

#####@> Merging...
tmp <- tmp01 %>%
    full_join(tmp02, by = "op_yr") %>%
    full_join(tmp03, by = "op_yr") %>%
    mutate(coverage01 = (Amostra / Total) * 100,
           coverage02 = (Amostra / Total2) * 100)

p00 <- ggplot(tmp, aes(x = op_yr, y = coverage01)) +
    geom_line() +
    geom_label(aes(label = paste(round(coverage01, 1), "%"))) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_x_continuous(limits = c(1998, 2018),
                       breaks = seq(1998, 2018, 2)) +
    geom_hline(yintercept = mean(tmp$coverage01), linetype = "dashed") +
    annotate(geom = "text", x = 2017.8, y = mean(tmp$coverage01) - 1,
             label = "Average Coverage") +
    labs(x = "Year", y = "Coverage (%)", caption = "Coverage R02") +
    my_theme()
p00

png("Coverage_R02_ver00.png", units = "cm", res = 200,
    w = 35, h = 20)
print(p00)
dev.off()

#####@> LognC Models...

####@> Defining the constant...
mn <- with(glmdat, 0.1 * mean(get(runsp)/hooks))

####@> Defining the weights per area...
wtt.all <- mk_wts(glmdat, wttype = "area")

####@> Running the model for no vessels included...

###@> Defining the formula...
form01 <- make_formula_IO(runsp, modtype = "logn", dohbf = TRUE,
                          addboat = FALSE, addcl = TRUE, nhbf = 1)

###@> Defining the labs for the models outputs...
modlab <- "lognC_novess_allyrs"
fname <- "BR_regY1_R02"

###@> Running the first model...
model <- glm(form01, data = glmdat, weights = wtt.all,
             family = "gaussian")

###@> Summary the results...
summarize_and_store(mod = model, dat = glmdat, fname, modlab,
                    dohbf = TRUE, keepd = TRUE)

####@> Running the model for vessels included...

###@> Defining the formula...
form02 <- make_formula_IO(runsp, modtype = "logn", dohbf = TRUE,
                          addboat = TRUE, addcl = TRUE, nhbf = 1)

###@> Defining the labs for the models outputs...
modlab <- "lognC_boat_allyrs"
fname <- "BR_regY1_R02"

###@> Running the first model...
model <- glm(form02, data = glmdat, weights = wtt.all,
             family = "gaussian")

###@> Summary the results...
summarize_and_store(mod = model, dat = glmdat, fname, modlab,
                    dohbf = TRUE, keepd = TRUE)

#####@> Delta-LogNormal Models...

####@> Running the model for no vessels included...

###@> Defining the labs for the models outputs...
modlab <- "dellog_novess_allyrs"
fname <- "BR_regY1_R02"

###@> Running the model...
do_deltalog(dat = glmdat, dohbf = TRUE, addboat = FALSE, addcl = TRUE,
            nhbf = 1, runsp = "alb", fname = fname, modlab = modlab,
            keepd = TRUE)

####@> Running the model for vessels included...

###@> Defining the labs for the models outputs...
modlab <- "dellog_boat_allyrs"
fname <- "BR_regY1_R02"

###@> Running the model...
do_deltalog(dat = glmdat, dohbf = TRUE, addboat = TRUE, addcl = TRUE,
            nhbf = 1, runsp = "alb", fname = fname, modlab = modlab,
            keepd = TRUE)

#####@> Region Y1 03...

####@> Loading datasets from cluster folder...
load("../../clustering/BR_regY1_3.RData")
jdat03 <- dataset; rm(dataset)

####@> watching the dataset for the proportions of catches of alb per
####@> vessel...
tab <- jdat03 %>%
    group_by(vessid) %>%
    summarise(alb = sum(alb, na.rm = TRUE),
              total = sum(Total, na.rm = TRUE)) %>%
    mutate(prop = alb/total) %>%
    arrange(-prop) %>%
    as.data.frame()
table(tab$prop == 0)

####@> id for boats that never caught a alb in life...
id <- tab$vessid[tab$prop == 0]

####@> removing the boats that never catched a alb in life...
jdat03.cut <- filter(jdat03, !vessid %in% id)

####@> Defining data for glm...
glmdat <- select_data_JointIO(jdat03.cut,
                              runreg = 3,
                              clk = clk_A,
                              minqtrs = 2,
                              runsp = "alb",
                              mt = "deltabin",
                              vars = c("vessid", "hooks", "yrqtr",
                                       "latlong", "hbf"),
                              oneflag = "BR",
                              maxqtrs = 500,
                              minvess = 5,
                              minll = 5,
                              minyrqtr = 50,
                              addcl = TRUE,
                              yrlims = c(1998, 2019),
                              cltype = "hcltrp",
                              addpca = NA,
                              samp = NA,
                              strsmp = NA)

#####@> Coverage in region 03...

####@> Prepdat3...
tmp01 <- prepdat3 %>%
    filter(regY1 == 3) %>%
    filter(op_yr %in% 1998:2018) %>%
    group_by(op_yr) %>%
    summarise(Total = n())

####@> Jdat03.cut...
tmp02 <- jdat03.cut %>%
    filter(regY1 == 3) %>%
    filter(op_yr %in% 1998:2018) %>%
    group_by(op_yr) %>%
    summarise(Total2 = n())

####@> Modelo...
tmp03 <- glmdat %>%
    mutate(op_yr = as.numeric(substr(yrqtr, 1, 4))) %>%
    filter(op_yr %in% 1998:2018) %>%
    group_by(op_yr) %>%
    summarise(Amostra = n())

#####@> Merging...
tmp <- tmp01 %>%
    full_join(tmp02, by = "op_yr") %>%
    full_join(tmp03, by = "op_yr") %>%
    mutate(coverage01 = (Amostra / Total) * 100,
           coverage02 = (Amostra / Total2) * 100)

p00 <- ggplot(tmp, aes(x = op_yr, y = coverage01)) +
    geom_line() +
    geom_label(aes(label = paste(round(coverage01, 1), "%"))) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_x_continuous(limits = c(1998, 2018),
                       breaks = seq(1998, 2018, 2)) +
    geom_hline(yintercept = mean(tmp$coverage01), linetype = "dashed") +
    annotate(geom = "text", x = 2017.8, y = mean(tmp$coverage01) - 1,
             label = "Average Coverage") +
    labs(x = "Year", y = "Coverage (%)", caption = "Coverage R03") +
    my_theme()
p00

png("Coverage_R03_ver00.png", units = "cm", res = 200,
    w = 35, h = 20)
print(p00)
dev.off()

#####@> LognC Models...

####@> Defining the constant...
mn <- with(glmdat, 0.1 * mean(get(runsp)/hooks))

####@> Defining the weights per area...
wtt.all <- mk_wts(glmdat, wttype = "area")

####@> Running the model for no vessels included...

###@> Defining the formula...
form01 <- make_formula_IO(runsp, modtype = "logn", dohbf = TRUE,
                          addboat = FALSE, addcl = TRUE, nhbf = 1)

###@> Defining the labs for the models outputs...
modlab <- "lognC_novess_allyrs"
fname <- "BR_regY1_R03"

###@> Running the first model...
model <- glm(form01, data = glmdat, weights = wtt.all,
             family = "gaussian")

###@> Summary the results...
summarize_and_store(mod = model, dat = glmdat, fname, modlab,
                    dohbf = TRUE, keepd = TRUE)

####@> Running the model for vessels included...

###@> Defining the formula...
form02 <- make_formula_IO(runsp, modtype = "logn", dohbf = TRUE,
                          addboat = TRUE, addcl = TRUE, nhbf = 1)

###@> Defining the labs for the models outputs...
modlab <- "lognC_boat_allyrs"
fname <- "BR_regY1_R03"

###@> Running the first model...
model <- glm(form02, data = glmdat, weights = wtt.all,
             family = "gaussian")

###@> Summary the results...
summarize_and_store(mod = model, dat = glmdat, fname, modlab,
                    dohbf = TRUE, keepd = TRUE)

#####@> Delta-LogNormal Models...

####@> Running the model for no vessels included...

###@> Defining the labs for the models outputs...
modlab <- "dellog_novess_allyrs"
fname <- "BR_regY1_R03"

###@> Running the model...
do_deltalog(dat = glmdat, dohbf = TRUE, addboat = FALSE, addcl = TRUE,
            nhbf = 1, runsp = "alb", fname = fname, modlab = modlab,
            keepd = TRUE)

####@> Running the model for vessels included...

###@> Defining the labs for the models outputs...
modlab <- "dellog_boat_allyrs"
fname <- "BR_regY1_R03"

###@> Running the model...
do_deltalog(dat = glmdat, dohbf = TRUE, addboat = TRUE, addcl = TRUE,
            nhbf = 1, runsp = "alb", fname = fname, modlab = modlab,
            keepd = TRUE)

#####@> Regions Integrated in One...

####@> Loading datasets from cluster folder...
load("../../clustering/BR_reg_All_Regions_Integrated_1All_Regions_Integrated.RData")
jdatAll <- dataset; rm(dataset)

####@> watching the dataset for the proportions of catches of alb per
####@> vessel...
tab <- jdatAll %>%
    group_by(vessid) %>%
    summarise(alb = sum(alb, na.rm = TRUE),
              total = sum(Total, na.rm = TRUE)) %>%
    mutate(prop = alb/total) %>%
    arrange(-prop) %>%
    as.data.frame()
table(tab$prop == 0)

####@> id for boats that never caught a alb in life...
id <- tab$vessid[tab$prop == 0]

####@> removing the boats that never catched a alb in life...
jdatAll.cut <- filter(jdatAll, !vessid %in% id)

####@> Defining data for glm...
glmdat <- select_data_JointIO(jdatAll.cut,
                              runreg = 1,
                              clk = clk_A,
                              minqtrs = 2,
                              runsp = "alb",
                              mt = "deltabin",
                              vars = c("vessid", "hooks", "yrqtr",
                                       "latlong", "hbf"),
                              oneflag = "BR",
                              maxqtrs = 500,
                              minvess = 5,
                              minll = 5,
                              minyrqtr = 50,
                              addcl = TRUE,
                              yrlims = c(1998, 2019),
                              cltype = "hcltrp",
                              addpca = NA,
                              samp = NA,
                              strsmp = NA)

#####@> Table for comparison of the coverage...
tab01 <- dat %>%
    group_by(op_yr) %>%
    summarise(alb = sum(alb, na.rm = TRUE)) %>%
    filter(op_yr %in% 1998:2018)

tab02 <- glmdat %>%
    mutate(op_yr = floor(as.numeric(as.character(yrqtr)))) %>%
    group_by(op_yr) %>%
    summarise(alb = sum(alb, na.rm = TRUE)) %>%
    filter(op_yr %in% 1998:2018)

tab03 <- merge(tab01, tab02, by = "op_yr")

tab03$prop <- with(tab03, alb.y / alb.x)

write.table(tab03, file = "BR_Coverage.csv", sep = ",", dec = ".",
            row.names = FALSE)

#####@> LognC Models...

####@> Defining the constant...
mn <- with(glmdat, 0.1 * mean(get(runsp)/hooks))

####@> Defining the weights per area...
wtt.all <- mk_wts(glmdat, wttype = "area")

####@> Running the model for no vessels included...

###@> Defining the formula...
form01 <- make_formula_IO(runsp, modtype = "logn", dohbf = TRUE,
                          addboat = FALSE, addcl = TRUE, nhbf = 1)

###@> Defining the labs for the models outputs...
modlab <- "lognC_novess_allyrs"
fname <- "BR_Reg_All_Integrated"

###@> Running the first model...
model <- glm(form01, data = glmdat, weights = wtt.all,
             family = "gaussian")

###@> Summary the results...
summarize_and_store(mod = model, dat = glmdat, fname, modlab,
                    dohbf = TRUE, keepd = TRUE)

####@> Running the model for vessels included...

###@> Defining the formula...
form02 <- make_formula_IO(runsp, modtype = "logn", dohbf = TRUE,
                          addboat = TRUE, addcl = TRUE, nhbf = 1)

###@> Defining the labs for the models outputs...
modlab <- "lognC_boat_allyrs"
fname <- "BR_Reg_All_Integrated"

###@> Running the first model...
model <- glm(form02, data = glmdat, weights = wtt.all,
             family = "gaussian")

###@> Summary the results...
summarize_and_store(mod = model, dat = glmdat, fname, modlab,
                    dohbf = TRUE, keepd = TRUE)

#####@> Delta-LogNormal Models...

####@> Running the model for no vessels included...

###@> Defining the labs for the models outputs...
modlab <- "dellog_novess_allyrs"
fname <- "BR_Reg_All_Integrated"

###@> Running the model...
do_deltalog(dat = glmdat, dohbf = TRUE, addboat = FALSE, addcl = TRUE,
            nhbf = 1, runsp = "alb", fname = fname, modlab = modlab,
            keepd = TRUE)

####@> Running the model for vessels included...

###@> Defining the labs for the models outputs...
modlab <- "dellog_boat_allyrs"
fname <- "BR_Reg_All_Integrated"

###@> Running the model...
do_deltalog(dat = glmdat, dohbf = TRUE, addboat = TRUE, addcl = TRUE,
            nhbf = 1, runsp = "alb", fname = fname, modlab = modlab,
            keepd = TRUE)

########################################################################
######@> Plot Results and Diagnostics...

######@> Creating a function to evaluate the results...
doplot_cpue <- function(a, vartype, mdti, regstr, runreg) {
    plot(a$yq, a$pr/mean(a$pr, na.rm = TRUE), xlab = "Year-quarter",
         ylab = "Relative CPUE", main = paste(vartype, mdti),
         type = "l", ylim = c(0, 3))
    points(a$yq, a$pr/mean(a$pr, na.rm = TRUE), cex = 0.7)
    points(a$yq, a$ul/mean(a$pr, na.rm = TRUE), pch = "-", col = 3)
    points(a$yq, a$ll/mean(a$pr, na.rm = TRUE), pch = "-", col = 2)
    mtext(paste0(regstr, " R", runreg), side = 3, outer = TRUE,
          line = -2)
}

######@> Creating a new directory to receive the outputs...
outdir <- "outputs"
dir.create(outdir)
setwd(outdir)

######@> Configurations to load and run diagnostics...

#####@> Model lognC - Region Y1 02 - No Vessels ID...
mdt <- "novess_allyrs"
mdti <- "1998 - present no vessid"
vartype <- "lognC"
regstr <- "regY1"
runreg <- "02"
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, "_R", runreg)

#####@> Loading the data...
load(paste0("../", fname, "_", modtype, "_predictions.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(as.character(nd$newdat$yrqtr)))
xx$pr1 <- switch(vartype, lognC = exp(nd$predresp$fit),
                 negbC = nd$predresp$fit)

####@> converting the predictions to index...
if(vartype == "lognC") {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
} else {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
}

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq), max(xx$yq), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 3:6])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

#####@> Model lognC - Region Y1 03 - No Vessels ID...
mdt <- "novess_allyrs"
mdti <- "1998 - present no vessid"
vartype <- "lognC"
regstr <- "regY1"
runreg <- "03"
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, "_R", runreg)

#####@> Loading the data...
load(paste0("../", fname, "_", modtype, "_predictions.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(as.character(nd$newdat$yrqtr)))
xx$pr1 <- switch(vartype, lognC = exp(nd$predresp$fit),
                 negbC = nd$predresp$fit)

####@> converting the predictions to index...
if(vartype == "lognC") {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
} else {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
}

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq), max(xx$yq), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 3:6])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

#####@> Model lognC - Region Y1 02 - Vessels ID...
mdt <- "boat_allyrs"
mdti <- "1998 - present with vessid"
vartype <- "lognC"
regstr <- "regY1"
runreg <- "02"
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, "_R", runreg)

#####@> Loading the data...
load(paste0("../", fname, "_", modtype, "_predictions.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(as.character(nd$newdat$yrqtr)))
xx$pr1 <- switch(vartype, lognC = exp(nd$predresp$fit),
                 negbC = nd$predresp$fit)

####@> converting the predictions to index...
if(vartype == "lognC") {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
} else {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
}

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq), max(xx$yq), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 3:6])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

#####@> Model lognC - Region Y1 03 - No Vessels ID...
mdt <- "boat_allyrs"
mdti <- "1998 - present with vessid"
vartype <- "lognC"
regstr <- "regY1"
runreg <- "03"
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, "_R", runreg)

#####@> Loading the data...
load(paste0("../", fname, "_", modtype, "_predictions.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(as.character(nd$newdat$yrqtr)))
xx$pr1 <- switch(vartype, lognC = exp(nd$predresp$fit),
                 negbC = nd$predresp$fit)

####@> converting the predictions to index...
if(vartype == "lognC") {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
} else {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
}

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq), max(xx$yq), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 3:6])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

#####@> Model Delta LogNormal - Region Y1 02 - No Vessels ID...
mdt <- "novess_allyrs"
mdti <- "1998 - present no vessid"
vartype <- "dellog"
regstr <- "regY1"
runreg <- "02"
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, "_R", runreg)

#####@> Loading the data...
load(paste0("../", fname, "_pos_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_bin_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_", modtype, "_indices.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(gsub("yrqtr", "", names(coefs.pos))))
xx <- data.frame(yq = as.numeric(as.character(ndpos$newdat$yrqtr)))
xx$pr <- pcoefs
xx$ln.cv <- ndpos$predterms$se.fit[, 1]
xx$ll <- exp(log(pcoefs) - 1.96 * ndpos$predterms$se.fit[,1])
xx$ul <- exp(log(pcoefs) + 1.96 * ndpos$predterms$se.fit[,1])

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq, na.rm = TRUE),
                         max(xx$yq, na.rm = TRUE), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 2:5])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

#####@> Model Delta LogNormal - Region Y1 03 - No Vessels ID...
mdt <- "novess_allyrs"
mdti <- "1998 - present no vessid"
vartype <- "dellog"
regstr <- "regY1"
runreg <- "03"
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, "_R", runreg)

#####@> Loading the data...
load(paste0("../", fname, "_pos_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_bin_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_", modtype, "_indices.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(gsub("yrqtr", "", names(coefs.pos))))
xx <- data.frame(yq = as.numeric(as.character(ndpos$newdat$yrqtr)))
xx$pr <- pcoefs
xx$ln.cv <- ndpos$predterms$se.fit[, 1]
xx$ll <- exp(log(pcoefs) - 1.96 * ndpos$predterms$se.fit[,1])
xx$ul <- exp(log(pcoefs) + 1.96 * ndpos$predterms$se.fit[,1])

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq, na.rm = TRUE),
                         max(xx$yq, na.rm = TRUE), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 2:5])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

#####@> Model Delta LogNormal - Region Y1 02 - Vessels ID...
mdt <- "boat_allyrs"
mdti <- "1998 - present with vessid"
vartype <- "dellog"
regstr <- "regY1"
runreg <- "02"
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, "_R", runreg)

#####@> Loading the data...
load(paste0("../", fname, "_pos_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_bin_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_", modtype, "_indices.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(gsub("yrqtr", "", names(coefs.pos))))
xx <- data.frame(yq = as.numeric(as.character(ndpos$newdat$yrqtr)))
xx$pr <- pcoefs
xx$ln.cv <- ndpos$predterms$se.fit[, 1]
xx$ll <- exp(log(pcoefs) - 1.96 * ndpos$predterms$se.fit[,1])
xx$ul <- exp(log(pcoefs) + 1.96 * ndpos$predterms$se.fit[,1])

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq, na.rm = TRUE),
                         max(xx$yq, na.rm = TRUE), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 2:5])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

#####@> Model Delta LogNormal - Region Y1 03 - Vessels ID...
mdt <- "boat_allyrs"
mdti <- "1998 - present with vessid"
vartype <- "dellog"
regstr <- "regY1"
runreg <- "03"
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, "_R", runreg)

#####@> Loading the data...
load(paste0("../", fname, "_pos_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_bin_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_", modtype, "_indices.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(gsub("yrqtr", "", names(coefs.pos))))
xx <- data.frame(yq = as.numeric(as.character(ndpos$newdat$yrqtr)))
xx$pr <- pcoefs
xx$ln.cv <- ndpos$predterms$se.fit[, 1]
xx$ll <- exp(log(pcoefs) - 1.96 * ndpos$predterms$se.fit[,1])
xx$ul <- exp(log(pcoefs) + 1.96 * ndpos$predterms$se.fit[,1])

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq, na.rm = TRUE),
                         max(xx$yq, na.rm = TRUE), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 2:5])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

#####@> Model lognC - Integrated Regions - No Vessels ID...
mdt <- "novess_allyrs"
mdti <- "1998 - present no vessid"
vartype <- "lognC"
regstr <- "Reg_All_Integrated"
runreg <- ""
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, runreg)

#####@> Loading the data...
load(paste0("../", fname, "_", modtype, "_predictions.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(as.character(nd$newdat$yrqtr)))
xx$pr1 <- switch(vartype, lognC = exp(nd$predresp$fit),
                 negbC = nd$predresp$fit)

####@> converting the predictions to index...
if(vartype == "lognC") {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
} else {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
}

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq), max(xx$yq), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 3:6])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

#####@> Model lognC - All Integreted Regions - Vessels ID...
mdt <- "boat_allyrs"
mdti <- "1998 - present with vessid"
vartype <- "lognC"
regstr <- "Reg_All_Integrated"
runreg <- ""
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, runreg)

#####@> Loading the data...
load(paste0("../", fname, "_", modtype, "_predictions.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(as.character(nd$newdat$yrqtr)))
xx$pr1 <- switch(vartype, lognC = exp(nd$predresp$fit),
                 negbC = nd$predresp$fit)

####@> converting the predictions to index...
if(vartype == "lognC") {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
} else {
    xx$pr <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant"))
    xx$cv <- nd$predterms$se.fit[, 1]
    xx$ll <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") -
                 1.96 * nd$predterms$se.fit[, 1])
    xx$ul <- exp(apply(nd$predterms$fit, 1, sum) +
                 attr(nd$predterms$fit, "constant") +
                 1.96 * nd$predterms$se.fit[, 1])
}

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq), max(xx$yq), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 3:6])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

#####@> Model Delta LogNormal - Integrated Regions - No Vessels ID...
mdt <- "novess_allyrs"
mdti <- "1998 - present no vessid"
vartype <- "dellog"
regstr <- "Reg_All_Integrated"
runreg <- ""
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, runreg)

#####@> Loading the data...
load(paste0("../", fname, "_pos_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_bin_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_", modtype, "_indices.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(gsub("yrqtr", "", names(coefs.pos))))
xx <- data.frame(yq = as.numeric(as.character(ndpos$newdat$yrqtr)))
xx$pr <- pcoefs
xx$ln.cv <- ndpos$predterms$se.fit[, 1]
xx$ll <- exp(log(pcoefs) - 1.96 * ndpos$predterms$se.fit[,1])
xx$ul <- exp(log(pcoefs) + 1.96 * ndpos$predterms$se.fit[,1])

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq, na.rm = TRUE),
                         max(xx$yq, na.rm = TRUE), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 2:5])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

#####@> Model Delta LogNormal - All Integrated Regions - Vessels ID...
mdt <- "boat_allyrs"
mdti <- "1998 - present with vessid"
vartype <- "dellog"
regstr <- "Reg_All_Integrated"
runreg <- ""
modtype <- paste(vartype, mdt, sep = "_")
fname <- paste0("BR_", regstr, runreg)

#####@> Loading the data...
load(paste0("../", fname, "_pos_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_bin_", modtype, "_predictions.RData"))
load(paste0("../", fname, "_", modtype, "_indices.RData"))

#####@> Preparing the prediction matrix to plot...

####@> converting the prediction to index...
xx <- data.frame(yq = as.numeric(gsub("yrqtr", "", names(coefs.pos))))
xx <- data.frame(yq = as.numeric(as.character(ndpos$newdat$yrqtr)))
xx$pr <- pcoefs
xx$ln.cv <- ndpos$predterms$se.fit[, 1]
xx$ll <- exp(log(pcoefs) - 1.96 * ndpos$predterms$se.fit[,1])
xx$ul <- exp(log(pcoefs) + 1.96 * ndpos$predterms$se.fit[,1])

####@> creating a new dataframe with all years and quarters...
a <- data.frame(yq = seq(min(xx$yq, na.rm = TRUE),
                         max(xx$yq, na.rm = TRUE), 0.25))
a <- cbind(yq = a$yq, xx[match(a$yq, xx$yq), 2:5])

####@> scaling the index by the mean...
a[, c(2, 4, 5)] <- a[, c(2, 4, 5)]/mean(a[, 2], na.rm = TRUE)

####@> convert year and quarter in factors...
a$yr <- as.factor(floor(a$yq))
a$qtr <- as.factor(a$yq - floor(a$yq))

####@> make a plot of the index...
doplot_cpue(a, vartype, mdti, regstr, runreg)
savePlot(file = paste(fname, mdt, vartype, "comp.png", sep = "_"),
         type = "png")
dev.off()

####@> make a new glm to estimate the average index per year...
mod <- glm(pr ~ yr + qtr, data = a)
nd <- data.frame(yr = sort(unique(a$yr[!is.na(a$pr)])),
                 qtr = levels(a$qtr)[2])
nd$pr <- predict(mod, newdat = nd, se.fit = FALSE)

####@> exporting the results...
write.csv(a, file = paste(fname, modtype, ".csv", sep = "_"))
write.csv(nd, file = paste(fname, modtype, "yr.csv", sep = "_"))

######@> Comparisons between indexes by regions and integration with
######@> boats - Delta LogNormal...

######@> Creating a tag for files...
tag <- c("BR_regY1_R02_dellog_boat_allyrs_yr.csv",
         "BR_regY1_R03_dellog_boat_allyrs_yr.csv",
         "BR_Reg_All_Integrated_dellog_boat_allyrs_yr.csv")

######@> Plot lines to compare...
plot(1998:2018, 1998:2018, type = "n", ylim = c(0, 3),
     xlim = c(1998, 2018), xlab = "Years", ylab = "Scaled index",
     axes = FALSE)
axis(1, at = 1998:2018, labels = 1998:2018)
axis(2, at = 0:3, labels = 0:3)
for(i in 1:3) {
    a <- read.csv(paste(tag[i]))
    lines(a$yr, a$pr, col = i, lwd = 2)
}
legend("topright", legend = c("Y1 - R02", "Y1 - R03", "All Regions"),
       col = 1:3, lwd = 2)
savePlot("ALB_DeltaLogN_Boat_Comparison_by_Regions_ver02.png", type = "png")

######@> Comparisons between indexes by regions and integration with
######@> no boats - Delta LogNormal...

######@> Creating a tag for files...
tag <- c("BR_regY1_R02_dellog_novess_allyrs_yr.csv",
         "BR_regY1_R03_dellog_novess_allyrs_yr.csv",
         "BR_Reg_All_Integrated_dellog_novess_allyrs_yr.csv")

######@> Plot lines to compare...
plot(1998:2018, 1998:2018, type = "n", ylim = c(0, 3),
     xlim = c(1998, 2018), xlab = "Years", ylab = "Scaled index",
     axes = FALSE)
axis(1, at = 1998:2018, labels = 1998:2018)
axis(2, at = 0:3, labels = 0:3)
for(i in 1:3) {
    a <- read.csv(paste(tag[i]))
    lines(a$yr, a$pr, col = i, lwd = 2)
}
legend("topright", legend = c("Y1 - R02", "Y1 - R03", "All Regions"),
       col = 1:3, lwd = 2)
savePlot("ALB_DeltaLogN_novess_Comparison_by_Regions.png", type = "png")

######@> Comparisons between indexes by regions and integration with
######@> boats - LognC...

######@> Creating a tag for files...
tag <- c("BR_regY1_R02_lognC_boat_allyrs_yr.csv",
         "BR_regY1_R03_lognC_boat_allyrs_yr.csv",
         "BR_Reg_All_Integrated_lognC_boat_allyrs_yr.csv")

######@> Plot lines to compare...
plot(1998:2018, 1998:2018, type = "n", ylim = c(0, 3),
     xlim = c(1998, 2018), xlab = "Years", ylab = "Scaled index",
     axes = FALSE)
axis(1, at = 1998:2018, labels = 1998:2018)
axis(2, at = 0:3, labels = 0:3)
for(i in 1:3) {
    a <- read.csv(paste(tag[i]))
    lines(a$yr, a$pr, col = i, lwd = 2)
}
legend("topright", legend = c("Y1 - R02", "Y1 - R03", "All Regions"),
       col = 1:3, lwd = 2)
savePlot("ALB_LognC_Boat_Comparison_by_Regions.png", type = "png")

######@> Comparisons between indexes by regions and integration with
######@> no boats - Delta LogNormal...

######@> Creating a tag for files...
tag <- c("BR_regY1_R02_lognC_novess_allyrs_yr.csv",
         "BR_regY1_R03_lognC_novess_allyrs_yr.csv",
         "BR_Reg_All_Integrated_lognC_novess_allyrs_yr.csv")

######@> Plot lines to compare...
plot(1998:2018, 1998:2018, type = "n", ylim = c(0, 3),
     xlim = c(1998, 2018), xlab = "Years", ylab = "Scaled index",
     axes = FALSE)
axis(1, at = 1998:2018, labels = 1998:2018)
axis(2, at = 0:3, labels = 0:3)
for(i in 1:3) {
    a <- read.csv(paste(tag[i]))
    lines(a$yr, a$pr, col = i, lwd = 2)
}
legend("topright", legend = c("Y1 - R02", "Y1 - R03", "All Regions"),
       col = 1:3, lwd = 2)
savePlot("ALB_LognC_novess_Comparison_by_Regions.png", type = "png")

######@> Influence analyses - LognC - Reg Y1 02 - No vessel id...

#####@> Loading the model results...
load("../BR_regY1_R02_lognC_novess_allyrs_model.RData")

#####@> Setup for the analyses...
runsp <- "alb"
mn <- with(mod$data, 0.1 * mean(get(runsp)/hooks))
assign("wtt.all", mk_wts(mod$data, wttype = "area"))
glmdat <- select_data_JointIO(jdat02.cut,
                              runreg = 2,
                              clk = clk_A,
                              minqtrs = 2,
                              runsp = "alb",
                              mt = "deltabin",
                              vars = c("vessid", "hooks", "yrqtr",
                                       "latlong", "hbf"),
                              oneflag = "BR",
                              maxqtrs = 500,
                              minvess = 5,
                              minll = 5,
                              minyrqtr = 50,
                              yrlims = c(1998, 2019),
                              addcl = TRUE,
                              cltype = "hcltrp",
                              addpca = NA,
                              samp = NA,
                              strsmp = NA)

#####@> Initializing the suport analyses to graphics...
infve <- Influence$new(mod)
infve$calc()

####@> Building graphics - standardise plot...
fname <- paste0("BR_Stanplot_LognC", runsp, "_regY1_R02_novess_allyrs")
infve$stanPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - step plot...
fname <- paste0("BR_Step_LognC", runsp, "_regY1_R02_novess_allyrs")
infve$stepPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - Influence plot...
fname <- paste0("BR_Influence_LognC", runsp, "_regY1_R02_novess_allyrs")
infve$influPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - latlong...
fname <- paste0("BR_Cdi_LatLong_LognC", runsp, "_regY1_R02_novess_allyrs")
infve$cdiPlot("latlong")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - Hooks...
fname <- paste0("BR_Cdi_Hooks_LognC", runsp, "_regY1_R02_novess_allyrs")
infve$cdiPlot("ns(hooks, 10)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - hpf...
fname <- paste0("BR_Cdi_Hbf_LognC", runsp, "_regY1_R02_novess_allyrs")
infve$cdiPlot("ns(hbf, 1)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - cluster...
fname <- paste0("BR_Cdi_Cluster_LognC", runsp, "_regY1_R02_novess_allyrs")
infve$cdiPlot("clust")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

#####@> Removing mod object...
rm(mod)

######@> Influence analyses - LognC - Reg Y1 02 - Boat ...

#####@> Loading the model results...
load("../BR_regY1_R02_lognC_boat_allyrs_model.RData")

#####@> Setup for the analyses...
runsp <- "alb"
mn <- with(mod$data, 0.1 * mean(get(runsp)/hooks))
assign("wtt.all", mk_wts(mod$data, wttype = "area"))
glmdat <- select_data_JointIO(jdat02.cut,
                              runreg = 2,
                              clk = clk_A,
                              minqtrs = 2,
                              runsp = "alb",
                              mt = "deltabin",
                              vars = c("vessid", "hooks", "yrqtr",
                                       "latlong", "hbf"),
                              oneflag = "BR",
                              maxqtrs = 500,
                              minvess = 5,
                              minll = 5,
                              minyrqtr = 50,
                              yrlims = c(1998, 2019),
                              addcl = TRUE,
                              cltype = "hcltrp",
                              addpca = NA,
                              samp = NA,
                              strsmp = NA)

#####@> Initializing the suport analyses to graphics...
infve <- Influence$new(mod)
infve$calc()

####@> Building graphics - standardise plot...
fname <- paste0("BR_Stanplot_LognC_", runsp, "_regY1_R02_boat_allyrs")
infve$stanPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - step plot...
fname <- paste0("BR_Step_LognC_", runsp, "_regY1_R02_boat_allyrs")
infve$stepPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - Influence plot...
fname <- paste0("BR_Influence_LognC_", runsp, "_regY1_R02_boat_allyrs")
infve$influPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - latlong...
fname <- paste0("BR_Cdi_LatLong_LognC_", runsp, "_regY1_R02_boat_allyrs")
infve$cdiPlot("latlong")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - Hooks...
fname <- paste0("BR_Cdi_Hooks_LognC_", runsp, "_regY1_R02_boat_allyrs")
infve$cdiPlot("ns(hooks, 10)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - hpf...
fname <- paste0("BR_Cdi_Hbf_LognC_", runsp, "_regY1_R02_boat_allyrs")
infve$cdiPlot("ns(hbf, 1)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - cluster...
fname <- paste0("BR_Cdi_Cluster_LognC_", runsp, "_regY1_R02_boat_allyrs")
infve$cdiPlot("clust")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

#####@> Removing the model file...
rm(mod)

######@> Influence analyses - LognC - Reg Y1 03 - No vessel id...

#####@> Loading the model results...
load("../BR_regY1_R03_lognC_novess_allyrs_model.RData", verbose = TRUE)

#####@> Setup for the analyses...
runsp <- "alb"
mn <- with(mod$data, 0.1 * mean(get(runsp)/hooks))
assign("wtt.all", mk_wts(mod$data, wttype = "area"))
glmdat <- select_data_JointIO(jdat03.cut,
                              runreg = 3,
                              clk = clk_A,
                              minqtrs = 2,
                              runsp = "alb",
                              mt = "deltabin",
                              vars = c("vessid", "hooks", "yrqtr",
                                       "latlong", "hbf"),
                              oneflag = "BR",
                              maxqtrs = 500,
                              minvess = 5,
                              minll = 5,
                              minyrqtr = 50,
                              yrlims = c(1998, 2019),
                              addcl = TRUE,
                              cltype = "hcltrp",
                              addpca = NA,
                              samp = NA,
                              strsmp = NA)

#####@> Initializing the suport analyses to graphics...
infve <- Influence$new(mod)
infve$calc()

####@> Building graphics - standardise plot...
fname <- paste0("BR_Stanplot_LognC_", runsp, "_regY1_R03_novess_allyrs")
infve$stanPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - step plot...
fname <- paste0("BR_Step_LognC_", runsp, "_regY1_R03_novess_allyrs")
infve$stepPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - Influence plot...
fname <- paste0("BR_Influence_LognC_", runsp, "_regY1_R03_novess_allyrs")
infve$influPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - latlong...
fname <- paste0("BR_Cdi_LatLong_LognC_", runsp, "_regY1_R03_novess_allyrs")
infve$cdiPlot("latlong")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - Hooks...
fname <- paste0("BR_Cdi_Hooks_LognC_", runsp, "_regY1_R03_novess_allyrs")
infve$cdiPlot("ns(hooks, 10)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - hpf...
fname <- paste0("BR_Cdi_Hbf_LognC_", runsp, "_regY1_R03_novess_allyrs")
infve$cdiPlot("ns(hbf, 1)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - cluster...
fname <- paste0("BR_Cdi_Cluster_LognC_", runsp, "_regY1_R03_novess_allyrs")
infve$cdiPlot("clust")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

#####@> Removing the model file...
rm(mod)

######@> Influence analyses - LognC - Reg Y1 03 - Boat ...

#####@> Loading the model results...
load("../BR_regY1_R03_lognC_boat_allyrs_model.RData")

#####@> Setup for the analyses...
runsp <- "alb"
mn <- with(mod$data, 0.1 * mean(get(runsp)/hooks))
assign("wtt.all", mk_wts(mod$data, wttype = "area"))
glmdat <- select_data_JointIO(jdat03.cut,
                              runreg = 3,
                              clk = clk_A,
                              minqtrs = 2,
                              runsp = "alb",
                              mt = "deltabin",
                              vars = c("vessid", "hooks", "yrqtr",
                                       "latlong", "hbf"),
                              oneflag = "BR",
                              maxqtrs = 500,
                              minvess = 5,
                              minll = 5,
                              minyrqtr = 50,
                              yrlims = c(1998, 2019),
                              addcl = TRUE,
                              cltype = "hcltrp",
                              addpca = NA,
                              samp = NA,
                              strsmp = NA)

#####@> Initializing the suport analyses to graphics...
infve <- Influence$new(mod)
infve$calc()

####@> Building graphics - standardise plot...
fname <- paste0("BR_Stanplot_LognC_", runsp, "_regY1_R03_boat_allyrs")
infve$stanPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - step plot...
fname <- paste0("BR_Step_LognC_", runsp, "_regY1_R03_boat_allyrs")
infve$stepPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - Influence plot...
fname <- paste0("BR_Influence_LognC_", runsp, "_regY1_R03_boat_allyrs")
infve$influPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - latlong...
fname <- paste0("BR_Cdi_LatLong_LognC_", runsp, "_regY1_R03_boat_allyrs")
infve$cdiPlot("latlong")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - Hooks...
fname <- paste0("BR_Cdi_Hooks_LognC_", runsp, "_regY1_R03_boat_allyrs")
infve$cdiPlot("ns(hooks, 10)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - hpf...
fname <- paste0("BR_Cdi_Hbf_LognC_", runsp, "_regY1_R03_boat_allyrs")
infve$cdiPlot("ns(hbf, 1)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - cluster...
fname <- paste0("BR_Cdi_Cluster_LognC_", runsp, "_regY1_R03_boat_allyrs")
infve$cdiPlot("clust")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

#####@> Removing the model file...
rm(mod)

######@> Influence analyses - LognC - Integrated Regions - No vessel id...

#####@> Loading the model results...
load("../BR_Reg_All_Integrated_lognC_boat_allyrs_model.RData",
     verbose = TRUE)

#####@> Setup for the analyses...
runsp <- "alb"
mn <- with(mod$data, 0.1 * mean(get(runsp)/hooks))
assign("wtt.all", mk_wts(mod$data, wttype = "area"))
glmdat <- select_data_JointIO(jdatAll.cut,
                              runreg = 1,
                              clk = clk_A,
                              minqtrs = 2,
                              runsp = "alb",
                              mt = "deltabin",
                              vars = c("vessid", "hooks", "yrqtr",
                                       "latlong", "hbf"),
                              oneflag = "BR",
                              maxqtrs = 500,
                              minvess = 5,
                              minll = 5,
                              minyrqtr = 50,
                              addcl = TRUE,
                              yrlims = c(1998, 2019),
                              cltype = "hcltrp",
                              addpca = NA,
                              samp = NA,
                              strsmp = NA)

#####@> Initializing the suport analyses to graphics...
infve <- Influence$new(mod)
infve$calc()

####@> Building graphics - standardise plot...
fname <- paste0("BR_Stanplot_LognC_", runsp, "All_Regions_novess_allyrs")
infve$stanPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - step plot...
fname <- paste0("BR_Step_LognC_", runsp, "_All_Regions_novess_allyrs")
infve$stepPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - Influence plot...
fname <- paste0("BR_Influence_LognC_", runsp, "_All_Regions_novess_allyrs")
infve$influPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - latlong...
fname <- paste0("BR_Cdi_LatLong_LognC_", runsp, "_All_Regions_novess_allyrs")
infve$cdiPlot("latlong")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - Hooks...
fname <- paste0("BR_Cdi_Hooks_LognC_", runsp, "_All_Regions_novess_allyrs")
infve$cdiPlot("ns(hooks, 10)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - hpf...
fname <- paste0("BR_Cdi_Hbf_LognC_", runsp, "_All_Regions_novess_allyrs")
infve$cdiPlot("ns(hbf, 1)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - cluster...
fname <- paste0("BR_Cdi_Cluster_LognC_", runsp, "_All_Regions_novess_allyrs")
infve$cdiPlot("clust")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

#####@> Removing the model file...
rm(mod)

######@> Influence analyses - LognC - All Regions - Boat ...

#####@> Loading the model results...
load("../BR_Reg_All_Integrated_lognC_boat_allyrs_model.RData")

#####@> Setup for the analyses...
runsp <- "alb"
mn <- with(mod$data, 0.1 * mean(get(runsp)/hooks))
assign("wtt.all", mk_wts(mod$data, wttype = "area"))
glmdat <- select_data_JointIO(jdatAll.cut,
                              runreg = 1,
                              clk = clk_A,
                              minqtrs = 2,
                              runsp = "alb",
                              mt = "deltabin",
                              vars = c("vessid", "hooks", "yrqtr",
                                       "latlong", "hbf"),
                              oneflag = "BR",
                              maxqtrs = 500,
                              minvess = 5,
                              minll = 5,
                              minyrqtr = 50,
                              yrlims = c(1998, 2019),
                              addcl = TRUE,
                              cltype = "hcltrp",
                              addpca = NA,
                              samp = NA,
                              strsmp = NA)

#####@> Initializing the suport analyses to graphics...
infve <- Influence$new(mod)
infve$calc()

####@> Building graphics - standardise plot...
fname <- paste0("BR_Stanplot_LognC_", runsp, "_All_Regions_boat_allyrs")
infve$stanPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - step plot...
fname <- paste0("BR_Step_LognC_", runsp, "_All_Regions_boat_allyrs")
infve$stepPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - Influence plot...
fname <- paste0("BR_Influence_LognC_", runsp, "_All_Regions_boat_allyrs")
infve$influPlot()
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - latlong...
fname <- paste0("BR_Cdi_LatLong_LognC_", runsp, "_All_Regions_boat_allyrs")
infve$cdiPlot("latlong")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - Hooks...
fname <- paste0("BR_Cdi_Hooks_LognC_", runsp, "_All_Regions_boat_allyrs")
infve$cdiPlot("ns(hooks, 10)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - hpf...
fname <- paste0("BR_Cdi_Hbf_LognC_", runsp, "_All_Regions_boat_allyrs")
infve$cdiPlot("ns(hbf, 1)")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

####@> Building graphics - cdi plot - cluster...
fname <- paste0("BR_Cdi_Cluster_LognC_", runsp, "_All_Regions_boat_allyrs")
infve$cdiPlot("clust")
savePlot(paste0(fname, ".png"), type = "png")
dev.off()

#####@> Removing the model file...
rm(mod)

#######@> Mapping clusters per decade - Region 02...

######@> Creating a new column - Decade...
a <- jdat02.cut
a$decade <- with(a, op_yr - op_yr %% 10)

######@> Aggregating data per decade...
tmp <- a %>%
    group_by(lon5, lat5, decade, hcltrp) %>%
    summarise(alb = sum(alb, na.rm = TRUE)) %>%
    as.data.frame()

#####@> Standardizing scale...
quantile(tmp$alb, probs = seq(0, 1, 0.01))

int.v <- round(quantile(tmp$alb, probs = seq(0, 1, l = 11)))
int.v[2] <- 1
int.c <- int.v[1:(length(int.v)-1)]
int.c2 <- int.v[2:length(int.v)]
tmp$disc <- cut(tmp$alb, breaks = int.v,
                include.lowest = TRUE, right = FALSE,
                labels = int.c, dig.lab = 0)

####@> Choosen the color pallete...
library(RColorBrewer)
gPal <- colorRampPalette(rev(brewer.pal(10, "RdYlGn")))
cor <- gPal(length(int.c))

#####@> mapping...
p00 <- ggplot(data = tmp, aes(x = lon5, y = lat5)) +
    geom_tile(aes(fill = disc), colour = "black") +
    geom_polygon(data = mm, aes(x = long, y = lat, group = group)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    facet_grid(hcltrp ~ decade) +
    ## scale_fill_manual(name = "Areas",
    ##                   values = c("red", "green", "blue")) +
    scale_fill_manual("Catches",
                      values = setNames(cor, int.c),
                      labels = rev(paste("[", int.c, ", ",
                                         int.c2, ")", sep = "")),
                      limits = int.c, breaks = rev(int.c)) +
    coord_fixed(xlim = c(-55, -15), ylim = c(-15, 10)) +
    xlab(expression(paste("Longitude ", "(", degree, ")"))) +
    ylab(expression(paste("Latitude ", "(", degree, ")"))) +
    my_theme()
p00

savePlot("Map_Clusters_Decade_Y1_R02.png", type = "png")

#######@> Mapping clusters per decade - Region 02...

######@> Creating a new column - Decade...
a <- jdat03.cut
a$decade <- with(a, op_yr - op_yr %% 10)

######@> Aggregating data per decade...
tmp <- a %>%
    group_by(lon5, lat5, decade, hcltrp) %>%
    summarise(alb = sum(alb, na.rm = TRUE)) %>%
    as.data.frame()

#####@> padronizando a escala...
quantile(tmp$alb, probs = seq(0, 1, l = 11))

int.v <- round(quantile(tmp$alb, probs = seq(0, 1, l = 11)))
int.v[2] <- 1
int.c <- int.v[1:(length(int.v)-1)]
int.c2 <- int.v[2:length(int.v)]
tmp$disc <- cut(tmp$alb, breaks = int.v,
                include.lowest = TRUE, right = FALSE,
                labels = int.c, dig.lab = 0)

####@> Escolher a paleta de cores...
gPal <- colorRampPalette(rev(brewer.pal(10, "RdYlGn")))
cor <- gPal(length(int.c))

#####@> mapping...
p00 <- ggplot(data = tmp, aes(x = lon5, y = lat5)) +
    geom_tile(aes(fill = disc), colour = "black") +
    geom_polygon(data = mm, aes(x = long, y = lat, group = group)) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    facet_grid(hcltrp ~ decade) +
    ## scale_fill_manual(name = "Areas",
    ##                   values = c("red", "green", "blue")) +
    scale_fill_manual("Catches",
                      values = setNames(cor, int.c),
                      labels = rev(paste("[", int.c, ", ",
                                         int.c2, ")", sep = "")),
                      limits = int.c, breaks = rev(int.c)) +
    coord_fixed(xlim = c(-55, -15), ylim = c(-40, -15)) +
    xlab(expression(paste("Longitude ", "(", degree, ")"))) +
    ylab(expression(paste("Latitude ", "(", degree, ")"))) +
    my_theme()
p00
savePlot("Map_Clusters_Decade_Y1_R03.png", type = "png")

## Vizualizando as producoes por grupos... Para isso é necessario transformar o
## data.frame e remover ano e cod...
temp <- dplyr::select(jdat02, hcltrp, yft, alb, bet, swo, sai, bum, bsh,
                      whm, sma)
temp.melt <- melt(temp, id.vars = "hcltrp")
col <- c("#000000", "#252525", rev(brewer.pal(9, "Purples")))

temp <- temp.melt %>%
    group_by(hcltrp, variable) %>%
    summarise(value = sum(value)) %>%
    mutate(P = value/sum(value))

p <- ggplot(temp, aes(x = factor(1), y = P, fill = factor(variable))) +
    geom_bar(stat = "identity") +
    coord_polar(theta = "y", direction = 2) +
    facet_wrap(~ hcltrp, ncol = 2) + xlab("") + ylab("") +
    labs(fill = "Species") +
    scale_fill_manual(values = col)
p

png("Pie-chart_species_by_group.png", width = 1080, height = 1080,
    res = 100)
print(p)
dev.off()

########################################################################
##
##                 Creative Commons License 4.0
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
