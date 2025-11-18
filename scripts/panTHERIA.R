#trying to get the pantheria database to work

# Install pak
install.packages("pak")
library(pak)

#with the pkg_install function from pak install the traitdata database from github
pkg_install("RS-eco/traitdata")

#load library
library(traitdata)
test <-data(pantheria)



# Section 2:  phylacine database ------------------------------------------

#tutorial from github https://github.com/MegaPast2Future/PHYLACINE_1.2#vignette 

#install.packages("pacman", repos="https://cloud.r-project.org")
# pacman::p_load(ggplot2,
#                dplyr,
#                stringr,
#                gridExtra,
#                viridisLite,
#                raster,
#                rasterVis,
#                rgdal,
#                maptools,
#                ape,
#                ggtree, update = F)

#install.packages("maptools")
# Source - https://stackoverflow.com/a
# Posted by margusl
# Retrieved 2025-11-11, License - CC BY-SA 4.0

install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")

library(maptools)
library(raster)
library(rasterVis)
library(rgdal)
library(virdisLite)


forest <- read.nexus("Data/Phylogenies/Complete_phylogeny.nex")
names(forest) <- NULL
set.seed(42)
forest <- forest[sample(1:1000, 30)]

# Load world map and subset to Australia
data(wrld_simpl)
australia <- wrld_simpl[wrld_simpl$NAME == "Australia", ]

# Load trait data. Remember to always use UTF-8 encoding with PHYLACINE files to avoid weird characters 
mam <- read.csv("Data/Traits/Trait_data.csv", fileEncoding = "UTF-8", stringsAsFactors = F)

# Set factor levels for IUCN status. "EP" is a new status we added to designate species that went extinct in prehistory like Diprotodon  
mam$IUCN.Status.1.2 <- factor(mam$IUCN.Status.1.2, levels=c("EP", "EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD"))

# Subset to species that are in Australian marsupial orders, over 20 kg, and over 90 % herbivorous
marsupial.orders <- c("Dasyuromorphia", "Peramelemorphia",
                      "Notoryctemorphia", "Diprotodontia")
marsupials <- mam[mam$Order.1.2 %in% marsupial.orders, ]
marsupials <- marsupials[marsupials$Mass.g > 20000, ]
marsupials <- marsupials[marsupials$Diet.Plant >= 90, ]

#range data example
maps.current <- paste0("C://Users//ameli//OneDrive//Documents//R_projects//fish_sleep//Data//Ranges//Current//", marsupials$Binomial.1.2, ".tif")
r.current <- stack(maps.current)
maps.pres.nat <- paste0("Data/Ranges/Present_natural/", marsupials$Binomial.1.2, ".tif")
r.pres.nat <- stack(maps.pres.nat)

#r.current <- "C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\fish_sleep\\Data\\Ranges\\Current\\Abditomys_latidens.tif"

australia <- spTransform(australia, crs(r.current))

# Crop range maps to just the extent of Australia for a cleaner plot
ext <- extent(australia)
ext[2] <- 15000000 # Reduce eastern extent
ext[3] <- -5200000 # Reduce southern extent
ext[4] <- -1370000 # Reduce northern extent
r.current <- crop(r.current, ext)
r.pres.nat <- crop(r.pres.nat, ext)

# Create a blank raster of the region
blank.r <- r.current[[1]]
blank.r[] <- NA
names(blank.r) <- NA

m.current <- matrix(NA, nrow=nrow(marsupials), ncol=dim(r.current)[1]*dim(r.current)[2])
rownames(m.current) <- marsupials$Binomial.1.2
for(i in 1:nrow(marsupials)) {
  m.current[i, ] <- getValues(r.current[[i]])
}
# Current species list
sp.current <- marsupials$Binomial.1.2[rowSums(m.current) > 0]

# Load all the present natural raster data into a matrix for faster handling
m.pres.nat <- matrix(NA, nrow=nrow(marsupials), ncol=dim(r.current)[1]*dim(r.current)[2])
rownames(m.pres.nat) <- marsupials$Binomial.1.2
for(i in 1:nrow(marsupials)) {
  m.pres.nat[i, ] <- getValues(r.pres.nat[[i]])
}
# Present natural species list
sp.pres.nat <- marsupials$Binomial.1.2[rowSums(m.pres.nat) > 0]

# Create rasters of taxonomic richness
current.div <- blank.r
current.div[] <- colSums(m.current)
names(current.div) <- "Current diversity"
pres.nat.div <- blank.r
pres.nat.div[] <- colSums(m.pres.nat)
names(pres.nat.div) <- "Present natural diversity"
div <- stack(current.div, pres.nat.div)
# Change 0 to NA, for nicer plotting
div[] <- ifelse(div[] == 0, NA, div[])

p.map <- 
  levelplot(div,
            layout = c(1,2),
            colorkey=list(
              space='left',                   
              labels=list(at = 1:max(div[], na.rm = T), font=4),
              axis.line=list(col = 'black'),
              width=0.75
            ),
            par.settings = list(
              strip.border = list(col='transparent'),
              strip.background = list(col='transparent'),
              axis.line = list(col = 'transparent')
            ),
            scales = list(draw = FALSE),            
            col.regions = viridis,                   
            at = 1:max(div[], na.rm = T),
            names.attr = str_replace_all(names(div), "\\.", " "))

# Add Australia polygon outline
p.map <- p.map + layer(sp.polygons(australia, col = "darkgrey"))


# Section 2: artiodactyla trait data --------------------------------------

#try for artiodactyla
mam <- read.csv("Data/Traits/Trait_data.csv", fileEncoding = "UTF-8", stringsAsFactors = F)
artio <- filter(mam, Order.1.2 == "Cetartiodactyla")

#contains data on terrestrial species and on cetaceans, 392 total sps
#info on taxonomy, enviro (terrestrial, marine, freshwater), mass, island endemicity, IUCN status, diet (plant, vertebrate, invertebrate)
artio <- artio[, c(1:8, 12, 17, 18, 20:22)]

#range data also exists in separate file

ggplot(artio, aes(x = Family.1.2, y = Mass.g)) + geom_point()

artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"))
colnames(artio) <- c("tips", "Order", "Family", "Genus", "Species", "Terrestrial" ,"Marine", "Freshwater", "Mass.g", "Island.Endemicity", "IUCN.Status", "Diet.Plant", "Diet.Vertebrate", "Diet.Invertebrate")
artio <- merge(artio, artio_full, by = "tips")


# Section 3: biogeography code from Bennie et al --------------------------

## R script to calculate duration of biologically-useful daylight, twilight and moonlight
## Jon Bennie May 2012
## using sections of code written by Tom Davies 2011
## Astronomical formulae derived from Meeus, J. (1988) Astronomical Formulae for Calculators (4th ed.) (Wilmann-Bell, Richmond, Virginia)

## load libraries

require(rgdal)

## set working directory for climate files to 'pathname'
setwd(pathname)

## read sample map from file in folder D:/CRU_TS_3 and extract latitude of grid cells

map <- readGDAL("cru_ts_3_00.1901.2006.tmp_2005_1.asc")
latitude <- coordinates(map)[,2]

## read in climate data

temperature <- read.table("CRU_ST_1961_1990_monthly averages.txt",header=TRUE)
maxtemp <- read.table("CRU_ST_1961_1990_max_monthly averages.txt",header=TRUE)
mintemp <- read.table("CRU_ST_1961_1990_min_monthly averages.txt",header=TRUE)
cloud <- read.table("CRU_ST_1961_1990_cloud_monthlyaverages.txt",header=TRUE)

## set number of days in each month
daysin <- array(c(31,28,31,30,31,30,31,31,30,31,30,31))

## set up arrays for output

hours_light <- array(0,c(12,259200))
hours_dark <- array(0,c(12,259200))
hours_twilight <- array(0,c(12,259200))
hours_moonlight <- array(0,c(12,259200))

## loop through months/days and calculate solar and lunar geometry

julianday <- 0
for (thismonth in 1:12)
{
  lastmonth <- thismonth-1
  if(lastmonth==0) lastmonth <- 12
  nextmonth <- thismonth+1
  if(nextmonth==13) nextmonth <- 1
  for (day in 1:daysin[thismonth])
  {
    ## solar equations to calculate duration of day, night and twilight
    
    julianday <- julianday+1
    declin <- (pi * 23.5 / 180) * cos(2 * pi * ((julianday - 171) / 365.25))
    angle <- -tan(pi*(latitude/180))*tan(declin)
    twilightangle <- angle - (pi * 12/180)/cos(pi*latitude/180)
    angle[angle>=1] <- 1
    angle[angle<=-1] <- -1
    twilightangle[twilightangle>=1] <- 1
    twilightangle[twilightangle<=-1] <- -1
    daylength <- acos(angle)/(pi*(7.5/180))
    nightlength <- 24-daylength
    twilightlength <- acos(twilightangle)/(pi*(7.5/180))-daylength
    
    ## interpolate daily mean, min and max temperatures from monthly means
    if(day<=daysin[thismonth]/2)
    {
      dailytemp <- temperature[,thismonth]-((daysin[thismonth]/2-day)/(daysin[thismonth]/2+daysin[lastmonth]/2))*(temperature[,thismonth]-temperature[,lastmonth])
      dailymax <- maxtemp[,thismonth]-((daysin[thismonth]/2-day)/(daysin[thismonth]/2+daysin[lastmonth]/2))*(maxtemp[,thismonth]-maxtemp[,lastmonth])
      dailymin <- mintemp[,thismonth]-((daysin[thismonth]/2-day)/(daysin[thismonth]/2+daysin[lastmonth]/2))*(mintemp[,thismonth]-mintemp[,lastmonth])
    } else {
      dailytemp <- temperature[,thismonth]+((day-daysin[thismonth]/2)/(daysin[thismonth]/2+daysin[nextmonth]/2))*(temperature[,nextmonth]-temperature[,thismonth])
      dailymax <- maxtemp[,thismonth]+((day-daysin[thismonth]/2)/(daysin[thismonth]/2+daysin[nextmonth]/2))*(maxtemp[,nextmonth]-maxtemp[,thismonth])
      dailymin <- mintemp[,thismonth]+((day-daysin[thismonth]/2)/(daysin[thismonth]/2+daysin[nextmonth]/2))*(mintemp[,nextmonth]-mintemp[,thismonth])
    }	
    
    ## adjust length of day, night and twilight for temperature limits
    ## days are not counted if 24 hr mean temp <= 0 deg C or max temp >= 35 deg C
    ## nights are not counted if min temp <= 0 deg C or mean temp >= 35 deg C
    ## morning twilight not counted if min temp <= 0 deg C or >= 35 deg C
    ## evening twilight not counted if mean temp <= 0 deg C or >= 35 deg C
    
    daylength[dailytemp<=0] <- 0
    daylength[dailymax>350] <- 0
    nightlength[dailymin<=0] <- 0
    nightlength[dailytemp>350] <- 0
    dawnmin <- twilightlength[which(dailymin<=0)]
    twilightlength[which(dailymin<=0)] <- 0.5*dawnmin
    twilightlength[dailytemp<=0] <- 0
    duskmax <- twilightlength[which(dailytemp>=350)]
    twilightlength[which(dailytemp>=350)] <- 0.5*duskmax
    twilightlength[dailymin>=350] <- 0
    hours_light[thismonth,] <- hours_light[thismonth,]+daylength
    hours_dark[thismonth,] <- hours_dark[thismonth,]+nightlength
    hours_twilight[thismonth,] <- hours_twilight[thismonth,]+twilightlength
    
    ## lunar calculations for moonlight
    
    d <- julianday+0.5
    N = as.double(125.1228-0.0529538083*d) %% 360
    i = as.double(5.1454)
    w = as.double(318.0634+0.1643573223*d) %%360
    a = as.double(60.2666)
    e = as.double(0.054900)
    M = as.double(115.3654+13.0649929509*d) %%360
    ecl = as.double(23.4393-3.563E-7*d)
    
    ## calculate phase of moon, illuminated fraction
    
    phase_angle <- 2*pi*((d-6)%%29.530588853)/29.530588853
    illuminated_fraction <- (1+cos(phase_angle))/2
    
    E = M + e*(180/pi) * sin(pi*M/180) * ( 1.0 + e * cos(pi*M/180) )
    
    xv = a * ( cos(pi*E/180) - e )
    yv = a * ( sqrt(1.0 - e*e) * sin(pi*E/180) )
    
    v = (atan2(yv, xv )*180/pi) %% 360
    r = sqrt( xv*xv + yv*yv )
    
    xh = r * ( cos(pi*N/180) * cos(pi*(v+w)/180) - sin(pi*N/180) * sin(pi*(v+w)/180) * cos(pi*i/180) )
    yh = r * ( sin(pi*N/180) * cos(pi*(v+w)/180) + cos(pi*N/180) * sin(pi*(v+w)/180) * cos(pi*i/180) )
    zh = r * ( sin(pi*(v+w)/180) * sin(pi*i/180) )
    
    xe = xh
    ye = yh * cos(pi*ecl/180)- zh * sin(pi*ecl/180)
    ze = yh * sin(pi*ecl/180)- zh * cos(pi*ecl/180)
    
    RA = atan2(ye,xe)
    declin = atan2(ze,sqrt(xe^2+ye^2))
    
    coslha <- -(sin(latitude*pi/180)*sin(declin))/(cos(latitude*pi/180)*cos(declin))
    coslha[coslha<=-1] <- -1
    coslha[coslha>=1] <- 1
    moonrise_sid <- RA-acos(coslha)
    moonset_sid <- RA+acos(coslha)
    
    w_sun = as.double(282.9404+4.70935E-5*d) %%360
    e_sun = as.double(0.016709-1.151E-9*d) %%360
    M_sun = as.double(356.0470+0.9856002585*d) %%360
    
    E_sun = M_sun + e_sun*(180/pi) * sin(pi*M_sun/180) * ( 1.0 + e_sun * cos(pi*M_sun/180) )
    
    xv_sun = cos(pi*E_sun/180) - e_sun
    yv_sun = sqrt(1.0 - e_sun^2) * sin(pi*E_sun/180)
    
    v_sun = (atan2(yv_sun, xv_sun )*180/pi) %% 360
    r_sun = sqrt( xv_sun^2 + yv_sun^2 )
    
    xe_sun = r_sun * cos((v_sun+w_sun)*pi/180)
    ye_sun = r_sun * sin((v_sun+w_sun)*pi/180) * cos(pi*ecl/180)
    ze_sun = r_sun * sin((v_sun+w_sun)*pi/180) * sin(pi*ecl/180)
    
    RA_sun  = atan2( ye_sun, xe_sun )
    declin_sun = atan2(ze_sun,sqrt(xe_sun^2+ye_sun^2))
    coslha_sun <- (sin(12*pi/180)-sin(latitude*pi/180)*sin(declin_sun))/(cos(latitude*pi/180)*cos(declin_sun))
    coslha_sun[coslha_sun<=-1] <- -1
    coslha_sun[coslha_sun>=1] <- 1
    
    ##times after sunset
    
    sunrise_sol <- 2*pi-2*acos(coslha_sun)
    sunset_sol <- 0*acos(coslha_sun)
    moonrise_sol <- moonrise_sid-RA_sun+acos(coslha_sun)
    moonset_sol <- moonset_sid-RA_sun+acos(coslha_sun)
    moonrise_sol <- (moonrise_sol)%%(2*pi)
    moonset_sol <- (moonset_sol)%%(2*pi)
    
    moonlightlength <- sunset_sol
    moonlightlength[moonrise_sol <= moonset_sol] <- pmin(sunrise_sol[moonrise_sol <= moonset_sol],moonset_sol[moonrise_sol <= moonset_sol])-pmin(sunrise_sol[moonrise_sol <= moonset_sol],moonrise_sol[moonrise_sol <= moonset_sol])
    moonlightlength[moonrise_sol > moonset_sol] <- pmin(sunrise_sol[moonrise_sol > moonset_sol],moonset_sol[moonrise_sol > moonset_sol]) + sunrise_sol[moonrise_sol > moonset_sol] - pmin(sunrise_sol[moonrise_sol > moonset_sol],moonrise_sol[moonrise_sol > moonset_sol])
    moonlightlength[illuminated_fraction<0.75] <- 0
    moonlightlength <- moonlightlength*24/(2*pi)
    moonlightlength <- moonlightlength*(1000-cloud[,thismonth])/1000
    moonlightlength[dailymin<=0] <- 0
    moonlightlength[dailymin>=350] <- 0
    hours_moonlight[thismonth,] <- hours_moonlight[thismonth,]+moonlightlength
  }
}

hours_light[is.na(hours_moonlight)] <- NA
hours_twilight[is.na(hours_moonlight)] <- NA

sum_moonlight <- apply(hours_moonlight,2,sum)
sum_twilight <- apply(hours_twilight,2,sum)
sum_daylight <- apply(hours_light,2,sum)
sum_darkness <- apply(hours_dark-(hours_moonlight+hours_twilight),2,sum)

ratio_moonlight <- sum_moonlight/(sum_moonlight+sum_twilight+sum_daylight+sum_darkness)
ratio_twilight <- sum_twilight/(sum_twilight+sum_moonlight+sum_daylight+sum_darkness)
ratio_daylight <- sum_daylight/(sum_twilight+sum_moonlight+sum_daylight+sum_darkness)

var_ratio_moon <- apply(hours_moonlight/(hours_twilight+hours_light+hours_moonlight+hours_dark),2,var,na.rm=T)
var_ratio_twi <- apply(hours_twilight/(hours_moonlight+hours_light+hours_twilight+hours_dark),2,var,na.rm=T)
var_ratio_day <- apply(hours_light/(hours_moonlight+hours_twilight+hours_light+hours_dark),2,var,na.rm=T)

## output maps of moonlight/twilight/daylight etc..

map$band1 <- ratio_moonlight
##spplot(map)
writeGDAL(map,"ratio_moonlight.tif",mvFlag=-999)

map$band1 <- ratio_twilight
##spplot(map)
writeGDAL(map,"ratio_twilight.tif",mvFlag=-999)

map$band1 <- ratio_daylight
##spplot(map)
writeGDAL(map,"ratio_daylight.tif",mvFlag=-999)

map$band1 <- var_ratio_moon
##spplot(map)
writeGDAL(map,"var_ratio_moon.tif",mvFlag=-999)

map$band1 <- var_ratio_twi
##spplot(map)
writeGDAL(map,"var_ratio_twi.tif",mvFlag=-999)

map$band1 <- var_ratio_day
##spplot(map)
writeGDAL(map,"var_ratio_day.tif",mvFlag=-999)
