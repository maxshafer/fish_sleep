library(stringr)
library(lubridate)
library(data.table)
library(rfishbase)
library(tidyr)
library(ggbiplot)
library(dplyr)

species1 <- read.csv("~/Downloads/GBIF_iNat_Amblyglyphidodon_curacao.csv", sep = "\t")
species2 <- read.csv("~/Downloads/GBIF_iNat_Porichthys_notatus.csv", sep = "\t")
species3 <- read.csv("~/Downloads/GBIF_iNat_Chaetodon_auriga.csv", sep = "\t")
species4 <- read.csv("~/Downloads/GBIF_iNat_Lepomis_gibbosus.csv", sep = "\t")

species <- rbind(species1, species2, species3, species4)

species$datetime <- str_replace(species$dateIdentified, "T", " ")

species <- species[species$datetime != "",]

head(as.POSIXct(species$datetime, format = "%Y-%m-%d %H:%M:%S"))

species$datetime <- as.POSIXct(species$datetime, format = "%Y-%m-%d %H:%M:%S")

day(species$datetime) <- 01
month(species$datetime) <- 01
year(species$datetime) <- 2024

ggplot(species, aes(x = datetime, fill = species)) + geom_density(alpha = 0.25)


## Load in partial full database?
## Full thing is 50Gb, but I only need 3 columns (unique ID, species, datetime)

## Want to keep 1,10,42
full_db <- fread("~/Downloads/GBIF_iNat_fulldatabase_2024-08-06.csv", drop = c(2:9,11:29,31:50))

fishbase_df <- load_taxa(collect = T, version = "21.06")
fishbase_df <- as.data.frame(fishbase_df)

## Take only fish
full_db <- full_db[full_db$species %in% fishbase_df$Species,]

full_db <- full_db[grepl("T", full_db$eventDate),]

iNat_fish <- table(full_db$species)

good_species <- names(iNat_fish)[iNat_fish > 50]

## This ends up being ~11k species

table(sleepy_fish$Species_name %in% good_species)
table(good_species %in% sleepy_fish$Species_name)

fish_db <- full_db[full_db$species %in% good_species,]

fish_db$datetime <- str_replace(fish_db$eventDate, "T", " ")

fish_db <- fish_db[fish_db$datetime != "",]
## pad out with seconds, if missing
index <- nchar(fish_db$datetime) != 19
fish_db$datetime[index] <- paste(fish_db$datetime[index], ":00", sep = "")

fish_db$datetime <- as.POSIXct(fish_db$datetime, format = "%Y-%m-%d %H:%M:%S")

day(fish_db$datetime) <- 01
month(fish_db$datetime) <- 01
year(fish_db$datetime) <- 2024


# ggplot(fish_db[fish_db$species == "Mola mola",], aes(x = datetime, fill = species)) + stat_bin()

# ggplot(fish_db[fish_db$species %in% unique(fish_db$species)[11:20],], aes(x = datetime, fill = species)) + geom_density(alpha = 0.25)

# ggplot(fish_db[fish_db$species %in% unique(fish_db$species)[11:20],], aes(x = datetime, fill = species)) + stat_bin() + facet_wrap(~species, scales = "free") + theme(legend.position = "none")




### seems like this works, and could be useful, but has the same interpretation problems as catch data. E.g. Mola mola looks nocturnal by this, but that's because during the day it is diving up and down hunting for food
### I would have to run some kind of stats, to determine if day vs night vs crepuscular periods have more or less sightings
### I could just ask if sightings are explained by time, and if the p-value is past cutoff, then determine which periods are higher?
### Or just do t-test between day and night (and dawn/dusk??)

### Alternative, I could do some kind of analysis like Annika did? With PCs?
### I'm sure I would find similar results, since the data is the same, but unclear how to make binary calls from this?

## OK, I also need to normalize for the # of observations per hour (which is very skewed)
# fish_db$fill <- 1
# ggplot(fish_db, aes(x = datetime, fill = fill)) + stat_bin()

## First convert to matrix?
## I think maybe I first need to summarize by half hour bin, then convert this (which will have a sum column)

fish_data <- fish_db %>% mutate(half_hour = floor_date(datetime, "30 minutes"))

summarised <- fish_data %>% count(species, half_hour)

summarised <- pivot_wider(summarised, names_from = species, values_from = n)
summarised <- summarised[!(is.na(summarised$half_hour)),]

row_names <- summarised$half_hour
summarised <- summarised[,-1]
summarised <- (summarised/rowSums(summarised,na.rm = T))*10000
rownames(summarised) <- row_names
summarised[is.na(summarised)] <- 0

## The above seems to work, and so does the below:
## PC2 is day vs night, and PC1 seems to be first half of day vs second (I suppose when people are awake and cataloging)

diel.pca <- prcomp(summarised, center = TRUE, scale. = TRUE)
biplot <- ggbiplot(diel.pca, choices = c(1,2), labels = row_names, var.axes = F) + theme_classic() 



### Can I convert back to long format? To plot and compare to before normalization?
summarised$time <- rownames(summarised)
summarised2 <- summarised %>% pivot_longer(!time, names_to = "species", values_to = "n")

summarised2$time <- as.POSIXct(summarised2$time, format = "%Y-%m-%d %H:%M:%S")

ggplot(summarised2[summarised2$species == "Acipenser fulvescens",], aes(x = time, fill = species, y = n)) + geom_point()

ggplot(summarised2[summarised2$species %in% unique(summarised2$species)[21:30],], aes(x = time, colour = species, y = n)) + geom_point() + facet_wrap(~species, scales = "free") + theme(legend.position = "none")

## Add column for day vs night vs dawn vs dusk (if only one of these is higher, it's crepuscular)
summarised2$diel_period <- NA
summarised2$diel_period <- ifelse(between(summarised2$time, as.POSIXct("2024-01-01 00:00:00 EST"), as.POSIXct("2024-01-01 06:00:00 EST")), "night", summarised2$diel_period)
summarised2$diel_period <- ifelse(between(summarised2$time, as.POSIXct("2024-01-01 06:00:00 EST"), as.POSIXct("2024-01-01 08:00:00 EST")), "dawn", summarised2$diel_period)
summarised2$diel_period <- ifelse(between(summarised2$time, as.POSIXct("2024-01-01 08:00:00 EST"), as.POSIXct("2024-01-01 18:00:00 EST")), "day", summarised2$diel_period)
summarised2$diel_period <- ifelse(between(summarised2$time, as.POSIXct("2024-01-01 18:00:00 EST"), as.POSIXct("2024-01-01 20:00:00 EST")), "dusk", summarised2$diel_period)
summarised2$diel_period <- ifelse(between(summarised2$time, as.POSIXct("2024-01-01 20:00:00 EST"), as.POSIXct("2024-01-01 24:00:00 EST")), "night", summarised2$diel_period)

ggplot(summarised2[summarised2$species == "Lepomis gibbosus",], aes(x = diel_period, fill = species, y = n)) + geom_boxplot()
ggplot(summarised2[summarised2$species %in% unique(summarised2$species)[21:30],], aes(x = diel_period, colour = species, y = n)) + geom_boxplot() + facet_wrap(~species, scales = "free") + theme(legend.position = "none")


fit = aov(n ~ diel_period, data = summarised2[summarised2$species == "Lepomis gibbosus",])
summary(fit)

library(ggfortify)
autoplot(fit)
