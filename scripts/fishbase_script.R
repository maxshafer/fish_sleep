library(RISmed)
library(rfishbase)


url <- 'https://docs.google.com/spreadsheets/d/18aNqHT73hX06cGRlf6oj7Y4TVKf6jd_Q5ojNIxm2rys/edit?usp=sharing'
sleepy_fish <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

species <- sleepy_fish$Species.name

## Load fishbase data

fishbase_species <- rfishbase::load_taxa(version = "19.04")

fishbase_diet <- diet(version = "19.04")
fishbase_ecology <- ecology(version = "19.04")
fishbase_morph <- morphometrics(version = "19.04")
fishbase_ecosystem <- ecosystem(version = "19.04")

fishbase_species_sub <- fishbase_species[fishbase_species$Species %in% species,]



species_sub <- as.data.frame(fishbase_species_sub)
species_sub$diel <- resolved_names$diel[match(species_sub$Species, resolved_names$unique_name)]

species_sub$DietTroph <- fishbase_ecology$DietTroph[match(species_sub$Species, fishbase_ecology$Species)]
species_sub$FeedingType <- fishbase_ecology$FeedingType[match(species_sub$Species, fishbase_ecology$Species)]


ggplot(species_sub, aes(fill = diel, x = FeedingType), alpha = 0.25) + geom_density(position = 'identity',bins = 20) + facet_wrap(~diel)





