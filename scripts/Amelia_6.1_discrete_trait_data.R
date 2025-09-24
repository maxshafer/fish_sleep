# Section 1: Churchill et al dive data, mass, length---------------------------------------------
church_pt2 <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt2.xlsx")
colnames(church_pt2) <- c("Species", "Dive_depth", "Log_dive_depth", "Dive_duration", "Log_dive_duration", "Mass", "Log_mass", "Total_body_length", "Log_body_length", "Assymetry_index", "Peak_frequency", "Log_frequency")

church_pt2$tips <- str_replace(church_pt2$Species, " ", "_")

#filter for species with dive data, leaves 62 species
church_pt2 <- church_pt2[!(is.na(church_pt2$Dive_depth)), c("Dive_depth", "Dive_duration", "Mass", "Total_body_length", "tips")]

#correct species spelling to match
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Cephalorhynchus_commersoni", replacement = "Cephalorhynchus_commersonii")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Berardius_bairdi", replacement = "Berardius_bairdii")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Sagmatias_australis", replacement = "Lagenorhynchus_australis")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Sagmatias_cruciger", replacement = "Lagenorhynchus_cruciger")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Sagmatias_obliquidens", replacement = "Lagenorhynchus_obliquidens")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Cephalorhynchus_heavisidi", replacement = "Cephalorhynchus_heavisidii")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Leucopleurus_acutus", replacement = "Lagenorhynchus_acutus")

#read in the additional dive data I collected
url <- 'https://docs.google.com/spreadsheets/d/1_0ZS_tbddOCckkcKn9H5HpVRDZty4jhkUU20Nc0YYQY/edit?usp=sharing'
dive_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

#average dive depth is very inconsistent (deep vs shallow dives, day vs night dives)
#and measured inconsistently (average of means across individuals, average of medians etc) so exclude for now
dive_full <- dive_full[!is.na(dive_full$Diel_Pattern),1:9]

#filter to rows with an alternative value
dive_full_alt <- dive_full[!is.na(dive_full$Alt_Max_1), ]
dive_full <- dive_full[is.na(dive_full$Alt_Max_1), ]

#pick the bigger value between max dive depth and alt dive depth
dive_full_alt[is.na(dive_full_alt)] <- 0
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_1)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_2)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_3)
dive_full_alt$Max_dive_depth_m <- pmax(dive_full_alt$Max_dive_depth_m, dive_full_alt$Alt_Max_4)

#combine with original values
dive_full <- rbind(dive_full, dive_full_alt)
dive_full$tips <- str_replace(dive_full$Species, " ", "_")

dive_full <- dive_full[, c("Max_dive_depth_m", "tips")]

church_pt2 <- merge(church_pt2, dive_full, all.x = TRUE, all.y = TRUE)
church_pt2[is.na(church_pt2)] <- 0
church_pt2$Final_dive_depth <- pmax(church_pt2$Max_dive_depth_m, church_pt2$Dive_depth)

#more data
#https://doi.org/10.1093/biolinnean/blaa161
Laeta_data <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Laeta_2020_dive_data.xlsx")
Laeta_data$tips <- str_replace(Laeta_data$Species, " ", "_")
Laeta_data <- Laeta_data[!is.na(Laeta_data$Species), c("Depth (m)", "tips")]
colnames(Laeta_data) <- c("Dive_depth_Laeta", "tips")
Laeta_data <- filter(Laeta_data, Dive_depth_Laeta != "-" )
Laeta_data$Dive_depth_Laeta <- as.integer(Laeta_data$Dive_depth_Laeta)

Dive_depth <- merge(church_pt2, Laeta_data, all.x = TRUE, all.y = TRUE)
Dive_depth[is.na(Dive_depth)] <- 0
Dive_depth$Dive_depth_final <- pmax(Dive_depth$Final_dive_depth, Dive_depth$Dive_depth_Laeta)

#save out only the relevant columns, leaves 65 species
Dive_depth <- Dive_depth %>% filter(Dive_depth_final != 0)
Dive_depth <- Dive_depth[, c("tips", "Mass", "Total_body_length", "Dive_depth_final")]
colnames(Dive_depth) <- c("tips", "Mass", "Body_length", "Dive_depth")
#four entries only have genus (Delphinus_sp, Lagenorhynchus_cruciger, Lissodelphis_sp, Neophocaena_sp,)

write.csv(Dive_depth, here("cetacean_dive_depth.csv"), row.names = FALSE)

# Section 2: Churchill et al, cetacean orbit ratio ---------------------------------

church_pt1 <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt1.xlsx")
colnames(church_pt1) <- c("Species", "Family", "Specimen_number", "Left_orbit_length", "Right_orbit_length", "Bizygomatic_width", "Average_orbit_length", "Orbit_ratio")

church_pt1$tips <- str_replace(church_pt1$Species, " ", "_")
church_pt1 <- church_pt1[!(is.na(church_pt1$Orbit_ratio)), c("Bizygomatic_width", "Average_orbit_length", "Orbit_ratio", "tips")]

#correct species spelling to match
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Cephalorhynchus_commersoni", replacement = "Cephalorhynchus_commersonii")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Berardius_bairdi", replacement = "Berardius_bairdii")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Sagmatias_australis", replacement = "Lagenorhynchus_australis")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Sagmatias_cruciger", replacement = "Lagenorhynchus_cruciger")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Sagmatias_obliquidens", replacement = "Lagenorhynchus_obliquidens")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Cephalorhynchus_heavisidi", replacement = "Cephalorhynchus_heavisidii")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Leucopleurus_acutus", replacement = "Lagenorhynchus_acutus")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Zygorhiza_kochi", replacement = "Zygorhiza_kochii")

#take the mean orbit size since there are some species with duplicates
church_pt1 <- church_pt1 %>% group_by(tips) %>% mutate(Orbit_ratio = mean(Orbit_ratio),  Bizygomatic_width = mean( Bizygomatic_width), Average_orbit_length = mean(Average_orbit_length))

#remove duplicates
church_pt1 <- church_pt1[!duplicated(church_pt1$tips),]

write.csv(church_pt1, here("cetacean_orbit_ratio.csv"), row.names = FALSE)

# Section 3: Baker et al artiodactyla orbit ratio -----------------------------------
artio_eyes <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Baker_2019.xlsx")
artio_eyes <- artio_eyes[2: nrow(artio_eyes),]
colnames(artio_eyes) <- c("tips", "Order", "Corneal_diameter", "Axial_length", "Activity_pattern", "Source")
artio_eyes <- filter(artio_eyes, Order == "Artiodactyla")
artio_eyes <- artio_eyes[artio_eyes$Corneal_diameter != "n/a", c("tips", "Corneal_diameter", "Axial_length") ]
artio_eyes$Axial_length <- as.numeric(artio_eyes$Axial_length)
artio_eyes$Corneal_diameter <- as.numeric(artio_eyes$Corneal_diameter)

#we can look at the ratio of corneal diameter to axial length
#From Hall et al: https://doi.org/10.1098/rspb.2012.2258
#The ratio of corneal diameter to axial length of the eye is a useful measure of relative sensitivity 
#and relative visual acuity that has been used in previous studies as a way to compare animals of disparate size
artio_eyes$Orbit_ratio <- artio_eyes$Corneal_diameter/artio_eyes$Axial_length

sleepy_artio <- read.csv(here("sleepy_artiodactyla_full.csv"))
sleepy_artio <- sleepy_artio[sleepy_artio$tips %in% artio_eyes$tips, c("Family", "tips", "Diel_Pattern")]
artio_eyes <- merge(sleepy_artio, artio_eyes, by = "tips")

write.csv(artio_eyes, here("artio_orbit_ratio.csv"), row.names = FALSE)

# Section 4: Coombs et al habitat, diet, dentition, echo ------------------------------------------------------

echo <- read_xlsx("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_discrete_traits/Coombs_et_al_2021.xlsx")
echo$tips <- str_replace(echo$`Museum ID`, " ", "_")
echo <- separate(echo, tips, into = c("tips", "museum_id"), sep = " ")
echo <- echo[!(is.na(echo$Echo)), c("Age", "Echo", "Diet", "Dentition", "FM", "Habitat", "tips")]

echo <- echo[!duplicated(echo$tips),]

write.csv(echo, here("cetacean_habitat_dentition_echo.csv"), row.names = FALSE)

# Section 5: Unknown source -----------------------------------------------

dive <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\resource.csv")
dive$tips <- dive$Taxon

write.csv(dive, here("cetacean_unknown_dive.csv"), row.names = FALSE)

# Section 6: Manger et al mass, brain size, sociality, longevity, feeding strategy-------------------------------------------------
Manger <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Manger_2013_discrete_cet.csv")
#body mass, encephalization quotient, group size, social dynamics, longevity, feeding strategy
#of these, body mass is probably the most relevant to activity patterns (ie are whales with large body size cathemeral)
colnames(Manger) <- c("Genus_species", "Male_mass", "Female_mass", "Average_body_mass", "Brain_mass", "Encephalization_quotient", "Longevity_days", "Sexual_maturity_days", "Group_size", "Group_size_range", "Group_social_dynamics", "Feeding_strategy")
Manger$tips <- str_replace(Manger$Genus_species, " ", "_")
Manger <- Manger[, c("Male_mass", "Female_mass", "Average_body_mass", "Brain_mass", "Encephalization_quotient", "Longevity_days", "Sexual_maturity_days", "Group_size", "Group_social_dynamics", "Feeding_strategy", "tips")]

Manger[Manger == ""]  <- NA

#delete the rows with no information (from 84 rows to 73)
#Manger <- Manger %>% filter(!is.na(Feeding_strategy) & !is.na(Group_size))

#replace feeding strategy codes with actual strings
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "SK", "Skim")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "SF", "Skim")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "LF", "Lunge")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "SwF", "Swallow")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "R", "Raptorial")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "/OCO", "")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "/CO", "")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "SuF", "Suction")

write.csv(Manger, here("cetacean_manger_et_al.csv"), row.names = FALSE)

# Section 7: Groot et al --------------------------------------------------

groot <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Groot_et_al_2023.xlsx")
#this data is coded so need to decode it with the original paper
#contains trait data that is in the other dataframes
#lifespan, length, mass, brain mass, EQ, age to reproduction, group size, gestation, sociality, group foraging, learned foraging, communication


#Section 8: Churchill qualitative variables
church_qual <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt2_page_2.xlsx")
church_qual <- as.data.frame(church_qual)
church_qual$tips <- str_replace(church_qual$Species, " ", "_")
church_qual <- church_qual[, c("Habitat", "Prey-capture", "tips")]
colnames(church_qual) <- c("Habitat_2", "Prey_capture", "tips")
write.csv(church_qual, here("churchill_habitat_prey_capture.csv"), row.names = FALSE)

# Section 9: Churchill et al dataframe  -----------------------------------
#to keep all data from one source can use just the churchill et al dataset 
dive_depth <- read.csv(here("cetacean_dive_depth.csv")) #65 species
orbit_ratio <- read.csv(here("cetacean_orbit_ratio.csv")) #99 species
habitat <- read.csv(here("churchill_habitat_prey_capture.csv")) #70 species

trait.data <- merge(dive_depth, orbit_ratio, by = "tips", all = TRUE)
trait.data <- merge(trait.data, habitat, all = TRUE)

colnames(trait.data) <- c("tips", "Body_mass", "Body_length", "Dive_depth", "Bizygomatic_width", "Average_orbit_length", "Orbit_ratio", "Habitat", "Prey_capture")

write.csv(trait.data, here("churchill_cetacean_dataset.csv"), row.names = FALSE)

# Section 11: One cetacean dataframe to rule them all -------------------------------
dive_depth <- read.csv(here("cetacean_dive_depth.csv")) #65 species
orbit_ratio <- read.csv(here("cetacean_orbit_ratio.csv")) #99 species
echo <- read.csv(here("cetacean_habitat_dentition_echo.csv")) #194 species
Manger <- read.csv(here("cetacean_manger_et_al.csv")) #84 species
habitat <- read.csv(here("churchill_habitat_prey_capture.csv")) #70 species

#115 species
trait.data <- merge(dive_depth, orbit_ratio, by = "tips", all = TRUE)

trait.data <- merge(trait.data, echo, all = TRUE)

trait.data <- merge(trait.data, habitat, all = TRUE)

trait.data <- merge(trait.data, Manger, all = TRUE)

#filter for only extant species
sleepy_artio <- read.csv(here("sleepy_artiodactyla_full.csv"))

#93 extant species with some data
trait.data <- trait.data %>% filter(tips %in% sleepy_artio$tips)

#46 species have complete data
#trait.data <- trait.data[complete.cases(trait.data),]

#combine the mass data into one column
trait.data.1 <- trait.data[, c("tips","Mass", "Average_body_mass")]
trait.data.1$Average_body_mass <- str_replace_all(trait.data.1$Average_body_mass, "\\,", "")
trait.data.1[is.na(trait.data.1)] <- 0
trait.data.1$Average_body_mass <- as.numeric(trait.data.1$Average_body_mass)

#Manger et al mass is in grams and Churchill mass is in KG. Convert Manger to kg by dividing by 1000
trait.data.1$Average_body_mass <- as.numeric(trait.data.1$Average_body_mass)/1000
#take the largest number as the mass
trait.data.1$Body_mass <- pmax(trait.data.1$Mass, trait.data.1$Average_body_mass)
trait.data.1[trait.data.1 == 0] <- NA
trait.data <- merge(trait.data, trait.data.1[, c("tips", "Body_mass")], all =TRUE)

#consolidate the habitat data
trait.data$Habitat_2 <- str_replace(trait.data$Habitat_2, "Freshwater", "riverine")
trait.data$Habitat_2 <- str_replace(trait.data$Habitat_2, "Nearshore", "coastal")
trait.data$Habitat_2 <- str_replace(trait.data$Habitat_2, "Offshore", "pelagic")

trait.data$Habitat_3 <- paste(trait.data$Habitat, trait.data$Habitat_2, sep = "")
trait.data$Habitat_3 <- str_replace_all(trait.data$Habitat_3, "NA", "")
trait.data$Habitat_3 <- str_replace(trait.data$Habitat_3, "riverineriverine", "riverine")
trait.data$Habitat_3 <- str_replace(trait.data$Habitat_3, "pelagicpelagic", "pelagic")
trait.data$Habitat_3 <- str_replace(trait.data$Habitat_3, "coastalcoastal", "coastal")
trait.data$Habitat_3 <- str_replace(trait.data$Habitat_3, "coastal-pelagiccoastal", "coastal-pelagic")
trait.data$Habitat_3 <- str_replace(trait.data$Habitat_3, "pelagiccoastal", "coastal-pelagic")
trait.data$Habitat <- trait.data$Habitat_3
trait.data[trait.data == ""] <- NA

#prey capture, FM and feeding strategy all have the same information so we can collapse them and resolve conflicts
#keep one column with more detailed information and another col
trait.data.1 <- trait.data[, c("tips", "Prey_capture", "Feeding_strategy", "FM")]
trait.data.1[is.na(trait.data.1)] <- "Unknown"

#four basic categories: skim filter, lunge filter, suction and raptorial
#define as many species into these categories, assign to the feeding method column

unique(trait.data.1$Feeding_strategy)
unique(trait.data.1$FM)
unique(trait.data.1$Prey_capture)

trait.data.1$FM <- str_to_title(trait.data.1$FM)
trait.data.1$FM <- str_replace(trait.data.1$FM, pattern = "Filter", replacement = "Skim")
trait.data.1$FM <- str_replace(trait.data.1$FM, pattern = "Biting", replacement = "Raptorial")

trait.data.1$Feeding_strategy <- str_replace(trait.data.1$Feeding_strategy, pattern = "Swallow", replacement = "Lunge")

trait.data.1$Prey_capture <- str_replace(trait.data.1$Prey_capture, pattern = "Skim", replacement = "Skim")
trait.data.1$Prey_capture <- str_replace(trait.data.1$Prey_capture, pattern = "Suction", replacement = "Suction")
trait.data.1$Prey_capture <- str_replace(trait.data.1$Prey_capture, pattern = "Engulfment", replacement = "Lunge")
trait.data.1$Prey_capture <- str_replace(trait.data.1$Prey_capture, pattern = "Snap", replacement = "Raptorial")
trait.data.1$Prey_capture <- str_replace(trait.data.1$Prey_capture, pattern = "Grip-and-Tear", replacement = "Raptorial")
trait.data.1$Prey_capture <- str_replace(trait.data.1$Prey_capture, pattern = "Ram", replacement = "Raptorial")

trait.data.1$all <- paste(trait.data.1$FM, trait.data.1$Prey_capture, trait.data.1$Feeding_strategy, sep = "")
trait.data.1$all <- str_replace_all(trait.data.1$all, "Unknown", "")
trait.data.1$all <- str_replace(trait.data.1$all, "RaptorialRaptorialRaptorial", "Raptorial")
trait.data.1$all <- str_replace(trait.data.1$all, "RaptorialRaptorial", "Raptorial")
trait.data.1$all <- str_replace(trait.data.1$all, "SuctionSuctionSuction", "Suction")
trait.data.1$all <- str_replace(trait.data.1$all, "SuctionSuction", "Suction")
trait.data.1$all <- str_replace(trait.data.1$all, "SkimSkimSkim", "Skim")
trait.data.1$all <- str_replace(trait.data.1$all, "SkimSkim", "Skim")
trait.data.1$all <- str_replace(trait.data.1$all, "SkimLungeLunge/Skim", "Skim/lunge")
trait.data.1$all <- str_replace(trait.data.1$all, "SkimLungeLunge", "Skim/lunge")
trait.data.1$all <- str_replace(trait.data.1$all, "SkimLunge", "Skim/lunge")
trait.data.1$all <- str_replace(trait.data.1$all, "RaptorialSuctionRaptorial", "Suction/raptorial")
trait.data.1$all <- str_replace(trait.data.1$all, "SuctionRaptorial", "Suction/raptorial")
trait.data.1$all <- str_replace(trait.data.1$all, "RaptorialSuction", "Suction/raptorial")

trait.data.1 <- trait.data.1[, c("tips", "all")]
colnames(trait.data.1) <- c("tips", "Feeding_method")

#merge with full dataset
trait.data <- merge(trait.data, trait.data.1)

#keep only the relevant columns
trait.data <- trait.data[, c("tips", "Body_length", "Body_mass", "Dive_depth", "Orbit_ratio", "Diet", "Habitat", "Prey_capture", "Brain_mass", "Feeding_method")]

#lastly add in the activity patterns
cetaceans_full <- read.csv(here("cetaceans_full.csv"))

cetaceans_full <- cetaceans_full[, c("Parvorder", "Family", "Diel_Pattern", "Confidence", "tips")]

trait.data <- merge(cetaceans_full, trait.data)
trait.data[trait.data == ""] <- NA

#add a colour for each of the families, easier to colour by later
trait.data$fam_colours <- "Unknown"
trait.data[trait.data$Family == "Ziphiidae", c("fam_colours")] <- "green"
trait.data[trait.data$Family == "Delphinidae", c("fam_colours")] <- "olivedrab"
trait.data[trait.data$Family == "Monodontidae", c("fam_colours")] <- "darkred"
trait.data[trait.data$Family == "Kogiidae", c("fam_colours")] <- "firebrick1"
trait.data[trait.data$Family == "Phocoenidae", c("fam_colours")] <- "lightslateblue"
trait.data[trait.data$Family == "Physeteridae", c("fam_colours")] <- "darkslategrey"
trait.data[trait.data$Family == "Neobalaenidae", c("fam_colours")] <- "orange2"
trait.data[trait.data$Family == "Eschrichtiidae", c("fam_colours")] <- "blue1"
trait.data[trait.data$Family == "Balaenidae", c("fam_colours")] <- "darkorchid1"
trait.data[trait.data$Family == "Iniiae", c("fam_colours")] <-"gold"
trait.data[trait.data$Family == "Lipotidae", c("fam_colours")] <-  "khaki3"
trait.data[trait.data$Family == "Balaenopteridae", c("fam_colours")] <- "pink"
trait.data[trait.data$Family == "Platanistidae", c("fam_colours")] <- "turquoise1"

trait.data <- trait.data[!is.na(trait.data$Diel_Pattern),]
trait.data <- trait.data[!is.na(trait.data$Dive_depth),]

ggplot(trait.data, aes(x = Diel_Pattern, y = Dive_depth)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(fill = Family), size = 3, width = 0.1, height = 0, colour = "black", pch = 21) + 
  scale_fill_discrete(values = trait.data$fam_colours) +
  labs(x = "Temporal activity pattern", y = "Dive_depth") + facet_wrap(~Parvorder)


write.csv(trait.data, here("cetacean_ecomorphology_dataset.csv"), row.names = FALSE)
