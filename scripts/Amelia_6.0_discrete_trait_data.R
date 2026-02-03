# Section 1: Churchill et al dive data, mass, length---------------------------------------------
#https://doi.org/10.1111/joa.13522 
church_pt2 <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt2.xlsx")
colnames(church_pt2) <- c("Species", "Dive_depth", "Log_dive_depth", "Dive_duration", "Log_dive_duration", "Mass", "Log_mass", "Total_body_length", "Log_body_length", "Assymetry_index", "Peak_frequency", "Log_frequency")

church_pt2$tips <- str_replace(church_pt2$Species, " ", "_")

#filter for species with dive data, leaves 62 species
church_pt2 <- church_pt2[!(is.na(church_pt2$Dive_depth)), c("Dive_depth", "Dive_duration", "Mass", "Total_body_length", "tips")]

#check for misspellings
church_pt2[!church_pt2$tips %in% mam.tree$tip.label,]

#correct species spelling to match
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Cephalorhynchus_commersoni", replacement = "Cephalorhynchus_commersonii")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Berardius_bairdi", replacement = "Berardius_bairdii")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Sagmatias_australis", replacement = "Lagenorhynchus_australis")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Sagmatias_cruciger", replacement = "Lagenorhynchus_cruciger")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Sagmatias_obliquidens", replacement = "Lagenorhynchus_obliquidens")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Cephalorhynchus_heavisidi", replacement = "Cephalorhynchus_heavisidii")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Leucopleurus_acutus", replacement = "Lagenorhynchus_acutus")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Sagmatias_obscurus", replacement = "Lagenorhynchus_obscurus")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Zygorhiza_kochi", replacement = "Zygorhiza_kochii")
church_pt2$tips <- str_replace(church_pt2$tips, pattern = "Mesoplodon_layardi", replacement = "Mesoplodon_layardii")

write.csv(church_pt2, here("Churchill_dive_body_length.csv"), row.names = FALSE)
# Section 2: Churchill et al, cetacean orbit ratio ---------------------------------
#https://doi.org/10.1111/joa.13522 
church_pt1 <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt1.xlsx")
colnames(church_pt1) <- c("Species", "Family", "Specimen_number", "Left_orbit_length", "Right_orbit_length", "Bizygomatic_width", "Average_orbit_length", "Orbit_ratio")

#only keep extant species
church_pt1 <- church_pt1 %>% filter(Family %in% c("Balaenidae", "Neobalaenidae", "Balaenopteridae", "Physeteridae", "Kogiidae", "Platanistidae", "Ziphiidae", "Lipotidae", "Pontoporiidae", "Iniidae", "Monodontidae", "Phocoenidae", "Delphinidae"))
#manually remove row 15 Zarhachis flagellator, row 65 Semirostrum ceruttii, row 16 Notocetus vanbenedeni, row 51 Parapontoporia sternbergi
church_pt1 <- church_pt1[-c(15, 16, 51, 65), ]

church_pt1$tips <- str_replace(church_pt1$Species, " ", "_")
church_pt1 <- church_pt1[!(is.na(church_pt1$Orbit_ratio)), c("Bizygomatic_width", "Average_orbit_length", "Orbit_ratio", "tips")]

#correct species spelling to match
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Cephalorhynchus_commersoni", replacement = "Cephalorhynchus_commersonii")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Berardius_bairdi", replacement = "Berardius_bairdii")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Sagmatias_australis", replacement = "Lagenorhynchus_australis")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Sagmatias_cruciger", replacement = "Lagenorhynchus_cruciger")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Sagmatias_obliquidens", replacement = "Lagenorhynchus_obliquidens")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Sagmatias_obscurus", replacement = "Lagenorhynchus_obscurus")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Cephalorhynchus_heavisidi", replacement = "Cephalorhynchus_heavisidii")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Leucopleurus_acutus", replacement = "Lagenorhynchus_acutus")
church_pt1$tips <- str_replace(church_pt1$tips, pattern = "Mesoplodon_layardi", replacement = "Mesoplodon_layardii")

#take the mean orbit size since there are some species with duplicates
church_pt1 <- church_pt1 %>% group_by(tips) %>% mutate(Orbit_ratio = mean(Orbit_ratio),  Bizygomatic_width = mean( Bizygomatic_width), Average_orbit_length = mean(Average_orbit_length))

#remove duplicates
church_pt1 <- church_pt1[!duplicated(church_pt1$tips), c("Orbit_ratio", "tips")]

write.csv(church_pt1, here("cetacean_orbit_ratio.csv"), row.names = FALSE)

# Section 3: Coombs et al habitat, diet, dentition, echo ------------------------------------------------------
#https://discovery.ucl.ac.uk/id/eprint/10135933/7/Coombs_10135933_thesis_revised.pdf
echo <- read_xlsx("C:/Users/ameli/OneDrive/Documents/R_projects/cetacean_discrete_traits/Coombs_et_al_2021.xlsx")
echo$tips <- str_replace(echo$`Museum ID`, " ", "_")
echo <- separate(echo, tips, into = c("tips", "museum_id"), sep = " ")

#filter for just extant species
echo <- echo %>% filter(Age == "Extant")

echo <- echo[!(is.na(echo$Echo)), c("Echo", "Diet", "Dentition", "FM", "Habitat", "tips")]

#remove duplicates because all information is identical
echo <- echo[duplicated(echo$tips),]

#check for misspellings
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
echo[!echo$tips %in% mam.tree$tip.label,]
echo$tips <- str_replace(echo$tips, pattern = "Kogia_simus", replacement = "Kogia_sima")

write.csv(echo, here("cetacean_habitat_dentition_echo.csv"), row.names = FALSE)

# Section 4: Travis Park et al, mass, divetype, diet, fm, habitat 2019 -----------------------------------------------
#https://doi.org/10.1186/s12862-019-1525-x

dive <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\resource.csv")
dive$tips <- dive$Taxon

#dive type (shallow (estimated max dive depth <100m), mid (estimated max dive depth ~500m), deep (estimated max dive depth ~1000m), very deep (estimated max dive depth>1000m))
#habitat 1 -all “riverine/nearshore” taxa are classed as “riverine” and “nearshore/oceanic” taxa are classed as “nearshore”
#habitat 2 - “riverine/nearshore” taxa are classed as “nearshore” and “nearshore/oceanic” taxa are classed as “oceanic”
#feeding1 - all “raptorial/suction” taxa are classed as “raptorial”
#feeding2 - all “raptorial/suction” taxa are classed as “suction”

dive <- dive[, c("tips", "Body.size", "Diet", "Divetype", "Feeding.behaviour", "Habitat")]

#check for misspellings
dive[!dive$tips %in% mam.tree$tip.label,]

#change terms to match with other datasets
dive$Habitat <- str_replace(dive$Habitat, "oceanic", "pelagic")
dive$Habitat <- str_replace(dive$Habitat, "nearshore", "coastal")

dive$Habitat <- tolower(dive$Habitat)
dive$Feeding.behaviour <- str_to_title(dive$Feeding.behaviour)

write.csv(dive, here("cetacean_Park_dive.csv"), row.names = FALSE)

# Section 5: Manger et al mass, brain size, sociality, longevity, feeding strategy-------------------------------------------------
#https://doi.org/10.1016/j.neuroscience.2013.07.041

Manger <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Manger_2013_discrete_cet.csv")
#body mass, encephalization quotient, group size, social dynamics, longevity, feeding strategy
#of these, body mass is probably the most relevant to activity patterns (ie are whales with large body size cathemeral)
colnames(Manger) <- c("Genus_species", "Male_mass", "Female_mass", "Average_body_mass", "Brain_mass", "Encephalization_quotient", "Longevity_days", "Sexual_maturity_days", "Group_size", "Group_size_range", "Group_social_dynamics", "Feeding_strategy")
Manger$tips <- str_replace(Manger$Genus_species, " ", "_")

Manger[Manger == ""]  <- NA

#replace feeding strategy codes with actual strings
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "SK", "Skim")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "SF", "Skim")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "LF", "Lunge")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "SwF", "Swallow")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "R", "Raptorial")
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "/OCO", "") #occasional cooperation in feeding
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "/CO", "") #cooperation in feeding
Manger$Feeding_strategy <- str_replace_all(Manger$Feeding_strategy, "SuF", "Suction")

Manger[!Manger$tips %in% mam.tree$tip.label, "tips"]

#fix spelling
Manger$tips <- str_replace(Manger$tips, pattern = "Cephalorhynchus_commersoni", replacement = "Cephalorhynchus_commersonii")
Manger$tips <- str_replace(Manger$tips, pattern = "Deplhinus_delphis", replacement = "Delphinus_delphis")
Manger$tips <- str_replace(Manger$tips, pattern = "Kogia_simus", replacement = "Kogia_sima")
Manger$tips <- str_replace(Manger$tips, pattern = "Mesoplodon_ginkodens", replacement = "Mesoplodon_ginkgodens")
Manger$tips <- str_replace(Manger$tips, pattern = "Lagenirhynchus_obliquidens", replacement = "Lagenorhynchus_obliquidens")

#Lagenorhynchus_australis has only female body mass data and no male, will use this for average body mass
Manger[Manger$tips == "Lagenorhynchus_australis", c("Average_body_mass")] <- Manger[Manger$tips == "Lagenorhynchus_australis", c("Female_mass")]

#delete the rows with no information (from 84 rows to 73)
Manger <- Manger %>% filter_at(vars(Feeding_strategy,Group_size, Average_body_mass),any_vars(!is.na(.)))

#keep only relevant columns
Manger <- Manger[, c("Average_body_mass", "Brain_mass", "Encephalization_quotient", "Longevity_days", "Sexual_maturity_days", "Group_size", "Group_social_dynamics", "Feeding_strategy", "tips")]

#drop row with Delphinus_capensis
Manger <- Manger[-c(18), ]

write.csv(Manger, here("cetacean_manger_et_al.csv"), row.names = FALSE)

# Section 6: Groot et al ignore for now --------------------------------------------------

groot <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Groot_et_al_2023.xlsx")
#this data is coded so need to decode it with the original paper
#contains trait data that is in the other dataframes
#lifespan, length, mass, brain mass, EQ, age to reproduction, group size, gestation, sociality, group foraging, learned foraging, communication


# Section 7: McCurry et al latitude ---------------------------------------
#https://doi.org/10.1093/biolinnean/blac128

McCurry <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\McCurry_2023.xlsx")

McCurry <- as.data.frame(McCurry[, c(6, 8:51)])

test <- McCurry %>% pivot_longer(cols = !`Absolute latitude`, names_to = "Species", values_to = "count")

test <- test %>% filter(count > 0)
test$Species <- str_replace_all(test$Species, pattern = "'", replacement = "")
test$Species <- str_replace_all(test$Species, pattern = " ", replacement = "_")
colnames(test) <- c("Absolute_latitude", "Species", "count")

test <- test %>% group_by(Species) %>% summarize(max_lat = max(Absolute_latitude), mean_lat = mean(Absolute_latitude), min_lat = min(Absolute_latitude))
names <- read.csv(here("cetaceans_full.csv"))

names <- names %>% separate(col = tips, into = c("Genus", "Species"), sep = "_")
test <- test %>% separate(col = Species, into = c("Genus", "Species"), sep = "_")

test <- merge(names, test, by = "Species", all.y =TRUE)
#two species have the species name attentuata, glacialis, australis and hectori. Drop the duplicates
latitude_df <- test[-c(5,8,9,11,31), ]

latitude_df %>% filter(!is.na(max_crep)) %>% ggplot(., aes(x = max_crep, y = max_lat)) + geom_boxplot() + stat_compare_means(method = "anova")

latitude_df$tips <- str_replace(latitude_df$Species_name, pattern = " ", replacement = "_")
latitude_df <- latitude_df[, c("tips", "max_lat", "mean_lat", "min_lat")]

#check if species names are spelled correctly
latitude_df[!latitude_df$tips %in% mam.tree$tip.label,]

#save out 
write.csv(latitude_df, here("cetacean_latitude_df.csv"), row.names = FALSE)

# Section 8: Churchill qualitative variables ------------------------------

church_qual <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Churchill_Baltz_2021_pt2_page_2.xlsx")
church_qual <- as.data.frame(church_qual)
church_qual$tips <- str_replace(church_qual$Species, " ", "_")

#check for misspellings
church_qual[!church_qual$tips %in% mam.tree$tip.label,]

church_qual$tips <- str_replace(church_qual$tips, pattern = "Cephalorhynchus_commersoni", replacement = "Cephalorhynchus_commersonii")
church_qual$tips <- str_replace(church_qual$tips, pattern = "Berardius_bairdi", replacement = "Berardius_bairdii")
church_qual$tips <- str_replace(church_qual$tips, pattern = "Sagmatias_australis", replacement = "Lagenorhynchus_australis")
church_qual$tips <- str_replace(church_qual$tips, pattern = "Sagmatias_cruciger", replacement = "Lagenorhynchus_cruciger")
church_qual$tips <- str_replace(church_qual$tips, pattern = "Sagmatias_obliquidens", replacement = "Lagenorhynchus_obliquidens")
church_qual$tips <- str_replace(church_qual$tips, pattern = "Sagmatias_obscurus", replacement = "Lagenorhynchus_obscurus")
church_qual$tips <- str_replace(church_qual$tips, pattern = "Cephalorhynchus_heavisidi", replacement = "Cephalorhynchus_heavisidii")
church_qual$tips <- str_replace(church_qual$tips, pattern = "Leucopleurus_acutus", replacement = "Lagenorhynchus_acutus")
church_qual$tips <- str_replace(church_qual$tips, pattern = "Zygorhiza_kochi", replacement = "Zygorhiza_kochii")
church_qual$tips <- str_replace(church_qual$tips, pattern = "Mesoplodon_layardi", replacement = "Mesoplodon_layardii")
church_qual <- church_qual[, c("Habitat", "Prey-capture", "tips")]
colnames(church_qual) <- c("Habitat_2", "Prey_capture", "tips")

write.csv(church_qual, here("churchill_habitat_prey_capture.csv"), row.names = FALSE)


# Section 8: Churchill et al full dataframe  -----------------------------------
#to keep all data from one source can use just the churchill et al dataset 
dive_depth <- read.csv(here("Churchill_dive_body_length.csv")) #28 species
orbit_ratio <- read.csv(here("cetacean_orbit_ratio.csv")) #70 species
habitat <- read.csv(here("churchill_habitat_prey_capture.csv")) #70 species

trait.data <- merge(dive_depth, orbit_ratio, by = "tips", all = TRUE)
trait.data <- merge(trait.data, habitat, all = TRUE)

colnames(trait.data) <- c("tips", "Body_mass", "Body_length", "Dive_depth", "Orbit_ratio", "Habitat", "Prey_capture")

write.csv(trait.data, here("churchill_cetacean_dataset.csv"), row.names = FALSE)

# Section 9: Chen et al cetacean data -------------------------------------
#https://doi.org/10.1111/gcb.16385
Chen <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Chen_2022.xlsx")

#Contains active region (AR: 1 = inland water, 2 = coastal waters, 3 = oceanic waters), range size in km2(RS)
#maximum school size (MSS), maximum dive depth in m (MDD), body weight in tonnes (MBW), maximum reproductive cycle (MRC).
# Each threat comprises four impact levels: no = 0, low = 1, moderate = 2 and high = 3.
Chen <- Chen[1:80, c("Species", "IUCN", "AR", "MDD", "MBW")]

#replace codes:1 = inland water, 2 = coastal waters, 3 = oceanic waters
#but use the same language as other databases (inland = riverine, oceanic = pelagic)
Chen$AR[which(Chen$AR == 1)] <- "riverine"
Chen$AR[which(Chen$AR == 2)] <- "coastal"
Chen$AR[which(Chen$AR == 3)] <- "pelagic"

Chen$IUCN[which(Chen$IUCN == 0)] <- "least_concern"
Chen$IUCN[which(Chen$IUCN == 1)] <- "near_threatened"
Chen$IUCN[which(Chen$IUCN == 2)] <- "vulnerable"
Chen$IUCN[which(Chen$IUCN == 3)] <- "endangered"
Chen$IUCN[which(Chen$IUCN == 4)] <- "critically_endangered"

#change columm names
Chen$tips <- str_replace(Chen$Species, pattern = " ", replacement = "_")

colnames(Chen) <- c("Species_name", "IUCN","Active_range", "max_dive_depth", "body_weight", "tips")
Chen <- Chen[, -1]

write.csv(Chen, here("Chen_cetacean_traits.csv"), row.names = FALSE)

# Section 10: Primary literature + Churchill + Laeta + Chen dive depth data --------------------------

#read in the additional dive data I collected
url <- 'https://docs.google.com/spreadsheets/d/1_0ZS_tbddOCckkcKn9H5HpVRDZty4jhkUU20Nc0YYQY/edit?usp=sharing'
dive_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)

#average dive depth is very inconsistent (deep vs shallow dives, day vs night dives)
#and measured inconsistently (average of means across individuals, average of medians etc) so exclude for now
dive_full <- dive_full[,1:9]

dive_full[is.na(dive_full)] <- 0
dive_full$Max_dive_depth_m <- pmax(dive_full$Max_dive_depth_m, dive_full$Alt_Max_1, dive_full$Alt_Max_2, dive_full$Alt_Max_3, dive_full$Alt_Max_4)
dive_full$tips <- str_replace(dive_full$Species, " ", "_")

church_pt2 <- read.csv(here("Churchill_dive_body_length.csv"))
dive.data <- merge(church_pt2[, c("tips", "Dive_depth")], dive_full[, c("Max_dive_depth_m", "tips")], all = TRUE, by = "tips")

#more data from Laeta et al
#https://doi.org/10.1093/biolinnean/blaa161
Laeta_data <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Laeta_2020_dive_data.xlsx")
Laeta_data$tips <- str_replace(Laeta_data$Species, " ", "_")

#check for misspellings
Laeta_data[!Laeta_data$tips %in% mam.tree$tip.label, "tips"]
Laeta_data <- Laeta_data[!Laeta_data$tips %in% c("Neophocaena_sp.", "Lissodelphis_sp.", "Tursiops_t. gephyreus", "Delphinus_sp."),] 

Laeta_data <- Laeta_data[!is.na(Laeta_data$Species), c("Depth (m)", "tips")]
colnames(Laeta_data) <- c("Dive_depth_Laeta", "tips")
Laeta_data <- filter(Laeta_data, Dive_depth_Laeta != "-" )
Laeta_data$Dive_depth_Laeta <- as.integer(Laeta_data$Dive_depth_Laeta)

dive.data <- merge(dive.data, Laeta_data, all = TRUE, by = "tips")

#combine Chen et al dive data with my dive data
Chen <- read.csv(here("Chen_cetacean_traits.csv"))

dive.data <- merge(dive.data, Chen[, c("tips", "max_dive_depth")], all = TRUE, by = "tips")
dive.data[is.na(dive.data)] <- 0

dive.data$Final_dive_depth <- pmax(dive.data$Dive_depth, dive.data$max_dive_depth, dive.data$Max_dive_depth_m, dive.data$Dive_depth_Laeta)

dive.data[dive.data == 0] <- NA

dive.data <- dive.data %>% mutate(., Mean_dive_depth = rowMeans(select(., 2:5), na.rm = TRUE))
dive.data[dive.data == "NaN"] <- NA

cetaceans_full <- read.csv(here("cetaceans_full.csv"))
cetaceans_full <- cetaceans_full[, c("Parvorder", "Family", "max_crep", "tips")]
dive.data <- merge(cetaceans_full, dive.data, by = "tips", all = TRUE)

#does anything associate with activity patterns...anything at all
dive.data %>% filter(!is.na(Final_dive_depth) & !is.na(max_crep)) %>% ggplot(., aes(x = max_crep, y = log(Final_dive_depth))) + geom_boxplot() + geom_point(aes(colour = Family)) + facet_wrap(~Parvorder) + stat_compare_means(method = "anova")

dive.data <- dive.data[, c("tips", "Final_dive_depth", "Mean_dive_depth")]
write.csv(dive.data, here("cetacean_dive_depth.csv"), row.names = FALSE)

# Section 8: One cetacean dataframe to rule them all -------------------------------
dive_depth <- read.csv(here("cetacean_dive_depth.csv")) #65 species
orbit_ratio <- read.csv(here("cetacean_orbit_ratio.csv")) #99 species
echo <- read.csv(here("cetacean_habitat_dentition_echo.csv")) #194 species
Manger <- read.csv(here("cetacean_manger_et_al.csv")) #84 species
habitat <- read.csv(here("churchill_habitat_prey_capture.csv")) #70 species
Parker <- read.csv(here("cetacean_Park_dive.csv")) #48 species

trait.data <- merge(dive_depth, orbit_ratio, by = "tips", all = TRUE)
trait.data <- merge(trait.data, echo, all = TRUE)
trait.data <- merge(trait.data, habitat, all = TRUE)
trait.data <- merge(trait.data, Manger, all = TRUE)

#add in new dataset
colnames(Parker) <- c("tips", "Body.size", "Diet1", "Divetype", "Feeding.behaviour", "Habitat_1")
trait.data <- merge(trait.data, Parker, all = TRUE)

#add in Chen et al data
Chen <- read.csv(here("Chen_cetacean_traits.csv"))
trait.data <- merge(trait.data, Chen, by = "tips", all = TRUE)

#combine the mass data into one column
trait.data.1 <- trait.data[, c("tips","Mass", "Average_body_mass", "Body.size", "body_weight")]
trait.data.1$Average_body_mass <- str_replace_all(trait.data.1$Average_body_mass, "\\,", "")
trait.data.1[is.na(trait.data.1)] <- 0
trait.data.1$Average_body_mass <- as.numeric(trait.data.1$Average_body_mass)

#Manger et al mass is in grams and Churchill and Parker et al mass is in KG.
#Convert Manger to kg by dividing by 1000
trait.data.1$Average_body_mass <- trait.data.1$Average_body_mass/1000
#Chen et al mass is in tonnes (?), convert to kg by multiplying by 1000
trait.data.1$body_weight <- trait.data.1$body_weight * 1000
#take the largest number as the mass
trait.data.1$Body_mass_kg <- pmax(trait.data.1$Mass, trait.data.1$Average_body_mass, trait.data.1$Body.size, trait.data.1$body_weight)
trait.data.1[trait.data.1 == 0] <- NA
trait.data.1$Average_body_mass_kg <- mean(trait.data.1$Mass, trait.data.1$Average_body_mass, trait.data.1$Body.size, trait.data.1$body_weight)
trait.data <- merge(trait.data, trait.data.1[, c("tips", "Body_mass_kg", "Average_body_mass_kg")], all =TRUE)

#consolidate the habitat data
trait.data$Habitat_2 <- str_replace(trait.data$Habitat_2, "Freshwater", "riverine")
trait.data$Habitat_2 <- str_replace(trait.data$Habitat_2, "Nearshore", "coastal")
trait.data$Habitat_2 <- str_replace(trait.data$Habitat_2, "Offshore", "pelagic")

trait.data$Habitat <- str_replace(trait.data$Habitat, "-", "/")

#combine habitat data
trait.data.1 <- trait.data[, c("tips", "Habitat", "Habitat_2", "Habitat_1", "Active_range")]
#take the mode
test <- trait.data.1 %>% rowwise() %>% mutate(Final_habitat = list(DescTools::Mode(c(Habitat, Habitat_2, Habitat_1, Active_range), na.rm = TRUE)))

#check for ties 
test$Final_habitat <- sapply(test$Final_habitat, paste, collapse = "/")
test$Final_habitat <- str_replace(test$Final_habitat, pattern = "coastal/coastal/pelagic", replacement = "coastal/pelagic")
test$Final_habitat <- str_replace(test$Final_habitat, pattern = "coastal/pelagic/pelagic", replacement = "coastal/pelagic")
test[test$tips == "Stenella_frontalis", "Final_habitat"] <- "coastal/pelagic"
test[test$tips == "Eubalaena_glacialis", "Final_habitat"] <- "coastal/pelagic"

#scans for species with a single datapoint (ie it is na in 3 of the 4 rows), replace with non NA value
test1 <- test[rowSums(is.na(test[, 2:5])) == 3, ] 
test1[is.na(test1)] <- ""
test1$Final_habitat <- paste(test1$Habitat, test1$Habitat_2, test1$Habitat_1, test1$Active_range, sep = "")
test[rowSums(is.na(test[, 2:5])) == 3, "Final_habitat"] <- test1[, "Final_habitat"]

test[test == "NA"] <- NA

#replace habitat column with trait.data
trait.data <- merge(trait.data, test[, c("tips", "Final_habitat")], by = "tips")

#prey capture, FM, feeding strategy and feeding.behaviour all have the same information so we can collapse them and resolve conflicts
#keep one column with more detailed information and another col
trait.data.1 <- trait.data[, c("tips", "Prey_capture", "Feeding_strategy", "FM", "Feeding.behaviour")]
#trait.data.1[is.na(trait.data.1)] <- "Unknown"

#four basic categories: skim filter, lunge filter, suction and raptorial
#biting and grip + tear are contained within raptorial, skim (continuous filter), lunge and suction feeding (both intermittent filter) are contained in filter feeding
#source: https://doi.org/10.1098/rspb.2017.1035

trait.data.1$FM <- str_to_title(trait.data.1$FM)
trait.data.1$FM <- str_replace(trait.data.1$FM, pattern = "Biting", replacement = "Raptorial")

#https://doi.org/10.1242/jeb.048157 second source to confirm blue whales do lunge feed (aka swallow feeding)
trait.data.1$Feeding_strategy <- str_replace(trait.data.1$Feeding_strategy, pattern = "Swallow", replacement = "Lunge")

trait.data.1$Prey_capture <- str_replace(trait.data.1$Prey_capture, pattern = "Engulfment", replacement = "Lunge")
trait.data.1$Prey_capture <- str_replace(trait.data.1$Prey_capture, pattern = "Snap", replacement = "Raptorial")
trait.data.1$Prey_capture <- str_replace(trait.data.1$Prey_capture, pattern = "Grip-and-Tear", replacement = "Raptorial")
trait.data.1$Prey_capture <- str_replace(trait.data.1$Prey_capture, pattern = "Ram", replacement = "Raptorial")

test <- trait.data.1 %>% rowwise() %>% mutate(Final_feeding = list(DescTools::Mode(c(Prey_capture, Feeding_strategy, FM, Feeding.behaviour), na.rm = TRUE)))

#check for ties 
test$Final_feeding <- sapply(test$Final_feeding, paste, collapse = "/")

#scans for species with a single datapoint (ie it is na in 3 of the 4 rows), replace with non NA value
test1 <- test[rowSums(is.na(test[, 2:5])) == 3, ] 
test1[is.na(test1)] <- ""
test1$Final_feeding <- paste(test1$FM, test1$Feeding_strategy, test1$Feeding.behaviour, test1$Prey_capture, sep = "")
test[rowSums(is.na(test[, 2:5])) == 3, "Final_feeding"] <- test1[, "Final_feeding"]

#ties between only two sources resolve manually
test[test$tips %in% c("Balaenoptera_musculus", "Megaptera_novaeangliae"), "Final_feeding"] <- "Lunge" #lunge is more informative than filter
test[test$tips == "Eubalaena_glacialis", "Final_feeding"] <- "Skim" #skim is more informative than filter
test[test$tips %in% c("Globicephala_macrorhynchus", "Hyperoodon_planifrons", "Neophocaena_asiaeorientalis", "Orcaella_heinsohni"), "Final_feeding"] <- "Raptorial/Suction"
test[test$tips == "Balaenoptera_borealis", "Final_feeding"] <- "Lunge"

test[test == "NA"] <- NA

#merge with full dataset
trait.data <- merge(trait.data, test[,c("tips", "Final_feeding")], by = "tips")

#consolidate the diet information
#diet is pretty complete so use diet 1 to supplement when a source is missing
trait.data[is.na(trait.data$Diet), c("Diet")] <- trait.data[is.na(trait.data$Diet), c("Diet1")]
trait.data$Diet <- str_replace(trait.data$Diet, "generalist", "cephalopods + fish")

#four species in the tree with no diet information, data added from NOAA fisheries species directory https://www.fisheries.noaa.gov/species-directory
#and from https://marinemammalscience.org/facts/balaenoptera-bonaerensis/#Prey 
trait.data[trait.data$tips == "Eubalaena_japonica", c("Diet")] <- "zooplankton + fish"
trait.data[trait.data$tips == "Stenella_frontalis", c("Diet")] <- "cephalopods + fish"
trait.data[trait.data$tips == "Balaenoptera_bonaerensis", c("Diet")] <- "zooplankton + fish"
trait.data[trait.data$tips == "Stenella_clymene", c("Diet")] <- "cephalopods + fish"

#dive type: shallow <100m, mid ~500m, deep ~1000m, very deep >1000m
trait.data.1 <- trait.data[, c("tips", "Final_dive_depth","Divetype")]
trait.data.1$Final_divetype <- "undetermined"

trait.data.1[trait.data.1$Final_dive_depth > 1000, c("Final_divetype")] <- "verydeep"
trait.data.1[trait.data.1$Final_dive_depth <= 1000, c("Final_divetype")] <- "deep"
trait.data.1[trait.data.1$Final_dive_depth < 500, c("Final_divetype")] <- "mid"
trait.data.1[trait.data.1$Final_dive_depth < 180, c("Final_divetype")] <- "shallow" #calls species that dive up to 166m shallow
trait.data.1[trait.data.1$Final_dive_depth == 0, c("Final_divetype")] <- "undetermined"

trait.data.1[trait.data.1 == 0] <- NA

trait.data <- merge(trait.data, trait.data.1[, c("tips", "Final_divetype")], by = "tips")
trait.data[trait.data$Final_dive_depth == 0, "Final_dive_depth"] <- NA

#keep only the relevant columns
trait.data <- trait.data[, c("tips", "Body_length", "Orbit_ratio", "Diet", "IUCN", "Body_mass_kg", "Final_habitat", "Final_feeding", "Final_dive_depth", "Final_divetype")]
colnames(trait.data) <- c("tips", "Body_length_m","Orbit_ratio", "Diet", "IUCN", "Body_mass_kg", "Habitat", "Feeding_method", "Dive_depth_m", "Divetype")

#add in latitude data for odontocetes
latitude_df <- read.csv(here("cetacean_latitude_df.csv"))

trait.data <- merge(trait.data, latitude_df, by = "tips", all = TRUE)

#lastly add in the activity patterns
cetaceans_full <- read.csv(here("cetaceans_full.csv"))
cetaceans_full <- cetaceans_full[, c("Parvorder", "Family", "Diel_Pattern", "max_crep", "Confidence", "tips")]

trait.data <- merge(cetaceans_full, trait.data, by = "tips", all = TRUE)
trait.data[trait.data == ""] <- NA

#add a colour for each of the families, easier to colour by later
trait.data$fam_colours <- "Unknown"
trait.data[trait.data$Family == "Balaenidae", c("fam_colours")] <- "darkorchid3"
trait.data[trait.data$Family == "Balaenopteridae", c("fam_colours")] <- "magenta"
trait.data[trait.data$Family == "Ziphiidae", c("fam_colours")] <- "green2"
trait.data[trait.data$Family == "Neobalaenidae", c("fam_colours")] <- "pink"
trait.data[trait.data$Family == "Delphinidae", c("fam_colours")] <- "grey"
trait.data[trait.data$Family == "Monodontidae", c("fam_colours")] <- "steelblue4"
trait.data[trait.data$Family == "Eschrichtiidae", c("fam_colours")] <- "yellow3"
trait.data[trait.data$Family == "Iniiae", c("fam_colours")] <-"dodgerblue2"
trait.data[trait.data$Family == "Kogiidae", c("fam_colours")] <- "lightslateblue"
trait.data[trait.data$Family == "Lipotidae", c("fam_colours")] <-  "gold"
trait.data[trait.data$Family == "Phocoenidae", c("fam_colours")] <- "orange2"
trait.data[trait.data$Family == "Physeteridae", c("fam_colours")] <- "green4"
trait.data[trait.data$Family == "Platanistidae", c("fam_colours")] <- "turquoise1"
trait.data[trait.data$Family == "Pontoporiidae", c("fam_colours")] <- "grey30"

write.csv(trait.data, here("cetacean_ecomorphology_dataset.csv"), row.names = FALSE)

# Section 9: eye mass vs body mass ---------------------------------------
eye_mass <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Burton_2006.xlsx")

eye_mass <- separate(eye_mass, col = SpeciesBodymassBrainmassEyemass, into = c("Genus", "Species", "Body_mass_g", "Brain_mass_g", "Eye_mass_g"), sep = " ")
eye_mass <- eye_mass %>% mutate(Eye_mass_g = as.numeric(Eye_mass_g), 
                                Body_mass_g = as.numeric(Body_mass_g)* 1000, 
                                Brain_mass_g = as.numeric(Brain_mass_g), 
                                relative_eye_mass = Eye_mass_g/Body_mass_g)

eye_mass$tips <- paste(eye_mass$Genus, eye_mass$Species, sep = "_")

#fix species names
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Megaptera_nodosa", replacement = "Megaptera_novaeangliae")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Phocaena_phocaena", replacement = "Phocoena_phocoena")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Cervus_axis", replacement = "Axis_axis")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Gazella_thomsonii", replacement = "Eudorcas_thomsonii")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Rangifer_arcticus", replacement = "Rangifer_tarandus")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Cercopithecus_aethiops", replacement = "Chlorocebus_aethiops")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Sciurus_hudsonicus", replacement = "Tamiasciurus_hudsonicus")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Vulpes_fulvus", replacement = "Vulpes_vulpes")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Thos_mesomelas", replacement = "Canis_mesomelas")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Felis_leo", replacement = "Panthera_leo")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Felis_oregonensis", replacement = "Puma_concolor")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Claviglis_saturatus", replacement = "Graphiurus_murinus") #maybe
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Felis_capensis", replacement = "Leptailurus_serval")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Felis_onca", replacement = "Panthera_onca")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Felis_pardus", replacement = "Panthera_pardus")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Mustela_articus", replacement = "Mustela_erminea") #maybe
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Pecari_angulatus", replacement = "Dicotyles_tajacu")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Tapirella_bairdii", replacement = "Tapirus_bairdii")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Thalarctos_maritimus", replacement = "Ursus_maritimus")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Macaca_rhesus", replacement = "Macaca_mulatta")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Rhinocerus_bicornis", replacement = "Diceros_bicornis")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Equus_caballus", replacement = "Equus_ferus")
eye_mass$tips <- str_replace(eye_mass$tips, pattern = "Citellus_parryii", replacement = "Spermophilus_citellus")

#compare eye size across orders
ggplot(eye_mass, aes(x = Order, y = log(relative_eye_mass), colour = Body_mass_g)) +geom_boxplot() + geom_jitter()

#compare to activity patterns
artio_df <- read.csv(here("Sleepy_artiodactyla_full.csv"))
artio_df <- artio_df[artio_df$tips %in% eye_mass$tips, c("tips", "max_crep")]
#16 of 84 species in artio database

mammals_df <- read.csv(here("Bennie_mam_data.csv"))
mammals_df <- mammals_df[mammals_df$tips %in% eye_mass$tips, c("tips", "max_crep")]
#only 41 of 84 species have activity pattern data from Bennie et al
#remove artios that I have data for (use my data)
mammals_df <- mammals_df[!mammals_df$tips %in% artio_df$tips, ]
#removes 10 artio sps

mammals_df <- rbind(mammals_df, artio_df)

eye_mass <- merge(eye_mass, mammals_df, all = TRUE)
eye_mass[eye_mass$tips == "Homo_sapiens", c("max_crep")] <- "diurnal"
eye_mass[eye_mass$tips == "Peromyscus_sp.", c("max_crep")] <- "nocturnal" #the entire genus is nocturnal

#optional add in missing species, most are domestic 
eye_mass[eye_mass$tips == "Felis_domesticus", c("max_crep")] <- "crepuscular" #https://doi.org/10.1007/s10530-017-1534-x
eye_mass[eye_mass$tips == "Canis_familiaris", c("max_crep")] <- "crepuscular" #https://doi.org/10.1016/j.applanim.2021.105449
eye_mass[eye_mass$tips == "Bos_taurus", c("max_crep")] <- "crepuscular" #https://doi.org/10.1371/journal.pone.0313086
eye_mass[eye_mass$tips == "Capra_hircus", c("max_crep")] <- "crepuscular" #https://doi.org/10.1139/z03-055

eye_mass <- eye_mass[!is.na(eye_mass$max_crep), ]

eye_mass %>% filter(max_crep %in% c("diurnal", "nocturnal")) %>% ggplot(., aes(x = log(Body_mass_g), y = log(Eye_mass_g))) + geom_point(aes(colour = max_crep)) +
  geom_smooth(method="lm", na.rm=T, se=F, formula=y~x, aes(color=max_crep)) 

ggplot(eye_mass, aes(x = log(Body_mass_g), y = log(Eye_mass_g))) + geom_point(aes(colour = max_crep)) +
  geom_smooth(method="lm", na.rm=T, se=F, formula=y~x, aes(color=max_crep)) 

ggplot(eye_mass, aes(x = max_crep, y = log(relative_eye_mass))) +geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(colour = tips), size = 2) +
  stat_compare_means(method = "anova") # + facet_wrap(~max_crep)

eye_mass%>% filter(Order %in% c("Ungulates (Artiodactyla)", "Cetacea"), max_crep %in% c("diurnal", "nocturnal", "crepusuclar", "cathemeral")) %>%
  ggplot(., aes(x = max_crep, y = log(relative_eye_mass))) +geom_boxplot(outlier.shape = NA) + geom_jitter(aes(colour = Order)) #+ facet_wrap(~Order)

#pinnipeds have no activity pattern data so we can remove them 
eye_mass%>% filter(max_crep %in% c("diurnal", "nocturnal", "crepuscular", "cathemeral")) %>%
  ggplot(., aes(x = max_crep, y = log(relative_eye_mass))) +geom_boxplot(aes(fill = max_crep), outlier.shape = NA) + geom_point(aes(colour = Order), size =2) +
  stat_compare_means(method = "anova") #+facet_wrap(~Order)

eye_mass%>% filter(Order %in% c("Ungulates (Artiodactyla)", "Cetacea"), max_crep %in% c("diurnal", "nocturnal", "crepusuclar", "cathemeral")) %>%
  ggplot(., aes(x = max_crep, y = log(relative_eye_mass))) +geom_boxplot(outlier.shape = NA) + geom_point(aes(colour = tips)) #+ facet_wrap(~Order)

ggplot(eye_mass, aes(x = max_crep, y = log(Eye_mass_g))) +geom_boxplot(outlier.shape = NA, aes(fill = max_crep)) + 
  geom_point(aes(colour = Order), size = 2) +
  stat_compare_means(method = "anova") # + facet_wrap(~max_crep)


# Section 10: Baker et al artiodactyla orbit ratio -----------------------------------
artio_eyes <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Baker_2019.xlsx")
artio_eyes <- artio_eyes[2: nrow(artio_eyes),]
colnames(artio_eyes) <- c("tips", "Order", "Corneal_diameter", "Axial_length", "Activity_pattern", "Source")
artio_eyes <- filter(artio_eyes, Order %in% c("Artiodactyla"))

artio_eyes <- artio_eyes[artio_eyes$Corneal_diameter != "n/a", c("tips", "Corneal_diameter", "Axial_length", "Activity_pattern", "Order") ]
artio_eyes$Axial_length <- as.numeric(artio_eyes$Axial_length)
artio_eyes$Corneal_diameter <- as.numeric(artio_eyes$Corneal_diameter)

#we can look at the ratio of corneal diameter to axial length
#From Hall et al: https://doi.org/10.1098/rspb.2012.2258
#The ratio of corneal diameter to axial length of the eye is a useful measure of relative sensitivity 
#and relative visual acuity that has been used in previous studies as a way to compare animals of disparate size
artio_eyes$Orbit_ratio <- artio_eyes$Corneal_diameter/artio_eyes$Axial_length

#drop unnecessary columns
artio_eyes <- artio_eyes[, c("tips", "Orbit_ratio")]

# Section 11: Artiodactyla pantheria data ---------------------------------

# library(pak)
# pkg_install("RS-eco/traitdata")
#load library and data
library(traitdata)
data(pantheria)

#filter for artoidactyla
pantheria <- pantheria %>% filter(Order == "Artiodactyla")
pantheria$tips <- paste(pantheria$Genus, pantheria$Species, sep = "_")

#filter for the trait data we're interested in
pantheria <- pantheria[, c("Order", "Family", "AdultBodyMass_g", "AdultHeadBodyLen_mm", "DietBreadth", "HabitatBreadth", "HomeRange_km2", "SocialGrpSize", "Terrestriality", "TrophicLevel", "GR_Area_km2", "GR_MidRangeLat_dd", "GR_MaxLat_dd", "GR_MinLat_dd", "GR_MaxLong_dd", "GR_MinLong_dd", "GR_MidRangeLong_dd", "tips")]

#add in diel pattern data
sleepy_artio <- read.csv(here("sleepy_artiodactyla_full.csv"))

#check for misspellings
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
pantheria[!pantheria$tips %in% mam.tree$tip.label, "tips"]

#fix alternative spellings
pantheria[pantheria$tips == "Hemitragus_hylocrius", "tips"] <- "Nilgiritragus_hylocrius"
pantheria[pantheria$tips == "Hemitragus_jayakari", "tips"] <- "Arabitragus_jayakari"
pantheria[pantheria$tips == "Hexaprotodon_liberiensis", "tips"] <- "Choeropsis_liberiensis"
pantheria[pantheria$tips == "Neotragus_moschatus", "tips"] <- "Nesotragus_moschatus"
pantheria[pantheria$tips == "Przewalskium_albirostris", "tips"] <- "Cervus_albirostris"
pantheria[pantheria$tips == "Pseudois_schaeferi", "tips"] <- "Pseudois_nayaur"
pantheria[pantheria$tips == "Saiga_borealis", "tips"] <- "Saiga_tatarica"
pantheria[pantheria$tips == "Sus_salvanius", "tips"] <- "Porcula_salvania"
pantheria[pantheria$tips == "Alces_americanus", "tips"] <- "Alces_alces"
pantheria[pantheria$tips == "Taurotragus_derbianus", "tips"] <- "Tragelaphus_derbianus"
pantheria[pantheria$tips == "Taurotragus_oryx", "tips"] <- "Tragelaphus_oryx"

#domesticated animals not in mam tree: bos frontalis, bos grunniens, camelus bactrianus, bos taurus, lama glama, ovis aries, capra hircus, bubalus bubalis
#C brookei only recently (?) upgraded to species status from C ogilbyi subspecies
#gazella arabica is an invalid species, synonymous with G erlangeri
#Damaliscus korrigum and Damaliscus superstes are subspecies of D lunatus
#Muntiacus puhoatensis poorly known, thought to be conspecific with M rooseveltorum
# Alcelaphus_caama and Alcelaphus_lichtensteinii are subspecies of Alcelaphus buselaphus
#Babyrousa_bolabatuensis only known from subfossil remains, may be a subspecies

sleepy_artio <- sleepy_artio[sleepy_artio$tips %in% pantheria$tips, c("tips", "max_crep")]
pantheria <- merge(sleepy_artio, pantheria, by = "tips", all = TRUE)

#add in eye size data
artio_eyes <- merge(pantheria, artio_eyes, by = "tips", all = TRUE)

#add to keep colours consistent for each family
artio_eyes$fam_colours <- "Unknown"
artio_eyes[artio_eyes$Family == "Cervidae", c("fam_colours")] <- "maroon"
artio_eyes[artio_eyes$Family == "Bovidae", c("fam_colours")] <- "red"
artio_eyes[artio_eyes$Family == "Camelidae", c("fam_colours")] <- "salmon1"
artio_eyes[artio_eyes$Family == "Rhinocerotidae", c("fam_colours")] <- "grey80"
artio_eyes[artio_eyes$Family == "Equidae", c("fam_colours")] <- "chocolate4"
artio_eyes[artio_eyes$Family == "Tragulidae", c("fam_colours")] <- "dimgrey"
artio_eyes[artio_eyes$Family == "Giraffidae", c("fam_colours")] <- "black"
artio_eyes[artio_eyes$Family == "Tayassuidae", c("fam_colours")] <- "moccasin"
artio_eyes[artio_eyes$Family == "Tapiridae", c("fam_colours")] <- "darkgoldenrod"
artio_eyes[artio_eyes$Family == "Suidae", c("fam_colours")] <- "brown4"
artio_eyes[artio_eyes$Family == "Moschidae", c("fam_colours")] <- "tan1"
artio_eyes[artio_eyes$Family == "Hippopotamidae", c("fam_colours")] <- "skyblue"
artio_eyes[artio_eyes$Family == "Antilocapridae", c("fam_colours")] <- "white"

#save out 
write.csv(artio_eyes, here("artiodactyla_ecomorphology_dataset.csv"), row.names = FALSE)
