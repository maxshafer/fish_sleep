# Section 1:Bennie dataframe --------

#read in the Bennie diel activity patterns
#from https://doi.org/10.1073/pnas.1216063110 

# fill in taxonomy
Bennie_mam_data <- read_excel(here("Bennie_diel_activity_data.xlsx"))
colnames(Bennie_mam_data) <- "SpeciesBehaviourReference"
Bennie_mam_data$SpeciesBehaviourReference <- str_replace(string = Bennie_mam_data$SpeciesBehaviourReference, pattern = " ", replacement  = "_")
Bennie_mam_data <- separate(Bennie_mam_data, col = SpeciesBehaviourReference, into = c("tips", "max_crep", "Reference"), sep = " ")
Bennie_mam_data$max_crep <- tolower(Bennie_mam_data$max_crep)
Bennie_mam_data$Species_name <- str_replace(string = Bennie_mam_data$tips, pattern = "_", replacement  = " ")
trait.data <- Bennie_mam_data[1:4477, c("tips", "max_crep", "Species_name") ]

resolved_names <- tnrs_match_names(names = trait.data$Species_name, context_name = "Vertebrates", do_approximate_matching = TRUE)
missing_names <- resolved_names[is.na(resolved_names$ott_id), ] #40 names not found
resolved_names <- resolved_names[!is.na(resolved_names$ott_id), ] #returns 4437 names 

get_rank <- function(tax_info, rank_name) {
  lineage <- tax_lineage(tax_info)[[1]]
  values <- lineage$name[lineage$rank == rank_name]
  if (length(values) == 0) return(NA_character_) else return(values[1])
}

#split it up because it takes so long lol
df1 <- resolved_names[1:1000,] %>% rowwise() %>% mutate(
  tax_info = list(taxonomy_taxon_info(ott_id, include_lineage = TRUE)),
  order = get_rank(tax_info, "order"),
  family = get_rank(tax_info, "family"),
  genus = get_rank(tax_info, "genus")) %>% ungroup() %>% select(-tax_info)

df2 <- resolved_names[1001:2000,] %>% rowwise() %>% mutate(
  tax_info = list(taxonomy_taxon_info(ott_id, include_lineage = TRUE)),
  order = get_rank(tax_info, "order"),
  family = get_rank(tax_info, "family"),
  genus = get_rank(tax_info, "genus")) %>% ungroup() %>% select(-tax_info)

df3 <- resolved_names[2001:3000,] %>% rowwise() %>% mutate(
  tax_info = list(taxonomy_taxon_info(ott_id, include_lineage = TRUE)),
  order = get_rank(tax_info, "order"),
  family = get_rank(tax_info, "family"),
  genus = get_rank(tax_info, "genus")) %>% ungroup() %>% select(-tax_info)

df4 <- resolved_names[3001:4437,] %>% rowwise() %>% mutate(
  tax_info = list(taxonomy_taxon_info(ott_id, include_lineage = TRUE)),
  order = get_rank(tax_info, "order"),
  family = get_rank(tax_info, "family"),
  genus = get_rank(tax_info, "genus")) %>% ungroup() %>% select(-tax_info)

df <- rbind(df1, df2, df3, df4)

df <- df[, c("search_string", "order", "family", "genus")]
df$search_string <- str_to_sentence(df$search_string)
colnames(df) <- c("Species_name","Order", "Family", "Genus")

#fill in missing info
df[df$Family %in% c("Aotidae", "Atelidae", "Cebidae", "Cercopithecidae", "Cheirogaleidae", "Cynocephalidae", "Daubentoniidae", "Galagidae", "Hominidae", "Hylobatidae", "Indriidae", "Lemuridae", "Lepilemuridae", "Lorisidae", "Pitheciidae", "Tarsiidae"), c("Order")] <- "Primates"
df[df$Genus %in% c("Microgale", "Tenrec", "Hemicentetes", "Oryzorictes", "Echinops", "Geogale", "Limnogale", "Setifer"), "Family"] <- "Tenrecidae"
df[df$Genus %in% c("Micropotamogale", "Potamogale"), "Family"] <- "Potamogalidae"
df[df$Family %in% c("Tenrecidae", "Potamogalidae"), "Order"] <- "Afrosoricida"

test2 <- cbind(trait.data, df)  
trait.data <- merge(trait.data, df, by = "Species_name", all = TRUE)

#40 (39?) species weren't found in the otl and so won't have taxonomic info
trait.data[is.na(trait.data$Genus), "Genus"] <- sub(" .*", "", trait.data[is.na(trait.data$Genus), "Species_name"])

#use existing taxonomic info to fill in those species by matching by genus
#this finds info for all but three
for(i in 1:nrow(trait.data)){
  if(is.na(trait.data[i, "Order"])){
    for(j in 1:nrow(trait.data)){
      if(trait.data[i, "Genus"] == trait.data[j, "Genus"] & !is.na(trait.data[j, "Order"])){
        trait.data[i, "Order"] <- trait.data[j, "Order"]
        trait.data[i, "Family"] <- trait.data[j, "Family"]
        break
      }
    }
  }
}

#species with no match, fill in manually (this is more than we need since setting approximate_match = TRUE found a lot of these already but keeping in case)
trait.data[trait.data$Genus %in% c("Smutsia", "Uromanis", "Phataginus"), "Family"] <- "Manidae"
trait.data[trait.data$Family %in% c("Manidae"), "Order"] <- "Pholidota"
trait.data[trait.data$Genus %in% c("Sphiggurus", "Echinoprocta"), "Family"] <- "Erethizontidae"
trait.data[trait.data$Genus %in% c("Loxodontomys", "Phaiomys"), "Family"] <- "Cricetidae"
trait.data[trait.data$Genus %in% c("Megadendromus"), "Family"] <- "Nesomyidae"
trait.data[trait.data$Family %in% c("Erethizontidae", "Cricetidae", "Nesomyidae"), "Order"] <- "Rodentia"
trait.data[trait.data$Genus %in% c("Pseudalopex", "Alopex"), "Family"] <- "Canidae"
trait.data[trait.data$Family %in% c("Canidae"), "Order"] <- "Carnivora"
trait.data[trait.data$Genus %in% c("Enchisthenes", "Lampronycteris", "Trinycteris"), "Family"] <- "Phyllostomidae"
trait.data[trait.data$Genus %in% c("Lissonycteris"), "Family"] <- "Pteropodidae"
trait.data[trait.data$Genus %in% c("Paracoelops"), "Family"] <- "Hipposideridae"
trait.data[trait.data$Genus %in% c("Vespadelus", "Bauerus"), "Family"] <- "Vespertilionidae"
trait.data[trait.data$Family %in% c("Vespertilionidae", "Phyllostomidae", "Pteropodidae", "Hipposideridae"), "Order"] <- "Chiroptera"
trait.data[trait.data$Genus %in% c("Cebuella"), "Family"] <- "Callitrichidae"
trait.data[trait.data$Genus %in% c("Oreonax"), "Family"] <- "Atelidae"
trait.data[trait.data$Family %in% c("Atelidae"), "Order"] <- "Primates"
trait.data[trait.data$Genus %in% c("Choeropsis"), "Family"] <- "Hippopotamidae"
trait.data[trait.data$Genus %in% c("Nesotragus", "Nilgiritragus"), "Family"] <- "Bovidae"
trait.data[trait.data$Family %in% c("Hippopotamidae" ,"Bovidae"), "Order"] <- "Artiodactyla"
trait.data[trait.data$Genus %in% c("Dactylonax"), "Family"] <- "Petauridae"
trait.data[trait.data$Family %in% c("Petauridae"), "Order"] <- "Diprotodontia"

#this species gets mislabeled as a gobi so fix it now
trait.data[trait.data$tips == "Tadarida_sarasinorum", "Family"] <- "Molossidae"
trait.data[trait.data$tips == "Tadarida_sarasinorum", "Order"] <- "Chiroptera"
  
#should have 4477 species
write.csv(trait.data, here("Bennie_mam_data.csv"), row.names = FALSE)

#trait.data <- read.csv(here("Bennie_mam_data.csv"))


# Section 2: Maor dataframe -----------------------------------------------
#read in the Maor diel activity patterns
#from https://doi.org/10.1038/s41559-017-0366-5 
maor_mam_data <- read_excel(here("Maor_diel_activity_data.xlsx"))
maor_mam_data <- maor_mam_data[17:nrow(maor_mam_data), 3:5]
colnames(maor_mam_data) <- c("Species", "Activity_pattern", "Reference")

#the issue with this df is that if they have an alternative activity they add it in a new row
duplicates1 <- maor_mam_data[duplicated(maor_mam_data$Species),]
#make another dataframe since some sps are repeated twice
duplicates2 <- duplicates1[duplicated(duplicates1$Species),] #1080 species have at least one alt diel pattern
duplicates1 <- duplicates1[!duplicated(duplicates1$Species),] #208 species have 2 alt diel patterns

#remove duplicates from Maor dataframe for now, we'll add alternative diel patterns back next
maor_mam_data <-maor_mam_data[!duplicated(maor_mam_data$Species),]
maor_mam_data$Species <- str_replace(string = maor_mam_data$Species, pattern = " ", replacement  = "_")

#format some entries to be max_crep (di/crep -> crep, noc/crep -> crep)
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Diurnal/ Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Nocturnal/Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Nocturnal /Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Diurnal or Cathemeral", "Cathemeral")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Diurnal/Crepuscular", "Crepuscular")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "Nocturnal - EXTINCT", "Nocturnal")
maor_mam_data$Activity_pattern <- str_replace_all(maor_mam_data$Activity_pattern, "EXTINCT", "")

#Don't add back alternative patterns for now
#save out Maor dataframe
write.csv(maor_full, here("maor_full.csv"), row.names  = FALSE)

# #add all the extra diel patterns, then sort the columns after since they're in a random order anyway
maor_mam_data <- merge(maor_mam_data, duplicates1, by='Species', all.x = TRUE, all.y = TRUE)
maor_full <- merge(maor_mam_data, duplicates2, by='Species', all.x = TRUE, all.y = TRUE)
maor_full <- maor_full[, c("Species", "Activity_pattern", "Activity_pattern.x", "Activity_pattern.y")]
maor_full <- relocate(maor_full, "Activity_pattern.x", .after = "Activity_pattern.y")
colnames(maor_full) <- c("Species", "alt_pattern_1", "alt_pattern_2", "Activity_pattern")
maor_full$alt_pattern_1 <- str_replace(maor_full$alt_pattern_1, pattern = "Nocturnal / Arrhythmic", replacement = "Nocturnal/Cathemeral")
maor_full$alt_pattern_2 <- str_replace(maor_full$alt_pattern_2, pattern = "Nocturnal / Arrhythmic", replacement = "Nocturnal/Cathemeral")
maor_full$alt_pattern_2 <- str_replace(maor_full$alt_pattern_2, pattern = "Ultradian", replacement = "Cathemeral")

#filter for just artiodactyls
artio_full <- read.csv(here("sleepy_artiodactyla_full.csv"))

maor_full$tips <- str_replace(maor_full$Species, pattern = " ", replacement = "_")
maor_full <- maor_full[maor_full$tips %in% artio_full$tips,] #leaves 200 species

#lazy solution
maor_full$Diel_pattern <- paste(maor_full$alt_pattern_1, maor_full$alt_pattern_2, maor_full$Activity_pattern, sep = "/")

maor_full$Diel_pattern <- str_replace_all(maor_full$Diel_pattern, pattern = "NA/", replacement = "")
maor_full$Diel_pattern <- str_replace_all(maor_full$Diel_pattern, pattern = "/NA", replacement = "")
maor_full$Diel_pattern <- str_replace_all(maor_full$Diel_pattern, pattern = "NA", replacement = "")

unique(maor_full$Diel_pattern)

maor_full$Diel_pattern <- str_replace_all(maor_full$Diel_pattern, pattern = "Nocturnal/Crepuscular/Nocturnal", replacement = "Nocturnal/Crepuscular")
maor_full$Diel_pattern <- str_replace_all(maor_full$Diel_pattern, pattern = "Crepuscular/Nocturnal", replacement = "Nocturnal/Crepuscular")
maor_full$Diel_pattern <- str_replace_all(maor_full$Diel_pattern, pattern = "Nocturnal/Diurnal", replacement = "Cathemeral")
maor_full$Diel_pattern <- str_replace_all(maor_full$Diel_pattern, pattern = "Diurnal/Crepuscular/Cathemeral", replacement = "Cathemeral/Crepuscular")
maor_full$Diel_pattern <- str_replace_all(maor_full$Diel_pattern, pattern = "Crepuscular/Cathemeral", replacement = "Cathemeral/Crepuscular")

maor_full <- maor_full[, c("Species", "Diel_pattern", "tips")]

#save out Maor dataframe
write.csv(maor_full, here("maor_artio_full.csv"), row.names  = FALSE)

# Section 3: How well do these sources agree? -----------------------------
diel_merge <- merge(Bennie_mam_data,maor_mam_data,by="Species")
colnames(diel_merge) <- c("Species", "Bennie_diel", "Bennie_source", "Maor_diel", "Maor_source")
diel_merge$Bennie_diel <- tolower(diel_merge$Bennie_diel)
diel_merge$Maor_diel <- tolower(diel_merge$Maor_diel)
#set default to idk
diel_merge$match <- "Idk"

for(i in 1:length(diel_merge$Species)){
  if(diel_merge[i, "Bennie_diel"] == diel_merge[i, "Maor_diel"]){
    diel_merge[i, "match"] <- "Yes"
  } else if(diel_merge[i, "Bennie_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal") & diel_merge[i, "Maor_diel"] %in% c("diurnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "diurnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Bennie_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")& diel_merge[i, "Maor_diel"] %in% c("nocturnal/crepuscular", "crepuscular", "cathemeral/crepuscular", "nocturnal")){
    diel_merge[i, "match"] <- "Approximate"
  } else if(diel_merge[i, "Bennie_diel"] %in% c("cathemeral", "cathemeral/crepuscular")& diel_merge[i, "Maor_diel"] %in% c("cathemeral", "cathemeral/crepuscular")){
    diel_merge[i, "match"] <- "Approximate"
  } else {
    diel_merge[i, "match"] <- "No"
  }
}

table(diel_merge$match)

#A surprising amount of consistency!
#170/2199 approximate matches, 244/2199 no matches, and 1785/2199 yes matches
#81% match!
diel_table <- diel_merge %>% count(match)
diel_table <- transform(diel_table, percent = (n/sum(diel_table$n)) * 100)

diel_table <- trait.data.all %>% count(Confidence)
diel_table <- transform(diel_table, percent = (n/sum(diel_table$n)) * 100)

ggplot(diel_table, aes(x="", y=n, fill=Confidence)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + theme_void()

#which species did not match
mismatch_species <- diel_merge %>% filter(match == "No")
table(mismatch_species$Bennie_diel, mismatch_species$Maor_diel)
#we should clean some of these entries up, could be a formatting error in why they don't match

#plot this on the mammal tree, which species tend to have conflicting data?
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))
trait.data <- mismatch_species[,1:2]
trait.data <- trait.data[trait.data$Species %in% mam.tree$tip.label,]
trpy_n_mam <- keep.tip(mam.tree, tip = trait.data$Species)
ggtree(trpy_n_mam, layout = "circular") + geom_tiplab(size = 3)

#seem to be distributed throughout the tree

#get taxonomic information from Maor et al
maor_mam_data <- read_excel(here("Maor_diel_activity_data.xlsx"))
maor_mam_data <- maor_mam_data[17:nrow(maor_mam_data), 1:3]
colnames(maor_mam_data) <- c("Order", "Family", "tips")
maor_mam_data$tips <- str_replace(maor_mam_data$tips, pattern = " ", replacement = "_")

#diel_merge currently has 2,199 species
diel_merge <- merge(diel_merge, maor_mam_data, by = "tips") #now is 3,129 species
diel_merge <- diel_merge[!duplicated(diel_merge$tips),] #remove duplicates, back to 2,199

#save out the final 
colnames(diel_merge) <- c("tips", "Bennie_diel", "Bennie_source", "Maor_diel", "Maor_source", "match", "Order", "Family")
#write.csv(diel_merge, here("sleepy_mammals_old.csv"))
