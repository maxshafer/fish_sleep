library(RISmed)
library(rfishbase)
library(googlesheets4)

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Euteleostomi_deil_patterns/")

PMID_list <- readRDS(file = "sleepy_fish_PMIDs.rds")


## OK, I need to design a search strategy that first searches Orders, then families, then genus, then species
## Combines all the PMIDs by order? There will be duplicates for sure
## Then, remove all the ones that we already have (from PMID scrapper)
## Need to make sure not to query the database too much

# First need to get a list of orders, families, genus, and species

fishbase_df <- load_taxa(collect = T, version = "19.04")
fishbase_df <- as.data.frame(fishbase_df)

# OK so now I have a dataframe with a row for each species, I should probably do nested lists, for order, family, genus, with vectors for each species?
# then I can do a several layered nested for loop? wow, that sounds like a bad idea.
fishbase <- fishbase_df
fishbase <- sapply(unique(fishbase$Order), function(x) fishbase[fishbase$Order == x,], simplify = F, USE.NAMES = T)
fishbase <- lapply(fishbase, function(orders) sapply(unique(orders$Family), function(family) orders[orders$Family == family,], simplify = F, USE.NAMES = T))
fishbase <- lapply(fishbase, function(orders) lapply(orders, function(family) sapply(unique(family$Genus), function(genus) family[family$Genus == genus,"Species"], simplify = F, USE.NAMES = T)))

# Now I have a nested list, with names, ending with a list of species for each genus grouped by family, then by order

# If I run the following for each search term, it will return a dataframe with PMIDs, which I can then use to eliminate doubles and ones we've already scoured
order_tally <- list()
x <- 1
for (fish_genera in names(fishbase)[x:length(names(fishbase))]) {
  # Sleep for 1 second to not overwhelm server
  Sys.sleep(0.25)
  
  family_tally <- list()
  y <- 1
  for (fg in names(fishbase[[fish_genera]])) {
    search_topic <- paste(paste("\"", fg, "\"", sep = ""), "AND (circadian OR diel OR diurnal OR nocturnal OR crepuscular)", sep = " ")
    search_query <- EUtilsSummary(search_topic, type="esearch", db="pubmed")
    records <- EUtilsGet(search_query@PMID, type = "efetch", db = "pubmed")
    if (!is.na(records@PMID)) {
      record_df <- data.frame('PMID' = PMID(records), 'DOI' = DOI(records), 'Title' = ArticleTitle(records), 'Abstract' = AbstractText(records), 'Year' = YearPubmed(records), 'Journal' = Title(records))
      family_tally[[y]] <- record_df
    } else {
      family_tally[[y]] <- "NA"
    }
    y <- y + 1
  }
  
  family_tally <- Reduce(rbind, family_tally[unlist(lapply(family_tally, function(x) is.data.frame(x)))])
  family_tally <- family_tally[!duplicated(family_tally$PMID),]
  
  # search the topic, fish_genera
  search_topic <- paste(paste("\"", fish_genera, "\"", sep = ""), "AND (sleep OR diel OR diurnal OR nocturnal OR crepuscular)", sep = " ")
  search_query <- EUtilsSummary(search_topic, type="esearch", db="pubmed")

  # return the df
  records <- EUtilsGet(search_query@PMID, type = "efetch", db = "pubmed")
  if (!is.na(records@PMID)) {
    record_df <- data.frame('PMID' = PMID(records), 'DOI' = DOI(records), 'Title' = ArticleTitle(records), 'Abstract' = AbstractText(records), 'Year' = YearPubmed(records), 'Journal' = Title(records))
    order_tally[[x]] <- record_df
  } else {
    order_tally[[x]] <- "NA"
  }
  order_tally[[x]] <- rbind(order_tally[[x]], family_tally)
  x <- x + 1
}


orders <- Reduce(rbind, order_tally[unlist(lapply(order_tally, function(x) is.data.frame(x)))])

orders <- orders[!duplicated(orders$PMID),]
orders <- orders[!(orders$PMID %in% PMID_list),]

saveRDS(orders, file = "PMID_results_orders.rds")

# Works well, give's me ~850 PMIDs for orders and family based searches

## Maybe run it for unique IDs in the fishbase_df

genus_tally <- list()
x <- 1
for (fish_genera in unique(fishbase_df$Genus)[x:length(unique(fishbase_df$Genus))]) {
  # Sleep for 1 second to not overwhelm server
  Sys.sleep(0.25)
  
  # search the topic, fish_genera
  search_topic <- paste(paste("\"", fish_genera, "\"", sep = ""), "AND (circadian OR diel OR diurnal OR nocturnal OR crepuscular)", sep = " ")
  search_query <- EUtilsSummary(search_topic, type="esearch", db="pubmed")
  
  # return the df
  records <- EUtilsGet(search_query@PMID, type = "efetch", db = "pubmed")
  if (!is.na(records@PMID)) {
    record_df <- data.frame('PMID' = PMID(records), 'DOI' = DOI(records), 'Title'  = ArticleTitle(records), 'Abstract' = AbstractText(records), 'Year' = YearPubmed(records), 'Journal' = Title(records))
    genus_tally[[x]] <- record_df
  } else {
    genus_tally[[x]] <- "NA"
  }
  x <- x + 1
}


genus <- Reduce(rbind, genus_tally[unlist(lapply(genus_tally, function(x) is.data.frame(x)))])

genus <- genus[!duplicated(genus$PMID),]
genus <- genus[!(genus$PMID %in% PMID_list),]

saveRDS(genus, file = "PMID_results_genus.rds")

# Maybe also remove entries where we already have species info? lapply a grep searching for each species in title or abstract?


# Include a tryCatch in case the server is busy or doesn't work (delays 2 secs, then moves on - genera which aren't searched are given an "ERROR" value)

# search_pubmed2 <- function(search_topic) {
  out <- tryCatch(
    expr = {
      
      records <- EUtilsGet(EUtilsSummary(search_topic, type="esearch", db="pubmed")@PMID, type = "efetch", db = "pubmed")
      if (!is.na(records@PMID)) {
        record_df <- data.frame('PMID' = PMID(records), 'DOI' = DOI(records), 'Title' = ArticleTitle(records), 'Abstract' = AbstractText(records), 'Year' = YearPubmed(records), 'Journal' = Title(records))
        species_tally[[x]] <- record_df
      } else {
        species_tally[[x]] <- "NA"
      }
      species_tally
    },
    error = function(e) {
      message("Server busy")
      print(e)
      return("ERROR")
      Sys.sleep(2)
    },
    warning = function(w) {
      message("Warning message happened")
      print(w)
    },
    finally = {
      message(paste("Finished fish #", x, sep = " "))
    }
  )
  return(out)
}

species_tally <- list()
x <- 1
# I'm starting again, but will make a new list to then append with the old, starting from the species I left off at

x <- grep(names(species_tally)[12804], unique(fishbase_df$Species))
for (fish_genera in unique(fishbase_df$Species)[x:length(unique(fishbase_df$Species))]) {
  # Sleep for 1 second to not overwhelm server
  #Sys.sleep(0.1)
  
  if (is.na(fish_genera) | nchar(fish_genera) < 5) {
    x <- x + 1
    fish_genera <- unique(fishbase_df$Species)[x]
  }
   
  # search the topic, fish_genera
  search_topic <- paste(paste("\"", fish_genera, "\"", sep = ""), "AND (sleep OR diel OR diurnal OR nocturnal OR crepuscular)", sep = " ")
  search_query <- EUtilsSummary(search_topic, type="esearch", db="pubmed")
  #search_query <- search_topic2(search_topic)
  
  # return the df
  records <- EUtilsGet(search_query@PMID, type = "efetch", db = "pubmed")
  if (!is.na(records@PMID)) {
    record_df <- data.frame('PMID' = PMID(records), 'DOI' = DOI(records), 'Title' = ArticleTitle(records), 'Abstract' = AbstractText(records), 'Year' = YearPubmed(records), 'Journal' = Title(records))
    species_tally[[x]] <- record_df
  } else {
    species_tally[[x]] <- "NA"
  }
  print(paste("Finished genera", x, sep = " "))
  x <- x + 1
}

# Saved in users/maxwellshafer

species_tally_1 <- readRDS("PMID_results_species_12804.rds")

saveRDS(species_tally, file = "PMID_results_species_34299_startfrom12804.rds")

species_tally <- c(species_tally_1, species_tally)

species <- Reduce(rbind, species_tally[unlist(lapply(species_tally, function(x) is.data.frame(x)))])

species <- species[!duplicated(species$PMID),]
species <- species[!(species$PMID %in% PMID_list),]

saveRDS(species, file = "PMID_results_species.rds")

## Read them back in and combine them together

orders <- readRDS("PMID_results_orders.rds")
genus <- readRDS("PMID_results_genus.rds")

combined <- Reduce(rbind, list(orders, genus, species))
combined <- combined[!duplicated(combined$PMID),]
combined <- combined[!(combined$PMID %in% PMID_list),]

combined$url <- paste("https://pubmed.ncbi.nlm.nih.gov/", combined$PMID, "/", sep = "")

# Mark entries we've already checked

checked <- read_sheet("https://docs.google.com/spreadsheets/d/1FAFU5U4ZA1CoX_0s_Eu20R5ehm9SPQ5XU1q6GCzC1Rw/edit#gid=72790608")

combined$checked <- checked$checked[match(combined$PMID, checked$PMID)]

# Maybe also remove entries where we already have species info? lapply a grep searching for each species in title or abstract?

saveRDS(combined, file = "combined_PMID_results.rds")

# Maybe prefilter for known words that won't match, like patient, mice, mouse, rat, rats

false.positives <- grep("patient", combined$Abstract, ignore.case = T)
false.positives <- c(false.positives, grep("patient", combined$Title, ignore.case = T))
false.positives <- c(false.positives, grep("human", combined$Abstract, ignore.case = T))
false.positives <- c(false.positives, grep("human", combined$Title, ignore.case = T))
false.positives <- c(false.positives, grep("children", combined$Abstract, ignore.case = T))
false.positives <- c(false.positives, grep("children", combined$Title, ignore.case = T))
false.positives <- c(false.positives, grep("disease", combined$Abstract, ignore.case = T))
false.positives <- c(false.positives, grep("disease", combined$Title, ignore.case = T))
false.positives <- c(false.positives, grep("disorder ", combined$Abstract, ignore.case = T))
false.positives <- c(false.positives, grep("disorder ", combined$Title, ignore.case = T))
false.positives <- c(false.positives, grep("apnoea", combined$Abstract, ignore.case = T))
false.positives <- c(false.positives, grep("apnoea", combined$Title, ignore.case = T))
false.positives <- c(false.positives, grep("epilepsy", combined$Abstract, ignore.case = T))
false.positives <- c(false.positives, grep("epilepsy", combined$Title, ignore.case = T))
false.positives <- c(false.positives, grep(" mice", combined$Abstract))
false.positives <- c(false.positives, grep(" mice", combined$Title))
false.positives <- c(false.positives, grep(" mouse", combined$Abstract))
false.positives <- c(false.positives, grep(" mouse", combined$Title))
false.positives <- c(false.positives, grep(" mammal", combined$Abstract))
false.positives <- c(false.positives, grep(" mammal", combined$Title))
false.positives <- c(false.positives, grep(" rat ", combined$Abstract))
false.positives <- c(false.positives, grep(" rat ", combined$Title))
false.positives <- c(false.positives, grep(" rats ", combined$Abstract))
false.positives <- c(false.positives, grep(" rats ", combined$Title))
false.positives <- c(false.positives, grep(" Drosophila", combined$Abstract))
false.positives <- c(false.positives, grep(" Drosophila", combined$Title))
false.positives <- c(false.positives, grep(" bird", combined$Abstract))
false.positives <- c(false.positives, grep(" bird", combined$Title))
false.positives <- c(false.positives, grep(" reptile", combined$Abstract))
false.positives <- c(false.positives, grep(" reptile", combined$Title))
false.positives <- c(false.positives, grep("REM", combined$Abstract))
false.positives <- c(false.positives, grep("REM", combined$Title))
false.positives <- c(false.positives, grep("EEG", combined$Abstract))
false.positives <- c(false.positives, grep("EEG", combined$Title))
false.positives <- c(false.positives, grep("education", combined$Abstract))
false.positives <- c(false.positives, grep("education", combined$Title))

length(unique(false.positives))

false.positives <- sort(unique(false.positives))

combined$checked[false.positives] <- "fp"
combined$checked[is.na(combined$checked)] <- ""

# Write to a gsheet

# gs4_create("PMID_output")
# sheet_write(combined, ss = "1tK0taHqjqEtaBtVWO_YhHGWArwqT7XZp7XwGuXZeNb8")

write.csv(combined, file = "combined_PMID_results.csv")

write.csv(combined[combined$checked == "",], file = "combined_PMID_results_filtered.csv")




