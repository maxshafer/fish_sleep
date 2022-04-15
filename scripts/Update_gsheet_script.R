library(googlesheets4)
library(rotl)

sleepy_fish <- read_sheet("https://docs.google.com/spreadsheets/d/18aNqHT73hX06cGRlf6oj7Y4TVKf6jd_Q5ojNIxm2rys/edit#gid=0")

species2 <- as.character(unique(sleepy_fish$`Species name`))

resolved_names <- tnrs_match_names(species2)
# resolved_names <- resolved_names[!(is.na(resolved_names$unique_name)),]
# resolved_names <- resolved_names[resolved_names$flags == "",]

# Maybe I should combine the two steps, lapply along the resolved names list, and only keep the 3 tax levels?
phylo_structure <- tax_lineage(taxonomy_taxon_info(ott_id(resolved_names), include_lineage = TRUE))
names <- resolved_names$unique_name[match(names(phylo_structure), resolved_names$ott_id)]

taxonomy <- lapply(phylo_structure, function(x) {
  tax <- list()
  tax[["genus"]] <- x[x$rank == "genus", "unique_name"][1]
  tax[["family"]] <- x[x$rank == "family", "unique_name"][1]
  tax[["order"]] <- x[x$rank == "order", "unique_name"][1]
  #tax <- x[x$rank == c("genus", "family", "order"), "unique_name"]
  return(tax)
})
names(taxonomy) <- names

# Make it into a dataframe?

taxonomy2 <- data.frame(unique_name = names(taxonomy), genus = unlist(lapply(taxonomy, function(x) x$genus)), family = unlist(lapply(taxonomy, function(x) x$family)), order = unlist(lapply(taxonomy, function(x) x$order)))
# OK so now I have the taxonomy for each species name (named by species), map it back and add it to resolved_names?

resolved_names$genus <- taxonomy2$genus[match(resolved_names$unique_name, taxonomy2$unique_name)]
resolved_names$family <- taxonomy2$family[match(resolved_names$unique_name, taxonomy2$unique_name)]
resolved_names$order <- taxonomy2$order[match(resolved_names$unique_name, taxonomy2$unique_name)]

resolved_names_sub <- resolved_names[,c(1,2,8,9,10)]

# Merge resolved_names with sleepy_fish?
sleepy_fish2 <- sleepy_fish

sleepy_fish2$`Species name` <- tolower(sleepy_fish2$`Species name`)
sleepy_fish2 <- merge(sleepy_fish2, resolved_names_sub, by.x = "Species name", by.y = "search_string")


gs4_create("Sleepy_fishies_2_test")
sheet_write(sleepy_fish2, ss = "1msI-PPf8ywelm01kBtleEJf97dNRaro9Jv-LNcsk87s")

googledrive::drive_get("Sleepy_fishies_2_test")
