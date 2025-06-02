library(rotl)
library(ggtree)
library(stringr)
library(scales)
library(gsheet)
library(ape)
library(patchwork)
library(ggpubr)
library(dplyr)
library(phytools)
library(rfishbase)
library(xlsx)
library(geiger)
library(here)
library(viridis)
library(tidyr)
library(dplyr)
library(readr)

setwd(here())

### LOAD IN OTHER DATA ### 

# # Fetch the taxonomic levels from fishbase for all of the species and make it into a dataframe
# # available_releases()
# # [1] "23.01" "23.05" "21.06" "19.04"
# 
# fishbase_df <- load_taxa(collect = T, version = "21.06")
# fishbase_df <- as.data.frame(fishbase_df)
# # fb_common <- common_names(version = "21.06")

url <- 'https://docs.google.com/spreadsheets/d/18aNqHT73hX06cGRlf6oj7Y4TVKf6jd_Q5ojNIxm2rys/edit#gid=0'
sleepy_fish <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
sleepy_fish$Diel_Pattern <- tolower(sleepy_fish$Diel_Pattern)


#######################################################################################################################################
############################## Test conversion to wide/long ###########################################################################

### Can I convert sleepy_fish from wide to long? By taking the columns representing confidence in different studies (and also columns with different comments and sources?)
### If they all become different rows (1 row for each of the confidence values, comments, and sources), I could then make it wide again
### To determine NEW_SPECIES and HIGHER_CONF, I could just compare it to the archived file from the biorxiv submission

test_fish_long <- sleepy_fish %>% pivot_longer(cols = X5:Source.13, names_to = "column", values_to = "value")

test_fish_long$column <- gsub("\\..*","",test_fish_long$column)

# Remove any row without a value
test_fish_long <- test_fish_long[test_fish_long$value != "",]

## Mutate to deal with crepuscularity
## This needs to be updated, so that new_diel is only if that source calls it crepuscular?
test_fish_long <- test_fish_long %>% group_by(Species_name) %>% mutate(crepuscular = any(grepl("crepuscular", tolower(value))), new_diel = ifelse(column %in% c("X1", "X2", "X3", "X4", "X5"), str_replace(tolower(value), "crepuscular/", ""), NA), new_confidence = ifelse("X5" %in% column, "E", ifelse("X4" %in% column, "D", ifelse("X3" %in% column, "C", ifelse("X1" %in% column, "A", ifelse("X2" %in% column, "B", "no_info"))))))
test_fish_long$new_diel <- str_replace(test_fish_long$new_diel, "crepuscular", "unclear")

## Modify so that the new diel reflects accurately the possibilities (for matching later)
index <- test_fish_long$column %in% c("X1", "X2", "X3", "X4", "X5")
test_fish_long$value <- ifelse(test_fish_long$value == "Crepuscular", "Crepuscular/unclear", test_fish_long$value)


## OK so now I need to write some kind of function
### Then just take the mode of everything?
### Or take mode of X5 + X4, if there isn't a mode (e.g. 1 di and 1 noc), then take mode of all?
### For above, first do X5+X4, then if not, add X3, then if not, add X1, then if not add X2, and if they don't have a clear mode, take mode of everything? And if nothing is mode, make unclear

which.max.simple=function(x,na.rm=TRUE,tie_value="NA"){
  if(na.rm)
  {
    x=x[!is.na(x)]
  }
  if(length(x)==0)
  {
    return(NA)
  }
  maxval=max(x)
  if(is.na(maxval))
  {
    return(NA)
  }
  if(sum(x %in% maxval) > 1)
  {
    # Ties exist, figure out what to do with them.
    if(tie_value=="NA")
    {
      return(NA)
    }
    
    if(tie_value=="random")
    {
      tie_postions=which(x==maxval)
      return(sample(tie_postions,size=1))
    }
    
    if(tie_value=="first")
    {
      tie_postions=which(x==maxval)
      return(tie_postions[1])
    }
    
  } else
  {
    return(which.max(x))
  }
}

tabulateFunc <- function(x) {
  
  if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("X5", "X4")], unique(x$new_diel))), tie_value = "NA"))) {
    if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("X5", "X4", "X3")], unique(x$new_diel))), tie_value = "NA"))) {
      if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("X5", "X4", "X3", "X1")], unique(x$new_diel))), tie_value = "NA"))) {
        if (is.na(which.max.simple(tabulate(match(x$new_diel[x$column %in% c("X5", "X4", "X3", "X1", "X2")], unique(x$new_diel))), tie_value = "NA"))) {
          activity_pattern <- "unclear"
        } else { activity_pattern <- "unclear"} 
      } else {
        activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("X5", "X4", "X3", "X1")], unique(x$new_diel))), tie_value = "NA")]
      }
    } else {
      activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("X5", "X4", "X3")], unique(x$new_diel))), tie_value = "NA")]
    }
  } else {
    activity_pattern <- unique(x$new_diel)[which.max.simple(tabulate(match(x$new_diel[x$column %in% c("X5", "X4")], unique(x$new_diel))), tie_value = "NA")]
  }

  return(activity_pattern)  
}

## This seems to work!!!! I'm not sure the best way to check it, maybe against sleepy_fish$Diel_Pattern?
test <- test_fish_long[test_fish_long$column %in% c("X1", "X2", "X3", "X4", "X5"),] %>% group_by(Species_name) %>% do(tabulated_diel_pattern = tabulateFunc(.)) %>% unnest()

## Will also need to infer the confidence for this call, I think it should probably just be the highest level for that species, which should be easy to get (using above)
test_fish_long$new_diel_tabulated <- NA
test_fish_long$new_diel_tabulated <- ifelse(test_fish_long$crepuscular, paste("crepuscular", test$tabulated_diel_pattern[match(test_fish_long$Species_name, test$Species_name)], sep = "/"), test$tabulated_diel_pattern[match(test_fish_long$Species_name, test$Species_name)])

## Convert back to wide?
## Need to keep in mind the multiple X1s?
sleepy_fish_wide <- test_fish_long %>% pivot_wider(-c(NEW, Common.name, Order, Genome, Assembly, Diel_Pattern, AltDiel.Pattern, Confidence, crepuscular, new_diel), names_from = column, values_from = value)
## Compare to biorxiv submission version

sleepy_fish_bx <- read_csv(file = here("sleepy_fish_database_local_biorxiv_20230519.csv"))
sleepy_fish_firstserach <- read_csv(file = here("sleepy_fish_database_local_2024-09-07.csv"))
table(unique(test_fish_long$Species_name) %in% sleepy_fish_bx$Species_name)
table(unique(test_fish_long$Species_name) %in% sleepy_fish_firstserach$Species_name)

## Everything seems to work, and comparisons to old data mostly agree (or are due to more info for those species)
## The format is now quite different, so I will need to rewrite the stuff above for comparing sources
## For example, all of the X5 examples are in a vector, in the column for each species


## Save out for input into next script

sleepy_fish_wide <- sleepy_fish_wide[,c("Species_name", "new_confidence", "new_diel_tabulated", "Method", "Clock", "X1", "X2", "X3", "X4", "X5", "Comments", "Source")]

saveRDS(sleepy_fish_wide, file = paste("sleepy_fish_wide", nrow(sleepy_fish_wide), "25-04-24.rds", sep = "_"))


#######################################################################################################################################
############################## New method for comparing sources #######################################################################
#######################################################################################################################################

compTwo <- function(comp1 = "comp1", comp2 = "comp2") {
  
  if(any(is.na(c(comp1, comp2)))) {
    return(NA)
  } else {
    comp1 <- str_split(comp1, "/")[[1]]
    comp2 <- str_split(comp2, "/")[[1]]
    
    if(any(comp1 %in% comp2)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
}

## If I take the long format version, I can separate by species, and then take each row and compare them
## and use the values from the columns to keep track of which confidence level

long_db <- test_fish_long[test_fish_long$column %in% c("X1", "X2", "X3", "X4", "X5"),]
long_db$value <- tolower(long_db$value)

species_list <- table(long_db$Species_name)
species_list <- names(species_list[species_list > 1])

output <- lapply(species_list, function(species) {
  
  df <- long_db[long_db$Species_name == species,]
  df$column <- make.unique(df$column)
  
  df_lists_comb <- expand(df, nesting(var = column, vector = value), nesting(var2 = column, vector2 = value), .name_repair = "universal")
  
  df_lists_comb <- df_lists_comb %>% filter(var != var2) %>% arrange(var, var2) %>% mutate(vars = paste0(var, ".", var2)) %>% select(contains("var"), everything())
  
  comparisons <- df_lists_comb %>% group_by(vars) %>% mutate(comp = compTwo(comp1 = vector, comp2 = vector2))
  # comparisons$var <- str_replace(comparisons$var, "\\.[0-50]", "")
  # comparisons$var2 <- str_replace(comparisons$var2, "\\.[0-50]", "")
  comparisons$var <- str_sub(comparisons$var, start = 1, end = 2)
  comparisons$var2 <- str_sub(comparisons$var2, start = 1, end = 2)
  
  comparisons$var_final <- paste(comparisons$var, comparisons$var2, sep = "_")
  
  return(comparisons[,c("comp","var_final")])
})


output <- Reduce(rbind, output)

# Great, this works, and matches roughly the results from before (~72% matching across all categories)
table <- table(output$var_final)
# prop.table(table, margin = 1)
table2 <- as.data.frame(prop.table(table(output$var_final, output$comp), margin = 1))
table2$Comp1 <- sapply(str_split(table2$Var1, "_"), `[`, 1)
table2$Comp2 <- sapply(str_split(table2$Var1, "_"), `[`, 2)
table2 <- table2[table2$Var2 == TRUE,]
table2$count <- table
plot <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = round(Freq, digits = 2))) + geom_tile() + geom_text() + scale_fill_viridis(limits = c(0,1))
plot <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = count)) + geom_tile() + geom_text() + scale_fill_viridis(limits = c(0,1))

plot + xlab("Category (1-5)") + ylab("Category(1-5)")




## This below was old, and trying to get URLs from DOI (not sure why I had to do this, since the entries are doi URLs)
## But this isn't needed anymore

# #######################################################################################################################################
# ############################## Get DOIs for conversion ################################################################################
# #######################################################################################################################################
# 
# # load reticulate so that I can load in the python module
# # then load in pydoi, and use it like pydoi$get_url
# library(reticulate)
# pydoi <- import("pydoi")
# 
# dois <- as.vector(as.matrix(sleepy_fish[,grep("Source.", colnames(sleepy_fish))]))
# dois <- dois[!(dois %in% "")]
# 
# links <- dois[!(grepl("doi.org", dois))]
# links <- links[grepl("https:", links)]
# 
# dois <- dois[grepl("doi.org", dois)]
# dois <- str_replace(dois, "https://doi.org/", replacement = "")
# 
# doi_links1 <- unlist(lapply(dois[1:1000], function(x) {
#   Sys.sleep(0.03)
#   url <- pydoi$get_url(x)
#   return(url)
#   }))
# doi_links2 <- unlist(lapply(dois[1001:2000], function(x) {
#   Sys.sleep(0.03)
#   url <- pydoi$get_url(x)
#   return(url)
# }))
# doi_links3 <- unlist(lapply(dois[2001:3000], function(x) {
#   Sys.sleep(0.03)
#   url <- pydoi$get_url(x)
#   return(url)
# }))
# doi_links4 <- unlist(lapply(dois[3001:4000], function(x) {
#   Sys.sleep(0.04)
#   url <- pydoi$get_url(x)
#   return(url)
# }))
# doi_links5 <- unlist(lapply(dois[4001:5000], function(x) {
#   Sys.sleep(0.03)
#   url <- pydoi$get_url(x)
#   return(url)
# }))
# doi_links6 <- unlist(lapply(dois[5001:6000], function(x) {
#   Sys.sleep(0.03)
#   url <- pydoi$get_url(x)
#   return(url)
# }))h
# doi_links7 <- unlist(lapply(dois[6001:length(dois)], function(x) {
#   Sys.sleep(0.03)
#   url <- pydoi$get_url(x)
#   return(url)
# }))
# 
# final_links <- c(links, doi_links1, doi_links2, doi_links3, doi_links4, doi_links5, doi_links6, doi_links7)
# # lol, most of these are replicates, could have saved a lot of time...
# final_links <- unique(final_links)
# 
# library(clipr)
# write_clip(final_links)
