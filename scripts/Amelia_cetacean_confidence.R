##Packages we will use ---------------------------------------------------
# For manipulating character strings
library(stringr)
# For controlling working directory
library(here)
# For plotting phylogenetic trees
library(ggplot2)
library(ggtree)
# For loading from google sheet
library(gsheet)
#saving dataframes
library(gridExtra)
#manipulating dataframes
library(dplyr)
#reading in excel sheets
library(readxl)
#open tree of life
library(rotl)
#adds timescale
library(deeptime)
#colours
library(RColorBrewer)
#apply two separate colour palettes
library(ggnewscale)
#more colours
library(pals)
#useful
library(tidyr)
#colours
library(viridis)

## Packages for phylogenetic analysis in R (overlapping uses)
## They aren't all used here, but you should have them all
library(ape) 
library(phytools)
library(geiger)
library(corHMM)
library(phangorn)

# Set the working directory and source the functions (not used yet)
setwd(here())

source("scripts/fish_sleep_functions.R")

# Section 1: Plot data on tree --------------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1eG_WIbhDzSv_g-PY90qpTMteESgPZZZt772g13v-H1o/edit?usp=sharing'
diel_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

#diel_full <- diel_full %>% filter(Parvorder == "Mysticeti")

diel_full$tips <- str_replace(diel_full$Species_name, pattern = " ", replacement = "_")

mam.tree <- readRDS(here("maxCladeCred_mammal_tree.rds"))

trait.data <- diel_full[diel_full$tips %in% mam.tree$tip.label, -c(2, 15:21)]
trpy_n <- keep.tip(mam.tree, tip = trait.data$tips)

custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "red", "black", "blue")
# custom.colours <- c("#dd8ae7",  "peachpuff2", "#FC8D62", "#66C2A5")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "Diel_Pattern_3", "New_Pattern")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = Diel_Pattern_3), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot +  new_scale_fill() + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x +2, y=y, fill = New_Pattern), inherit.aes = FALSE, colour = "transparent") 
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

custom.colours <- c("#dd8ae7", "#FC8D62", "#fbbe30", "#66C2A5", "#A6D854", "red", "black", "blue")
diel.plot <- ggtree(trpy_n, layout = "circular") %<+% trait.data[,c("tips", "New_Pattern")]
diel.plot <- diel.plot + geom_tile(data = diel.plot$data[1:length(trpy_n$tip.label),], aes(x=x, y=y, fill = New_Pattern), inherit.aes = FALSE, colour = "transparent") + scale_fill_manual(values = custom.colours, name = "Temporal activity pattern")
diel.plot <- diel.plot + geom_tiplab(size = 2, offset = 3)
diel.plot

png("C:/Users/ameli/OneDrive/Documents/R_projects/Amelia_data_diel_plots/new_cetacean_diel.png", width=20,height=15,units="cm",res=1200)
diel.plot
dev.off()

# Section 2: How well does the data agree -------------------------------
#take only the columns we're interested in
diel_full <- diel_full[, c(1,3,4,8:15)]

#convert strings to lowercase
diel_full$Conf_1 <- tolower(diel_full$Conf_1)
diel_full$Conf_2_daytime <- tolower(diel_full$Conf_2_daytime)
diel_full$Conf_2_nighttime <- tolower(diel_full$Conf_2_nighttime)
diel_full$Conf_3_PAM <- tolower(diel_full$Conf_3_PAM)
diel_full$Conf_3_Stomach_bycatch <- tolower(diel_full$Conf_3_Stomach_bycatch)
diel_full$Conf_4 <- tolower(diel_full$Conf_4)
diel_full$Conf_5 <- tolower(diel_full$Conf_5)

#separate all entries into each category into a separate column
diel_full <- separate(data = diel_full, col = Conf_1, into = c("Conf1.1", "Conf1.2", "Conf1.3", "Conf1.4", "Conf1.5"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_2_daytime, into = c("Conf2.1", "Conf2.2", "Conf2.3", "Conf2.4", "Conf2.5"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_2_nighttime, into = c("Conf2N.1", "Conf2N.2"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_3_PAM, into = c("Conf3.1", "Conf3.2", "Conf3.3", "Conf3.4", "Conf3.5", "Conf3.6", "Conf3.7", "Conf3.8"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_3_Stomach_bycatch, into = c("Conf3ByS.1", "Conf3ByS.2"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_4, into = c("Conf4.1", "Conf4.2", "Conf4.3", "Conf4.4", "Conf4.5", "Conf4.6", "Conf4.7", "Conf4.8"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_5, into = c("Conf5.1", "Conf5.2", "Conf5.3", "Conf5.4"), sep = ",")

#how well do entries in the same category agree?

#filter for category 1 data
diel_full_C1 <- diel_full[, c("Conf1.1", "Conf1.2", "Conf1.3", "Conf1.4", "Conf1.5")]
#filter out rows with only one entry (nothing to compare)
diel_full_C1 <- diel_full_C1[!is.na(diel_full_C1$Conf1.2),]

#five species with multiple level 1 evidence

#str_trim removes whitespace from entry
diel_full_C1$match <- diel_full_C1$Conf1.1 == str_trim(diel_full_C1$Conf1.2)



# Section 3: Compare data sources using Max's code ------------------------
url <- 'https://docs.google.com/spreadsheets/d/1eG_WIbhDzSv_g-PY90qpTMteESgPZZZt772g13v-H1o/edit?usp=sharing'
diel_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
diel_full <- diel_full[, c(1,3,4,8:15)]
diel_full$Conf_1 <- tolower(diel_full$Conf_1)
diel_full$Conf_2_daytime <- tolower(diel_full$Conf_2_daytime)
diel_full$Conf_2_nighttime <- tolower(diel_full$Conf_2_nighttime)
diel_full$Conf_3_PAM <- tolower(diel_full$Conf_3_PAM)
diel_full$Conf_3_Stomach_bycatch <- tolower(diel_full$Conf_3_Stomach_bycatch)
diel_full$Conf_4 <- tolower(diel_full$Conf_4)
diel_full$Conf_5 <- tolower(diel_full$Conf_5)
diel_full <- separate(data = diel_full, col = Conf_1, into = c("Conf1.1", "Conf1.2", "Conf1.3", "Conf1.4", "Conf1.5"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_2_daytime, into = c("Conf2.1", "Conf2.2", "Conf2.3", "Conf2.4", "Conf2.5"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_2_nighttime, into = c("Conf2N.1", "Conf2N.2"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_3_PAM, into = c("Conf3.1", "Conf3.2", "Conf3.3", "Conf3.4", "Conf3.5", "Conf3.6", "Conf3.7", "Conf3.8"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_3_Stomach_bycatch, into = c("Conf3ByS.1", "Conf3ByS.2"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_4, into = c("Conf4.1", "Conf4.2", "Conf4.3", "Conf4.4", "Conf4.5", "Conf4.6", "Conf4.7", "Conf4.8"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_5, into = c("Conf5.1", "Conf5.2", "Conf5.3", "Conf5.4"), sep = ",")

#convert into long format
diel_full_long <- diel_full %>% pivot_longer(cols = Conf1.1:Conf5.4, names_to = "column", values_to = "value")

#remove whitespace
diel_full_long$value <- str_trim(diel_full_long$value)

#remove rows with empty values
diel_full_long <- diel_full_long[diel_full_long$value != "",]
diel_full_long <- diel_full_long[!is.na(diel_full_long$value),]

#interested in how this breaks down between odontocetes and mysticetes
#diel_full_long <- filter(diel_full_long, Parvorder == "Odontoceti")
#diel_full_long <- filter(diel_full_long, Parvorder == "Mysticeti")

#for now convert annoying values, REVISIT THIS 
unique(diel_full_long$value)
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "reverse-dvm", replacement = "unclear")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "dvm/nocturnal", replacement = "nocturnal")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "dvm", replacement = "unclear")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "tbd", replacement = "unclear")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "cathemeral-variable", replacement = "cathemeral")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "cathemeral-invariate", replacement = "cathemeral")

#calling these nocturnal, diurnal, cathemeral (removing unclear as an option)
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear-nocturnal", replacement = "nocturnal")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear-diurnal", replacement = "diurnal")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear-cathemeral", replacement = "cathemeral")

#calling these nocturnal, diurnal, crepuscular, etc may switch to cathemeral later
#idk what to do about this one
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "weakly-nocturnal/crepuscular", replacement = "nocturnal/crepuscular")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "weakly-diurnal", replacement = "diurnal/cathemeral")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "weakly-crepuscular", replacement = "crepuscular/cathemeral")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "weakly-nocturnal", replacement = "nocturnal/cathemeral")

#optional: remove unclear as an option since it gives no new information
diel_full_long[diel_full_long == "unclear"] <- NA
diel_full_long <- diel_full_long[!is.na(diel_full_long$value),]

#use below for four state data, will also need to make a call on partially cathemeral species
# diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "nocturnal/crepuscular", replacement = "nocturnal")
# diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "diurnal/crepuscular", replacement = "diurnal")
# diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#get a list of all the species with more than one source (should be most of them)
species_list <- table(diel_full_long$Species_name)
#should be 76 species
species_list <- names(species_list[species_list > 1])


#function Max wrote for comparing entries
compTwo <- function(comp1 = "comp1", comp2 = "comp2") {
  
  if(any(is.na(c(comp1, comp2)))) {
    return(NA)
  } else {
    #splits any entries with a backslash into two components (ie nocturnal/crepuscular into nocturnal and crepuscular)
    comp1 <- str_split(comp1, "/")[[1]]
    comp2 <- str_split(comp2, "/")[[1]]
    #then compares if any of the components match
    if(any(comp1 %in% comp2)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
}


#apply this function across all species with multiple entries
output <- lapply(species_list, function(species) {
  
  #filter for one species at a time
  df <- diel_full_long[diel_full_long$Species_name == species,]
  #rename the column names to be unique for every entry for this species (ie for multiple column 2 entries column 2.1, 2.2 etc)
  df$column <- make.unique(df$column)
  
  #converts the dataframe so it compares every entry with each other (ie for A,B,C A-A, A-B, A-C, B-A, B-B, B-C, etc)
  df_lists_comb <- expand(df, nesting(var = column, vector = value), nesting(var2 = column, vector2 = value), .name_repair = "universal")
  
  #??? idk
  df_lists_comb <- df_lists_comb %>% filter(var != var2) %>% arrange(var, var2) %>% mutate(vars = paste0(var, ".", var2)) %>% select(contains("var"), everything())
  
  #evaluates the activity patterns for each of these sources and returns if they agree or not (TRUE or FALSE)
  comparisons <- df_lists_comb %>% group_by(vars) %>% mutate(comp = compTwo(comp1 = vector, comp2 = vector2))
  #manipulate the strings for both variable names to revert them back to the original name (back to column 2 from col 2.1)
  comparisons$var <- str_sub(comparisons$var, start = 1, end = 5)
  comparisons$var2 <- str_sub(comparisons$var2, start = 1, end = 5)
  
  #create a column returning the comparison being made (ie col2-col2, col1-col2, etc)
  comparisons$var_final <- paste(comparisons$var, comparisons$var2, sep = "_")
  
  #return just the comparison result column (TRUE or FALSE match) and the comparison being made (ie col1 vs col1)
  return(comparisons[,c("comp","var_final")])
})

#combine this list of results 
output <- Reduce(rbind, output)

table <- table(output$var_final)
# prop.table(table, margin = 1)
table2 <- as.data.frame(prop.table(table(output$var_final, output$comp), margin = 1))
table2$Comp1 <- sapply(str_split(table2$Var1, "_"), `[`, 1)
table2$Comp2 <- sapply(str_split(table2$Var1, "_"), `[`, 2)
table2 <- table2[table2$Var2 == TRUE,]
table2$count <- table
plot_freq <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = round(Freq, digits = 2))) + geom_tile() + geom_text() + scale_fill_viridis(limits = c(0,1))
plot_freq
plot_count <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = count)) + geom_tile() + geom_text() + scale_fill_viridis(limits = c(0,1))
plot_count

#things this doesn't account for: species with some diurnal calls, some nocturnal calls, some cathemeral calls that are categorized as cathemeral-variable


# #Section 4: Cetacean confidence sankey ----------------------------------
library(ggsankey)
url <- 'https://docs.google.com/spreadsheets/d/1eG_WIbhDzSv_g-PY90qpTMteESgPZZZt772g13v-H1o/edit?usp=sharing'
diel_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
diel_full <- diel_full[, c(1,3,4,8:15)]
diel_full$Conf_1 <- tolower(diel_full$Conf_1)
diel_full$Conf_2_daytime <- tolower(diel_full$Conf_2_daytime)
diel_full$Conf_2_nighttime <- tolower(diel_full$Conf_2_nighttime)
diel_full$Conf_3_PAM <- tolower(diel_full$Conf_3_PAM)
diel_full$Conf_3_Stomach_bycatch <- tolower(diel_full$Conf_3_Stomach_bycatch)
diel_full$Conf_4 <- tolower(diel_full$Conf_4)
diel_full$Conf_5 <- tolower(diel_full$Conf_5)
diel_full <- separate(data = diel_full, col = Conf_1, into = c("Conf1.1", "Conf1.2", "Conf1.3", "Conf1.4", "Conf1.5"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_2_daytime, into = c("Conf2.1", "Conf2.2", "Conf2.3", "Conf2.4", "Conf2.5"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_2_nighttime, into = c("Conf2N.1", "Conf2N.2"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_3_PAM, into = c("Conf3.1", "Conf3.2", "Conf3.3", "Conf3.4", "Conf3.5", "Conf3.6", "Conf3.7", "Conf3.8"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_3_Stomach_bycatch, into = c("Conf3ByS.1", "Conf3ByS.2"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_4, into = c("Conf4.1", "Conf4.2", "Conf4.3", "Conf4.4", "Conf4.5", "Conf4.6", "Conf4.7", "Conf4.8"), sep = ",")
diel_full <- separate(data = diel_full, col = Conf_5, into = c("Conf5.1", "Conf5.2", "Conf5.3", "Conf5.4"), sep = ",")

diel_full <- diel_full %>% mutate_all(na_if,"")
diel_full <- diel_full %>% mutate_all(na_if," ")

#use below for all columns
#df <- diel_full %>% make_long(3:36)

#use below for custom column selection
#df <- diel_full %>% make_long(Conf1.1, Conf2.1, Conf3.1, Conf4.1, Conf5.1)

#for all the columns with level 5 data, how does it agree with other data sources?
#df <- diel_full %>% make_long(Conf3.1, Conf3.2, Conf3.3, Conf3.4, Conf3.5)
#df <- diel_full %>% make_long(Conf5.1, Conf5.2, Conf5.3, Conf5.4)
#df <- diel_full %>% make_long(Conf4.1, Conf4.2, Conf4.3, Conf4.4, Conf4.5)
df <- diel_full %>% make_long(Conf3.1, Conf4.1, Conf5.1)

df$node <- tolower(df$node)
df$node <- str_trim(df$node)

df$node <- str_replace_all(df$node, pattern = "reverse-dvm", replacement = "unclear")
df$node <- str_replace_all(df$node, pattern = "dvm/nocturnal", replacement = "nocturnal")
df$node <- str_replace_all(df$node, pattern = "dvm", replacement = "unclear")
df$node <- str_replace_all(df$node, pattern = "tbd", replacement = "unclear")
df$node <- str_replace_all(df$node, pattern = "cathemeral-variable", replacement = "cathemeral")
df$node <- str_replace_all(df$node, pattern = "cathemeral-invariate", replacement = "cathemeral")
df$node <- str_replace_all(df$node, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#calling these nocturnal, diurnal, cathemeral (removing unclear as an option)
df$node <- str_replace_all(df$node, pattern = "unclear-nocturnal", replacement = "nocturnal")
df$node <- str_replace_all(df$node, pattern = "unclear-diurnal", replacement = "diurnal")
df$node <- str_replace_all(df$node, pattern = "unclear-cathemeral", replacement = "cathemeral")

#calling these nocturnal, diurnal, crepuscular, etc may switch to cathemeral later
#idk what to do about this one
df$node <- str_replace_all(df$node, pattern = "weakly-nocturnal/crepuscular", replacement = "nocturnal/crepuscular")
df$node <- str_replace_all(df$node, pattern = "weakly-diurnal", replacement = "diurnal/cathemeral")
df$node <- str_replace_all(df$node, pattern = "weakly-crepuscular", replacement = "crepuscular")
df$node <- str_replace_all(df$node, pattern = "weakly-nocturnal", replacement = "nocturnal/cathemeral")

df$next_node <- tolower(df$next_node)
df$next_node <- str_trim(df$next_node)
df$next_node <- str_replace_all(df$next_node, pattern = "reverse-dvm", replacement = "unclear")
df$next_node <- str_replace_all(df$next_node, pattern = "dvm/nocturnal", replacement = "nocturnal")
df$next_node <- str_replace_all(df$next_node, pattern = "dvm", replacement = "unclear")
df$next_node <- str_replace_all(df$next_node, pattern = "tbd", replacement = "unclear")
df$next_node <- str_replace_all(df$next_node, pattern = "cathemeral-variable", replacement = "cathemeral")
df$next_node <- str_replace_all(df$next_node, pattern = "cathemeral-invariate", replacement = "cathemeral")
df$next_node <- str_replace_all(df$next_node, pattern = "cathemeral/crepuscular", replacement = "crepuscular")

#calling these nocturnal, diurnal, cathemeral (removing unclear as an option)
df$next_node <- str_replace_all(df$next_node, pattern = "unclear-nocturnal", replacement = "nocturnal")
df$next_node <- str_replace_all(df$next_node, pattern = "unclear-diurnal", replacement = "diurnal")
df$next_node <- str_replace_all(df$next_node, pattern = "unclear-cathemeral", replacement = "cathemeral")

#calling these nocturnal, diurnal, crepuscular, etc may switch to cathemeral later
#idk what to do about this one
df$next_node <- str_replace_all(df$next_node, pattern = "weakly-nocturnal/crepuscular", replacement = "nocturnal/crepuscular")
df$next_node <- str_replace_all(df$next_node, pattern = "weakly-diurnal", replacement = "diurnal/cathemeral")
df$next_node <- str_replace_all(df$next_node, pattern = "weakly-crepuscular", replacement = "crepuscular")
df$next_node <- str_replace_all(df$next_node, pattern = "weakly-nocturnal", replacement = "nocturnal/cathemeral")
df$next_node <- str_replace_all(df$next_node, pattern = "weakly-nocturnal", replacement = "nocturnal/cathemeral")


df[df == "unclear"] <- NA
#remove entries where both the starting and ending node are NA
df <- df %>% filter_at(vars(node, next_node), any_vars(!is.na(.)))

#optional: for four states
df$node <- str_replace_all(df$node, pattern = "nocturnal/crepuscular", replacement = "nocturnal")
df$node <- str_replace_all(df$node, pattern = "diurnal/crepuscular", replacement = "diurnal")
df$next_node <- str_replace_all(df$next_node, pattern = "nocturnal/crepuscular", replacement = "nocturnal")
df$next_node <- str_replace_all(df$next_node, pattern = "diurnal/crepuscular", replacement = "diurnal")
df$node <- str_replace_all(df$node, pattern = "diurnal/cathemeral", replacement = "cathemeral")
df$node <- str_replace_all(df$node, pattern = "nocturnal/cathemeral", replacement = "cathemeral")
#not sure whether to call these cathemeral or di/noc
df$next_node <- str_replace_all(df$next_node, pattern = "diurnal/cathemeral", replacement = "cathemeral")
df$next_node <- str_replace_all(df$next_node, pattern = "nocturnal/cathemeral", replacement = "cathemeral")

ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node))) + geom_sankey() + theme_sankey(base_size = 16)



# Section 5: Artiodactyl   ------------------------------------------------
url <- 'https://docs.google.com/spreadsheets/d/1JGC7NZE_S36-IgUWpXBYyl2sgnBHb40DGnwPg2_F40M/edit?gid=562902012#gid=562902012'
diel_full <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
diel_full <- diel_full[!is.na(diel_full$Confidence_primary_source), c(1, 2, 6:13)]
diel_full <- diel_full %>% filter(Confidence_secondary_source != "")

diel_full <- diel_full %>% pivot_longer(cols = c(Confidence_primary_source, Confidence_secondary_source, Confidence_tertiary_source, Confidence_4th_source), names_to = "column", values_to = "values")

diel_full$Confidence_level <- paste(diel_full$column, diel_full$values, sep = "_")

diel_full_1 <- diel_full %>% filter(column == "Confidence_primary_source") %>% pivot_wider(names_from = Confidence_level, values_from = Diel_Pattern_primary) 
diel_full_1 <- diel_full_1[, c(1, 2, 8:12)]
diel_full_2 <- diel_full %>% filter(column == "Confidence_secondary_source") %>% pivot_wider(names_from = Confidence_level, values_from = Diel_Pattern_secondary)
diel_full_2 <- diel_full_2[, c(8:11)]
diel_full_3 <- diel_full %>% filter(column == "Confidence_tertiary_source") %>% pivot_wider(names_from = Confidence_level, values_from = Diel_Pattern_tertiary)
diel_full_3 <- diel_full_3[, c(8:12)]
diel_full_4 <- diel_full %>% filter(column == "Confidence_4th_source") %>% pivot_wider(names_from = Confidence_level, values_from = Diel_pattern_4th_source)
diel_full_4 <- diel_full_4[, c(8:11)]

diel_full <- cbind(diel_full_1, diel_full_2, diel_full_3, diel_full_4)
#rename the columns to the confidence level, R will add .1 for every duplicate so they will end up with unique identifiers
colnames(diel_full) <- c("Species_name", "Family", "Conf1", "Conf2", "Conf3", "Conf4", "Conf5", "Conf1", "Conf3", "Conf2", "Conf2", "tertiary_source_NA", "Conf3", "Conf4", "Conf2", "Conf1", "fourth_source_NA", "Conf3", "Conf2", "Conf4")
#drop the two columns with only NA values ("tertiary_source_NA", "fourth_source_NA")
diel_full <- diel_full[, -c(12, 17)]

#convert into long format
diel_full_long <- diel_full %>% pivot_longer(cols = c(3:18), names_to = "column", values_to = "value")

#replace any unconventional strings
unique(diel_full_long$value)
diel_full_long$value <- tolower(diel_full_long$value)
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "weakly-diurnal", replacement = "diurnal/cathemeral")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "weakly-nocturnal", replacement = "nocturnal/cathemeral")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "unclear-diurnal", replacement = "diurnal")
#check what these are later
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "diurnal*", replacement = "diurnal")
diel_full_long$value <- str_replace_all(diel_full_long$value, pattern = "nocturnal*", replacement = "nocturnal")

diel_full_long[diel_full_long == "unclear"] <- NA
diel_full_long[diel_full_long == ""] <- NA
diel_full_long <- diel_full_long[!is.na(diel_full_long$value),]

species_list <- table(diel_full_long$Species_name)
#should be 116 species
species_list <- names(species_list[species_list > 1])

#function Max wrote for comparing entries
compTwo <- function(comp1 = "comp1", comp2 = "comp2") {
  
  if(any(is.na(c(comp1, comp2)))) {
    return(NA)
  } else {
    #splits any entries with a backslash into two components (ie nocturnal/crepuscular into nocturnal and crepuscular)
    comp1 <- str_split(comp1, "/")[[1]]
    comp2 <- str_split(comp2, "/")[[1]]
    #then compares if any of the components match
    if(any(comp1 %in% comp2)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
}


#apply this function across all species with multiple entries
output <- lapply(species_list, function(species) {
  
  #filter for one species at a time
  df <- diel_full_long[diel_full_long$Species_name == species,]
  #rename the column names to be unique for every entry for this species (ie for multiple column 2 entries column 2.1, 2.2 etc)
  df$column <- make.unique(df$column)
  
  #converts the dataframe so it compares every entry with each other (ie for A,B,C A-A, A-B, A-C, B-A, B-B, B-C, etc)
  df_lists_comb <- expand(df, nesting(var = column, vector = value), nesting(var2 = column, vector2 = value), .name_repair = "universal")
  
  #??? idk
  df_lists_comb <- df_lists_comb %>% filter(var != var2) %>% arrange(var, var2) %>% mutate(vars = paste0(var, ".", var2)) %>% select(contains("var"), everything())
  
  #evaluates the activity patterns for each of these sources and returns if they agree or not (TRUE or FALSE)
  comparisons <- df_lists_comb %>% group_by(vars) %>% mutate(comp = compTwo(comp1 = vector, comp2 = vector2))
  #manipulate the strings for both variable names to revert them back to the original name (back to column 2 from col 2.1)
  comparisons$var <- str_sub(comparisons$var, start = 1, end = 5)
  comparisons$var2 <- str_sub(comparisons$var2, start = 1, end = 5)
  
  #create a column returning the comparison being made (ie col2-col2, col1-col2, etc)
  comparisons$var_final <- paste(comparisons$var, comparisons$var2, sep = "-")
  
  #return just the comparison result column (TRUE or FALSE match) and the comparison being made (ie col1 vs col1)
  return(comparisons[,c("comp","var_final")])
})

#combine this list of results 
output <- Reduce(rbind, output)

table <- table(output$var_final)
# prop.table(table, margin = 1)
table2 <- as.data.frame(prop.table(table(output$var_final, output$comp), margin = 1))
table2$Comp1 <- sapply(str_split(table2$Var1, "-"), `[`, 1)
table2$Comp2 <- sapply(str_split(table2$Var1, "-"), `[`, 2)
table2 <- table2[table2$Var2 == TRUE,]
table2$count <- table
plot_freq <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = round(Freq, digits = 2))) + geom_tile() + geom_text() + scale_fill_viridis(limits = c(0,1))
plot_freq
plot_count <- ggplot(table2, aes(x = Comp1, y = Comp2, fill = Freq, label = count)) + geom_tile() + geom_text() + scale_fill_viridis(limits = c(0,1))
plot_count
