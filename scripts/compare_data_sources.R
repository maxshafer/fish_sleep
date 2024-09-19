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

setwd(here())

### LOAD IN OTHER DATA ### 

# Fetch the taxonomic levels from fishbase for all of the species and make it into a dataframe
# available_releases()
# [1] "23.01" "23.05" "21.06" "19.04"

fishbase_df <- load_taxa(collect = T, version = "21.06")
fishbase_df <- as.data.frame(fishbase_df)

url <- 'https://docs.google.com/spreadsheets/d/18aNqHT73hX06cGRlf6oj7Y4TVKf6jd_Q5ojNIxm2rys/edit#gid=0'
sleepy_fish <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
sleepy_fish$Diel_Pattern <- tolower(sleepy_fish$Diel_Pattern)

### I want to make new columns for every possible pairwise comparison of columns 9-22 (or possible columns that might exist)
### I also need to figure out how to allow "Crepuscular" and "Nocturnal" to match with "Crepuscular/Nocturnal", etc, and from both directions

data <- sleepy_fish[!(sleepy_fish$X5. == "" & sleepy_fish$X5 == "" & sleepy_fish$X4. == "" & sleepy_fish$X4.. == "" & sleepy_fish$X4 == "" & sleepy_fish$X3.. == "" & sleepy_fish$X3. == "" & sleepy_fish$X3 == "" & sleepy_fish$X2. == "" & sleepy_fish$X2.. == "" & sleepy_fish$X2 == "" & sleepy_fish$X1. == "" & sleepy_fish$X1.. == "" & sleepy_fish$X1 == ""),]

data <- data[,c("Species_name","Diel_Pattern","X5", "X5.", "X5..","X4", "X4.", "X4..","X3", "X3.", "X3..","X2", "X2.", "X2..","X1", "X1.", "X1..")]

# OK this works and returns what I had planned
data[data == ""] <- NA

## I could likely make a new function that compares two values, and outputs TRUE/FALSE based on my own logic rules
## based on many ifelse statements

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

## These are the comparisons that need to be run, like above
## They should become the new columns
comparisons <- expand.grid(colnames(data)[3:ncol(data)], colnames(data)[3:ncol(data)])
comparisons <- comparisons[!(comparisons$Var1 == comparisons$Var2),]


## This seems to work
test <- apply(comparisons, 1, function(x) {
  
  out <- apply(data, 1, function(y) compTwo(comp1 = y[as.character(x[1])], comp2 = y[as.character(x[2])]))
  
  #out <- data[,as.character(x[1])] == data[,as.character(x[2])]
  return(out)
  
})

## This works, but doesn't take into account the different #s involved (just averages averages)

df <- data.frame(names_1 = comparisons$Var1, names_2 = comparisons$Var2, overlap = apply(test, 2, function(x) mean(x, na.rm = TRUE)))

df$trues <- colSums(test*1, na.rm = T)
df$falses <- colSums(!(test)*1, na.rm = T)
## Maybe remove "." and then combine by same rows?

df$names_1 <- str_replace_all(df$names_1, "[.]", "")

df$names_2 <- str_replace_all(df$names_2, "[.]", "")

df$names <- paste(df$names_1, df$names_2, sep = "_")

df <- df %>% group_by(names) %>% mutate(agree = sum(trues, na.rm = T), disagree = sum(falses, na.rm = T))
df$average <- df$agree/(df$agree + df$disagree)

df2 <- unique(df[,c(1,2,6,9)])

ggplot(df2, aes(x = names_1, y = names_2, fill = average, label = average)) + geom_tile() + scale_fill_continuous(limits=c(0, 1)) + geom_text()


## Count the number of sources (and histogram of species per source)

sources <- c(sleepy_fish[,grep("Source", colnames(sleepy_fish))])
sources <- unlist(sources)
sources <- sources[sources != ""]

hist(log(table(sources)))
