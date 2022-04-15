library(rotl)
library(ggtree)
library(stringr)
library(scales)
library(gsheet)
library(ape)
library(patchwork)
library(ggpubr)
library(dplyr)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Sleepy_fishes/")

## Make a script to extract DOIs or PMIDs, to make a list of all PMIDs that we already have in the list (for filtering the machine-derived PMIDs)
# Load in the file, to get the urls, which contain either the PMID, or the DOI

url <- 'https://docs.google.com/spreadsheets/d/18aNqHT73hX06cGRlf6oj7Y4TVKf6jd_Q5ojNIxm2rys/edit?usp=sharing'

sleepy_fish <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)


## OK so the IDs or doi's come in a few flavours:
# Can be a #######, with "/"'s on either side, can also be PIM#####, but both have "ncbi" in the url
# Need to search all columns, maybe an apply on just that subset? go by row, run this for each entry (won't find if it doesn't exist)

subset <- sleepy_fish[,13:21]

# This works well, and extracts 201 ids for 201 instances of "ncbi" in the google sheet
# It also returns a structure that includes empty spots, so I can easily match it with other approachs to find DOIs

PMIDs <- apply(subset, 1, function(x) {
  line <- x
  # if its just a number, it looks like it doesn't have "pmc" or "PMC", and the number is the last entry separate by "/"
  index.ncbi <- grepl("ncbi", line)
  line <- line[index.ncbi]
  index.pmc <- grepl("PMC|pmc", line)
  index.nonpmc <- !(index.pmc)
  line.str <- lapply(line, function(x) unlist(strsplit(x, split = "/"))) # this gives me a list of vectors (str split)
  # extract the PMIDs for those that match the index
  line.str.pmc <- unlist(line.str[index.pmc])
  ids.pmc <- line.str.pmc[grep("PMC", line.str.pmc)]
  line.str.nonpmc <- line.str[index.nonpmc]
  ids.non.pmc <- lapply(line.str.nonpmc, function(x) x[length(x)]) # this takes the last entry, but sometimes it's not the number
  ids <- unname(c(unlist(ids.pmc), unlist(ids.non.pmc)))
  return(ids)
})



## OK can I find doi's now? What do they look like?
## Below works well at pulling dois, but it also gets a few extra that aren't dois - looks pretty good otherwise

DOIs <- apply(subset, 1, function(x) {
  line <- x
  
  index.doi <- !(grepl("ncbi", line))
  line <- line[index.doi]
  index.doi <- grepl("10.", line)
  line <- line[index.doi]
  
  # looks like it is always the last two entries separated by slashs, with anything past a # extra and can be removed
  # two layer strsplit?
  
  line.str <- lapply(line, function(x) unlist(strsplit(x, split = "/")))
  #line.str <- line.str[!(unlist(lapply(line.str, function(x) identical(x, character(0)))))]
  
  dois <- unlist(lapply(line.str, function(x) paste(x[length(x)-1], x[length(x)], sep = "/")))
  if(any(grep("#", dois))) {
    dois <- strsplit(dois, split = "#")
    dois <- unlist(lapply(dois, function(x) x[1]))
  }
  if(any(grep("pdf", dois))) {
    dois <- strsplit(dois, split = "pdf")
    dois <- unlist(lapply(dois, function(x) x[1]))
  }
  dois <- unname(dois)
  return(dois)
})


# How many of the fish entries do I have either a PMID or a doi for? Or how many are missing one or the other
# Actually a pretty surprising number, though how many of the DOIs aren't really DOIs, I do not know

test <- lapply(DOIs, function(x) is.null(x))
length(test[!(unlist(test))])
# [1] 936
test <- lapply(PMIDs, function(x) is.null(x))
length(test[!(unlist(test))])
# [1] 146

## I can make a two vectors, and then when I pull the pubmed results, I can also return the DOIs, and check the PMIDs and the dois against this list?

## dois can be made activatable by pasting "https://doi.org/" in front of them

saveRDS(PMIDs, file = "sleepy_fish_PMIDs.rds")
saveRDS(DOIs, file = "sleepy_fish_DOIs.rds")

