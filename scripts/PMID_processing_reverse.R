library(RISmed)
library(rfishbase)
library(googlesheets4)

setwd("/Volumes/BZ/Scientific Data/RG-AS04-Data01/Euteleostomi_deil_patterns/")

PMID_list <- readRDS(file = "sleepy_fish_PMIDs.rds")



fishbase_df <- load_taxa(collect = T, version = "19.04")
fishbase_df <- as.data.frame(fishbase_df)

`%rin%` = function (pattern, list) {
  vapply(pattern, function (p) any(grepl(p, list)), logical(1L), USE.NAMES = FALSE)
}



search_topics <- list("diel", "diurnal", "nocturnal", "crepuscular", "cathemeral")
search_query <- lapply(search_topics, function(x) EUtilsSummary(x, type="esearch", db="pubmed", retmax = 99999))
counts <- lapply(search_query, function(x) x@count)
record_df <- list()
for (i in 1:length(search_query)) {
  if (counts[[i]] < 1000){
    records <- EUtilsGet(search_query[[i]], type = "efetch", db = "pubmed")
    record_df[[i]] <- data.frame('PMID' = PMID(records), 'DOI' = DOI(records), 'Title' = ArticleTitle(records), 'Abstract' = AbstractText(records), 'Year' = YearPubmed(records), 'Journal' = Title(records))
  }
  if (counts[[i]] > 1000) {
    bins <- c(1, (1:floor(counts[[i]]/1000))*1000, counts[[i]])
    records <- list()
    record_df_2 <- list()
    for (j in 1:(length(bins)-1)) {
      records[[j]] <- EUtilsGet(search_query[[i]]@PMID[bins[j]:bins[j+1]], type = "efetch", db = "pubmed")
      record_df_2[[j]] <- data.frame('PMID' = PMID(records[[j]]), 'DOI' = DOI(records[[j]]), 'Title' = ArticleTitle(records[[j]]), 'Abstract' = AbstractText(records[[j]]), 'Year' = YearPubmed(records[[j]]), 'Journal' = Title(records[[j]]))
    }
    record_df[[i]] <- Reduce(rbind, record_df_2)
  }
}

search_topic <- "diel"
search_query <- EUtilsSummary(search_topic, type = "esearch", db = "pubmed", retmax = 99999)
records <- EUtilsGet(search_query@PMID[2000:3000], type = "efetch", db = "pubmed")
record_df <- data.frame('PMID' = PMID(records), 'DOI' = DOI(records), 'Title' = ArticleTitle(records), 'Abstract' = AbstractText(records), 'Year' = YearPubmed(records), 'Journal' = Title(records))


## Wowow this works! and doesn't take forever!
## So just need to run the above on small enough subsets that it works (diel is ok, at 10k, but not diurnal at ~30k)

test <- str_c(unique(fishbase_df$Species), collapse = "|")

test2 <- record_df %>% filter(str_detect(record_df$Title,test))
