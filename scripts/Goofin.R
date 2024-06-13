#how well does data from different sources (categories) agree with each other

#import fish database
fish_df <- read.csv(here("sleepy_fish_database_local.csv"))

#look for exact matches
df$match <- "Idk"
  
for(i in 1:length(fish_df$Species_name)){
  if(fish_df[i, "column_1"] == fish_df[i, "Column_2"]){
    fish_df[i, "match"] <- "Yes"
  } else if(fish_df[i, "column_1"] %in% c("crepuscular/diurnal", "crepuscular", "diurnal") & fish_df[i, "column_2"] %in% c("crepuscular/diurnal", "crepuscular", "diurnal")){
    fish_df[i, "match"] <- "Approximate"
  } else if(fish_df[i, "column_1"] %in% c("crepuscular/nocturnal", "crepuscular", "nocturnal") & fish_df[i, "column_2"] %in% c("crepuscular/nocturnal", "crepuscular", "nocturnal")){
    fish_df[i, "match"] <- "Approximate"
  } else {
    fish_df[i, "match"] <- "No"
  }
}

#allow for approximate matches
fish <- c("crepuscular/diurnal", "crepuscular", "diurnal")
a <- "crepuscular"
b <- "diurnal"

a %in% fish & b %in% fish
