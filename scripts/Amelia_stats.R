library(phangorn)
library(stringr)
library(here)
library(gsheet)
library(dplyr)
library(phytools)
library(geiger)
library(dplyr)
library(readxl)
library(tidyr)
library(parsedate)
library(lubridate)
library(ggplot2)
library(forcats)
library(ggforce)

#we want to see if the difference between the mean likelihood of one model is greater than another

#ie is the ARD model statistically more likely than the bridge_only model

# Compute the analysis of variance
res.aov <- aov(weight ~ group, data = my_data)
# Summary of the analysis
summary(res.aov)