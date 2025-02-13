#trying to get the pantheria database to work

# Install pak
install.packages("pak")
library(pak)

#with the pkg_install function from pak install the traitdata database from github
pkg_install("RS-eco/traitdata")

#load library
library(traitdata)
data(pantheria)



