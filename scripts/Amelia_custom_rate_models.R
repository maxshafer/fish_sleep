# Section 1.0 Cetacean cathemeral dead end-------------------------------

#create a custom matrix to test if cathemerality is a dead end 
# #ie rate from cathemeral into other diel patterns is zero

# generic_ratemat <- getStateMat4Dat(cor_model_ard3$data)
# test_ratemat <- generic_ratemat$rate.mat
# test_ratemat
# 
# custom_rate_matrix <- dropStateMatPars(test_ratemat, c(1,6))
# custom_rate_matrix
# 
# custom_ard <- corHMM(phy = cor_model_ard3$phy, data = cor_model_ard3$data, rate.cat = 1, rate.mat = custom_rate_matrix)

#matrix2 <- matrix(c(NA, 1L, 2L, NA, NA, 3L, NA, 4L, NA), nrow =3, ncol = 3)
#ace_deadend_ard <- ace(trait.vector3, trpy_n3, model = matrix2, type = "discrete")

#likelihoods <- rbind(likelihoods, c("ace_deadend_ard", ace_deadend_ard$loglik, "diurnal, nocturnal, and cathemeral as a dead-end"))


# Section 2.0 Cetacean cathemeral bridge ----------------------------------
head(cetaceans_full)


# Cetacean discrete traits  -----------------------------------------------
###discrete traits to model with
discrete_cet <- read.csv("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\resource.csv")
more_cet_traits <- read_xlsx("C:\\Users\\ameli\\OneDrive\\Documents\\R_projects\\cetacean_discrete_traits\\Coombs_et_al_2022.xlsx")
more_cet_traits <- more_cet_traits[!(is.na(more_cet_traits$No.)),]
