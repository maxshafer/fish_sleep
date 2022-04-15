library(ggplot2)
library(ggalt)
library(rworldmap)
library(rnaturalearth)

df <- fishbase_distribution[fishbaase_distribution$Species %in% "Siganus fuscescens",c("LatDeg", "N_S", "LongDeg", "E_W")]
df$LatDeg <- ifelse(df$N_S == "S", df$LatDeg*-1, df$LatDeg)
df$LongDeg <- ifelse(df$N_S == "S", df$LongDeg*-1, df$LongDeg)

ggplot(data = df) + geom_polygon(data = df, aes(y = LatDeg, x = LongDeg))




df2 <- fishbase_ecosystem[fishbase_ecosystem$Species %in% "Siganus fuscescens",c("NorthernLat", "NrangeNS", "SouthernLat", "SrangeNS", "WesternLat", "WrangeEW", "EasternLat", "ErangeEW")]





## Combine fishbase_ecosystem and fishbase_distrution?


test <- fishbase_ecosystem[1:10000,]
test <- fishbase_ecosystem[fishbase_ecosystem$Species %in% resolved_names$unique_name,]

df <- test[,c("Species", "Salinity", "NorthernLat", "NrangeNS", "SouthernLat", "SrangeNS", "WesternLat", "WrangeEW", "EasternLat", "ErangeEW")]
df <- df[!(is.na(df$NorthernLat)),]
colnames(df) <- c("Species", "Salinity", "Lat1", "Lat1_sign", "Lat2", "Lat2_sign", "Long1", "Long1_sign", "Long2", "Long2_sign")

df$Lat1_sign <- ifelse(df$Lat1_sign %in% c("n", "N"), 1, -1)
df$Lat2_sign <- ifelse(df$Lat2_sign %in% c("n", "N"), 1, -1)
df$Long1_sign <- ifelse(df$Long1_sign %in% c("w", "W"), -1, 1)
df$Long2_sign <- ifelse(df$Long2_sign %in% c("w", "W"), -1, 1)

df$Lat1 <- as.numeric(df$Lat1)*df$Lat1_sign
df$Lat2 <- as.numeric(df$Lat2)*df$Lat2_sign
df$Long1 <- as.numeric(df$Long1)*df$Long1_sign
df$Long2 <- as.numeric(df$Long2)*df$Long2_sign

df2 <- data.frame(Species = c(df$Species, df$Species), Salinity = c(df$Salinity, df$Salinity), Lat = c(df$Lat1, df$Lat2), Long = c(df$Long1, df$Long2))

df2$Diel <- resolved_names$diel2[match(df2$Species, resolved_names$unique_name)]
df2$Family <- resolved_names$family[match(df2$Species, resolved_names$unique_name)]

test.plot <- ggplot(df2, aes(y = Lat, x = Long)) + geom_point(aes(color = Diel), show.legend = TRUE) + scale_colour_manual(values = c("yellow", "red", "blue", "black", "black")) + ylim(c(-90,90)) + xlim(c(-180, 180)) + theme_classic() #+ geom_encircle(aes(fill = Species), s_shape = 1, expand = 0, alpha = 0.1, color = "grey20", show.legend = FALSE)

world <- ne_countries(scale = "medium", returnclass = "sf")

world <- ggplot(data = world) +  geom_sf() +  labs( x = "Longitude", y = "Latitude") +  ggtitle("World map", subtitle = paste0("(", length(unique(world$admin)), " countries)"))

salinity.plot <- world + geom_point(data = df2, aes(y = Lat, x = Long, color = Salinity), show.legend = TRUE, size = 0.5) + scale_colour_manual(values = c("yellow", "black", "blue")) + ylim(c(-90,90)) + xlim(c(-180, 180))

diel.plot <- world + geom_point(data = df2, aes(y = Lat, x = Long, color = Diel), show.legend = TRUE, size = 0.5) + scale_colour_manual(values = c("yellow", "red", "blue", "black")) + ylim(c(-90,90)) + xlim(c(-180, 180))

# family.plot <- world + geom_point(data = df2, aes(y = Lat, x = Long, color = Family), show.legend = FALSE) + ylim(c(-90,90)) + xlim(c(-180, 180))

# species.plot <- world + geom_point(data = df2, aes(y = Lat, x = Long, color = Species), show.legend = FALSE) + ylim(c(-90,90)) + xlim(c(-180, 180))

world + salinity.plot + diel.plot + plot_layout(nrow = 3)




ggplot(df2, aes(y = Lat, x = Long)) + stat_density2d(geom = "polygon", aes(fill = Diel), alpha = 0.1)



+ geom_point(aes(color = Diel), show.legend = TRUE) + scale_colour_manual(values = c("yellow", "red", "blue", "black", "black")) + ylim(c(-90,90)) + xlim(c(-180, 180)) + theme_classic()








