
library("sp")
library("raster")
library("maptools")
library("rgdal")
library("dismo")
library(geodata)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(gridExtra)
#install.packages("viridis")
library(viridis)
#install.packages("scico")
library(scico)

#These are observations I downloaded off iNaturalist
library(readr)
GlobalRhizoObs <- read_csv("~/GlobalRhizostomeDistribution/GlobalRhizoObs.csv")
#View(GlobalRhizoObs)
summary(GlobalRhizoObs)

#This will remove any observations with broken coordinate data
obs_data <- GlobalRhizoObs[!is.na(GlobalRhizoObs$latitude), ]
obs_data <- obs_data[!is.na(obs_data$longitude), ]
summary(obs_data)

#We can visualize the data with pretty graphs here
#Compare longitudes, lattitudes, etc.

bar <- ggplot(obs_data, aes(x = fct_infreq(taxon_genus_name), fill = taxon_suborder_name))+
  geom_bar()+
  scale_fill_manual(values = c("#1E88E5","#D81B60"))+
  scale_y_log10(name = "Observation #")+
  xlab("Genus")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_text())
bar
box <- ggplot(obs_data, aes(x = taxon_suborder_name, y = latitude,  fill = taxon_suborder_name))+
  geom_violin(alpha = 0.5)+
  geom_boxplot(alpha = 0.5, size = 1.5, )+
  scale_fill_manual(values = c("#1E88E5","#D81B60"))+
  scale_color_manual(values = c("#1E88E5","#D81B60"))+
  xlab("Order")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90), 
        legend.position = "none")
box
box2 <- ggplot(obs_data, aes(x = fct_infreq(taxon_genus_name), y = latitude, fill = taxon_suborder_name, color = taxon_suborder_name))+
  geom_boxplot(alpha = 0.5, size = 1.5)+
  scale_fill_manual(values = c("#1E88E5","#D81B60"))+
  scale_color_manual(values = c("#1E88E5","#D81B60"))+
  theme_classic()+
  xlab("Genus")+
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.y = element_blank(), legend.position = "none")
box2

#This allows us to arrange all the graphs on the same panel
library(gridExtra)
boxes <- grid.arrange(box, box2, nrow = 1)
boxes
grid.arrange(bar, boxes, nrow = 2)

#Raw visualization of distributions before any filtering
ggplot(data = obs_data, aes(longitude, latitude, color = taxon_genus_name))+
         geom_point(alpha = 0.3)+
  facet_grid(rows = vars(taxon_suborder_name))+
  theme_classic()
#Looks like there are some unique distributions

##Data curation for extrapolation curve
lowresobs<-obs_data
lowresobs$latitude <- round(lowresobs$latitude, digits = -1)
head(lowresobs)
lowresobs$longitude <- round(lowresobs$longitude, digits = -1)

abundbylowrescoordinate <- lowresobs %>%
  group_by(latitude, longitude, taxon_suborder_name, taxon_genus_name) %>%
  summarize(observations = n())
head(abundbylowrescoordinate)

genus.wide.coordinate <- pivot_wider(abundbylowrescoordinate, names_from = taxon_genus_name, values_from = observations) 
head(genus.wide.coordinate)
class(genus.wide.coordinate)
genus.wide.coordinate <- as.data.frame(genus.wide.coordinate)

genus.wide.coordinate <- dplyr::select(genus.wide.coordinate, -c(latitude, longitude, taxon_suborder_name))
head(genus.wide.coordinate)

genus.wide.coordinate[is.na(genus.wide.coordinate)] <- 0
head(genus.wide.coordinate)

#install.packages("iNEXT")
library(iNEXT)
out <- iNEXT(genus.wide.coordinate, q = 0, datatype = "abundance")
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 22
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
summary(out)

#Visualization rarefaction curves
ggiNEXT(out)+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors)+
  ylab("Coordinate diversity")+
  theme_classic()
#While this is helpful, just looking at the curves makes it a bit arbitrary

### Use SC to look at saturation?
Extrap <- out[["iNextEst"]][["coverage_based"]]
head(Extrap)

ggplot(Extrap, aes(Assemblage, SC, fill = Method, color = Method))+
  geom_hline(yintercept = 0.95, size = 2, color = "#56b4e9ff")+
  geom_boxplot(alpha = 0.5, size = 1, outlier.alpha = 0)+
  scale_fill_manual(values = c("#5e5e5eff","#3a3a3aff","#000000ff"))+
  scale_color_manual(values = c("#5e5e5eff","#3a3a3aff","#000000ff"))+
  theme_classic()+
  xlab("Genus")+
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.y = element_blank())

#If we set our desired sampling as >95% after rarefaction,
#Suggests we have sufficient sampling for:
#Phyllorhiza, Rhizostoma, Cassiopea, Catostylus, 
#Rhopilema, Cotylorhiza, Pseudorhiza, Stomolophus

#Let's make genus-specific subsets
Cassiopea <- subset(obs_data, taxon_genus_name == c("Cassiopea"))
Catostylus <- subset(obs_data, taxon_genus_name == c("Catostylus"))
Cotylorhiza <- subset(obs_data, taxon_genus_name == c("Cotylorhiza"))
Eupilema <- subset(obs_data, taxon_genus_name == c("Eupilema"))
Lobonema <- subset(obs_data, taxon_genus_name == c("Lobonema"))
Lychnorhiza <- subset(obs_data, taxon_genus_name == c("Lychnorhiza"))
Pseudorhiza <- subset(obs_data, taxon_genus_name == c("Pseudorhiza"))
Rhizostoma <- subset(obs_data, taxon_genus_name == c("Rhizostoma"))
Rhopilema <- subset(obs_data, taxon_genus_name == c("Rhopilema"))
Stomolophus <- subset(obs_data, taxon_genus_name == c("Stomolophus"))

#Now we have unique datesets for each genus,
#But I would like to make a filtered observation database too
filtered_obs_data <- rbind(Cassiopea, Catostylus)
filtered_obs_data <- rbind(filtered_obs_data, Cotylorhiza)
filtered_obs_data <- rbind(filtered_obs_data, Eupilema)
filtered_obs_data <- rbind(filtered_obs_data, Lobonema)
filtered_obs_data <- rbind(filtered_obs_data, Lychnorhiza)
filtered_obs_data <- rbind(filtered_obs_data, Pseudorhiza)
filtered_obs_data <- rbind(filtered_obs_data, Rhizostoma)
filtered_obs_data <- rbind(filtered_obs_data, Rhopilema)
filtered_obs_data <- rbind(filtered_obs_data, Stomolophus)
head(filtered_obs_data)
#View(filtered_obs_data)

#Let's look at what we have left
box3 <- ggplot(filtered_obs_data, aes(x = fct_infreq(taxon_genus_name), y = latitude, fill = taxon_suborder_name, color = taxon_suborder_name))+
  geom_boxplot(alpha = 0.5, size = 1.5)+
  geom_violin(alpha = 0.5, linewidth = 1)+
  #geom_jitter(alpha = 0.3, color = "black")+
  scale_fill_manual(values = c("#1E88E5","#D81B60"))+
  scale_color_manual(values = c("#1E88E5","#D81B60"))+
  theme_classic()+
  xlab("Genus")+
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")
box4 <- ggplot(filtered_obs_data, aes(x = fct_infreq(taxon_genus_name), y = longitude, fill = taxon_suborder_name, color = taxon_suborder_name))+
  geom_boxplot(alpha = 0.5, size = 1.5)+
  geom_violin(alpha = 0.5, linewidth = 1)+
  #geom_jitter(alpha = 0.3, color = "black")+
  scale_fill_manual(values = c("#1E88E5","#D81B60"))+
  scale_color_manual(values = c("#1E88E5","#D81B60"))+
  theme_classic()+
  xlab("Genus")+
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.y = element_blank(), legend.position = "none")+
  coord_flip()
grid.arrange(box3, box4)

#Which distributions are the same? Different?
#install.packages("FSA")
library(FSA)
#install.packages("rcompanion")
library(rcompanion)

#Suborder latitudinal comparison
kruskal.test(latitude ~ taxon_suborder_name, data = filtered_obs_data)
#Genus latitudinal comparison
kruskal.test(latitude ~ taxon_genus_name, data = filtered_obs_data)
#Post-hoc comparison of genera
latstats <- dunnTest(latitude ~ taxon_genus_name, data = filtered_obs_data, method = "bh")
latstats
latstats = latstats$res
cldList(comparison = latstats$Comparison, p.value = latstats$P.adj, threshold = 0.05)
#There is some clumping, but mostly genus specific

#Let's check the longitude for genera
kruskal.test(longitude ~ taxon_genus_name, data = filtered_obs_data)
longstats <- dunnTest(longitude ~ taxon_genus_name, data = filtered_obs_data, method = "bh")
longstats
longstats = longstats$res
cldList(comparison = longstats$Comparison, p.value = longstats$P.adj, threshold = 0.05)
#Seemingly sporadic and probably driven by geographic barriers

#I'd like to see where the residual rhizostomes are
OtherRhizo <- subset(filtered_obs_data, taxon_genus_name != c("Cassiopea")) 
OtherRhizo <- subset(OtherRhizo, taxon_genus_name != c("Catostylus")) 
OtherRhizo <- subset(OtherRhizo, taxon_genus_name != c("Cotylorhiza"))
OtherRhizo <- subset(OtherRhizo, taxon_genus_name != c("Eupilema")) 
OtherRhizo <- subset(OtherRhizo, taxon_genus_name != c("Lobonema"))
OtherRhizo <- subset(OtherRhizo, taxon_genus_name != c("Lychnorhiza")) 
OtherRhizo <- subset(OtherRhizo, taxon_genus_name != c("Pseudorhiza")) 
OtherRhizo <- subset(OtherRhizo, taxon_genus_name != c("Rhizostoma")) 
OtherRhizo <- subset(OtherRhizo, taxon_genus_name != c("Rhopilema")) 
OtherRhizo <- subset(OtherRhizo, taxon_genus_name != c("Stomolophus"))

#Now we can start making our map. Based on all observations, this sets our map boundary
max_lat <- ceiling(max(obs_data$latitude))
min_lat <- floor(min(obs_data$latitude))
max_lon <- ceiling(max(obs_data$longitude))
min_lon <- floor(min(obs_data$longitude))
geographic_extent <- extent(x = c(min_lon, max_lon, min_lat, max_lat))

# Load the data to use for our base map
data(wrld_simpl)
#install.packages("scales")
library(scales)

#Now you can plot our groups of interest like this
plot(wrld_simpl, 
     xlim = c(min_lon, max_lon),
     ylim = c(-60, 80),
     axes = TRUE,
     col = "grey95")
points(x =OtherRhizo$longitude, 
       y =OtherRhizo$latitude, 
       main = "OtherDactyl",
       col = alpha("black", 0.1),
       #col = "#004D40", 
       pch = 20, 
       cex = 2)
points(x = Rhopilema$longitude, 
       y = Rhopilema$latitude,
       col = alpha("#D55E00", 0.01),
       pch = 20, 
       cex = 2)
points(x = Stomolophus$longitude, 
       y = Stomolophus$latitude,
       col = alpha("#0072B2", 0.01), 
       pch = 20, 
       cex = 2)
points(x = Rhopilema$longitude, 
       y = Rhopilema$latitude,
       col = alpha("#D55E00", 0.01),
       pch = 20, 
       cex = 2)
points(x = Stomolophus$longitude, 
       y = Stomolophus$latitude,
       col = alpha("#0072B2", 0.02), 
       pch = 20, 
       cex = 2)
points(x = Rhopilema$longitude, 
       y = Rhopilema$latitude,
       col = alpha("#D55E00", 0.01),
       pch = 20, 
       cex = 2)
points(x = Cassiopea$longitude, 
       y = Cassiopea$latitude,
       col = alpha("#D81B60", 0.06),
       pch = 20, 
       cex = 2)
points(x = Catostylus$longitude, 
       y = Catostylus$latitude,
       col = alpha("#009E73", 0.03), 
       pch = 20, 
       cex = 2)
points(x = Eupilema$longitude, 
       y = Eupilema$latitude,
       col = alpha("#1E88E5", 0.07),
       pch = 20, 
       cex = 2)
points(x = Lobonema$longitude, 
       y = Lobonema$latitude,
       col = alpha("#004D40", 0.3),
       pch = 20, 
       cex = 2)
points(x = Lychnorhiza$longitude, 
       y = Lychnorhiza$latitude,
       col = alpha("#CC79A7", 0.1),
       pch = 20, 
       cex = 2)
points(x = Pseudorhiza$longitude, 
       y = Pseudorhiza$latitude,
       col = alpha("#E69F00", 0.07),
       pch = 20, 
       cex = 2)
points(x = Rhizostoma$longitude, 
       y = Rhizostoma$latitude, 
       main = "Genera with saturated distribution reports",
       col = alpha("#56B4E9", 0.03),
       #col = "#004D40", 
       pch = 20, 
       cex = 2)
points(x = Cotylorhiza$longitude, 
       y = Cotylorhiza$latitude,
       col = alpha("#F0E442", 0.03), 
       pch = 20, 
       cex = 2)

legend("bottomleft", legend=c("Cassiopea", "Catostylus", "Cotylorhiza", "Eupilema", "Lobonema", "Lychnorhiza", 
                              "Pseudorhiza", "Rhizostoma", "Rhopilema", "Stomolophus", "Other Rhizostomae"),
       fill=c("#D81B60", "#009E73", "#F0E442", "#1E88E5", "#004D40", "#CC79A7", "#E69F00", "#56B4E9","#D55E00","#0072B2","#000000"),
        cex=1)
box()

#Latitudinal and longitudinal revisualization as density curves
#install.packages("ggridges")
library(ggridges)
latdensity <- ggplot(filtered_obs_data, aes(latitude, taxon_genus_name, fill = taxon_genus_name, color = taxon_genus_name))+
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_manual(values = c("#D81B60", "#009E73", "#F0E442", "#1E88E5", "#004D40", "#CC79A7", "#E69F00", "#56B4E9","#D55E00","#0072B2"))+
  scale_color_manual(values = c("#D81B60", "#009E73", "#F0E442", "#1E88E5", "#004D40", "#CC79A7", "#E69F00", "#56B4E9","#D55E00","#0072B2"))+
  #scale_x_continuous(limits = c(-60,80))+
  theme_minimal()+
  theme(legend.position = "bottom")
longdensity <- ggplot(filtered_obs_data, aes(longitude, taxon_genus_name, fill = taxon_genus_name, color = taxon_genus_name))+
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_manual(values = c("#D81B60", "#009E73", "#F0E442", "#1E88E5", "#004D40", "#CC79A7", "#E69F00", "#56B4E9","#D55E00","#0072B2"))+
  scale_color_manual(values = c("#D81B60", "#009E73", "#F0E442", "#1E88E5", "#004D40", "#CC79A7", "#E69F00", "#56B4E9","#D55E00","#0072B2"))+
  #scale_x_continuous(limits = c(-60,80))+
  theme_minimal()+
  theme(legend.position = "bottom")
grid.arrange(latdensity, longdensity)


###
#Let's integrate temperature data
###


#Combine all temperature data into a single database
#Only keep 0m, 5m, 10m, 20m
#install.packages("fs")
library(fs)
files <- fs::dir_ls("GlobalRhizostomeDistribution/EnvironmentalData/Temperature/", glob="*.csv")
df <- read_csv(files, id="path")
head(df)
#install.packages("dplyr")
library(dplyr)

df <- dplyr::select(df, -c(path))
head(df)
df2 <- df %>% 
  rename(temp0 = 4, temp5 = 5, temp10 = 6, temp20 = 8)
head(df2)
df2<-dplyr::select(df2, c(Month, Latitude, Longitude, temp0, temp5,temp10, temp20))
head(df2)
ggplot(df2, aes(Longitude, Latitude, fill = temp0))+
  geom_tile()


obs_data3 <- tidyr::separate(obs_data, 'observed_on',
                       into = c('Month', 'Day', 'Year'),
                       sep= '/')
head(obs_data3)
df2 <- df2 %>% 
  rename("longitude" = "Longitude",
         "latitude" = "Latitude")
temp <- df2
head(temp)
#After all that, we have ea workable temperature database... 
#partitioned by monthly means at lat and long coordinates

###
###Let's incorporate more environmental data.
###

#Salinity
(files <- fs::dir_ls("GlobalRhizostomeDistribution/EnvironmentalData/Salinity/", 
                     glob="*.csv"))
salinity <- read_csv(files, id="path")
head(salinity)

salinity <- dplyr::select(salinity, -c(path))
sal2 <- salinity %>% 
  rename(sal0 = 4, sal5 = 5, sal10 = 6, sal20 = 8, sal30 = 10)
head(sal2)
sal2<-dplyr::select(sal2, c(Month, latitude, longitude, sal0, sal5, sal10, sal20))
head(sal2)
ggplot(sal2, aes(longitude, latitude, fill = sal5))+
  geom_tile()

sal2 <- sal2 %>% 
  group_by(Month, latitude, longitude) %>% 
  summarize(sal0=mean(sal0, na.rm = T), 
            sal5=mean(sal5, na.rm = T), 
            sal10=mean(sal10, na.rm = T), 
            sal20=mean(sal20, na.rm = T))
head(sal2)

#DO
(files <- fs::dir_ls("GlobalRhizostomeDistribution/EnvironmentalData/DO/", 
                     glob="*.csv"))
DO <- read_csv(files, id="path")
head(DO)

DO <- dplyr::select(DO, -c(path))
DO2 <- DO %>% 
  rename(DO0 = 4, DO5 = 5, DO10 = 6, DO20 = 8, DO30 = 10, DO50 = 14)
head(DO2)
DO2<-dplyr::select(DO2, c(Month, latitude, longitude, DO0, DO5, DO10, DO20))
head(DO2)
ggplot(DO2, aes(longitude, latitude, fill = DO5))+
  geom_tile()

DO2 <- DO2 %>% 
  group_by(Month, latitude, longitude) %>% 
  summarize(DO0=mean(DO0, na.rm = T), 
            DO5=mean(DO5, na.rm = T), 
            DO10=mean(DO10, na.rm = T), 
            DO20=mean(DO20, na.rm = T))
head(DO2)

##Nitrate
(files <- fs::dir_ls("GlobalRhizostomeDistribution/EnvironmentalData/Nitrate/", 
                     glob="*.csv"))
Nitrate <- read_csv(files, id="path")
head(Nitrate)

Nitrate <- dplyr::select(Nitrate, -c(path))
NO2 <- Nitrate %>% 
  rename(NO0 = 4, NO5 = 5, NO10 = 6, NO20 = 8, NO30 = 10, NO50 = 14)
head(NO2)
NO2<-dplyr::select(NO2, c(Month, latitude, longitude, NO0, NO5, NO10, NO20))
head(NO2)
ggplot(NO2, aes(longitude, latitude, fill = NO5))+
  geom_tile()

NO2 <- NO2 %>% 
  group_by(Month, latitude, longitude) %>% 
  summarize(NO0=mean(NO0, na.rm = T), 
            NO5=mean(NO5, na.rm = T), 
            NO10=mean(NO10, na.rm = T), 
            NO20=mean(NO20, na.rm = T))
head(NO2)

#OSat
(files <- fs::dir_ls("GlobalRhizostomeDistribution/EnvironmentalData/OSat/", 
                     glob="*.csv"))
OSat <- read_csv(files, id="path")
head(OSat)

OSat <- dplyr::select(OSat, -c(path))
OSat2 <- OSat %>% 
  rename(OSat0 = 4, OSat5 = 5, OSat10 = 6, OSat20 = 8)
head(OSat2)
OSat2<-dplyr::select(OSat2, c(Month, latitude, longitude, OSat0, OSat5, OSat10, OSat20))
head(OSat2)
ggplot(OSat2, aes(longitude, latitude, fill = OSat5))+
  geom_tile()

OSat2 <- OSat2 %>% 
  group_by(Month, latitude, longitude) %>% 
  summarize(OSat0=mean(OSat0, na.rm = T), 
            OSat5=mean(OSat5, na.rm = T), 
            OSat10=mean(OSat10, na.rm = T), 
            OSat20=mean(OSat20, na.rm = T))
head(OSat2)

#Phosphate
(files <- fs::dir_ls("GlobalRhizostomeDistribution/EnvironmentalData/Phosphate/", 
                     glob="*.csv"))
PO <- read_csv(files, id="path")
head(PO)

PO <- dplyr::select(PO, -c(path))
PO2 <- PO %>% 
  rename(PO0 = 4, PO5 = 5, PO10 = 6, PO20 = 8)
head(PO2)
PO2<-dplyr::select(PO2, c(Month, latitude, longitude, PO0, PO5, PO10, PO20))
head(PO2)
ggplot(PO2, aes(longitude, latitude, fill = PO5))+
  geom_tile()

PO2 <- PO2 %>% 
  group_by(Month, latitude, longitude) %>% 
  summarize(PO0=mean(PO0, na.rm = T), 
            PO5=mean(PO5, na.rm = T), 
            PO10=mean(PO10, na.rm = T), 
            PO20=mean(PO20, na.rm = T))
head(PO2)

#Silicate
(files <- fs::dir_ls("GlobalRhizostomeDistribution/EnvironmentalData/Silicate/", 
                     glob="*.csv"))
Silicate <- read_csv(files, id="path")
head(Silicate)

Silicate <- dplyr::select(Silicate, -c(path))
Sil2 <- Silicate %>% 
  rename(Sil0 = 4, Sil5 = 5, Sil10 = 6, Sil20 = 8)
head(Sil2)
Sil2<-dplyr::select(Sil2, c(Month, latitude, longitude, Sil0, Sil5, Sil10, Sil20))
head(Sil2)
ggplot(Sil2, aes(longitude, latitude, fill = Sil5))+
  geom_tile()

Sil2 <- Sil2 %>% 
  group_by(Month, latitude, longitude) %>% 
  summarize(Sil0=mean(Sil0, na.rm = T), 
            Sil5=mean(Sil5, na.rm = T), 
            Sil10=mean(Sil10, na.rm = T), 
            Sil20=mean(Sil20, na.rm = T))
head(Sil2)

sal2 <- as.data.frame(sal2)
temp <- as.data.frame(temp)

head(temp)
head(sal2)

#Merge Environmental data
Wombo<-merge(temp, sal2, by=c("Month","longitude","latitude"))
Wombo<-merge(Wombo, DO2, by=c("Month","longitude","latitude"))
Wombo<-merge(Wombo, NO2, by=c("Month","longitude","latitude"))
Wombo<-merge(Wombo, OSat2, by=c("Month","longitude","latitude"))
Wombo<-merge(Wombo, PO2, by=c("Month","longitude","latitude"))
Wombo<-merge(Wombo, Sil2, by=c("Month","longitude","latitude"))

ggplot(Wombo, aes(longitude, latitude, fill = temp0))+
  geom_tile()
#We lost some resolution, due to sparse sampling in some variables, 
#but this will be corrected when join the dataset to jellyfish observations

#install.packages("powerjoin")
library(powerjoin)

head(filtered_obs_data$latitude)
filtered_obs_data <- tidyr::separate(filtered_obs_data, 'observed_on',
                             into = c('Month', 'Day', 'Year'),
                             sep= '/')
head(filtered_obs_data)
class(Wombo)
class(filtered_obs_data)
filtered_obs_data <- as.data.frame(filtered_obs_data)
selected_obs_data <- dplyr::select(filtered_obs_data, 
                            c(id, Month, Day, Year, 
                                                 quality_grade, latitude, longitude, 
                                                 species_guess, scientific_name, 
                              common_name, taxon_order_name, taxon_suborder_name, 
                              taxon_family_name, taxon_genus_name, taxon_species_name))

head(selected_obs_data)
nrow(selected_obs_data)
#7808
library(stats)
#assign environment 2 degrees each way
tol <- 2
NRC2 <- power_left_join(selected_obs_data, Wombo,
                        by = c(~ sqrt((.x$latitude - .y$latitude)^2 + (.x$longitude - .y$longitude)^2) < tol))

#Pair up months
NRC2 = filter(NRC2, Month.x == Month.y)
NRCclean2 <- NRC2 %>%
  group_by(id, Month.x, latitude.x, longitude.x, taxon_suborder_name, taxon_genus_name) %>%
  summarize(temp0=mean(temp0, na.rm = T), 
            temp5=mean(temp5, na.rm = T), 
            temp10=mean(temp10, na.rm = T), 
            temp20=mean(temp20, na.rm = T),
            sal0=mean(sal0, na.rm = T), 
            sal5=mean(sal5, na.rm = T), 
            sal10=mean(sal10, na.rm = T), 
            sal20=mean(sal20, na.rm = T),
            DO0=mean(DO0, na.rm = T), 
            DO5=mean(DO5, na.rm = T), 
            DO10=mean(DO10, na.rm = T), 
            DO20=mean(DO20, na.rm = T),
            NO0=mean(NO0, na.rm = T), 
            NO5=mean(NO5, na.rm = T), 
            NO10=mean(NO10, na.rm = T), 
            NO20=mean(NO20, na.rm = T),
            OSat0=mean(OSat0, na.rm = T), 
            OSat5=mean(OSat5, na.rm = T), 
            OSat10=mean(OSat10, na.rm = T), 
            OSat20=mean(OSat20, na.rm = T),
            PO0=mean(PO0, na.rm = T), 
            PO5=mean(PO5, na.rm = T), 
            PO10=mean(PO10, na.rm = T), 
            PO20=mean(PO20, na.rm = T),
            Sil0=mean(Sil0, na.rm = T), 
            Sil5=mean(Sil5, na.rm = T), 
            Sil10=mean(Sil10, na.rm = T), 
            Sil20=mean(Sil20, na.rm = T)
  )
NRCclean2 <- na.omit(NRCclean2)
nrow(NRCclean2)
#4662
#Lost 3146 observations at this coordinate resolution (7808 - 4662)

#I don't think we should go more resolved but let's check
#Anything over 2 degrees would also be too course and inaccurate

#Assign environment by 1 degree each way
tol <- 1
NRC1 <- power_left_join(selected_obs_data, Wombo,
                        by = c(~ sqrt((.x$latitude - .y$latitude)^2 + (.x$longitude - .y$longitude)^2) < tol))

NRC1 = filter(NRC1, Month.x == Month.y)
NRCclean1 <- NRC1 %>%
  group_by(id, Month.x, latitude.x, longitude.x, taxon_suborder_name, taxon_genus_name) %>%
  summarize(temp0=mean(temp0, na.rm = T), 
            temp5=mean(temp5, na.rm = T), 
            temp10=mean(temp10, na.rm = T), 
            temp20=mean(temp20, na.rm = T),
            sal0=mean(sal0, na.rm = T), 
            sal5=mean(sal5, na.rm = T), 
            sal10=mean(sal10, na.rm = T), 
            sal20=mean(sal20, na.rm = T),
            DO0=mean(DO0, na.rm = T), 
            DO5=mean(DO5, na.rm = T), 
            DO10=mean(DO10, na.rm = T), 
            DO20=mean(DO20, na.rm = T),
            NO0=mean(NO0, na.rm = T), 
            NO5=mean(NO5, na.rm = T), 
            NO10=mean(NO10, na.rm = T), 
            NO20=mean(NO20, na.rm = T),
            OSat0=mean(OSat0, na.rm = T), 
            OSat5=mean(OSat5, na.rm = T), 
            OSat10=mean(OSat10, na.rm = T), 
            OSat20=mean(OSat20, na.rm = T),
            PO0=mean(PO0, na.rm = T), 
            PO5=mean(PO5, na.rm = T), 
            PO10=mean(PO10, na.rm = T), 
            PO20=mean(PO20, na.rm = T),
            Sil0=mean(Sil0, na.rm = T), 
            Sil5=mean(Sil5, na.rm = T), 
            Sil10=mean(Sil10, na.rm = T), 
            Sil20=mean(Sil20, na.rm = T)
  )
NRCclean1 <- na.omit(NRCclean1)
nrow(NRCclean1)
#Left with 1962 observations;
#We lost too many with this.
#I will be continuing with 2 degree dataframe, 
#but let's visualize to make sure we don't bias too heavily
head(NRCclean1)

NRCunrounded1 <- NRCclean1
NRCunrounded2 <- NRCclean2
#See if filtering changed distributions
lat_initial <- ggplot(selected_obs_data, aes(x = latitude, y = fct_infreq(taxon_genus_name), fill = after_stat(x))
) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Lat", option = "D") +
  facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") + 
  theme_minimal()
lat_initial

lat_filtered1 <- ggplot(NRCunrounded1, aes(x = latitude.x, y = fct_infreq(taxon_genus_name), fill = after_stat(x))
) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Lat", option = "D") +
  facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") + 
  theme_minimal()
lat_filtered2 <- ggplot(NRCunrounded2, aes(x = latitude.x, y = fct_infreq(taxon_genus_name), fill = after_stat(x)))+
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Lat", option = "D") +
  facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") + 
  theme_minimal()
grid.arrange(lat_initial, lat_filtered1, lat_filtered2)

lon_initial <- ggplot(selected_obs_data, aes(x = longitude, y = fct_infreq(taxon_genus_name), fill = after_stat(x))
) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Lat", option = "D") +
  facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") + 
  theme_minimal()
lon_filtered1 <- ggplot(NRCunrounded1, aes(x = longitude.x, y = fct_infreq(taxon_genus_name), fill = after_stat(x))
) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Lat", option = "D") +
  facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") + 
  theme_minimal()
lon_filtered2 <- ggplot(NRCunrounded2, aes(x = longitude.x, y = fct_infreq(taxon_genus_name), fill = after_stat(x))
) +
  geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Lat", option = "D") +
  facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside") + 
  theme_minimal()
grid.arrange(lon_initial, lon_filtered1, lon_filtered2)

#subset data for some statistical tests
latog <- dplyr::select(selected_obs_data, c(id, Month, taxon_suborder_name, taxon_genus_name, latitude, longitude))
lat1 <- dplyr::select(NRCunrounded1, c(id, taxon_genus_name, latitude.x, longitude.x))
lat2 <- dplyr::select(NRCunrounded2, c(id, taxon_genus_name, latitude.x, longitude.x))

latog$round <- 0
lat1$round <- 1
lat2$round <- 2

latog$Month.x <- latog$Month
latog$latitude.x <- latog$latitude
latog$longitude.x <- latog$longitude
CoordinateComparison <- rbind(lat1, lat2)
latog <- dplyr::select(latog, -c(Month, longitude, latitude))
CoordinateComparison <- rbind(latog, CoordinateComparison)

ks.test(lat1$latitude.x, latog$latitude.x)
ks.test(lat2$latitude.x, latog$latitude.x)

ks.test(lat1$longitude.x, latog$longitude.x)
ks.test(lat2$longitude.x, latog$longitude.x)

#While significantly different, the effect size of the 2 degree filter is comparatively small to the 1 degree filter
#Combine this with the fact that 2 degrees is the most coarse level we should go... 
#I'm comfortable with 2 degrees.
#Citizen science and publically available data will need to improve for better resolution

###
### Genera correlations
###

NRCclean1$latitude.x <- round(NRCclean1$latitude.x, digits = 1)
NRCclean1$longitude.x <- round(NRCclean1$longitude.x, digits = 1)

abund1 <- NRCclean1 %>%
  group_by(latitude.x, longitude.x, taxon_suborder_name, taxon_genus_name, temp0, sal0, DO0, NO0, PO0, OSat0, Sil0) %>%
  summarize(observations = n())
abund1
tempdataabund1<-dplyr::select(abund1, c(temp0, sal0, DO0, NO0, PO0, OSat0, Sil0))
class(tempdataabund1)
tempdataabund1 <- as.data.frame(tempdataabund1)
tempdataabund1<-dplyr::select(tempdataabund1, c(temp0, sal0, DO0, NO0, PO0, OSat0, Sil0))
head(tempdataabund1)

genus.wide1 <- pivot_wider(abund1, names_from = taxon_genus_name, values_from = observations) 
head(genus.wide1)
class(genus.wide1)
genus.wide1 <- as.data.frame(genus.wide1)
genus.wide1<-dplyr::select(genus.wide1, -c(latitude.x, longitude.x, taxon_suborder_name, temp0, sal0, DO0, NO0, PO0, OSat0, Sil0))
head(genus.wide1)
genus.wide1[is.na(genus.wide1)] <- 0
head(genus.wide1)

#install.packages("mvabund")
#install.packages("ecoCopula")
library(mvabund)
library(ecoCopula)
spiderAbund1 = mvabund(genus.wide1)
head(spiderAbund1)

spider_glmInt1 = manyglm(spiderAbund1~1,data=tempdataabund1)

ord_spiderInt1 = cord(spider_glmInt1)
plot(ord_spiderInt1, biplot = TRUE) #for a biplot
#Stomolophus and Rhizostoma most different at 1 degree rounding

###Now with 2 degree rounding
###
### Genera correlations
###

NRCclean2$latitude.x <- round(NRCclean2$latitude.x, digits = 1)
NRCclean2$longitude.x <- round(NRCclean2$longitude.x, digits = 1)

abund2 <- NRCclean2 %>%
  group_by(latitude.x, longitude.x, taxon_suborder_name, taxon_genus_name, temp0, sal0, DO0, NO0, PO0, OSat0, Sil0) %>%
  summarize(observations = n())
abund2 <- na.omit(abund2)

genus.wide2 <- pivot_wider(abund2, names_from = taxon_genus_name, values_from = observations) 
head(genus.wide2)
class(genus.wide2)
genus.wide2 <- as.data.frame(genus.wide2)

tempdataabund2<-dplyr::select(genus.wide2, c(temp0, sal0, DO0, NO0, PO0, OSat0, Sil0))
nrow(tempdataabund2)

genus.wide2<-dplyr::select(genus.wide2, -c(latitude.x, longitude.x, taxon_suborder_name, temp0, sal0, DO0, NO0, PO0, OSat0, Sil0))
head(genus.wide2)
genus.wide2[is.na(genus.wide2)] <- 0
head(genus.wide2)
spiderAbund2 = mvabund(genus.wide2)
#Niche model for 2 degree rounding
spider_glmInt2 = manyglm(spiderAbund2~1,data=tempdataabund2)
ord_spiderInt2 = cord(spider_glmInt2)
plot(ord_spiderInt2, biplot = TRUE) #for a biplot

par(mfrow=c(1,2))
plot(ord_spiderInt1, biplot = TRUE) #for a biplot
plot(ord_spiderInt2, biplot = TRUE) #for a biplot
#SAME NICHE RESULTS INDEPENDENT OF JOIN RESOLUTION

###
#Statistical different distributions after filtering
#However, the 1 and 2 lat filtered comparisons produced identical glms
#Random forest models were also minimally different even at 5 degrees coordinate rounding (Not in this code)
#Therefore, I would move forward, to conclusions. I like the 2 degree one for a balance between sampling, and resolution
#This however is a bit arbitrary and I must move forward with caution
###

###
#Random forest model
###
#This will tell us which environmental paramaters are most influential

#install.packages("dplyr")
library(dplyr)
#install.packages("randomForest")
library(randomForest)
#install.packages("caret")
library(caret)
citation("caret")
packageVersion("caret")
#install.packages("e1071")
library(e1071)
citation("e1071")
packageVersion("e1071")
#install.packages("caTools")    # For Logistic regression 
library(caTools)
citation("caTools")
packageVersion("caTools")
par(mfrow=c(1,1))

#Convert data to observation numbers, this will later be converted to binary metrics
randomforest<- NRCunrounded2 %>%
  group_by(Month.x, latitude.x, longitude.x, taxon_suborder_name, taxon_genus_name, 
           temp0, temp5, temp10, temp20,  
           sal0, sal5, sal10, sal20,   
           DO0, DO5, DO10, DO20, 
           NO0, NO5, NO10, NO20, 
           OSat0, OSat5, OSat10, OSat20, 
           PO0, PO5, PO10, PO20, 
           Sil0, Sil5, Sil10, Sil20) %>%
  summarize(observations = n())


randomforest <- pivot_wider(randomforest, names_from = taxon_genus_name, values_from = observations) 
randomforest <- as.data.frame(randomforest)
randomforest <- dplyr::select(randomforest, -c(Month.x, latitude.x, longitude.x, taxon_suborder_name))
randomforest[is.na(randomforest)] <- 0
head(randomforest)

#subset into genera for seperate models

#Just gonna start with pseudorhiza
PseudoOccur <- dplyr::select(randomforest, -c(Catostylus, Lychnorhiza, Rhizostoma, Eupilema, Cassiopea, Stomolophus,Rhopilema, Cotylorhiza, Lobonema))
PseudoOccur$Pseudorhiza[PseudoOccur$Pseudorhiza > 0] <- 'Yes'
PseudoOccur$Pseudorhiza[PseudoOccur$Pseudorhiza == 0] <- 'No'
head(PseudoOccur)

split <- sample.split(PseudoOccur, SplitRatio = 0.8) 
split
data_train <- subset(PseudoOccur, split == "TRUE")
dim(data_train)
head(data_train)
data_test <- subset(PseudoOccur, split == "FALSE")
dim(data_test)
head(data_test)
head(data_train)

PseudoOccur$Pseudorhiza <- as.factor(PseudoOccur$Pseudorhiza)    
data_train$Pseudorhiza <- as.factor(data_train$Pseudorhiza)
data_test$Pseudorhiza <- as.factor(data_test$Pseudorhiza)

model <- randomForest(Pseudorhiza~.,data= data_train, mtry = 5)
importance(model)
varImpPlot(model)
pseudomodel <- model
pred_test <- predict(model, newdata = data_test, type= "class")
pred_test

pseudoconfusion <- confusionMatrix(table(pred_test,data_test$Pseudorhiza))
pseudoconfusion
#Catostylus model
CatoOccur <- dplyr::select(randomforest, -c(Pseudorhiza, Lychnorhiza, Rhizostoma, Eupilema, Cassiopea, Stomolophus,Rhopilema, Cotylorhiza, Lobonema))
head(CatoOccur)
CatoOccur$Catostylus[CatoOccur$Catostylus > 0] <- 'Yes'
CatoOccur$Catostylus[CatoOccur$Catostylus == 0] <- 'No'
#View(CatoOccur)

split <- sample.split(CatoOccur, SplitRatio = 0.8) 
split
data_train <- subset(CatoOccur, split == "TRUE")
dim(data_train)
head(data_train)
data_test <- subset(CatoOccur, split == "FALSE")
dim(data_test)
head(data_test)

CatoOccur$Catostylus <- as.factor(CatoOccur$Catostylus)    
data_train$Catostylus <- as.factor(data_train$Catostylus)
data_test$Catostylus <- as.factor(data_test$Catostylus)

head(data_train)
model <- randomForest(Catostylus~., mtry = 5, data= data_train)
importance(model)
varImpPlot(model)
catomodel <- model
pred_test <- predict(model, newdata = data_test, type= "class")

catoconfusion <- confusionMatrix(table(pred_test,data_test$Catostylus))

#Lychnorhiza model
LychnorhizaOccur <- dplyr::select(randomforest, -c(Pseudorhiza, Catostylus, Rhizostoma, Eupilema, Cassiopea, Stomolophus,Rhopilema, Cotylorhiza, Lobonema))
head(LychnorhizaOccur)
LychnorhizaOccur$Lychnorhiza[LychnorhizaOccur$Lychnorhiza > 0] <- 'Yes'
LychnorhizaOccur$Lychnorhiza[LychnorhizaOccur$Lychnorhiza == 0] <- 'No'
head(LychnorhizaOccur)

split <- sample.split(LychnorhizaOccur, SplitRatio = 0.8) 
split
data_train <- subset(LychnorhizaOccur, split == "TRUE")
dim(data_train)
head(data_train)
data_test <- subset(LychnorhizaOccur, split == "FALSE")
dim(data_test)
head(data_test)

LychnorhizaOccur$Lychnorhiza <- as.factor(LychnorhizaOccur$Lychnorhiza)    
data_train$Lychnorhiza <- as.factor(data_train$Lychnorhiza)
data_test$Lychnorhiza <- as.factor(data_test$Lychnorhiza)

head(data_train)
model <- randomForest(Lychnorhiza~., mtry = 5, data= data_train)
importance(model)
varImpPlot(model)
Lychnorhizamodel<-model
pred_test <- predict(model, newdata = data_test, type= "class")

Lychnorhizaconfusion <- confusionMatrix(table(pred_test,data_test$Lychnorhiza))
Lychnorhizaconfusion

LychnorhizaOccur <- dplyr::select(randomforest, -c(Pseudorhiza, Catostylus, Rhizostoma, Eupilema, Cassiopea, Stomolophus,Rhopilema, Cotylorhiza, Lobonema))
head(LychnorhizaOccur)
LychnorhizaOccur$Lychnorhiza[LychnorhizaOccur$Lychnorhiza > 0] <- 'Yes'
LychnorhizaOccur$Lychnorhiza[LychnorhizaOccur$Lychnorhiza == 0] <- 'No'
head(LychnorhizaOccur)

split <- sample.split(LychnorhizaOccur, SplitRatio = 0.8) 
split
data_train <- subset(LychnorhizaOccur, split == "TRUE")
dim(data_train)
head(data_train)
data_test <- subset(LychnorhizaOccur, split == "FALSE")
dim(data_test)
head(data_test)

LychnorhizaOccur$Lychnorhiza <- as.factor(LychnorhizaOccur$Lychnorhiza)    
data_train$Lychnorhiza <- as.factor(data_train$Lychnorhiza)
data_test$Lychnorhiza <- as.factor(data_test$Lychnorhiza)

head(data_train)

model <- randomForest(Lychnorhiza~., mtry = 6, data= data_train)
importance(model)
varImpPlot(model)
Lychnorhizamodel<-model
pred_test <- predict(model, newdata = data_test, type= "class")

Lychnorhizaconfusion <- confusionMatrix(table(pred_test,data_test$Lychnorhiza))
Lychnorhizaconfusion

#Rhizostoma model
RhizOccur <- dplyr::select(randomforest, -c(Pseudorhiza, Catostylus, Lychnorhiza, Eupilema, Cassiopea, Stomolophus,Rhopilema, Cotylorhiza, Lobonema))
head(RhizOccur)
RhizOccur$Rhizostoma[RhizOccur$Rhizostoma > 0] <- 'Yes'
RhizOccur$Rhizostoma[RhizOccur$Rhizostoma == 0] <- 'No'
head(RhizOccur)

split <- sample.split(RhizOccur, SplitRatio = 0.8) 
split
data_train <- subset(RhizOccur, split == "TRUE")
dim(data_train)
head(data_train)
data_test <- subset(RhizOccur, split == "FALSE")
dim(data_test)
head(data_test)

RhizOccur$Rhizostoma <- as.factor(RhizOccur$Rhizostoma)    
data_train$Rhizostoma <- as.factor(data_train$Rhizostoma)
data_test$Rhizostoma <- as.factor(data_test$Rhizostoma)

head(data_train)
model <- randomForest(Rhizostoma~., mtry = 7, data= data_train)
importance(model)
varImpPlot(model)
rhizomodel <- model
pred_test <- predict(model, newdata = data_test, type= "class")

rhizoconfusion <- confusionMatrix(table(pred_test,data_test$Rhizostoma))
rhizoconfusion

#Eupilema
EupilemaOccur <- dplyr::select(randomforest, -c(Pseudorhiza, Catostylus, Lychnorhiza, Rhizostoma, Cassiopea, Stomolophus, Rhopilema, Cotylorhiza, Lobonema))
head(EupilemaOccur)
EupilemaOccur$Eupilema[EupilemaOccur$Eupilema > 0] <- 'Yes'
EupilemaOccur$Eupilema[EupilemaOccur$Eupilema == 0] <- 'No'
head(EupilemaOccur)

split <- sample.split(EupilemaOccur, SplitRatio = 0.8) 
split
data_train <- subset(EupilemaOccur, split == "TRUE")
dim(data_train)
head(data_train)
data_test <- subset(EupilemaOccur, split == "FALSE")
dim(data_test)
head(data_test)

EupilemaOccur$Eupilema <- as.factor(EupilemaOccur$Eupilema)    
data_train$Eupilema <- as.factor(data_train$Eupilema)
data_test$Eupilema <- as.factor(data_test$Eupilema)

head(data_train)
model <- randomForest(Eupilema~., mtry = 5, data= data_train)
importance(model)
varImpPlot(model)
Eupilemamodel<-model
pred_test <- predict(model, newdata = data_test, type= "class")

Eupilemaconfusion <- confusionMatrix(table(pred_test,data_test$Eupilema))
Eupilemaconfusion

#Cassiopea model
CassOccur <- dplyr::select(randomforest, -c(Pseudorhiza, Catostylus, Lychnorhiza, Rhizostoma, Eupilema, Stomolophus,Rhopilema, Cotylorhiza, Lobonema))
head(CassOccur)
CassOccur$Cassiopea[CassOccur$Cassiopea > 0] <- 'Yes'
CassOccur$Cassiopea[CassOccur$Cassiopea == 0] <- 'No'
head(CassOccur)

split <- sample.split(CassOccur, SplitRatio = 0.8) 
split
data_train <- subset(CassOccur, split == "TRUE")
dim(data_train)
head(data_train)
data_test <- subset(CassOccur, split == "FALSE")
dim(data_test)
head(data_test)

CassOccur$Cassiopea <- as.factor(CassOccur$Cassiopea)    
data_train$Cassiopea <- as.factor(data_train$Cassiopea)
data_test$Cassiopea <- as.factor(data_test$Cassiopea)

head(data_train)
model <- randomForest(Cassiopea~., mtry = 6, data= data_train)
importance(model)
varImpPlot(model)
cassiomodel <- model
pred_test <- predict(model, newdata = data_test, type= "class")

cassioconfusion <- confusionMatrix(table(pred_test,data_test$Cassiopea))
cassioconfusion

#Stomolophus model
StomOccur <- dplyr::select(randomforest, -c(Pseudorhiza, Catostylus, Lychnorhiza, Rhizostoma, Eupilema, Cassiopea, Rhopilema, Cotylorhiza, Lobonema))
head(StomOccur)
StomOccur$Stomolophus[StomOccur$Stomolophus > 0] <- 'Yes'
StomOccur$Stomolophus[StomOccur$Stomolophus == 0] <- 'No'
head(StomOccur)

split <- sample.split(StomOccur, SplitRatio = 0.8) 
split
data_train <- subset(StomOccur, split == "TRUE")
dim(data_train)
head(data_train)
data_test <- subset(StomOccur, split == "FALSE")
dim(data_test)
head(data_test)

StomOccur$Stomolophus <- as.factor(StomOccur$Stomolophus)    
data_train$Stomolophus <- as.factor(data_train$Stomolophus)
data_test$Stomolophus <- as.factor(data_test$Stomolophus)

head(data_train)

model <- randomForest(Stomolophus~., mtry = 6, data= data_train)
importance(model)
varImpPlot(model)
stommodel <- model
pred_test <- predict(model, newdata = data_test, type= "class")

stomconfusion <- confusionMatrix(table(pred_test,data_test$Stomolophus))
stomconfusion

#Rhopilema model
RhopOccur <- dplyr::select(randomforest, -c(Pseudorhiza, Catostylus, Lychnorhiza, Rhizostoma, Eupilema, Cassiopea, Stomolophus, Cotylorhiza, Lobonema))
head(RhopOccur)
RhopOccur$Rhopilema[RhopOccur$Rhopilema > 0] <- 'Yes'
RhopOccur$Rhopilema[RhopOccur$Rhopilema == 0] <- 'No'
head(RhopOccur)

split <- sample.split(RhopOccur, SplitRatio = 0.8) 
split
data_train <- subset(RhopOccur, split == "TRUE")
dim(data_train)
head(data_train)
data_test <- subset(RhopOccur, split == "FALSE")
dim(data_test)
head(data_test)

RhopOccur$Rhopilema <- as.factor(RhopOccur$Rhopilema)    
data_train$Rhopilema <- as.factor(data_train$Rhopilema)
data_test$Rhopilema <- as.factor(data_test$Rhopilema)

model <- randomForest(Rhopilema~., mtry = 5, data= data_train)
importance(model)
varImpPlot(model)
rhopmodel <- model
pred_test <- predict(model, newdata = data_test, type= "class")

rhopconfusion <- confusionMatrix(table(pred_test,data_test$Rhopilema))
rhopconfusion

#Cotylorhiza model
CotylOccur <- dplyr::select(randomforest, -c(Pseudorhiza, Catostylus, Lychnorhiza, Rhizostoma, Eupilema, Cassiopea, Stomolophus,Rhopilema, Lobonema))
head(CotylOccur)
CotylOccur$Cotylorhiza[CotylOccur$Cotylorhiza > 0] <- 'Yes'
CotylOccur$Cotylorhiza[CotylOccur$Cotylorhiza == 0] <- 'No'
head(CotylOccur)

split <- sample.split(CotylOccur, SplitRatio = 0.8) 
split
data_train <- subset(CotylOccur, split == "TRUE")
dim(data_train)
head(data_train)
data_test <- subset(CotylOccur, split == "FALSE")
dim(data_test)
head(data_test)

CotylOccur$Cotylorhiza <- as.factor(CotylOccur$Cotylorhiza)    
data_train$Cotylorhiza <- as.factor(data_train$Cotylorhiza)
data_test$Cotylorhiza <- as.factor(data_test$Cotylorhiza)

head(data_train)

model <- randomForest(Cotylorhiza~., mtry = 5, data= data_train)
importance(model)
varImpPlot(model)
cotylmodel <- model
pred_test <- predict(model, newdata = data_test, type= "class")

cotylconfusion <- confusionMatrix(table(pred_test,data_test$Cotylorhiza))
cotylconfusion

#Lobonema
LobonemaOccur <- dplyr::select(randomforest, -c(Pseudorhiza, Catostylus, Rhizostoma, Cassiopea, Stomolophus,Rhopilema, Cotylorhiza, Lychnorhiza, Eupilema))
head(LobonemaOccur)
LobonemaOccur$Lobonema[LobonemaOccur$Lobonema > 0] <- 'Yes'
LobonemaOccur$Lobonema[LobonemaOccur$Lobonema == 0] <- 'No'
head(LobonemaOccur)

split <- sample.split(LobonemaOccur, SplitRatio = 0.8) 
split
data_train <- subset(LobonemaOccur, split == "TRUE")
dim(data_train)
head(data_train)
data_test <- subset(LobonemaOccur, split == "FALSE")
dim(data_test)
head(data_test)

LobonemaOccur$Lobonema <- as.factor(LobonemaOccur$Lobonema)    
data_train$Lobonema <- as.factor(data_train$Lobonema)
data_test$Lobonema <- as.factor(data_test$Lobonema)

model <- randomForest(Lobonema~., mtry = 5, data= data_train)
importance(model)
varImpPlot(model)
Lobonemamodel<-model
pred_test <- predict(model, newdata = data_test, type= "class")

Lobonemaconfusion <- confusionMatrix(table(pred_test,data_test$Lobonema))
Lobonemaconfusion

#Final models
par(mfrow=c(1,1))
rhizomodel
rhizoconfusion

rhizoplot<-varImpPlot(rhizomodel)
rhizoplot <- as.data.frame(rhizoplot)
rhizoplot$varnames <- rownames(rhizoplot) # row names to column
rownames(rhizoplot) <- NULL 
head(rhizoplot)
rhizoplot <- cbind(Genus = c("Rhizostoma"), rhizoplot)

rhopmodel
rhopconfusion

rhopplot<-varImpPlot(rhopmodel)
rhopplot <- as.data.frame(rhopplot)
rhopplot$varnames <- rownames(rhopplot) # row names to column
rownames(rhopplot) <- NULL 
head(rhopplot)
rhopplot <- cbind(Genus = c("Rhopilema"), rhopplot)

cassiomodel
cassioconfusion

cassioplot<-varImpPlot(cassiomodel)
cassioplot <- as.data.frame(cassioplot)
cassioplot$varnames <- rownames(cassioplot) # row names to column
rownames(cassioplot) <- NULL 
head(cassioplot)
cassioplot <- cbind(Genus = c("Cassiopea"), cassioplot)

pseudomodel
pseudoconfusion

pseudoplot<-varImpPlot(pseudomodel)
pseudoplot <- as.data.frame(pseudoplot)
pseudoplot$varnames <- rownames(pseudoplot) # row names to column
rownames(pseudoplot) <- NULL 
head(pseudoplot)
pseudoplot <- cbind(Genus = c("Pseudorhiza"), pseudoplot)

stommodel
stomconfusion

stomplot<-varImpPlot(stommodel)
stomplot <- as.data.frame(stomplot)
stomplot$varnames <- rownames(stomplot) # row names to column
rownames(stomplot) <- NULL 
head(stomplot)
stomplot <- cbind(Genus = c("Stomolophus"), stomplot)

catomodel
catoconfusion

catoplot<-varImpPlot(catomodel)
catoplot <- as.data.frame(catoplot)
catoplot$varnames <- rownames(catoplot) # row names to column
rownames(catoplot) <- NULL 
head(catoplot)
catoplot <- cbind(Genus = c("Catostylus"), catoplot)

cotylmodel
cotylconfusion

cotylplot<-varImpPlot(cotylmodel)
cotylplot <- as.data.frame(cotylplot)
cotylplot$varnames <- rownames(cotylplot) # row names to column
rownames(cotylplot) <- NULL 
head(cotylplot)
cotylplot <- cbind(Genus = c("Cotylorhiza"), cotylplot)


Eupilemamodel
Eupilemaconfusion

Eupilemaplot<-varImpPlot(Eupilemamodel)
Eupilemaplot <- as.data.frame(Eupilemaplot)
Eupilemaplot$varnames <- rownames(Eupilemaplot) # row names to column
rownames(Eupilemaplot) <- NULL 
head(Eupilemaplot)
Eupilemaplot <- cbind(Genus = c("Eupilema"), Eupilemaplot)

Lychnorhizamodel
Lychnorhizaconfusion

Lychnorhizaplot<-varImpPlot(Lychnorhizamodel)
Lychnorhizaplot <- as.data.frame(Lychnorhizaplot)
Lychnorhizaplot$varnames <- rownames(Lychnorhizaplot) # row names to column
rownames(Lychnorhizaplot) <- NULL 
head(Lychnorhizaplot)
Lychnorhizaplot <- cbind(Genus = c("Lychnorhiza"), Lychnorhizaplot)

Lobonemamodel
Lobonemaconfusion

Lobonemaplot<-varImpPlot(Lobonemamodel)
Lobonemaplot <- as.data.frame(Lobonemaplot)
Lobonemaplot$varnames <- rownames(Lobonemaplot) # row names to column
rownames(Lobonemaplot) <- NULL 
head(Lobonemaplot)
Lobonemaplot <- cbind(Genus = c("Lobonema"), Lobonemaplot)

head(Lobonemaplot)

#Let's combine all random forest models now and compare how each genus (the model) 
#is influenced by different variables
IMP <- rbind(pseudoplot, rhopplot)
IMP <- rbind(IMP, catoplot)
IMP <- rbind(IMP, rhizoplot)
IMP <- rbind(IMP, stomplot)
IMP <- rbind(IMP, cassioplot)
IMP <- rbind(IMP, cotylplot)
IMP <- rbind(IMP, Lobonemaplot)
IMP <- rbind(IMP, Lychnorhizaplot)
IMP <- rbind(IMP, Eupilemaplot)

IMP$varnames <- factor(IMP$varnames, levels = c("sal0", "sal5", "sal10", "sal20", "sal50",
                                                "temp0", "temp5", "temp10", "temp20", "temp50",
                                                "Sil0", "Sil5", "Sil10", "Sil20", "Sil50",
                                                "DO0", "DO5", "DO10", "DO20", "DO50",
                                                "OSat0", "OSat5", "OSat10", "OSat20", "OSat50",
                                                "NO0","NO5","NO10","NO20","NO50",
                                                "PO0","PO5","PO10","PO20","PO50"))
IMP$Genus <- factor(IMP$Genus, levels = c("Cassiopea", "Catostylus", "Cotylorhiza", 
                                          "Eupilema", "Lobonema", "Lychnorhiza", 
                                          "Pseudorhiza", "Rhizostoma", "Rhopilema", 
                                          "Stomolophus"))
head(IMP)

IMP %>%
  ggplot(aes(Genus, weight = MeanDecreaseGini, fill = varnames))+
    geom_bar(alpha = 0.8, position = "fill")+
    scale_fill_manual(
      values = c("temp0" = "#D55E00","temp5"= "#B34200", "temp10"="#922500","temp20"= "#730000",
                 "DO0"="#56B4E9", "DO5" = "#1F8FC2","DO10"= "#006B9C", "DO20" = "#004A77",
                 "Sil0" = "#F0E442",  "Sil5" = "#BBB300", "Sil10" = "#878400",  "Sil20"="#575800",
                 "sal0" = "#CC79A7","sal5" = "#AA5987", "sal10" = "#883B68", "sal20" = "#681C4B",
                 "OSat0" = "#A2E5FF", "OSat5" = "#81C7FF", "OSat10" = "#60AAEF", "OSat20" = "#3D8DD0",
                 "NO0" = "#009E73", "NO5" = "#008E7D", "NO10" = "#007C7F", "NO20" = "#006B79",
                 "PO0" = "#E69F00", "PO5" = "#BA7A00", "PO10" = "#905800", "PO20" = "#6A3700"))+
   scale_x_discrete(limits = rev)+
    scale_y_reverse()+
    coord_flip()+
    ylab("Relative Importance to Model")+
    theme_classic()

NRCunrounded2$taxon_genus_name <- factor(NRCunrounded2$taxon_genus_name, 
                                         levels = c("Cassiopea", "Catostylus", "Cotylorhiza", 
                                                     "Eupilema", "Lobonema", "Lychnorhiza", 
                                                     "Pseudorhiza", "Rhizostoma", "Rhopilema", 
                                                     "Stomolophus"))
###Now that we know what environmental variables help us predict jellyfish occurence,
#Let's look at the underlying data
tempplot <- ggplot(NRCunrounded2, aes(y = taxon_genus_name)) +
  #geom_density_ridges(aes(x = temp50), fill = "#570000", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = temp20), fill = "#730000", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = temp10), fill = "#922500", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = temp5), fill = "#B34200", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = temp0), fill = "#D55E00", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  #facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  scale_y_discrete(limits = rev)+
  theme_minimal()+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 1))
tempplot

#Genus temp comparison
kruskal.test(temp0 ~ taxon_genus_name, data = NRCunrounded2)
#Post-hoc comparison of genera
tempstats <- dunnTest(temp0 ~ taxon_genus_name, data = NRCunrounded2, method = "bh")
tempstats
tempstats = tempstats$res
cldList(comparison = tempstats$Comparison, p.value = tempstats$P.adj, threshold = 0.05)
#There is some clumping, but mostly genus specific

salplot <- ggplot(NRCunrounded2, aes(y = taxon_genus_name)) +
  #geom_density_ridges(aes(x = sal50), fill = "#48002F", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = sal20), fill = "#681C4B", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = sal10), fill = "#883B68", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = sal5), fill = "#AA5987", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = sal0), fill = "#CC79A7", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  #facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(25,41))+
  scale_y_discrete(limits = rev)+
  theme_minimal()+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 1))
salplot

#Genus sal comparison
kruskal.test(sal0 ~ taxon_genus_name, data = NRCunrounded2)
#Post-hoc comparison of genera
salstats <- dunnTest(sal0 ~ taxon_genus_name, data = NRCunrounded2, method = "bh")
salstats
salstats = salstats$res
cldList(comparison = salstats$Comparison, p.value = salstats$P.adj, threshold = 0.05)

DOplot <- ggplot(NRCunrounded2, aes(y = taxon_genus_name)) +
  #geom_density_ridges(aes(x = DO50), fill = "#002B55", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = DO20), fill = "#004A77", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = DO10), fill = "#006B9C", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = DO5), fill = "#1F8FC2", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = DO0), fill = "#56B4E9", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  #facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(75,340))+
  scale_y_discrete(limits = rev)+
  theme_minimal()+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 1))
DOplot

#Genus sal comparison
kruskal.test(DO0 ~ taxon_genus_name, data = NRCunrounded2)
#Post-hoc comparison of genera
DOstats <- dunnTest(DO0 ~ taxon_genus_name, data = NRCunrounded2, method = "bh")
DOstats
DOstats = DOstats$res
cldList(comparison = DOstats$Comparison, p.value = DOstats$P.adj, threshold = 0.05)


NOplot <- ggplot(NRCunrounded2, aes(y =taxon_genus_name)) +
  #geom_density_ridges(aes(x = NO50), fill = "#2F4858", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = NO20), fill = "#006B79", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = NO10), fill = "#007C7F", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = NO5), fill = "#008E7D", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = NO0), fill = "#009E73", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  #facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(-0.5,10))+
  scale_y_discrete(limits = rev)+
  theme_minimal()+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 1))
NOplot

#Genus NO comparison
kruskal.test(NO0 ~ taxon_genus_name, data = NRCunrounded2)
#Post-hoc comparison of genera
NOstats <- dunnTest(NO0 ~ taxon_genus_name, data = NRCunrounded2, method = "bh")
NOstats
NOstats = NOstats$res
cldList(comparison = NOstats$Comparison, p.value = NOstats$P.adj, threshold = 0.05)

POplot <- ggplot(NRCunrounded2, aes(y = taxon_genus_name)) +
  #geom_density_ridges(aes(x = PO50), fill = "#4A1700", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = PO20), fill = "#6A3700", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = PO10), fill = "#905800", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = PO5), fill = "#BA7A00", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = PO0), fill = "#E69F00", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  #facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(-0.1, 1))+
  scale_y_discrete(limits = rev)+
  theme_minimal()+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 1))
POplot

#Genus PO comparison
kruskal.test(PO0 ~ taxon_genus_name, data = NRCunrounded2)
#Post-hoc comparison of genera
POstats <- dunnTest(PO0 ~ taxon_genus_name, data = NRCunrounded2, method = "bh")
POstats
POstats = POstats$res
cldList(comparison = POstats$Comparison, p.value = POstats$P.adj, threshold = 0.05)

OSatplot <- ggplot(NRCunrounded2, aes(y = taxon_genus_name)) +
  #geom_density_ridges(aes(x = OSat50), fill = "#0072B2", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = OSat20), fill = "#3D8DD0", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = OSat10), fill = "#60AAEF", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = OSat5), fill = "#81C7FF", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = OSat0), fill = "#A2E5FF", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  #facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(85,116))+
  scale_y_discrete(limits = rev)+
  theme_minimal()+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 1))
OSatplot

#Genus OSat comparison
kruskal.test(OSat0 ~ taxon_genus_name, data = NRCunrounded2)
#Post-hoc comparison of genera
OSatstats <- dunnTest(OSat0 ~ taxon_genus_name, data = NRCunrounded2, method = "bh")
OSatstats
OSatstats = OSatstats$res
cldList(comparison = OSatstats$Comparison, p.value = OSatstats$P.adj, threshold = 0.05)

Silplot <- ggplot(NRCunrounded2, aes(y = taxon_genus_name)) +
  #geom_density_ridges(aes(x = Sil50), fill = "#342F00", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = Sil20), fill = "#575800", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = Sil10), fill = "#878400", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  geom_density_ridges(aes(x = Sil5), fill = "#BBB300", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01)+
  geom_density_ridges(aes(x = Sil0), fill = "#F0E442", alpha = 0.8, color = NA, scale = 1, size = 0.3, rel_min_height = 0.01) +
  #facet_grid(rows = vars(taxon_suborder_name), switch = "y", scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0,15))+
  scale_y_discrete(limits = rev)+
  theme_minimal()+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 1))
Silplot

#Genus OSat comparison
kruskal.test(Sil0 ~ taxon_genus_name, data = NRCunrounded2)
#Post-hoc comparison of genera
Silstats <- dunnTest(Sil0 ~ taxon_genus_name, data = NRCunrounded2, method = "bh")
Silstats
Silstats = Silstats$res
cldList(comparison = Silstats$Comparison, p.value = Silstats$P.adj, threshold = 0.05)

DOplot2 <- DOplot +                    # Modify margins of second ggplot2 plot
  theme(plot.margin = unit(c(5.5, 0, 5.5, 0), "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank()
        )
Silplot2 <- Silplot+                    # Modify margins of second ggplot2 plot
  theme(plot.margin = unit(c(5.5, 0, 5.5, 0), "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank()
  )
Tempplot2 <- tempplot +                    # Modify margins of second ggplot2 plot
  theme(plot.margin = unit(c(5.5, 0, 5.5, 0), "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank()
  )
OSatplot2 <- OSatplot +                    # Modify margins of second ggplot2 plot
  theme(plot.margin = unit(c(5.5, 0, 5.5, 0), "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank()
  )
NOplot2 <- NOplot +                    # Modify margins of second ggplot2 plot
  theme(plot.margin = unit(c(5.5, 0, 5.5, 0), "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank()
  )
POplot2 <- POplot +                    # Modify margins of second ggplot2 plot
  theme(plot.margin = unit(c(5.5, 0, 5.5, 0), "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank()
  )
grid.arrange(salplot, Tempplot2,  Silplot2, DOplot2, 
             OSatplot2, NOplot2, POplot2, nrow = 1, widths = c(1.5,1,1,1,1,1,1))
#This is informative

###
#Take these values and backplot them onto maps based on random forest models
###

library(readxl)
X1_Coordinates_List <- read_excel("C:/Users/colin/Downloads/1_Coordinates_List.xlsx", 
                                  col_types = c("text", "text", "text", 
                                                "text", "text", "numeric", "numeric", 
                                                "text", "text", "text", "text", "text", 
                                                "text", "text", "text", "text"))
View(X1_Coordinates_List)

#Take the values from NRCunrounded2 and back plot them onto original global data
#Base it on strongest predictors from randomforest
#Rhopilema
varImpPlot(stommodel)
abline(v=140)
abline(v=240)

#Sal0, 5, 10, 20
#PO20

StomolophusYes <- subset(StomOccur, Stomolophus == "Yes")
head(StomolophusYes)
#PO20
res<-quantile(StomolophusYes$DO20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
DO2 <- as.data.frame(DO2)
StomolophusDORangeOuter <- subset(DO2, DO20 >  res[2] & DO20 < res[6])
StomolophusDORangeInner <- subset(DO2, DO20 >  res[3] & DO20 < res[5])

#sal20
res<-quantile(StomolophusYes$sal20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
sal2 <- as.data.frame(sal2)
head(sal2)
StomolophusSalRangeOuter <- subset(sal2, sal20 >  res[2] & sal20 < res[6])
StomolophusSalRangeInner <- subset(sal2, sal20 >  res[3] & sal20 < res[5])

res<-quantile(StomolophusYes$sal10, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
StomolophusSalRangeOuter <- subset(StomolophusSalRangeOuter, sal20 >  res[2] & sal20 < res[6])
StomolophusSalRangeInner <- subset(StomolophusSalRangeInner, sal10 >  res[3] & sal10 < res[5])

#temp10
res<-quantile(StomolophusYes$temp10, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res

StomolophusTempRangeOuter <- subset(temp, temp10 >  res[2] & temp10 < res[6])
StomolophusTempRangeInner <- subset(temp, temp10 >  res[3] & temp10 < res[5])

#Combine
StomolophusSpaceOuter <- merge(StomolophusDORangeOuter, StomolophusSalRangeOuter, by = c("latitude", "longitude"))
StomolophusSpaceOuter <- merge(StomolophusSpaceOuter, StomolophusTempRangeOuter, by = c("latitude", "longitude"))

StomolophusSpaceOuter <- StomolophusSpaceOuter %>%
  group_by(latitude, longitude) %>%
  summarize(DO20 = mean(DO20), sal20 = mean(sal20), sal10 = mean(sal10), temp10 = mean(temp10))
head(StomolophusSpaceOuter)

StomolophusSpaceInner <- merge(StomolophusDORangeInner, StomolophusSalRangeInner, by = c("latitude", "longitude"))
StomolophusSpaceInner <- merge(StomolophusSpaceInner, StomolophusTempRangeInner, by = c("latitude", "longitude"))

StomolophusSpaceInner <- StomolophusSpaceInner %>%
  group_by(latitude, longitude) %>%
  summarize(DO20 = mean(DO20), sal20 = mean(sal20), sal10 = mean(sal10), temp = mean(temp10))
head(StomolophusSpaceInner)

StomolophusSpace <- rbind(StomolophusSpaceInner, StomolophusSpaceOuter)

StomoLit <- subset(X1_Coordinates_List, Genus == "Stomolophus")

StomolophusDistribution<- ggplot() +
  geom_density_2d_filled(aes(x = StomolophusSpace$longitude, y = StomolophusSpace$latitude, 
                             fill = stat(level)), h = 7, adjust = 1/4)+
  scico::scale_fill_scico_d(palette = "grayC")+
  geom_polygon(data = wrld_simpl[wrld_simpl@data$UN!="10",], 
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "grey95")+
  geom_point(aes(x = Stomolophus$longitude, y = Stomolophus$latitude), color = "#D81B60", shape =15, size = 3, alpha = 0.8)+
  geom_point(aes(x = StomoLit$longitude, y = StomoLit$latitude), color = "#1E88E5", shape = 16, size = 3, alpha = 0.8)+
  coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat))+
  theme_classic()+
  theme(legend.position = "none")
StomolophusDistribution

#Rhopilema
varImpPlot(rhopmodel)
abline(v=14)

#Sal0, 5, 10, 20
#PO20

RhopilemaYes <- subset(RhopOccur, Rhopilema == "Yes")
head(RhopilemaYes)
#PO20
res<-quantile(RhopilemaYes$PO20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
PO2 <- as.data.frame(PO2)
RhopilemaPORangeOuter <- subset(PO2, PO20 >  res[2] & PO20 < res[6])
RhopilemaPORangeInner <- subset(PO2, PO20 >  res[3] & PO20 < res[5])

#sal20
res<-quantile(RhopilemaYes$sal20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
sal2 <- as.data.frame(sal2)
head(sal2)
RhopilemaSalRangeOuter <- subset(sal2, sal20 >  res[2] & sal20 < res[6])
RhopilemaSalRangeInner <- subset(sal2, sal20 >  res[3] & sal20 < res[5])

res<-quantile(RhopilemaYes$sal10, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
RhopilemaSalRangeOuter <- subset(RhopilemaSalRangeOuter, sal20 >  res[2] & sal20 < res[6])
RhopilemaSalRangeInner <- subset(RhopilemaSalRangeInner, sal10 >  res[3] & sal10 < res[5])
res<-quantile(RhopilemaYes$sal5, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
RhopilemaSalRangeOuter <- subset(RhopilemaSalRangeOuter, sal5 >  res[2] & sal5 < res[6])
RhopilemaSalRangeInner <- subset(RhopilemaSalRangeInner, sal5 >  res[3] & sal5 < res[5])
res<-quantile(RhopilemaYes$sal0, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
RhopilemaSalRangeOuter <- subset(RhopilemaSalRangeOuter, sal0 >  res[2] & sal0 < res[6])
RhopilemaSalRangeInner <- subset(RhopilemaSalRangeInner, sal0 >  res[3] & sal0 < res[5])

#Combine
RhopilemaSpaceOuter <- merge(RhopilemaPORangeOuter, RhopilemaSalRangeOuter, by = c("latitude", "longitude"))

RhopilemaSpaceOuter <- RhopilemaSpaceOuter %>%
  group_by(latitude, longitude) %>%
  summarize(PO20 = mean(PO20), sal20 = mean(sal20), sal10 = mean(sal10), sal5 = mean(sal5), sal0 =mean(sal0))
head(RhopilemaSpaceOuter)

RhopilemaSpaceInner <- merge(RhopilemaPORangeInner, RhopilemaSalRangeInner, by = c("latitude", "longitude"))

RhopilemaSpaceInner <- RhopilemaSpaceInner %>%
  group_by(latitude, longitude) %>%
  summarize(PO20 = mean(PO20), sal20 = mean(sal20), sal10 = mean(sal10), sal5 = mean(sal5), sal0 =mean(sal0))
head(RhopilemaSpaceInner)

RhopilemaSpace <- rbind(RhopilemaSpaceInner, RhopilemaSpaceOuter)

RhopLit <- subset(X1_Coordinates_List, Genus == "Rhopilema")

RhopilemaDistribution<- ggplot() +
  geom_density_2d_filled(aes(x = RhopilemaSpace$longitude, y =RhopilemaSpace$latitude, 
                             fill = stat(level)), h = 7, adjust = 1/4)+
  scico::scale_fill_scico_d(palette = "grayC")+
  geom_polygon(data = wrld_simpl[wrld_simpl@data$UN!="10",], 
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "grey95")+
  geom_point(aes(x = Rhopilema$longitude, y = Rhopilema$latitude), color = "#D81B60", shape =15, size = 3, alpha = 0.8)+
  geom_point(aes(x = RhopLit$longitude, y = RhopLit$latitude), color = "#1E88E5", shape = 16, size = 3, alpha = 0.8)+
  coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat))+
  theme_classic()+
  theme(legend.position = "none")
RhopilemaDistribution

#Cassiopea
varImpPlot(cassiomodel)
abline(v=20)
abline(v=10)

#Sil20
#sal20, 0, 5, 10

CassiopeaYes <- subset(CassOccur, Cassiopea == "Yes")
head(CassiopeaYes)
#temp20
res<-quantile(CassiopeaYes$temp20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
temp <- as.data.frame(temp)
temp
CassiopeaTempRangeOuter <- subset(temp, temp20 >  res[2] & temp20 < res[6])
CassiopeaTempRangeInner <- subset(temp, temp20 >  res[3] & temp20 < res[5])

res<-quantile(CassiopeaYes$temp10, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
CassiopeaTempRangeOuter <- subset(CassiopeaTempRangeOuter, temp10 >  res[2] & temp10 < res[6])
CassiopeaTempRangeInner <- subset(CassiopeaTempRangeInner, temp10 >  res[3] & temp10 < res[5])
res<-quantile(CassiopeaYes$temp5, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
CassiopeaTempRangeOuter <- subset(CassiopeaTempRangeOuter, temp5 >  res[2] & temp5 < res[6])
CassiopeaTempRangeInner <- subset(CassiopeaTempRangeInner, temp5 >  res[3] & temp5 < res[5])
res<-quantile(CassiopeaYes$temp0, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
CassiopeaTempRangeOuter <- subset(CassiopeaTempRangeOuter, temp0 >  res[2] & temp0 < res[6])
CassiopeaTempRangeInner <- subset(CassiopeaTempRangeInner, temp0 >  res[3] & temp0 < res[5])

#sal20
res<-quantile(CassiopeaYes$sal20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
sal2 <- as.data.frame(sal2)
head(sal2)
CassiopeaSalRangeOuter <- subset(sal2, sal20 >  res[2] & sal20 < res[6])
CassiopeaSalRangeInner <- subset(sal2, sal20 >  res[3] & sal20 < res[5])

res<-quantile(CassiopeaYes$sal10, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
CassiopeaSalRangeOuter <- subset(CassiopeaSalRangeOuter, sal20 >  res[2] & sal20 < res[6])
CassiopeaSalRangeInner <- subset(CassiopeaSalRangeInner, sal10 >  res[3] & sal10 < res[5])
res<-quantile(CassiopeaYes$sal5, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
CassiopeaSalRangeOuter <- subset(CassiopeaSalRangeOuter, sal5 >  res[2] & sal5 < res[6])
CassiopeaSalRangeInner <- subset(CassiopeaSalRangeInner, sal5 >  res[3] & sal5 < res[5])
res<-quantile(CassiopeaYes$sal0, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
CassiopeaSalRangeOuter <- subset(CassiopeaSalRangeOuter, sal0 >  res[2] & sal0 < res[6])
CassiopeaSalRangeInner <- subset(CassiopeaSalRangeInner, sal0 >  res[3] & sal0 < res[5])

#Combine
CassiopeaSpaceOuter <- merge(CassiopeaTempRangeOuter, CassiopeaSalRangeOuter, by = c("latitude", "longitude"))

CassiopeaSpaceOuter <- CassiopeaSpaceOuter %>%
  group_by(latitude, longitude) %>%
  summarize(temp20 = mean(temp20), temp10 = mean(temp10), temp5 = mean(temp5), temp0 = mean(temp0), sal20 = mean(sal20), sal10 = mean(sal10), sal5 = mean(sal5), sal0 =mean(sal0))
head(CassiopeaSpaceOuter)

CassiopeaSpaceInner <- merge(CassiopeaTempRangeInner, CassiopeaSalRangeInner, by = c("latitude", "longitude"))

CassiopeaSpaceInner <- CassiopeaSpaceInner %>%
  group_by(latitude, longitude) %>%
  summarize(temp20 = mean(temp20), temp10 = mean(temp10), temp5 = mean(temp5), temp0 = mean(temp0), sal20 = mean(sal20), sal10 = mean(sal10), sal5 = mean(sal5), sal0 =mean(sal0))
head(CassiopeaSpaceInner)

CassiopeaSpace <- rbind(CassiopeaSpaceInner, CassiopeaSpaceOuter)

LitCassio <- subset(X1_Coordinates_List, Genus == "Cassiopea")

CassiopeaDistribution<- ggplot() +
  geom_density_2d_filled(aes(x =CassiopeaSpace$longitude, y =CassiopeaSpace$latitude, 
                             fill = stat(level)), h = 7, adjust = 1/4)+
  scico::scale_fill_scico_d(palette = "grayC")+
  geom_polygon(data = wrld_simpl[wrld_simpl@data$UN!="10",], 
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "grey95")+
  geom_point(aes(x = Cassiopea$longitude, y = Cassiopea$latitude), color = "#D81B60", shape =15, size = 3, alpha = 0.8)+
  geom_point(aes(x = LitCassio$longitude, y = LitCassio$latitude), color = "#1E88E5", shape = 16, size = 3, alpha = 0.8)+
  coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat))+
  theme_classic()+
  theme(legend.position = "none")
CassiopeaDistribution

#Catostylus
varImpPlot(catomodel)
abline(v=30)
#Sil20
#sal20, 0, 5, 10

CatostylusYes <- subset(CatoOccur, Catostylus == "Yes")

#Sil20
res<-quantile(CatostylusYes$Sil20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
Sil2 <- as.data.frame(Sil2)
CatostylusSilRangeOuter <- subset(Sil2, Sil20 >  res[2] & Sil20 < res[6])
CatostylusSilRangeInner <- subset(Sil2, Sil20 >  res[3] & Sil20 < res[5])

#sal20
res<-quantile(CatostylusYes$sal20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
sal2 <- as.data.frame(sal2)
head(sal2)
CatostylusSalRangeOuter <- subset(sal2, sal20 >  res[2] & sal20 < res[6])
CatostylusSalRangeInner <- subset(sal2, sal20 >  res[3] & sal20 < res[5])
#sal0
res<-quantile(CatostylusYes$sal0, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))

CatostylusSalRangeOuter <- subset(CatostylusSalRangeOuter, sal0 >  res[2] & sal0 < res[6])
CatostylusSalRangeInner <- subset(CatostylusSalRangeInner, sal0 >  res[3] & sal0 < res[5])
#sal5
res<-quantile(CatostylusYes$sal5, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))

CatostylusSalRangeOuter <- subset(CatostylusSalRangeOuter, sal5 >  res[2] & sal5 < res[6])
CatostylusSalRangeInner <- subset(CatostylusSalRangeInner, sal5 >  res[3] & sal5 < res[5])
#sal10
res<-quantile(CatostylusYes$sal10, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))

CatostylusSalRangeOuter <- subset(CatostylusSalRangeOuter, sal10 >  res[2] & sal10 < res[6])
CatostylusSalRangeInner <- subset(CatostylusSalRangeInner, sal10 >  res[3] & sal10 < res[5])
#Combine
CatostylusSpaceOuter <- merge(CatostylusSilRangeOuter, CatostylusSalRangeOuter, by = c("latitude", "longitude"))

CatostylusSpaceOuter <- CatostylusSpaceOuter %>%
  group_by(latitude, longitude) %>%
  summarize(Sil20 = mean(Sil20), sal20 = mean(sal20), sal0 = mean(sal0), sal5 = mean(sal5), sal10 = mean(sal10))
head(CatostylusSpaceOuter)

CatostylusSpaceInner <- merge(CatostylusSilRangeInner, CatostylusSalRangeInner, by = c("latitude", "longitude"))

CatostylusSpaceInner <- CatostylusSpaceInner %>%
  group_by(latitude, longitude) %>%
  summarize(Sil20 = mean(Sil20), sal20 = mean(sal20), sal0 = mean(sal0), sal5 = mean(sal5), sal10 = mean(sal10))
head(CatostylusSpaceInner)

CatoSpace <- rbind(CatostylusSpaceInner, CatostylusSpaceOuter)
LitCato <- subset(X1_Coordinates_List, Genus == "Catostylus")

CatostylusDistribution<- ggplot() +
  geom_density_2d_filled(aes(x =CatoSpace$longitude, y =CatoSpace$latitude, 
                             fill = stat(level)), h = 7, adjust = 1/4)+
  scico::scale_fill_scico_d(palette = "grayC")+
  geom_polygon(data = wrld_simpl[wrld_simpl@data$UN!="10",], 
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "grey95")+
  geom_point(aes(x = Catostylus$longitude, y = Catostylus$latitude), color = "#D81B60", shape =15, size = 3, alpha = 0.8)+
  geom_point(aes(x = LitCato$longitude, y = LitCato$latitude), color = "#1E88E5", shape = 16, size = 3, alpha = 0.8)+
  coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat))+
  theme_classic()+
  theme(legend.position = "none")
CatostylusDistribution

#Eupilema
varImpPlot(Eupilemamodel)
abline(v=2.5)
#sal0
#NO20

#sal0
EupilemaYes <- subset(EupilemaOccur, Eupilema == "Yes")

res<-quantile(EupilemaYes$sal0, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
sal2 <- as.data.frame(sal2)
head(sal2)
EupilemaSalRangeOuter <- subset(sal2, sal0 >  res[2] & sal0 < res[6])
EupilemaSalRangeInner <- subset(sal2, sal0 >  res[3] & sal0 < res[5])

#sal20
res<-quantile(EupilemaYes$sal20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))

EupilemaSalRangeOuter <- subset(EupilemaSalRangeOuter, sal20 >  res[2] & sal20 < res[6])
EupilemaSalRangeInner <- subset(EupilemaSalRangeInner, sal20 >  res[3] & sal20 < res[5])
#sal5
res<-quantile(EupilemaYes$sal5, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))

EupilemaSalRangeOuter <- subset(EupilemaSalRangeOuter, sal5 >  res[2] & sal5 < res[6])
EupilemaSalRangeInner <- subset(EupilemaSalRangeInner, sal5 >  res[3] & sal5 < res[5])
#sal10
res<-quantile(EupilemaYes$sal10, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))

EupilemaSalRangeOuter <- subset(EupilemaSalRangeOuter, sal10 >  res[2] & sal10 < res[6])
EupilemaSalRangeInner <- subset(EupilemaSalRangeInner, sal10 >  res[3] & sal10 < res[5])

#NO20
res<-quantile(EupilemaYes$NO20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
NO2 <- as.data.frame(NO2)
EupilemaNORangeOuter <- subset(NO2, NO20 >  res[2] & NO20 < res[6])
EupilemaNORangeInner <- subset(NO2, NO20 >  res[3] & NO20 < res[5])

#Combine
EupilemaSpaceOuter <- merge(EupilemaNORangeOuter, EupilemaSalRangeOuter, by = c("latitude", "longitude"))

EupilemaSpaceOuter <- EupilemaSpaceOuter %>%
  group_by(latitude, longitude) %>%
  summarize(NO20 = mean(NO20), sal0 = mean(sal0), sal5 = mean(sal5), sal10 = mean(sal10), sal20 = mean(sal20))
head(EupilemaSpaceOuter)

EupilemaSpaceInner <- merge(EupilemaNORangeInner, EupilemaSalRangeInner, by = c("latitude", "longitude"))

EupilemaSpaceInner <- EupilemaSpaceInner %>%
  group_by(latitude, longitude) %>%
  summarize(NO20 = mean(NO20), sal0 = mean(sal0), sal5 = mean(sal5), sal10 = mean(sal10), sal20 = mean(sal20))
head(EupilemaSpaceInner)

EupilemaSpace <- rbind(EupilemaSpaceInner, EupilemaSpaceOuter)

EupLit <- subset(X1_Coordinates_List, Genus == "Eupilema")

EupilemaDistribution<- ggplot() +
  geom_density_2d_filled(aes(x =EupilemaSpace$longitude, y =EupilemaSpace$latitude, 
                             fill = stat(level)), h = 7, adjust = 1/4)+
  scico::scale_fill_scico_d(palette = "grayC")+
  geom_polygon(data = wrld_simpl[wrld_simpl@data$UN!="10",], 
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "grey95")+
  geom_point(aes(x = Eupilema$longitude, y = Eupilema$latitude), color = "#D81B60", shape =15, size = 3, alpha = 0.8)+
  geom_point(aes(x = EupLit$longitude, y = EupLit$latitude), color = "#1E88E5", shape = 16, size = 3, alpha = 0.8)+
  coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat))+
  theme_classic()+
  theme(legend.position = "none")
EupilemaDistribution

#Lobonema
varImpPlot(Lobonemamodel)
abline(v=0.5)
abline(v=0.7)

#Temp5
head(LobonemaOccur)
nrow(LobonemaOccur)
LobonemaYes <- subset(LobonemaOccur, Lobonema == "Yes")
head(LobonemaYes)

res<-quantile(LobonemaYes$temp5, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
temp <- as.data.frame(temp)
LobonemaTempRangeOuter <- subset(temp, temp5 >  res[2] & temp5 < res[6])
LobonemaTempRangeInner <- subset(temp, temp5 >  res[3] & temp5 < res[5])
res<-quantile(LobonemaYes$temp10, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))

LobonemaTempRangeOuter <- subset(LobonemaTempRangeOuter, temp10 >  res[2] & temp10 < res[6])
LobonemaTempRangeInner <- subset(LobonemaTempRangeInner, temp10 >  res[3] & temp10 < res[5])
res<-quantile(LobonemaYes$temp0, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))

LobonemaTempRangeOuter <- subset(LobonemaTempRangeOuter, temp0 >  res[2] & temp0 < res[6])
LobonemaTempRangeInner <- subset(LobonemaTempRangeInner, temp0 >  res[3] & temp0 < res[5])
#Sil0
res<-quantile(LobonemaYes$Sil0, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
Sil2 <- as.data.frame(Sil2)
LobonemaSilRangeOuter <- subset(Sil2, Sil0 >  res[2] & Sil0 < res[6])
LobonemaSilRangeInner <- subset(Sil2, Sil0 >  res[3] & Sil0 < res[5])

#Only maintain coordinates that fulfill both requirements
LobonemaSpaceOuter <- merge(LobonemaTempRangeOuter, LobonemaSilRangeOuter, by = c("latitude", "longitude"))
head(LobonemaSpaceOuter)
LobonemaSpaceOuter <- LobonemaSpaceOuter %>%
  group_by(latitude, longitude) %>%
  summarize(temp5 = mean(temp5), Sil0 = mean(Sil0))

LobonemaSpaceInner <- merge(LobonemaTempRangeInner, LobonemaSilRangeInner, by = c("latitude", "longitude"))
head(LobonemaSpaceInner)
LobonemaSpaceInner <- LobonemaSpaceInner %>%
  group_by(latitude, longitude) %>%
  summarize(temp5 = mean(temp5), Sil0 = mean(Sil0))

LobonemaSpace <- rbind(LobonemaSpaceInner, LobonemaSpaceOuter)

LobLit <- subset(X1_Coordinates_List, Genus == "Lobonema")

LobonemaDistribution<- ggplot() +
  geom_density_2d_filled(aes(x =LobonemaSpace$longitude, y =LobonemaSpace$latitude, 
                             fill = stat(level)), h = 7, adjust = 1/4)+
  scico::scale_fill_scico_d(palette = "grayC")+
  geom_polygon(data = wrld_simpl[wrld_simpl@data$UN!="10",], 
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "grey95")+
  geom_point(aes(x = Lobonema$longitude, y = Lobonema$latitude), color = "#D81B60", shape =15, size = 3, alpha = 0.8)+
  geom_point(aes(x = LobLit$longitude, y = LobLit$latitude), color = "#1E88E5", shape = 16, size = 3, alpha = 0.8)+
  coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat))+
  theme_classic()+
  theme(legend.position = "none")
LobonemaDistribution
#Some hotspots of interest here...

#Lychnorhiza
varImpPlot(Lychnorhizamodel)
abline(v=5)
abline(v=3)

#Lychnorhiza = PO20
LychnorhizaYes <- subset(LychnorhizaOccur, Lychnorhiza == "Yes")

res<-quantile(LychnorhizaYes$PO20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
PO2 <- as.data.frame(PO2)
head(PO2)
LychnorhizaPORangeOuter <- subset(PO2, PO20 >  res[2] & PO20 < res[6])
LychnorhizaPORangeInner <- subset(PO2, PO20 >  res[3] & PO20 < res[5])

#Lychnorhiza = Sil10
res<-quantile(LychnorhizaYes$Sil10, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
Sil2 <- as.data.frame(Sil2)
LychnorhizaSilRangeOuter <- subset(Sil2, Sil10 >  res[2] & Sil10 < res[6])
LychnorhizaSilRangeInner <- subset(Sil2, Sil10 >  res[3] & Sil10 < res[5])

#Lychnorhiza = Sal0
res<-quantile(LychnorhizaYes$sal0, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
sal2 <- as.data.frame(sal2)
LychnorhizaSalRangeOuter <- subset(sal2, sal0 >  res[2] & sal0 < res[6])
LychnorhizaSalRangeInner <- subset(sal2, sal0 >  res[3] & sal0 < res[5])

#Combine
LychnoSpaceOuter <- merge(LychnorhizaPORangeOuter, LychnorhizaSilRangeOuter, by = c("latitude", "longitude"))
LychnoSpaceOuter <- merge(LychnoSpaceOuter, LychnorhizaSalRangeOuter, by = c("latitude", "longitude"))

LychnoSpaceOuter <- LychnoSpaceOuter %>%
  group_by(latitude, longitude) %>%
  summarize(PO20 = mean(PO20), Sil10 = mean(Sil10), sal0 = mean(sal0))
head(LychnoSpaceOuter)

LychnoSpaceInner <- merge(LychnorhizaPORangeInner, LychnorhizaSilRangeInner, by = c("latitude", "longitude"))
LychnoSpaceInner <- merge(LychnoSpaceInner, LychnorhizaSalRangeInner, by = c("latitude", "longitude"))

LychnoSpaceInner <- LychnoSpaceInner %>%
  group_by(latitude, longitude) %>%
  summarize(PO20 = mean(PO20), Sil10 = mean(Sil10), sal0 = mean(sal0))
head(LychnoSpaceInner)

LitLychno <- subset(X1_Coordinates_List, Genus == "Lychnorhiza")
head(LitLychno)

LychnoDistribution<- ggplot() +
  geom_density_2d_filled(aes(x =LychnoSpaceOuter$longitude, y =LychnoSpaceOuter$latitude, 
                             fill = stat(level)), h = 7, adjust = 1/4)+
  scico::scale_fill_scico_d(palette = "grayC")+
  geom_polygon(data = wrld_simpl[wrld_simpl@data$UN!="10",], 
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "grey95")+
  geom_point(aes(x = Lychnorhiza$longitude, y = Lychnorhiza$latitude), color = "#D81B60", shape =15, size = 3, alpha = 0.8)+
  geom_point(aes(x = LitLychno$longitude, y = LitLychno$latitude), color = "#1E88E5", shape = 16, size = 3, alpha = 0.8)+
  coord_cartesian(xlim = c(min_lon, max_lon))+
  coord_cartesian(ylim = c(min_lat, max_lat))+
  theme_classic()+
  theme(legend.position = "none")
LychnoDistribution


#Pseudorhiza
varImpPlot(pseudomodel)
abline(v=2.7)
#Sil20, Sal10

#Sil20,5,0,10
PseudoYes <- subset(PseudoOccur, Pseudorhiza == "Yes")
head(PseudoYes)

res<-quantile(PseudoYes$Sil20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
Sil2 <- as.data.frame(Sil2)
PseudoSilRangeOuter <- subset(Sil2, Sil20 >  res[2] & Sil20 < res[6])
PseudoSilRangeInner <- subset(Sil2, Sil20 >  res[3] & Sil20 < res[5])

res<-quantile(PseudoYes$Sil5, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
PseudoSilRangeOuter <- subset(PseudoSilRangeOuter, Sil5 >  res[2] & Sil5 < res[6])
PseudoSilRangeInner <- subset(PseudoSilRangeInner, Sil5 >  res[3] & Sil5 < res[5])

res<-quantile(PseudoYes$Sil0, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
PseudoSilRangeOuter <- subset(PseudoSilRangeOuter, Sil0 >  res[2] & Sil0 < res[6])
PseudoSilRangeInner <- subset(PseudoSilRangeInner, Sil0 >  res[3] & Sil0 < res[5])

res<-quantile(PseudoYes$Sil10, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
PseudoSilRangeOuter <- subset(PseudoSilRangeOuter, Sil10 >  res[2] & Sil10 < res[6])
PseudoSilRangeInner <- subset(PseudoSilRangeInner, Sil10 >  res[3] & Sil10 < res[5])

#sal 10,0,20,5
res<-quantile(PseudoYes$sal10, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
sal2 <- as.data.frame(sal2)
PseudoSalRangeOuter <- subset(sal2, sal10 >  res[2] & sal10 < res[6])
PseudoSalRangeInner <- subset(sal2, sal10 >  res[3] & sal10 < res[5])

res<-quantile(PseudoYes$sal0, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))

PseudoSalRangeOuter <- subset(PseudoSalRangeOuter, sal0 >  res[2] & sal0 < res[6])
PseudoSalRangeInner <- subset(PseudoSalRangeInner, sal0 >  res[3] & sal0 < res[5])

res<-quantile(PseudoYes$sal20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))

PseudoSalRangeOuter <- subset(PseudoSalRangeOuter, sal20 >  res[2] & sal20 < res[6])
PseudoSalRangeInner <- subset(PseudoSalRangeInner, sal20 >  res[3] & sal20 < res[5])

res<-quantile(PseudoYes$sal5, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))

PseudoSalRangeOuter <- subset(PseudoSalRangeOuter, sal5 >  res[2] & sal5 < res[6])
PseudoSalRangeInner <- subset(PseudoSalRangeInner, sal5 >  res[3] & sal5 < res[5])

PseudoSpaceOuter <- merge(PseudoSilRangeOuter, PseudoSalRangeOuter, by = c("latitude", "longitude"))
PseudoSpaceOuter <- PseudoSpaceOuter %>%
  group_by(latitude, longitude) %>%
  summarize(Sil20= mean(Sil20), Sil5= mean(Sil5),Sil0= mean(Sil0),Sil10= mean(Sil10),
            sal10 = mean(sal10),sal0 = mean(sal0),sal20 = mean(sal20),sal5 = mean(sal5))
head(PseudoSpaceOuter)

PseudoSpaceInner <- merge(PseudoSilRangeInner, PseudoSalRangeInner, by = c("latitude", "longitude"))
PseudoSpaceInner <- PseudoSpaceInner %>%
  group_by(latitude, longitude) %>%
  summarize(Sil20= mean(Sil20), Sil5= mean(Sil5),Sil0= mean(Sil0),Sil10= mean(Sil10),
            sal10 = mean(sal10),sal0 = mean(sal0),sal20 = mean(sal20),sal5 = mean(sal5))
head(PseudoSpaceInner)

PseudoSpace <- rbind(PseudoSpaceInner, PseudoSpaceOuter)

PseudoLit <- subset(X1_Coordinates_List, Genus == "Pseudorhiza")

PseudoDistribution<- ggplot() +
  geom_density_2d_filled(aes(x =PseudoSpace$longitude, y = PseudoSpace$latitude, 
                             fill = stat(level)), h = 7, adjust = 1/4)+
  scico::scale_fill_scico_d(palette = "grayC")+
  geom_polygon(data = wrld_simpl[wrld_simpl@data$UN!="10",], 
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "grey95")+
  geom_point(aes(x = Pseudorhiza$longitude, y = Pseudorhiza$latitude), color = "#D81B60", shape =15, size = 3, alpha = 0.8)+
  geom_point(aes(x = PseudoLit$longitude, y = PseudoLit$latitude), color = "#1E88E5", shape = 16, size = 3, alpha = 0.8)+
  coord_cartesian(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat))+
  theme_classic()+
  theme(legend.position = "none")
PseudoDistribution

###
#Rhizostoma
varImpPlot(rhizomodel)
abline(v=60)
abline(v=125)
#temp20, DO20

head(RhizOccur)
nrow(RhizOccur)
RhizYes <- subset(RhizOccur, Rhizostoma == "Yes")
head(RhizYes)

#temp20
res<-quantile(RhizYes$temp20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
RhizTempRangeOuter <- subset(temp, temp20 >  res[2] & temp20 < res[6])
nrow(RhizTempRangeOuter)
RhizTempRangeInner <- subset(temp, temp20 >  res[3] & temp20 < res[5])
nrow(RhizTempRangeInner)

#DO20
res<-quantile(RhizYes$DO20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
RhizDORangeOuter <- subset(DO2, DO20 > res[2] & DO20 < res[6])
nrow(RhizDORangeOuter)
RhizDORangeInner <- subset(DO2, DO20 > res[3] & DO20 < res[5])
nrow(RhizDORangeInner)

#OuterRange
RhizSpaceOuter <- merge(RhizTempRangeOuter, RhizDORangeOuter, by = c("latitude", "longitude"))
nrow(RhizSpaceOuter)
RhizSpaceOuter <- RhizSpaceOuter %>%
  group_by(latitude, longitude) %>%
  summarize(temp20 = mean(temp20), DO20= mean(DO20))

#InnerRange
RhizSpaceInner <- merge(RhizTempRangeInner, RhizDORangeInner, by = c("latitude", "longitude"))
RhizSpaceInner <- RhizSpaceInner %>%
  group_by(latitude, longitude) %>%
  summarize(temp20 = mean(temp20), DO20= mean(DO20))

RhizSpace <- rbind(RhizSpaceInner, RhizSpaceOuter)

RhizLit <- subset(X1_Coordinates_List, Genus == "Rhizostoma")

RhizDistribution<- ggplot() +
  geom_density_2d_filled(aes(x =RhizSpace$longitude, y =RhizSpace$latitude, 
                             fill = stat(level)), h = 7, adjust = 1/4)+
  scico::scale_fill_scico_d(palette = "grayC")+
  geom_polygon(data = wrld_simpl[wrld_simpl@data$UN!="10",], 
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "grey95")+
  geom_point(aes(x = Rhizostoma$longitude, y = Rhizostoma$latitude), color = "#D81B60", shape =15, size = 2, alpha = 0.8)+
  geom_point(aes(x = RhizLit$longitude, y = RhizLit$latitude), color = "#1E88E5", shape = 16, size = 2, alpha = 0.8)+
  coord_cartesian(xlim = c(min_lon, max_lon))+
  coord_cartesian(ylim = c(min_lat, max_lat))+
  theme_classic()+
  theme(legend.position = "none")
RhizDistribution

###
#Cotylorhiza
varImpPlot(cotylmodel)
abline(v=30)
abline(v=70)
#sal20

head(CotylOccur)
nrow(CotylOccur)
CotylYes <- subset(CotylOccur, Cotylorhiza == "Yes")
head(CotylYes)

#sal20
res<-quantile(CotylYes$sal20, probs = c(0, 0.05, 0.25,0.5,0.75, 0.95, 1))
res
CotylSalRangeOuter <- subset(sal2, sal20 >  res[2] & sal20 < res[6])
nrow(CotylSalRangeOuter)
CotylSalRangeInner <- subset(sal2, sal20 >  res[3] & sal20 < res[5])
nrow(CotylSalRangeInner)

CotylSpace <- rbind(CotylSalRangeInner, CotylSalRangeOuter)

CotylLit <- subset(X1_Coordinates_List, Genus == "Cotylorhiza")

CotylDistribution<- ggplot() +
  geom_density_2d_filled(aes(x =CotylSpace$longitude, y =CotylSpace$latitude, 
                             fill = stat(level)), h = 7, adjust = 1/4)+
  scico::scale_fill_scico_d(palette = "grayC")+
  geom_polygon(data = wrld_simpl[wrld_simpl@data$UN!="10",], 
               aes(x = long, y = lat, group = group),
               colour = "black", fill = "grey95")+
  geom_point(aes(x = Cotylorhiza$longitude, y = Cotylorhiza$latitude), color = "#D81B60", shape =15, size = 1, alpha = 1)+
  geom_point(aes(x = CotylLit$longitude, y = CotylLit$latitude), color = "#1E88E5", shape = 16, size = 1, alpha = 1)+
  coord_cartesian(xlim = c(min_lon, max_lon))+
  coord_cartesian(ylim = c(min_lat, max_lat))+
  theme_classic()+
  theme(legend.position = "none")
CotylDistribution

