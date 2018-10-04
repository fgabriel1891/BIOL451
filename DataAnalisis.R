library(vegan)
library(stats)
####################################################################
# Load data
# Add weight data
dataPeso = read.csv("DATA/Weights.csv", header = T)
head(dataPeso)
tail(dataPeso)
# Add canopy data 
canopy = read.csv("DATA/CanopyApp.csv", header = T)
head(canopy)
canopy$X = NULL
head(canopy)
sort(canopy$Site)
# Add species data
species = read.csv("DATA/BiodiversityData2OpenRefine.csv", header = T)
species2 = read.csv("DATA/BiodiversityData2.csv", header = T)


sort(unique(species$Morfospecies))
head(species)
###############

# Summarize data
summary(dataPeso) # Resumenes 
summary(canopy)
summary(species)

# Exploratory analysis 
plot(dataPeso$WetWeight~dataPeso$LitterType)
plot(dataPeso$DryWeight~dataPeso$LitterType)
plot(canopy$Site, canopy$Percentage)


# Calculate total weights per sample 

WetWeight = aggregate(dataPeso$WetWeight, by =  list(dataPeso$Sample), FUN = sum)
DryWeight = aggregate(dataPeso$DryWeight, list(dataPeso$Sample), FUN = sum)

# Make a new summary dataPeso 

Weights =  cbind(WetWeight, DryWeight$x)
names(Weights)
names(Weights) = c("Sample", "WetWeight", "DryWeight")

# Infer humidity content per sample as the difference of wet vs dry weight
Weights$WeightDiff = Weights$WetWeight - Weights$DryWeight
# Calculate how the sample deviates from the mean (to standardize a measurement of relative humidity)
Weights$HumStd = Weights$WeightDiff / max(Weights$WeightDiff)
# Add canopy data 
Weights$canopy = canopy$Percentage[match(Weights$Sample,canopy$Site)]


## Match  and calculate relative weights back to original dataPeso
dataPeso$pesoTotDry = Weights$DryWeight[match(dataPeso$Sample, Weights$Sample)]
# Divide the dry weight vs the total dry weight 
dataPeso$pesoDryProp = dataPeso$DryWeight/dataPeso$pesoTotDry


## Calculate environmental variable Matrix with proportional dry weight (make dummy variables)
envData = data.frame(reshape2::acast(data = dataPeso, Sample~LitterType, sum, value.var = "pesoDryProp"))

## Add humidity and canopy op data to envmatrix
envData$HumStd = Weights$HumStd # Is better to use std humidity or raw humidity
envData$Canopy = Weights$canopy 
envData2 = envData
envData$Canopy = NULL

### Looking at the species occurrences 

## Make species occurrence matrix
names = names(species)
names(species) = c("Sample", names[-1])
names(species)
specTable = reshape2::acast(Sample~Morfospecies, sum, value.var = c("Abundance") , data = species)

dim(specTable) # dimensions

##
pdf("figs/Figure1.pdf")
par(las = 1, 
    mar=c(5.1, 5.1, 4.1, 8.1), # mar argument of par() means the margins 
    xpd=TRUE) # parametros 
paletteBarPlot = c("#c1883f", "#cc554f", "#80a33f", "#95c3ff", "#ccb54e", "#92e473", "#7a7f55", "#1f2a0a")
barplot(t(envData), # data to plot
        col = paletteBarPlot, # "col" argument means a vector of colors with the same length of the rownames from the data 
        xlab = "site",
        ylab = "litter dry relative abundance",
        cex.lab = 1.5) 
legend("topright", 
       rownames(t(envData)), # rownames =  
       fill = paletteBarPlot, # fill = 
       inset=c(-0.32,0))
dev.off() #


#############

## Diversity index of litter

envData3 = envData[1:8]

simpsonDiv = vegan::diversity(envData3, index = "shannon") 
specRich = specnumber(specTable)

relData = data.frame( simpsonDiv,specRich)
relData$abundance = rowSums(specTable)


divRel = data.frame(simpsonDiv, specRich)
divRel = divRel[order(divRel$simpsonDiv),]


# lets make a regression 
lm1 = lm(specRich~simpsonDiv)
lm2 = lm(relData$abundance~relData$simpsonDiv)

pdf("figs/litterComplexityRel.pdf")
par(las = 1, 
    mar=c(5.1, 5.1, 2.1, 2.1))
plot(relData$specRich~relData$simpsonDiv,
     xlab ="litter complexity (Simpson)",
     ylab = "number of species",
     col = "#73b34a", 
     pch = 16, 
     cex.lab = 2,
     cex = log(relData$abundance)
     )
points(relData$specRich~relData$simpsonDiv,
       cex = log(relData$abundance))
dev.off()

summary(lm1) # species richness
summary(lm2) # species abundance

## Scatterplot of humidity vs diversity 
relData$hum = envData$HumStd
# make regressions 
humSpeRich = lm(relData$specRich~relData$hum)
summary(humSpeRich)
humAbudRel = lm(relData$abundance~relData$hum)
summary(humAbudRel)

pdf("figs/humidityRel.pdf")
par(las = 1, 
    mar=c(5.1, 5.1, 2.1, 2.1))
plot(relData$specRich~relData$hum,
     xlab = "relative humidity",
     ylab = "number of species",
     col = "#73b34a", 
     pch = 16, 
     cex.lab = 2,
     cex = log(relData$abundance)
)
points(relData$specRich~relData$hum,
       cex = log(relData$abundance))
dev.off()


## Scatterplot of Canopy vs diversity 
relData$canopy = canopy[match(rownames(relData),canopy$Site),]$Percentage

CanoSpRic = lm(relData$specRich~relData$canopy)
summary(CanoSpRic)
CanoSpAbun = lm(relData$abundance~relData$canopy)
summary(CanoSpAbun)
# make regression 
pdf("figs/Canopy.pdf")
par(las = 1, 
    mar=c(5.1, 5.1, 2.1, 2.1))
plot(relData$specRich~relData$canopy,
     xlab = "relative humidity",
     ylab = "number of species",
     col = "#73b34a", 
     pch = 16, 
     cex.lab = 2,
     cex = log(relData$abundance)
)
points(relData$specRich~relData$canopy,
       cex = log(relData$abundance))
dev.off()

##################

## Ordinations 
## Make constrained ordination 
ccaMatrix = vegan::rda(specTable, envData)
plot(ccaMatrix)
## Make permutation test
anova(ccaMatrix)
# Observe loadings and asssociated data
summary(ccaMatrix)

## What if we transform the abundance data so we can see variations better. 
# Logarithm (check Log1p function)
specTableLog = log1p(specTable)
## Make constrained ordination 
ccaMatrix2 = vegan::rda(specTableLog, envData)
plot(ccaMatrix2)
anova(ccaMatrix2)
summary(ccaMatrix2)


# erasing the p13 

envData2 = envData2[!rownames(envData2) == "p13",]
specTable2 = specTable[!rownames(specTable) == "p13",]


ccaMatrix3 = vegan::rda(log1p(specTable2), envData2)
plot(ccaMatrix3)
anova(ccaMatrix3)
summary(ccaMatrix3)

#############

