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

##
paletteBarPlot = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666')
barplot(t(envData), col = paletteBarPlot)





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

## Diversity index of litter

envData3 = envData[1:8]

simpsonDiv = vegan::diversity(envData3, index = "shannon") 
specRich = specnumber(specTable)


divRel = data.frame(simpsonDiv, specRich)
divRel = divRel[order(divRel$simpsonDiv),]


# lets make a regression 
lm1 = lm(specRich~simpsonDiv)

par(las=1)
plot(specRich~simpsonDiv,
     xlab ="Litter Complexity (Simpson)",
     ylab = "Number of Species",
     col = "#73b34a", 
     pch = 16, 
     cex = 2
     )
abline(lm1)
anova(lm1)
