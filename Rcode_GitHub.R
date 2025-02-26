library(ggplot2)
library(purrr)
library(tidyr)
library(tidyverse)
library(terra)
library(sf)
library(data.table)
library(grid)
library(corrplot)
library(bibtex)
library(readxl)
library(leaflegend)
library(kableExtra)
library(ape)
library(vegan)
library(corrr)
library(factoextra)
library(gridExtra)
library(ggfortify)
library(sdmTMB)
library(landscapemetrics)
library(FNN)
library(geodata)
library(predicts)
library(ENMTools)
library(jtools)
library(RColorBrewer)
library(flextable)
library(dismo)
library(gbm)
library(caret)
library(pdp)

#### Raster data ####

# Details on sources for all raster layers found in Table 1 of the journal article
full_stack = rast("~/full_stack.tif")

#### create points ####
set.seed(3888)

# create random sample of points within study area
bgsamp1 = st_sample(surveyarea, 10000, "random") 

# pull out coordinate columns, x (longitude) first, then y (latitude) from data
species_obs_car <- species_obs_car %>% dplyr::mutate(lon = st_coordinates(.)[,1], 
                                                     lat = st_coordinates(.)[,2])
presence = species_obs_car[, c("lon", "lat")]
presence$pa = 1

#convert background data to df
absence = st_as_sf(bgsamp1)
absence <- absence %>% dplyr::mutate(lon = st_coordinates(.)[,1],
                                     lat = st_coordinates(.)[,2])
absence$pa = 0
colnames(absence)[1] = "geometry" 
st_geometry(absence) = "geometry"

#join data into single df
all_points = rbind(presence, absence)

## extract to points
extract = terra::extract(full_stack, all_points[,c("lon", "lat")], ID = FALSE)
points1 = cbind(all_points, extract)
#identify columns that are lon and lat
drop_cols = which(colnames(points1) %in% c("lon", "lat"))
#remove geographic coordinates from df
points1 = points1[, -drop_cols]
points1$geometry = NULL
points1 = na.omit(points1)

#### training and testing data ####
fold = folds(points1, k = 5, by = points1$pa)
table(fold)  
testing = points1[fold == 1,]  
training = points1[fold != 1,]

#### models ####

# model repeated for all combinations of tree complexity, learning rate, and bag fraction listed in Table A2 of the journal article Supplemental data.
brt.tc5.lr005.75 = gbm.step(data = training, gbm.x = c("MCMT", "MSP", "bFFP", "eFFP", "elevation", "ruggedness", "ndvi", "ndmi", "nleaf", "Water", "ThicketSwamp", "ConiferousSwamp", "DeciduousSwamp", "OpenBog", "TreedBog", "OpenFen", "TreedFen", "Coniferous", "Sparse", "Deciduous", "NonTreed", "Anthropogenic", "Disturbance", "snowdepth"), gbm.y = 1, family = "bernoulli", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.75)
save(brt.tc5.lr005.75, file = "~/brt.tc5.lr005.75.RData")

#### summary table ####

brt.model.summaries = as.data.frame(matrix(nrow = 13, ncol = 1))
colnames(brt.model.summaries) = "flight.brt.tc5.lr005.75"
rownames(brt.model.summaries) = c("Tree complexity", "Learning rate", "Bag fraction", "Number of trees", "Total mean deviance", "Residual mean deviance", "RMSE", "Rsquared", "MAE", "tPos", "tNeg", "fPos", "fNeg")

# code below repeated for all combinations of tree complexity, learning rate, and bag fraction listed in Table A2 of the journal article Supplemental data.
flight.brt.tc5.lr005.75.e = dismo::evaluate(subset(points1, pa == "1"), subset(points1, pa == "0"), flight.brt.tc5.lr005.75, type = "response")

brt.model.summaries$flight.brt.tc5.lr005.75 = c(flight.brt.tc5.lr005.75$gbm.call$tree.complexity, flight.brt.tc5.lr005.75$gbm.call$learning.rate, flight.brt.tc5.lr005.75$gbm.call$bag.fraction, flight.brt.tc5.lr005.75$n.trees, flight.brt.tc5.lr005.75$self.statistics$mean.null, flight.brt.tc5.lr005.75$self.statistics$mean.resid, (postResample(pred = predict(flight.brt.tc5.lr005.75, testing), obs = testing$pa)), mean(flight.brt.tc5.lr005.e@TPR), mean(flight.brt.tc5.lr005.e@TNR), mean(flight.brt.tc5.lr005.e@FPR), mean(flight.brt.tc5.lr005.e@FNR))

### final evaluation table
brt.model.summaries = as.data.frame(t(brt.model.summaries))

brt.model.summaries2 = brt.model.summaries

brt.model.summaries$modelID = rownames(brt.model.summaries)

get_sens <- function(data){
  TP <- data$tPos
  FN <- data$fNeg
  as.vector(TP / (TP + FN))
}

brt.model.summaries$sens = get_sens(brt.model.summaries)

get_spec <- function(data){
  TN <- data$tNeg
  FP <- data$fNeg
  as.vector(TN / (TN + FP))
}

brt.model.summaries$spec = get_spec(brt.model.summaries)

#### ROC ####

topmodel = brt.tc5.lr005.75
eval = dismo::evaluate(subset(points1, pa == "1"), subset(points1, pa == "0"), topmodel, type = "response")
plot(eval, "ROC")
summary.gbm(topmodel)

#### partial dependence plots ####
gbm.plot(topmodel, common.scale = T, n.plots = 23, smooth = T,  rug = F, plot.layout = c(3,3), write.title = F)

#### interactions ####
find.int = gbm.interactions(topmodel)
find.int$rank.list
gbm.perspec2(topmodel, 2, 1, theta = 225, phi = 25, perspective = T)


#### threshold determination and spatial plotting ####
Method = factor('electric', levels = levels(training$Method))
add = data.frame(Method)
p = predict(full_stack, topmodel, const = add, n.trees = topmodel$gbm.call$best.trees, type = "response")

threshold = threshold(dismo::evaluate(subset(points1, pa == "1"), subset(points1, pa == "0"), topmodel, type = "response"), stat = "spec_sens")
plot(p > threshold)