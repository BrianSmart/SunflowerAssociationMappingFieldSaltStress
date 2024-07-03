library(dplyr)
library(ggplot2)
library(scales)
library(magrittr)
library(gstat)
library(sp)

#Load training and test data. Training data is the soil samples and test data will be used to estimate salinity for each plot
train <- read.csv("krig_train_data.csv", stringsAsFactors = F)
test <- read.csv("krig_test_data.csv", stringsAsFactors = F)

#Kriging is performed for 2 years at 2 soil depths, 0 to 6 cm and 6 to 24 cm. 

#### 2016 ####
# Subset for 2016
train <- train[train$Year == 2016,]
test <- test[test$Year == 2016,]

coordinates(train) <- ~x+y
coordinates(test) <- ~x+y
# Create variogram using pre-determined values for cutoff and width
vario2016_0t6 <- variogram(salinity0_6 ~ x+y, cutoff = 350, data=train, width = 36)
vario2016_6t24 <- variogram(salinity6_24 ~ x+y, cutoff = 350, data=train, width = 34)
# Fit model using pre-determined values for model, nugget, and range
fit2016_0t6 <- fit.variogram(vario2016_0t6, model = vgm(psill = NA, model = 'Pen', nugget = 0, range = 225))
fit2016_6t24 <- fit.variogram(vario2016_6t24, model = vgm(psill = NA, model = 'Pen', nugget = 0, range = 225))
# Estimate salinity for individual plots
krig_0t6 <- krige(salinity0_6 ~ x+y, train, test, model = fit2016_0t6)
krig_0t6 <- as.data.frame(krig_0t6)
krig_6t24 <- krige(salinity6_24 ~ x+y, train, test, model = fit2016_6t24)
krig_6t24 <- as.data.frame(krig_6t24)

#### 2017 ####
train <- read.csv("krig_train_data.csv", stringsAsFactors = F)
test <- read.csv("krig_test_data.csv", stringsAsFactors = F)

#Two subsets are required for 2017, one for Rep 3 and the other Rep 1 & 2
trainR3 <- train[train$Year == '2017W',]
testR3 <- test[test$Year == '2017W',]
trainR1_2 <- train[train$Year == '2017H',]
testR1_2 <- test[test$Year == '2017H',]

coordinates(trainR3) <- ~x+y
coordinates(testR3) <- ~x+y
coordinates(trainR1_2) <- ~x+y
coordinates(testR1_2) <- ~x+y

# Create variogram using pre-determined values for cutoff and width
vario2017R3_0t6 <- variogram(salinity0_6 ~ x+y, cutoff = 300, data=trainR3, width = 22)
vario2017R3_6t24 <- variogram(salinity6_24 ~ x+y, cutoff = 300, data=trainR3, width = 24)
vario2017R1_2_0t6 <- variogram(salinity0_6 ~ x+y, cutoff = 250, data=trainR1_2, width = 14)
vario2017R1_2_6t24 <- variogram(salinity6_24 ~ x+y, cutoff = 350, data=trainR1_2, width = 12)
# Fit model using pre-determined values for model, nugget, and range
fit2017R3_0t6 <- fit.variogram(vario2017R3_0t6, model = vgm(psill = NA, model = 'Pen', nugget = 0, range = 150))
fit2017R3_6t24 <- fit.variogram(vario2017R3_6t24, model = vgm(psill = NA, model = 'Pen', nugget = 0, range = 175))
fit2017R1_2_0t6 <- fit.variogram(vario2017R1_2_0t6, model = vgm(psill = NA, model = 'Cir', nugget = 0, range = 175))
fit2017R1_2_6t24 <- fit.variogram(vario2017R1_2_6t24, model = vgm(psill = NA, model = 'Cir', nugget = 0, range = 200))
# Estimate salinity for individual plots
krig2017_R3_0t6 <- krige(salinity0_6 ~ x+y, trainR3, testR3, model = fit2017R3_0t6)
krig2017_R3_0t6 <- as.data.frame(krig2017_R3_0t6)
krig2017_R3_6t24 <- krige(salinity6_24 ~ x+y, trainR3, testR3, model = fit2017R3_6t24)
krig2017_R3_6t24 <- as.data.frame(krig2017_R3_6t24)
krig2017_R1_2_0t6 <- krige(salinity0_6 ~ x+y, trainR1_2, testR1_2, model = fit2017R1_2_0t6)
krig2017_R1_2_0t6 <- as.data.frame(krig2017_R1_2_0t6)
krig2017_R1_2_6t24 <- krige(salinity6_24 ~ x+y, trainR1_2, testR1_2, model = fit2017R1_2_6t24)
krig2017_R1_2_6t24 <- as.data.frame(krig2017_R1_2_6t24)

