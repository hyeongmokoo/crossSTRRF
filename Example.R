rm(list=ls(all=TRUE))    # Clean objects from workspace
library(maptools)
library(foreign)
library(spatstat)        # Key library for spatial point pattern analysis        
library(smacpod)         # Relative risk kernel densities based on statstat

source("Functions_CrossSTRRF")
s.area <- readShapePoly("Data\\city_of_plano_Buffer.shp") ##move to stkde_new folder
sample.dbf <- as.data.frame(read.dbf("Data\\Plano_Assualt_R.dbf"))
plano.owin <- as.owin(s.area)

t.size <- 1400 ##Temporal cell size.

###Bandwidth calculation with Scott Rules 
sub.sample <- subset(sample.dbf, (adjT >= 0)&(adjT < 182700))
sample.n <- nrow(sub.sample)
bw.t <- (sd(sub.sample$adjT)*(sample.n^(-1/6)))
sample.ppp <- ppp(sub.sample$x, sub.sample$y, window = plano.owin)
res.bw.scott <- bw.scott(sample.ppp)

##the cross-STRRF calculation
res.cross.strrf <- cross.strrf(sample.dbf$x, sample.dbf$y, sample.dbf$adjT, res.bw.scott, bw.t, 1400, s.area)

##the cross-STRRF calculation with a significance test
res.cross.strrf <- cross.strrf(sample.dbf$x, sample.dbf$y, sample.dbf$adjT, res.bw.scott, bw.t, 1400, s.area, nsim=999)

