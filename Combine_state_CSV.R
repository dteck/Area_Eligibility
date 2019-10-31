library(tigris)
library(leaflet)
library(leaflet.extras)
library(rgeos)
library(readxl)
library(RColorBrewer)
library(tidycensus)
library(rgdal)
library(rgeos)
library(htmltools)
library(spatstat)
library(spdep)
library(arrangements)
library(tidyverse)

setwd(choose.dir())
workdir<-getwd()
files <- list.files(path = workdir, pattern = "*.csv", full.names = T)

states<-data.frame()
for(i in files){
  print(i)
  a<-read.csv(header = TRUE, file = i)
  states<-rbind.data.frame(states,a)
  rm(a)
}

write.csv(states, "FY20_all.csv", row.names = FALSE)





