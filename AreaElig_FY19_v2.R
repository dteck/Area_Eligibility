############################
# Author: Mark Richards    #
############################

#load libraries--
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
library(rgeos)

#create direcotories and download/unzip files
setwd(choose.dir()) #asks user to choose Working dir
wd<-getwd() #gets copy of WD path
USDA<-paste(wd,"/USDA", sep="") #path to save USDA file
dir.create(file.path(USDA)) #create dir
setwd(file.path(USDA)) 
  download.file('http://data-cacfp-sfsp.opendata.arcgis.com/datasets/7a39c0ae571149c5a09c5790f0ff99f9_0.zip', method = 'auto',
                destfile ='FY19.zip', quiet = TRUE) #download file from USDA
unzip("./FY19.zip") #unzip file
#------

#read in shapefile and remove extra data
USDA_SHP<-readOGR('FY19_CACFP_SFSP_Census_Eligibility_ACS2012_2016.shp')
USDA_SHP<-USDA_SHP[,c(3,7,11,12,18,19,20,22,23)]
USDA_SHP$STATEFP<-as.numeric(as.character(USDA_SHP$STATEFP))
#------

#Separate shp into states by FIPS code
StateData<-list()
for (i in seq(1,56,1)){ 
  if (i %in% c(3,7,14,43,52)){
    print(paste('Skip non-state:', i))
  } else {
    print(i)
    StateData[[i]]<-subset.data.frame(USDA_SHP,USDA_SHP$STATEFP == i)
  }
}
setwd(wd)
saveRDS(StateData, "StateDataUSDA.rds")
rm(USDA_SHP)
#------

#Build Lists of neighbors for each blockgroup in a state
Neighbors<-list()
for (i in seq(1,56,1)){
  if (i %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", i))
  } else{
    print(paste("Finding Neighbors for:", i))
    Neighbors[[i]]<-poly2nb(StateData[[i]], queen = TRUE)
  }
}
saveRDS(Neighbors, "NeighborsUSDA.rds")
#------


#extract state data from USDA shape files
USDAData<-list()
for (i in seq(1,56,1)){
  if (i %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", i))
  } else{
    print(paste("Extracting Data for:", i))
    USDAData[[i]]<-data.frame(StateData[[i]])
    USDAData[[i]]$TotPovUniv<-as.numeric(as.character(USDAData[[i]]$TotPovUniv))
    USDAData[[i]]$Num18BG<-as.numeric(as.character(USDAData[[i]]$Num18BG))
    USDAData[[i]]$TotPovUn_1<-as.numeric(as.character(USDAData[[i]]$TotPovUn_1))
    USDAData[[i]]$Num12BG<-as.numeric(as.character(USDAData[[i]]$Num12BG))
  }
}
saveRDS(USDAData, "USDAData.rds")
rm(StateData)
#------



#function to calc combos for under 18
USDA_calc<- function(StateFIPS){ #function calc combos for under 18
  dfcombo<-data.frame() #build df to hold results
  for(i in seq(1,length(USDAData[[StateFIPS]]$GEOID),1)){ #run through all group blocks
    if (USDAData[[StateFIPS]]$ELIGFY19[i]=="No"){
      #print(paste(i, length(neighbor[[i]])))
      if (length(Neighbors[[StateFIPS]][[i]]) >1){ #if the # of neighboring BGs is more than 1
        comb<-arrangements::combinations(Neighbors[[StateFIPS]][[i]],2, replace = FALSE)
        for (n in seq(1,length(comb[,1]),1)){
          print(paste('FIPS:',StateFIPS,'Block Group:', i, "of", length(USDAData[[StateFIPS]]$GEOID),"Calculation", n, "of", length(comb[,1])))
          dfcomb2<-data.frame(
            GEOID=USDAData[[StateFIPS]]$GEOID[i], #bg in question id
            ID18=USDAData[[StateFIPS]]$Num18BG[i], #bg in question Identified under 18
            Tot18=USDAData[[StateFIPS]]$TotPovUniv[i], #bg in question Total under 18
            ID12=USDAData[[StateFIPS]]$Num12BG[i], #bg in question Identified under 12
            Tot12=USDAData[[StateFIPS]]$TotPovUn_1[i], #bg in question Total under 12
            firstGEOID18=USDAData[[StateFIPS]]$GEOID[comb[n,1]], #first calc block id
            B1_ID18=USDAData[[StateFIPS]]$Num18BG[comb[n,1]], #first calc block Identified under 18
            B1_Tot18=USDAData[[StateFIPS]]$TotPovUniv[comb[n,1]], #first calc block Total under 18
            B1_ID12=USDAData[[StateFIPS]]$Num12BG[comb[n,1]], #first calc block Identified under 12
            B1_Tot12=USDAData[[StateFIPS]]$TotPovUn_1[comb[n,1]],#first calc block Total under 12
            secondGEOID18=USDAData[[StateFIPS]]$GEOID[comb[n,2]], #second calc block id
            B2_ID18=USDAData[[StateFIPS]]$Num18BG[comb[n,2]], #second calc block Identified under 18
            B2_Tot18=USDAData[[StateFIPS]]$TotPovUniv[comb[n,2]], #second calc block Total under 18
            B2_ID12=USDAData[[StateFIPS]]$Num12BG[comb[n,2]], #second calc block Identified under 12
            B2_Tot12=USDAData[[StateFIPS]]$TotPovUn_1[comb[n,2]])#second calc block Total under 12
          dfcomb2$weightnum18<-dfcomb2$ID18+dfcomb2$B1_ID18+dfcomb2$B2_ID18 #weight identified under 18
          dfcomb2$weightpov18<-dfcomb2$Tot18+dfcomb2$B1_Tot18+dfcomb2$B2_Tot18 #weight total under 18
          dfcomb2$weightPerc18<-(dfcomb2$weightnum18/dfcomb2$weightpov18)*100 #weight percent under 18
          dfcomb2$weightnum12<-dfcomb2$ID12+dfcomb2$B1_ID12+dfcomb2$B2_ID12 #weight identified under 12
          dfcomb2$weightpov12<-dfcomb2$Tot12+dfcomb2$B1_Tot12+dfcomb2$B2_Tot12 #weight total under 12
          dfcomb2$weightPerc12<-(dfcomb2$weightnum12/dfcomb2$weightpov12)*100 # weight percent under 12
          dfcombo<-rbind(dfcombo, dfcomb2)
          rm(dfcomb2)
        }
      } else { #if the # of neighboring BGs is not greater than 1
        if (Neighbors[[StateFIPS]][[i]][1]==0){
          dfcomb2<-data.frame(GEOID=USDAData[[StateFIPS]]$GEOID[i],
                              ID18=USDAData[[StateFIPS]]$Num18BG[i],
                              Tot18=USDAData[[StateFIPS]]$TotPovUniv[i],
                              ID12=USDAData[[StateFIPS]]$Num12BG[i],
                              Tot12=USDAData[[StateFIPS]]$TotPovUn_1[i],
                              firstGEOID18=as.factor(0),
                              B1_ID18=0,
                              B1_Tot18=0,
                              B1_ID12=0,
                              B1_Tot12=0,
                              secondGEOID18=as.factor(0),
                              B2_ID18=0,
                              B2_Tot18=0,
                              B2_ID12=0,
                              B2_Tot12=0,
                              weightnum18=0,
                              weightpov18=0,
                              weightPerc18=0,
                              weightnum12=0,
                              weightpov12=0,
                              weightPerc12=0)
          dfcombo<-rbind(dfcombo, dfcomb2)
          rm(dfcomb2)
        } else {
          dfcomb2<-data.frame(
            GEOID=USDAData[[StateFIPS]]$GEOID[i],
            ID18=USDAData[[StateFIPS]]$Num18BG[i],
            Tot18=USDAData[[StateFIPS]]$TotPovUniv[i],
            ID12=USDAData[[StateFIPS]]$Num12BG[i],
            Tot12=USDAData[[StateFIPS]]$TotPovUn_1[i],
            firstGEOID18=USDAData[[StateFIPS]]$GEOID[Neighbors[[StateFIPS]][[i]][1]], #neighbor[[i]][1]
            B1_ID18=USDAData[[StateFIPS]]$Num18BG[Neighbors[[StateFIPS]][[i]][1]], 
            B1_Tot18=USDAData[[StateFIPS]]$TotPovUniv[Neighbors[[StateFIPS]][[i]][1]], 
            B1_ID12=USDAData[[StateFIPS]]$Num12BG[Neighbors[[StateFIPS]][[i]][1]],
            B1_Tot12=USDAData[[StateFIPS]]$TotPovUn_1[Neighbors[[StateFIPS]][[i]][1]],
            secondGEOID18=as.factor(0),
            B2_ID18=0,
            B2_Tot18=0,
            B2_ID12=0,
            B2_Tot12=0)
          dfcomb2$weightnum18<-dfcomb2$ID18+dfcomb2$B1_ID18 
          dfcomb2$weightpov18<-dfcomb2$Tot18+dfcomb2$B1_Tot18
          dfcomb2$weightPerc18<-(dfcomb2$weightnum18/dfcomb2$weightpov18)*100 
          dfcomb2$weightnum12<-dfcomb2$ID12+dfcomb2$B1_ID12
          dfcomb2$weightpov12<-dfcomb2$Tot12+dfcomb2$B1_Tot12
          dfcomb2$weightPerc12<-(dfcomb2$weightnum12/dfcomb2$weightpov12)*100
          dfcombo<-rbind(dfcombo, dfcomb2)
          rm(dfcomb2)
        }
      }
    }
    print(paste('FIPS:',StateFIPS,'Block Group:', i, "of", length(USDAData[[StateFIPS]]$GEOID)))
  }
  return(dfcombo)
}



StateDataAll<<-list()
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    print(paste("Calculating Averages for:", r))
    StateDataAll[[r]]<-USDA_calc(r)
    saveRDS(StateDataAll[[r]],paste("NeighborCacl_",r,".rds",sep = ""))
  }
}
saveRDS(StateDataAll, "StateDataAll.rds")

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------


#keep only one result above 50%
Over50<-list()
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    print(paste('Keeping one result per GEOID for:',r))
    Over50[[r]]<-subset.data.frame(StateDataAll[[r]], StateDataAll[[r]]$weightPerc18 >= 50 | StateDataAll[[r]]$weightPerc12 >=50)
    Over50[[r]]<-Over50[[r]][!duplicated(Over50[[r]]$GEOID),]
  }
}
saveRDS(Over50, "Over50.rds")

#strip original geoid from 0ver 50 list
Over50strip<-list()
for (i in seq(1,56,1)){ 
  if (i %in% c(3,7,14,43,52)){
    print(paste('Skip non-state:', i))
  } else {
    print(i)
    Over50strip[[i]]<-Over50[[i]][,c(1,6:21)]
  }
}


#Create list of state outlines
StatesSHP<-states(cb=TRUE) #State shapefiles
StateOutline<-list() #list to hold state oulines
for(i in seq(1,56,1)){ #subset state shapes to indivigual states
  print(i)
  StateOutline[[i]]<-subset.data.frame(StatesSHP, as.numeric(StatesSHP$STATEFP)==i)
}
rm(StatesSHP) #drop unified state shapefile
#------

#create list of counties in each state
CountyOutline<-list() #create list to hold county shapes
for (i in seq(1,56,1)){ #pull county shapes for each state by FIPS code
  if (i %in% c(03,07,14,43,52)){ # list of non states if pulled it grabs *ALL* data, so exclude (159mb down to 29mb)
    print(paste('Do not pull County:', i))
    CountyOutline[[i]]<-0
  } else {
    print(i)
    CountyOutline[[i]]<-counties(i,cb=TRUE)
  }
}
#------

#Create list of block groups in each state
BGOutline<-list() #create list to hold Block Group shapes
for (i in seq(1,56,1)){ #pull county shapes for each state by FIPS code
  if (i %in% c(03,07,14,43,52)){ # list of non states, exclude
    print(paste('Do not pull BG:', i))
    CountyOutline[[i]]<-0
  } else {
    print(i)
    BGOutline[[i]]<-block_groups(i,cb=TRUE)
  }
}
#-----
saveRDS(BGOutline,"BGOutline.rds")



#strip extra data from shapefiles
for (i in seq(1,56,1)){ 
  if (i %in% c(3,7,14,43,52)){
    print(paste('Skip non-state:', i))
  } else {
    print(i)
    BGOutline[[i]]@data<-BGOutline[[i]]@data[,c(6,9,10)]
  }
}

#strip extra data from USDAData
USDADatastrip<-list()
for (i in seq(1,56,1)){ 
  if (i %in% c(3,7,14,43,52)){
    print(paste('Skip non-state:', i))
  } else {
    print(i)
    USDADatastrip[[i]]<-USDAData[[i]][,c(2,5:9)]
    colnames(USDADatastrip[[i]])<-c("GEOID","ELIGFY19","Tot18","ID18","Tot12","ID12")
  }
}


#geojoin shapefile and usda data
BGJoin<-list()
for (i in seq(1,56,1)){ 
  if (i %in% c(3,7,14,43,52)){
    print(paste('Skip non-state:', i))
  } else {
    print(i)
    BGJoin[[i]]<-geo_join(BGOutline[[i]], USDADatastrip[[i]], by_sp='GEOID', by_df='GEOID', how='left')
    BGJoin[[i]]$GEOID.1<-NULL
  }
}





#geojoin shapefile over50
for (i in seq(1,56,1)){ 
  if (i %in% c(3,7,14,43,52)){
    print(paste('Skip non-state:', i))
  } else {
    print(i)
    BGJoin[[i]]<-geo_join(BGJoin[[i]], Over50strip[[i]], by_sp='GEOID', by_df='GEOID', how='left')
    BGJoin[[i]]$GEOID.1<-NULL
  }
}

#convert NA's to zeros for use in operator subsetting
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    print(paste('cleaning data for:', r))
    BGJoin[[r]]@data$weightnum18[is.na(BGJoin[[r]]@data$weightnum18)]<-0
    BGJoin[[r]]@data$weightpov18[is.na(BGJoin[[r]]@data$weightpov18)]<-0
    BGJoin[[r]]@data$weightPerc18[is.na(BGJoin[[r]]@data$weightPerc18)]<-0
    BGJoin[[r]]@data$weightnum12[is.na(BGJoin[[r]]@data$weightnum12)]<-0
    BGJoin[[r]]@data$weightpov12[is.na(BGJoin[[r]]@data$weightpov12)]<-0
    BGJoin[[r]]@data$weightPerc12[is.na(BGJoin[[r]]@data$weightPerc12)]<-0
  }
}
#------ 


#separeate the results into yes/no/calc so they can be mapped separately
BGSubsets<-list()
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    YNC<-list()
    print(paste('Subsetting Yes/No/Calc for:', r))
    YNC[[1]]<-subset.data.frame(BGJoin[[r]],BGJoin[[r]]$ELIGFY19=="Yes")
    YNC[[2]]<-subset.data.frame(BGJoin[[r]],BGJoin[[r]]$ELIGFY19=="No" & BGJoin[[r]]$weightPerc18 < 50 & BGJoin[[r]]$weightPerc12 < 50)
    YNC[[3]]<-subset.data.frame(BGJoin[[r]],BGJoin[[r]]$ELIGFY19=="No" & BGJoin[[r]]$weightPerc18 >=50 | BGJoin[[r]]$weightPerc12 >=50)
    BGSubsets[[r]]<-YNC
    rm(YNC)
  }
}
#------


#Build leaflets for each state
BGLeaflets<-list()
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    print(paste('Building Leaflet for:', r))
    BGLeaflets[[r]] <- leaflet() %>% enableTileCaching() %>%
      addProviderTiles(providers$OpenStreetMap.BlackAndWhite, group = "Grey")
    
    BGLeaflets[[r]]<- BGLeaflets[[r]] %>% addPolygons(data=BGSubsets[[r]][[1]],
                                                      weight = 1,
                                                      fill = TRUE,
                                                      fillOpacity = 0.3,
                                                      fillColor = 'red',
                                                      stroke = TRUE,
                                                      color = "black",
                                                      group = 'plotT',
                                                      label = BGSubsets[[r]][[1]]$GEOID, #paste("FY19Elig:", BGSubsets[[r]][[1]]$ELIGFY19),
                                                      labelOptions = labelOptions(textsize = "16px"),
                                                      highlight = highlightOptions(
                                                        weight = 5,
                                                        color = "red",
                                                        fillOpacity = 0.7,
                                                        bringToFront = TRUE),
                                                      popup = paste('<b>FY19Elig</b>: Eligable<br>',
                                                                    '<b>GEOID:</b>', BGSubsets[[r]][[1]]$GEOID, '<br>',
                                                                    '<b>Identified Under 18:</b>', BGSubsets[[r]][[1]]$ID18, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[1]]$Tot18, '<br>',
                                                                    '<b>Identified Under 12:</b>', BGSubsets[[r]][[1]]$ID12, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[1]]$Tot12, '<br>'))
    BGLeaflets[[r]]<- BGLeaflets[[r]] %>% addPolygons(data=BGSubsets[[r]][[2]],
                                                      weight = 1,
                                                      fill = TRUE,
                                                      fillOpacity = 0.3,
                                                      fillColor = 'blue',
                                                      stroke = TRUE,
                                                      color = "black",
                                                      group = 'plotT',
                                                      label = BGSubsets[[r]][[2]]$GEOID, #paste("FY19Elig:", BGSubsets[[r]][[2]]$ELIGFY19),
                                                      labelOptions = labelOptions(textsize = "16px"),
                                                      highlight = highlightOptions(
                                                        weight = 5,
                                                        color = "red",
                                                        fillOpacity = 0.7,
                                                        bringToFront = TRUE),
                                                      popup = paste('<b>FY19Elig</b>: Not Eligable<br>',
                                                                    '<b>GEOID:</b>', BGSubsets[[r]][[2]]$GEOID, '<br>',
                                                                    '<b>Identified Under 18:</b>', BGSubsets[[r]][[2]]$ID18, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[2]]$Tot18, '<br>',
                                                                    '<b>Identified Under 12:</b>', BGSubsets[[r]][[2]]$ID12, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[2]]$Tot12, '<br>'))
    BGLeaflets[[r]]<- BGLeaflets[[r]] %>% addPolygons(data=BGSubsets[[r]][[3]],
                                                      weight = 1,
                                                      fill = TRUE,
                                                      fillOpacity = 0.3,
                                                      fillColor = 'yellow',
                                                      stroke = TRUE,
                                                      color = "black",
                                                      group = 'plotT',
                                                      label = BGSubsets[[r]][[3]]$GEOID, #paste("FY19Elig: Calc Eligable"),
                                                      labelOptions = labelOptions(textsize = "16px"),
                                                      highlight = highlightOptions(
                                                        weight = 5,
                                                        color = "red",
                                                        fillOpacity = 0.7,
                                                        bringToFront = TRUE),
                                                      popup = paste('<b>FY19Elig</b>: Calculation Eligable<br>',
                                                                    '<b>GEOID:</b>', BGSubsets[[r]][[3]]$GEOID, '<br>',
                                                                    '<b>Identified Under 18:</b>', BGSubsets[[r]][[3]]$ID18, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[3]]$Tot18, '<br>',
                                                                    '<b>Identified Under 12:</b>', BGSubsets[[r]][[3]]$ID12, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[3]]$Tot12, '<br><hr>',
                                                                    '<b><center>Calculation</center></b><hr>',
                                                                    '<table border="1">
                                                                      <tr>
                                                                      <th style="background-color:#c5d9d5;" align="center">GEOID</th>
                                                                      <th style="background-color:#c5d9d5;" align="center">Under 18</th>
                                                                      <th style="background-color:#c5d9d5;" align="center">18 pop</th>
                                                                      <th style="background-color:#c5d9d5;" align="center">Under 12</th>
                                                                      <th style="background-color:#c5d9d5;" align="center">12 pop</th>
                                                                      </tr>
                                                                      <tbody>
                                                                      <tr>
                                                                      <td align="center">',BGSubsets[[r]][[3]]$firstGEOID18,'</td>
                                                                      <td align="center">',BGSubsets[[r]][[3]]$B1_ID18,'</td>
                                                                      <td align="center">',BGSubsets[[r]][[3]]$B1_Tot18,'</td>
                                                                      <td align="center">',BGSubsets[[r]][[3]]$B1_ID12,'</td>
                                                                      <td align="center">',BGSubsets[[r]][[3]]$B1_Tot12,'</td>
                                                                      </tr>
                                                                      <tr>
                                                                      <td align="center">', BGSubsets[[r]][[3]]$secondGEOID18, '</td>
                                                                      <td align="center">', BGSubsets[[r]][[3]]$B2_ID18, '</td>
                                                                      <td align="center">', BGSubsets[[r]][[3]]$B2_Tot18, '</td>
                                                                      <td align="center">', BGSubsets[[r]][[3]]$B2_ID12, '</td>
                                                                      <td align="center">', BGSubsets[[r]][[3]]$B2_Tot12, '</td>
                                                                      </tr>
                                                                      <tr>
                                                                      <td style="background-color:#f7daad;" align="center">Totals</td>
                                                                      <td style="background-color:#f7daad;" align="center">', BGSubsets[[r]][[3]]$ID18+BGSubsets[[r]][[3]]$B1_ID18+BGSubsets[[r]][[3]]$B2_ID18, '</td>
                                                                      <td style="background-color:#f7daad;" align="center">', BGSubsets[[r]][[3]]$Tot18+BGSubsets[[r]][[3]]$B1_Tot18+BGSubsets[[r]][[3]]$B2_Tot18, '</td>
                                                                      <td style="background-color:#f7daad;" align="center">', BGSubsets[[r]][[3]]$ID12+BGSubsets[[r]][[3]]$B1_ID12+BGSubsets[[r]][[3]]$B2_ID12, '</td>
                                                                      <td style="background-color:#f7daad;" align="center">', BGSubsets[[r]][[3]]$Tot12+BGSubsets[[r]][[3]]$B1_Tot12+BGSubsets[[r]][[3]]$B2_Tot12, '</td>
                                                                      </tr>
                                                                      </tbody>
                                                                      </table>',
                                                                    '<b>Weighted % under 18:</b>', round(BGSubsets[[r]][[3]]$weightPerc18,2), '<br>',
                                                                    '<b>Weighted % under 12:</b>', round(BGSubsets[[r]][[3]]$weightPerc12,2), '<br>'))
    BGLeaflets[[r]] <- BGLeaflets[[r]] %>% addLegend(position = "bottomright", 
                                                     labels =  c("USDA Eligible","Calculated Eligible", "Ineligible"), 
                                                     colors =  c('red','yellow','blue'),
                                                     opacity = 0.5,
                                                     group = 'BlockG')
    
    BGLeaflets[[r]] <- BGLeaflets[[r]] %>% addSearchOSM(options = searchOptions(autoCollapse = TRUE, minLength = 2))
  }
}
#------
saveRDS(BGLeaflets, "BGLeaflets.rds") #save leaflet maps to prevent having to rerun




#export all data into folders
SumAll<-data.frame()
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    print(paste("Saving Map for:", r))
    mainDir <- getwd()
    subDir <- paste(mainDir,"/",gsub(" ","_",StateOutline[[r]]$NAME),sep="")
    if (file.exists(subDir)){
      setwd(file.path(subDir))
    } else {
      dir.create(file.path(subDir))
      setwd(file.path(subDir))
      
    }
    name<-paste(StateOutline[[r]]$NAME,"_FY19.html", sep="")
    name<-gsub(" ","_",name)
    htmlwidgets::saveWidget(BGLeaflets[[r]], name, selfcontained = TRUE)
    write.csv(BGJoin[[r]],paste(StateOutline[[r]]$NAME,".csv", sep = ""), row.names = FALSE)
    saveRDS(BGSubsets[[r]][[1]], "Elig_Yes.rds")
    saveRDS(BGSubsets[[r]][[2]], "Elig_No.rds")
    saveRDS(BGSubsets[[r]][[3]], "Elig_Calc.rds")
    Summary<-data.frame(
      Name=StateOutline[[r]]$NAME, #name of area
      AreasYes=length(BGSubsets[[r]][[1]]$GEOID), # number of USDA areas
      Under18Yes=sum(BGSubsets[[r]][[1]]$ID18), # Number of USDA identified under 18
      Under12Yes=sum(BGSubsets[[r]][[1]]$ID12), # Number of USDA identified under 12
      
      AreasCalc=length(BGSubsets[[r]][[3]]$GEOID), # number of calc eligable areas
      Under18Calc=sum(BGSubsets[[r]][[3]]$ID18), # Number of calc identified under 18
      Under12Calc=sum(BGSubsets[[r]][[3]]$ID12), # Number of calc identified under 12
      
      AreasDeltaPerc=(((length(BGSubsets[[r]][[1]]$GEOID)+length(BGSubsets[[r]][[3]]$GEOID))-length(BGSubsets[[r]][[1]]$GEOID))/length(BGSubsets[[r]][[1]]$GEOID))*100, # %change num of areas
      Under18DeltaPerc=(((sum(BGSubsets[[r]][[1]]$ID18)+sum(BGSubsets[[r]][[3]]$ID18))-sum(BGSubsets[[r]][[1]]$ID18))/sum(BGSubsets[[r]][[1]]$ID18))*100, # %change Ident under 18
      Under12DeltaPerc=(((sum(BGSubsets[[r]][[1]]$ID12)+sum(BGSubsets[[r]][[3]]$ID12))-sum(BGSubsets[[r]][[1]]$ID12))/sum(BGSubsets[[r]][[1]]$ID12))*100, # %change Ident under 12
      
      AreasNo=length(BGSubsets[[r]][[2]]$GEOID), #num of not eligable
      Under18No=sum(BGSubsets[[r]][[2]]$ID18) #num of under 18 still not served
    )
    SumAll<-rbind.data.frame(SumAll,Summary)
    write.csv(Summary,paste(StateOutline[[r]]$NAME,"_Summary.csv", sep = ""), row.names = FALSE)
    setwd(mainDir)
  }
}
#------
write.csv(SumAll, "Overall_Summary.csv", row.names = FALSE)




# test<-data.frame()
# #Create list of block groups in each state
# BGOutline<-list() #create list to hold Block Group shapes
# for (i in seq(1,56,1)){ #pull county shapes for each state by FIPS code
#   if (i %in% c(03,07,14,43,52)){ # list of non states, exclude
#     print(paste('Do not pull BG:', i))
#   } else {
#     print(i)
#     test<-rbind(test,Over50[[i]])
#   }
# }


#merge Combined=USDAData and Over50

#merge combined and BGoutline
#build maps
#build metrics








